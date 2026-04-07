"""ACF/DLW Production Function Estimator with Analytical Standard Errors.

Two-step estimator following Ackerberg, Caves & Frazer (2015, Econometrica):
  Stage 1: OLS polynomial -> phi-hat (purge measurement error)
  Stage 2: GMM minimizing E[xi*Z]=0 (recover PF parameters)

Follows the coding conventions of pyblp (Conlon & Gortmaker, 2020):
  - Configuration classes for optimization
  - Structured result objects
  - Module-level options
  - NumPy-style docstrings
  - Type aliases

References
----------
ACF (2015): Identification Properties of Recent PF Estimators, Econometrica.
DLW (2012): Markups and Firm-Level Export Status, AER.
ACH (2012): Asymptotic Theory for Nonparametric Two-Stage Estimators, REStat.
DGM (2026): Hitchhiker's Guide to Markup Estimation, Econometrica.
CWDL (2015): Reallocation and Technology, AER (survival correction).
CWDL (2020): PF Estimation with Measurement Error in Capital (IV correction).
ADL (2024): Production Function Identification Under Imperfect Competition, CEPR DP.
ADL (2024b): QEC Model, AER P&P.
Hall (2018): New Evidence on the Markup of Prices over Marginal Costs, JEP.

Author: Marek Chadim (Yale, Tobin Center)
"""

from __future__ import annotations

import time
import warnings
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy.optimize import minimize

# type aliases (pyblp convention)
Array = Any
Options = Dict[str, Any]


# ========================================================================== #
#  Module-level options (pyblp convention: mutable module attributes)
# ========================================================================== #

class _Options:
    """Global configuration for ACF estimation."""
    digits: int = 6
    verbose: bool = True
    verbose_output: Callable = print
    dtype: type = np.float64
    finite_differences_epsilon: float = np.sqrt(np.finfo(np.float64).eps)
    pseudo_inverses: bool = True
    singular_tol: float = 1 / np.finfo(np.float64).eps
    flush_output: bool = False

options = _Options()


def _output(message: str) -> None:
    """Print a message if verbosity is on."""
    if options.verbose:
        options.verbose_output(message)
        if options.flush_output:
            import sys
            sys.stdout.flush()


# ========================================================================== #
#  Exceptions (pyblp convention: domain-specific errors)
# ========================================================================== #

class ACFError(Exception):
    """Base exception for ACF estimation errors."""

class ConvergenceError(ACFError):
    """Raised when the GMM optimizer fails to converge."""

class DataError(ACFError):
    """Raised when the input data has issues."""

class SingularMatrixError(ACFError):
    """Raised when a matrix inversion fails."""


# ========================================================================== #
#  Configuration classes (pyblp convention: Optimization, Formulation)
# ========================================================================== #

@dataclass
class Optimization:
    r"""Configuration for the GMM optimization routine.

    Parameters
    ----------
    method : str
        Optimization method. Supported:
            - ``'nelder-mead'`` : Derivative-free simplex (default, recommended
              for ACF). Uses a large initial simplex for global exploration
              followed by polishing restarts.
            - ``'bfgs'`` : Quasi-Newton with numerical gradients. Used as a
              polishing step after Nelder-Mead.
            - ``'nm+bfgs'`` : Hybrid Nelder-Mead then BFGS polish.
            - ``'fsolve'`` : Root-finding for just-identified systems (DGM).
              Solves m(β)=0 directly using :func:`scipy.optimize.fsolve`.
              More robust than minimization on flat criterion surfaces.
            - ``'basin_hopping'`` : Global optimizer with random perturbations
              (:func:`scipy.optimize.basinhopping`). Wraps NM local search.
            - ``'powell'`` : Derivative-free Powell's method.
    method_options : dict, optional
        Options passed to :func:`scipy.optimize.minimize`.
    simplex_delta : float
        Initial simplex perturbation for Nelder-Mead (default 0.1). The
        initial simplex is constructed by perturbing each starting parameter
        by this amount. Larger values explore more of the parameter space.
    n_restarts : int
        Number of Nelder-Mead polishing restarts with small simplex
        (delta=1e-5). Default 3.
    maxiter : int
        Maximum iterations per optimization pass (default 10000).

    Examples
    --------
    Default (Nelder-Mead with BFGS polish)::

        >>> opt = Optimization()

    Pure BFGS::

        >>> opt = Optimization(method='bfgs')

    Custom simplex::

        >>> opt = Optimization(simplex_delta=0.5, n_restarts=5)
    """
    method: str = 'nm+bfgs'
    method_options: Optional[Options] = None
    simplex_delta: float = 0.1
    n_restarts: int = 3
    maxiter: int = 10000


@dataclass
class Formulation:
    r"""Configuration for the production function specification.

    Parameters
    ----------
    spec : str
        ``'cd'`` for Cobb-Douglas, ``'tl'`` for translog.
    poly_order : int
        Polynomial order for first-stage OLS (default 3).
    ar_order : int
        Polynomial order of the first-order Markov transition g(ω_{t-1}).
        ω_t = g(ω_{t-1}) + ξ_t where g is approximated by a polynomial:
        1 = linear (g = a + ρω), 2 = quadratic (+ γω²),
        3 = cubic (+ δω³). Default 1 (linear), following DGM (2026)
        and CWDL (2015). DGM offer quadratic as a robustness option.
    pp_in_markov : bool
        Include lagged procurement dummy in the Markov process (default True).
    pp_interactions : bool
        Interact first-stage polynomial with pp_dummy (default True).
    year_fe : bool
        Include year fixed effects in first stage (default True).
    nace2_fe : bool
        Include nace2 fixed effects in first stage (default False).
    first_stage_controls : list of str
        Additional control variables in first-stage OLS (default []).
        DGM (2026) include firm market share (``salessharefirm``);
        CWDL (2015) include employment (``lne``). These enter linearly,
        not in the polynomial expansion.
    """
    spec: str = 'tl'
    poly_order: int = 3
    ar_order: int = 1
    pp_in_markov: bool = True
    pp_interactions: bool = True
    year_fe: bool = True
    nace2_fe: bool = False
    variable_input: str = 'cogs'
    """Column name for the variable (flexible) input. Default ``'cogs'``.

    Can be set to ``'ii'`` (intermediate inputs), ``'cogs'`` (cost of goods
    sold), or any log-deflated expenditure column present in the data.
    The chosen column is mapped internally to ``'cogs'`` so the rest of the
    estimator works unchanged. The expenditure share α̂ = exp(v)/exp(go)
    uses this variable, affecting markup levels.
    """
    additional_inputs: List[str] = field(default_factory=list)
    """Additional production function inputs beyond k and the variable input.

    E.g., ``['o']`` for overhead/services, ``['o', 'le']`` for overhead
    and labor. These enter the first-stage polynomial, the second-stage
    GMM (as additional beta parameters), and the instrument set (lagged).
    Following DGM (2026) 4-input framework: s = f(v, k, m, o).

    Only supported for ``spec='cd'`` currently. Each additional input adds
    one beta parameter (CD) and one instrument (lagged value).
    """
    first_stage_controls: List[str] = field(default_factory=list)
    overidentify: bool = False
    """Add deeper lags (L.k, L2.cogs) as extra instruments for translog.

    Creates overidentification that improves finite-sample identification
    of linear translog terms (Kim, Luo & Su 2019, JAE). Enables Hansen J
    test. Only affects translog specification.
    """
    optimal_instruments: str = ''
    """Chamberlain (1987) / Amemiya (1977) optimal instruments strategy.

    ``''`` (empty): disabled (default).
    ``'replace'``: replace Z with K_beta projected Jacobians.
        Produces exactly K_beta optimal instruments Z* = E[∂ξ/∂β | Z],
        estimated via sieve projection of the Jacobian onto all exogenous
        variables. The system stays just-identified but with the
        asymptotically efficient combination of instrument information.
        This is the approach used by pyblp (Conlon & Gortmaker 2020).
        Uses fsolve (root-finding) for exact solution.
    ``'augment'``: augment Z with projected Jacobians.
        Stacks K_beta Chamberlain columns alongside standard Z, creating
        overidentification. Hansen J is available but will reject if the
        Markov process is misspecified, since the Jacobian is contaminated.

    .. warning::
        Empirically, both modes perform poorly for ACF estimation because
        the standard instrument set (k, L.cogs) is too small for the sieve
        projection to add value. Chamberlain optimal instruments help when
        many instruments are compressed into K_beta efficient combinations
        (as in BLP with many excluded instruments). For ACF, use oligopoly
        instruments (ADL 2024) as standard instruments in Z instead.

    Step 1: Estimate with standard instruments to get consistent β̂.
    Step 2: Compute observation-level Jacobian ∂ξ/∂β, project onto
    sieve basis of Z, and re-estimate with optimal instruments.
    """
    optimal_instruments_sieve_order: int = 2
    """Polynomial order for the sieve basis in optimal instrument
    construction. Higher orders capture more nonlinearity in E[∂ξ/∂β|Z]
    but risk overfitting. Default 2 (quadratic in Z).
    """


@dataclass
class CWDLExtensions:
    r"""Configuration for CWDL extensions to ACF.

    Parameters
    ----------
    iv_capital : bool
        Use lagged investment as IV for capital (CWDL 2020). Requires
        ``'investment'`` column. Default False.
    survival_correction : bool
        Include predicted survival probability in Markov process
        (CWDL 2015, AER). Requires ``'survival'`` column. Default False.
    weighting : str
        GMM weighting matrix: ``'optimal'`` for inv(Z'Z), ``'identity'``
        for I. Default ``'optimal'``.
    markov_controls : list of str
        Additional variables to include in the Markov process. Lagged
        values will be used. Default [].
    markov_interactions : bool
        Interact Markov controls with lagged omega (CWDL 2015). Default
        False.
    """
    iv_capital: bool = False
    survival_correction: bool = False
    weighting: str = 'optimal'
    markov_controls: List[str] = field(default_factory=list)
    markov_interactions: bool = False


@dataclass
class ImperfectCompetition:
    r"""Configuration for imperfect competition in PF estimation (ADL 2024).

    Implements the sufficient statistic approach from Ackerberg & De Loecker
    (2024, CEPR DP 19640). Under imperfect competition, input demand depends
    on competitors' actions:

    .. math::

        V_{jt} = V(\omega_{jt}, F_{jt}, Q_{-jt}, Z_t)

    where :math:`Q_{-jt}` is the sufficient statistic for competition. For
    homogeneous goods Cournot, :math:`Q_{-jt} = \sum_{k \neq j} Q_{kt}`.
    For the QEC model (heterogeneous output, ADL 2024b AER P&P), deflated
    revenues :math:`\tilde{R}_{-jt}` replace quantities.

    The sufficient statistic enters:
    1. First-stage polynomial (control function inversion)
    2. Markov process for productivity (competitive environment is a state)
    3. Oligopoly instruments in the second-stage moment conditions

    Parameters
    ----------
    enabled : bool
        Enable imperfect competition correction. Default False.
    sufficient_statistic : str
        Column name for the sufficient statistic (e.g., ``'comp_output'``
        for :math:`Q_{-jt}` or ``'comp_revenue'`` for :math:`\tilde{R}_{-jt}`).
        Must be present in data.
    oligopoly_instruments : list of str
        Column names for oligopoly instruments. Typically summary statistics
        of competitors: ``['comp_k_mean', 'comp_omega_lag_mean']`` following
        ADL Section 4.2 (use moments to avoid many-instruments bias).
    in_first_stage : bool
        Include sufficient statistic in first-stage polynomial. Default True.
        This corrects the control function inversion (ADL eq. 34).
    in_markov : bool
        Include sufficient statistic in the Markov process. Default True.
        Competition affects the evolution of productivity.
    ss_interactions : bool
        Interact sufficient statistic with inputs in first-stage polynomial.
        Default True (needed for flexible control function approximation).
    """
    enabled: bool = False
    sufficient_statistic: str = ''
    oligopoly_instruments: List[str] = field(default_factory=list)
    in_first_stage: bool = True
    in_markov: bool = True
    ss_interactions: bool = True


# ========================================================================== #
#  Result classes (pyblp convention: structured immutable results)
# ========================================================================== #

@dataclass
class OptimizationDiagnostics:
    """Diagnostics for each optimization start.

    Attributes
    ----------
    criterion : float
        GMM criterion value at the solution.
    converged : bool
        Whether the optimizer converged.
    n_iterations : int
        Number of function evaluations.
    message : str
        Optimizer status message.
    method : str
        Which method produced this result.
    """
    criterion: float
    converged: bool
    n_iterations: int
    message: str
    method: str


@dataclass
class ACFResults:
    """Results from ACF production function estimation.

    Attributes
    ----------
    betas : Array
        Estimated production function parameters.
    beta_names : list of str
        Names corresponding to each element of ``betas``.
    se : Array
        Analytical standard errors (ACH 2012, clustered by firm).
    vcov : Array
        Variance-covariance matrix of ``betas``.
    spec : str
        Specification (``'cd'`` or ``'tl'``).
    first_stage_r2 : float
        R-squared from first-stage OLS.
    gmm_criterion : float
        GMM objective value at the solution.
    n_obs : int
        Number of observations in GMM sample.
    n_clusters : int
        Number of firm clusters.
    data : pd.DataFrame
        Estimation sample with markups, SEs, omega, etc.
    optimization_diagnostics : list of OptimizationDiagnostics
        Per-start diagnostics for convergence assessment.
    elapsed_time : float
        Total estimation time in seconds.
    """
    betas: Array
    beta_names: List[str]
    se: Array
    vcov: Array
    spec: str
    first_stage_r2: float
    gmm_criterion: float
    n_obs: int
    n_clusters: int
    data: pd.DataFrame
    optimization_diagnostics: List[OptimizationDiagnostics]
    elapsed_time: float
    hansen_j: Optional[float] = None
    hansen_j_pvalue: Optional[float] = None
    n_overid: int = 0

    def __repr__(self) -> str:
        lines = [
            f"ACFResults(spec='{self.spec}', N={self.n_obs}, "
            f"N_clusters={self.n_clusters})",
            f"  First-stage R² = {self.first_stage_r2:.4f}",
            f"  GMM criterion  = {self.gmm_criterion:.6e}",
            f"  Elapsed time   = {self.elapsed_time:.1f}s",
            "",
            f"  {'Param':8s}  {'Coef':>10s}  {'SE':>10s}  {'t':>8s}",
            f"  {'-'*40}",
        ]
        for name, b, s in zip(self.beta_names, self.betas, self.se):
            t = b / s if s > 0 else np.inf
            lines.append(f"  {name:8s}  {b:10.{options.digits}f}  "
                         f"{s:10.{options.digits}f}  {t:8.2f}")
        mu = self.data['markup']
        valid = np.isfinite(mu) & (mu > 0)
        mv = mu[valid]
        lines.extend([
            "",
            f"  Markups: mean={mv.mean():.4f}, sd={mv.std():.4f}, "
            f"median={mv.median():.4f}",
            f"  p10={mv.quantile(0.1):.4f}, p90={mv.quantile(0.9):.4f}, "
            f"frac<1={100*(mv < 1).mean():.1f}%",
        ])
        return "\n".join(lines)

    @property
    def markup_stats(self) -> Dict[str, float]:
        """Summary statistics for markups."""
        mu = self.data['markup']
        valid = np.isfinite(mu) & (mu > 0)
        mv = mu[valid]
        return {
            'mean': float(mv.mean()), 'sd': float(mv.std()),
            'p10': float(mv.quantile(0.1)), 'p25': float(mv.quantile(0.25)),
            'p50': float(mv.quantile(0.5)), 'p75': float(mv.quantile(0.75)),
            'p90': float(mv.quantile(0.9)), 'frac_below_1': float((mv < 1).mean()),
            'mean_se': float(self.data.loc[valid, 'se_markup'].mean()),
        }


# ========================================================================== #
#  Core estimator
# ========================================================================== #

class ACFEstimator:
    r"""ACF/DLW two-step production function estimator.

    Estimates a production function via the Ackerberg, Caves & Frazer (2015)
    control function approach with GMM. Computes firm-level markups and
    analytical standard errors following Ackerberg, Chen & Hahn (2012).

    Parameters
    ----------
    data : pd.DataFrame
        Panel data with columns: ``id``, ``year``, ``go``, ``k``, ``cogs``,
        ``pp_dummy``, ``nace2``. Variables ``go``, ``k``, ``cogs`` should be
        in logs.
    formulation : Formulation, optional
        Production function specification. Default: translog.
    optimization : Optimization, optional
        Optimizer configuration. Default: NM+BFGS hybrid.
    extensions : CWDLExtensions, optional
        CWDL extensions (IV capital, survival correction, etc.).
    imperfect_competition : ImperfectCompetition, optional
        ADL (2024) imperfect competition configuration.
    n_starts : int
        Number of starting points for multi-start optimization (default 3).

    Returns
    -------
    ACFResults
        Call :meth:`solve` to run the estimation and return results.

    Examples
    --------
    Basic translog estimation::

        >>> est = ACFEstimator(data)
        >>> results = est.solve()
        >>> print(results)

    Cobb-Douglas with CWDL extensions::

        >>> est = ACFEstimator(
        ...     data,
        ...     formulation=Formulation(spec='cd'),
        ...     extensions=CWDLExtensions(iv_capital=True),
        ... )
        >>> results = est.solve()
    """

    def __init__(
        self,
        data: pd.DataFrame,
        formulation: Optional[Formulation] = None,
        optimization: Optional[Optimization] = None,
        extensions: Optional[CWDLExtensions] = None,
        imperfect_competition: Optional[ImperfectCompetition] = None,
        n_starts: int = 3,
    ) -> None:
        self._raw_data = data.copy()
        self._formulation = formulation or Formulation()
        self._optimization = optimization or Optimization()
        self._extensions = extensions or CWDLExtensions()
        self._ic = imperfect_competition or ImperfectCompetition()
        self._n_starts = n_starts

        # variable input mapping: user-specified column → internal 'cogs'
        vi = self._formulation.variable_input
        self._variable_input_name = vi
        if vi != 'cogs':
            if vi not in data.columns:
                raise DataError(f"Variable input '{vi}' not in data columns")
            self._raw_data['cogs'] = self._raw_data[vi]
            _output(f"  Variable input: '{vi}' → mapped to internal 'cogs'")

        # validate required columns
        required = ['id', 'year', 'go', 'k', 'cogs']
        for ai in self._formulation.additional_inputs:
            required.append(ai)
        missing = [c for c in required if c not in data.columns]
        if missing:
            raise DataError(f"Missing required columns: {missing}")

        # validate imperfect competition columns
        if self._ic.enabled:
            if not self._ic.sufficient_statistic:
                raise DataError("ImperfectCompetition.sufficient_statistic "
                                "must be specified when enabled=True")
            ic_cols = [self._ic.sufficient_statistic] + \
                self._ic.oligopoly_instruments
            ic_missing = [c for c in ic_cols if c not in data.columns]
            if ic_missing:
                raise DataError(f"Missing imperfect competition columns: "
                                f"{ic_missing}")

    def solve(self) -> ACFResults:
        """Run the full estimation pipeline and return results.

        Returns
        -------
        ACFResults
            Production function estimates, markups, and diagnostics.
        """
        t0 = time.time()

        ic_tag = " [ADL IC]" if self._ic.enabled else ""
        _output(f"ACF Estimation: spec={self._formulation.spec.upper()}, "
                f"method={self._optimization.method}{ic_tag}")

        # Stage 1: first-stage OLS
        self._prepare_data()
        r2 = self._first_stage()

        # Prepare GMM matrices
        self._prepare_gmm()
        _output(f"  GMM sample: N={self._N}, N_clusters={self._N_clusters}")

        # Stage 2: multi-start GMM
        betas, criterion, diagnostics = self._multi_start_gmm()

        # Chamberlain (1987) optimal instruments: two-step procedure
        oi_mode = self._formulation.optimal_instruments
        if oi_mode in ('replace', 'augment'):
            _output(f"\n  --- Chamberlain Optimal Instruments "
                    f"(mode='{oi_mode}') ---")
            Z_proj = self._construct_optimal_instruments(betas)

            # save standard Z for diagnostics
            self._Z_standard = self._Z.copy()

            if oi_mode == 'replace':
                self._Z = Z_proj
                _output(f"  Replaced Z: {Z_proj.shape[1]} optimal instruments "
                        f"(= K_beta, just-identified)")
            else:  # augment
                self._Z = np.hstack([self._Z, Z_proj])
                n_new = self._Z.shape[1]
                _output(f"  Augmented Z: {n_new} instruments "
                        f"(standard + {Z_proj.shape[1]} Chamberlain)")

            # update weighting matrix for new Z
            K_z = self._Z.shape[1]
            if self._extensions.weighting == 'identity':
                self._W = np.eye(K_z, dtype=options.dtype)
            else:
                ZtZ = self._Z.T @ self._Z
                if options.pseudo_inverses:
                    self._W = np.linalg.pinv(ZtZ)
                else:
                    self._W = np.linalg.inv(ZtZ)

            # re-estimate with optimal instruments
            # For replace mode (just-identified), use fsolve for exact root
            if oi_mode == 'replace':
                saved_optimization = self._optimization
                self._optimization = Optimization(method='fsolve')
                # use standard-Z solution as starting point
                self._beta_init = betas.copy()
            betas_opt, criterion_opt, diag_opt = self._multi_start_gmm()
            if oi_mode == 'replace':
                self._optimization = saved_optimization
            _output(f"  Optimal instruments GMM Q: {criterion_opt:.6e} "
                    f"(standard: {criterion:.6e})")
            betas = betas_opt
            criterion = criterion_opt
            diagnostics = diagnostics + diag_opt

        # Analytical SEs (ACH 2012)
        vcov, se = self._analytical_vcov(betas)

        # Markups and delta-method SEs
        data_out = self._compute_markups(betas, vcov)

        # Hansen J test (overidentification)
        n_moments = self._Z.shape[1]
        n_params = len(betas)
        n_overid = n_moments - n_params
        hansen_j = None
        hansen_j_pvalue = None
        if n_overid > 0:
            hansen_j, hansen_j_pvalue = self._hansen_j_test(betas)
            _output(f"\n  Hansen J test: χ²({n_overid}) = {hansen_j:.3f}, "
                    f"p = {hansen_j_pvalue:.3f}")

        elapsed = time.time() - t0

        results = ACFResults(
            betas=betas,
            beta_names=self._beta_names,
            se=se,
            vcov=vcov,
            spec=self._formulation.spec,
            first_stage_r2=r2,
            gmm_criterion=criterion,
            n_obs=self._N,
            n_clusters=self._N_clusters,
            data=data_out,
            optimization_diagnostics=diagnostics,
            elapsed_time=elapsed,
            hansen_j=hansen_j,
            hansen_j_pvalue=hansen_j_pvalue,
            n_overid=n_overid,
        )

        _output(f"\n{results}")
        return results

    # ---------------------------------------------------------------------- #
    #  Stage 1: Data preparation and first-stage OLS
    # ---------------------------------------------------------------------- #

    def _prepare_data(self) -> None:
        """Generate polynomial terms for first-stage OLS."""
        df = self._raw_data.sort_values(['id', 'year']).reset_index(drop=True)
        M = self._formulation.poly_order

        for i in range(2, M + 1):
            df[f'k{i}'] = df['k'] ** i
            df[f'cogs{i}'] = df['cogs'] ** i
        for i in range(1, M + 1):
            for j in range(1, M + 1):
                if i + j <= M and not (i == 1 and j == 0) \
                        and not (i == 0 and j == 1):
                    df[f'k{i}cogs{j}'] = df['k'] ** i * df['cogs'] ** j

        # polynomial terms for additional inputs (CD only, linear in GMM
        # but full polynomial in first stage for control function)
        for ai in self._formulation.additional_inputs:
            if ai not in df.columns:
                continue
            for i in range(2, M + 1):
                df[f'{ai}{i}'] = df[ai] ** i
            # cross terms with k and cogs
            for base in ['k', 'cogs']:
                for pk in range(1, M):
                    for pa in range(1, M):
                        if pk + pa <= M:
                            cname = f'{base}{pk}{ai}{pa}'
                            df[cname] = df[base] ** pk * df[ai] ** pa
            # cross terms between additional inputs
            for ai2 in self._formulation.additional_inputs:
                if ai2 <= ai or ai2 not in df.columns:
                    continue
                for p1 in range(1, M):
                    for p2 in range(1, M):
                        if p1 + p2 <= M:
                            cname = f'{ai}{p1}{ai2}{p2}'
                            df[cname] = df[ai] ** p1 * df[ai2] ** p2

        self._full_data = df

    def _first_stage(self) -> float:
        """OLS: go ~ polynomial(k, cogs) [x pp_dummy] + FEs -> phi-hat."""
        df = self._full_data.copy()
        M = self._formulation.poly_order
        form = self._formulation

        # build regressor list
        poly_vars = ['k', 'cogs']
        for i in range(2, M + 1):
            poly_vars.extend([f'k{i}', f'cogs{i}'])
        for i in range(1, M + 1):
            for j in range(1, M + 1):
                if i + j <= M:
                    name = f'k{i}cogs{j}'
                    if name in df.columns and name not in poly_vars:
                        poly_vars.append(name)

        # additional inputs and their polynomial/cross terms
        for ai in form.additional_inputs:
            if ai not in df.columns:
                continue
            poly_vars.append(ai)
            for i in range(2, M + 1):
                pname = f'{ai}{i}'
                if pname in df.columns:
                    poly_vars.append(pname)
            # cross terms
            for base in ['k', 'cogs']:
                for pk in range(1, M):
                    for pa in range(1, M):
                        if pk + pa <= M:
                            cname = f'{base}{pk}{ai}{pa}'
                            if cname in df.columns and cname not in poly_vars:
                                poly_vars.append(cname)
            for ai2 in form.additional_inputs:
                if ai2 <= ai or ai2 not in df.columns:
                    continue
                for p1 in range(1, M):
                    for p2 in range(1, M):
                        if p1 + p2 <= M:
                            cname = f'{ai}{p1}{ai2}{p2}'
                            if cname in df.columns and cname not in poly_vars:
                                poly_vars.append(cname)

        # imperfect competition: add sufficient statistic to polynomial
        # ADL (2024) eq. 34: q_jt = v^{-1}(v_jt, f_jt, Q_{-jt}, z_t) + eps
        ic = self._ic
        if ic.enabled and ic.in_first_stage:
            ss = ic.sufficient_statistic
            if ss in df.columns:
                poly_vars.append(ss)
                # higher powers of ss
                for p in range(2, M + 1):
                    pname = f'{ss}{p}'
                    df[pname] = df[ss] ** p
                    poly_vars.append(pname)
                # interactions: ss x k, ss x cogs (and higher if M > 2)
                if ic.ss_interactions:
                    for base in ['k', 'cogs']:
                        for pk in range(1, M):
                            for ps in range(1, M):
                                if pk + ps <= M:
                                    iname = f'{base}{pk}_{ss}{ps}'
                                    df[iname] = df[base] ** pk * df[ss] ** ps
                                    if iname not in poly_vars:
                                        poly_vars.append(iname)
                _output(f"  IC first stage: added {ss} + interactions "
                        f"({len(poly_vars)} poly vars)")

        rhs_vars = ['_const'] + poly_vars[:]
        df['_const'] = 1.0

        # procurement interactions
        if form.pp_interactions and 'pp_dummy' in df.columns:
            for v in poly_vars:
                iname = f'{v}_x_pp'
                df[iname] = df[v] * df['pp_dummy']
                rhs_vars.append(iname)

        # year FE
        if form.year_fe:
            years = sorted(df['year'].unique())
            for y in years[1:]:
                dname = f'_yr{y}'
                df[dname] = (df['year'] == y).astype(options.dtype)
                rhs_vars.append(dname)

        # nace2 FE
        if form.nace2_fe and 'nace2' in df.columns:
            naces = sorted(df['nace2'].unique())
            for n in naces[1:]:
                dname = f'_nace{n}'
                df[dname] = (df['nace2'] == n).astype(options.dtype)
                rhs_vars.append(dname)

        # year x nace2 interaction FE
        if form.year_fe and form.nace2_fe and 'nace2' in df.columns:
            years = sorted(df['year'].unique())
            naces = sorted(df['nace2'].unique())
            for y in years[1:]:
                for n in naces[1:]:
                    dname = f'_yr{y}_nace{n}'
                    df[dname] = ((df['year'] == y) & (df['nace2'] == n)
                                 ).astype(options.dtype)
                    rhs_vars.append(dname)

        # additional first-stage controls (DGM: salessharefirm; CWDL: lne)
        for ctrl in form.first_stage_controls:
            if ctrl in df.columns:
                rhs_vars.append(ctrl)

        # estimate
        mask = df[['go'] + rhs_vars].notna().all(axis=1)
        Y = df.loc[mask, 'go'].values.astype(options.dtype)
        X = df.loc[mask, rhs_vars].values.astype(options.dtype)

        betas_ols = np.linalg.lstsq(X, Y, rcond=None)[0]
        phi = X @ betas_ols
        epsilon = Y - phi

        df['phi'] = np.nan
        df.loc[mask, 'phi'] = phi
        df['epsilon'] = np.nan
        df.loc[mask, 'epsilon'] = epsilon
        self._full_data = df

        r2 = 1 - np.var(epsilon) / np.var(Y)
        _output(f"  First-stage R² = {r2:.4f}")
        return float(r2)

    # ---------------------------------------------------------------------- #
    #  GMM matrix construction
    # ---------------------------------------------------------------------- #

    def _prepare_gmm(self) -> None:
        """Build X, X_lag, Z, W matrices for second-stage GMM."""
        df = self._full_data.copy()
        ext = self._extensions
        form = self._formulation

        # additional inputs list (used throughout)
        ai_list = [ai for ai in form.additional_inputs if ai in df.columns]

        # corrected expenditure share
        df['y_c'] = df['go'] - df['epsilon']
        df['Y_c'] = np.exp(df['y_c'])
        df['alphahat'] = np.exp(df['cogs']) / df['Y_c']

        # lags (L1 and optionally L2 for overidentification)
        df = df.sort_values(['id', 'year'])
        lag_vars = ['phi', 'pp_dummy', 'k', 'cogs']
        if ext.iv_capital and 'investment' in df.columns:
            lag_vars.append('investment')
        for mc in ext.markov_controls:
            if mc not in lag_vars and mc in df.columns:
                lag_vars.append(mc)
        # imperfect competition: lag sufficient statistic and oligo instruments
        ic = self._ic
        if ic.enabled:
            ss = ic.sufficient_statistic
            if ss not in lag_vars and ss in df.columns:
                lag_vars.append(ss)
            for oi in ic.oligopoly_instruments:
                if oi not in lag_vars and oi in df.columns:
                    lag_vars.append(oi)
        for v in lag_vars:
            df[f'L{v}'] = df.groupby('id')[v].shift(1)

        # deeper lags for overidentification (Kim, Luo & Su 2019)
        if form.overidentify and form.spec == 'tl':
            df['L2cogs'] = df.groupby('id')['cogs'].shift(2)

        if ext.survival_correction and 'survival' in df.columns:
            df['Lsurvival'] = df.groupby('id')['survival'].shift(1)

        # lags for additional inputs
        for ai in form.additional_inputs:
            if ai in df.columns and ai not in lag_vars:
                lag_vars.append(ai)
                df[f'L{ai}'] = df.groupby('id')[ai].shift(1)

        # translog terms (always generated, used only if spec='tl')
        df['k2'] = df['k'] ** 2
        df['cogs2'] = df['cogs'] ** 2
        df['kcogs'] = df['k'] * df['cogs']
        df['Lk2'] = df['Lk'] ** 2
        df['Lcogs2'] = df['Lcogs'] ** 2
        df['LkLcogs'] = df['Lk'] * df['Lcogs']
        df['kLcogs'] = df['k'] * df['Lcogs']
        df['const'] = 1.0

        if ext.iv_capital and 'Linvestment' in df.columns:
            df['Linv2'] = df['Linvestment'] ** 2
            df['LinvLcogs'] = df['Linvestment'] * df['Lcogs']

        # drop missing
        needed = ['phi', 'Lphi', 'k', 'cogs', 'Lk', 'Lcogs',
                  'alphahat', 'id', 'year']
        for ai in ai_list:
            needed.extend([ai, f'L{ai}'])
        if form.overidentify and form.spec == 'tl' and 'L2cogs' in df.columns:
            needed.append('L2cogs')
        if form.pp_in_markov:
            needed.append('Lpp_dummy')
        if ext.iv_capital and 'Linvestment' in df.columns:
            needed.append('Linvestment')
        if ext.survival_correction and 'Lsurvival' in df.columns:
            needed.append('Lsurvival')
        for mc in ext.markov_controls:
            lmc = f'L{mc}'
            if lmc in df.columns:
                needed.append(lmc)
        # imperfect competition: need current and lagged ss + oligo instruments
        if ic.enabled:
            ss = ic.sufficient_statistic
            if ss in df.columns:
                needed.append(ss)
            lss = f'L{ss}'
            if lss in df.columns:
                needed.append(lss)
            for oi in ic.oligopoly_instruments:
                if oi in df.columns:
                    needed.append(oi)
        mask = df[needed].notna().all(axis=1)
        self._est_data = df[mask].copy().reset_index(drop=True)
        self._N = len(self._est_data)

        # variable definitions by spec
        if form.spec == 'cd':
            self._beta_names = ['const', 'k', 'cogs'] + ai_list
            x_vars = ['const', 'k', 'cogs'] + ai_list
            x_lag_vars = ['const', 'Lk', 'Lcogs'] + [f'L{ai}' for ai in ai_list]
            z_vars = ['const', 'k', 'Lcogs'] + [f'L{ai}' for ai in ai_list]
            if ai_list:
                _output(f"  Additional inputs (CD): {ai_list} "
                        f"({len(self._beta_names)} params, "
                        f"{len(z_vars)} instruments)")
        else:
            if ai_list:
                raise NotImplementedError(
                    "additional_inputs not yet supported for translog. "
                    "Use spec='cd' or extend _prepare_gmm for TL cross terms.")
            self._beta_names = ['const', 'k', 'cogs', 'k2', 'cogs2', 'kcogs']
            x_vars = ['const', 'k', 'cogs', 'k2', 'cogs2', 'kcogs']
            x_lag_vars = ['const', 'Lk', 'Lcogs', 'Lk2', 'Lcogs2', 'LkLcogs']
            z_vars = ['const', 'k', 'Lcogs', 'k2', 'Lcogs2', 'kLcogs']
            # Overidentification: add Lk and L2cogs (Kim, Luo & Su 2019)
            if form.overidentify:
                z_extra = ['Lk']
                if 'L2cogs' in self._est_data.columns:
                    z_extra.append('L2cogs')
                z_vars = z_vars + z_extra
                _output(f"  Overidentified TL: {len(z_vars)} instruments "
                        f"for {len(x_vars)} params "
                        f"(+{len(z_extra)}: {z_extra})")

        # IV capital: replace k instruments with lagged investment
        if ext.iv_capital and 'Linvestment' in self._est_data.columns:
            if form.spec == 'cd':
                z_vars = ['const', 'Linvestment', 'Lcogs']
            else:
                z_vars = ['const', 'Linvestment', 'Lcogs',
                          'Linv2', 'Lcogs2', 'LinvLcogs']
            _output("  IV capital: using lagged investment as instrument")

        # ADL (2024) oligopoly instruments: f_{-jt}, omega_{-jt-1}
        # These shift residual demand and are orthogonal to xi_jt by timing.
        if ic.enabled and ic.oligopoly_instruments:
            z_oligo = [oi for oi in ic.oligopoly_instruments
                       if oi in self._est_data.columns]
            if z_oligo:
                z_vars = z_vars + z_oligo
                _output(f"  ADL oligopoly instruments: +{len(z_oligo)} "
                        f"({z_oligo})")

        d = self._est_data
        self._X = d[x_vars].values.astype(options.dtype)
        self._X_lag = d[x_lag_vars].values.astype(options.dtype)
        self._Z = d[z_vars].values.astype(options.dtype)
        self._PHI = d['phi'].values.astype(options.dtype)
        self._PHI_LAG = d['Lphi'].values.astype(options.dtype)
        self._C = d['const'].values.astype(options.dtype).reshape(-1, 1)

        # Markov controls
        self._markov_extras: List[Array] = []
        if form.pp_in_markov and 'Lpp_dummy' in d.columns:
            self._markov_extras.append(
                d['Lpp_dummy'].values.astype(options.dtype).reshape(-1, 1))
        if ext.survival_correction and 'Lsurvival' in d.columns:
            self._markov_extras.append(
                d['Lsurvival'].values.astype(options.dtype).reshape(-1, 1))
        for mc in ext.markov_controls:
            lmc = f'L{mc}'
            if lmc in d.columns:
                self._markov_extras.append(
                    d[lmc].values.astype(options.dtype).reshape(-1, 1))

        # ADL (2024): sufficient statistic in Markov process
        # Competition is a state variable affecting productivity evolution
        if ic.enabled and ic.in_markov:
            lss = f'L{ic.sufficient_statistic}'
            if lss in d.columns:
                self._markov_extras.append(
                    d[lss].values.astype(options.dtype).reshape(-1, 1))
                _output(f"  ADL Markov: added L.{ic.sufficient_statistic}")

        # weighting matrix
        if ext.weighting == 'identity':
            self._W = np.eye(self._Z.shape[1], dtype=options.dtype)
        else:
            ZtZ = self._Z.T @ self._Z
            cond = np.linalg.cond(ZtZ)
            if cond > options.singular_tol:
                warnings.warn(f"Z'Z condition number ({cond:.2e}) exceeds "
                              f"tolerance ({options.singular_tol:.2e})")
            if options.pseudo_inverses:
                self._W = np.linalg.pinv(ZtZ)
            else:
                self._W = np.linalg.inv(ZtZ)

        # starting values from OLS: phi ~ X
        self._beta_init = np.linalg.lstsq(self._X, self._PHI, rcond=None)[0]

        # cluster info
        self._cluster_ids = d['id'].values
        self._unique_clusters = np.unique(self._cluster_ids)
        self._N_clusters = len(self._unique_clusters)

    # ---------------------------------------------------------------------- #
    #  GMM criterion and innovation
    # ---------------------------------------------------------------------- #

    def _compute_xi(self, betas: Array) -> Array:
        """Compute productivity innovation xi for given betas.

        First-order Markov: omega = g(omega_lag, controls) + xi,
        where g is a polynomial in omega_lag (order set by ar_order)
        with controls optionally including PP_lag, survival, etc.
        """
        betas = np.asarray(betas, dtype=options.dtype).ravel()
        omega = self._PHI - self._X @ betas
        omega_lag = self._PHI_LAG - self._X_lag @ betas

        omega_lag_col = omega_lag.reshape(-1, 1)
        omega_lag_pol = [self._C, omega_lag_col]
        for p in range(2, self._formulation.ar_order + 1):
            omega_lag_pol.append((omega_lag ** p).reshape(-1, 1))

        for ctrl in self._markov_extras:
            omega_lag_pol.append(ctrl)
            if self._extensions.markov_interactions:
                omega_lag_pol.append(ctrl * omega_lag_col)

        omega_lag_pol = np.hstack(omega_lag_pol)
        g_b = np.linalg.lstsq(omega_lag_pol, omega, rcond=None)[0]
        return omega - omega_lag_pol @ g_b

    def _gmm_criterion(self, betas: Array) -> float:
        """GMM objective: (Z'xi)' W (Z'xi)."""
        xi = self._compute_xi(betas)
        moments = self._Z.T @ xi
        return float(moments.T @ self._W @ moments)

    # ---------------------------------------------------------------------- #
    #  Multi-start optimization
    # ---------------------------------------------------------------------- #

    def _run_single_optimization(
        self, x0: Array, label: str
    ) -> Tuple[Array, float, OptimizationDiagnostics]:
        """Run one optimization from starting values x0."""
        opt = self._optimization
        K = len(x0)

        if opt.method in ('nelder-mead', 'nm+bfgs'):
            # Phase 1: large simplex for global exploration
            simplex = np.zeros((K + 1, K), dtype=options.dtype)
            simplex[0] = x0
            for i in range(K):
                simplex[i + 1] = x0.copy()
                simplex[i + 1, i] += opt.simplex_delta

            res = minimize(
                self._gmm_criterion, x0, method='Nelder-Mead',
                options={
                    'xatol': 1e-10, 'fatol': 1e-10,
                    'maxiter': opt.maxiter, 'maxfev': opt.maxiter * 10,
                    'initial_simplex': simplex, 'adaptive': False,
                    **(opt.method_options or {}),
                },
            )
            p = res.x
            nfev = res.nfev

            # Phase 2: polishing restarts (small simplex)
            for _ in range(opt.n_restarts):
                simplex_small = np.zeros((K + 1, K), dtype=options.dtype)
                simplex_small[0] = p
                for i in range(K):
                    simplex_small[i + 1] = p.copy()
                    simplex_small[i + 1, i] += 1e-5

                res_polish = minimize(
                    self._gmm_criterion, p, method='Nelder-Mead',
                    options={
                        'xatol': 1e-10, 'fatol': 1e-10,
                        'maxiter': opt.maxiter, 'maxfev': opt.maxiter * 10,
                        'initial_simplex': simplex_small, 'adaptive': False,
                    },
                )
                nfev += res_polish.nfev
                if res_polish.fun < self._gmm_criterion(p):
                    p = res_polish.x

            crit_nm = self._gmm_criterion(p)
            method_used = 'nelder-mead'

            # Phase 3: BFGS polish (if nm+bfgs)
            if opt.method == 'nm+bfgs':
                try:
                    res_bfgs = minimize(
                        self._gmm_criterion, p, method='BFGS',
                        options={'gtol': 1e-10, 'maxiter': opt.maxiter},
                    )
                    nfev += res_bfgs.nfev
                    if res_bfgs.success and res_bfgs.fun < crit_nm:
                        p = res_bfgs.x
                        method_used = 'nm+bfgs'
                        _output(f"    {label}: BFGS improved "
                                f"{crit_nm:.8e} -> {res_bfgs.fun:.8e}")
                    else:
                        _output(f"    {label}: BFGS did not improve "
                                f"(NM: {crit_nm:.8e})")
                except Exception:
                    _output(f"    {label}: BFGS failed, keeping NM")

        elif opt.method == 'bfgs':
            res = minimize(
                self._gmm_criterion, x0, method='BFGS',
                options={
                    'gtol': 1e-10, 'maxiter': opt.maxiter,
                    **(opt.method_options or {}),
                },
            )
            p = res.x
            nfev = res.nfev
            method_used = 'bfgs'

        elif opt.method == 'powell':
            res = minimize(
                self._gmm_criterion, x0, method='Powell',
                options={
                    'ftol': 1e-10, 'maxiter': opt.maxiter,
                    **(opt.method_options or {}),
                },
            )
            p = res.x
            nfev = res.nfev
            method_used = 'powell'

        elif opt.method == 'fsolve':
            # Root-finding: solve m(β)=0 directly (DGM approach).
            # Only valid for just-identified systems.
            from scipy.optimize import fsolve as _fsolve

            def _moments_vec(betas):
                xi = self._compute_xi(betas)
                return (self._Z.T @ xi / self._N).ravel()

            n_moments = self._Z.shape[1]
            n_params = len(x0)
            if n_moments != n_params:
                raise ACFError(
                    f"fsolve requires just-identified system "
                    f"({n_moments} moments != {n_params} params). "
                    f"Use 'nm+bfgs' for overidentified systems.")
            p, info, ier, mesg = _fsolve(
                _moments_vec, x0, full_output=True, maxfev=opt.maxiter * 10)
            nfev = info['nfev']
            if ier != 1:
                _output(f"    {label}: fsolve warning: {mesg.strip()}")
            method_used = 'fsolve'

        elif opt.method == 'basin_hopping':
            # Global optimizer: random perturbation + local NM (DGM approach)
            from scipy.optimize import basinhopping

            minimizer_kwargs = {
                'method': 'Nelder-Mead',
                'options': {'xatol': 1e-10, 'fatol': 1e-10,
                            'maxiter': opt.maxiter},
            }
            res = basinhopping(
                self._gmm_criterion, x0, minimizer_kwargs=minimizer_kwargs,
                niter=100, T=1.0, stepsize=opt.simplex_delta,
                disp=False, seed=42)
            p = res.x
            nfev = res.nfev if hasattr(res, 'nfev') else -1
            method_used = 'basin_hopping'

        else:
            raise ACFError(f"Unknown optimization method: {opt.method}")

        final_crit = self._gmm_criterion(p)
        converged = final_crit < 1.0  # heuristic: reasonable criterion

        diag = OptimizationDiagnostics(
            criterion=final_crit,
            converged=converged,
            n_iterations=nfev,
            message=f"{label}: {method_used}",
            method=method_used,
        )

        return p, final_crit, diag

    def _multi_start_gmm(
        self,
    ) -> Tuple[Array, float, List[OptimizationDiagnostics]]:
        """Multi-start GMM optimization.

        Uses multiple starting points to escape local optima:
          1. OLS starting values
          2. Perturbed OLS (0.5x scale)
          3. Zero-intercept start (for translog)
        """
        diagnostics = []
        best_betas = None
        best_crit = np.inf

        # Start 1: OLS
        starts = [
            (self._beta_init.copy(), "Start 1 (OLS)"),
        ]

        # Start 2: perturbed OLS (scale non-constant params by 0.5)
        x0_half = self._beta_init.copy()
        x0_half[1:] *= 0.5
        starts.append((x0_half, "Start 2 (0.5*OLS)"))

        # Start 3: perturbed OLS (scale by 1.5)
        if self._n_starts >= 3:
            x0_15 = self._beta_init.copy()
            x0_15[1:] *= 1.5
            starts.append((x0_15, "Start 3 (1.5*OLS)"))

        for x0, label in starts[:self._n_starts]:
            _output(f"  {label}...")
            p, crit, diag = self._run_single_optimization(x0, label)
            diagnostics.append(diag)
            _output(f"    criterion = {crit:.6e}")
            if crit < best_crit:
                best_crit = crit
                best_betas = p.copy()

        # report convergence
        crits = [d.criterion for d in diagnostics]
        crit_range = max(crits) - min(crits)
        if crit_range > 0.01:
            _output(f"  WARNING: criteria range {crit_range:.6f} "
                    f"-- local optima detected")
        else:
            _output(f"  OK: all starts within {crit_range:.6f} of best")

        _output(f"  Best criterion: {best_crit:.6e}")
        for name, b in zip(self._beta_names, best_betas):
            _output(f"    {name:8s} = {b:10.{options.digits}f}")

        return best_betas, best_crit, diagnostics

    # ---------------------------------------------------------------------- #
    #  Analytical standard errors (ACH 2012)
    # ---------------------------------------------------------------------- #

    def _moment_vector(self, betas: Array) -> Array:
        """(1/N) Z'xi -- the moment condition vector."""
        xi = self._compute_xi(betas)
        return self._Z.T @ xi / self._N

    def _numerical_jacobian(self, betas: Array) -> Array:
        """Central-difference Jacobian of moment vector w.r.t. betas.

        Uses epsilon from module options (default: sqrt of machine epsilon).
        """
        h = options.finite_differences_epsilon
        K = len(betas)
        m0 = self._moment_vector(betas)
        G = np.zeros((len(m0), K), dtype=options.dtype)
        for j in range(K):
            bp = betas.copy()
            bm = betas.copy()
            bp[j] += h
            bm[j] -= h
            G[:, j] = (self._moment_vector(bp) - self._moment_vector(bm)) / (2 * h)
        return G

    def _analytical_vcov(self, betas: Array) -> Tuple[Array, Array]:
        """GMM sandwich VCV with firm-level clustering.

        V = (1/N) inv(G'WG) G'W S W G inv(G'WG)

        where G is the Jacobian, W is the weighting matrix, and S is the
        clustered meat matrix with small-sample correction.

        Returns
        -------
        vcov : Array
            Variance-covariance matrix (K x K).
        se : Array
            Standard errors (K,).
        """
        G = self._numerical_jacobian(betas)
        xi = self._compute_xi(betas)
        N = self._N

        # clustered meat
        Z_xi = self._Z * xi.reshape(-1, 1)
        K_z = self._Z.shape[1]
        S = np.zeros((K_z, K_z), dtype=options.dtype)
        for c in self._unique_clusters:
            sel = self._cluster_ids == c
            mc = Z_xi[sel].sum(axis=0)
            S += np.outer(mc, mc)
        N_c = self._N_clusters
        S = S / N * (N_c / (N_c - 1))

        # sandwich
        GWG = G.T @ self._W @ G
        cond = np.linalg.cond(GWG)
        if cond > options.singular_tol:
            warnings.warn(f"G'WG condition number ({cond:.2e}) exceeds "
                          f"tolerance -- SEs may be unreliable")
        if options.pseudo_inverses:
            GWG_inv = np.linalg.pinv(GWG)
        else:
            GWG_inv = np.linalg.inv(GWG)

        V = GWG_inv @ G.T @ self._W @ S @ self._W @ G @ GWG_inv / N
        se = np.sqrt(np.maximum(np.diag(V), 0))

        _output(f"\n  Analytical SEs (ACH 2012, {N_c} clusters):")
        for name, b, s in zip(self._beta_names, betas, se):
            _output(f"    {name:8s}  coef={b:10.{options.digits}f}  "
                    f"se={s:10.{options.digits}f}")

        return V, se

    def _hansen_j_test(self, betas: Array) -> Tuple[float, float]:
        """Hansen J test for overidentifying restrictions.

        J = N * m(β)' W_opt m(β) ~ χ²(n_overid)

        Uses the efficient (optimal) weighting matrix W = inv(S) where S
        is the clustered covariance of moment conditions.

        Returns
        -------
        j_stat : float
            Hansen J statistic.
        p_value : float
            p-value from chi-squared distribution.
        """
        from scipy.stats import chi2

        xi = self._compute_xi(betas)
        N = self._N
        n_moments = self._Z.shape[1]
        n_params = len(betas)
        n_overid = n_moments - n_params
        if n_overid <= 0:
            return 0.0, 1.0

        # efficient W = inv(S) where S is clustered moment covariance
        Z_xi = self._Z * xi.reshape(-1, 1)
        K_z = n_moments
        S = np.zeros((K_z, K_z), dtype=options.dtype)
        for c in self._unique_clusters:
            sel = self._cluster_ids == c
            mc = Z_xi[sel].sum(axis=0)
            S += np.outer(mc, mc)
        S = S / N

        if options.pseudo_inverses:
            S_inv = np.linalg.pinv(S)
        else:
            S_inv = np.linalg.inv(S)

        m = self._Z.T @ xi / N
        j_stat = float(N * m.T @ S_inv @ m)
        p_value = float(1 - chi2.cdf(j_stat, n_overid))
        return j_stat, p_value

    # ---------------------------------------------------------------------- #
    #  Chamberlain (1987) optimal instruments
    # ---------------------------------------------------------------------- #

    def _observation_jacobian(self, betas: Array) -> Array:
        r"""Analytical observation-level Jacobian :math:`\partial\xi_i/\partial\beta`.

        For ACF, :math:`\xi(\beta) = \omega - g(\omega_{lag})` where
        :math:`\omega = \phi - X\beta` and
        :math:`\omega_{lag} = \phi_{lag} - X_{lag}\beta`.

        The derivative is:

        .. math::

            \frac{\partial\xi_i}{\partial\beta} = -X_i + g'(\omega_{lag,i}) X_{lag,i}

        where :math:`g'` is the derivative of the Markov polynomial.

        Parameters
        ----------
        betas : Array
            Production function parameters.

        Returns
        -------
        J : Array
            N x K matrix of observation-level Jacobians.
        """
        betas = np.asarray(betas, dtype=options.dtype).ravel()
        omega = self._PHI - self._X @ betas
        omega_lag = self._PHI_LAG - self._X_lag @ betas

        # estimate Markov process g(omega_lag) to get coefficients
        omega_lag_col = omega_lag.reshape(-1, 1)
        omega_lag_pol = [self._C, omega_lag_col]
        for p in range(2, self._formulation.ar_order + 1):
            omega_lag_pol.append((omega_lag ** p).reshape(-1, 1))
        for ctrl in self._markov_extras:
            omega_lag_pol.append(ctrl)
            if self._extensions.markov_interactions:
                omega_lag_pol.append(ctrl * omega_lag_col)
        omega_lag_pol = np.hstack(omega_lag_pol)
        gamma = np.linalg.lstsq(omega_lag_pol, omega, rcond=None)[0]

        # g'(omega_lag) = gamma[1] + 2*gamma[2]*omega_lag + ...
        g_prime = np.full(self._N, gamma[1], dtype=options.dtype)
        for p in range(2, self._formulation.ar_order + 1):
            g_prime += p * gamma[p] * omega_lag ** (p - 1)

        # J_i = -X_i + g'(omega_lag_i) * X_lag_i
        J = -self._X + g_prime.reshape(-1, 1) * self._X_lag
        return J

    def _construct_optimal_instruments(self, betas: Array) -> Array:
        r"""Construct Chamberlain (1987) / Amemiya (1977) optimal instruments.

        The optimal instrument for each observation is:

        .. math::

            Z^*_i = E\left[\frac{\partial\xi_i}{\partial\beta}
            \mid Z_i\right]

        (up to a scaling by the inverse error variance, which is absorbed
        by the weighting matrix). The conditional expectation is approximated
        via sieve regression: project each column of the Jacobian J = ∂ξ/∂β
        onto a polynomial basis in Z.

        This produces exactly K_beta columns — the same dimensionality as
        the parameter vector. When used in 'replace' mode, the system is
        just-identified but with the asymptotically efficient combination
        of instrument information (Chamberlain 1987; Conlon & Gortmaker 2020).

        Parameters
        ----------
        betas : Array
            Initial consistent parameter estimates.

        Returns
        -------
        Z_opt : Array
            N x K_beta optimal instrument matrix.
        """
        J = self._observation_jacobian(betas)
        K = J.shape[1]

        # build expanded sieve basis from ALL current instruments
        # (includes standard ACF instruments + any oligopoly instruments)
        sieve_order = self._formulation.optimal_instruments_sieve_order
        Z_base = self._Z  # N x n_z

        # normalize Z for numerical stability in polynomial expansion
        Z_mean = Z_base.mean(axis=0)
        Z_std = Z_base.std(axis=0)
        Z_std[Z_std < 1e-10] = 1.0
        Z_norm = (Z_base - Z_mean) / Z_std

        # build polynomial sieve: columns of Z, their squares, cross-products
        sieve_cols = [Z_norm]
        if sieve_order >= 2:
            n_z = Z_norm.shape[1]
            for j in range(n_z):
                if Z_std[j] > 1e-10:
                    sieve_cols.append((Z_norm[:, j] ** 2).reshape(-1, 1))
            for j in range(n_z):
                for l in range(j + 1, n_z):
                    if Z_std[j] > 1e-10 and Z_std[l] > 1e-10:
                        sieve_cols.append(
                            (Z_norm[:, j] * Z_norm[:, l]).reshape(-1, 1))
        if sieve_order >= 3:
            for j in range(n_z):
                if Z_std[j] > 1e-10:
                    sieve_cols.append((Z_norm[:, j] ** 3).reshape(-1, 1))

        Z_sieve = np.hstack(sieve_cols)
        n_sieve = Z_sieve.shape[1]
        _output(f"  Optimal instruments: sieve basis with {n_sieve} terms "
                f"from {Z_base.shape[1]} base instruments (order {sieve_order})")

        # project each column of J onto sieve basis: E[J_k | Z] ≈ Z_sieve @ δ_k
        # ridge regression for stability (small penalty avoids overfitting)
        ridge_lambda = 1e-6
        ZtZ_sieve = Z_sieve.T @ Z_sieve
        reg = ZtZ_sieve + ridge_lambda * np.eye(n_sieve, dtype=options.dtype)

        Z_opt = np.zeros((self._N, K), dtype=options.dtype)
        for k in range(K):
            delta = np.linalg.solve(reg, Z_sieve.T @ J[:, k])
            Z_opt[:, k] = Z_sieve @ delta

        # check rank of optimal instrument matrix
        rank = np.linalg.matrix_rank(Z_opt, tol=1e-8)
        _output(f"  Optimal instruments: rank = {rank} / {K} "
                f"(mode='{self._formulation.optimal_instruments}')")
        if rank < K:
            warnings.warn(f"Optimal instrument matrix is rank-deficient "
                          f"({rank} < {K}). Results may be unreliable.")

        return Z_opt

    # ---------------------------------------------------------------------- #
    #  Markups and delta-method SEs
    # ---------------------------------------------------------------------- #

    def _compute_markups(self, betas: Array, vcov: Array) -> pd.DataFrame:
        """Compute mu = theta^V / alpha^V and delta-method SEs.

        For CD: theta = beta_cogs (constant).
        For TL: theta = beta_cogs + 2*beta_cogs2*cogs + beta_kcogs*k.
        """
        df = self._est_data.copy()
        alpha = df['alphahat'].values.astype(options.dtype)

        if self._formulation.spec == 'cd':
            theta = np.full(len(df), betas[2], dtype=options.dtype)
            grad = np.zeros((len(df), len(betas)), dtype=options.dtype)
            grad[:, 2] = 1.0
        else:
            cogs_v = df['cogs'].values.astype(options.dtype)
            k_v = df['k'].values.astype(options.dtype)
            theta = betas[2] + 2 * betas[4] * cogs_v + betas[5] * k_v
            grad = np.zeros((len(df), len(betas)), dtype=options.dtype)
            grad[:, 2] = 1.0
            grad[:, 4] = 2 * cogs_v
            grad[:, 5] = k_v

        markup = theta / alpha

        # delta method: se_theta_i = sqrt(grad_i' V grad_i)
        # vectorized: se_theta = sqrt(diag(grad @ V @ grad.T))
        # but N can be large, so compute row-by-row via einsum
        gVg = np.einsum('ij,jk,ik->i', grad, vcov, grad)
        se_theta = np.sqrt(np.maximum(gVg, 0))
        se_markup = se_theta / np.abs(alpha)

        df['theta'] = theta
        df['markup'] = markup
        df['se_theta'] = se_theta
        df['se_markup'] = se_markup
        df['omega'] = self._PHI - self._X @ betas

        return df


# ========================================================================== #
#  Industry-by-industry estimation
# ========================================================================== #

def estimate_by_industry(
    data: pd.DataFrame,
    specs: Sequence[str] = ('cd', 'tl'),
    formulation_kwargs: Optional[Options] = None,
    optimization_kwargs: Optional[Options] = None,
    extensions_kwargs: Optional[Options] = None,
    ic_kwargs: Optional[Options] = None,
    n_starts: int = 3,
) -> Tuple[Dict[Tuple, ACFResults], pd.DataFrame, pd.DataFrame]:
    """Estimate ACF separately by nace2, combine results.

    Parameters
    ----------
    data : pd.DataFrame
        Full panel with id, year, go, k, cogs, pp_dummy, nace2.
    specs : sequence of str
        Which specifications to estimate (default: both CD and TL).
    formulation_kwargs : dict, optional
        Keyword arguments for :class:`Formulation`.
    optimization_kwargs : dict, optional
        Keyword arguments for :class:`Optimization`.
    extensions_kwargs : dict, optional
        Keyword arguments for :class:`CWDLExtensions`.
    ic_kwargs : dict, optional
        Keyword arguments for :class:`ImperfectCompetition`.
    n_starts : int
        Number of optimization starting points (default 3).

    Returns
    -------
    results : dict
        ``{(nace2, spec): ACFResults}`` for each industry x specification.
    coefficients : pd.DataFrame
        Summary table of PF parameters and SEs.
    markups : pd.DataFrame
        Combined firm-level markups from all industries.
    """
    results = {}
    coeff_rows = []
    markup_frames = []

    formulation_kwargs = formulation_kwargs or {}
    optimization_kwargs = optimization_kwargs or {}
    extensions_kwargs = extensions_kwargs or {}
    ic_kwargs = ic_kwargs or {}

    industries = sorted(data['nace2'].unique())

    for nace in industries:
        df_nace = data[data['nace2'] == nace].copy()
        _output(f"\n{'=' * 60}")
        _output(f"  NACE2 = {nace} (N_firms={df_nace['id'].nunique()}, "
                f"N_obs={len(df_nace)})")
        _output(f"{'=' * 60}")

        for spec in specs:
            _output(f"\n--- {spec.upper()} ---")

            form = Formulation(spec=spec, nace2_fe=False, **formulation_kwargs)
            opt = Optimization(**optimization_kwargs)
            ext = CWDLExtensions(**extensions_kwargs)
            ic = ImperfectCompetition(**ic_kwargs) if ic_kwargs else None

            est = ACFEstimator(df_nace, formulation=form, optimization=opt,
                               extensions=ext,
                               imperfect_competition=ic,
                               n_starts=n_starts)
            res = est.solve()
            results[(nace, spec)] = res

            # coefficient row
            row: Dict[str, Any] = {
                'nace2': nace, 'spec': spec, 'N_obs': res.n_obs,
            }
            for name, b, s in zip(res.beta_names, res.betas, res.se):
                row[name] = b
                row[f'se_{name}'] = s
            row['gmm_criterion'] = res.gmm_criterion
            coeff_rows.append(row)

            # markups
            mu_df = res.data[
                ['id', 'year', 'nace2', 'markup', 'se_markup',
                 'theta', 'omega', 'alphahat', 'pp_dummy', 'k', 'cogs']
            ].copy()
            mu_df['spec'] = spec
            markup_frames.append(mu_df)

    coefficients = pd.DataFrame(coeff_rows)
    markups = pd.concat(markup_frames, ignore_index=True)

    return results, coefficients, markups


# ========================================================================== #
#  Hall (2018) distributional decomposition
# ========================================================================== #

def hall_decomposition(
    markups: Array, label: str = '',
) -> Optional[Dict[str, float]]:
    """Hall (2018) decomposition of markup variance.

    Decomposes observed markup variance into true heterogeneity and
    estimation noise. Assumes mu-1 = V + e where V ~ LogNormal(delta, sigma^2)
    and e ~ N(0, gamma^2). Matches first 3 moments of (mu-1).

    Parameters
    ----------
    markups : array-like
        Firm-level markup estimates.
    label : str
        Label for printing.

    Returns
    -------
    dict or None
        Keys: gamma, delta, sigma, mean_corrected, sd_true, sd_total,
        noise_share, n_used. Returns None if fewer than 10 observations.
    """
    mu = np.asarray(markups, dtype=options.dtype)
    mu_minus1 = mu - 1.0
    sel = mu_minus1 > 0
    if sel.sum() < 10:
        if label:
            _output(f"  {label}: fewer than 10 obs with mu > 1, skipping")
        return None

    x = mu_minus1[sel]
    M1, M2, M3 = x.mean(), (x ** 2).mean(), (x ** 3).mean()

    def objective(params: Array) -> float:
        gamma, delta, sigma = params
        pred_M1 = np.exp(delta + 0.5 * sigma ** 2)
        pred_M2 = gamma ** 2 + np.exp(2 * delta + 2 * sigma ** 2)
        pred_M3 = (np.exp(3 * delta + 4.5 * sigma ** 2)
                    + 3 * gamma ** 2 * np.exp(delta + 0.5 * sigma ** 2))
        return (pred_M1 - M1) ** 2 + (pred_M2 - M2) ** 2 + (pred_M3 - M3) ** 2

    res = minimize(objective, x0=[0.1, 0.1, 0.5], method='Nelder-Mead',
                   options={'xatol': 1e-16, 'fatol': 1e-16, 'maxiter': 10000})

    gamma, delta, sigma = res.x
    var_V = (np.exp(sigma ** 2) - 1) * np.exp(2 * delta + sigma ** 2)
    sd_V = np.sqrt(max(var_V, 0))
    mean_V = np.exp(delta + sigma ** 2 / 2)
    denom = var_V + gamma ** 2
    noise_share = gamma ** 2 / denom if denom > 0 else 0.0

    result = {
        'gamma': float(gamma), 'delta': float(delta), 'sigma': float(sigma),
        'mean_corrected': float(1 + mean_V), 'sd_true': float(sd_V),
        'sd_total': float(np.sqrt(max(denom, 0))),
        'noise_share': float(noise_share), 'n_used': int(sel.sum()),
    }

    if label:
        _output(f"  {label}: N={sel.sum()}, corrected mean="
                f"{result['mean_corrected']:.4f}, true SD={sd_V:.4f}, "
                f"noise share={100*noise_share:.1f}%")

    return result


# ========================================================================== #
#  Main entry point
# ========================================================================== #

if __name__ == '__main__':
    import sys

    # Prefer rebuilt data (has mktshare, additional variables)
    from pathlib import Path as _Path
    _script_dir = _Path(__file__).resolve().parent
    _input_dir = _script_dir.parent / 'input'
    rebuilt_path = str(_input_dir / 'data_rebuilt.dta')
    orig_path = str(_input_dir / 'data.dta')
    import os
    data_path = rebuilt_path if os.path.exists(rebuilt_path) else orig_path
    print(f'Loading data from {data_path}')
    df = pd.read_stata(data_path)

    required = ['id', 'year', 'go', 'k', 'cogs', 'pp_dummy', 'nace2']
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f'Missing columns: {missing}')
        print(f'Available: {sorted(df.columns.tolist())}')
        sys.exit(1)

    # Construct market share if not present (DGM first-stage control)
    if 'mktshare' not in df.columns:
        df['rGO'] = np.exp(df['go'])
        total = df.groupby(['nace2', 'year'])['rGO'].transform('sum')
        df['mktshare'] = df['rGO'] / total
        print(f'  Constructed mktshare: {df["mktshare"].notna().sum()} obs')

    print(f'Data: {len(df)} obs, {df["id"].nunique()} firms, '
          f'years {df["year"].min()}-{df["year"].max()}')
    print(f'Industries: {sorted(df["nace2"].unique())}')

    # estimate by industry (DGM-consistent: include market share in first stage)
    results, coefficients, markups = estimate_by_industry(
        df,
        formulation_kwargs={'first_stage_controls': ['mktshare']},
    )

    # summary
    print('\n' + '=' * 70)
    print('  COEFFICIENT SUMMARY')
    print('=' * 70)
    print(coefficients.to_string(index=False, float_format='{:.6f}'.format))

    # procurement premium: OLS with controls (matching Stata)
    print('\n' + '=' * 70)
    print('  PROCUREMENT PREMIUM')
    print('=' * 70)
    for spec in ['cd', 'tl']:
        mu_spec = markups[markups['spec'] == spec].copy()
        mu_spec = mu_spec[mu_spec['markup'] > 0].copy()
        mu_spec['lmu'] = np.log(mu_spec['markup'])

        # raw difference
        pp0 = mu_spec[mu_spec['pp_dummy'] == 0]['lmu']
        pp1 = mu_spec[mu_spec['pp_dummy'] == 1]['lmu']
        raw_diff = pp1.mean() - pp0.mean()

        # OLS: lmu ~ pp_dummy + k + cogs + year*nace2 FE, clustered by id
        rhs = mu_spec[['pp_dummy', 'k', 'cogs']].copy()
        rhs['_const'] = 1.0
        # year x nace2 FE
        yr_nace = mu_spec['year'].astype(str) + '_' + mu_spec['nace2'].astype(str)
        for val in sorted(yr_nace.unique())[1:]:
            rhs[f'_fe_{val}'] = (yr_nace == val).astype(float)
        X = rhs.values.astype(np.float64)
        y = mu_spec['lmu'].values.astype(np.float64)
        mask = np.isfinite(X).all(axis=1) & np.isfinite(y)
        X, y = X[mask], y[mask]
        bhat = np.linalg.lstsq(X, y, rcond=None)[0]
        resid = y - X @ bhat
        pp_coef = bhat[0]

        # clustered SE on pp_dummy coefficient
        ids = mu_spec.loc[mu_spec.index[mask], 'id'].values
        XtXinv = np.linalg.inv(X.T @ X)
        meat = np.zeros((X.shape[1], X.shape[1]))
        for cid in np.unique(ids):
            idx = ids == cid
            score = X[idx].T @ resid[idx].reshape(-1, 1)
            meat += score @ score.T
        N, K = X.shape
        G = len(np.unique(ids))
        dfc = G / (G - 1) * (N - 1) / (N - K)
        vcov = dfc * XtXinv @ meat @ XtXinv
        pp_se = np.sqrt(vcov[0, 0])

        print(f'  {spec.upper()}: raw diff = {raw_diff:.4f}, '
              f'OLS(controls) = {pp_coef:.4f} (SE = {pp_se:.4f})')

    # Hall decomposition
    print('\n' + '=' * 70)
    print('  HALL (2018) DECOMPOSITION')
    print('=' * 70)
    for spec in ['cd', 'tl']:
        mu_spec = markups[markups['spec'] == spec]
        firm_means = mu_spec.groupby('id')['markup'].mean()
        hall_decomposition(firm_means.values, label=f'{spec.upper()} all')

    # save
    _output_dir = _script_dir.parent / 'output'
    out_path = str(_output_dir / 'acf_python_markups.csv')
    markups.to_csv(out_path, index=False)
    print(f'\nMarkups saved to {out_path}')

    coeff_path = str(_output_dir / 'acf_python_coefficients.csv')
    coefficients.to_csv(coeff_path, index=False)
    print(f'Coefficients saved to {coeff_path}')

    # ================================================================== #
    #  OVERIDENTIFIED TRANSLOG (Kim, Luo & Su 2019)
    # ================================================================== #
    print('\n' + '=' * 70)
    print('  OVERIDENTIFIED TRANSLOG (deeper lags: Lk, L2cogs)')
    print('=' * 70)

    form_overid = Formulation(spec='tl', overidentify=True)
    for nace in sorted(df['nace2'].unique()):
        df_n = df[df['nace2'] == nace].copy()
        print(f'\n--- NACE2 = {nace} ---')
        try:
            est = ACFEstimator(
                data=df_n,
                formulation=form_overid,
                optimization=Optimization(method='nm+bfgs'),
            )
            res = est.solve()
            print(res)
        except Exception as e:
            print(f'  Failed: {e}')

    # ================================================================== #
    #  FSOLVE on just-identified TL (DGM approach)
    # ================================================================== #
    print('\n' + '=' * 70)
    print('  FSOLVE: just-identified TL (root-finding, DGM)')
    print('=' * 70)

    form_tl = Formulation(spec='tl', overidentify=False)
    for nace in sorted(df['nace2'].unique()):
        df_n = df[df['nace2'] == nace].copy()
        print(f'\n--- NACE2 = {nace} ---')
        try:
            est = ACFEstimator(
                data=df_n,
                formulation=form_tl,
                optimization=Optimization(method='fsolve'),
            )
            res = est.solve()
            print(res)
        except Exception as e:
            print(f'  Failed: {e}')

    # ================================================================== #
    #  CHAMBERLAIN OPTIMAL INSTRUMENTS
    # ================================================================== #
    print('\n' + '=' * 70)
    print('  CHAMBERLAIN (1987) OPTIMAL INSTRUMENTS')
    print('=' * 70)

    for spec in ['cd', 'tl']:
        for oi_mode in ['replace', 'augment']:
            form_opt = Formulation(spec=spec, optimal_instruments=oi_mode,
                                   optimal_instruments_sieve_order=2)
            for nace in sorted(df['nace2'].unique()):
                df_n = df[df['nace2'] == nace].copy()
                print(f'\n--- {spec.upper()} NACE2 = {nace} '
                      f'(Chamberlain {oi_mode}) ---')
                try:
                    est = ACFEstimator(
                        data=df_n,
                        formulation=form_opt,
                        optimization=Optimization(method='nm+bfgs'),
                    )
                    res = est.solve()
                    print(res)
                except Exception as e:
                    print(f'  Failed: {e}')

    # ================================================================== #
    #  BLP-STYLE INSTRUMENTS & ADL IMPERFECT COMPETITION
    # ================================================================== #
    print('\n' + '=' * 70)
    print('  BLP-STYLE INSTRUMENTS & ADL (2024) IMPERFECT COMPETITION')
    print('=' * 70)

    # ------------------------------------------------------------------
    #  Construct BLP instruments from year x nace2 markets
    # ------------------------------------------------------------------
    # Market = year x nace2 (all firms, not just procurement).
    # For each firm j in market t:
    #   Z_jt^rival(x) = sum of x_kt for all k != j in market t
    # Following pyblp.build_blp_instruments logic (no own-firm/rival
    # distinction since Czech construction firms are single-establishment).
    #
    # Characteristics used:
    #   k (capital), go (output), cogs (intermediate inputs)
    # Instruments constructed:
    #   comp_k_sum, comp_go_sum, comp_cogs_sum  — BLP sums of rivals
    #   comp_k_mean, comp_cogs_mean             — BLP means (scaled)
    #   comp_n, comp_n_log                      — number of competitors
    #   comp_k_sd                               — dispersion of rival k
    # Sufficient statistic (QEC, ADL 2024b):
    #   comp_go_lse — leave-one-out log-sum-exp of rivals' output
    # ------------------------------------------------------------------

    print('\nConstructing BLP instruments (market = year x nace2)...')

    df_ic = df.copy()
    df_ic['mkt_id'] = (df_ic['year'].astype(str) + '_'
                       + df_ic['nace2'].astype(int).astype(str))

    # market aggregates
    chars = ['k', 'go', 'cogs']
    mkt = df_ic.groupby('mkt_id').agg(
        mkt_n=('id', 'count'),
        **{f'mkt_{c}_sum': (c, 'sum') for c in chars},
        **{f'mkt_{c}_mean': (c, 'mean') for c in chars},
        mkt_k_sd=('k', 'std'),
        mkt_go_lse=('go', lambda x: np.log(np.sum(np.exp(x)))),
    ).reset_index()
    df_ic = df_ic.merge(mkt, on='mkt_id', how='left')

    # leave-one-out: subtract own characteristics
    for c in chars:
        df_ic[f'comp_{c}_sum'] = df_ic[f'mkt_{c}_sum'] - df_ic[c]
    df_ic['comp_n'] = df_ic['mkt_n'] - 1
    for c in chars:
        df_ic[f'comp_{c}_mean'] = np.where(
            df_ic['comp_n'] > 0,
            df_ic[f'comp_{c}_sum'] / df_ic['comp_n'],
            np.nan)
    df_ic['comp_n_log'] = np.log(df_ic['comp_n'].clip(lower=1))
    df_ic['comp_k_sd'] = df_ic['mkt_k_sd']  # market-level, no LOO needed

    # sufficient statistic: leave-one-out log-sum-exp of output
    df_ic['comp_go_lse'] = np.where(
        df_ic['comp_n'] > 0,
        np.log(np.clip(
            np.exp(df_ic['mkt_go_lse']) - np.exp(df_ic['go']),
            1e-10, None)),
        np.nan)

    # competitors' lagged productivity (ADL Theorem 6: optimal IV uses
    # theta_f * f_{k1} + rho * omega_{k0}). Approximate omega with Solow
    # residual (go - beta_k*k - beta_cogs*cogs ≈ go - k - cogs as proxy).
    df_ic['omega_proxy'] = df_ic['go'] - 0.05 * df_ic['k'] - 0.95 * df_ic['cogs']
    df_ic['L_omega_proxy'] = df_ic.groupby('id')['omega_proxy'].shift(1)
    mkt_omega = df_ic.groupby('mkt_id').agg(
        mkt_Lomega_sum=('L_omega_proxy', 'sum'),
        mkt_Lomega_n=('L_omega_proxy', 'count'),
    ).reset_index()
    df_ic = df_ic.merge(mkt_omega, on='mkt_id', how='left')
    df_ic['comp_Lomega_mean'] = np.where(
        df_ic['comp_n'] > 0,
        (df_ic['mkt_Lomega_sum'] - df_ic['L_omega_proxy'].fillna(0))
        / df_ic['comp_n'],
        np.nan)

    # report coverage
    has_comp = df_ic['comp_go_lse'].notna() & (df_ic['comp_n'] > 0)
    print(f'  Firms with competitors: {has_comp.sum()} / {len(df_ic)} '
          f'({100 * has_comp.mean():.1f}%)')
    for nace in sorted(df_ic['nace2'].unique()):
        mask_n = (df_ic['nace2'] == nace) & has_comp
        print(f'    NACE {nace}: {mask_n.sum()} obs, '
              f'avg N_rivals = {df_ic.loc[mask_n, "comp_n"].mean():.1f}')

    # ------------------------------------------------------------------
    #  Instrument sets to compare
    # ------------------------------------------------------------------
    # (a) Baseline: comp_k_mean + comp_n_log (2 instruments)
    # (b) BLP sums: comp_k_sum + comp_cogs_sum + comp_n_log (3 instruments)
    # (c) ADL optimal: comp_k_mean + comp_Lomega_mean (Theorem 6: f_{-j} + omega_{-j,t-1})
    # (d) ADL full: comp_k_mean + comp_Lomega_mean + comp_n_log (3 instruments)
    iv_sets = {
        'baseline': ['comp_k_mean', 'comp_n_log'],
        'blp_sums': ['comp_k_sum', 'comp_cogs_sum', 'comp_n_log'],
        'adl_optimal': ['comp_k_mean', 'comp_Lomega_mean'],
        'adl_full': ['comp_k_mean', 'comp_Lomega_mean', 'comp_n_log'],
    }

    # ------------------------------------------------------------------
    #  Estimate: standard vs IC-corrected across instrument sets
    # ------------------------------------------------------------------
    all_results = []
    for spec in ['cd']:
        for nace in sorted(df_ic['nace2'].unique()):
            df_n = df_ic[(df_ic['nace2'] == nace) & has_comp].copy()
            if len(df_n) < 100:
                continue

            print(f'\n{"="*60}')
            print(f'  {spec.upper()} NACE {int(nace)} '
                  f'(N={len(df_n)}, firms={df_n["id"].nunique()})')
            print(f'{"="*60}')

            form = Formulation(spec=spec, pp_in_markov=True,
                              pp_interactions=True, year_fe=True)
            opt = Optimization(method='nm+bfgs')

            # standard (no IC)
            try:
                est_std = ACFEstimator(df_n, formulation=form,
                                       optimization=opt, n_starts=2)
                res_std = est_std.solve()
                print(res_std)
            except Exception as e:
                print(f'  Standard failed: {e}')
                continue

            # IC-corrected with each instrument set
            for iv_name, iv_cols in iv_sets.items():
                print(f'\n  --- IC + {iv_name} ---')
                try:
                    ic_config = ImperfectCompetition(
                        enabled=True,
                        sufficient_statistic='comp_go_lse',
                        oligopoly_instruments=iv_cols,
                    )
                    est_ic = ACFEstimator(df_n, formulation=form,
                                          optimization=opt,
                                          imperfect_competition=ic_config,
                                          n_starts=2)
                    res_ic = est_ic.solve()
                    print(res_ic)
                except Exception as e:
                    print(f'  IC ({iv_name}) failed: {e}')
                    continue

                # compare
                print(f'\n  {iv_name}: '
                      f'k={res_ic.betas[1]:.4f} '
                      f'(std: {res_std.betas[1]:.4f}), '
                      f'cogs={res_ic.betas[2]:.4f} '
                      f'(std: {res_std.betas[2]:.4f})')
                mu_std = res_std.data['markup']
                mu_ic = res_ic.data['markup']
                print(f'  Mean markup: IC={mu_ic.mean():.4f} '
                      f'vs std={mu_std.mean():.4f} '
                      f'(diff={mu_ic.mean()-mu_std.mean():+.4f})')

                # procurement premium
                for label, res in [('std', res_std), ('ic', res_ic)]:
                    d = res.data
                    mu_v = d[d['markup'] > 0].copy()
                    mu_v['lmu'] = np.log(mu_v['markup'])
                    pp0 = mu_v[mu_v['pp_dummy'] == 0]['lmu']
                    pp1 = mu_v[mu_v['pp_dummy'] == 1]['lmu']
                    if len(pp1) > 0 and len(pp0) > 0:
                        print(f'  PP premium ({label}): '
                              f'{pp1.mean()-pp0.mean():.4f}')
                if res_ic.hansen_j is not None:
                    print(f'  Hansen J: chi2({res_ic.n_overid})='
                          f'{res_ic.hansen_j:.3f}, '
                          f'p={res_ic.hansen_j_pvalue:.3f}')

                all_results.append({
                    'spec': spec, 'nace': int(nace), 'iv': iv_name,
                    'k': res_ic.betas[1], 'cogs': res_ic.betas[2],
                    'mu_mean': mu_ic.mean(),
                    'premium': (pp1.mean() - pp0.mean()
                                if len(pp1) > 0 and len(pp0) > 0
                                else np.nan),
                    'hansen_j': res_ic.hansen_j,
                    'hansen_p': res_ic.hansen_j_pvalue,
                    'n_overid': res_ic.n_overid,
                })

    # summary table
    if all_results:
        print('\n' + '=' * 70)
        print('  BLP INSTRUMENTS SUMMARY')
        print('=' * 70)
        print(f'  {"NACE":<6} {"IV set":<12} {"k":>8} {"cogs":>8} '
              f'{"mu":>8} {"prem":>8} {"J":>8} {"p":>6} {"n_oi":>5}')
        print(f'  {"-"*6} {"-"*12} {"-"*8} {"-"*8} '
              f'{"-"*8} {"-"*8} {"-"*8} {"-"*6} {"-"*5}')
        for r in all_results:
            j_str = f'{r["hansen_j"]:.2f}' if r['hansen_j'] is not None else ''
            p_str = f'{r["hansen_p"]:.3f}' if r['hansen_p'] is not None else ''
            n_str = str(r['n_overid']) if r['n_overid'] else ''
            print(f'  {r["nace"]:<6} {r["iv"]:<12} {r["k"]:>8.4f} '
                  f'{r["cogs"]:>8.4f} {r["mu_mean"]:>8.4f} '
                  f'{r["premium"]:>8.4f} {j_str:>8} {p_str:>6} '
                  f'{n_str:>5}')

    # ================================================================== #
    #  VARIABLE INPUT COMPARISON: COGS vs II vs W+II
    # ================================================================== #
    print('\n' + '=' * 70)
    print('  VARIABLE INPUT COMPARISON: COGS vs II')
    print('=' * 70)

    # Check ii is available
    if 'ii' not in df.columns:
        print('  SKIPPED: ii (log intermediate inputs) not in data')
    else:
        vi_results = []
        for vi_name in ['cogs', 'ii']:
            for nace in sorted(df['nace2'].unique()):
                df_n = df[df['nace2'] == nace].copy()

                form = Formulation(spec='cd', variable_input=vi_name)
                opt = Optimization(method='nm+bfgs')
                try:
                    est = ACFEstimator(df_n, formulation=form,
                                       optimization=opt, n_starts=3)
                    res = est.solve()
                except Exception as e:
                    print(f'  {vi_name} NACE {int(nace)}: failed ({e})')
                    continue

                # procurement premium
                d_r = res.data
                mu_v = d_r[d_r['markup'] > 0].copy()
                mu_v['lmu'] = np.log(mu_v['markup'])
                pp0 = mu_v[mu_v['pp_dummy'] == 0]['lmu']
                pp1 = mu_v[mu_v['pp_dummy'] == 1]['lmu']
                prem = pp1.mean() - pp0.mean() if len(pp1) > 0 else np.nan

                vi_results.append({
                    'vi': vi_name, 'nace': int(nace),
                    'k': res.betas[1], 'vi_coef': res.betas[2],
                    'mu_mean': d_r['markup'].mean(),
                    'mu_med': d_r['markup'].median(),
                    'mu_sd': d_r['markup'].std(),
                    'premium': prem,
                    'N': res.n_obs,
                    'share_med': d_r['alphahat'].median(),
                })
                print(f'  {vi_name:>5} NACE {int(nace)}: '
                      f'k={res.betas[1]:.4f}, '
                      f'vi={res.betas[2]:.4f}, '
                      f'mu={d_r["markup"].mean():.3f} '
                      f'(med={d_r["markup"].median():.3f}), '
                      f'share={d_r["alphahat"].median():.3f}, '
                      f'prem={prem:.4f}')

        if vi_results:
            print(f'\n  {"VI":<6} {"NACE":<6} {"k":>7} {"β_vi":>7} '
                  f'{"α̂_med":>7} {"μ_mean":>7} {"μ_med":>7} '
                  f'{"prem":>7} {"N":>6}')
            print(f'  {"-"*6} {"-"*6} {"-"*7} {"-"*7} '
                  f'{"-"*7} {"-"*7} {"-"*7} {"-"*7} {"-"*6}')
            for r in vi_results:
                print(f'  {r["vi"]:<6} {r["nace"]:<6} '
                      f'{r["k"]:>7.4f} {r["vi_coef"]:>7.4f} '
                      f'{r["share_med"]:>7.3f} {r["mu_mean"]:>7.3f} '
                      f'{r["mu_med"]:>7.3f} {r["premium"]:>7.4f} '
                      f'{r["N"]:>6}')

        # --- CWDL base spec with ii ---
        print('\n  --- CWDL base (survival + pp_dummy) with ii ---')
        for nace in sorted(df['nace2'].unique()):
            df_n = df[df['nace2'] == nace].copy()
            form = Formulation(spec='cd', variable_input='ii')
            opt = Optimization(method='nm+bfgs')
            ext = CWDLExtensions(survival_correction=True)
            try:
                est = ACFEstimator(df_n, formulation=form,
                                   optimization=opt, extensions=ext,
                                   n_starts=3)
                res = est.solve()
            except Exception as e:
                print(f'  ii CWDL NACE {int(nace)}: failed ({e})')
                continue

            d_r = res.data
            mu_v = d_r[d_r['markup'] > 0].copy()
            mu_v['lmu'] = np.log(mu_v['markup'])
            pp0 = mu_v[mu_v['pp_dummy'] == 0]['lmu']
            pp1 = mu_v[mu_v['pp_dummy'] == 1]['lmu']
            prem = pp1.mean() - pp0.mean() if len(pp1) > 0 else np.nan
            print(f'  ii CWDL NACE {int(nace)}: '
                  f'k={res.betas[1]:.4f}, ii={res.betas[2]:.4f}, '
                  f'mu={d_r["markup"].mean():.3f} '
                  f'(med={d_r["markup"].median():.3f}), '
                  f'share={d_r["alphahat"].median():.3f}, '
                  f'prem={prem:.4f}')
