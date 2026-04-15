"""
Misspecification Diagnostics (Andrews, Chen, and Tecchio 2026).

Part 1: Hansen J-statistic as a weight-sensitivity bound.
  Under local misspecification, sqrt(J) bounds the number of standard errors
  by which the estimate can be moved by varying the GMM weight matrix
  (Proposition 5.2).

Part 2: Misspecification-robust standard errors for the procurement premium
  via non-recentered cluster bootstrap (pp.18-19). The non-recentered
  bootstrap is valid under misspecification, while the analytical SE
  assumes correct specification.

References
----------
Andrews, Chen & Tecchio (2026): arXiv:2508.13076v5.

Author: Marek Chadim (Yale, Tobin Center)
"""

from __future__ import annotations

import sys
import time
from typing import Optional
import numpy as np
import pandas as pd
import statsmodels.api as sm
from pathlib import Path
from scipy.stats import chi2

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR / 'lib'))

from acf_estimator import (
    ACFEstimator, Formulation, Optimization, CWDLExtensions,
    options as acf_options,
)

INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
acf_options.verbose = False

B_BOOT = 999
RNG = np.random.default_rng(42)


# ================================================================== #
#  Part 1: J-statistic for ACF specifications
# ================================================================== #

def compute_j_diagnostic(est: ACFEstimator, beta_hat: np.ndarray,
                         spec_label: str) -> dict:
    """Compute Hansen J and sqrt(J) weight-hacking bound.

    Parameters
    ----------
    est : ACFEstimator
        Fitted estimator (after solve()).
    beta_hat : array
        Estimated betas.
    spec_label : str
        Label for display and output.

    Returns
    -------
    dict with J, sqrt_J, p_value, df, n_moments, n_params.
    """
    g = est._moment_vector(beta_hat)
    M = len(g)
    K = len(beta_hat)
    df = M - K
    N = est._N

    # Clustered variance of moment vector
    xi = est._compute_xi(beta_hat)
    Z_xi = est._Z * xi.reshape(-1, 1)
    Sigma = np.zeros((M, M))
    for c in est._unique_clusters:
        sel = est._cluster_ids == c
        mc = Z_xi[sel].sum(axis=0)
        Sigma += np.outer(mc, mc)
    Sigma /= N

    # J = N * g' Sigma^{-1} g
    try:
        J = float(N * g @ np.linalg.solve(Sigma, g))
    except np.linalg.LinAlgError:
        J = float(N * g @ np.linalg.lstsq(Sigma, g, rcond=None)[0])

    sqrt_J = np.sqrt(max(J, 0.0))

    if df > 0:
        p_value = float(1 - chi2.cdf(J, df))
    else:
        p_value = np.nan

    print(f'\n  {spec_label}:')
    print(f'    Moments={M}, Params={K}, df={df}')
    print(f'    J = {J:.4f}, sqrt(J) = {sqrt_J:.4f}, '
          f'p = {p_value:.4f}' if np.isfinite(p_value)
          else f'    J = {J:.4f}, sqrt(J) = {sqrt_J:.4f}, p = NA (just-identified)')

    return {
        'Specification': spec_label,
        'J': J, 'sqrt_J': sqrt_J, 'p_value': p_value,
        'df': df, 'n_moments': M, 'n_params': K,
    }


# ================================================================== #
#  Part 2: Non-recentered cluster bootstrap for premium SE
# ================================================================== #

def ols_premium(df: pd.DataFrame) -> tuple[float, float]:
    """Run OLS of log(markup_A) on pp_dummy + year/nace2 FEs, clustered.

    Returns (coefficient, clustered SE).
    """
    sub = df[['id', 'year', 'nace2', 'ln_mu', 'D']].dropna().copy()
    if len(sub) < 30:
        return np.nan, np.nan

    yr_dum = pd.get_dummies(sub['year'], prefix='yr', drop_first=True, dtype=float)
    n_dum = pd.get_dummies(sub['nace2'], prefix='n', drop_first=True, dtype=float)
    X = pd.concat([sub[['D']], yr_dum, n_dum], axis=1)
    X = sm.add_constant(X)
    y = sub['ln_mu']

    res = sm.OLS(y, X).fit(cov_type='cluster',
                           cov_kwds={'groups': sub['id']})
    return float(res.params['D']), float(res.bse['D'])


def bootstrap_premium_se(df: pd.DataFrame, B: int = B_BOOT,
                          theta_hat: Optional[float] = None,
                          pivotal_ci: bool = True) -> dict:
    """Non-recentered cluster bootstrap for premium SE.

    Resamples firms (clusters) with replacement, reconstructs panel,
    re-runs OLS. The non-recentered bootstrap captures both sampling
    uncertainty AND misspecification bias (Andrews, Chen, Tecchio p.18-19).

    When ``theta_hat`` is provided AND ``pivotal_ci=True`` (the default),
    returns Conlon's "Better Way" pivotal CI ``[2*theta_hat - q97.5,
    2*theta_hat - q2.5]`` (bootstrap.tex lines 75-82) plus a
    bias-corrected point estimate. Otherwise returns the legacy
    percentile CI ``[q2.5, q97.5]`` (backward-compatible).
    """
    firms = df['id'].unique()
    N_firms = len(firms)
    beta_boot = np.full(B, np.nan)

    t0 = time.time()
    for b in range(B):
        # Resample firms with replacement
        boot_firms = RNG.choice(firms, size=N_firms, replace=True)
        # Reconstruct panel: keep all years for each resampled firm
        # Handle duplicates by assigning new ids
        pieces = []
        for i, fid in enumerate(boot_firms):
            chunk = df.loc[df['id'] == fid].copy()
            chunk['boot_id'] = i  # unique id for clustering
            pieces.append(chunk)
        boot_df = pd.concat(pieces, ignore_index=True)
        boot_df['id'] = boot_df['boot_id']

        try:
            beta_b, _ = ols_premium(boot_df)
            beta_boot[b] = beta_b
        except Exception:
            pass

    elapsed = time.time() - t0
    valid = beta_boot[np.isfinite(beta_boot)]
    boot_se = float(np.std(valid, ddof=1))
    q_lo = float(np.percentile(valid, 2.5))
    q_hi = float(np.percentile(valid, 97.5))
    pct_valid = 100 * len(valid) / B

    if theta_hat is not None and pivotal_ci:
        ci_lo = float(2.0 * theta_hat - q_hi)
        ci_hi = float(2.0 * theta_hat - q_lo)
        ci_method = "pivotal"
    else:
        ci_lo = q_lo
        ci_hi = q_hi
        ci_method = "percentile"

    print(f'  Bootstrap: {len(valid)}/{B} valid ({pct_valid:.0f}%), '
          f'{elapsed:.1f}s')
    print(f'  Boot SE = {boot_se:.4f}, 95% CI ({ci_method}) '
          f'= [{ci_lo:.4f}, {ci_hi:.4f}]')

    out = {
        'boot_se': boot_se, 'ci_lo': ci_lo, 'ci_hi': ci_hi,
        'n_valid': len(valid), 'elapsed': elapsed,
    }
    if theta_hat is not None:
        out['bias'] = float(np.mean(valid) - theta_hat)
        out['bias_corrected'] = float(2.0 * theta_hat - np.mean(valid))
    return out


# ================================================================== #
#  Main
# ================================================================== #

def main():
    print('=== Misspecification Diagnostics (Andrews, Chen, Tecchio 2026) ===')

    # ---- Load data ---- #
    data_path = INPUT_DIR / 'data_rebuilt.dta'
    if not data_path.exists():
        data_path = INPUT_DIR / 'data.dta'
    df = pd.read_stata(str(data_path))
    if 'ico' in df.columns and 'id' not in df.columns:
        df = df.rename(columns={'ico': 'id'})
    df = df.dropna(subset=['go', 'k', 'cogs', 'pp_dummy', 'year', 'nace2']).copy()
    print(f'Data: N={len(df)}, firms={df["id"].nunique()}')

    # ---- Part 1: J-statistic diagnostics ---- #
    print('\n--- Part 1: Hansen J / sqrt(J) Weight-Sensitivity Bound ---')

    j_results = []

    # Spec A: translog, pp in Markov, NACE 41
    df_41 = df[df['nace2'] == 41].copy()
    form_a = Formulation(spec='tl', overidentify=True,
                         pp_in_markov=True, pp_interactions=True,
                         year_fe=True, nace2_fe=False)
    est_a = ACFEstimator(
        data=df_41, formulation=form_a,
        optimization=Optimization(method='nm+bfgs'),
        extensions=CWDLExtensions(survival_correction=True),
        n_starts=3,
    )
    res_a = est_a.solve()
    j_results.append(compute_j_diagnostic(est_a, res_a.betas,
                                          'Spec A (TL, NACE 41)'))

    # Spec E: CD base, pp in Markov, NACE 41
    form_e = Formulation(spec='cd', pp_in_markov=True, pp_interactions=True,
                         year_fe=True, nace2_fe=False)
    est_e = ACFEstimator(
        data=df_41, formulation=form_e,
        optimization=Optimization(method='nm+bfgs'),
        extensions=CWDLExtensions(survival_correction=True),
        n_starts=3,
    )
    res_e = est_e.solve()
    j_results.append(compute_j_diagnostic(est_e, res_e.betas,
                                          'Spec E (CD, NACE 41)'))

    # Pooled construction (all NACE, with nace2 FE)
    form_pool = Formulation(spec='tl', overidentify=True,
                            pp_in_markov=True, pp_interactions=True,
                            year_fe=True, nace2_fe=True)
    est_pool = ACFEstimator(
        data=df, formulation=form_pool,
        optimization=Optimization(method='nm+bfgs'),
        extensions=CWDLExtensions(survival_correction=True),
        n_starts=3,
    )
    res_pool = est_pool.solve()
    j_results.append(compute_j_diagnostic(est_pool, res_pool.betas,
                                          'Pooled (CD, all NACE)'))

    j_df = pd.DataFrame(j_results)

    # ---- Part 2: Non-recentered cluster bootstrap ---- #
    print('\n--- Part 2: Non-Recentered Cluster Bootstrap ---')

    mk = pd.read_stata(str(OUTPUT_DIR / 'data' / 'paper_markups.dta'))
    if hasattr(mk['year'].iloc[0], 'year'):
        mk['year'] = mk['year'].dt.year
    mk['ln_mu'] = np.log(mk['markup_A'])
    mk['D'] = mk['pp_dummy'].astype(int)
    mk = mk.dropna(subset=['ln_mu', 'D'])
    print(f'Premium panel: N={len(mk)}, firms={mk["id"].nunique()}')

    # Analytical clustered SE
    beta_hat, se_analytical = ols_premium(mk)
    print(f'  Analytical: beta_pp = {beta_hat:.4f}, SE = {se_analytical:.4f}')

    # Bootstrap SE — pass theta_hat to enable Conlon's pivotal CI
    boot = bootstrap_premium_se(mk, B=B_BOOT, theta_hat=beta_hat)
    se_ratio = boot['boot_se'] / se_analytical if se_analytical > 0 else np.nan

    print(f'\n  SE comparison:')
    print(f'    Analytical (clustered): {se_analytical:.4f}')
    print(f'    Bootstrap (non-recentered): {boot["boot_se"]:.4f}')
    print(f'    Ratio (boot/analytical): {se_ratio:.3f}')
    print(f'    95% CI: [{boot["ci_lo"]:.4f}, {boot["ci_hi"]:.4f}]')

    # Add SE info to j_df rows for table output
    se_rows = []
    for _, row in j_df.iterrows():
        se_rows.append({
            **row.to_dict(),
            'se_analytical': se_analytical,
            'se_bootstrap': boot['boot_se'],
            'se_ratio': se_ratio,
        })
    combined = pd.DataFrame(se_rows)

    # ---- Save CSV ---- #
    combined.to_csv(OUTPUT_DIR / 'misspecification_diagnostics.csv', index=False)
    print(f'\nSaved: misspecification_diagnostics.csv')

    # ---- LaTeX table ---- #
    (OUTPUT_DIR / 'tables').mkdir(parents=True, exist_ok=True)

    def fmt(x, d=2):
        if pd.isna(x) or not np.isfinite(x):
            return '---'
        return f'{x:.{d}f}'

    tex = [
        r'\begin{table}[htbp]\centering',
        r'\caption{Misspecification Diagnostics (Andrews, Chen, and Tecchio 2026)}'
        r'\label{tab:misspec}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lcccccc}',
        r'\toprule',
        r'& \multicolumn{3}{c}{J-statistic} '
        r'& \multicolumn{3}{c}{Premium SE} \\',
        r'\cmidrule(lr){2-4}\cmidrule(lr){5-7}',
        r'Specification & $J$ & $\sqrt{J}$ & $p$-value '
        r'& Analytical & Bootstrap & Ratio \\',
        r'\midrule',
    ]
    for _, r in combined.iterrows():
        label = r['Specification'].replace('_', r'\_')
        tex.append(
            f'{label} & {fmt(r["J"])} & {fmt(r["sqrt_J"])} '
            f'& {fmt(r["p_value"])} '
            f'& {fmt(r["se_analytical"], 4)} '
            f'& {fmt(r["se_bootstrap"], 4)} '
            f'& {fmt(r["se_ratio"])} \\\\'
        )
    tex += [
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r'\item \textit{Notes:} $J$ is Hansen\textquotesingle s overidentification '
        r'statistic. Under local misspecification, $\sqrt{J}$ bounds the number '
        r'of standard errors by which the estimate can be moved by varying the '
        r'GMM weight matrix (Andrews, Chen, and Tecchio 2026, Proposition 5.2). '
        r'Bootstrap SE: 500 cluster-level replications without moment recentering '
        r'(valid under misspecification). Analytical SE: firm-clustered sandwich '
        r'(valid under correct specification only).',
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}',
    ]
    tex_path = OUTPUT_DIR / 'tables' / 'misspecification_diagnostics.tex'
    with open(tex_path, 'w') as f:
        f.write('\n'.join(tex))
    print(f'Saved: tables/misspecification_diagnostics.tex')

    print('\n=== Done ===')


if __name__ == '__main__':
    main()
