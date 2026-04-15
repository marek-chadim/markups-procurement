"""dml_core.py — shared infrastructure for Double/Debiased Machine Learning robustness.

Reference implementation adapted from Chernozhukov, Hansen, Kallus, Spindler, and
Syrgkanis (2022), *Applied Causal Inference Powered by ML and AI*, notebook AC1
(`python-sensitivity-analysis-with-sensmakr-and-debiased-ml.ipynb`), cell 19 —
the canonical `dml()` function for partially linear model estimation with
cross-fitting and clustered standard errors.

Used by:
    dml_premium.py            Partially linear DML for the 14% premium
    dml_sensitivity.py        Chernozhukov-Cinelli-Hazlett-Kuchibhotla-Robins 2019 sensitivity
    dml_cate.py               DML CATE by NACE and reform era
    dml_iv.py                 DML partially linear IV model
    dml_strong_exclusion.py   Orthogonal moments for ABGRS strong exclusion

The partially linear model is

    Y = α · D + g(X) + U,   E[U | D, X] = 0

and the Neyman-orthogonal score (Chernozhukov et al. 2018 Econometrics Journal)

    ψ(W; α, g, m) = (Y - g(X)) · (D - m(X)) - α · (D - m(X))²

where g(X) = E[Y|X] and m(X) = E[D|X] are cross-fitted ML nuisances. The
point estimate is the residual-on-residual OLS

    α̂ = Σ (Y_i - ĝ(X_i))(D_i - m̂(X_i)) / Σ (D_i - m̂(X_i))²

with standard errors clustered by firm via statsmodels.

Panel treatment (Chernozhukov-Newey-Singh 2022): Y, D, and continuous X are
pre-demeaned by firm before cross-fitting. Year and NACE fixed effects are
included as one-hot dummies inside X.

References:
    - Chernozhukov et al. 2018 Econometrics Journal 21(1): C1–C68
    - Chernozhukov-Newey-Singh 2022 JoE
    - MetricsMLNotebooks/AC1/python-sensitivity-analysis-with-sensmakr-and-debiased-ml.ipynb
"""

from __future__ import annotations
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LassoCV, LogisticRegressionCV
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import statsmodels.formula.api as smf

# ---- Paths -----------------------------------------------------------------
# This module lives in 2_analysis/source/lib/; climb 2 levels to reach
# 2_analysis/, then into input/ or output/. Previously .parent was enough
# because the module lived at source/ root; after the 2026-04-15 lib/
# refactor it's one level deeper.
LIB_DIR = Path(__file__).resolve().parent
MODULE_DIR = LIB_DIR.parent.parent
INPUT_DIR = MODULE_DIR / "input"
OUTPUT_DIR = MODULE_DIR / "output"
TAB_DIR = OUTPUT_DIR / "tables"
FIG_DIR = OUTPUT_DIR / "figures"
DAT_DIR = OUTPUT_DIR / "data"
for d in (TAB_DIR, FIG_DIR, DAT_DIR):
    d.mkdir(parents=True, exist_ok=True)

PAPER_MARKUPS = OUTPUT_DIR / "data" / "paper_markups.dta"
DATA_REBUILT = INPUT_DIR / "data_rebuilt.dta"

# ---- Defaults --------------------------------------------------------------
DEFAULT_SEED = 42
DEFAULT_NFOLDS = 5

# Covariate set for nuisance X matrix. Two critical exclusions:
#
# 1. Variables that are mechanically a function of pp_dummy (pp_share,
#    pp_stock_*, pp_dummy_L1, single_bid_share, avg_bids, hhi_revenue, entry) —
#    these would let the propensity nuisance m(X)=E[D|X] near-perfectly predict
#    D, zeroing out the orthogonal score.
#
# 2. Production inputs k and cogs — log(markup_A) is derived from the ACF formula
#    μ = θ_cogs · (go/cogs), so log(markup) is a mechanical function of cogs.
#    Including cogs in X lets the outcome nuisance g(X)=E[Y|X] absorb the
#    treatment effect non-parametrically.
#
# The safe covariates below are pre-treatment firm characteristics that are
# NEITHER mechanical functions of the treatment NOR of the ACF markup formula:
# market structure (mktshare) and ownership (foreign). omega_A is available
# separately for the productivity-decomposition spec (see dml_premium.py).
BASE_COVARIATES = [
    "mktshare",           # firm sales / industry sales — exogenous market structure
    "foreign",            # foreign ownership dummy — exogenous
]

# Productivity decomposition spec: adds omega_A as a confounder. Replicates the
# paper's finding that controlling for ω̂_it reduces the premium to ~3.2%
# (line 341 of markups_procurement.tex).
PRODUCTIVITY_COVARIATES = BASE_COVARIATES + ["omega_A"]


# ---------------------------------------------------------------------------
# Data loading and preparation
# ---------------------------------------------------------------------------

def load_merged_panel() -> pd.DataFrame:
    """Merge paper_markups.dta with data_rebuilt.dta on (id, year).

    Returns a single panel DataFrame with firm-year covariates for DML.
    Inner join on non-null markup_A; expected N ≈ 7,666.
    """
    print(f"[load] {PAPER_MARKUPS.name}")
    mk = pd.read_stata(PAPER_MARKUPS)
    mk = mk.loc[mk["markup_A"].notna()].copy()
    mk["id"] = mk["id"].astype(int)
    mk["year"] = mk["year"].astype(int)

    print(f"[load] {DATA_REBUILT.name}")
    rb = pd.read_stata(DATA_REBUILT)
    rb["id"] = rb["id"].astype(int)
    rb["year"] = rb["year"].astype(int)

    # Drop duplicate column names from the rebuilt panel that already exist in mk
    drop_cols = [c for c in ("nace2", "pp_dummy", "k", "cogs", "go")
                 if c in rb.columns and c in mk.columns]
    rb_merge = rb.drop(columns=drop_cols)

    merged = mk.merge(rb_merge, on=["id", "year"], how="inner",
                       validate="one_to_one")
    merged["log_mu"] = np.log(merged["markup_A"])
    print(f"  merged panel: {len(merged):,} firm-year observations, "
          f"{merged['id'].nunique():,} firms, "
          f"years {merged['year'].min()}–{merged['year'].max()}")
    return merged


def construct_X(df: pd.DataFrame, covariates: list[str] | None = None,
                 fe_policy: str = "demean") -> pd.DataFrame:
    """Build the nuisance matrix X for DML.

    Parameters
    ----------
    df : DataFrame with the merged panel
    covariates : list of column names to include (defaults to BASE_COVARIATES)
    fe_policy : 'dummy' for one-hot year/NACE dummies (default); 'demean' for
                within-firm demeaning of all continuous covariates

    Returns
    -------
    DataFrame of shape (N, p) with numeric columns only. NaNs replaced with
    column median for numerical stability.
    """
    if covariates is None:
        covariates = BASE_COVARIATES
    # De-duplicate and keep columns that exist
    cols = []
    for c in covariates:
        if c not in cols and c in df.columns:
            cols.append(c)

    X = df[cols].copy()
    for c in X.columns:
        X[c] = pd.to_numeric(X[c], errors="coerce")
        med = X[c].median()
        X[c] = X[c].fillna(med if pd.notna(med) else 0.0)

    # One-hot year dummies
    year_d = pd.get_dummies(df["year"], prefix="year", drop_first=True).astype(float)
    X = pd.concat([X, year_d], axis=1)

    # One-hot NACE dummies
    nace_d = pd.get_dummies(df["nace2"], prefix="nace", drop_first=True).astype(float)
    X = pd.concat([X, nace_d], axis=1)

    if fe_policy == "demean":
        # Within-firm demean continuous covariates (keep dummies as-is).
        # Chernozhukov-Newey-Singh 2022 panel DML with firm FE.
        cont_cols = [c for c in cols if c in X.columns]
        X_demean = X.copy()
        firm = df["id"].values
        for c in cont_cols:
            means = pd.Series(X[c].values).groupby(firm).transform("mean").values
            X_demean[c] = X[c].values - means
        X = X_demean

    return X.astype(float)


def firm_demean(series: pd.Series, firm: pd.Series) -> np.ndarray:
    """Within-firm demean a series. Reusable helper."""
    s = pd.Series(np.asarray(series, dtype=float))
    means = s.groupby(firm.values).transform("mean").values
    return s.values - means


# ---------------------------------------------------------------------------
# Cross-fitted nuisance estimation (Chernozhukov et al 2018)
# ---------------------------------------------------------------------------

def cross_fit(X: np.ndarray | pd.DataFrame, target: np.ndarray,
              model, *, n_folds: int = DEFAULT_NFOLDS,
              classifier: bool = False, seed: int = DEFAULT_SEED) -> np.ndarray:
    """Cross-fitted out-of-fold predictions for one ML model.

    Mirrors MetricsMLNotebooks/AC1 cell 19 `dml()` function's cross-fitting step.
    For binary treatment nuisance use classifier=True to return predict_proba[:,1].
    """
    cv = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    X = np.asarray(X, dtype=float)
    target = np.asarray(target)
    if classifier:
        y_hat = cross_val_predict(model, X, target, cv=cv,
                                    method="predict_proba", n_jobs=-1)[:, 1]
    else:
        y_hat = cross_val_predict(model, X, target, cv=cv, n_jobs=-1)
    return np.asarray(y_hat)


# ---------------------------------------------------------------------------
# Orthogonal-score estimator (partially linear model)
# ---------------------------------------------------------------------------

def plr_orthogonal(y: np.ndarray, d: np.ndarray, g_hat: np.ndarray,
                    m_hat: np.ndarray, cluster: np.ndarray) -> dict:
    """Compute α̂ from the Neyman-orthogonal score with cluster SEs.

    Parameters
    ----------
    y, d : outcome and treatment (numpy arrays)
    g_hat : cross-fitted predictions of E[Y|X]
    m_hat : cross-fitted predictions of E[D|X]
    cluster : firm ids for clustered SEs

    Returns
    -------
    dict with: point, stderr, n, t, ci_lo, ci_hi, resY, resD
    """
    resY = y - g_hat
    resD = d - m_hat
    df = pd.DataFrame({"resY": resY, "resD": resD, "cluster": cluster})
    ols = smf.ols("resY ~ 1 + resD", data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df["cluster"]}
    )
    point = float(ols.params["resD"])
    stderr = float(ols.bse["resD"])
    ci = ols.conf_int().loc["resD"].tolist()
    return dict(
        point=point,
        stderr=stderr,
        n=int(len(y)),
        t=point / stderr if stderr > 0 else float("nan"),
        ci_lo=float(ci[0]),
        ci_hi=float(ci[1]),
        resY=resY,
        resD=resD,
        ols=ols,
    )


# ---------------------------------------------------------------------------
# Orthogonal IV score (partially linear IV model)
# ---------------------------------------------------------------------------

def plr_iv_orthogonal(y: np.ndarray, d: np.ndarray, z: np.ndarray,
                        g_hat: np.ndarray, m_hat_d: np.ndarray,
                        m_hat_z: np.ndarray, cluster: np.ndarray) -> dict:
    """Compute α̂_IV for the partially linear IV model with cluster SEs.

    Neyman orthogonal IV score: α̂ = Σ(Y - g̃)(Z - m̃_Z) / Σ(D - m̃_D)(Z - m̃_Z)
    Implemented as 2SLS of resY on resD instrumented by resZ.
    """
    resY = y - g_hat
    resD = d - m_hat_d
    resZ = z - m_hat_z
    # First-stage partial R²
    first_stage = np.polyfit(resZ, resD, 1)
    resD_hat = first_stage[0] * resZ + first_stage[1]
    ss_res = np.sum((resD - resD_hat) ** 2)
    ss_tot = np.sum((resD - resD.mean()) ** 2)
    partial_r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # IV estimator: α̂ = cov(resY, resZ) / cov(resD, resZ)
    cov_yz = np.mean((resY - resY.mean()) * (resZ - resZ.mean()))
    cov_dz = np.mean((resD - resD.mean()) * (resZ - resZ.mean()))
    point = cov_yz / cov_dz if abs(cov_dz) > 1e-12 else float("nan")

    # Clustered SE via influence function: ψ = (Y - g̃ - α(D - m̃_D))(Z - m̃_Z) / cov_dz
    psi = (resY - point * resD) * resZ / cov_dz
    df = pd.DataFrame({"psi": psi, "cluster": cluster})
    # Clustered variance estimator: sum of squared cluster means
    cluster_sums = df.groupby("cluster")["psi"].sum()
    var_cluster = (cluster_sums ** 2).sum() / len(y) ** 2
    stderr = float(np.sqrt(var_cluster))

    return dict(
        point=float(point),
        stderr=stderr,
        partial_r2=float(partial_r2),
        first_stage_coef=float(first_stage[0]),
        n=int(len(y)),
        t=point / stderr if stderr > 0 else float("nan"),
        ci_lo=point - 1.96 * stderr,
        ci_hi=point + 1.96 * stderr,
        resY=resY,
        resD=resD,
        resZ=resZ,
    )


# ---------------------------------------------------------------------------
# Weak-IV-robust Anderson-Rubin confidence set for PLR-IV
# ---------------------------------------------------------------------------

def ar_confidence_set(resY: np.ndarray, resD: np.ndarray, resZ: np.ndarray,
                       cluster: np.ndarray, alpha: float = 0.05,
                       beta_grid: np.ndarray | None = None) -> dict:
    """Cluster-robust Anderson-Rubin confidence set for the single-IV PLR-IV
    model, inverted over a grid of candidate coefficients.

    For each candidate $\\beta_0$ the moment vector is
    $m_i(\\beta_0) = (resY_i - \\beta_0 \\cdot resD_i) \\cdot resZ_i$,
    and the cluster-robust AR statistic is
    $AR(\\beta_0) = \\bar{m}(\\beta_0)^2 / \\widehat{Var}_{cluster}(\\bar{m})$
    where the cluster variance is computed as the sum of squared within-cluster
    moment totals divided by $N^2$. Under $H_0: \\beta = \\beta_0$ the statistic
    is $\\chi^2_1$-distributed regardless of instrument strength, so the set
    $\\{\\beta_0 : AR(\\beta_0) \\le \\chi^2_{1, 1-\\alpha}\\}$ is a valid
    $(1-\\alpha)$ confidence set (Anderson and Rubin 1949; Stock and Wright 2000).

    Returns a dict with keys: `ar_lo`, `ar_hi` (the grid endpoints of the
    accepted region), `ar_empty` (True if no grid point is accepted),
    `ar_min_stat` (the minimum AR stat over the grid, should be near zero
    at the IV point estimate), and `ar_n_accepted` (count of accepted
    grid points).
    """
    from scipy.stats import chi2
    if beta_grid is None:
        beta_grid = np.linspace(-0.5, 1.5, 401)
    crit = chi2.ppf(1 - alpha, df=1)

    resY = np.asarray(resY)
    resD = np.asarray(resD)
    resZ = np.asarray(resZ)
    cluster = np.asarray(cluster)
    n = len(resY)

    accepted = []
    stats = []
    for b in beta_grid:
        m = (resY - b * resD) * resZ
        df = pd.DataFrame({"m": m, "cluster": cluster})
        cluster_sums = df.groupby("cluster")["m"].sum().values
        var_cluster = float((cluster_sums ** 2).sum() / (n ** 2))
        mean_sq = float(m.mean() ** 2)
        ar_stat = mean_sq / var_cluster if var_cluster > 0 else float("inf")
        stats.append(ar_stat)
        if ar_stat < crit:
            accepted.append(float(b))

    stats = np.array(stats)
    if not accepted:
        return dict(
            ar_lo=float("nan"), ar_hi=float("nan"), ar_empty=True,
            ar_min_stat=float(stats.min()) if len(stats) else float("nan"),
            ar_n_accepted=0,
        )
    return dict(
        ar_lo=min(accepted), ar_hi=max(accepted), ar_empty=False,
        ar_min_stat=float(stats.min()),
        ar_n_accepted=len(accepted),
    )


# ---------------------------------------------------------------------------
# Sensitivity analysis helpers (Chernozhukov-Cinelli-Hazlett-Kuchibhotla-Robins 2019)
# ---------------------------------------------------------------------------

def sensitivity_bias(point: float, resY: np.ndarray, resD: np.ndarray,
                      r2_yc: float, r2_dc: float) -> float:
    """Compute the absolute bias from a hypothetical confounder with
    partial R² of (r2_yc, r2_dc) in the (outcome, treatment) equations.

    Formula (AC1 notebook cell 17, Chernozhukov et al 2019):
        Bias² = [R²_YC × R²_DC / (1 - R²_DC)] × (var(epsilon) / var(resD))
    where epsilon is the PLR residual after the final residual-on-residual OLS.
    """
    kappa = (r2_yc * r2_dc) / (1 - r2_dc) if r2_dc < 1 else float("inf")
    # epsilon residual of resY regressed on resD
    alpha = np.sum(resY * resD) / np.sum(resD ** 2)
    epsilon = resY - alpha * resD
    variance_ratio = np.mean(epsilon ** 2) / np.mean(resD ** 2)
    return float(np.sqrt(kappa * variance_ratio))


def robustness_value(point: float, stderr: float, resY: np.ndarray,
                       resD: np.ndarray, q: float = 0.0,
                       alpha: float = 0.05) -> float:
    """Robustness value (RV) — minimum equal partial R² in both equations
    at which the estimate loses significance at level α.

    Cinelli-Hazlett (2020) formula generalized to DML (Chernozhukov et al 2019).
    We solve for the RV that makes |bias| = |point - q · stderr · 1.96|.

    Parameters
    ----------
    point : point estimate
    stderr : standard error
    resY, resD : residualized outcome and treatment
    q : fraction of the stderr to leave — q=0 means bias flips sign
    alpha : significance level for the critical value
    """
    from scipy.optimize import brentq
    target_bias = abs(point) - q * stderr * 1.96
    if target_bias <= 0:
        return 0.0

    def bias_equal(r):
        if r <= 0 or r >= 1:
            return -target_bias
        return sensitivity_bias(point, resY, resD, r, r) - target_bias

    try:
        rv = brentq(bias_equal, 1e-6, 0.9999, xtol=1e-6)
    except ValueError:
        # Bias never reaches target at any R² ∈ (0, 1) — very robust
        rv = 1.0
    return float(rv)


# ---------------------------------------------------------------------------
# ML estimator factories
# ---------------------------------------------------------------------------

def make_outcome_estimators(seed: int = DEFAULT_SEED) -> dict:
    """Return the three default outcome-nuisance ML models for E[Y|X]."""
    return {
        "Lasso": make_pipeline(StandardScaler(), LassoCV(cv=5, random_state=seed, max_iter=10_000)),
        "RF":    RandomForestRegressor(n_estimators=500, max_depth=8,
                                         min_samples_leaf=5, random_state=seed,
                                         n_jobs=-1),
        "GB":    GradientBoostingRegressor(n_estimators=500, max_depth=4,
                                             learning_rate=0.05, random_state=seed),
    }


def make_treatment_estimators(seed: int = DEFAULT_SEED,
                                 binary: bool = True) -> dict:
    """Return the three default treatment-nuisance ML models for E[D|X].

    Binary treatment → classifier with predict_proba.
    Continuous treatment → regressor.
    """
    if binary:
        return {
            "Lasso": make_pipeline(StandardScaler(),
                                     LogisticRegressionCV(cv=5, random_state=seed,
                                                            max_iter=10_000,
                                                            penalty="l2")),
            "RF":    RandomForestClassifier(n_estimators=500, max_depth=8,
                                              min_samples_leaf=5, random_state=seed,
                                              n_jobs=-1),
            "GB":    GradientBoostingClassifier(n_estimators=500, max_depth=4,
                                                  learning_rate=0.05, random_state=seed),
        }
    return make_outcome_estimators(seed)


# ---------------------------------------------------------------------------
# Firm-level cluster bootstrap
# ---------------------------------------------------------------------------

def cluster_bootstrap(point_fn: Callable, firms: np.ndarray,
                        n_rep: int = 999, seed: int = DEFAULT_SEED,
                        theta_hat: Optional[float] = None,
                        pivotal_ci: bool = True) -> dict:
    """Block bootstrap by firm cluster.

    Parameters
    ----------
    point_fn : Callable
        Takes an index array and returns a point estimate.
    firms : np.ndarray
        Cluster ids (typically firm ids).
    n_rep : int
        Bootstrap replications. Default 999 (Conlon bootstrap.tex convention
        for percentile CIs; pivotal CIs are less sensitive to small B).
    seed : int
        RNG seed.
    theta_hat : Optional[float]
        Full-sample point estimate. When provided AND ``pivotal_ci=True``
        (the default), pivotal CIs are returned along with bias-corrected
        point estimate. When ``None``, falls back to percentile CIs
        (backward-compatible with pre-2026-04-15 behavior).
    pivotal_ci : bool
        If True (default) and ``theta_hat`` is provided, return
        Conlon's "Better Way" CI ``[2*theta_hat - q97.5, 2*theta_hat - q2.5]``
        (bootstrap.tex lines 75-82). If False, return naive percentile CI
        ``[q2.5, q97.5]`` even when ``theta_hat`` is provided.

    Returns
    -------
    dict with keys:
        se : sample SD of bootstrap reps (with ddof=1)
        ci_lo, ci_hi : 95% CI (pivotal if ``theta_hat`` given and
            ``pivotal_ci=True``, percentile otherwise)
        n_rep : number of valid replicates
        bias : bootstrap bias estimate ``mean(reps) - theta_hat``
            (only if ``theta_hat`` is provided)
        bias_corrected : ``2*theta_hat - mean(reps)`` (only if
            ``theta_hat`` is provided; Conlon bootstrap.tex line 61)
    """
    rng = np.random.default_rng(seed)
    uniq_firms = np.unique(firms)
    firm_to_idx = {f: np.where(firms == f)[0] for f in uniq_firms}

    reps = []
    for _ in range(n_rep):
        boot_firms = rng.choice(uniq_firms, size=len(uniq_firms), replace=True)
        boot_idx = np.concatenate([firm_to_idx[f] for f in boot_firms])
        try:
            reps.append(point_fn(boot_idx))
        except Exception:
            continue
    reps = np.asarray(reps)
    n_valid = int(len(reps))

    # Conlon "When Does It Fail?" frame (bootstrap.tex lines 102-108): if too
    # many reps are silently dropped, the empirical distribution may not
    # approximate the true sampling distribution.
    if n_rep > 0 and n_valid / n_rep < 0.95:
        print(f"  WARNING: cluster_bootstrap dropped "
              f"{n_rep - n_valid}/{n_rep} reps "
              f"({100*(1 - n_valid/n_rep):.0f}%)")

    q_lo = float(np.quantile(reps, 0.025))
    q_hi = float(np.quantile(reps, 0.975))

    out = dict(
        se=float(np.std(reps, ddof=1)),
        n_rep=n_valid,
    )

    if theta_hat is not None and pivotal_ci:
        # Conlon bootstrap.tex lines 75-82: pivotal CI carries automatic
        # bias correction and the higher-order Edgeworth refinement that's
        # the reason to prefer bootstrap over delta method.
        out["ci_lo"] = float(2.0 * theta_hat - q_hi)
        out["ci_hi"] = float(2.0 * theta_hat - q_lo)
    else:
        # Naive percentile CI (Conlon's "obvious way"): backward-compatible
        # path for callers that don't supply theta_hat.
        out["ci_lo"] = q_lo
        out["ci_hi"] = q_hi

    if theta_hat is not None:
        out["bias"] = float(reps.mean() - theta_hat)
        out["bias_corrected"] = float(2.0 * theta_hat - reps.mean())

    return out


# ---------------------------------------------------------------------------
# Utility: stars for significance
# ---------------------------------------------------------------------------

def stars(p: float) -> str:
    if p < 0.01:
        return "***"
    if p < 0.05:
        return "**"
    if p < 0.10:
        return "*"
    return ""


if __name__ == "__main__":
    # Smoke test: verify load + cross-fit works end-to-end
    print("[dml_core] smoke test")
    df = load_merged_panel()
    X = construct_X(df, fe_policy="demean")
    print(f"  X shape: {X.shape}, fe_policy=demean (within-firm)")
    # Demean y and d by firm as well — panel DML with firm FE
    y = firm_demean(df["log_mu"], df["id"])
    d = firm_demean(df["pp_dummy"].astype(float), df["id"])
    firm = df["id"].values

    outcome_models = make_outcome_estimators()
    g_hat = cross_fit(X, y, outcome_models["Lasso"])
    print(f"  g_hat stats: mean={g_hat.mean():.4f}, sd={g_hat.std():.4f}, "
          f"corr(y,g_hat)={np.corrcoef(y, g_hat)[0,1]:.4f}")
    treatment_models = make_treatment_estimators(binary=False)
    m_hat = cross_fit(X, d, treatment_models["Lasso"], classifier=False)
    print(f"  m_hat stats: mean={m_hat.mean():.4f}, sd={m_hat.std():.4f}")

    res = plr_orthogonal(y, d, g_hat, m_hat, firm)
    print(f"  DML-Lasso premium (demeaned): α̂ = {res['point']:.4f} "
          f"(SE {res['stderr']:.4f}, t = {res['t']:.2f}, "
          f"95% CI [{res['ci_lo']:.4f}, {res['ci_hi']:.4f}])")
    print("  Target (OLS firm+year FE baseline): 0.1324")
    print("[dml_core] smoke test OK")
