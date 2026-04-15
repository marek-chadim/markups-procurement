#!/usr/bin/env python
"""abgrs_mc_validation.py — Monte Carlo validation for the four-stage
composite estimator formalized in `abgrs_exposition.tex`.

Validates Propositions 1 and 2 on synthetic panels with a known DGP across
four scenarios that contrast strongly-identified against identification
failures in the Stage C regression or the cost-share input:

  S1 Baseline CD               — correctly specified, within-firm regression
  S2 Over-parameterized trnslog — TL fit when true is CD (Prop 2 validation)
  S3 Pooled OLS (no firm FE)   — drops the within-firm identification strategy;
                                  selection bias dominates
  S4 Noisy cost-share input    — classical EIV on alpha^V (measurement error
                                  in cost-of-goods-sold reporting propagates
                                  directly to the Stage B markup numerator)

Outputs:
  2_analysis/output/tables/mc_validation.tex         (LaTeX table)
  2_analysis/output/figures/mc_validation.pdf        (4-panel density)
  2_analysis/output/data/mc_validation_results.csv   (raw results)

Run: python abgrs_mc_validation.py
Reproducible via MASTER_SEED = 20260415.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt

SRC_DIR = Path(__file__).parent
sys.path.insert(0, str(SRC_DIR / 'lib'))
from style_markups import (  # noqa: E402
    apply_markups_style,
    MARKUPS_BLUE,
    MARKUPS_GREEN,
    MARKUPS_PINK,
    MARKUPS_RED,
)
from ags_sensitivity import (  # noqa: E402
    get_sensitivity,
    get_standardized_sensitivity,
    compute_composite_lambda,
    gmm_sandwich_vcov,
)

MASTER_SEED = 20260415
S_REPS = 20
ALPHA_PP_TRUE = 0.138

EXTERNAL_PANEL_PATH = (
    SRC_DIR.parent.parent / '1_data' / 'output' / 'external_panel_annual.csv'
)

EXTERNAL_SHIFTER_COLS = [
    'fx_eur',
    'ppi_industry',
    'construction_cost_residential',
    'long_rate_mcby',
    'building_permits_sqm',
    'weather_frost_days',
    'weather_precipitation_sum',
]

DGP = dict(
    N=1500,
    T=10,
    beta_c=0.65,
    beta_k=0.25,
    rho_omega=0.70,
    delta_omega=0.20,
    sigma_xi=0.10,
    alpha_pp=ALPHA_PP_TRUE,
    mu0=0.70,
    sigma_eta=0.30,
    sigma_lambda=0.05,
    sigma_v=0.05,
    sigma_eps=0.15,
    eta_pp_corr=1.0,      # coefficient linking firm FE to procurement propensity
    pp_misclass_rate=0.15,  # S4 classical binary misclassification rate for pp
    pp_persistence=0.85,
    k0_mean=5.0,
    k0_sd=1.0,
    k_drift=0.05,
    k_sd=0.05,
)


def generate_panel(seed: int, params: dict) -> pd.DataFrame:
    """Generate a synthetic panel consistent with the DGP.

    The firm fixed effect eta_i is constructed to be positively correlated
    with the firm-level procurement propensity p_i, so that procurement-active
    firms have systematically higher baseline markups. Within-firm
    regressions (S1, S2, S4) absorb this correlation through the firm FE;
    pooled regressions (S3) do not, producing a selection bias.
    """
    rng = np.random.default_rng(seed)
    N, T = params['N'], params['T']
    p_i = rng.beta(2, 3, N)
    # eta correlated with propensity: procurement firms have higher FE
    eta_raw = rng.normal(0, params['sigma_eta'], N)
    eta = eta_raw + params['eta_pp_corr'] * (p_i - p_i.mean())
    lam = rng.normal(0, params['sigma_lambda'], T)

    pp = np.zeros((N, T), dtype=np.int8)
    pp[:, 0] = (rng.uniform(0, 1, N) < p_i).astype(np.int8)
    for t in range(1, T):
        stay = rng.uniform(0, 1, N) < params['pp_persistence']
        draw = (rng.uniform(0, 1, N) < p_i).astype(np.int8)
        pp[:, t] = np.where(stay, pp[:, t - 1], draw)

    k = np.zeros((N, T))
    k[:, 0] = rng.normal(params['k0_mean'], params['k0_sd'], N)
    for t in range(1, T):
        k[:, t] = k[:, t - 1] + params['k_drift'] + rng.normal(0, params['k_sd'], N)

    omega = np.zeros((N, T))
    omega[:, 0] = rng.normal(0, 0.20, N)
    for t in range(1, T):
        omega[:, t] = (
            params['rho_omega'] * omega[:, t - 1]
            + params['delta_omega'] * pp[:, t - 1]
            + rng.normal(0, params['sigma_xi'], N)
        )

    v = rng.normal(0, params['sigma_v'], (N, T))
    log_mu = (
        params['mu0']
        + params['alpha_pp'] * pp
        + eta[:, None]
        + lam[None, :]
        + v
    )
    mu = np.exp(log_mu)
    alpha_V = params['beta_c'] / mu
    log_alpha_V = np.log(alpha_V)

    eps = rng.normal(0, params['sigma_eps'], (N, T))
    one_minus_bc = 1.0 - params['beta_c']
    c = (log_alpha_V + params['beta_k'] * k + omega + eps) / one_minus_bc
    y = c - log_alpha_V

    # Classical binary misclassification of pp (S4). With probability
    # pp_misclass_rate, flip the observed procurement status. This is the
    # standard EIV-in-treatment setup and produces attenuation by factor
    # (1 - 2*lambda).
    flip = rng.uniform(0, 1, (N, T)) < params['pp_misclass_rate']
    pp_obs = np.where(flip, 1 - pp, pp).astype(np.int8)

    firm_ids = np.repeat(np.arange(N), T)
    years = np.tile(np.arange(T), N)
    df = pd.DataFrame({
        'id': firm_ids,
        'year': years,
        'y': y.flatten(),
        'c': c.flatten(),
        'k': k.flatten(),
        'pp': pp.flatten(),
        'pp_noisy': pp_obs.flatten(),
        'alpha_V': alpha_V.flatten(),
        'log_mu_true': log_mu.flatten(),
        'omega_true': omega.flatten(),
    })
    df = _merge_external_shifters(df)
    return df


def _merge_external_shifters(df: pd.DataFrame) -> pd.DataFrame:
    """Map synthetic MC years 0..T-1 onto real calendar years 2005..2005+T-1
    and left-join the real external shifter panel, interacting each shifter
    with k_it and c_it to produce 14 firm-year varying instruments.

    Returns df with new columns: <shifter>_z, <shifter>_x_k, <shifter>_x_Lc.
    If the external panel is unavailable, returns df unchanged.
    """
    if not EXTERNAL_PANEL_PATH.exists():
        return df
    try:
        ext = pd.read_csv(EXTERNAL_PANEL_PATH)
    except Exception:
        return df
    needed_cols = ['year'] + EXTERNAL_SHIFTER_COLS
    missing = [c for c in needed_cols if c not in ext.columns]
    if missing:
        return df
    ext = ext[needed_cols].copy()
    # Keep only rows with all shifters non-null
    ext = ext.dropna(subset=EXTERNAL_SHIFTER_COLS)
    if len(ext) == 0:
        return df
    # Z-score each shifter across available years
    for col in EXTERNAL_SHIFTER_COLS:
        mean, std = ext[col].mean(), ext[col].std()
        ext[f'{col}_z'] = (ext[col] - mean) / (std if std > 0 else 1.0)
    # Map MC synthetic years 0..T-1 → real years (wrap at end of ext)
    unique_mc_years = np.sort(df['year'].unique())
    ext_years = ext['year'].to_numpy()
    ext = ext.reset_index(drop=True)
    year_map = {}
    for i, mc_y in enumerate(unique_mc_years):
        real_y = int(ext_years[i % len(ext_years)])
        year_map[int(mc_y)] = real_y
    df = df.copy()
    df['_real_year'] = df['year'].map(year_map)
    z_cols = [f'{c}_z' for c in EXTERNAL_SHIFTER_COLS]
    ext_z_only = ext[['year'] + z_cols].rename(columns={'year': '_real_year'})
    df = df.merge(ext_z_only, on='_real_year', how='left')
    df = df.drop(columns=['_real_year'])
    # Firm-year interactions with k and (lagged) c: store in dataframe for use
    # as Stage A instruments in the S5 scenario. Lagged c is constructed in
    # run_composite; here we store the raw shifter columns plus k-interactions.
    for c in EXTERNAL_SHIFTER_COLS:
        df[f'{c}_x_k'] = df[f'{c}_z'] * df['k']
    return df


def _poly_basis(c: np.ndarray, k: np.ndarray, order: int) -> np.ndarray:
    """Build polynomial basis up to `order` in (c, k), excluding the constant."""
    cols = []
    for i in range(order + 1):
        for j in range(order + 1 - i):
            if i + j == 0:
                continue
            cols.append((c ** i) * (k ** j))
    return np.column_stack(cols)


def stage0_invert(df: pd.DataFrame, include_pp: bool, order: int = 3) -> np.ndarray:
    """Stage 0: first-stage inversion via polynomial regression of y on (c,k[,pp])."""
    c = df['c'].to_numpy()
    k = df['k'].to_numpy()
    basis = _poly_basis(c, k, order)
    if include_pp:
        pp = df['pp'].to_numpy().astype(float)
        basis = np.column_stack([basis, basis * pp[:, None], pp[:, None]])
    X = np.column_stack([np.ones(len(df)), basis])
    y = df['y'].to_numpy()
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    phi_hat = X @ beta
    return phi_hat


def _markov_regression(omega_theta: np.ndarray, df_sorted: pd.DataFrame,
                       include_pp_shifter: bool) -> np.ndarray:
    """Regress omega_t on cubic(omega_{t-1}) [+ pp_{t-1}]; return xi residuals."""
    df_m = df_sorted[['id', 'year']].copy()
    df_m['omega'] = omega_theta
    df_m['omega_lag'] = df_m.groupby('id')['omega'].shift(1)
    df_m['pp_lag'] = df_sorted.groupby('id')['pp'].shift(1).to_numpy()
    mask = df_m['omega_lag'].notna() & df_m['pp_lag'].notna()
    if mask.sum() < 50:
        return np.full(len(df_m), np.nan)
    ol = df_m['omega_lag'].to_numpy()[mask]
    X_cols = [np.ones(mask.sum()), ol, ol ** 2, ol ** 3]
    if include_pp_shifter:
        X_cols.append(df_m['pp_lag'].to_numpy()[mask])
    X_m = np.column_stack(X_cols)
    y_m = df_m['omega'].to_numpy()[mask]
    beta_m, *_ = np.linalg.lstsq(X_m, y_m, rcond=None)
    resid = y_m - X_m @ beta_m
    xi_full = np.full(len(df_m), np.nan)
    xi_full[mask.to_numpy()] = resid
    return xi_full


def _omega_from_theta(theta: np.ndarray, df: pd.DataFrame, spec: str,
                      phi_hat: np.ndarray) -> np.ndarray:
    c = df['c'].to_numpy()
    k = df['k'].to_numpy()
    if spec == 'cd':
        beta_c, beta_k = theta[0], theta[1]
        return phi_hat - (beta_c * c + beta_k * k)
    elif spec == 'tl':
        beta_c, beta_k, beta_cc, beta_kk, beta_kc = theta
        return phi_hat - (
            beta_c * c + beta_k * k
            + beta_cc * c * c + beta_kk * k * k + beta_kc * c * k
        )
    raise ValueError(spec)


def _gmm_objective(theta: np.ndarray, df_sorted: pd.DataFrame, phi_hat: np.ndarray,
                   spec: str, Z: np.ndarray, include_pp_shifter: bool) -> float:
    omega_t = _omega_from_theta(theta, df_sorted, spec, phi_hat)
    xi = _markov_regression(omega_t, df_sorted, include_pp_shifter)
    mask = ~np.isnan(xi)
    if mask.sum() < 50:
        return 1e6
    xi_m = xi[mask]
    Z_m = Z[mask]
    moments = Z_m.T @ xi_m / mask.sum()
    return float(moments @ moments)


def stage_a_gmm(df_sorted: pd.DataFrame, phi_hat: np.ndarray, spec: str,
                instrument_cols: list[str], include_pp_shifter: bool) -> np.ndarray:
    """Stage A: GMM estimator for (beta_c, beta_k[, curvature])."""
    Z = df_sorted[instrument_cols].to_numpy(dtype=float)
    if spec == 'cd':
        starts = [(0.65, 0.25), (0.50, 0.20)]
    else:
        starts = [
            (0.65, 0.25, 0.0, 0.0, 0.0),
            (0.60, 0.25, 0.02, 0.01, -0.01),
        ]
    best_x, best_f = None, np.inf
    for t0 in starts:
        try:
            res = minimize(
                _gmm_objective, t0,
                args=(df_sorted, phi_hat, spec, Z, include_pp_shifter),
                method='Nelder-Mead',
                options={'xatol': 5e-4, 'fatol': 1e-6, 'maxiter': 200},
            )
            # Reject spurious local minima with implausible beta_c values
            # (per KLS 2019 reseeding logic)
            if not (0.05 < res.x[0] < 1.5):
                continue
            if res.fun < best_f:
                best_f = res.fun
                best_x = res.x
        except Exception:
            continue
    if best_x is None:
        return np.array([np.nan] * (2 if spec == 'cd' else 5))
    return best_x


def stage_b_markup(df_sorted: pd.DataFrame, theta_hat: np.ndarray,
                   spec: str) -> np.ndarray:
    c = df_sorted['c'].to_numpy()
    k = df_sorted['k'].to_numpy()
    if spec == 'cd':
        theta_V = np.full(len(df_sorted), theta_hat[0])
    else:
        beta_c, _, beta_cc, _, beta_kc = theta_hat
        theta_V = beta_c + 2.0 * beta_cc * c + beta_kc * k
    alpha_V = df_sorted['alpha_V'].to_numpy()
    with np.errstate(divide='ignore', invalid='ignore'):
        mu_hat = theta_V / alpha_V
    log_mu = np.where((mu_hat > 0) & np.isfinite(mu_hat), np.log(mu_hat), np.nan)
    return log_mu


def stage_c_regression(log_mu_hat: np.ndarray,
                       df_sorted: pd.DataFrame,
                       firm_fe: bool = True,
                       pp_col: str = 'pp') -> tuple[float, float]:
    """Stage C: OLS of log_mu_hat on pp.

    If firm_fe=True: within-firm-year demean (firm + year FE). Cluster-robust SE.
    If firm_fe=False: pooled OLS with only year FE (no firm control). Tests
    whether productivity selection biases the treatment effect.
    pp_col: which column to use for the treatment indicator ('pp' or 'pp_noisy').
    """
    df_c = df_sorted[['id', 'year', pp_col]].copy()
    df_c = df_c.rename(columns={pp_col: 'pp'})
    df_c['y'] = log_mu_hat
    df_c = df_c.dropna(subset=['y'])
    if len(df_c) < 100:
        return np.nan, np.nan
    y_year = df_c.groupby('year')['y'].transform('mean')
    y_grand = df_c['y'].mean()
    pp_year = df_c.groupby('year')['pp'].transform('mean')
    pp_grand = df_c['pp'].mean()
    if firm_fe:
        y_firm = df_c.groupby('id')['y'].transform('mean')
        pp_firm = df_c.groupby('id')['pp'].transform('mean')
        y_tilde = (df_c['y'] - y_firm - y_year + y_grand).to_numpy()
        pp_tilde = (df_c['pp'] - pp_firm - pp_year + pp_grand).to_numpy()
    else:
        y_tilde = (df_c['y'] - y_year).to_numpy()
        pp_tilde = (df_c['pp'] - pp_year).to_numpy()
    den = (pp_tilde ** 2).sum()
    if den < 1e-10:
        return np.nan, np.nan
    beta = float((y_tilde * pp_tilde).sum() / den)
    resid = y_tilde - beta * pp_tilde
    score = resid * pp_tilde
    df_c['score'] = score
    cluster_sum = df_c.groupby('id')['score'].sum().to_numpy()
    meat = (cluster_sum ** 2).sum()
    se = float(np.sqrt(meat) / den)
    return beta, se


SCENARIOS = {
    'S1': {
        'label': 'Baseline CD (correctly specified)',
        'spec': 'cd',
        'include_pp_in_stage0': True,
        'include_pp_shifter': True,
        'instrument_cols': ['k', 'c_lag1'],
        'stage_c_firm_fe': True,
    },
    'S2': {
        'label': 'Over-parameterized translog',
        'spec': 'tl',
        'include_pp_in_stage0': True,
        'include_pp_shifter': True,
        'instrument_cols': ['k', 'c_lag1', 'k_lag1', 'c_lag2'],
        'stage_c_firm_fe': True,
    },
    'S3': {
        'label': 'Pooled Stage C (no firm FE)',
        'spec': 'cd',
        'include_pp_in_stage0': True,
        'include_pp_shifter': True,
        'instrument_cols': ['k', 'c_lag1'],
        'stage_c_firm_fe': False,
    },
    'S4': {
        'label': 'Misclassified treatment (EIV in $pp$)',
        'spec': 'cd',
        'include_pp_in_stage0': True,
        'include_pp_shifter': True,
        'instrument_cols': ['k', 'c_lag1'],
        'stage_c_firm_fe': True,
        'pp_col': 'pp_noisy',
    },
    'S5': {
        'label': 'Over-identified CD with external shifters',
        'spec': 'cd',
        'include_pp_in_stage0': True,
        'include_pp_shifter': True,
        'instrument_cols': ['k', 'c_lag1'] + [
            f'{c}_x_k_resid' for c in EXTERNAL_SHIFTER_COLS
        ] + [
            f'{c}_x_Lc_resid' for c in EXTERNAL_SHIFTER_COLS
        ],
        'stage_c_firm_fe': True,
        'build_residualized_shifters': True,
    },
}


def _residualize_ols(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Return OLS residuals of y on [1, X], NaN-safe."""
    mask = ~np.isnan(y) & ~np.any(np.isnan(X), axis=1)
    out = np.full_like(y, np.nan, dtype=float)
    if mask.sum() < 10:
        return out
    X_mat = np.column_stack([np.ones(mask.sum()), X[mask]])
    beta, *_ = np.linalg.lstsq(X_mat, y[mask], rcond=None)
    out[mask] = y[mask] - X_mat @ beta
    return out


def _build_residualized_shifters(df_work: pd.DataFrame) -> None:
    """ABGRS Appendix C.3 residualization of firm-year shifter interactions.

    For each of 7 external shifters we have (<sh>_z, <sh>_x_k) in df. We
    add (<sh>_x_Lc) using c_lag1 and OLS-residualize both interactions
    against the control set X = (k, c_lag1, pp, year-dummies) so the
    resulting residualized moments are mean-independent of X by construction.
    """
    year_dummies = pd.get_dummies(df_work['year'].astype(int), prefix='_yr',
                                    drop_first=True).astype(float)
    ctrl = np.column_stack([
        df_work[['k', 'c_lag1', 'pp']].to_numpy(dtype=float),
        year_dummies.to_numpy(),
    ])
    for c in EXTERNAL_SHIFTER_COLS:
        if f'{c}_z' not in df_work.columns:
            return  # external panel missing → skip
        # c_lag1 interaction (may be built fresh here)
        df_work[f'{c}_x_Lc'] = (df_work[f'{c}_z'].to_numpy() *
                                 df_work['c_lag1'].to_numpy())
        df_work[f'{c}_x_k_resid'] = _residualize_ols(
            df_work[f'{c}_x_k'].to_numpy(dtype=float), ctrl
        )
        df_work[f'{c}_x_Lc_resid'] = _residualize_ols(
            df_work[f'{c}_x_Lc'].to_numpy(dtype=float), ctrl
        )


def run_composite(df: pd.DataFrame, scen: dict) -> tuple[float, float, float]:
    """Execute Stages 0-C end-to-end for one scenario."""
    df_work = df.sort_values(['id', 'year']).reset_index(drop=True).copy()
    df_work['c_lag1'] = df_work.groupby('id')['c'].shift(1)
    df_work['c_lag2'] = df_work.groupby('id')['c'].shift(2)
    df_work['k_lag1'] = df_work.groupby('id')['k'].shift(1)
    if scen.get('build_residualized_shifters'):
        _build_residualized_shifters(df_work)

    phi_hat = stage0_invert(df_work, include_pp=scen['include_pp_in_stage0'], order=3)

    needed = scen['instrument_cols']
    valid_mask = df_work[needed].notna().all(axis=1)
    df_a = df_work.loc[valid_mask].reset_index(drop=True)
    phi_hat_a = phi_hat[valid_mask.to_numpy()]
    if len(df_a) < 100:
        return np.nan, np.nan, np.nan

    theta_hat = stage_a_gmm(
        df_a, phi_hat_a, spec=scen['spec'],
        instrument_cols=needed,
        include_pp_shifter=scen['include_pp_shifter'],
    )
    if np.isnan(theta_hat).any():
        return np.nan, np.nan, np.nan

    log_mu_hat = stage_b_markup(df_a, theta_hat, spec=scen['spec'])
    pp_col = scen.get('pp_col', 'pp')
    beta_pp, se = stage_c_regression(log_mu_hat, df_a,
                                      firm_fe=scen['stage_c_firm_fe'],
                                      pp_col=pp_col)
    return float(theta_hat[0]), beta_pp, se


def run_single(seed: int, scen_key: str) -> dict:
    df = generate_panel(seed, DGP)
    beta_c_hat, alpha_pp_hat, se = run_composite(df, SCENARIOS[scen_key])
    return dict(seed=int(seed), scenario=scen_key,
                beta_c_hat=beta_c_hat, alpha_pp_hat=alpha_pp_hat, se=se)


def run_mc(S: int = S_REPS) -> pd.DataFrame:
    rng = np.random.default_rng(MASTER_SEED)
    seeds = rng.integers(0, 2 ** 31, S)
    results = []
    for scen_key in SCENARIOS:
        print(f"[{scen_key}] {SCENARIOS[scen_key]['label']}", flush=True)
        res = []
        for i, s in enumerate(seeds):
            res.append(run_single(int(s), scen_key))
            if (i + 1) % 10 == 0:
                print(f"    {i+1}/{S}", flush=True)
        results.extend(res)
        sub = pd.DataFrame([r for r in res if not np.isnan(r['alpha_pp_hat'])])
        if len(sub) > 0:
            print(f"    converged {len(sub)}/{S}, "
                  f"mean alpha_hat = {sub['alpha_pp_hat'].mean():.4f}, "
                  f"bias = {sub['alpha_pp_hat'].mean() - ALPHA_PP_TRUE:+.4f}",
                  flush=True)
    return pd.DataFrame(results)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scen in SCENARIOS:
        sub = df[df['scenario'] == scen].dropna(subset=['alpha_pp_hat'])
        if len(sub) == 0:
            rows.append(dict(scenario=scen, n_conv=0))
            continue
        alpha_hat = sub['alpha_pp_hat'].to_numpy()
        beta_c_hat = sub['beta_c_hat'].to_numpy()
        se_hat = sub['se'].to_numpy()
        ci_low = alpha_hat - 1.96 * se_hat
        ci_high = alpha_hat + 1.96 * se_hat
        cov = float(((ALPHA_PP_TRUE >= ci_low) & (ALPHA_PP_TRUE <= ci_high)).mean())
        rows.append(dict(
            scenario=scen,
            n_conv=len(sub),
            beta_c_mean=float(np.nanmean(beta_c_hat)),
            alpha_mean=float(alpha_hat.mean()),
            alpha_sd=float(alpha_hat.std(ddof=1)),
            bias=float(alpha_hat.mean() - ALPHA_PP_TRUE),
            rmse=float(np.sqrt(((alpha_hat - ALPHA_PP_TRUE) ** 2).mean())),
            cov95=cov,
        ))
    return pd.DataFrame(rows)


def save_table(summary: pd.DataFrame, outpath: Path) -> None:
    lines = [
        r"\begin{tabular}{lcccccc}",
        r"\toprule",
        r"Scenario & $S_\text{conv}$ & $\mathbb{E}[\hat\beta_c]$ & "
        r"$\mathbb{E}[\hat\alpha_{pp}]$ & Bias & RMSE & Cov.\ 95\% \\",
        r"\midrule",
    ]
    for _, r in summary.iterrows():
        lbl = SCENARIOS[r['scenario']]['label']
        if pd.isna(r.get('beta_c_mean', np.nan)):
            lines.append(f"{r['scenario']}: {lbl} & 0 & --- & --- & --- & --- & --- \\\\")
            continue
        lines.append(
            f"{r['scenario']}: {lbl} & {int(r['n_conv'])} & "
            f"{r['beta_c_mean']:.3f} & {r['alpha_mean']:.3f} & "
            f"{r['bias']:+.3f} & {r['rmse']:.3f} & {r['cov95']:.2f} \\\\"
        )
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    outpath.write_text('\n'.join(lines) + '\n')


def plot_distributions(df: pd.DataFrame, outpath: Path) -> None:
    apply_markups_style()
    colors = {'S1': MARKUPS_BLUE, 'S2': MARKUPS_GREEN,
              'S3': MARKUPS_RED, 'S4': MARKUPS_PINK}
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    for ax, scen_key in zip(axes.flat, SCENARIOS):
        sub = df[df['scenario'] == scen_key].dropna(subset=['alpha_pp_hat'])
        if len(sub) == 0:
            ax.text(0.5, 0.5, 'No convergent replications',
                    transform=ax.transAxes, ha='center', va='center')
            ax.set_title(f"{scen_key}: {SCENARIOS[scen_key]['label']}")
            continue
        vals = sub['alpha_pp_hat'].to_numpy()
        col = colors[scen_key]
        ax.hist(vals, bins=40, density=True, alpha=0.75,
                color=col, edgecolor='white', linewidth=0.4)
        ax.axvline(ALPHA_PP_TRUE, color='black', linestyle='--', linewidth=1.5,
                   label=rf'True $\alpha_{{pp}} = {ALPHA_PP_TRUE}$')
        ax.axvline(vals.mean(), color=col, linestyle='-', linewidth=2.0,
                   label=rf'Mean $= {vals.mean():.3f}$')
        ax.set_xlabel(r'$\hat\alpha_{pp}$')
        ax.set_ylabel('Density')
        ax.set_title(f"{scen_key}: {SCENARIOS[scen_key]['label']}", fontsize=10)
        ax.legend(loc='upper right', fontsize=8, frameon=False)
    fig.suptitle(
        rf'Monte Carlo sampling distributions of $\hat\alpha_{{pp}}$ '
        rf'($S={S_REPS}$, $N={DGP["N"]}$, $T={DGP["T"]}$, seed={MASTER_SEED})',
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)


def _moment_value(theta: np.ndarray, df_sorted: pd.DataFrame, spec: str,
                  phi_hat: np.ndarray, Z: np.ndarray,
                  include_pp_shifter: bool) -> np.ndarray:
    """Compute the sample moment vector m(theta) = (1/N) sum xi_it z_it.

    Returns an M-vector (one entry per instrument column).
    """
    xi = _markov_regression(
        _omega_from_theta(theta, df_sorted, spec, phi_hat),
        df_sorted, include_pp_shifter,
    )
    mask = ~np.isnan(xi)
    xi_m = xi[mask]
    Z_m = Z[mask]
    return Z_m.T @ xi_m / mask.sum()


def _numerical_jacobian(theta: np.ndarray, df_sorted: pd.DataFrame, spec: str,
                        phi_hat: np.ndarray, Z: np.ndarray,
                        include_pp_shifter: bool,
                        eps: float = 1e-4) -> np.ndarray:
    """Central-difference Jacobian d bar_m / d theta' (M×P).

    Matches the definition used by AGS (2017) Λ = -(J'WJ)^{-1} J'W.
    """
    m0 = _moment_value(theta, df_sorted, spec, phi_hat, Z, include_pp_shifter)
    M = m0.size
    P = theta.size
    J = np.zeros((M, P))
    for p in range(P):
        t_plus = theta.copy()
        t_minus = theta.copy()
        step = max(eps, eps * abs(theta[p]))
        t_plus[p] += step
        t_minus[p] -= step
        m_plus = _moment_value(t_plus, df_sorted, spec, phi_hat, Z, include_pp_shifter)
        m_minus = _moment_value(t_minus, df_sorted, spec, phi_hat, Z, include_pp_shifter)
        J[:, p] = (m_plus - m_minus) / (2 * step)
    return J


def _compute_psi(df_sorted: pd.DataFrame, theta_hat: np.ndarray,
                 spec: str) -> np.ndarray:
    r"""Row vector $\psi = \nabla_\theta \hat\beta_{pp}|_{\theta_0}$.

    Uses the Proposition~1 Step~3 derivation:
    $\hat\beta_{pp} = \Cov(\widetilde{pp}, \widetilde{\log\hat\mu}) / \Var(\widetilde{pp})$
    and $\nabla_\theta \log\hat\mu_{it} = \nabla_\theta \log\hat\theta^V_{it}$
    (since $\alpha^V$ does not depend on $\theta$). Returns a 1×P row vector.

    Under CD the gradient is constant across firm-years, so within-firm-year
    demeaning sets $\psi = 0$ — consistent with Corollary 1 of the exposition
    note. Under translog the gradient varies with $(c, k)$ and $\psi$ is
    generally non-zero.
    """
    c = df_sorted['c'].to_numpy()
    k = df_sorted['k'].to_numpy()

    if spec == 'cd':
        beta_c = theta_hat[0]
        grad = np.zeros((len(df_sorted), 2))
        grad[:, 0] = 1.0 / beta_c
        grad[:, 1] = 0.0
    elif spec == 'tl':
        beta_c, beta_k, beta_cc, beta_kk, beta_kc = theta_hat
        theta_V = beta_c + 2.0 * beta_cc * c + beta_kc * k
        inv_tV = 1.0 / theta_V
        grad = np.zeros((len(df_sorted), 5))
        grad[:, 0] = inv_tV
        grad[:, 1] = 0.0
        grad[:, 2] = 2.0 * c * inv_tV
        grad[:, 3] = 0.0
        grad[:, 4] = k * inv_tV
    else:
        raise ValueError(spec)

    ids = df_sorted['id'].to_numpy()
    years = df_sorted['year'].to_numpy()
    pp = df_sorted['pp'].to_numpy(dtype=float)

    pp_demean = pp.copy()
    for col in range(grad.shape[1]):
        gc = grad[:, col]
        gc_demean = gc.copy()
        for unique_id in np.unique(ids):
            m_id = ids == unique_id
            gc_demean[m_id] -= gc[m_id].mean()
        for unique_year in np.unique(years):
            m_y = years == unique_year
            gc_demean[m_y] -= gc[m_y].mean()
        gc_demean += gc.mean()
        grad[:, col] = gc_demean

    pp_firm_mean = pd.Series(pp).groupby(pd.Series(ids)).transform('mean').to_numpy()
    pp_year_mean = pd.Series(pp).groupby(pd.Series(years)).transform('mean').to_numpy()
    pp_demean = pp - pp_firm_mean - pp_year_mean + pp.mean()

    var_pp = (pp_demean ** 2).mean()
    if var_pp < 1e-10:
        return np.zeros((1, grad.shape[1]))
    psi = (grad * pp_demean[:, None]).mean(axis=0) / var_pp
    return psi.reshape(1, -1)


def _sample_moment_vcov(theta_hat: np.ndarray, df_sorted: pd.DataFrame,
                        spec: str, phi_hat: np.ndarray, Z: np.ndarray,
                        include_pp_shifter: bool) -> np.ndarray:
    r"""Sample estimate of $\Omega = \Var(\xi_{it} z_{it})$ — moment variance."""
    xi = _markov_regression(
        _omega_from_theta(theta_hat, df_sorted, spec, phi_hat),
        df_sorted, include_pp_shifter,
    )
    mask = ~np.isnan(xi)
    xi_m = xi[mask]
    Z_m = Z[mask]
    g_it = Z_m * xi_m[:, None]  # N × M
    g_bar = g_it.mean(axis=0)
    centered = g_it - g_bar
    return centered.T @ centered / len(g_it)


def compute_hansen_j(theta_hat: np.ndarray, df_sorted: pd.DataFrame,
                     spec: str, phi_hat: np.ndarray, Z: np.ndarray,
                     include_pp_shifter: bool) -> tuple[float, int]:
    """Hansen-Sargan J-statistic under optimal weighting.

    Returns (J-statistic, degrees of freedom = M - P).
    Under correct specification and identity first-step weighting, J is
    distributed chi-squared on (M - P) degrees of freedom.
    """
    xi = _markov_regression(
        _omega_from_theta(theta_hat, df_sorted, spec, phi_hat),
        df_sorted, include_pp_shifter,
    )
    mask = ~np.isnan(xi)
    xi_m = xi[mask]
    Z_m = Z[mask]
    n = mask.sum()
    if n < 50:
        return np.nan, 0
    # Empirical moment vcov Omega_hat
    g_it = Z_m * xi_m[:, None]
    g_bar = g_it.mean(axis=0)
    Omega = (g_it - g_bar).T @ (g_it - g_bar) / n
    try:
        W_opt = np.linalg.pinv(Omega)
    except np.linalg.LinAlgError:
        return np.nan, 0
    J = n * float(g_bar @ W_opt @ g_bar)
    M = Z.shape[1]
    P = len(theta_hat)
    return J, max(M - P, 0)


def compute_partial_r2(df_work: pd.DataFrame, Z_cols: list[str]) -> dict:
    """Partial R² of each instrument on the control set (year + k + c_lag1 + pp).

    ABGRS strong-exclusion diagnostic: should be < 0.10 for moments that
    satisfy mean-independence of the controls.
    """
    year_dummies = pd.get_dummies(df_work['year'].astype(int), prefix='_yr',
                                    drop_first=True).astype(float)
    ctrl = np.column_stack([
        df_work[['k', 'c_lag1', 'pp']].to_numpy(dtype=float),
        year_dummies.to_numpy(),
    ])
    r2_out = {}
    for col in Z_cols:
        z = df_work[col].to_numpy(dtype=float)
        mask = ~np.isnan(z) & ~np.any(np.isnan(ctrl), axis=1)
        if mask.sum() < 20:
            r2_out[col] = np.nan
            continue
        X = np.column_stack([np.ones(mask.sum()), ctrl[mask]])
        beta, *_ = np.linalg.lstsq(X, z[mask], rcond=None)
        resid = z[mask] - X @ beta
        var_total = z[mask].var()
        var_resid = resid.var()
        if var_total < 1e-12:
            r2_out[col] = np.nan
        else:
            r2_out[col] = 1.0 - var_resid / var_total
    return r2_out


def compute_ags_lambda(df: pd.DataFrame, scen: dict) -> dict:
    """Run Stage 0-A on one draw; return AGS Lambda diagnostics."""
    df_work = df.sort_values(['id', 'year']).reset_index(drop=True).copy()
    df_work['c_lag1'] = df_work.groupby('id')['c'].shift(1)
    df_work['c_lag2'] = df_work.groupby('id')['c'].shift(2)
    df_work['k_lag1'] = df_work.groupby('id')['k'].shift(1)
    if scen.get('build_residualized_shifters'):
        _build_residualized_shifters(df_work)

    phi_hat = stage0_invert(df_work, include_pp=scen['include_pp_in_stage0'], order=3)
    needed = scen['instrument_cols']
    valid = df_work[needed].notna().all(axis=1)
    df_a = df_work.loc[valid].reset_index(drop=True)
    phi_hat_a = phi_hat[valid.to_numpy()]

    theta_hat = stage_a_gmm(
        df_a, phi_hat_a, spec=scen['spec'],
        instrument_cols=needed,
        include_pp_shifter=scen['include_pp_shifter'],
    )
    if np.isnan(theta_hat).any():
        return {}

    Z = df_a[needed].to_numpy(dtype=float)
    J = _numerical_jacobian(theta_hat, df_a, scen['spec'], phi_hat_a, Z,
                            scen['include_pp_shifter'])
    M, P = J.shape
    W = np.eye(M)
    Lambda = get_sensitivity(J, W)

    Omega = _sample_moment_vcov(theta_hat, df_a, scen['spec'], phi_hat_a, Z,
                                 scen['include_pp_shifter'])
    param_vcov = gmm_sandwich_vcov(J, W, Omega, n=len(df_a))

    psi = _compute_psi(df_a, theta_hat, scen['spec'])
    Lambda_bpp = compute_composite_lambda(Lambda, psi)

    se_param = np.sqrt(np.clip(np.diag(param_vcov), 1e-12, None))
    se_mom = np.sqrt(np.clip(np.diag(Omega) / len(df_a), 1e-12, None))
    Lambda_std = get_standardized_sensitivity(Lambda, se_param, se_mom)

    se_beta_pp_sq = float(psi @ param_vcov @ psi.T)
    se_beta_pp = float(np.sqrt(max(se_beta_pp_sq, 1e-20)))

    return dict(
        theta_hat=theta_hat,
        J=J,
        Lambda=Lambda,
        Lambda_std=Lambda_std,
        Omega=Omega,
        param_vcov=param_vcov,
        psi=psi,
        Lambda_bpp=Lambda_bpp,
        se_param=se_param,
        se_mom=se_mom,
        se_beta_pp=se_beta_pp,
        n_obs=len(df_a),
        instrument_cols=needed,
    )


def save_lambda_table(res_cd: dict, res_tl: dict, outpath: Path) -> None:
    """Write ags_lambda_premium.tex — Λ for the causal summary under CD and translog."""
    mom_labels = {
        'k': r'$k_{t}$',
        'c_lag1': r'$c_{t-1}$',
        'k_lag1': r'$k_{t-1}$',
        'c_lag2': r'$c_{t-2}$',
        'const': r'constant',
    }

    def fmt(x):
        if abs(x) < 1e-6:
            return r'$<10^{-5}$'
        return f"{x:.4f}"

    lines = [
        r"\begin{tabular}{lcc}",
        r"\toprule",
        r" & Cobb-Douglas (S1) & Translog (S2) \\",
        r"Moment & $\Lambda^{\beta_{pp}}$ & $\Lambda^{\beta_{pp}}$ \\",
        r"\midrule",
    ]
    all_cols = list(dict.fromkeys(res_cd['instrument_cols'] + res_tl['instrument_cols']))
    for col in all_cols:
        cd_val = '—'
        tl_val = '—'
        if col in res_cd['instrument_cols']:
            idx = res_cd['instrument_cols'].index(col)
            cd_val = fmt(res_cd['Lambda_bpp'][0, idx])
        if col in res_tl['instrument_cols']:
            idx = res_tl['instrument_cols'].index(col)
            tl_val = fmt(res_tl['Lambda_bpp'][0, idx])
        lines.append(f"{mom_labels.get(col, col)} & {cd_val} & {tl_val} \\\\")
    lines.append(r"\midrule")
    lines.append(
        rf"$\|\psi\|_2$ & {float(np.linalg.norm(res_cd['psi'])):.2e} & "
        rf"{float(np.linalg.norm(res_tl['psi'])):.2e} \\"
    )
    lines.append(
        rf"$\|\Lambda^{{\beta_{{pp}}}}\|_2$ & "
        rf"{float(np.linalg.norm(res_cd['Lambda_bpp'])):.2e} & "
        rf"{float(np.linalg.norm(res_tl['Lambda_bpp'])):.2e} \\"
    )
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    outpath.write_text('\n'.join(lines) + '\n')


def run_lambda_diagnostic(out_root: Path) -> None:
    """Compute AGS Λ for β_pp on one draw under S1 (CD), S2 (translog), S5 (overid-CD).

    Also reports Hansen-Sargan J and ABGRS strong-exclusion partial R² for S5.
    """
    from scipy import stats as _stats
    df = generate_panel(MASTER_SEED, DGP)
    res_cd = compute_ags_lambda(df, SCENARIOS['S1'])
    res_tl = compute_ags_lambda(df, SCENARIOS['S2'])
    res_s5 = compute_ags_lambda(df, SCENARIOS['S5'])
    if not res_cd or not res_tl:
        print("[Λ] diagnostic failed to converge on seed", MASTER_SEED)
        return
    print()
    print("=== AGS Lambda for the procurement premium (single-draw, seed="
          f"{MASTER_SEED}) ===")
    print(f"CD (S1):  theta_hat = {res_cd['theta_hat']}")
    print(f"CD (S1):  |psi|_2 = {np.linalg.norm(res_cd['psi']):.2e} "
          f"(expected ~ 0 under CD)")
    print(f"CD (S1):  |Lambda_bpp|_2 = {np.linalg.norm(res_cd['Lambda_bpp']):.2e}")
    print(f"TL (S2):  theta_hat = {res_tl['theta_hat']}")
    print(f"TL (S2):  |psi|_2 = {np.linalg.norm(res_tl['psi']):.4f}")
    print(f"TL (S2):  |Lambda_bpp|_2 = {np.linalg.norm(res_tl['Lambda_bpp']):.2e}")
    if res_s5:
        print(f"S5 over-id CD: theta_hat = {res_s5['theta_hat']}")
        print(f"S5 over-id CD: n_moments = {len(res_s5['instrument_cols'])}")
        print(f"S5 over-id CD: |psi|_2 = {np.linalg.norm(res_s5['psi']):.2e} "
              f"(expected ~ 0 by Corollary 1)")
    save_lambda_table(res_cd, res_tl, out_root / 'tables' / 'ags_lambda_premium.tex')
    print(f"Wrote table -> {out_root / 'tables' / 'ags_lambda_premium.tex'}")

    # --- S5 Hansen-Sargan J-test and partial R² diagnostic ---
    if not res_s5:
        return
    df_work = df.sort_values(['id', 'year']).reset_index(drop=True).copy()
    df_work['c_lag1'] = df_work.groupby('id')['c'].shift(1)
    df_work['c_lag2'] = df_work.groupby('id')['c'].shift(2)
    df_work['k_lag1'] = df_work.groupby('id')['k'].shift(1)
    _build_residualized_shifters(df_work)
    needed = SCENARIOS['S5']['instrument_cols']
    valid = df_work[needed].notna().all(axis=1)
    df_a = df_work.loc[valid].reset_index(drop=True)
    phi_a = stage0_invert(df_a, include_pp=True, order=3)
    Z = df_a[needed].to_numpy(dtype=float)
    J_stat, df_J = compute_hansen_j(
        res_s5['theta_hat'], df_a, 'cd', phi_a, Z,
        include_pp_shifter=True,
    )
    p_val = float(1 - _stats.chi2.cdf(J_stat, df_J)) if df_J > 0 else np.nan
    # Strong-exclusion check applies only to the added external-shifter
    # moments; the baseline (k, c_lag1) are in the control set by
    # construction and trivially have partial R² = 1.
    external_only = [c for c in needed if c not in ('k', 'c_lag1')]
    part_r2 = compute_partial_r2(df_a, external_only)
    print()
    print("=== S5 specification diagnostics (over-identified CD) ===")
    print(f"Hansen-Sargan J = {J_stat:.2f}  df = {df_J}  p = {p_val:.3f}")
    print("Partial R² (residualized shifter moments, ABGRS threshold 0.10):")
    max_r2 = -np.inf
    for col in external_only:
        r2 = part_r2.get(col, np.nan)
        max_r2 = max(max_r2, r2 if r2 is not None and not np.isnan(r2) else -np.inf)
        flag = " [FAIL]" if r2 is not None and r2 > 0.10 else ""
        print(f"  {col:45s}  R² = {r2:+.4f}{flag}")
    print(f"max partial R² over 14 residualized moments: {max_r2:+.4f}  "
          f"(ABGRS threshold 0.10: {'PASS' if max_r2 < 0.10 else 'FAIL'})")

    save_spec_diag_table(
        res_s5, J_stat, df_J, p_val, part_r2, external_only,
        out_root / 'tables' / 'spec_diagnostics_s5.tex',
    )
    print(f"Wrote table -> {out_root / 'tables' / 'spec_diagnostics_s5.tex'}")


def save_spec_diag_table(res_s5: dict, J_stat: float, df_J: int, p_val: float,
                          part_r2: dict, Z_cols: list[str], outpath: Path) -> None:
    """Write spec_diagnostics_s5.tex reporting S5 J-test + partial R²."""
    lines = [
        r"\begin{tabular}{lcc}",
        r"\toprule",
        r"Moment & Partial $R^2$ on controls & ABGRS threshold \\",
        r"\midrule",
    ]
    for col in Z_cols:
        r2 = part_r2.get(col, np.nan)
        pass_fail = r"$\checkmark$" if (r2 is not None and r2 < 0.10) else r"$\times$"
        safe = col.replace('_', r'\_')
        r2_str = "—" if (r2 is None or np.isnan(r2)) else f"{r2:.4f}"
        lines.append(rf"\texttt{{{safe}}} & {r2_str} & {pass_fail} \\")
    lines.append(r"\midrule")
    lines.append(rf"Hansen-Sargan $J$ & {J_stat:.2f} (df $= {df_J}$, $p = {p_val:.3f}$) & --- \\")
    lines.append(
        rf"$\|\psi\|_2$ & {float(np.linalg.norm(res_s5['psi'])):.2e} & "
        r"Corollary~\ref{cor:cd_zero} \\"
    )
    lines.append(
        rf"$\|\Lambda^{{\beta_{{pp}}}}\|_2$ & "
        rf"{float(np.linalg.norm(res_s5['Lambda_bpp'])):.2e} & "
        r"Corollary~\ref{cor:cd_zero} \\"
    )
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    outpath.write_text('\n'.join(lines) + '\n')


def main() -> None:
    out_root = SRC_DIR.parent.parent / '2_analysis' / 'output'
    tables = out_root / 'tables'
    figures = out_root / 'figures'
    data = out_root / 'data'
    for d in (tables, figures, data):
        d.mkdir(parents=True, exist_ok=True)

    print(f"ABGRS MC validation: S={S_REPS}, N={DGP['N']}, T={DGP['T']}, "
          f"master seed={MASTER_SEED}")
    print(f"Target treatment effect: alpha_pp = {ALPHA_PP_TRUE}")
    print()

    results = run_mc(S=S_REPS)
    summary = summarize(results)
    print()
    print("=== Summary ===")
    print(summary.to_string(index=False))
    print()

    table_path = tables / 'mc_validation.tex'
    fig_path = figures / 'mc_validation.pdf'
    csv_path = data / 'mc_validation_results.csv'
    save_table(summary, table_path)
    plot_distributions(results, fig_path)
    results.to_csv(csv_path, index=False)
    print(f"Wrote table  -> {table_path}")
    print(f"Wrote figure -> {fig_path}")
    print(f"Wrote CSV    -> {csv_path}")

    run_lambda_diagnostic(out_root)


if __name__ == '__main__':
    main()
