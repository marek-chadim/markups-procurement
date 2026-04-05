"""ACF production function estimation on the Orbis all-industry panel.

Runs Cobb-Douglas ACF estimation by industry on the Orbis panel to get
proper markup-based procurement premiums. Covers top 5 procurement
industries (NACE 41, 43, 46, 71, 62), pooled construction (41+42+43),
and top-5 pooled.

Author: Marek Chadim (Yale, Tobin Center)
"""

import sys
import os
import gc
import time
import warnings
import numpy as np
import pandas as pd

from pathlib import Path

from acf_estimator import (
    ACFEstimator, Formulation, Optimization, CWDLExtensions,
    options as acf_options
)

warnings.filterwarnings('ignore')
acf_options.verbose = False

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
DATA_PATH = str(INPUT_DIR / 'orbis_panel.dta')
OUT_PATH = str(OUTPUT_DIR / 'tables' / 'orbis_acf_results.csv')

COLS = ['ico', 'year', 'go', 'k', 'cogs', 'pp_dummy', 'nace2', 'mktshare']


def compute_premium(res_data):
    """log(mu) ~ pp_dummy + year FE + k + cogs. Returns (coef, HC1 SE, R2, N)."""
    mu = res_data[res_data['markup'] > 0].copy()
    mu['lmu'] = np.log(mu['markup'])
    mu = mu.dropna(subset=['lmu', 'pp_dummy', 'k', 'cogs'])
    if mu['pp_dummy'].nunique() < 2:
        return np.nan, np.nan, np.nan, len(mu)
    X = pd.DataFrame({
        'const': 1.0, 'pp': mu['pp_dummy'].values,
        'k': mu['k'].values, 'cogs': mu['cogs'].values
    }, index=mu.index)
    years = mu['year'].astype(str)
    for val in sorted(years.unique())[1:]:
        X[f'yr_{val}'] = (years == val).astype(float)
    y = mu['lmu'].values
    Xm = X.values.astype(np.float64)
    bhat = np.linalg.lstsq(Xm, y, rcond=None)[0]
    resid = y - Xm @ bhat
    N = len(y)
    k = Xm.shape[1]
    XtX_inv = np.linalg.inv(Xm.T @ Xm)
    # Efficient HC1: X'diag(e^2)X = (X*e)' (X*e), avoids N x N matrix
    Xe = Xm * resid[:, None]
    meat = Xe.T @ Xe
    V = XtX_inv @ meat @ XtX_inv * N / (N - k)
    return bhat[1], np.sqrt(V[1, 1]), 1 - np.sum(resid**2) / np.sum((y - y.mean())**2), N


def run_one(df, label, nace_str, use_nace_fe=False):
    """Estimate ACF CD and compute procurement premium for one group."""
    form = Formulation(
        spec='cd', pp_in_markov=True, pp_interactions=True,
        year_fe=True, nace2_fe=use_nace_fe,
        first_stage_controls=['mktshare'],
    )
    t0 = time.time()
    est = ACFEstimator(
        data=df, formulation=form,
        optimization=Optimization(method='nm+bfgs'),
        extensions=CWDLExtensions(survival_correction=False),
        n_starts=3,
    )
    res = est.solve()
    elapsed = time.time() - t0

    bd = dict(zip(res.beta_names, res.betas))
    sd = dict(zip(res.beta_names, res.se))
    mu = res.data['markup'][(res.data['markup'] > 0) & np.isfinite(res.data['markup'])]
    prem, prem_se, r2, n_reg = compute_premium(res.data)

    row = {
        'industry': label, 'nace2': nace_str,
        'N_obs': res.n_obs, 'N_firms': res.n_clusters,
        'b_k': bd.get('k', np.nan), 'se_b_k': sd.get('k', np.nan),
        'b_cogs': bd.get('cogs', np.nan), 'se_b_cogs': sd.get('cogs', np.nan),
        'gmm_obj': float(res.gmm_criterion),
        'first_stage_r2': float(res.first_stage_r2),
        'markup_mean': float(mu.mean()), 'markup_sd': float(mu.std()),
        'markup_p10': float(mu.quantile(0.1)), 'markup_p50': float(mu.quantile(0.5)),
        'markup_p90': float(mu.quantile(0.9)),
        'frac_below_1': float((mu < 1).mean()),
        'pp_premium': prem, 'pp_premium_se': prem_se,
        'pp_premium_r2': r2, 'pp_premium_N': n_reg,
        'elapsed_s': elapsed,
    }
    print(f'  {label}: N={res.n_obs}, b_k={row["b_k"]:.4f}({row["se_b_k"]:.4f}), '
          f'b_cogs={row["b_cogs"]:.4f}({row["se_b_cogs"]:.4f}), '
          f'mean_mu={row["markup_mean"]:.4f}, '
          f'premium={prem:.4f}({prem_se:.4f}), {elapsed:.0f}s')
    return row


def main():
    print("Loading Orbis panel...")
    # Use parquet cache if available (much faster than .dta for 1.3M rows)
    PARQUET_CACHE = str(INPUT_DIR / 'orbis_subset.parquet')
    NACES_NEEDED = {41, 42, 43, 46, 71, 62}
    if os.path.exists(PARQUET_CACHE):
        df_all = pd.read_parquet(PARQUET_CACHE)
        df_all = df_all.rename(columns={'ico': 'id'})
        print(f"Loaded from parquet cache: {len(df_all)} obs")
    else:
        df_all = pd.read_stata(DATA_PATH, columns=COLS)
        df_all = df_all[df_all['nace2'].isin(NACES_NEEDED)]
        df_all.to_parquet(PARQUET_CACHE, index=False)
        df_all = df_all.rename(columns={'ico': 'id'})
        print(f"Loaded from Stata, cached: {len(df_all)} obs")

    results = []

    # Individual industries
    for nace, label in [
        (41, 'Buildings'),
        (43, 'Specialized construction'),
        (46, 'Wholesale trade'),
        (71, 'Architecture & engineering'),
        (62, 'Computer programming'),
    ]:
        print(f'\nNACE {nace}...')
        row = run_one(df_all[df_all['nace2'] == nace].copy(), label, str(nace))
        results.append(row)
        gc.collect()

    # Pooled construction
    print('\nConstruction pooled...')
    row = run_one(
        df_all[df_all['nace2'].isin([41, 42, 43])].copy(),
        'Construction (41+42+43)', '41+42+43', use_nace_fe=True)
    results.append(row)
    gc.collect()

    # Top-5 PP pooled
    print('\nTop-5 PP pooled...')
    row = run_one(
        df_all[df_all['nace2'].isin([41, 43, 46, 71, 62])].copy(),
        'Top-5 PP industries pooled', 'pooled', use_nace_fe=True)
    results.append(row)
    gc.collect()

    # Save
    out = pd.DataFrame(results)
    out.to_csv(OUT_PATH, index=False, float_format='%.6f')
    print(f'\nSaved to {OUT_PATH}')

    # Summary
    print(f'\n{"="*90}')
    print(f'  SUMMARY: Orbis ACF Estimation Results')
    print(f'{"="*90}')
    fmt = '{:<35s} {:>7s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}'
    print(fmt.format('Industry', 'N', 'b_k', 'b_cogs', 'Mean mu', 'Premium', 'SE'))
    print(f'  {"-"*82}')
    for r in results:
        print(f'  {r["industry"]:<33s} {r["N_obs"]:>7d} {r["b_k"]:>8.4f} '
              f'{r["b_cogs"]:>8.4f} {r["markup_mean"]:>8.4f} '
              f'{r["pp_premium"]:>8.4f} {r["pp_premium_se"]:>8.4f}')


if __name__ == '__main__':
    main()
