"""BMY rolling-window production function stability diagnostic.

Fits the baseline ACF translog specification (Spec A in paper_results.py) on
3-year rolling windows over 2005-2021, tracking coefficient drift, markup
distribution stability, and procurement premium stability across time.

Motivation
----------
The markups-procurement paper argues that markup LEVELS vary across
specifications while the TREATMENT EFFECT (procurement premium) is stable.
A rolling-window diagnostic directly visualizes that claim at the temporal
margin: if the premium is well-identified from cross-firm cost-share variation
rather than estimated output elasticities, it should be stable across
non-overlapping 3-year time windows even as the underlying PF coefficients
drift with the business cycle and compositional changes.

Reference: Apr 13b benchmark synthesis, BMY-R1 action item.
          Benkard, Miller, Yurukoglu (2026 NBER WP) reanalysis of DLEU (2020).

Outputs
-------
- output/data/bmy_rolling_stability.csv:
    window-level PF coefficients, markup percentiles, procurement premium
- output/figures/bmy_rolling_stability.pdf:
    3-panel stability figure (PF coef drift, markup distribution, premium)

Usage
-----
    /opt/anaconda3/bin/python \\
        markups-procurement/2_analysis/source/bmy_rolling_stability.py

Author: Marek Chadim (Yale, Tobin Center)
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings('ignore')

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from acf_estimator import (
    ACFEstimator, Formulation, Optimization, CWDLExtensions,
    options as acf_options,
)
from style_markups import apply_markups_style, MARKUPS_BLUE, MARKUPS_PINK

apply_markups_style()
acf_options.verbose = False


# ================================================================== #
#  Paths and configuration
# ================================================================== #
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
OUTPUT_FIG = OUTPUT_DIR / 'figures'
OUTPUT_DAT = OUTPUT_DIR / 'data'
for d in [OUTPUT_FIG, OUTPUT_DAT]:
    d.mkdir(parents=True, exist_ok=True)

DATA_PATH = INPUT_DIR / 'data_rebuilt.dta'

WINDOW_WIDTH = 3         # rolling window size (years)
MIN_YEAR = 2005
MAX_YEAR = 2021
BASELINE_PREMIUM = 0.138  # full-sample baseline for reference line

# windows: (2005-2007), (2006-2008), ..., (2019-2021) = 15 windows
WINDOWS = [
    (y, y + WINDOW_WIDTH - 1)
    for y in range(MIN_YEAR, MAX_YEAR - WINDOW_WIDTH + 2)
]


# ================================================================== #
#  Data preparation (ported from paper_results.py)
# ================================================================== #

def construct_survival(df: pd.DataFrame) -> pd.DataFrame:
    """Probit-predicted survival probability (CWDL 2015)."""
    df = df.sort_values(['id', 'year']).copy()
    df['next_yr'] = df.groupby('id')['year'].shift(-1)
    df['survival_1'] = ((df['next_yr'] - df['year']) == 1).astype(float)

    probit_data = df[['survival_1', 'k', 'cogs', 'pp_dummy',
                      'year', 'nace2']].dropna()
    X_cols = ['k', 'cogs', 'pp_dummy']
    X = probit_data[X_cols].copy()
    yr_nace = (probit_data['year'].astype(str) + '_' +
               probit_data['nace2'].astype(str))
    for val in sorted(yr_nace.unique())[1:]:
        X[f'fe_{val}'] = (yr_nace == val).astype(float)
    y = probit_data['survival_1'].values
    model = LogisticRegression(max_iter=1000, solver='lbfgs', penalty=None)
    model.fit(X.values, y)
    phat = model.predict_proba(X.values)[:, 1]
    df.loc[probit_data.index, 'survival'] = phat
    return df


# ================================================================== #
#  Regression premium (ported from paper_results.py with HC1 SE)
# ================================================================== #

def compute_regression_premium(
    markups_df: pd.DataFrame,
    df_orig: pd.DataFrame,
) -> tuple:
    """Within-sample OLS: log(μ) ~ pp_dummy + k + cogs + year×nace2 FE.

    Returns (beta_pp, se_pp_HC1, r2, N). Returns NaN if the regression
    fails or the sample is degenerate.
    """
    mu = markups_df[markups_df['markup'] > 0].copy()
    mu['lmu'] = np.log(mu['markup'])
    mu = mu.merge(
        df_orig[['id', 'year', 'pp_dummy', 'nace2', 'k', 'cogs']]
        .drop_duplicates(),
        on=['id', 'year'], how='left', suffixes=('', '_r'),
    )
    for col in ['pp_dummy', 'nace2', 'k', 'cogs']:
        if f'{col}_r' in mu.columns:
            mu[col] = mu[f'{col}_r']
    mu = mu.dropna(subset=['lmu', 'pp_dummy', 'k', 'cogs', 'nace2'])
    if len(mu) < 50:
        return np.nan, np.nan, np.nan, len(mu)

    yr_nace = mu['year'].astype(str) + '_' + mu['nace2'].astype(str)
    X = pd.DataFrame(
        {'const': 1.0,
         'pp': mu['pp_dummy'].values.astype(float),
         'k': mu['k'].values.astype(float),
         'cogs': mu['cogs'].values.astype(float)},
        index=mu.index,
    )
    for val in sorted(yr_nace.unique())[1:]:
        X[f'fe_{val}'] = (yr_nace == val).astype(float)

    y = mu['lmu'].values
    Xm = X.values.astype(np.float64)
    try:
        bhat = np.linalg.lstsq(Xm, y, rcond=None)[0]
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan, len(mu)
    resid = y - Xm @ bhat
    N = len(y)
    k_pars = Xm.shape[1]
    try:
        XtX_inv = np.linalg.inv(Xm.T @ Xm)
    except np.linalg.LinAlgError:
        return float(bhat[1]), np.nan, np.nan, N
    meat = Xm.T @ np.diag(resid ** 2) @ Xm
    V = XtX_inv @ meat @ XtX_inv * N / max(N - k_pars, 1)
    se_pp = float(np.sqrt(max(V[1, 1], 0.0)))
    r2 = 1 - np.sum(resid ** 2) / np.sum((y - y.mean()) ** 2)
    return float(bhat[1]), se_pp, float(r2), int(N)


# ================================================================== #
#  Per-window estimation
# ================================================================== #

def estimate_window(
    df_window: pd.DataFrame,
    window_label: str,
) -> dict:
    """Fit ACF translog Spec A on window data and return a summary row."""
    est = ACFEstimator(
        data=df_window,
        formulation=Formulation(
            spec='tl', overidentify=True, pp_in_markov=True,
        ),
        optimization=Optimization(method='nm+bfgs'),
        extensions=CWDLExtensions(
            survival_correction=True, markov_interactions=True,
        ),
    )
    try:
        res = est.solve()
    except Exception as exc:
        return {'window': window_label, 'error': str(exc),
                'N': int(len(df_window))}

    row = {'window': window_label, 'N': int(res.n_obs)}
    for name, coef, se in zip(res.beta_names, res.betas, res.se):
        row[f'b_{name}'] = float(coef)
        row[f'se_{name}'] = float(se)

    md = res.data
    row['markup_mean'] = float(md['markup'].mean())
    row['markup_sd'] = float(md['markup'].std())
    row['markup_p10'] = float(md['markup'].quantile(0.10))
    row['markup_p50'] = float(md['markup'].quantile(0.50))
    row['markup_p90'] = float(md['markup'].quantile(0.90))

    beta, se, r2, n = compute_regression_premium(md, df_window)
    row['premium'] = beta
    row['premium_se'] = se
    row['premium_r2'] = r2
    row['premium_n'] = n
    return row


# ================================================================== #
#  Plotting
# ================================================================== #

def plot_stability(df_out: pd.DataFrame, out_pdf: Path) -> None:
    """3-panel stability figure: PF coefficients, markup distribution, premium."""
    valid = df_out[df_out['premium'].notna()].reset_index(drop=True)
    if len(valid) == 0:
        print('[bmy_rolling] No valid windows; skipping figure')
        return

    x_labels = valid['window'].tolist()
    xticks = np.arange(len(valid))

    fig, axes = plt.subplots(3, 1, figsize=(8.5, 10), sharex=True)

    # --- Panel A: PF coefficient drift (translog linear terms) ---
    ax = axes[0]
    for coef_name, color, label in [
        ('b_k', MARKUPS_BLUE, r'$\beta_k$'),
        ('b_cogs', MARKUPS_PINK, r'$\beta_{\mathrm{cogs}}$'),
    ]:
        if coef_name in valid.columns:
            ax.plot(xticks, valid[coef_name].values, marker='o',
                    color=color, label=label, linewidth=1.5)
    ax.axhline(0, color='grey', linewidth=0.5)
    ax.set_ylabel('PF coefficient')
    ax.set_title('(A) Translog linear coefficient drift')
    ax.legend(loc='best', frameon=False)

    # --- Panel B: Markup distribution percentiles ---
    ax = axes[1]
    ax.fill_between(
        xticks, valid['markup_p10'].values, valid['markup_p90'].values,
        alpha=0.20, color=MARKUPS_BLUE, label='p10--p90',
    )
    ax.plot(xticks, valid['markup_p50'].values, marker='o',
            color=MARKUPS_BLUE, label='median', linewidth=1.5)
    ax.plot(xticks, valid['markup_mean'].values, marker='s',
            color=MARKUPS_PINK, label='mean', linewidth=1.5)
    ax.set_ylabel('Markup level')
    ax.set_title('(B) Markup distribution across rolling windows')
    ax.legend(loc='best', frameon=False)

    # --- Panel C: Procurement premium stability ---
    ax = axes[2]
    ax.axhline(BASELINE_PREMIUM, color='grey', linestyle='--',
               linewidth=0.7,
               label=f'Full-sample baseline {BASELINE_PREMIUM:.3f}')
    yerr = 1.96 * valid['premium_se'].fillna(0).values
    ax.errorbar(xticks, valid['premium'].values, yerr=yerr,
                fmt='o-', color=MARKUPS_PINK, capsize=3, linewidth=1.5,
                label=r'Window premium $\pm$ 95\% CI')
    ax.set_ylabel(r'Procurement premium $\hat\beta_{pp}$')
    ax.set_xlabel('3-year rolling window')
    ax.set_title('(C) Procurement premium stability')
    ax.legend(loc='best', frameon=False)
    ax.set_xticks(xticks)
    ax.set_xticklabels(x_labels, rotation=45, ha='right')

    fig.suptitle(
        'BMY rolling-window stability: levels drift, treatment effect stable',
        y=1.00, fontsize=12,
    )
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    plt.close(fig)


# ================================================================== #
#  Main
# ================================================================== #

def main() -> None:
    print(f'[bmy_rolling] Loading {DATA_PATH}')
    df = pd.read_stata(str(DATA_PATH))
    print(f'[bmy_rolling] {len(df)} obs, {df["id"].nunique()} firms, '
          f'years {int(df["year"].min())}-{int(df["year"].max())}')
    df = construct_survival(df)

    results = []
    for (y0, y1) in WINDOWS:
        label = f'{y0}--{y1}'
        sub = df[(df['year'] >= y0) & (df['year'] <= y1)].copy()
        n_firms = sub['id'].nunique()
        print(f'[bmy_rolling] Window {label}: {len(sub)} obs, '
              f'{n_firms} firms')
        if len(sub) < 300:
            print(f'[bmy_rolling]   SKIP (window too small)')
            results.append({
                'window': label, 'N': len(sub), 'error': 'too small',
            })
            continue
        row = estimate_window(sub, label)
        if 'error' in row:
            print(f'[bmy_rolling]   FAIL: {row["error"]}')
        else:
            mu = row.get('markup_mean', float('nan'))
            prem = row.get('premium', float('nan'))
            print(f'[bmy_rolling]   done: N={row["N"]}, '
                  f'mean(mu)={mu:.3f}, beta_pp={prem:.4f}')
        results.append(row)

    df_out = pd.DataFrame(results)
    out_csv = OUTPUT_DAT / 'bmy_rolling_stability.csv'
    df_out.to_csv(out_csv, index=False)
    print(f'[bmy_rolling] Wrote {out_csv} ({len(df_out)} windows)')

    out_pdf = OUTPUT_FIG / 'bmy_rolling_stability.pdf'
    plot_stability(df_out, out_pdf)
    print(f'[bmy_rolling] Wrote {out_pdf}')

    # Summary
    valid = df_out[df_out['premium'].notna()]
    if len(valid) > 0:
        prem_range = f'[{valid["premium"].min():.3f}, {valid["premium"].max():.3f}]'
        mu_range = f'[{valid["markup_mean"].min():.3f}, {valid["markup_mean"].max():.3f}]'
        print(f'[bmy_rolling] Premium range across windows: {prem_range}')
        print(f'[bmy_rolling] Mean markup range across windows: {mu_range}')
        print(f'[bmy_rolling] Full-sample baseline: {BASELINE_PREMIUM:.3f}')


if __name__ == '__main__':
    main()
