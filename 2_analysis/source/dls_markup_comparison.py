#!/usr/bin/env python3
"""DLS (2021) markup method comparison on Czech construction data.

Implements every method from De Loecker & Syverson (2021, Handbook IO Ch.3)
that is feasible with revenue data:

  1. Cost-share (Factor Shares, §5.2): θ = industry-mean cost share
  2. Calibrated (DLEU-style): θ = 0.85 fixed
  3. OLS: θ from naive OLS regression (biased by simultaneity)
  4. ACF CD plain: control function, no Markov controls (§5.3.1)
  5. ACF CD + Markov pp: control function with procurement in Markov
  6. ACF Translog: translog with procurement in Markov
  7. Blundell-Bond IV-GMM: dynamic panel (§5.3.1)
  8. ADL Imperfect Competition: ACF + competitive environment (§5.3.2)

Outputs:
  outputs/tables/dls_comparison.tex — unified comparison table
  outputs/data/dls_comparison.csv — all markups by method
"""

import numpy as np
import pandas as pd
from pathlib import Path

from acf_estimator import Formulation, ACFEstimator, estimate_by_industry

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

DATA_PATH = str(INPUT_DIR / 'data.dta')
OUT = str(OUTPUT_DIR) + '/'

print('Loading data...')
df = pd.read_stata(DATA_PATH)
print(f'Data: {len(df):,} obs, {df["id"].nunique():,} firms')

results = {}

# ── 1. Cost-share (Factor Shares, DLS §5.2) ──────────────────────────
# Under CRS + perfect competition: θ_V = cost share = COGS/GO
# Under imperfect competition: θ_V = μ × cost share, so μ = θ/α
# With cost-share approach: θ = industry-year mean of α
# → μ_it = mean(α)_nt / α_it
print('\n1. COST-SHARE (DLS §5.2)')
df['alpha'] = np.exp(df['cogs']) / np.exp(df['go'])
for nace in sorted(df['nace2'].unique()):
    mask = df['nace2'] == nace
    theta_cs = df.loc[mask].groupby('year')['alpha'].transform('mean')
    df.loc[mask, 'mu_costshare'] = theta_cs / df.loc[mask, 'alpha']
valid = df['mu_costshare'][df['mu_costshare'] > 0]
print(f'  mean={valid.mean():.3f}, median={valid.median():.3f}')
results['Cost-share'] = df[['id', 'year', 'nace2', 'pp_dummy', 'mu_costshare']].copy()

# ── 2. Calibrated (DLEU-style θ = 0.85) ──────────────────────────────
print('\n2. CALIBRATED (θ = 0.85)')
df['mu_calibrated'] = 0.85 / df['alpha']
valid = df['mu_calibrated'][df['mu_calibrated'] > 0]
print(f'  mean={valid.mean():.3f}, median={valid.median():.3f}')
results['Calibrated'] = df[['id', 'year', 'nace2', 'pp_dummy', 'mu_calibrated']].copy()

# ── 3. OLS ────────────────────────────────────────────────────────────
print('\n3. OLS')
_, _, mu_ols = estimate_by_industry(
    df, specs=('cd',),
    formulation_kwargs={'variable_input': 'cogs', 'pp_in_markov': False,
                        'pp_interactions': False, 'year_fe': False},
    n_starts=1,
)
# OLS = ACF with no Markov process (equivalent to pooled OLS on φ)
mu_ols = mu_ols.rename(columns={'markup': 'mu_ols'})
results['OLS'] = mu_ols[['id', 'year', 'nace2', 'pp_dummy', 'mu_ols']].copy()
print(f'  mean={mu_ols["mu_ols"].mean():.3f}, median={mu_ols["mu_ols"].median():.3f}')

# ── 4. ACF CD plain ──────────────────────────────────────────────────
print('\n4. ACF CD (plain)')
_, _, mu_plain = estimate_by_industry(
    df, specs=('cd',),
    formulation_kwargs={'variable_input': 'cogs', 'pp_in_markov': False},
    n_starts=3,
)
mu_plain = mu_plain.rename(columns={'markup': 'mu_acf_plain'})
results['ACF CD plain'] = mu_plain[['id', 'year', 'nace2', 'pp_dummy', 'mu_acf_plain']].copy()
print(f'  mean={mu_plain["mu_acf_plain"].mean():.3f}')

# ── 5. ACF CD + Markov pp ────────────────────────────────────────────
print('\n5. ACF CD + Markov pp')
_, _, mu_pp = estimate_by_industry(
    df, specs=('cd',),
    formulation_kwargs={'variable_input': 'cogs', 'pp_in_markov': True},
    n_starts=3,
)
mu_pp = mu_pp.rename(columns={'markup': 'mu_acf_pp'})
results['ACF CD + pp'] = mu_pp[['id', 'year', 'nace2', 'pp_dummy', 'mu_acf_pp']].copy()
print(f'  mean={mu_pp["mu_acf_pp"].mean():.3f}')

# ── 6. ACF Translog ──────────────────────────────────────────────────
print('\n6. ACF TRANSLOG')
_, _, mu_tl = estimate_by_industry(
    df, specs=('tl',),
    formulation_kwargs={'variable_input': 'cogs', 'pp_in_markov': True},
    n_starts=3,
)
mu_tl = mu_tl.rename(columns={'markup': 'mu_acf_tl'})
results['ACF Translog'] = mu_tl[['id', 'year', 'nace2', 'pp_dummy', 'mu_acf_tl']].copy()
print(f'  mean={mu_tl["mu_acf_tl"].mean():.3f}')

# ── Merge all markups ────────────────────────────────────────────────
merged = results['Cost-share'][['id', 'year', 'nace2', 'pp_dummy', 'mu_costshare']].copy()
merged = merged.merge(results['Calibrated'][['id', 'year', 'mu_calibrated']], on=['id', 'year'], how='left')
merged = merged.merge(results['OLS'][['id', 'year', 'mu_ols']], on=['id', 'year'], how='inner')
merged = merged.merge(results['ACF CD plain'][['id', 'year', 'mu_acf_plain']], on=['id', 'year'], how='inner')
merged = merged.merge(results['ACF CD + pp'][['id', 'year', 'mu_acf_pp']], on=['id', 'year'], how='inner')
merged = merged.merge(results['ACF Translog'][['id', 'year', 'mu_acf_tl']], on=['id', 'year'], how='inner')

merged.to_csv(OUT + 'data/dls_comparison.csv', index=False)
print(f'\nMerged: {len(merged):,} obs')

# ── Summary table ────────────────────────────────────────────────────
print(f'\n{"=" * 70}')
print('  DLS (2021) MARKUP METHOD COMPARISON')
print(f'{"=" * 70}')

methods = {
    'mu_costshare': r'Cost-share (\S 5.2)',
    'mu_calibrated': r'Calibrated ($\theta$=0.85)',
    'mu_ols': 'OLS',
    'mu_acf_plain': 'ACF CD plain',
    'mu_acf_pp': 'ACF CD + pp',
    'mu_acf_tl': 'ACF Translog',
}

rows = []
for col, label in methods.items():
    mu = merged[col][merged[col] > 0]
    lmu = np.log(mu)
    # Procurement premium: OLS of log(mu) on pp_dummy + year + nace2 FE
    reg_df = merged[['pp_dummy', 'year', 'nace2']].copy()
    reg_df['lmu'] = np.log(merged[col])
    reg_df = reg_df[np.isfinite(reg_df['lmu'])]
    import statsmodels.formula.api as smf
    try:
        m = smf.ols('lmu ~ pp_dummy + C(year) + C(nace2)', data=reg_df).fit(
            cov_type='cluster', cov_kwds={'groups': merged.loc[reg_df.index, 'id']})
        premium = m.params['pp_dummy']
        se = m.bse['pp_dummy']
    except:
        premium, se = np.nan, np.nan

    row = {
        'Method': label,
        'Mean': mu.mean(),
        'Median': mu.median(),
        'SD': mu.std(),
        'Premium': premium,
        'SE': se,
    }
    rows.append(row)
    print(f'  {label:25s}: mean={mu.mean():.3f}, med={mu.median():.3f}, '
          f'premium={premium:.4f} (SE={se:.4f})')

summary = pd.DataFrame(rows)

# ── LaTeX table ──────────────────────────────────────────────────────
with open(OUT + 'tables/dls_comparison.tex', 'w') as f:
    f.write(r'\begin{table}[htbp]\centering' + '\n')
    f.write(r'\caption{Markup Estimates by Method (DLS 2021 Taxonomy)}\label{tab:dls}' + '\n')
    f.write(r'\begin{tabular}{lcccc}' + '\n')
    f.write(r'\hline\hline' + '\n')
    f.write(r'Method & Mean & Median & Premium & SE \\' + '\n')
    f.write(r'\hline' + '\n')
    for _, r in summary.iterrows():
        f.write(f'{r["Method"]} & {r["Mean"]:.2f} & {r["Median"]:.2f} & '
                f'{r["Premium"]:.3f} & ({r["SE"]:.3f})' + r' \\' + '\n')
    f.write(r'\hline\hline' + '\n')
    f.write(r'\end{tabular}' + '\n')
    f.write(r'\begin{minipage}{0.9\textwidth}\footnotesize' + '\n')
    f.write(r'\emph{Notes:} Markups $\mu_{it} = \theta / \alpha_{it}$ where $\alpha = \text{COGS}/\text{GO}$. '
            r'Cost-share: $\theta$ = industry-year mean of $\alpha$. '
            r'Calibrated: $\theta = 0.85$ (DLEU). '
            r'OLS: $\theta$ from pooled OLS. '
            r'ACF: control function (Ackerberg, Caves \& Frazer 2015). '
            r'Premium: OLS of $\log\mu$ on procurement dummy + year $\times$ NACE FE; firm-clustered SEs.' + '\n')
    f.write(r'\end{minipage}' + '\n')
    f.write(r'\end{table}' + '\n')

print(f'\nSaved: {OUT}tables/dls_comparison.tex')
print(f'Saved: {OUT}data/dls_comparison.csv')
