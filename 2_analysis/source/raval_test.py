#!/usr/bin/env python3
"""Raval (2023) overidentification test on Czech construction data.

Estimates markups using four flexible inputs and tests whether they agree:
  - Labor (W): wages
  - Materials (COGS): cost of goods sold
  - Materials net (II-W): intermediate inputs minus wages
  - Composite (W+COGS): labor + materials bundled

Under Hicks-neutral productivity, all should give the same markup.
Raval (2023) rejects this across 7 datasets.
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from pathlib import Path

from acf_estimator import (
    Formulation, ACFEstimator, estimate_by_industry
)

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

DATA_PATH = str(INPUT_DIR / 'data.dta')
OUT = str(OUTPUT_DIR) + '/'

# ── Load data ──────────────────────────────────────────────────────────
print('Loading data...')
df = pd.read_stata(DATA_PATH)

# Construct additional variables
df['mat_net'] = np.log(np.exp(df['ii']) - np.exp(df['w']))  # II - W
df.loc[~np.isfinite(df['mat_net']), 'mat_net'] = np.nan
df['composite'] = np.log(np.exp(df['w']) + np.exp(df['cogs']))  # W + COGS

for v in ['w', 'cogs', 'mat_net', 'composite']:
    n = df[v].notna().sum()
    print(f'  {v}: {n:,} non-missing')

print(f'Data: {len(df):,} obs, {df["id"].nunique():,} firms')

# ── Estimate markups with each variable input (2-input CD) ─────────────
INPUTS = {
    'w': 'Labor (W)',
    'cogs': 'Materials (COGS)',
    'composite': 'Composite (W+COGS)',
}

all_mu = {}
for vi, label in INPUTS.items():
    print(f'\n{"=" * 60}')
    print(f'  {label}: variable_input = {vi}')
    print(f'{"=" * 60}')

    results, coefficients, markups = estimate_by_industry(
        df, specs=('cd',),
        formulation_kwargs={'variable_input': vi, 'pp_in_markov': True},
        n_starts=3,
    )
    markups = markups.rename(columns={'markup': f'mu_{vi}'})
    all_mu[vi] = markups[['id', 'year', 'nace2', 'pp_dummy', f'mu_{vi}']].copy()

    med = markups.loc[markups[f'mu_{vi}'] > 0, f'mu_{vi}'].median()
    print(f'  → median markup = {med:.3f}')

# ── Merge ──────────────────────────────────────────────────────────────
merged = all_mu['cogs'].copy()
for vi in ['w', 'composite']:
    merged = merged.merge(all_mu[vi][['id', 'year', f'mu_{vi}']],
                          on=['id', 'year'], how='inner')

for vi in INPUTS:
    merged = merged[merged[f'mu_{vi}'] > 0]
    merged[f'lmu_{vi}'] = np.log(merged[f'mu_{vi}'])
merged = merged.dropna()

print(f'\nMerged: {len(merged):,} obs with all 4 markups positive')
merged.to_csv(OUT + 'data/raval_test_markups.csv', index=False)

# ── Diagnostic 1: Dispersion ──────────────────────────────────────────
print(f'\n{"=" * 60}')
print('  DISPERSION (Raval Table II)')
print(f'{"=" * 60}')

disp_rows = []
for vi, label in INPUTS.items():
    mu = merged[f'mu_{vi}']
    p10, p25, p50, p75, p90 = mu.quantile([.10, .25, .50, .75, .90])
    row = {'Input': label, 'Mean': mu.mean(), 'Median': p50,
           '90/50': p90/p50, '90/10': p90/p10, 'Frac<1': (mu<1).mean()}
    disp_rows.append(row)
    print(f'  {label:22s}: mean={mu.mean():.3f} med={p50:.3f} '
          f'90/50={p90/p50:.2f} frac<1={row["Frac<1"]:.3f}')
disp_df = pd.DataFrame(disp_rows)

# ── Diagnostic 2: Cross-Input Correlation ─────────────────────────────
print(f'\n{"=" * 60}')
print('  CROSS-INPUT CORRELATION (Raval Table III)')
print(f'{"=" * 60}')

correl_rows = []
pairs = [('cogs', 'w'), ('composite', 'w'), ('composite', 'cogs')]

for vi_y, vi_x in pairs:
    reg = merged[['id', 'year', 'nace2', f'lmu_{vi_y}', f'lmu_{vi_x}']].copy()
    reg = reg.replace([np.inf, -np.inf], np.nan).dropna()
    reg['yc'] = pd.Categorical(reg['year'])
    reg['nc'] = pd.Categorical(reg['nace2'])
    try:
        m = smf.ols(f'lmu_{vi_y} ~ lmu_{vi_x} + C(yc) + C(nc)', data=reg).fit(
            cov_type='cluster', cov_kwds={'groups': reg['id']})
        b, se = m.params[f'lmu_{vi_x}'], m.bse[f'lmu_{vi_x}']
        correl_rows.append({'Dep': INPUTS[vi_y], 'Indep': INPUTS[vi_x],
                            'beta': b, 'se': se, 'N': int(m.nobs)})
        print(f'  {INPUTS[vi_y]:22s} on {INPUTS[vi_x]:22s}: '
              f'beta={b:.4f} (SE={se:.4f})')
    except Exception as e:
        print(f'  FAILED: {e}')
correl_df = pd.DataFrame(correl_rows)

# ── Diagnostic 3: Time Trends ─────────────────────────────────────────
print(f'\n{"=" * 60}')
print('  TIME TRENDS (Raval Figure 2)')
print(f'{"=" * 60}')

fig, ax = plt.subplots(figsize=(10, 6))
colors = {'w': '#b2182b', 'cogs': '#2166ac', 'composite': '#4daf4a'}
styles = {'w': '-', 'cogs': '--', 'composite': '-.'}
years = sorted(merged['year'].unique())

for vi, label in INPUTS.items():
    reg = merged[['year', 'nace2', f'lmu_{vi}']].replace([np.inf, -np.inf], np.nan).dropna()
    reg['yc'] = pd.Categorical(reg['year'])
    reg['nc'] = pd.Categorical(reg['nace2'])
    m = smf.ols(f'lmu_{vi} ~ C(yc) + C(nc)', data=reg).fit()
    fe = [0.0] + [m.params.get(f'C(yc)[T.{yr}]', np.nan) * 100 for yr in years[1:]]
    ax.plot(years, fe, color=colors[vi], linestyle=styles[vi], linewidth=2, label=label)
    print(f'  {label}: {years[0]}-{years[-1]} change = {fe[-1]:.1f}%')

ax.axhline(0, color='gray', linestyle=':', linewidth=0.5)
ax.set_xlabel('Year')
ax.set_ylabel('Percent Change')
ax.legend()
plt.tight_layout()
plt.savefig(OUT + 'figures/raval_time_trends.pdf', bbox_inches='tight')
print(f'  Saved: raval_time_trends.pdf')

# ── LaTeX table ───────────────────────────────────────────────────────
with open(OUT + 'tables/raval_overid.tex', 'w') as f:
    f.write(r'\begin{table}[htbp]\centering' + '\n')
    f.write(r'\caption{Overidentification Test: Variable Input Choice}\label{tab:raval}' + '\n')
    f.write(r'\begin{tabular}{lccc}' + '\n')
    f.write(r'\hline\hline' + '\n')
    f.write(r' & Labor & COGS & Composite \\' + '\n')
    f.write(r'\hline' + '\n')
    f.write(r'\multicolumn{4}{l}{\emph{Panel A: Markup Dispersion}} \\[3pt]' + '\n')
    order = [INPUTS[v] for v in ['w', 'cogs', 'composite']]
    for stat in ['Mean', 'Median', '90/50', '90/10', 'Frac<1']:
        vals = [disp_df.loc[disp_df['Input']==lab, stat].values[0] for lab in order]
        fmt = '{:.3f}' if stat == 'Frac<1' else '{:.2f}'
        f.write(f'{stat} & ' + ' & '.join(fmt.format(v) for v in vals) + r' \\' + '\n')
    f.write(r'\\' + '\n')
    f.write(r'\multicolumn{4}{l}{\emph{Panel B: Cross-Input Correlation ($\beta$, null = 1)}} \\[3pt]' + '\n')
    for _, r in correl_df.iterrows():
        f.write(f'{r["Dep"]} on {r["Indep"]} & '
                f'\\multicolumn{{2}}{{c}}{{{r["beta"]:.3f} ({r["se"]:.3f})}} & '
                f'\\multicolumn{{2}}{{c}}{{{r["N"]:,d}}}' + r' \\' + '\n')
    f.write(r'\hline\hline' + '\n')
    f.write(r'\end{tabular}' + '\n')
    f.write(r'\begin{minipage}{0.9\textwidth}\footnotesize' + '\n')
    f.write(r'\emph{Notes:} Cobb-Douglas production functions estimated by NACE 2-digit. '
            r'Labor = wages; COGS = cost of goods sold; '
            r'Composite = wages + COGS. Panel B: '
            r'$\log\mu^Y_{it} = \alpha + \beta\log\mu^X_{it} + \gamma_t + \delta_n + \varepsilon_{it}$; '
            r'firm-clustered SEs. Under Hicks neutrality, $\beta = 1$.' + '\n')
    f.write(r'\end{minipage}' + '\n')
    f.write(r'\end{table}' + '\n')
print(f'  Saved: raval_overid.tex')

# ── Summary ───────────────────────────────────────────────────────────
print(f'\n{"=" * 60}')
print(f'  SUMMARY: N = {len(merged):,}')
for _, r in disp_df.iterrows():
    print(f'    {r["Input"]:22s}: 90/50 = {r["90/50"]:.2f}')
print(f'  Correlations (null: beta=1):')
for _, r in correl_df.iterrows():
    sign = '+' if r['beta'] > 0 else ''
    print(f'    {r["Dep"]:22s} on {r["Indep"]:22s}: {sign}{r["beta"]:.3f}')
