#!/usr/bin/env python3
"""DLS Table 2 replication: markups across production technologies.

Reproduces De Loecker & Scott (2025, Table 2) format on Czech construction
data. Reports sales-weighted mean markups for each combination of:
  - Functional form: Cobb-Douglas vs Translog
  - Variable input: COGS (materials) vs W (labor)
  - Output measure: Gross output vs Value added

Standard errors on theta (output elasticity) via ACH (2012) analytical GMM
sandwich estimator with firm-level clustering.

Author: Marek Chadim (Yale, Tobin Center)
"""

import numpy as np
import pandas as pd
from pathlib import Path

from acf_estimator import (ACFEstimator, Formulation, Optimization,
                           CWDLExtensions, options)

options.verbose = False

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

# ── Load data ────────────────────────────────────────────────────────────
print('Loading data...')
df_base = pd.read_stata(str(INPUT_DIR / 'data.dta'))
df_base = df_base.dropna(subset=['go', 'k', 'cogs', 'w', 'pp_dummy',
                                  'year', 'nace2']).copy()
df_base['sales'] = np.exp(df_base['go'])

# Construct VA if missing: va = log(exp(go) - exp(cogs))
if 'va' not in df_base.columns or df_base['va'].isna().all():
    go_level = np.exp(df_base['go'])
    cogs_level = np.exp(df_base['cogs'])
    va_level = go_level - cogs_level
    df_base['va'] = np.where(va_level > 0, np.log(va_level), np.nan)
    print('  Constructed VA = log(GO - COGS)')

print(f'Data: {len(df_base):,} obs, {df_base["id"].nunique():,} firms')

# ── Specifications ───────────────────────────────────────────────────────
# (label, func_form, output_col, variable_input, output_type)
specs = [
    ('CD, Variable COGS',      'cd', 'go', 'cogs', 'Gross output'),
    ('CD, Variable W',         'cd', 'go', 'w',    'Gross output'),
    ('TL, Variable COGS',      'tl', 'go', 'cogs', 'Gross output'),
    ('TL, Variable W',         'tl', 'go', 'w',    'Gross output'),
    ('CD, Variable W (VA)',    'cd', 'va', 'w',    'Value added'),
    ('TL, Variable W (VA)',    'tl', 'va', 'w',    'Value added'),
]

# ── Estimation loop ──────────────────────────────────────────────────────
results = []
for label, func, out_col, vi, out_type in specs:
    print(f'\n--- {label} ---')
    df = df_base.copy()

    # For VA specs, replace go with va
    if out_col == 'va':
        df = df.dropna(subset=['va']).copy()
        df['go'] = df['va']  # ACFEstimator uses 'go' internally

    form = Formulation(spec=func, variable_input=vi, pp_in_markov=True,
                       pp_interactions=True, year_fe=True, nace2_fe=True)
    try:
        est = ACFEstimator(df, formulation=form,
                           optimization=Optimization(method='nm+bfgs'),
                           extensions=CWDLExtensions(survival_correction=True),
                           n_starts=3)
        res = est.solve()

        # Get markup and compute sales-weighted mean
        md = res.data[res.data['markup'] > 0].copy()
        md = md.merge(df_base[['id', 'year', 'sales']].drop_duplicates(
            ['id', 'year']), on=['id', 'year'], how='left')
        md['sw'] = md['sales'] / md['sales'].sum()
        mu_sw = (md['sw'] * md['markup']).sum()
        mu_mean = md['markup'].mean()
        mu_median = md['markup'].median()

        # Extract theta (output elasticity of variable input)
        # variable_input is mapped to 'cogs' internally
        vi_key = 'cogs'
        vi_idx = res.beta_names.index(vi_key)
        theta = res.betas[vi_idx]
        se_theta = res.se[vi_idx]

        # Markov rho
        rho = (res.markov_coefs[1]
               if res.markov_coefs is not None and len(res.markov_coefs) > 1
               else np.nan)

        results.append({
            'label': label, 'func_form': func, 'output': out_type,
            'variable_input': vi,
            'theta': theta, 'se_theta': se_theta,
            'mu_sw': mu_sw, 'mu_mean': mu_mean, 'mu_median': mu_median,
            'rho': rho, 'N': res.n_obs,
        })
        print(f'  theta={theta:.4f} (SE {se_theta:.4f}), '
              f'mu_sw={mu_sw:.3f}, mu_mean={mu_mean:.3f}')
    except Exception as e:
        print(f'  FAILED: {e}')
        results.append({
            'label': label, 'func_form': func, 'output': out_type,
            'variable_input': vi,
            'theta': np.nan, 'se_theta': np.nan,
            'mu_sw': np.nan, 'mu_mean': np.nan, 'mu_median': np.nan,
            'rho': np.nan, 'N': 0,
        })

# ── Summary table ────────────────────────────────────────────────────────
print(f'\n{"=" * 70}')
print('  DLS TABLE 2: MARKUPS ACROSS PRODUCTION TECHNOLOGIES')
print(f'{"=" * 70}')

rdf = pd.DataFrame(results)
for _, r in rdf.iterrows():
    print(f'  {r["label"]:25s}: theta={r["theta"]:.4f} (SE {r["se_theta"]:.4f}), '
          f'mu_sw={r["mu_sw"]:.3f}, N={r["N"]}')

# ── Save CSV ─────────────────────────────────────────────────────────────
csv_path = str(OUTPUT_DIR / 'data' / 'dls_table2.csv')
rdf.to_csv(csv_path, index=False)
print(f'\nSaved: {csv_path}')

# ── LaTeX table (DLS Table 2 format) ─────────────────────────────────────
# Build lookup: (func_form, variable_input, output) -> (mu_sw, se_theta)
lookup = {}
for _, r in rdf.iterrows():
    key = (r['func_form'], r['variable_input'], r['output'])
    lookup[key] = (r['mu_sw'], r['se_theta'])


def cell(func, vi, out_type):
    """Format a table cell as 'X.XX (X.XX)' or '---'."""
    key = (func, vi, out_type)
    if key not in lookup or np.isnan(lookup[key][0]):
        return '---'
    mu, se = lookup[key]
    return f'{mu:.2f} ({se:.2f})'


tex_path = str(OUTPUT_DIR / 'tables' / 'dls_table2.tex')
with open(tex_path, 'w') as f:
    f.write(r'\begin{table}[htbp]\centering' + '\n')
    f.write(r'\caption{Markups: Production Technologies '
            r'(Czech Construction, 2007--2021)}\label{tab:dls_table2}' + '\n')
    f.write(r'\begin{threeparttable}' + '\n')
    f.write(r'\begin{tabular}{llcc}' + '\n')
    f.write(r'\toprule' + '\n')
    f.write(r'& & Gross output & Value added \\' + '\n')
    f.write(r'\midrule' + '\n')
    # Cobb-Douglas
    f.write(r'\multirow{2}{*}{Cobb-Douglas}' + '\n')
    f.write(f'& Variable COGS & {cell("cd", "cogs", "Gross output")} '
            f'& --- \\\\\n')
    f.write(f'& Variable W & {cell("cd", "w", "Gross output")} '
            f'& {cell("cd", "w", "Value added")} \\\\[4pt]\n')
    # Translog
    f.write(r'\multirow{2}{*}{Translog}' + '\n')
    f.write(f'& Variable COGS & {cell("tl", "cogs", "Gross output")} '
            f'& --- \\\\\n')
    f.write(f'& Variable W & {cell("tl", "w", "Gross output")} '
            f'& {cell("tl", "w", "Value added")} \\\\\n')
    f.write(r'\bottomrule' + '\n')
    f.write(r'\end{tabular}' + '\n')
    f.write(r'\begin{tablenotes}\footnotesize' + '\n')
    f.write(r'\item \textit{Notes:} Sales-weighted mean markups across all '
            r'Czech construction firms (NACE 41--43), pooled 2007--2021. '
            r'Standard errors (ACH 2012 analytical) in parentheses refer to '
            r'the output elasticity $\hat{\theta}^V$. '
            r'Markup $\mu = \hat{\theta}^V / \alpha^V$ where $\alpha^V$ is '
            r'the expenditure share of the variable input. '
            r'Variable COGS: cost of goods sold ($\approx 50\%$ of GO). '
            r'Variable W: wage bill. '
            r'Value added specifications net out intermediate inputs from '
            r'output. '
            r'Replicates \citet{DeLoeckerScott2025} Table~2 format on '
            r'Czech data.' + '\n')
    f.write(r'\end{tablenotes}' + '\n')
    f.write(r'\end{threeparttable}' + '\n')
    f.write(r'\end{table}' + '\n')

print(f'Saved: {tex_path}')
print('\nDone.')
