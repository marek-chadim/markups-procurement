"""
Generate publication-quality summary statistics table in LaTeX.

Input:  2_analysis/input/data_rebuilt.dta
Output: 2_analysis/output/tables/summary_stats.tex
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'


def fmt(x, decimals=3):
    """Format number with commas for thousands and fixed decimals."""
    if pd.isna(x):
        return ''
    if decimals == 0:
        return f'{x:,.0f}'
    return f'{x:,.{decimals}f}'


def compute_panel_a(df):
    """Compute Panel A: Full Sample summary statistics."""
    # Define variables: (display_name, series, conditional_on_pp)
    df['sales_cogs_ratio'] = np.exp(df['go']) / np.exp(df['cogs'])

    var_specs = [
        ('Log sales (GO)',       df['go'],                False),
        ('Log COGS',             df['cogs'],              False),
        ('Log capital (K)',      df['k'],                 False),
        ('Employment',           df['empl_mid'],          False),
        ('Sales/COGS ratio',    df['sales_cogs_ratio'],   False),
        ('Market share',         df['mktshare'],          False),
        ('Procurement dummy',    df['pp_dummy'],          False),
        ('Number of contracts',  df['n_contracts'],       True),
        ('Average bidders',      df['avg_bids'],          True),
        ('Single-bid share',     df['single_bid_share'],  True),
    ]

    rows = []
    for display, series, cond_pp in var_specs:
        if cond_pp:
            s = series[df['pp_dummy'] == 1].dropna()
        else:
            s = series.dropna()
        if len(s) == 0:
            continue
        rows.append({
            'Variable': display,
            'Mean': s.mean(),
            'SD': s.std(),
            'P10': np.percentile(s, 10),
            'Median': np.percentile(s, 50),
            'P90': np.percentile(s, 90),
            'N': len(s),
        })
    return rows


def compute_panel_b(df):
    """Compute Panel B: Treated vs Control comparison."""
    prod_vars = [
        ('Log sales (GO)',  'go'),
        ('Log COGS',        'cogs'),
        ('Log capital (K)', 'k'),
        ('Employment',      'empl_mid'),
        ('Market share',    'mktshare'),
    ]

    treated = df[df['pp_dummy'] == 1]
    control = df[df['pp_dummy'] == 0]

    rows = []
    for display, col in prod_vars:
        t = treated[col].dropna()
        c = control[col].dropna()
        if len(t) == 0 or len(c) == 0:
            continue
        diff = t.mean() - c.mean()
        tstat, _ = stats.ttest_ind(t, c, equal_var=False)
        rows.append({
            'Variable': display,
            'PP1_Mean': t.mean(),
            'PP0_Mean': c.mean(),
            'Diff': diff,
            'tstat': tstat,
        })

    # NACE2 shares
    for nace in sorted(df['nace2'].dropna().unique()):
        nace_int = int(nace)
        t_share = (treated['nace2'] == nace).mean()
        c_share = (control['nace2'] == nace).mean()
        diff = t_share - c_share
        # Two-proportion z-test
        n1, n2 = len(treated), len(control)
        p_pool = (t_share * n1 + c_share * n2) / (n1 + n2)
        se = np.sqrt(p_pool * (1 - p_pool) * (1/n1 + 1/n2)) if p_pool > 0 else np.nan
        tstat = diff / se if se and se > 0 else np.nan
        rows.append({
            'Variable': f'NACE {nace_int} share',
            'PP1_Mean': t_share,
            'PP0_Mean': c_share,
            'Diff': diff,
            'tstat': tstat,
        })

    return rows


def build_latex(panel_a, panel_b, n_firms, n_obs, year_min, year_max):
    """Build the full LaTeX table string."""
    lines = []
    lines.append(r'\begin{table}[htbp]')
    lines.append(r'\centering')
    lines.append(r'\begin{threeparttable}')
    lines.append(r'\caption{Summary Statistics}')
    lines.append(r'\label{tab:summary}')

    # --- Panel A ---
    lines.append(r'\begin{tabular}{l' + 'r' * 6 + '}')
    lines.append(r'\toprule')
    lines.append(r'\multicolumn{7}{l}{\textit{Panel A: Full Sample}} \\')
    lines.append(r'\midrule')
    lines.append(r'Variable & Mean & SD & P10 & Median & P90 & N \\')
    lines.append(r'\midrule')

    for row in panel_a:
        # Use 0 decimals for N, 3 for rest
        is_count = row['Variable'] in ('Number of contracts',)
        is_dummy = row['Variable'] in ('Procurement dummy',)
        dec = 0 if is_count else 3
        n_str = fmt(row['N'], 0)
        line = (f"{row['Variable']} & {fmt(row['Mean'], dec)} & {fmt(row['SD'], dec)} & "
                f"{fmt(row['P10'], dec)} & {fmt(row['Median'], dec)} & "
                f"{fmt(row['P90'], dec)} & {n_str} \\\\")
        lines.append(line)

    lines.append(r'\midrule')

    # --- Panel B ---
    lines.append(r'\multicolumn{7}{l}{\textit{Panel B: Treated vs.\ Control}} \\')
    lines.append(r'\midrule')
    # Redefine column headers for Panel B (reuse same tabular)
    lines.append(r'Variable & PP$=$1 & PP$=$0 & Diff. & $t$-stat & & \\')
    lines.append(r'\midrule')

    for row in panel_b:
        is_share = 'share' in row['Variable'].lower()
        dec = 3
        t_str = fmt(row['tstat'], 2) if not pd.isna(row['tstat']) else ''
        line = (f"{row['Variable']} & {fmt(row['PP1_Mean'], dec)} & "
                f"{fmt(row['PP0_Mean'], dec)} & {fmt(row['Diff'], dec)} & "
                f"{t_str} & & \\\\")
        lines.append(line)

    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular}')

    # Notes
    lines.append(r'\begin{tablenotes}')
    lines.append(r'\small')
    lines.append(
        r'\item \textit{Notes:} '
        f'Sample period: {year_min}--{year_max}. '
        f'{n_firms:,} firms, {n_obs:,} firm-year observations. '
        r'Data sources: MagnusWeb (firm financials) and Datlab (procurement register). '
        r'Panel A reports full-sample moments; procurement variables conditioned on PP$=$1. '
        r'Panel B reports Welch $t$-tests for production variables and $z$-tests for NACE shares.'
    )
    lines.append(r'\end{tablenotes}')
    lines.append(r'\end{threeparttable}')
    lines.append(r'\end{table}')

    return '\n'.join(lines)


def main():
    df = pd.read_stata(str(INPUT_DIR / 'data_rebuilt.dta'))

    # Basic sample info
    n_obs = len(df)
    n_firms = df['id'].nunique()
    year_min = int(df['year'].min())
    year_max = int(df['year'].max())

    # Guard: check for columns that may not exist
    for col in ['pp_share']:
        if col not in df.columns:
            df[col] = np.nan

    panel_a = compute_panel_a(df)
    panel_b = compute_panel_b(df)
    latex = build_latex(panel_a, panel_b, n_firms, n_obs, year_min, year_max)

    os.makedirs(str(OUTPUT_DIR / 'tables'), exist_ok=True)
    with open(OUTPUT_DIR / 'tables' / 'summary_stats.tex', 'w') as f:
        f.write(latex)
    print(f'Saved: {OUTPUT_DIR / "tables" / "summary_stats.tex"}')
    print(f'Sample: {n_firms:,} firms, {n_obs:,} obs, {year_min}-{year_max}')


if __name__ == '__main__':
    main()
