"""BMY-style markup analysis for Czech construction data.

Adapts Benkard, Miller & Yurukoglu (2026) "The Rise of Market Power: Comment"
methodology to Czech construction procurement data. Implements:
  A. Sample restriction sensitivity (original vs rebuilt)
  B. Missing data characterization
  C. Markup decomposition (within-firm, reallocation, net entry)
  D. Markup percentiles over time
  E. Variable input sensitivity (COGS vs II)
  F. Procurement subgroup decomposition

References
----------
BMY (2026): The Rise of Market Power and the Macroeconomic Implications: Comment.
DLEU (2020): The Rise of Market Power, QJE.
DLW (2012): Markups and Firm-Level Export Status, AER.

Author: Marek Chadim (Yale, Tobin Center)
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# paths
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
OUTPUT_FIG = OUTPUT_DIR / 'figures'
OUTPUT_TAB = OUTPUT_DIR / 'tables'
OUTPUT_DAT = OUTPUT_DIR / 'data'

for d in [OUTPUT_FIG, OUTPUT_TAB, OUTPUT_DAT]:
    d.mkdir(parents=True, exist_ok=True)

# import ACF estimator
from acf_estimator import (ACFEstimator, Formulation, Optimization,
                           estimate_by_industry, options as acf_options)

acf_options.verbose = False


# ========================================================================== #
#  Data Loading
# ========================================================================== #

def load_data(path: Optional[str] = None) -> pd.DataFrame:
    """Load rebuilt Czech construction data."""
    if path is None:
        path = str(INPUT_DIR / 'data_rebuilt.dta')
    df = pd.read_stata(path)
    return df


def load_original_ids() -> set:
    """Load firm-year IDs from original (pre-rebuild) dataset."""
    orig_path = INPUT_DIR / 'data.dta'
    if orig_path.exists():
        orig = pd.read_stata(str(orig_path))
        return set(zip(orig['id'].astype(int), orig['year'].astype(int)))
    return set()


# ========================================================================== #
#  A. Sample Restriction Sensitivity (BMY Figure 1 analog)
# ========================================================================== #

def compute_markups_by_sample(df: pd.DataFrame) -> pd.DataFrame:
    """Estimate ACF markups on original vs rebuilt samples, COGS vs II."""
    results_list = []

    for var_input, label in [('cogs', 'COGS'), ('ii', 'II')]:
        formulation_kw = {'variable_input': var_input} if var_input == 'ii' else {}

        for sample_name, mask_series in [('rebuilt', pd.Series(True, index=df.index)),
                                          ('original', df['_in_original'])]:
            sub = df.loc[mask_series].copy()
            if len(sub) < 100:
                continue

            try:
                print(f"    Estimating {sample_name}/{label} (N={len(sub)})...")
                _, _, markups = estimate_by_industry(
                    sub, specs=('cd',),
                    formulation_kwargs=formulation_kw,
                )
                # merge go back from sub for sales computation
                markups = markups.merge(
                    sub[['id', 'year', 'go', 'pp_dummy']].drop_duplicates(),
                    on=['id', 'year'], how='left', suffixes=('', '_orig')
                )
                markups['sample'] = sample_name
                markups['var_input'] = label
                markups['sales'] = np.exp(markups['go'])
                results_list.append(markups)
            except Exception as e:
                print(f"  Warning: {sample_name}/{label} failed: {e}")

    if not results_list:
        return pd.DataFrame()
    return pd.concat(results_list, ignore_index=True)


def aggregate_markups(mu: pd.DataFrame) -> pd.DataFrame:
    """Sales-weighted aggregate markup by year, sample, var_input."""
    rows = []
    for (sample, vi, year), g in mu.groupby(['sample', 'var_input', 'year']):
        w = g['sales'] / g['sales'].sum()
        agg_mu = (w * g['markup']).sum()
        rows.append({'sample': sample, 'var_input': vi, 'year': year,
                     'markup_agg': agg_mu, 'n_firms': len(g)})
    return pd.DataFrame(rows)


def figure_sample_sensitivity(agg: pd.DataFrame):
    """BMY Figure 1 analog: markup trends by sample/variable input."""
    fig, ax = plt.subplots(figsize=(8, 5))

    styles = {
        ('rebuilt', 'COGS'): ('solid', 'black', 'Rebuilt (COGS)'),
        ('original', 'COGS'): ('--', 'red', 'Original (COGS)'),
        ('rebuilt', 'II'): ('-.', 'blue', 'Rebuilt (II)'),
    }

    for (sample, vi), style in styles.items():
        sub = agg[(agg['sample'] == sample) & (agg['var_input'] == vi)]
        if len(sub) > 0:
            ax.plot(sub['year'], sub['markup_agg'], linestyle=style[0],
                    color=style[1], label=style[2], linewidth=2)

    ax.set_xlabel('')
    ax.set_ylabel('Sales-Weighted Markup')
    ax.legend(frameon=False, loc='upper left')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title('Sample Restriction Sensitivity')
    fig.tight_layout()
    fig.savefig(OUTPUT_FIG / 'bmy_sample_sensitivity.pdf', dpi=150)
    plt.close(fig)
    print(f"  Saved: {OUTPUT_FIG / 'bmy_sample_sensitivity.pdf'}")


# ========================================================================== #
#  C. Markup Decomposition (BMY Figure 3 analog)
# ========================================================================== #

def markup_decomposition(
    mu: pd.DataFrame,
    base_year: int = 2006,
    label: str = 'all',
) -> pd.DataFrame:
    """Olley-Pakes style decomposition following BMY Markup_Figures.do lines 168-218.

    Decomposes Δμ_agg into:
      - within-firm: Σ(Δμ_i × L.share_i)
      - reallocation: Σ(Δshare_i × L.demeaned_μ_i) + Σ(Δshare_i × Δμ_i)
      - net entry: residual

    Parameters
    ----------
    mu : DataFrame with id, year, markup, sales
    base_year : initialization year
    label : label for this decomposition
    """
    panel = mu[['id', 'year', 'markup', 'sales']].dropna().copy()
    panel = panel.sort_values(['id', 'year'])

    years = sorted(panel['year'].unique())
    years = [y for y in years if y >= base_year]

    decomp_rows = []
    for year in years:
        curr = panel[panel['year'] == year].copy()
        prev = panel[panel['year'] == year - 1].copy()

        if len(curr) == 0 or len(prev) == 0:
            continue

        # current aggregates
        curr['share'] = curr['sales'] / curr['sales'].sum()
        mu_agg = (curr['share'] * curr['markup']).sum()

        # previous aggregates
        prev['share'] = prev['sales'] / prev['sales'].sum()
        mu_agg_prev = (prev['share'] * prev['markup']).sum()

        # merge for continuing firms
        cont = pd.merge(
            prev[['id', 'markup', 'share']].rename(columns={'markup': 'mu_prev', 'share': 'share_prev'}),
            curr[['id', 'markup', 'share']].rename(columns={'markup': 'mu_curr', 'share': 'share_curr'}),
            on='id', how='inner'
        )

        if len(cont) == 0:
            continue

        cont['d_mu'] = cont['mu_curr'] - cont['mu_prev']
        cont['d_share'] = cont['share_curr'] - cont['share_prev']
        cont['mu_prev_demean'] = cont['mu_prev'] - mu_agg_prev

        # within-firm: Σ(Δμ_i × L.share_i)
        within = (cont['d_mu'] * cont['share_prev']).sum()
        # reallocation: Σ(Δshare_i × L.demeaned_μ_i) + cross-term Σ(Δshare_i × Δμ_i)
        realloc_pure = (cont['d_share'] * cont['mu_prev_demean']).sum()
        cross = (cont['d_share'] * cont['d_mu']).sum()
        reallocation = realloc_pure + cross

        d_mu_agg = mu_agg - mu_agg_prev
        net_entry = d_mu_agg - within - reallocation

        decomp_rows.append({
            'year': year,
            'mu_agg': mu_agg,
            'd_mu_agg': d_mu_agg,
            'within': within,
            'reallocation': reallocation,
            'net_entry': net_entry,
            'n_continuing': len(cont),
            'n_total': len(curr),
            'label': label,
        })

    decomp = pd.DataFrame(decomp_rows)

    # cumulate from base year
    if len(decomp) > 0:
        base_mu = panel[panel['year'] == base_year]
        if len(base_mu) > 0:
            base_mu_val = (base_mu['sales'] / base_mu['sales'].sum() * base_mu['markup']).sum()
        else:
            base_mu_val = decomp['mu_agg'].iloc[0]

        decomp['within_cum'] = decomp['within'].cumsum() + base_mu_val
        decomp['realloc_cum'] = decomp['reallocation'].cumsum() + base_mu_val
        decomp['entry_cum'] = decomp['net_entry'].cumsum() + base_mu_val

    return decomp


def figure_decomposition(decomp: pd.DataFrame, filename: str = 'bmy_decomposition.pdf',
                         title: str = 'Decomposition of Markup Growth'):
    """BMY Figure 3 analog: within, reallocation, net entry."""
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(decomp['year'], decomp['mu_agg'], 'k-.',
            linewidth=2.5, label='Markup (aggregate)')
    ax.plot(decomp['year'], decomp['within_cum'], 'r--',
            linewidth=2, label='Within-firm')
    ax.plot(decomp['year'], decomp['realloc_cum'], 'b:',
            linewidth=2, label='Reallocation')
    ax.plot(decomp['year'], decomp['entry_cum'], 'g-',
            linewidth=1.5, label='Net Entry')

    ax.set_ylabel('Sales-Weighted Markup')
    ax.set_xlabel('')
    ax.legend(frameon=False, loc='upper left')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(OUTPUT_FIG / filename, dpi=150)
    plt.close(fig)
    print(f"  Saved: {OUTPUT_FIG / filename}")


# ========================================================================== #
#  D. Markup Percentiles (BMY Figure 4 analog)
# ========================================================================== #

def figure_percentiles(mu: pd.DataFrame, filename: str = 'bmy_percentiles.pdf'):
    """BMY Figure 4 analog: p10/p25/p50/p75/p90 over time."""
    years = sorted(mu['year'].unique())
    pcts = {p: [] for p in [10, 25, 50, 75, 90]}
    means = []

    for year in years:
        sub = mu[mu['year'] == year]
        w = sub['sales'] / sub['sales'].sum()
        # weighted percentiles (approximate via sorting)
        sorted_idx = sub['markup'].argsort()
        sorted_mu = sub['markup'].iloc[sorted_idx].values
        sorted_w = w.iloc[sorted_idx].values
        cum_w = np.cumsum(sorted_w)

        for p in pcts:
            idx = np.searchsorted(cum_w, p / 100)
            idx = min(idx, len(sorted_mu) - 1)
            pcts[p].append(sorted_mu[idx])
        means.append((w * sub['markup']).sum())

    fig, ax = plt.subplots(figsize=(8, 5))
    styles = {90: '--', 75: ':', 50: '-.', 25: ':', 10: '--'}
    for p, vals in pcts.items():
        ax.plot(years, vals, 'r' + styles[p], linewidth=1.5, label=f'P{p}')
    ax.plot(years, means, 'r-', linewidth=2.5, label='Average')

    ax.set_ylabel('Markup')
    ax.set_xlabel('')
    ax.legend(frameon=False, loc='upper left', ncol=2)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title('Markup Percentiles Over Time')
    fig.tight_layout()
    fig.savefig(OUTPUT_FIG / filename, dpi=150)
    plt.close(fig)
    print(f"  Saved: {OUTPUT_FIG / filename}")


# ========================================================================== #
#  E. Variable Input Sensitivity (BMY F&I analog)
# ========================================================================== #

def figure_variable_input(agg: pd.DataFrame, filename: str = 'bmy_variable_input.pdf'):
    """Side-by-side COGS vs II markup trends."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # left: aggregate markup trends
    for vi, style, color in [('COGS', '-', 'black'), ('II', '--', 'blue')]:
        sub = agg[(agg['sample'] == 'rebuilt') & (agg['var_input'] == vi)]
        if len(sub) > 0:
            ax1.plot(sub['year'], sub['markup_agg'], linestyle=style,
                     color=color, label=vi, linewidth=2)

    ax1.set_ylabel('Sales-Weighted Markup')
    ax1.set_title('Aggregate Markup by Variable Input')
    ax1.legend(frameon=False)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    # right: procurement premium by variable input (annual)
    # compute annual PP premium as ratio of means
    ax2.set_title('Procurement Premium by Variable Input')
    ax2.set_ylabel('Premium (log points)')
    ax2.axhline(0, color='gray', linewidth=0.5)
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    fig.tight_layout()
    fig.savefig(OUTPUT_FIG / filename, dpi=150)
    plt.close(fig)
    print(f"  Saved: {OUTPUT_FIG / filename}")


# ========================================================================== #
#  B. Missing Data Characterization (BMY Table 1 analog)
# ========================================================================== #

def table_sample_comparison(df: pd.DataFrame) -> pd.DataFrame:
    """Characterize original-only vs rebuilt-only observations."""
    df = df.copy()
    df['ln_sales'] = df['go']
    df['ln_cogs'] = df['cogs']
    df['s_g'] = np.exp(df['go']) / np.exp(df['cogs'])

    # summary by sample flag
    rows = []
    for flag in ['both', 'rebuilt_only', 'original_only']:
        sub = df[df['_sample_flag'] == flag]
        if len(sub) == 0:
            continue
        rows.append({
            'Sample': flag,
            'N': len(sub),
            'N_firms': sub['id'].nunique(),
            'mean_ln_sales': sub['ln_sales'].mean(),
            'mean_ln_cogs': sub['ln_cogs'].mean(),
            'mean_s_g': sub['s_g'].mean(),
            'mean_pp': sub['pp_dummy'].mean() if 'pp_dummy' in sub else np.nan,
        })

    result = pd.DataFrame(rows)
    # save as LaTeX
    latex = result.to_latex(index=False, float_format='%.3f')
    (OUTPUT_TAB / 'bmy_sample_comparison.tex').write_text(latex)
    print(f"  Saved: {OUTPUT_TAB / 'bmy_sample_comparison.tex'}")
    return result


# ========================================================================== #
#  Main Pipeline
# ========================================================================== #

def main():
    print('=' * 70)
    print('  BMY-STYLE ANALYSIS: CZECH CONSTRUCTION MARKUPS')
    print('=' * 70)

    # load data
    print('\n[1] Loading data...')
    df = load_data()
    print(f"  Rebuilt sample: {len(df)} obs, {df['id'].nunique()} firms")

    # tag original sample
    orig_ids = load_original_ids()
    if orig_ids:
        df['_in_original'] = list(zip(df['id'].astype(int), df['year'].astype(int)))
        df['_in_original'] = df['_in_original'].isin(orig_ids)
        df['_sample_flag'] = 'both'
        df.loc[~df['_in_original'], '_sample_flag'] = 'rebuilt_only'
        n_orig = df['_in_original'].sum()
        print(f"  Original sample: {n_orig} obs matched")
        print(f"  Rebuilt-only: {(~df['_in_original']).sum()} obs")
    else:
        df['_in_original'] = True
        df['_sample_flag'] = 'both'
        print("  Warning: original data not found, skipping sample comparison")

    # A. Sample sensitivity — estimate markups under different samples/inputs
    print('\n[2] Estimating markups by sample and variable input...')
    mu_all = compute_markups_by_sample(df)

    if len(mu_all) > 0:
        agg = aggregate_markups(mu_all)
        figure_sample_sensitivity(agg)

        # E. Variable input sensitivity
        print('\n[5] Variable input sensitivity figure...')
        figure_variable_input(agg)
    else:
        print("  Warning: markup estimation failed, using fallback")
        mu_all = pd.DataFrame()
        agg = pd.DataFrame()

    # C. Decomposition (use rebuilt/COGS markups)
    print('\n[3] Markup decomposition...')
    if len(mu_all) > 0 and 'sample' in mu_all.columns:
        mu_decomp = mu_all[(mu_all['sample'] == 'rebuilt') & (mu_all['var_input'] == 'COGS')]
    else:
        mu_decomp = pd.DataFrame()
    if len(mu_decomp) > 0:
        decomp = markup_decomposition(mu_decomp, base_year=2006, label='all')
        if len(decomp) > 0:
            figure_decomposition(decomp, 'bmy_decomposition.pdf',
                                 'Decomposition of Markup Growth — Czech Construction')

            # F. Procurement subgroup decomposition
            print('\n[6] Procurement subgroup decomposition...')
            for pp_val, pp_label in [(0, 'non_procurement'), (1, 'procurement')]:
                mu_pp = mu_decomp[mu_decomp['pp_dummy'] == pp_val]
                if len(mu_pp) > 100:
                    decomp_pp = markup_decomposition(mu_pp, base_year=2006, label=pp_label)
                    if len(decomp_pp) > 0:
                        figure_decomposition(
                            decomp_pp,
                            f'bmy_decomposition_{pp_label}.pdf',
                            f'Decomposition — {pp_label.replace("_", " ").title()} Firms'
                        )

    # D. Percentiles
    print('\n[4] Markup percentiles...')
    if len(mu_decomp) > 0:
        figure_percentiles(mu_decomp, 'bmy_percentiles.pdf')

    # B. Sample comparison table
    print('\n[7] Sample comparison table...')
    table_sample_comparison(df)

    # save combined markup data
    if len(mu_all) > 0:
        mu_all.to_stata(str(OUTPUT_DAT / 'bmy_czech_markups.dta'), write_index=False)
        print(f"\n  Saved: {OUTPUT_DAT / 'bmy_czech_markups.dta'}")

    print('\n' + '=' * 70)
    print('  BMY ANALYSIS COMPLETE')
    print('=' * 70)


if __name__ == '__main__':
    main()
