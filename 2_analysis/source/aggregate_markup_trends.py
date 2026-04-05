"""
Aggregate Markup Trends for Czech Construction.

Computes sales-weighted aggregate markups following De Loecker, Eeckhout, and
Unger (2020, QJE) "The Rise of Market Power and the Macroeconomic Implications",
with Olley-Pakes decomposition of year-over-year changes into within, between,
cross, and net-entry components.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'


def main():
    # Load data
    df = pd.read_stata(str(INPUT_DIR / 'data.dta'))
    mk = pd.read_stata(str(OUTPUT_DIR / 'data' / 'paper_markups.dta'))

    # Handle year dtype
    for d in [df, mk]:
        if hasattr(d['year'].iloc[0], 'year'):
            d['year'] = d['year'].dt.year
        d['year'] = d['year'].astype(int)

    # Merge: need sales (exp(go)) and markup_A
    panel = df[['id', 'year', 'nace2', 'go']].merge(
        mk[['id', 'year', 'markup_A']], on=['id', 'year'])
    panel = panel[(panel['markup_A'] > 0) & panel['markup_A'].notna()].copy()
    panel['sales'] = np.exp(panel['go'])

    # Compute sales-weighted aggregate markup by year
    results = []
    for year, grp in panel.groupby('year'):
        total_sales = grp['sales'].sum()
        grp = grp.copy()
        grp['share'] = grp['sales'] / total_sales
        m_sw = (grp['share'] * grp['markup_A']).sum()  # sales-weighted
        m_mean = grp['markup_A'].mean()  # simple mean
        # Sales-weighted median
        sorted_grp = grp.sort_values('markup_A')
        sorted_grp['cum_share'] = sorted_grp['share'].cumsum()
        median_idx = (sorted_grp['cum_share'] >= 0.5).idxmax()
        m_sw_median = sorted_grp.loc[median_idx, 'markup_A']
        results.append({
            'year': year, 'M_sw': m_sw, 'M_mean': m_mean, 'M_sw_median': m_sw_median,
            'N': len(grp), 'total_sales': total_sales
        })

    agg = pd.DataFrame(results).sort_values('year')
    print('Aggregate Markup Trends (Czech Construction):')
    print(agg[['year', 'M_sw', 'M_mean', 'M_sw_median', 'N']].to_string(index=False))

    # Olley-Pakes decomposition
    # Panel form: markup_it, share_it
    panel['sales_share'] = panel.groupby('year')['sales'].transform(lambda s: s / s.sum())

    decomp = []
    years = sorted(panel['year'].unique())
    for t in years[1:]:
        curr = panel[panel['year'] == t].set_index('id')
        prev = panel[panel['year'] == (t - 1)].set_index('id')
        incumbents = curr.index.intersection(prev.index)
        entrants = curr.index.difference(prev.index)
        exiters = prev.index.difference(curr.index)

        if len(incumbents) == 0:
            continue
        c = curr.loc[incumbents]
        p = prev.loc[incumbents]

        dM_within = (p['sales_share'] * (c['markup_A'] - p['markup_A'])).sum()
        dM_between = ((c['sales_share'] - p['sales_share']) * p['markup_A']).sum()
        dM_cross = ((c['sales_share'] - p['sales_share']) * (c['markup_A'] - p['markup_A'])).sum()

        # Net entry
        dM_entry = (curr.loc[entrants, 'sales_share'] * curr.loc[entrants, 'markup_A']).sum() if len(entrants) > 0 else 0
        dM_exit = -(prev.loc[exiters, 'sales_share'] * prev.loc[exiters, 'markup_A']).sum() if len(exiters) > 0 else 0

        dM_total = agg.loc[agg['year'] == t, 'M_sw'].iloc[0] - agg.loc[agg['year'] == (t - 1), 'M_sw'].iloc[0]

        decomp.append({
            'year': t, 'dM_total': dM_total,
            'within': dM_within, 'between': dM_between, 'cross': dM_cross,
            'entry': dM_entry, 'exit': dM_exit,
            'N_incumbents': len(incumbents), 'N_entrants': len(entrants), 'N_exiters': len(exiters)
        })

    decomp_df = pd.DataFrame(decomp)
    print('\nOlley-Pakes Decomposition (annual changes):')
    print(decomp_df[['year', 'dM_total', 'within', 'between', 'cross', 'entry', 'exit']].to_string(index=False))

    # Save outputs
    agg.to_csv(OUTPUT_DIR / 'aggregate_markup_trends.csv', index=False)
    decomp_df.to_csv(OUTPUT_DIR / 'aggregate_markup_decomposition.csv', index=False)

    # Figure: aggregate trends
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(agg['year'], agg['M_sw'], 'o-', label='Sales-weighted mean', linewidth=2)
    ax.plot(agg['year'], agg['M_mean'], 's--', label='Simple mean', linewidth=1.5, alpha=0.7)
    ax.plot(agg['year'], agg['M_sw_median'], '^:', label='Sales-weighted median', linewidth=1.5, alpha=0.7)
    ax.axvline(2012, color='red', linestyle='--', alpha=0.5, label='Act 55/2012')
    ax.axvline(2016, color='green', linestyle='--', alpha=0.5, label='Act 134/2016')
    ax.set_xlabel('Year')
    ax.set_ylabel('Aggregate Markup')
    ax.set_title('Aggregate Markup Trends, Czech Construction')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'figures' / 'aggregate_markup_trends.pdf', dpi=200, bbox_inches='tight')
    print(f'\nSaved: figures/aggregate_markup_trends.pdf')

    # LaTeX table
    tex = [
        r'\begin{table}[htbp]\centering',
        r'\caption{Aggregate Markup Trends, Czech Construction}\label{tab:aggregate_trends}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lcccc}',
        r'\toprule',
        r'Year & Sales-weighted & Simple mean & SW Median & N firms \\',
        r'\midrule'
    ]
    for _, r in agg.iterrows():
        tex.append(f'{int(r["year"])} & {r["M_sw"]:.3f} & {r["M_mean"]:.3f} & {r["M_sw_median"]:.3f} & {int(r["N"]):,} \\\\')
    tex += [
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r'\item \textit{Notes:} Sales-weighted aggregate markup following De Loecker, Eeckhout, and Unger (2020). $M_t^{sw} = \sum_i (\text{sales}_{it}/\sum_j \text{sales}_{jt}) \times \mu_{it}^A$ where $\mu^A$ is the ACF Cobb-Douglas markup with procurement in the Markov transition. Construction panel: CZ-NACE F (41-43).',
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}'
    ]
    with open(OUTPUT_DIR / 'tables' / 'aggregate_markup_trends.tex', 'w') as f:
        f.write('\n'.join(tex))
    print('Saved: tables/aggregate_markup_trends.tex')


if __name__ == '__main__':
    main()
