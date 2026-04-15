"""
premium_timeseries.py

Anchor time-series figure: procurement markup premium by year (2006-2021),
with Czech procurement reform dates marked. For each year, regresses
log(markup_A) on pp_dummy plus NACE 2-digit dummies and input controls
(k, cogs). Standard errors clustered at the firm level.

Inputs:
    input/data_rebuilt.dta        (reference panel)
    output/paper_markups.dta      (markup_A by id-year plus pp_dummy, k, cogs, nace2)

Outputs:
    output/premium_timeseries.csv
    output/figures/premium_timeseries.pdf
"""

from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
FIG_DIR = OUTPUT_DIR / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)

YEARS = list(range(2006, 2022))
REFORMS = [
    (2012, 'Act 55/2012 (single-bid ban)', '#d62728'),
    (2016, 'Act 134/2016 (MEAT criteria)', '#2ca02c'),
    (2017, 'Register of Contracts',         '#9467bd'),
]


def estimate_year(df_year):
    """OLS of log(markup_A) on pp_dummy, NACE2 FE, k, cogs. Cluster at firm."""
    d = df_year.dropna(subset=['markup_A', 'pp_dummy', 'k', 'cogs', 'nace2', 'id']).copy()
    d = d[d['markup_A'] > 0]
    if len(d) < 20 or d['pp_dummy'].nunique() < 2:
        return None
    d['log_markup'] = np.log(d['markup_A'])
    nace_dum = pd.get_dummies(d['nace2'].astype(int), prefix='nace2', drop_first=True).astype(float)
    X = pd.concat([d[['pp_dummy', 'k', 'cogs']].astype(float), nace_dum], axis=1)
    X = sm.add_constant(X, has_constant='add')
    y = d['log_markup'].astype(float)
    groups = d['id'].astype(int).values
    model = sm.OLS(y, X, missing='drop')
    res = model.fit(cov_type='cluster', cov_kwds={'groups': groups})
    coef = float(res.params['pp_dummy'])
    se = float(res.bse['pp_dummy'])
    n_obs = int(res.nobs)
    n_firms = int(pd.Series(groups).nunique())
    return coef, se, n_obs, n_firms


def main():
    print('Loading paper_markups.dta ...')
    df = pd.read_stata(OUTPUT_DIR / 'paper_markups.dta')
    # Confirm the reference panel is present
    _ = pd.read_stata(INPUT_DIR / 'data_rebuilt.dta', columns=['id', 'year', 'pp_dummy', 'nace2'])

    rows = []
    for y in YEARS:
        sub = df[df['year'] == y]
        out = estimate_year(sub)
        if out is None:
            rows.append({'year': y, 'premium': np.nan, 'se': np.nan,
                         'ci_lo': np.nan, 'ci_hi': np.nan, 'N': 0, 'N_firms': 0})
            continue
        coef, se, n_obs, n_firms = out
        rows.append({
            'year': y,
            'premium': coef,
            'se': se,
            'ci_lo': coef - 1.96 * se,
            'ci_hi': coef + 1.96 * se,
            'N': n_obs,
            'N_firms': n_firms,
        })

    out_df = pd.DataFrame(rows)
    csv_path = OUTPUT_DIR / 'premium_timeseries.csv'
    out_df.to_csv(csv_path, index=False)

    # Print clean table
    print('\nProcurement Markup Premium by Year')
    print('-' * 64)
    print('{:>6} {:>10} {:>10} {:>10} {:>10} {:>7} {:>7}'.format(
        'Year', 'Premium', 'SE', 'CI_lo', 'CI_hi', 'N', 'Firms'))
    print('-' * 64)
    for _, r in out_df.iterrows():
        if np.isnan(r['premium']):
            print('{:>6} {:>10} {:>10} {:>10} {:>10} {:>7} {:>7}'.format(
                int(r['year']), 'NA', 'NA', 'NA', 'NA', int(r['N']), int(r['N_firms'])))
        else:
            print('{:>6d} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f} {:>7d} {:>7d}'.format(
                int(r['year']), r['premium'], r['se'],
                r['ci_lo'], r['ci_hi'], int(r['N']), int(r['N_firms'])))
    print('-' * 64)
    print('Saved: ' + str(csv_path))

    # Figure
    valid = out_df.dropna(subset=['premium'])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax.axhline(0.0, color='black', linewidth=0.7, alpha=0.6)

    for year_r, label, color in REFORMS:
        ax.axvline(year_r, color=color, linestyle='--', linewidth=1.3,
                   alpha=0.8, label=label)

    ax.fill_between(valid['year'], valid['ci_lo'], valid['ci_hi'],
                    alpha=0.2, color='#1f77b4', label='95% CI')
    ax.plot(valid['year'], valid['premium'],
            marker='o', linewidth=1.8, color='#1f77b4',
            markersize=5, label='Point estimate')

    ax.set_xlabel('Year')
    ax.set_ylabel('Procurement markup premium (log markup difference)')
    ax.set_title('Procurement Markup Premium by Year')
    ax.set_xticks(YEARS)
    ax.set_xticklabels([str(y) for y in YEARS], rotation=45, ha='right')
    ax.legend(loc='best', fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig_path = FIG_DIR / 'premium_timeseries.pdf'
    fig.savefig(fig_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('Saved: ' + str(fig_path))


if __name__ == '__main__':
    main()
