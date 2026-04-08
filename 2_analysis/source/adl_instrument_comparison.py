#!/usr/bin/env python3
"""ADL (2024) instrument comparison on Czech construction data.

Compares oligopoly instrument sets for production function estimation under
imperfect competition, following ADL (2024, CEPR DP 19640) Section 5.

Instrument rows (ADL Table analogs):
  Row 1: Competitors' fixed inputs (f_{-j}) — leave-one-out mean of capital
  Row 2: Own lagged variable input (standard ACF, no oligopoly instruments)
  Row 3: Competitors' fixed inputs + lagged productivity

Each row is estimated with increasing polynomial complexity of instruments.

Outputs:
  outputs/tables/adl_instrument_comparison.tex — LaTeX table
  outputs/data/adl_instrument_comparison.csv   — raw results
"""

import numpy as np
import pandas as pd
from pathlib import Path

from acf_estimator import (
    ACFEstimator, Formulation, Optimization, CWDLExtensions,
    ImperfectCompetition, options,
)

options.verbose = False

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

# ── Load and prepare data ────────────────────────────────────────────────

print('Loading data...')
df = pd.read_stata(str(INPUT_DIR / 'data.dta'))
df = df.dropna(subset=['go', 'k', 'cogs', 'pp_dummy', 'year', 'nace2']).copy()
print(f'Data: {len(df):,} obs, {df["id"].nunique():,} firms')

# Market definition: year x nace2
df['mkt_id'] = df['year'].astype(str) + '_' + df['nace2'].astype(int).astype(str)

# Market aggregates
chars = ['k', 'go', 'cogs']
mkt = df.groupby('mkt_id').agg(
    mkt_n=('id', 'count'),
    **{f'mkt_{c}_sum': (c, 'sum') for c in chars},
    mkt_go_lse=('go', lambda x: np.log(np.sum(np.exp(x)))),
).reset_index()
df = df.merge(mkt, on='mkt_id', how='left')

# Leave-one-out competitor variables
for c in chars:
    df[f'comp_{c}_sum'] = df[f'mkt_{c}_sum'] - df[c]
df['comp_n'] = df['mkt_n'] - 1
df['comp_k_mean'] = np.where(
    df['comp_n'] > 0, df['comp_k_sum'] / df['comp_n'], np.nan)
df['comp_cogs_mean'] = np.where(
    df['comp_n'] > 0, df['comp_cogs_sum'] / df['comp_n'], np.nan)

# Sufficient statistic: log-sum-exp of competitors' output
df['comp_go_lse'] = np.where(
    df['comp_n'] > 0,
    np.log(np.clip(np.exp(df['mkt_go_lse']) - np.exp(df['go']), 1e-10, None)),
    np.nan)

# Competitors' lagged productivity (omega proxy = go - 0.05*k - 0.95*cogs)
df['omega_proxy'] = df['go'] - 0.05 * df['k'] - 0.95 * df['cogs']
df = df.sort_values(['id', 'year'])
df['L_omega_proxy'] = df.groupby('id')['omega_proxy'].shift(1)
mkt_omega = df.groupby('mkt_id').agg(
    mkt_Lomega_sum=('L_omega_proxy', 'sum'),
).reset_index()
df = df.merge(mkt_omega, on='mkt_id', how='left')
df['comp_Lomega_mean'] = np.where(
    df['comp_n'] > 0,
    (df['mkt_Lomega_sum'] - df['L_omega_proxy'].fillna(0)) / df['comp_n'],
    np.nan)

# Lagged cogs for interactions
df['L_cogs'] = df.groupby('id')['cogs'].shift(1)

# Polynomial terms — competitors' fixed inputs
df['comp_k_mean_sq'] = df['comp_k_mean'] ** 2
df['comp_k_mean_x_k'] = df['comp_k_mean'] * df['k']
df['comp_k_mean_x_Lcogs'] = df['comp_k_mean'] * df['L_cogs']

# Polynomial terms — competitors' lagged productivity
df['comp_Lomega_sq'] = df['comp_Lomega_mean'] ** 2
df['comp_Lomega_x_k'] = df['comp_Lomega_mean'] * df['k']
df['comp_Lomega_x_Lcogs'] = df['comp_Lomega_mean'] * df['L_cogs']

# Cross term
df['comp_k_x_Lomega'] = df['comp_k_mean'] * df['comp_Lomega_mean']

# ── Row 4: Tender-level competition (N bidders from Datlab) ─────────────
# avg_bids is observed only for procurement firms; fill 0 for non-procurement
# (valid as interaction: avg_bids * pp_dummy captures within-procurement variation)
df['avg_bids_filled'] = df['avg_bids'].fillna(0)
df['avg_bids_sq'] = df['avg_bids_filled'] ** 2
df['avg_bids_x_k'] = df['avg_bids_filled'] * df['k']
df['avg_bids_x_Lcogs'] = df['avg_bids_filled'] * df['L_cogs']
df['avg_bids_x_compk'] = df['avg_bids_filled'] * df['comp_k_mean']
df['single_bid_filled'] = df['single_bid_share'].fillna(0)
df['n_contracts_log'] = np.log(df['n_contracts'].clip(lower=1))

print(f'Tender competition coverage: avg_bids non-null = '
      f'{df["avg_bids"].notna().sum():,} / {len(df):,}')

# Filter to firms with competitors
df = df[df['comp_go_lse'].notna() & (df['comp_n'] > 0)].copy()
print(f'After filtering: {len(df):,} obs, {df["id"].nunique():,} firms')

# ── Instrument set definitions ───────────────────────────────────────────

instrument_configs = {
    # Row 1: Competitors' fixed inputs only
    'f_k (linear)': ['comp_k_mean'],
    'f_k (quadratic)': ['comp_k_mean', 'comp_k_mean_sq'],
    'f_k (interactions)': ['comp_k_mean', 'comp_k_mean_x_k',
                           'comp_k_mean_x_Lcogs'],
    'f_k (full 2nd)': ['comp_k_mean', 'comp_k_mean_sq',
                        'comp_k_mean_x_k', 'comp_k_mean_x_Lcogs'],

    # Row 2: Standard ACF (sufficient statistic in first stage, no oligo IVs)
    'v_lag (standard)': [],

    # Row 3: Competitors' fixed inputs + lagged productivity
    'f_k+w_k (linear)': ['comp_k_mean', 'comp_Lomega_mean'],
    'f_k+w_k (quadratic)': ['comp_k_mean', 'comp_Lomega_mean',
                             'comp_k_mean_sq', 'comp_Lomega_sq'],
    'f_k+w_k (interactions)': ['comp_k_mean', 'comp_Lomega_mean',
                                'comp_k_mean_x_k', 'comp_Lomega_x_k'],
    'f_k+w_k (full 2nd)': ['comp_k_mean', 'comp_Lomega_mean',
                            'comp_k_mean_sq', 'comp_Lomega_sq',
                            'comp_k_mean_x_k', 'comp_Lomega_x_k',
                            'comp_k_x_Lomega'],

    # Row 4: Tender-level competition (N bidders from Datlab procurement register)
    'N_bids (linear)': ['avg_bids_filled'],
    'N_bids (quadratic)': ['avg_bids_filled', 'avg_bids_sq'],
    'N_bids (interactions)': ['avg_bids_filled', 'avg_bids_x_k',
                              'avg_bids_x_Lcogs'],
    'N_bids (full 2nd)': ['avg_bids_filled', 'avg_bids_sq',
                           'avg_bids_x_k', 'avg_bids_x_Lcogs'],

    # Row 5: Combined — accounting + tender competition
    'f_k+N_bids (linear)': ['comp_k_mean', 'avg_bids_filled'],
    'f_k+N_bids+w_k (full)': ['comp_k_mean', 'comp_Lomega_mean',
                               'avg_bids_filled', 'single_bid_filled',
                               'n_contracts_log'],
}

# ── Estimation loop ─────────────────────────────────────────────────────

results = []
for iv_label, iv_cols in instrument_configs.items():
    print(f'\n--- {iv_label} ---')
    valid_cols = [c for c in iv_cols
                  if c in df.columns and df[c].notna().sum() > 100]

    ic = ImperfectCompetition(
        enabled=True,
        sufficient_statistic='comp_go_lse',
        oligopoly_instruments=valid_cols,
    )
    form = Formulation(spec='cd', pp_in_markov=True, pp_interactions=True,
                       year_fe=True, nace2_fe=True)

    try:
        est = ACFEstimator(
            df, formulation=form,
            optimization=Optimization(method='nm+bfgs'),
            extensions=CWDLExtensions(survival_correction=True),
            imperfect_competition=ic, n_starts=3,
        )
        res = est.solve()

        # Production function parameters
        cogs_idx = res.beta_names.index('cogs')
        theta_v = res.betas[cogs_idx]
        se_v = res.se[cogs_idx]
        k_idx = res.beta_names.index('k')
        theta_k = res.betas[k_idx]
        se_k = res.se[k_idx]

        # Markov parameters
        markov_rho = (res.markov_coefs[1]
                      if res.markov_coefs is not None else np.nan)
        markov_pp = (res.markov_coefs[2]
                     if res.markov_coefs is not None
                     and len(res.markov_coefs) > 2
                     else np.nan)

        # Premium: raw log markup difference
        d = res.data[res.data['markup'] > 0].copy()
        d['lmu'] = np.log(d['markup'])
        pp0 = d[d['pp_dummy'] == 0]['lmu'].mean()
        pp1 = d[d['pp_dummy'] == 1]['lmu'].mean()
        premium = pp1 - pp0

        # Hansen J
        hj = res.hansen_j if res.hansen_j is not None else np.nan
        hp = res.hansen_j_pvalue if res.hansen_j_pvalue is not None else np.nan

        results.append({
            'instruments': iv_label,
            'n_oligo_iv': len(valid_cols),
            'theta_cogs': theta_v, 'se_cogs': se_v,
            'theta_k': theta_k, 'se_k': se_k,
            'rho': markov_rho, 'gamma_pp': markov_pp,
            'premium': premium,
            'hansen_j': hj, 'hansen_p': hp,
            'n_overid': res.n_overid, 'N': res.n_obs,
        })
        print(f'  theta_cogs={theta_v:.4f} (SE {se_v:.4f}), '
              f'premium={premium:.4f}, J={hj:.3f} p={hp:.3f}')
    except Exception as e:
        print(f'  FAILED: {e}')
        results.append({
            'instruments': iv_label, 'n_oligo_iv': len(valid_cols),
            'theta_cogs': np.nan, 'se_cogs': np.nan,
            'theta_k': np.nan, 'se_k': np.nan,
            'rho': np.nan, 'gamma_pp': np.nan,
            'premium': np.nan,
            'hansen_j': np.nan, 'hansen_p': np.nan,
            'n_overid': 0, 'N': 0,
        })

# ── Output ───────────────────────────────────────────────────────────────

rdf = pd.DataFrame(results)
print('\n\n=== ADL Instrument Comparison ===')
print(rdf.to_string(index=False))

# Save CSV
(OUTPUT_DIR / 'data').mkdir(parents=True, exist_ok=True)
rdf.to_csv(OUTPUT_DIR / 'data' / 'adl_instrument_comparison.csv', index=False)
print(f'\nSaved: {OUTPUT_DIR / "data" / "adl_instrument_comparison.csv"}')

# ── LaTeX table ──────────────────────────────────────────────────────────


def _fmt(x, d=4):
    if pd.isna(x):
        return '--'
    return f'{x:.{d}f}'


row1_labels = {
    'f_k (linear)': 'Linear',
    'f_k (quadratic)': 'Quadratic',
    'f_k (interactions)': 'Interactions',
    'f_k (full 2nd)': 'Full 2nd order',
}
row2_labels = {
    'v_lag (standard)': 'Standard ACF',
}
row3_labels = {
    'f_k+w_k (linear)': 'Linear',
    'f_k+w_k (quadratic)': 'Quadratic',
    'f_k+w_k (interactions)': 'Interactions',
    'f_k+w_k (full 2nd)': 'Full 2nd order',
}

lines = [
    r'\begin{table}[htbp]\centering',
    r'\caption{ADL (2024) Instrument Comparison: '
    r'Czech Construction Data}\label{tab:adl_instruments}',
    r'\begin{threeparttable}',
    r'\begin{tabular}{lccccccc}',
    r'\toprule',
    r'Instruments & $\hat{\theta}^V$ & SE & $\hat{\rho}$ '
    r'& $\hat{\gamma}_{pp}$ & Premium & $J$ & $p$ \\',
    r'\midrule',
    r"\multicolumn{8}{l}{\textit{Row 1: Competitors' "
    r"fixed inputs ($f_{-j}$)}} \\",
]

for key, label in row1_labels.items():
    r = rdf[rdf['instruments'] == key]
    if r.empty:
        continue
    r = r.iloc[0]
    lines.append(
        f'{label} & {_fmt(r["theta_cogs"])} & {_fmt(r["se_cogs"])} '
        f'& {_fmt(r["rho"], 3)} & {_fmt(r["gamma_pp"], 3)} '
        f'& {_fmt(r["premium"], 3)} & {_fmt(r["hansen_j"], 2)} '
        f'& {_fmt(r["hansen_p"], 3)} \\\\'
    )

lines.append(r'\addlinespace')
lines.append(
    r'\multicolumn{8}{l}{\textit{Row 2: Own lagged input '
    r'(standard ACF)}} \\'
)

for key, label in row2_labels.items():
    r = rdf[rdf['instruments'] == key]
    if r.empty:
        continue
    r = r.iloc[0]
    lines.append(
        f'{label} & {_fmt(r["theta_cogs"])} & {_fmt(r["se_cogs"])} '
        f'& {_fmt(r["rho"], 3)} & {_fmt(r["gamma_pp"], 3)} '
        f'& {_fmt(r["premium"], 3)} & {_fmt(r["hansen_j"], 2)} '
        f'& {_fmt(r["hansen_p"], 3)} \\\\'
    )

lines.append(r'\addlinespace')
lines.append(
    r"\multicolumn{8}{l}{\textit{Row 3: Competitors' "
    r"fixed inputs $+$ lagged productivity}} \\"
)

for key, label in row3_labels.items():
    r = rdf[rdf['instruments'] == key]
    if r.empty:
        continue
    r = r.iloc[0]
    lines.append(
        f'{label} & {_fmt(r["theta_cogs"])} & {_fmt(r["se_cogs"])} '
        f'& {_fmt(r["rho"], 3)} & {_fmt(r["gamma_pp"], 3)} '
        f'& {_fmt(r["premium"], 3)} & {_fmt(r["hansen_j"], 2)} '
        f'& {_fmt(r["hansen_p"], 3)} \\\\'
    )

row4_labels = {
    'N_bids (linear)': 'Linear',
    'N_bids (quadratic)': 'Quadratic',
    'N_bids (interactions)': 'Interactions',
    'N_bids (full 2nd)': 'Full 2nd order',
}
row5_labels = {
    'f_k+N_bids (linear)': 'Linear',
    'f_k+N_bids+w_k (full)': 'Full combined',
}

lines.append(r'\addlinespace')
lines.append(
    r'\multicolumn{8}{l}{\textit{Row 4: Tender-level competition '
    r'($N$ bidders, Datlab)}} \\'
)
for key, label in row4_labels.items():
    r = rdf[rdf['instruments'] == key]
    if r.empty:
        continue
    r = r.iloc[0]
    lines.append(
        f'{label} & {_fmt(r["theta_cogs"])} & {_fmt(r["se_cogs"])} '
        f'& {_fmt(r["rho"], 3)} & {_fmt(r["gamma_pp"], 3)} '
        f'& {_fmt(r["premium"], 3)} & {_fmt(r["hansen_j"], 2)} '
        f'& {_fmt(r["hansen_p"], 3)} \\\\'
    )

lines.append(r'\addlinespace')
lines.append(
    r'\multicolumn{8}{l}{\textit{Row 5: Combined --- accounting '
    r'$+$ tender competition}} \\'
)
for key, label in row5_labels.items():
    r = rdf[rdf['instruments'] == key]
    if r.empty:
        continue
    r = r.iloc[0]
    lines.append(
        f'{label} & {_fmt(r["theta_cogs"])} & {_fmt(r["se_cogs"])} '
        f'& {_fmt(r["rho"], 3)} & {_fmt(r["gamma_pp"], 3)} '
        f'& {_fmt(r["premium"], 3)} & {_fmt(r["hansen_j"], 2)} '
        f'& {_fmt(r["hansen_p"], 3)} \\\\'
    )

lines.extend([
    r'\bottomrule',
    r'\end{tabular}',
    r'\begin{tablenotes}\footnotesize',
    r'\item \textit{Notes:} ACF Cobb-Douglas with $pp_{it-1}$ in Markov, '
    r'year $\times$ NACE FE, pooled construction (NACE 41--43). '
    r'Sufficient statistic (comp.\ deflated revenue) included in first '
    r'stage for all specifications. Oligopoly instruments added to the '
    r'standard ACF instrument set. Row 4 uses average number of bidders '
    r'from the Datlab procurement register (filled with 0 for non-procurement '
    r'firms). Row 5 combines accounting-based and tender-level instruments. '
    r'$\rho$: productivity persistence. '
    r'$\gamma_{pp}$: procurement effect on productivity dynamics. Premium: raw '
    r'log markup difference. $J$: '
    r'Hansen overidentification test ($p$-value).',
    r'\end{tablenotes}',
    r'\end{threeparttable}',
    r'\end{table}',
])

tex = '\n'.join(lines) + '\n'
(OUTPUT_DIR / 'tables').mkdir(parents=True, exist_ok=True)
tex_path = OUTPUT_DIR / 'tables' / 'adl_instrument_comparison.tex'
with open(tex_path, 'w') as f:
    f.write(tex)
print(f'Saved: {tex_path}')

print('\nDone.')
