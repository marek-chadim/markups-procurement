#!/usr/bin/env python3
"""Instrument comparison for translog ACF on Czech construction data.

Main specification: translog with Kim-Luo-Su (2019) overidentification,
comparing seven instrument-set rows that combine internal (ADL 2024)
accounting moments with external (ABGRS 2025) exogenous shifters.

Instrument rows:
  Row 1-3: Internal accounting instruments (ADL 2024 §5)
  Row 4-5: Tender-level competition (Datlab N bidders)
  Row 6:   External raw shifters (year-level, not residualized)
  Row 7:   Combined internal + external raw
  Row 8:   External interactions RESIDUALIZED against controls
           (ABGRS 2025 §IV.D direct procedure for strong exclusion)
  Row 9:   Row 8 pool + Chamberlain (1987) optimal-instrument sieve
           compression to K_beta efficient combinations

External shifters from `1_data/output/external_panel_annual.csv`:
CZK/EUR, industry PPI, construction cost index, 10Y bond yield,
building permits, frost days, annual precipitation.

Outputs:
  outputs/tables/adl_instrument_comparison.tex — LaTeX table
  outputs/data/adl_instrument_comparison.csv   — raw results
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
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

# ── Row 6/7: External exogenous shifters (ABGRS strong exclusion) ─────
# Built by 1_data/source/build_external_panel.py from Eurostat, CNB FX,
# and Open-Meteo daily weather. Year-level shifters that are strictly
# outside the firm's productivity process, hence mean-independent of the
# ACF innovation xi_jt conditional on market×year controls.
EXTERNAL_CSV = (SCRIPT_DIR.parent.parent / '1_data' / 'output' /
                'external_panel_annual.csv')
ext_z: list[str] = []
ext_interactions: list[str] = []
if EXTERNAL_CSV.exists():
    ext = pd.read_csv(EXTERNAL_CSV)
    ext['year'] = ext['year'].astype(int)
    # Curated 7-channel shifter set (all full 2005-2023 coverage)
    ext_cols = [
        'fx_eur',                       # imported-materials cost
        'ppi_industry',                 # upstream PPI (NACE B-E)
        'construction_cost_residential',  # sector-specific material cost
        'long_rate_mcby',               # financing (10Y gov yield)
        'building_permits_sqm',         # downstream demand
        'weather_frost_days',           # outdoor activity restriction
        'weather_precipitation_sum',    # supply disruption
    ]
    ext = ext[['year'] + ext_cols]
    df['year'] = df['year'].astype(int)
    df = df.merge(ext, on='year', how='left')
    # Standardize to mean 0, std 1 for numerical stability in GMM
    for c in ext_cols:
        mu, sd = df[c].mean(), df[c].std()
        if sd and sd > 0:
            df[f'{c}_z'] = (df[c] - mu) / sd
            ext_z.append(f'{c}_z')
    # Interactions with k for firm-year variation (4 strongest channels)
    for c in ext_z[:4]:
        df[f'{c}_x_k'] = df[c] * df['k']
        ext_interactions.append(f'{c}_x_k')
    print(f'External shifters merged: {len(ext_z)} linear z-cols, '
          f'{len(ext_interactions)} k-interactions')
else:
    print(f'WARN: {EXTERNAL_CSV} not found — external rows will be empty')

# Filter to firms with competitors
df = df[df['comp_go_lse'].notna() & (df['comp_n'] > 0)].copy()
print(f'After filtering: {len(df):,} obs, {df["id"].nunique():,} firms')

# ── Row 8/9: ABGRS §IV.D direct procedure — residualize instruments ───
# Strong exclusion (ABGRS 2025, Definition 4) requires the excluded
# moment components to be mean-independent of the included controls X.
# Year-level shifters alone would be wiped out by year FE, so we form
# firm-year interactions (z_t × k_it, z_t × L_cogs_it) and then
# residualize those against the full control set. The OLS residual of
# (interaction) on (controls) has mean zero conditional on X by
# construction, satisfying the mean-independence requirement.
ext_resid_cols: list[str] = []
if ext_z:
    # L_cogs is generated inside the ACFEstimator but we need it here
    if 'L_cogs' not in df.columns:
        df = df.sort_values(['id', 'year'])
        df['L_cogs'] = df.groupby('id')['cogs'].shift(1)

    # Build the control matrix X_i: translog-basis inputs (linear only;
    # quadratic terms are what we're instrumenting for, so we can't
    # residualize against them), pp_dummy, year FE, NACE FE.
    year_dummies = pd.get_dummies(
        df['year'].astype(int), prefix='yr', drop_first=True).astype(float)
    nace_dummies = pd.get_dummies(
        df['nace2'].astype(int), prefix='nace', drop_first=True).astype(float)
    ctrl_df = pd.concat(
        [df[['k', 'L_cogs', 'pp_dummy']].astype(float).reset_index(drop=True),
         year_dummies.reset_index(drop=True),
         nace_dummies.reset_index(drop=True)],
        axis=1)
    ctrl_df.index = df.index

    def _residualize(series: pd.Series, X: pd.DataFrame) -> np.ndarray:
        """OLS residuals of series on [1, X]. NaN-safe; missing rows stay NaN."""
        y = series.values.astype(float)
        Xm = X.values.astype(float)
        mask = ~np.isnan(y) & ~np.any(np.isnan(Xm), axis=1)
        if mask.sum() < 10:
            return np.full_like(y, np.nan)
        y0 = y[mask]
        X0 = np.column_stack([np.ones(mask.sum()), Xm[mask]])
        try:
            beta, *_ = np.linalg.lstsq(X0, y0, rcond=None)
        except np.linalg.LinAlgError:
            return np.full_like(y, np.nan)
        out = np.full_like(y, np.nan)
        out[mask] = y0 - X0 @ beta
        return out

    # For each external shifter (z-scored), form z*k and z*L_cogs, then
    # residualize against the control matrix. These residuals are the
    # ABGRS-compliant strongly-excluded instruments.
    for c in ext_z:
        for interact_with, suffix in [('k', 'k'), ('L_cogs', 'Lcogs')]:
            raw = df[c] * df[interact_with]
            resid_col = f'{c}_x_{suffix}_resid'
            df[resid_col] = _residualize(raw, ctrl_df)
            ext_resid_cols.append(resid_col)

    n_ctrls = ctrl_df.shape[1]
    print(f'ABGRS residualized instruments: {len(ext_resid_cols)} cols '
          f'(7 shifters x 2 interactions), orthogonal to {n_ctrls} controls '
          f'(k, L_cogs, pp_dummy, year FE, NACE FE)')

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

    # Row 6: External exogenous shifters (ABGRS strong exclusion)
    'ext (linear)': list(ext_z),
    'ext (linear + k interactions)': list(ext_z) + list(ext_interactions),

    # Row 7: Combined — internal accounting + external shifters
    'f_k+w_k+ext (full ABGRS)':
        ['comp_k_mean', 'comp_Lomega_mean',
         'comp_k_mean_sq', 'comp_Lomega_sq',
         'comp_k_mean_x_k', 'comp_Lomega_x_k'] + list(ext_z),

    # Row 8: ABGRS strong exclusion — residualized interactions
    'ext_resid (ABGRS direct)': list(ext_resid_cols),

    # Row 9: ABGRS + Chamberlain (1987) optimal-instrument sieve
    'ext_resid + Chamberlain (optimal)': list(ext_resid_cols),
}

# Row 9 uses the same IV pool as Row 8 but triggers the Chamberlain
# `optimal_instruments='replace'` flag in Formulation, which sieve-
# projects the Jacobian onto the instrument space and returns K_beta
# optimal combinations. This is the efficient-instrument analog of
# the ABGRS direct procedure.
CHAMBERLAIN_LABELS = {'ext_resid + Chamberlain (optimal)'}

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
    form_kwargs = dict(spec='tl', overidentify=True,
                       pp_in_markov=True, pp_interactions=True,
                       year_fe=True, nace2_fe=True)
    if iv_label in CHAMBERLAIN_LABELS:
        form_kwargs['optimal_instruments'] = 'replace'
    form = Formulation(**form_kwargs)

    try:
        est = ACFEstimator(
            df, formulation=form,
            optimization=Optimization(method='nm+bfgs'),
            extensions=CWDLExtensions(survival_correction=True),
            imperfect_competition=ic, n_starts=3,
        )
        res = est.solve()

        # Production function parameters (linear terms)
        cogs_idx = res.beta_names.index('cogs')
        theta_v = res.betas[cogs_idx]
        se_v = res.se[cogs_idx]
        k_idx = res.beta_names.index('k')
        theta_k = res.betas[k_idx]
        se_k = res.se[k_idx]

        # For translog, report mean firm-year output elasticity of cogs:
        #   ∂ log y / ∂ log cogs = θ_c + 2·θ_cc·cogs + θ_kc·k
        # This is the primitive that enters the markup formula μ = θ/α.
        if form.spec == 'tl' and 'cogs2' in res.beta_names:
            try:
                c2_idx = res.beta_names.index('cogs2')
                kc_idx = res.beta_names.index('kcogs')
                d_all = res.data
                elast_c = (res.betas[cogs_idx]
                           + 2 * res.betas[c2_idx] * d_all['cogs']
                           + res.betas[kc_idx] * d_all['k'])
                theta_v_mean = float(elast_c.mean())
            except Exception:
                theta_v_mean = float(theta_v)
        else:
            theta_v_mean = float(theta_v)

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
            'theta_cogs_mean': theta_v_mean,
            'theta_k': theta_k, 'se_k': se_k,
            'rho': markov_rho, 'gamma_pp': markov_pp,
            'premium': premium,
            'hansen_j': hj, 'hansen_p': hp,
            'n_overid': res.n_overid, 'N': res.n_obs,
        })
        print(f'  theta_cogs={theta_v:.4f} (SE {se_v:.4f}), '
              f'mean elast={theta_v_mean:.4f}, '
              f'premium={premium:.4f}, J={hj:.3f} p={hp:.3f}')
    except Exception as e:
        print(f'  FAILED: {e}')
        results.append({
            'instruments': iv_label, 'n_oligo_iv': len(valid_cols),
            'theta_cogs': np.nan, 'se_cogs': np.nan,
            'theta_cogs_mean': np.nan,
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
row6_labels = {
    'ext (linear)': 'Linear z-scored',
    'ext (linear + k interactions)': r'Linear $+$ $k$ interactions',
}
row7_labels = {
    'f_k+w_k+ext (full ABGRS)': 'Combined full',
}
row8_labels = {
    'ext_resid (ABGRS direct)': 'Residualized (ABGRS direct, 14 IVs)',
}
row9_labels = {
    'ext_resid + Chamberlain (optimal)':
        r'Residualized $+$ Chamberlain sieve (optimal)',
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

lines.append(r'\addlinespace')
lines.append(
    r'\multicolumn{8}{l}{\textit{Row 6: External exogenous shifters '
    r'(ABGRS strong exclusion)}} \\'
)
for key, label in row6_labels.items():
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
    r'\multicolumn{8}{l}{\textit{Row 7: Internal accounting $+$ external '
    r'shifters (raw, not residualized)}} \\'
)
for key, label in row7_labels.items():
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
    r'\multicolumn{8}{l}{\textit{Row 8: ABGRS direct procedure --- '
    r'external interactions residualized against $X_i$}} \\'
)
for key, label in row8_labels.items():
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
    r'\multicolumn{8}{l}{\textit{Row 9: Chamberlain (1987) optimal '
    r'instruments --- Row 8 pool sieve-compressed to $K_\beta$}} \\'
)
for key, label in row9_labels.items():
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
    r'\item \textit{Notes:} ACF \textbf{translog} ($Y = \theta_k k + '
    r'\theta_c c + \theta_{kk} k^2 + \theta_{cc} c^2 + \theta_{kc} kc + '
    r'\omega$) with Kim, Luo \& Su (2019) overidentification lags '
    r'($L.k, L^2.c$), $pp_{it-1}$ in Markov, year $\times$ NACE FE, pooled '
    r'construction (NACE 41--43). Column $\hat{\theta}^V$ is the linear '
    r'cogs coefficient $\theta_c$; the firm-year output elasticity is '
    r'$\theta_c + 2\theta_{cc} c_{it} + \theta_{kc} k_{it}$ and is used to '
    r'compute the firm-specific markup $\mu_{it} = \theta^V_{it}/\alpha_{it}$. '
    r'Sufficient statistic (competitors'' deflated revenue) included in '
    r'first stage for all specifications. '
    r'Rows 1--5 use internal accounting instruments (ADL 2024 \S5). '
    r'Row 4 uses average number of bidders from the Datlab procurement '
    r'register. Row 6 uses seven external exogenous shifters (EUR/CZK, '
    r'industry PPI, construction cost index, 10Y bond yield, building '
    r'permits, frost days, precipitation), z-scored and merged on year. '
    r'Row 8 implements the \emph{direct procedure} of ABGRS (2025, '
    r'\S IV.D): each interaction $z_t \times k_{it}$ and $z_t \times '
    r'L.c_{it}$ is OLS-residualized against the controls $X_i = (k, L.c, '
    r'pp, \text{year FE}, \text{NACE FE})$ so that the moment function is '
    r'mean-independent of $X_i$ by construction, enforcing ABGRS strong '
    r'exclusion (Definition 4) and hence approximate causal consistency. '
    r'Row 9 passes the same Row 8 pool through the Chamberlain (1987) '
    r'optimal-instrument sieve via \texttt{optimal\_instruments=''replace''}: '
    r'the 14 residualized instruments are compressed to the $K_\beta = 6$ '
    r'efficient combinations $Z^* = \widehat{E}[\partial\xi/\partial\beta|Z]$. '
    r'$\rho$: productivity persistence. $\gamma_{pp}$: procurement effect '
    r'on productivity dynamics. Premium: mean log-markup difference. $J$: '
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
