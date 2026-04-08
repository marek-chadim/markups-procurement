"""Favoritism vs Productivity Decomposition of the Procurement Markup Premium.

Decomposes the ~14% procurement markup premium into:
  (A) Competition heterogeneity: split by single-bid share
  (B) Reform timing: pre- vs post-2012 single-bid ban (Act 55/2012)
  (C) Calibrated decomposition: total - Baranek & Titl (2024) favoritism
  (D) Oster (2019) selection robustness bounds

References
----------
Baranek & Titl (2024): Cost of Favoritism in Public Procurement, JLE 67(2).
Oster (2019): Unobservable Selection and Coefficient Stability, JBES 37(2).

Author: Marek Chadim (Yale, Tobin Center)
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings

warnings.filterwarnings('ignore')

# ================================================================== #
#  Paths
# ================================================================== #

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

DATA_PATH = INPUT_DIR / 'data.dta'
MARKUPS_PATH = OUTPUT_DIR / 'data' / 'paper_markups.dta'
TABLE_DIR = OUTPUT_DIR / 'tables'
OSTER_CSV = TABLE_DIR / 'oster_bounds.csv'
TENDERS_CSV = SCRIPT_DIR.parents[1] / '1_data' / 'input' / 'datlab' / 'master_tender_analytics.csv'


# ================================================================== #
#  Data Loading
# ================================================================== #

def load_merged() -> pd.DataFrame:
    """Load firm-year panel merged with markups."""
    df = pd.read_stata(str(DATA_PATH))
    mk = pd.read_stata(str(MARKUPS_PATH))

    # Normalize year types
    if hasattr(mk['year'].iloc[0], 'year'):
        mk['year'] = mk['year'].dt.year

    # Merge markups onto panel
    merge_cols = ['id', 'year', 'markup_A']
    if 'nace2' in mk.columns:
        merge_cols.append('nace2')
    mk_sub = mk[merge_cols].copy()
    mk_sub = mk_sub.rename(columns={'nace2': 'nace2_mk'})

    merged = df.merge(mk_sub, on=['id', 'year'], how='inner')
    merged = merged.dropna(subset=['markup_A'])
    merged = merged[merged['markup_A'] > 0].copy()
    merged['ln_mu'] = np.log(merged['markup_A'])

    print(f'Merged panel: {len(merged):,} obs, {merged["id"].nunique():,} firms, '
          f'{int(merged["year"].min())}-{int(merged["year"].max())}')
    return merged


# ================================================================== #
#  Regression Helper
# ================================================================== #

def run_premium_regression(df: pd.DataFrame) -> dict:
    """Run: log(markup_A) ~ pp_dummy + k + cogs + year x nace2 FE, cluster(id).

    Returns dict with premium, se, n_obs, n_firms.
    """
    sub = df[['id', 'year', 'ln_mu', 'pp_dummy', 'k', 'cogs', 'nace2']].dropna().copy()
    if len(sub) < 20:
        return {'premium': np.nan, 'se': np.nan, 'n_obs': len(sub),
                'n_firms': sub['id'].nunique()}

    # Year x NACE fixed effects
    yr_nace = sub['year'].astype(int).astype(str) + '_' + sub['nace2'].astype(int).astype(str)
    fe_dummies = pd.get_dummies(yr_nace, prefix='yn', drop_first=True, dtype=float)

    X = pd.concat([
        sub[['pp_dummy', 'k', 'cogs']].reset_index(drop=True),
        fe_dummies.reset_index(drop=True),
    ], axis=1)
    X = sm.add_constant(X)
    y = sub['ln_mu'].reset_index(drop=True)
    groups = sub['id'].reset_index(drop=True)

    res = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': groups})

    return {
        'premium': res.params['pp_dummy'],
        'se': res.bse['pp_dummy'],
        'n_obs': int(res.nobs),
        'n_firms': sub['id'].nunique(),
    }


# ================================================================== #
#  Test 1: Competition Heterogeneity
# ================================================================== #

def test_competition_heterogeneity(df: pd.DataFrame) -> list[dict]:
    """Split procurement firms by competition intensity."""
    print('\n=== Panel A: Competition Heterogeneity ===')
    results = []

    # Full sample
    full = run_premium_regression(df)
    full['label'] = 'Full sample'
    results.append(full)
    print(f'  Full sample: premium = {full["premium"]:.3f} '
          f'(SE {full["se"]:.3f}), N = {full["n_obs"]:,}, firms = {full["n_firms"]:,}')

    # Determine competition measure: prefer single_bid_share, fallback to avg_bids
    pp_firms = df[df['pp_dummy'] == 1].copy()

    has_sbs = pp_firms['single_bid_share'].notna().sum() > 0.5 * len(pp_firms)
    has_ab = pp_firms['avg_bids'].notna().sum() > 0.5 * len(pp_firms)

    if has_sbs:
        # Compute firm-level mean single_bid_share
        firm_comp = pp_firms.groupby('id')['single_bid_share'].mean()
        med = firm_comp.median()
        low_comp_ids = set(firm_comp[firm_comp > med].index)
        high_comp_ids = set(firm_comp[firm_comp <= med].index)
        comp_var = 'single\\_bid\\_share'
        print(f'  Competition variable: single_bid_share (median = {med:.3f})')
    elif has_ab:
        firm_comp = pp_firms.groupby('id')['avg_bids'].mean()
        med = firm_comp.median()
        low_comp_ids = set(firm_comp[firm_comp < med].index)
        high_comp_ids = set(firm_comp[firm_comp >= med].index)
        comp_var = 'avg\\_bids'
        print(f'  Competition variable: avg_bids (median = {med:.2f})')
    else:
        print('  WARNING: No competition variable available. Skipping split.')
        results.append({'label': 'Low competition', 'premium': np.nan,
                        'se': np.nan, 'n_obs': 0, 'n_firms': 0})
        results.append({'label': 'High competition', 'premium': np.nan,
                        'se': np.nan, 'n_obs': 0, 'n_firms': 0})
        return results

    # Low competition subsample: pp firms with low competition + all non-pp firms
    df_low = df[(df['id'].isin(low_comp_ids)) | (df['pp_dummy'] == 0)].copy()
    low = run_premium_regression(df_low)
    low['label'] = 'Low competition'
    results.append(low)
    print(f'  Low competition:  premium = {low["premium"]:.3f} '
          f'(SE {low["se"]:.3f}), N = {low["n_obs"]:,}')

    # High competition subsample: pp firms with high competition + all non-pp firms
    df_high = df[(df['id'].isin(high_comp_ids)) | (df['pp_dummy'] == 0)].copy()
    high = run_premium_regression(df_high)
    high['label'] = 'High competition'
    results.append(high)
    print(f'  High competition: premium = {high["premium"]:.3f} '
          f'(SE {high["se"]:.3f}), N = {high["n_obs"]:,}')

    return results


# ================================================================== #
#  Test 2: Reform Timing
# ================================================================== #

def test_reform_timing(df: pd.DataFrame) -> list[dict]:
    """Split sample pre/post 2012 single-bid ban."""
    print('\n=== Panel B: Reform Timing ===')
    results = []

    df_pre = df[df['year'] < 2012].copy()
    pre = run_premium_regression(df_pre)
    pre['label'] = 'Pre-2012'
    results.append(pre)
    print(f'  Pre-2012:  premium = {pre["premium"]:.3f} '
          f'(SE {pre["se"]:.3f}), N = {pre["n_obs"]:,}, firms = {pre["n_firms"]:,}')

    df_post = df[df['year'] >= 2012].copy()
    post = run_premium_regression(df_post)
    post['label'] = 'Post-2012'
    results.append(post)
    print(f'  Post-2012: premium = {post["premium"]:.3f} '
          f'(SE {post["se"]:.3f}), N = {post["n_obs"]:,}, firms = {post["n_firms"]:,}')

    return results


# ================================================================== #
#  Test 3: Calibrated Decomposition
# ================================================================== #

def calibrated_decomposition(total_premium: float) -> list[dict]:
    """Accounting decomposition: total = favoritism + residual."""
    print('\n=== Panel C: Calibrated Decomposition ===')

    bt_favoritism = 0.060  # Baranek & Titl (2024 JLE): +6.2% overpricing
    residual = total_premium - bt_favoritism

    results = [
        {'label': 'Total markup premium', 'value': total_premium},
        {'label': 'Favoritism channel', 'value': bt_favoritism},
        {'label': 'Residual (productivity / market power)', 'value': residual},
    ]

    print(f'  Total premium:     {total_premium:.3f}')
    print(f'  Favoritism (B&T):  {bt_favoritism:.3f}')
    print(f'  Residual:          {residual:.3f}')

    return results


# ================================================================== #
#  Test 4: Oster Bounds
# ================================================================== #

def load_oster_bounds() -> dict:
    """Load Oster bounds from CSV or return defaults."""
    print('\n=== Panel D: Oster (2019) Selection Robustness ===')

    if OSTER_CSV.exists():
        oster = pd.read_csv(OSTER_CSV)
        delta_star = oster['delta_star'].iloc[0]
        beta_star = oster['beta_star'].iloc[0]
        print(f'  Loaded from {OSTER_CSV.name}')
    else:
        # Default values from panel_treatment.py
        delta_star = -6.05
        beta_star = 0.115
        print(f'  Using default values (oster_bounds.csv not found)')

    print(f'  delta* = {delta_star:.2f}')
    print(f'  beta*  = {beta_star:.3f}')
    if delta_star < 0:
        print(f'  -> delta* < 0: unobservables would need to work in the OPPOSITE')
        print(f'     direction to explain away the premium. Very strong evidence.')

    return {'delta_star': delta_star, 'beta_star': beta_star}


# ================================================================== #
#  Test 5: Policy Reform DiD
# ================================================================== #

def run_did_regression(df: pd.DataFrame, post_year: int, reform_label: str) -> dict:
    """DiD: log(markup) ~ pp_dummy + post + pp_dummy*post + k + cogs + year x nace2 FE.

    Returns the DiD coefficient (pp_dummy x post interaction) with firm-clustered SE.
    """
    sub = df[['id', 'year', 'ln_mu', 'pp_dummy', 'k', 'cogs', 'nace2']].dropna().copy()
    if len(sub) < 20:
        return {'label': reform_label, 'did': np.nan, 'se': np.nan,
                'n_obs': len(sub), 'n_firms': sub['id'].nunique()}

    post = (sub['year'].astype(int) >= post_year).astype(float)
    pp = sub['pp_dummy'].astype(float).reset_index(drop=True)
    post = post.reset_index(drop=True)
    interaction = (pp * post).rename('pp_x_post')

    # Year x NACE fixed effects (absorbs post main effect via year FE)
    yr_nace = sub['year'].astype(int).astype(str) + '_' + sub['nace2'].astype(int).astype(str)
    fe_dummies = pd.get_dummies(yr_nace, prefix='yn', drop_first=True, dtype=float)

    X = pd.concat([
        pp.rename('pp_dummy'),
        interaction,
        sub[['k', 'cogs']].reset_index(drop=True),
        fe_dummies.reset_index(drop=True),
    ], axis=1)
    X = sm.add_constant(X)
    y = sub['ln_mu'].reset_index(drop=True)
    groups = sub['id'].reset_index(drop=True)

    res = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': groups})

    return {
        'label': reform_label,
        'did': res.params['pp_x_post'],
        'se': res.bse['pp_x_post'],
        'n_obs': int(res.nobs),
        'n_firms': sub['id'].nunique(),
    }


def test_policy_reform_did(df: pd.DataFrame) -> list[dict]:
    """Run DiD specifications for Act 55/2012 and Act 134/2016 reforms."""
    print('\n=== Panel E: Policy Reform Difference-in-Differences ===')
    results = []

    r1 = run_did_regression(df, post_year=2012, reform_label='Act 55/2012 (single-bid ban)')
    results.append(r1)
    print(f'  {r1["label"]}: beta_3 = {r1["did"]:.3f} '
          f'(SE {r1["se"]:.3f}), N = {r1["n_obs"]:,}, firms = {r1["n_firms"]:,}')

    r2 = run_did_regression(df, post_year=2016, reform_label='Act 134/2016 (MEAT criteria)')
    results.append(r2)
    print(f'  {r2["label"]}: beta_3 = {r2["did"]:.3f} '
          f'(SE {r2["se"]:.3f}), N = {r2["n_obs"]:,}, firms = {r2["n_firms"]:,}')

    return results


# ================================================================== #
#  Test 6: Rel_Price Direct Test (Engineer Cost Estimates)
# ================================================================== #

def build_firm_year_rel_price() -> pd.DataFrame:
    """Build firm-year average Rel_Price from contract-level tender data.

    Rel_Price = bid_final_price / lot_estimated_price, trimmed to (0.2, 5)
    following Baranek & Titl (2024) methodology.
    """
    if not TENDERS_CSV.exists():
        print(f'  WARNING: {TENDERS_CSV} not found. Skipping Rel_Price analysis.')
        return pd.DataFrame()

    tenders = pd.read_csv(str(TENDERS_CSV), low_memory=False)
    tenders = tenders.dropna(subset=['lot_estimated_price', 'bid_final_price'])
    tenders = tenders[
        (tenders['lot_estimated_price'] > 0) & (tenders['bid_final_price'] > 0)
    ].copy()
    tenders['rel_price'] = tenders['bid_final_price'] / tenders['lot_estimated_price']

    # Trim outliers (B&T methodology)
    tenders = tenders[(tenders['rel_price'] > 0.2) & (tenders['rel_price'] < 5)]

    # Aggregate to firm-year
    tenders['year'] = pd.to_numeric(tenders['year'], errors='coerce')
    tenders['bidder_id'] = tenders['bidder_id'].astype(str)
    firm_year = tenders.groupby(['bidder_id', 'year']).agg(
        mean_rel_price=('rel_price', 'mean'),
        n_contracts_est=('rel_price', 'count'),
    ).reset_index()
    firm_year = firm_year.rename(columns={'bidder_id': 'id_str'})

    print(f'  Tender data: {len(tenders):,} contracts -> '
          f'{len(firm_year):,} firm-year Rel_Price observations')
    return firm_year


def test_rel_price(df: pd.DataFrame) -> list[dict]:
    """Panel F: Direct test using engineer cost estimates.

    Three regressions on the subsample with Rel_Price data:
    (1) Baseline premium (subsample comparison)
    (2) Premium controlling for Rel_Price (mediation test)
    (3) Within procurement firms: markup ~ Rel_Price (correlation test)
    """
    print('\n=== Panel F: Engineer Cost Estimate Test (Rel_Price) ===')
    results = []

    firm_year_rp = build_firm_year_rel_price()
    if firm_year_rp.empty:
        for label in ['Subsample baseline', 'Controlling for Rel\\_Price',
                       'Within-PP: markup $\\sim$ Rel\\_Price']:
            results.append({'label': label, 'premium': np.nan, 'se': np.nan,
                            'n_obs': 0, 'n_firms': 0})
        return results

    # Merge Rel_Price onto the main panel
    df_rp = df.copy()
    df_rp['id_str'] = df_rp['id'].astype(int).astype(str)
    df_rp['year_int'] = df_rp['year'].astype(int)
    firm_year_rp['year_int'] = firm_year_rp['year'].astype(int)

    df_rp = df_rp.merge(
        firm_year_rp[['id_str', 'year_int', 'mean_rel_price', 'n_contracts_est']],
        on=['id_str', 'year_int'], how='left',
    )

    # Subsample: firm-years where at least one party has Rel_Price
    # (procurement firms with estimates + all non-PP firms as control)
    has_rp = df_rp['mean_rel_price'].notna()
    pp_with_rp = df_rp[has_rp & (df_rp['pp_dummy'] == 1)]
    pp_ids_with_rp = set(pp_with_rp['id'].unique())

    # Subsample: PP firms that have Rel_Price + all non-PP firms
    sub = df_rp[(df_rp['id'].isin(pp_ids_with_rp)) | (df_rp['pp_dummy'] == 0)].copy()
    # Fill Rel_Price = 0 for non-PP firms (no contracts, no overpricing)
    sub['mean_rel_price'] = sub['mean_rel_price'].fillna(0)

    print(f'  PP firms with Rel_Price: {len(pp_ids_with_rp):,}')
    print(f'  Subsample: {len(sub):,} obs ({sub["id"].nunique():,} firms)')

    # --- (1) Subsample baseline: same spec as full sample ---
    sub_base = sub[['id', 'year', 'ln_mu', 'pp_dummy', 'k', 'cogs', 'nace2']].dropna()
    if len(sub_base) < 20:
        for label in ['Subsample baseline', 'Controlling for Rel\\_Price',
                       'Within-PP: markup $\\sim$ Rel\\_Price']:
            results.append({'label': label, 'premium': np.nan, 'se': np.nan,
                            'n_obs': 0, 'n_firms': 0})
        return results

    yr_nace = sub_base['year'].astype(int).astype(str) + '_' + sub_base['nace2'].astype(int).astype(str)
    fe = pd.get_dummies(yr_nace, prefix='yn', drop_first=True, dtype=float)

    X1 = pd.concat([sub_base[['pp_dummy', 'k', 'cogs']].reset_index(drop=True),
                     fe.reset_index(drop=True)], axis=1)
    X1 = sm.add_constant(X1)
    y = sub_base['ln_mu'].reset_index(drop=True)
    groups = sub_base['id'].reset_index(drop=True)

    res1 = sm.OLS(y, X1).fit(cov_type='cluster', cov_kwds={'groups': groups})
    r1 = {'label': 'Subsample baseline',
           'premium': res1.params['pp_dummy'], 'se': res1.bse['pp_dummy'],
           'n_obs': int(res1.nobs), 'n_firms': sub_base['id'].nunique()}
    results.append(r1)
    print(f'  (1) Subsample baseline: {r1["premium"]:.3f} (SE {r1["se"]:.3f})')

    # --- (2) Controlling for Rel_Price ---
    sub2 = sub[['id', 'year', 'ln_mu', 'pp_dummy', 'k', 'cogs', 'nace2',
                 'mean_rel_price']].dropna(subset=['ln_mu', 'pp_dummy', 'k', 'cogs', 'nace2'])
    yr_nace2 = sub2['year'].astype(int).astype(str) + '_' + sub2['nace2'].astype(int).astype(str)
    fe2 = pd.get_dummies(yr_nace2, prefix='yn', drop_first=True, dtype=float)

    X2 = pd.concat([sub2[['pp_dummy', 'mean_rel_price', 'k', 'cogs']].reset_index(drop=True),
                     fe2.reset_index(drop=True)], axis=1)
    X2 = sm.add_constant(X2)
    y2 = sub2['ln_mu'].reset_index(drop=True)
    groups2 = sub2['id'].reset_index(drop=True)

    res2 = sm.OLS(y2, X2).fit(cov_type='cluster', cov_kwds={'groups': groups2})
    r2 = {'label': 'Controlling for Rel\\_Price',
           'premium': res2.params['pp_dummy'], 'se': res2.bse['pp_dummy'],
           'n_obs': int(res2.nobs), 'n_firms': sub2['id'].nunique(),
           'rp_coef': res2.params['mean_rel_price'],
           'rp_se': res2.bse['mean_rel_price']}
    results.append(r2)
    print(f'  (2) Controlling for RP: {r2["premium"]:.3f} (SE {r2["se"]:.3f}), '
          f'Rel_Price coef = {r2["rp_coef"]:.3f} (SE {r2["rp_se"]:.3f})')

    # --- (3) Within procurement firms: markup ~ Rel_Price ---
    pp_only = sub2[sub2['pp_dummy'] == 1].copy()
    pp_only = pp_only[pp_only['mean_rel_price'] > 0]  # actual Rel_Price, not filled zeros
    if len(pp_only) < 20:
        results.append({'label': 'Within-PP: markup $\\sim$ Rel\\_Price',
                        'premium': np.nan, 'se': np.nan, 'n_obs': 0, 'n_firms': 0})
        return results

    yr_nace3 = pp_only['year'].astype(int).astype(str) + '_' + pp_only['nace2'].astype(int).astype(str)
    fe3 = pd.get_dummies(yr_nace3, prefix='yn', drop_first=True, dtype=float)

    X3 = pd.concat([pp_only[['mean_rel_price', 'k', 'cogs']].reset_index(drop=True),
                     fe3.reset_index(drop=True)], axis=1)
    X3 = sm.add_constant(X3)
    y3 = pp_only['ln_mu'].reset_index(drop=True)
    groups3 = pp_only['id'].reset_index(drop=True)

    res3 = sm.OLS(y3, X3).fit(cov_type='cluster', cov_kwds={'groups': groups3})
    r3 = {'label': 'Within-PP: markup $\\sim$ Rel\\_Price',
           'premium': res3.params['mean_rel_price'], 'se': res3.bse['mean_rel_price'],
           'n_obs': int(res3.nobs), 'n_firms': pp_only['id'].nunique()}
    results.append(r3)
    print(f'  (3) Within-PP correlation: {r3["premium"]:.3f} (SE {r3["se"]:.3f})')

    return results


# ================================================================== #
#  LaTeX Table
# ================================================================== #

def write_latex_table(
    comp_results: list[dict],
    reform_results: list[dict],
    decomp_results: list[dict],
    oster: dict,
    did_results: list[dict],
    rp_results: list[dict],
    outpath: Path,
) -> None:
    """Write booktabs LaTeX table."""

    def fmt_coef(x):
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return ''
        return f'{x:.3f}'

    def fmt_se(x):
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return ''
        return f'({x:.3f})'

    def fmt_int(x):
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return ''
        return f'{int(x):,}'

    lines = []
    lines.append(r'\begin{table}[htbp]\centering')
    lines.append(r'\caption{Markup Premium Decomposition: Favoritism vs Productivity}\label{tab:favoritism}')
    lines.append(r'\begin{threeparttable}')
    lines.append(r'\begin{tabular}{lcccc}')
    lines.append(r'\toprule')
    lines.append(r'Sample / Channel & Premium & SE & $N$ & Firms \\')
    lines.append(r'\midrule')

    # Panel A
    lines.append(r'\multicolumn{5}{l}{\textit{Panel A: Competition Heterogeneity}} \\')
    for r in comp_results:
        label = r['label']
        if label.lower().startswith('low'):
            label = r'Low competition (single-bid $>$ median)'
        elif label.lower().startswith('high'):
            label = r'High competition (single-bid $\leq$ median)'
        lines.append(
            f'{label} & {fmt_coef(r["premium"])} & {fmt_se(r["se"])} '
            f'& {fmt_int(r["n_obs"])} & {fmt_int(r["n_firms"])} \\\\'
        )
    lines.append(r'\midrule')

    # Panel B
    lines.append(r'\multicolumn{5}{l}{\textit{Panel B: Reform Timing}} \\')
    for r in reform_results:
        label = r['label']
        if 'pre' in label.lower():
            label = 'Pre-2012 (before single-bid ban)'
        elif 'post' in label.lower():
            label = 'Post-2012 (after single-bid ban)'
        lines.append(
            f'{label} & {fmt_coef(r["premium"])} & {fmt_se(r["se"])} '
            f'& {fmt_int(r["n_obs"])} & {fmt_int(r["n_firms"])} \\\\'
        )
    lines.append(r'\midrule')

    # Panel C
    lines.append(r'\multicolumn{5}{l}{\textit{Panel C: Calibrated Decomposition}} \\')
    for r in decomp_results:
        label = r['label']
        if 'favoritism' in label.lower():
            label = r"Favoritism channel (Baran\'{e}k \& Titl 2024)"
        lines.append(f'{label} & {fmt_coef(r["value"])} & & & \\\\')
    lines.append(r'\midrule')

    # Panel D
    lines.append(r'\multicolumn{5}{l}{\textit{Panel D: Oster (2019) Selection Robustness}} \\')
    lines.append(
        f'$\\delta^*$ for $\\beta = 0$ & '
        f'${oster["delta_star"]:.2f}$ & & & \\\\'
    )
    lines.append(
        f'Bias-adjusted $\\beta^*$ ($\\delta = 1$) & '
        f'{oster["beta_star"]:.3f} & & & \\\\'
    )
    lines.append(r'\midrule')

    # Panel E
    lines.append(r'\multicolumn{5}{l}{\textit{Panel E: Policy Reform DiDs}} \\')
    for r in did_results:
        lines.append(
            f'{r["label"]}: $\\beta_3$ & {fmt_coef(r["did"])} & {fmt_se(r["se"])} '
            f'& {fmt_int(r["n_obs"])} & {fmt_int(r["n_firms"])} \\\\'
        )
    lines.append(r'\midrule')

    # Panel F
    lines.append(r'\multicolumn{5}{l}{\textit{Panel F: Engineer Cost Estimate Test}} \\')
    for r in rp_results:
        coef = r.get('premium', r.get('did', np.nan))
        lines.append(
            f'{r["label"]} & {fmt_coef(coef)} & {fmt_se(r["se"])} '
            f'& {fmt_int(r["n_obs"])} & {fmt_int(r["n_firms"])} \\\\'
        )

    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular}')
    lines.append(r'\begin{tablenotes}\footnotesize')
    lines.append(
        r"\item \textit{Notes:} Dependent variable is $\log(\hat{\mu}_{jt})$ "
        r"from baseline ACF Cobb-Douglas specification. "
        r"Controls: $\log K_{jt}$, $\log \mathrm{COGS}_{jt}$, year $\times$ NACE-2 fixed effects. "
        r"Standard errors clustered by firm. "
        r"Panel~A splits procurement firms by median single-bid share (higher share = less competitive). "
        r"Panel~B uses 2012 as cutoff (Act 55/2012 single-bid ban). "
        r"Panel~C decomposes the total premium using the +6\% overpricing estimate from "
        r"Baran\'{e}k \& Titl (2024, \textit{JLE}). "
        r"Panel~D reports Oster (2019) bounds: $\delta^* < 0$ means unobservables would need "
        r"to work in the opposite direction to observables to explain away the premium. "
        r"Panel~E reports difference-in-differences estimates from "
        r"$\log(\hat{\mu}_{jt}) = \beta_0 + \beta_1\,\mathrm{PP}_{jt} + \beta_2\,\mathrm{Post}_t "
        r"+ \beta_3\,(\mathrm{PP}_{jt} \times \mathrm{Post}_t) + \gamma\,k_{jt} + \delta\,\mathrm{cogs}_{jt} "
        r"+ \alpha_{t,s} + \varepsilon_{jt}$, "
        r"where $\mathrm{Post}_t$ equals one after the reform year; the DiD coefficient $\beta_3$ "
        r"captures the change in the procurement markup premium around each reform. "
        r"Panel~F uses engineering cost estimates from the procurement register "
        r"($\mathrm{Rel\_Price}_{jt} = \overline{\mathrm{price}/\mathrm{estimate}}$, "
        r"trimmed to $(0.2, 5)$ following Baran\'{e}k \& Titl 2024). "
        r"Row~1 re-estimates the baseline on the subsample with estimates. "
        r"Row~2 adds $\mathrm{Rel\_Price}$ as a control (mediation test). "
        r"Row~3 reports the within-procurement-firms correlation between markups and Rel\_Price."
    )
    lines.append(r'\end{tablenotes}')
    lines.append(r'\end{threeparttable}')
    lines.append(r'\end{table}')

    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text('\n'.join(lines), encoding='utf-8')
    print(f'\nSaved LaTeX table: {outpath}')


# ================================================================== #
#  Main
# ================================================================== #

if __name__ == '__main__':
    df = load_merged()

    # Panel A: Competition heterogeneity
    comp_results = test_competition_heterogeneity(df)

    # Panel B: Reform timing
    reform_results = test_reform_timing(df)

    # Panel C: Calibrated decomposition (use full-sample premium)
    total_premium = comp_results[0]['premium']  # full sample
    decomp_results = calibrated_decomposition(total_premium)

    # Panel D: Oster bounds
    oster = load_oster_bounds()

    # Panel E: Policy reform DiDs
    did_results = test_policy_reform_did(df)

    # Panel F: Rel_Price direct test
    rp_results = test_rel_price(df)

    # Write LaTeX table
    write_latex_table(
        comp_results=comp_results,
        reform_results=reform_results,
        decomp_results=decomp_results,
        oster=oster,
        did_results=did_results,
        rp_results=rp_results,
        outpath=TABLE_DIR / 'favoritism_decomposition.tex',
    )

    # Print summary
    print('\n' + '=' * 60)
    print('SUMMARY')
    print('=' * 60)
    print(f'  Full sample premium:       {comp_results[0]["premium"]:.3f} '
          f'(SE {comp_results[0]["se"]:.3f})')
    print(f'  Low competition premium:   {comp_results[1]["premium"]:.3f} '
          f'(SE {comp_results[1]["se"]:.3f})')
    print(f'  High competition premium:  {comp_results[2]["premium"]:.3f} '
          f'(SE {comp_results[2]["se"]:.3f})')
    print(f'  Pre-2012 premium:          {reform_results[0]["premium"]:.3f} '
          f'(SE {reform_results[0]["se"]:.3f})')
    print(f'  Post-2012 premium:         {reform_results[1]["premium"]:.3f} '
          f'(SE {reform_results[1]["se"]:.3f})')
    print(f'  Favoritism channel (B&T):  {decomp_results[1]["value"]:.3f}')
    print(f'  Residual channel:          {decomp_results[2]["value"]:.3f}')
    print(f'  Oster delta*:              {oster["delta_star"]:.2f}')
    print(f'  Oster beta*:               {oster["beta_star"]:.3f}')
    print(f'  DiD Act 55/2012 beta_3:    {did_results[0]["did"]:.3f} '
          f'(SE {did_results[0]["se"]:.3f})')
    print(f'  DiD Act 134/2016 beta_3:   {did_results[1]["did"]:.3f} '
          f'(SE {did_results[1]["se"]:.3f})')
    if rp_results and not np.isnan(rp_results[0].get('premium', np.nan)):
        print(f'  RP subsample baseline:     {rp_results[0]["premium"]:.3f} '
              f'(SE {rp_results[0]["se"]:.3f})')
        print(f'  RP-controlled premium:     {rp_results[1]["premium"]:.3f} '
              f'(SE {rp_results[1]["se"]:.3f})')
        if 'rp_coef' in rp_results[1]:
            print(f'  Rel_Price coefficient:     {rp_results[1]["rp_coef"]:.3f} '
                  f'(SE {rp_results[1]["rp_se"]:.3f})')
        print(f'  Within-PP RP correlation:  {rp_results[2]["premium"]:.3f} '
              f'(SE {rp_results[2]["se"]:.3f})')
