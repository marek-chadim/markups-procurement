"""Rebuild estimation data from raw CSVs with proper merge.

Replaces the m:m merge in Create_Data.do with a correct pipeline:
  1. Load raw CSVs (financial, ratios, selections, deflators, tenders)
  2. Deduplicate each dataset on (ID, Year) BEFORE merging
  3. Merge 1:1 on (ID, Year) — no cartesian products
  4. Apply same cleaning rules as Create_Data.do
  5. Compare against existing magnus.dta to quantify differences

Conflict resolution strategy for duplicates with different values:
  - Financial: keep the row with median sales (avoid extremes from
    consolidated vs. unconsolidated reporting)
  - Ratios: keep the row with cm closest to industry median (~40-50%)

Author: Marek Chadim (Yale, Tobin Center)
Date: March 2026
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore', category=FutureWarning)

# ========================================================================== #
#  Paths
# ========================================================================== #

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

CSV_DIR = INPUT_DIR / 'magnusweb'
OUT_DIR = OUTPUT_DIR


def header(title: str) -> None:
    print(f'\n{"=" * 70}')
    print(f'  {title}')
    print(f'{"=" * 70}\n')


# ========================================================================== #
#  Step 1: Load raw CSVs
# ========================================================================== #

def load_financial() -> pd.DataFrame:
    """Load and concatenate financial CSV files."""
    dfs = [pd.read_csv(CSV_DIR / f'financial{i}.csv') for i in range(1, 7)]
    df = pd.concat(dfs, ignore_index=True)
    df = df.rename(columns={
        'ID': 'id', 'Year': 'year',
        'C - Costs': 'costs',
        'FA - Fixed assets': 'assets',
        'Sal - Sales, Outputs': 'sales',
    })
    print(f'  Financial: {len(df):,} rows from 6 CSVs')
    return df


def load_ratios() -> pd.DataFrame:
    """Load and concatenate ratios CSV files."""
    dfs = [pd.read_csv(CSV_DIR / f'ratios{i}.csv') for i in range(1, 4)]
    df = pd.concat(dfs, ignore_index=True)
    df = df.rename(columns={
        'ID': 'id', 'Year': 'year',
        'WVA - Wages / Value added': 'wva',
        'WS - Wages / Sales': 'ws',
        'LP - Labour productivity': 'lp',
        'CM III - Contribution margin': 'cm',
    })
    print(f'  Ratios: {len(df):,} rows from 3 CSVs')
    return df


def _empl_cat_to_midpoint(s) -> float:
    """Convert employment category string to numeric midpoint.

    E.g., '10 - 19 employees' → 14.5, '250 - 499 employees' → 374.5.
    """
    if pd.isna(s):
        return np.nan
    s = str(s).strip()
    if '-' not in s:
        return np.nan
    parts = s.split('-')
    try:
        lo = int(parts[0].strip().replace(' ', '').replace(',', ''))
        hi_str = parts[1].strip().split()[0].replace(' ', '').replace(',', '')
        hi = int(hi_str)
        return (lo + hi) / 2
    except (ValueError, IndexError):
        return np.nan


def load_selections() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load selections: firm-level attributes + time-varying employment panel.

    Returns
    -------
    firm_info : DataFrame
        One row per firm: id, nace, empl_num, legal_form, inst_sector.
    empl_panel : DataFrame
        Long panel: id, year, empl_mid (midpoint of employment category).
    """
    dfs = [pd.read_csv(CSV_DIR / f'selections{i}.csv', encoding='latin1')
           for i in range(1, 6)]
    df = pd.concat(dfs, ignore_index=True)
    # First column has encoding issue in name
    col0 = df.columns[0]
    df = df.rename(columns={col0: 'id'})
    df['id'] = pd.to_numeric(df['id'], errors='coerce')
    df = df.dropna(subset=['id'])

    # Identify columns
    nace_col = [c for c in df.columns if 'NACE' in c.upper()][0]
    empl_col = 'Number of employees'
    lf_col = 'Legal form'
    inst_col = 'Institutional sectors (ESA 2010)'

    # --- Firm-level attributes (one row per firm) ---
    keep_cols = ['id', nace_col, empl_col]
    if lf_col in df.columns:
        keep_cols.append(lf_col)
    if inst_col in df.columns:
        keep_cols.append(inst_col)
    firm_info = df[keep_cols].copy()
    firm_info = firm_info.rename(columns={
        nace_col: 'nace', empl_col: 'empl_num',
        lf_col: 'legal_form', inst_col: 'inst_sector',
    })
    firm_info['nace'] = pd.to_numeric(firm_info['nace'], errors='coerce')
    firm_info['empl_num'] = pd.to_numeric(firm_info['empl_num'], errors='coerce')
    # Foreign ownership dummy: 1 if foreign-controlled
    if 'inst_sector' in firm_info.columns:
        firm_info['foreign'] = firm_info['inst_sector'].str.contains(
            'foreign', case=False, na=False).astype(int)
    firm_info = firm_info.drop_duplicates(subset=['id'], keep='first')
    print(f'  Selections (firm-level): {len(firm_info):,} firms')

    # --- Employment panel (long format, one row per firm-year) ---
    empl_yr_cols = {c: int(c[:4]) for c in df.columns
                    if 'classification' in c.lower() and c[0].isdigit()}
    rows = []
    for _, row in df.drop_duplicates(subset=['id']).iterrows():
        fid = row['id']
        for col, yr in empl_yr_cols.items():
            mid = _empl_cat_to_midpoint(row[col])
            if not np.isnan(mid):
                rows.append({'id': fid, 'year': yr, 'empl_mid': mid})
    empl_panel = pd.DataFrame(rows)
    if len(empl_panel) > 0:
        empl_panel = empl_panel.drop_duplicates(subset=['id', 'year'])
    print(f'  Selections (empl panel): {len(empl_panel):,} firm-year obs '
          f'({empl_panel["id"].nunique():,} firms, '
          f'{len(empl_yr_cols)} years)')

    return firm_info, empl_panel


def load_deflators() -> pd.DataFrame:
    """Load deflators."""
    df = pd.read_csv(CSV_DIR / 'deflators.csv')
    df.columns = df.columns.str.lower()
    print(f'  Deflators: {len(df):,} rows')
    return df


def load_tenders() -> pd.DataFrame:
    """Load tender data and construct firm-year competition measures.

    Variables constructed (per firm-year):
      pp_sales       : total procurement revenue (sum of winning bids)
      n_contracts    : number of contracts won
      avg_bids       : average number of bidders on won contracts
      max_bids       : max bidders on any single contract
      single_bid_share : fraction of contracts with only 1 bidder
      avg_discount   : mean (est_price - bid) / est_price
      hhi_revenue    : Herfindahl of contract values (revenue concentration)
    """
    df = pd.read_csv(INPUT_DIR / 'datlab' / 'master_tender_analytics.csv',
                     low_memory=False)
    df = df.rename(columns={'bidder_id': 'id', 'bid_final_price': 'pp_sales'})

    # Drop foreign bidders (id length == 2)
    df['id'] = df['id'].astype(str)
    df = df[df['id'].str.len() != 2]
    df['id'] = pd.to_numeric(df['id'], errors='coerce')
    df = df.dropna(subset=['id'])
    df['id'] = df['id'].astype(int)
    n_tender_rows = len(df)

    # --- FIX 1: Use signature date year where available, fall back to year ---
    df['sig_year'] = pd.to_datetime(
        df['contract_signature_date'], errors='coerce'
    ).dt.year
    n_sig = df['sig_year'].notna().sum()
    n_reassigned = ((df['sig_year'].notna()) &
                    (df['sig_year'] != df['year'])).sum()
    df['year_orig'] = df['year']
    df['year'] = df['sig_year'].fillna(df['year'])
    df = df.dropna(subset=['year'])
    df['year'] = df['year'].astype(int)
    print(f'  Year fix: {n_sig:,}/{n_tender_rows:,} have signature date, '
          f'{n_reassigned:,} reassigned to signature year')

    # --- FIX 2: bids_count == 0 → NaN (impossible: winner must bid) ---
    n_zero_bids = (df['bids_count'] == 0).sum()
    df.loc[df['bids_count'] == 0, 'bids_count'] = np.nan
    print(f'  bids_count == 0 → NaN: {n_zero_bids}')

    # Compute bid discount where estimated price is available
    df['discount'] = np.where(
        df['lot_estimated_price'].notna() & (df['lot_estimated_price'] > 0),
        (df['lot_estimated_price'] - df['pp_sales']) / df['lot_estimated_price'],
        np.nan
    )
    # Cap discount to [-1, 1] range (winsorize outliers)
    df['discount'] = df['discount'].clip(-1, 1)

    # Single-bidder indicator
    df['single_bid'] = (df['bids_count'] == 1).astype(float)
    df.loc[df['bids_count'].isna(), 'single_bid'] = np.nan

    # --- FIX 3: n_contracts counts only rows with non-null bid price ---
    df['has_price'] = df['pp_sales'].notna().astype(int)

    # Collapse to firm-year with competition measures
    agg = df.groupby(['id', 'year']).agg(
        pp_sales=('pp_sales', 'sum'),
        n_contracts=('has_price', 'sum'),  # FIX: count priced contracts only
        avg_bids=('bids_count', 'mean'),
        max_bids=('bids_count', 'max'),
        single_bid_share=('single_bid', 'mean'),
        avg_discount=('discount', 'mean'),
    ).reset_index()

    # HHI of contract values within firm-year (revenue concentration)
    def _hhi(group):
        vals = group['pp_sales'].dropna()
        total = vals.sum()
        if total <= 0 or len(vals) <= 1:
            return np.nan
        shares = vals / total
        return (shares ** 2).sum()

    hhi = df.groupby(['id', 'year']).apply(_hhi, include_groups=False).reset_index()
    hhi.columns = ['id', 'year', 'hhi_revenue']
    agg = agg.merge(hhi, on=['id', 'year'], how='left')

    print(f'  Tenders: {n_tender_rows:,} contract rows → {len(agg):,} firm-year obs')
    print(f'    avg_bids coverage:    {agg["avg_bids"].notna().mean():.0%}')
    print(f'    avg_discount coverage: {agg["avg_discount"].notna().mean():.0%}')
    print(f'    single_bid_share coverage: {agg["single_bid_share"].notna().mean():.0%}')
    return agg


# ========================================================================== #
#  Step 2: Deduplicate
# ========================================================================== #

def deduplicate_financial(df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate financial data on (id, year).

    Strategy: for conflicts (different values for same id-year),
    keep the row with sales closest to the group median.
    This avoids picking consolidated vs. unconsolidated extremes.
    """
    n_before = len(df)

    # Drop exact duplicates first
    df = df.drop_duplicates()
    n_exact = n_before - len(df)

    # Find remaining id-year duplicates
    dup_mask = df.duplicated(subset=['id', 'year'], keep=False)
    n_conflict_rows = dup_mask.sum()

    if n_conflict_rows > 0:
        clean = df[~dup_mask].copy()
        dups = df[dup_mask].copy()

        # For each group, keep the row with sales closest to group median
        resolved = []
        for (id_, yr), grp in dups.groupby(['id', 'year']):
            if grp['sales'].notna().any():
                med = grp['sales'].median()
                idx = (grp['sales'] - med).abs().idxmin()
                resolved.append(grp.loc[[idx]])
            else:
                # No sales info; just keep first
                resolved.append(grp.head(1))

        resolved = pd.concat(resolved, ignore_index=True)
        df = pd.concat([clean, resolved], ignore_index=True)

    n_conflict_groups = n_conflict_rows // 2  # approximate
    print(f'  Financial dedup: {n_exact:,} exact dups, '
          f'~{n_conflict_groups} conflicting groups resolved (median-sales rule)')
    print(f'  {n_before:,} → {len(df):,} rows')
    return df


def deduplicate_ratios(df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate ratios data on (id, year).

    Strategy: for conflicts, keep the row with cm closest to 40
    (typical contribution margin on 0-100 scale for construction).
    """
    n_before = len(df)
    df = df.drop_duplicates()
    n_exact = n_before - len(df)

    dup_mask = df.duplicated(subset=['id', 'year'], keep=False)
    n_conflict_rows = dup_mask.sum()

    if n_conflict_rows > 0:
        clean = df[~dup_mask].copy()
        dups = df[dup_mask].copy()

        resolved = []
        for (id_, yr), grp in dups.groupby(['id', 'year']):
            if grp['cm'].notna().any():
                # Keep row with cm closest to typical value (~40%)
                idx = (grp['cm'] - 40).abs().idxmin()
                resolved.append(grp.loc[[idx]])
            else:
                resolved.append(grp.head(1))

        resolved = pd.concat(resolved, ignore_index=True)
        df = pd.concat([clean, resolved], ignore_index=True)

    n_conflict_groups = n_conflict_rows // 2
    print(f'  Ratios dedup: {n_exact:,} exact dups, '
          f'~{n_conflict_groups} conflicting groups resolved (median-cm rule)')
    print(f'  {n_before:,} → {len(df):,} rows')
    return df


# ========================================================================== #
#  Step 3: Merge (1:1)
# ========================================================================== #

def merge_data(financial: pd.DataFrame, ratios: pd.DataFrame,
               firm_info: pd.DataFrame, empl_panel: pd.DataFrame,
               deflators: pd.DataFrame,
               tenders: pd.DataFrame) -> pd.DataFrame:
    """Merge all datasets using proper 1:1 and m:1 merges."""

    # Verify uniqueness before merging
    assert not financial.duplicated(subset=['id', 'year']).any(), \
        'Financial still has duplicates!'
    assert not ratios.duplicated(subset=['id', 'year']).any(), \
        'Ratios still has duplicates!'

    # 1:1 merge financial + ratios on (id, year)
    df = financial.merge(ratios, on=['id', 'year'], how='outer',
                         indicator='_merge_fr')
    n_both = (df['_merge_fr'] == 'both').sum()
    n_fin_only = (df['_merge_fr'] == 'left_only').sum()
    n_rat_only = (df['_merge_fr'] == 'right_only').sum()
    print(f'  Financial × Ratios (1:1):')
    print(f'    Both: {n_both:,}, Financial only: {n_fin_only:,}, '
          f'Ratios only: {n_rat_only:,}')
    df = df.drop(columns=['_merge_fr'])

    # m:1 merge with firm-level selections on id
    df = df.merge(firm_info, on='id', how='left')
    print(f'  + Firm info (m:1): {len(df):,} rows')

    # 1:1 merge with employment panel on (id, year)
    if len(empl_panel) > 0:
        df = df.merge(empl_panel, on=['id', 'year'], how='left')
        n_empl = df['empl_mid'].notna().sum()
        print(f'  + Employment panel (1:1): {n_empl:,} matched '
              f'({n_empl / len(df):.0%})')

    # Verify no duplicates introduced
    n_dups_pre_defl = df.duplicated(subset=['id', 'year']).sum()
    if n_dups_pre_defl > 0:
        print(f'  WARNING: {n_dups_pre_defl} id-year duplicates after selections merge')
        df = df.drop_duplicates(subset=['id', 'year'], keep='first')

    # Generate nace2 before deflator merge
    df['nace2'] = (df['nace'] // 10000).astype('Int64')

    # m:1 merge with deflators on (year, nace2)
    n_pre = len(df)
    df = df.merge(deflators, on=['year', 'nace2'], how='left',
                  indicator='_merge_defl')
    n_matched = (df['_merge_defl'] == 'both').sum()
    n_unmatched = (df['_merge_defl'] == 'left_only').sum()
    print(f'  + Deflators (m:1): {n_matched:,} matched, {n_unmatched:,} unmatched')
    df = df.drop(columns=['_merge_defl'])

    # Check for duplicates after deflator merge
    n_dups_post = df.duplicated(subset=['id', 'year']).sum()
    if n_dups_post > 0:
        print(f'  WARNING: {n_dups_post} duplicates after deflator merge — dropping')
        df = df.drop_duplicates(subset=['id', 'year'], keep='first')
    else:
        print(f'  No duplicates after deflator merge (FIXED vs original: 2,348)')

    # 1:1 merge with tenders (now includes competition measures)
    df = df.merge(tenders, on=['id', 'year'], how='left')
    df['pp_sales'] = df['pp_sales'].fillna(0)
    df['pp_dummy'] = (df['pp_sales'] > 0).astype(int)
    df['pp_share'] = df['pp_sales'] / df['sales']
    df.loc[df['pp_share'] > 1, 'pp_share'] = 1
    # Fill competition measures with 0/NaN for non-PP firms
    df['n_contracts'] = df['n_contracts'].fillna(0).astype(int)
    n_with_comp = df.loc[df['pp_dummy'] == 1, 'avg_bids'].notna().sum()
    n_pp = (df['pp_dummy'] == 1).sum()
    print(f'  + Tenders (1:1): {len(df):,} rows')
    print(f'    Competition measures for {n_with_comp:,}/{n_pp:,} PP obs '
          f'({n_with_comp / max(n_pp, 1):.0%})')

    return df


# ========================================================================== #
#  Step 4: Apply cleaning (replicating Create_Data.do)
# ========================================================================== #

def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    """Apply same cleaning rules as Create_Data.do."""

    header('Step 4: Cleaning (replicating Create_Data.do)')
    n0 = len(df)

    # Correct sales
    n_neg_sales = (df['sales'] < 0).sum()
    df.loc[df['sales'] < 0, 'sales'] = np.nan
    print(f'  Sales < 0 → missing: {n_neg_sales}')

    # FIX 5: costs == 0 with sales > 0 → NaN (data error, not dormant firm)
    n_zero_costs = ((df['costs'] == 0) & (df['sales'] > 0)).sum()
    df.loc[(df['costs'] == 0) & (df['sales'] > 0), 'costs'] = np.nan
    print(f'  Costs == 0 with sales > 0 → missing: {n_zero_costs}')

    # Correct ws (wage share) — THE 100x CORRECTION
    n_ws_corr = ((df['ws'] < 1) & (df['ws'] > 0)).sum()
    df.loc[(df['ws'] < 1) & (df['ws'] > 0), 'ws'] = \
        df.loc[(df['ws'] < 1) & (df['ws'] > 0), 'ws'] * 100
    n_ws_neg = (df['ws'] < 0).sum()
    df.loc[df['ws'] < 0, 'ws'] = np.nan
    n_ws_high = (df['ws'] > 100).sum()
    df.loc[df['ws'] > 100, 'ws'] = np.nan
    df['ws'] = df['ws'] / 100  # convert to share
    print(f'  ws: {n_ws_corr} corrected (×100), {n_ws_neg} neg→NaN, {n_ws_high} >100→NaN')

    # Correct wva
    n_wva_corr = ((df['wva'] < 1) & (df['wva'] > 0)).sum()
    df.loc[(df['wva'] < 1) & (df['wva'] > 0), 'wva'] = \
        df.loc[(df['wva'] < 1) & (df['wva'] > 0), 'wva'] * 100
    df.loc[df['wva'] < 0, 'wva'] = np.nan
    df.loc[df['wva'] > 500, 'wva'] = np.nan
    df['wva'] = df['wva'] / 100
    print(f'  wva: {n_wva_corr} corrected (×100)')

    # Correct cm — NO upper bound in MSc (only cm < 0)
    n_cm_corr = ((df['cm'] < 1) & (df['cm'] > 0)).sum()
    df.loc[(df['cm'] < 1) & (df['cm'] > 0), 'cm'] = \
        df.loc[(df['cm'] < 1) & (df['cm'] > 0), 'cm'] * 100
    n_cm_neg = (df['cm'] < 0).sum()
    df.loc[df['cm'] < 0, 'cm'] = np.nan
    # FIX 4: cm > 90 implausible for construction (>90% contribution margin)
    n_cm_high = (df['cm'] > 90).sum()
    df.loc[df['cm'] > 90, 'cm'] = np.nan
    df['cm'] = df['cm'] / 100
    print(f'  cm: {n_cm_corr} corrected (×100), {n_cm_neg} neg→NaN, '
          f'{n_cm_high} >90→NaN')

    # cs = costs / sales
    df['cs'] = df['costs'] / df['sales']
    df.loc[df['cs'] < 0, 'cs'] = np.nan
    df.loc[df['cs'] > 10, 'cs'] = np.nan

    # iis = cs - ws (intermediate inputs share)
    df['iis'] = df['cs'] - df['ws']
    df.loc[df['iis'] < 0, 'iis'] = np.nan

    # COGSS = 1 - cm (COGS share)
    df['COGSS'] = 1 - df['cm']
    df.loc[df['COGSS'] < 0, 'COGSS'] = np.nan
    df.loc[df['COGSS'] > 10, 'COGSS'] = np.nan

    # Correct lp
    df.loc[df['lp'] > 1e8, 'lp'] = df.loc[df['lp'] > 1e8, 'lp'] / 1000

    # Generate variables
    df['GO'] = df['sales']
    df['W'] = df['ws'] * df['sales']
    df['II'] = df['iis'] * df['sales']
    df['COGS'] = df['COGSS'] * df['sales']
    df['VA'] = df['GO'] - df['II']
    df['L'] = np.where((df['VA'] / df['lp']) > 0, df['VA'] / df['lp'], np.nan)
    df['K'] = df['assets']

    # Correct negatives
    for var in ['GO', 'COGS', 'II', 'W', 'K']:
        n_neg = (df[var] < 0).sum()
        df.loc[df[var] < 0, var] = np.nan
        if n_neg > 0:
            print(f'  {var} < 0 → missing: {n_neg}')

    # Deflate
    df['rGO'] = df['GO'] / df['deflatorprdp']
    df['rVA'] = df['VA'] / df['deflatorvalp']
    df['rII'] = df['II'] / df['deflatorintp']
    df['rW'] = df['W'] / df['deflatorcpi']
    df['rK'] = df['K'] / df['deflatorgfcp']
    df['rCOGS'] = df['COGS'] / df['deflatorintp']

    # Log variables
    for (raw, log) in [('rGO', 'go'), ('rW', 'w'), ('rII', 'ii'),
                        ('rVA', 'va'), ('rK', 'k'), ('rCOGS', 'cogs')]:
        df[log] = np.log(df[raw])
    df['l'] = np.log(df['L'])

    # ---- NEW: DGM-style additional variables ----

    # Overhead/services: O = II - COGS (non-production non-wage costs)
    # DGM analog: o = ln(autachaR) = "other purchases"
    df['O'] = df['II'] - df['COGS']
    df.loc[df['O'] <= 0, 'O'] = np.nan
    df['rO'] = df['O'] / df['deflatorintp']
    df['o'] = np.log(df['rO'])
    n_o = df['o'].notna().sum()
    print(f'  Overhead (o = ln(II - COGS)): {n_o:,} non-missing '
          f'({n_o / len(df):.0%})')

    # Employment from panel (empl_mid): log employment
    if 'empl_mid' in df.columns:
        df['le'] = np.log(df['empl_mid'])
        n_le = df['le'].notna().sum()
        print(f'  Log employment (le): {n_le:,} non-missing '
              f'({n_le / len(df):.0%})')

    # Market share: firm revenue / nace2-year total (DGM: ms5d)
    df['mktshare'] = np.nan
    valid = df['rGO'].notna() & df['nace2'].notna()
    total_rev = df.loc[valid].groupby(['nace2', 'year'])['rGO'].transform('sum')
    df.loc[valid, 'mktshare'] = df.loc[valid, 'rGO'] / total_rev
    n_ms = df['mktshare'].notna().sum()
    print(f'  Market share (mktshare): {n_ms:,} non-missing '
          f'({n_ms / len(df):.0%})')

    # Sales share within nace2-year (DGM first-stage control)
    # This is the same as mktshare but named to match acf_estimator.py convention
    df['salessharefirm'] = df['mktshare']

    # Panel structure: require consecutive years
    df = df.sort_values(['id', 'year'])
    df['lag_year'] = df.groupby('id')['year'].shift(1)
    df['lag_id'] = df.groupby('id')['id'].shift(1)
    df['lagExists'] = ((df['year'] == df['lag_year'] + 1) &
                        (df['id'] == df['lag_id'])).astype(int)
    n_pre_lag = len(df)
    df = df[df['lagExists'] == 1].copy()
    print(f'  Consecutive-year filter: {n_pre_lag:,} → {len(df):,}')

    # Entry/exit
    df['entry'] = (df.groupby('id').cumcount() == 0).astype(int)
    max_year = df.groupby('id')['year'].transform('max')
    df['exit'] = ((df['year'] == max_year) & (df['year'] < 2021)).astype(int)

    # Drop temp columns
    df = df.drop(columns=['lag_year', 'lag_id', 'lagExists'], errors='ignore')

    # ---- Procurement indicators ----
    # Three distinct roles:
    #
    # A. MARKOV STATE (for PF estimation: ω_t = g(ω_{t-1}, state_{t-1}))
    #    The firm's information set at t-1 when choosing inputs. Should
    #    capture accumulated procurement experience/capacity, not a single
    #    contract outcome. Stock variables work best here.
    #
    # B. TREATMENT INDICATOR (for markup comparison: μ(PP=1) vs μ(PP=0))
    #    Should separate the *decision* to participate (endogenous) from
    #    the *outcome* of winning (more plausibly exogenous). Revenue
    #    recognition lag means the production we observe at t reflects
    #    contracts won at t-1 or t-2.
    #
    # C. COMPETITION ENVIRONMENT (for ADL sufficient statistic)
    #    avg_bids, single_bid_share — already constructed in load_tenders().

    df = df.sort_values(['id', 'year'])

    # --- A. Markov state variables (lagged stock = firm's info at t-1) ---

    # Lags of pp_sales for rolling windows
    df['pp_sales_L1'] = df.groupby('id')['pp_sales'].shift(1).fillna(0)
    df['pp_sales_L2'] = df.groupby('id')['pp_sales'].shift(2).fillna(0)

    # Rolling procurement revenue stock (2-year and 3-year windows)
    # These capture the firm's "procurement pipeline" — ongoing contracts.
    df['pp_stock_2y'] = df['pp_sales'] + df['pp_sales_L1']
    df['pp_stock_3y'] = df['pp_stock_2y'] + df['pp_sales_L2']

    # Cumulative procurement revenue since start of panel
    # Captures total procurement EXPERIENCE — doesn't depreciate.
    # Analogous to human capital: organizational capacity, compliance
    # infrastructure, network effects with contracting authorities.
    # Better Markov state than rolling window because a firm with 50M
    # cumulative PP has different capabilities than one with 2M, even
    # if both had pp_dummy=1 last year.
    df['pp_cumul_revenue'] = df.groupby('id')['pp_sales'].cumsum()

    # Log versions (0 → enters as 0 via indicator interaction)
    for var in ['pp_stock_2y', 'pp_stock_3y', 'pp_cumul_revenue']:
        pos_var = f'{var}_pos'
        log_var = f'l{var}'
        df[pos_var] = (df[var] > 0).astype(int)
        df[log_var] = np.where(df[var] > 0, np.log(df[var]), 0.0)

    # Active procurement dummy (rolling window definitions)
    df['pp_active_2y'] = (df['pp_stock_2y'] > 0).astype(int)
    df['pp_active_3y'] = (df['pp_stock_3y'] > 0).astype(int)

    # --- B. Treatment indicators ---

    # Lagged dummies (binary: won anything in year t-k?)
    df['pp_dummy_L1'] = df.groupby('id')['pp_dummy'].shift(1).fillna(0).astype(int)
    df['pp_dummy_L2'] = df.groupby('id')['pp_dummy'].shift(2).fillna(0).astype(int)

    # Ever-PP window dummies (won anything in rolling window?)
    df['pp_ever_2y'] = ((df['pp_dummy'] == 1) |
                         (df['pp_dummy_L1'] == 1)).astype(int)
    df['pp_ever_3y'] = ((df['pp_dummy'] == 1) |
                         (df['pp_dummy_L1'] == 1) |
                         (df['pp_dummy_L2'] == 1)).astype(int)

    # First procurement entry (for event-study / staggered DiD)
    first_pp = df.loc[df['pp_dummy'] == 1].groupby('id')['year'].min()
    first_pp.name = 'pp_entry_year'
    df = df.merge(first_pp, on='id', how='left')
    df['pp_years_since_entry'] = df['year'] - df['pp_entry_year']
    # Never-PP firms get NaN for entry year and years_since_entry

    # Procurement intensity (share of revenue, lagged for treatment timing)
    df['pp_share_L1'] = df.groupby('id')['pp_share'].shift(1).fillna(0)

    # Cumulative procurement experience (total contracts won up to t)
    df['pp_cumul_contracts'] = df.groupby('id')['n_contracts'].cumsum()

    n_pp = (df['pp_dummy'] == 1).sum()
    n_active2 = (df['pp_active_2y'] == 1).sum()
    n_active3 = (df['pp_active_3y'] == 1).sum()
    n_ever_entered = df['pp_entry_year'].notna().sum()
    print(f'  Procurement indicators:')
    print(f'    pp_dummy=1 (current year):  {n_pp:,}')
    print(f'    pp_active_2y (t or t-1):    {n_active2:,}')
    print(f'    pp_active_3y (t, t-1, t-2): {n_active3:,}')
    print(f'    Ever entered procurement:   {n_ever_entered:,}')

    print(f'  Final magnus equivalent: {len(df):,} obs, {df["id"].nunique():,} firms')
    return df


def winsorize_and_trim(df: pd.DataFrame, by: str = 'nace2',
                       pct: float = 0.02) -> pd.DataFrame:
    """DGM-style winsorization + outlier trimming on new variables.

    Two layers:
      1. Winsorize log variables at `pct` tails by `by`-group (DGM 2026: 2%)
      2. Trim implausible observations (O/GO > 0.90, etc.)

    Parameters
    ----------
    df : DataFrame
        Must contain log variables and nace2.
    by : str
        Group variable for winsorization (default: 'nace2').
    pct : float
        Fraction to winsorize at each tail (default: 0.02 = 2%).

    Returns
    -------
    DataFrame with winsorized log vars (originals preserved as {var}_raw).
    """
    header('Step 5a: Winsorization & outlier trimming')
    n0 = len(df)

    # --- DGM-style winsorization on log variables ---
    log_vars = ['go', 'k', 'cogs', 'ii', 'o', 'le', 'va', 'w', 'l']
    log_vars = [v for v in log_vars if v in df.columns]

    n_winsorized = {}
    for var in log_vars:
        # Save raw version
        df[f'{var}_raw'] = df[var].copy()
        n_wins = 0
        for grp in df[by].dropna().unique():
            mask = (df[by] == grp) & df[var].notna()
            vals = df.loc[mask, var]
            if len(vals) < 20:
                continue
            lo = vals.quantile(pct)
            hi = vals.quantile(1 - pct)
            clip_lo = mask & (df[var] < lo)
            clip_hi = mask & (df[var] > hi)
            n_wins += clip_lo.sum() + clip_hi.sum()
            df.loc[clip_lo, var] = lo
            df.loc[clip_hi, var] = hi
        n_winsorized[var] = n_wins

    total_wins = sum(n_winsorized.values())
    print(f'  Winsorized {total_wins:,} values ({pct:.0%} tails by {by}):')
    for var, n in n_winsorized.items():
        if n > 0:
            print(f'    {var:6s}: {n:,}')

    # --- Trim implausible observations ---
    # O/GO > 0.90: overhead exceeding 90% of revenue
    if 'O' in df.columns and 'GO' in df.columns:
        ogo = df['O'] / df['GO']
        bad_ogo = (ogo > 0.90) & ogo.notna()
        n_bad_ogo = bad_ogo.sum()
        df = df[~bad_ogo].copy()
        print(f'  Trimmed O/GO > 0.90: {n_bad_ogo}')

    # n_contracts at p99 cap (keep obs, just cap the value)
    if 'n_contracts' in df.columns:
        p99 = df.loc[df['n_contracts'] > 0, 'n_contracts'].quantile(0.99)
        n_capped = (df['n_contracts'] > p99).sum()
        df.loc[df['n_contracts'] > p99, 'n_contracts'] = int(p99)
        print(f'  Capped n_contracts > {int(p99)} (p99): {n_capped}')

    # avg_discount < -0.50: likely data errors
    if 'avg_discount' in df.columns:
        bad_disc = (df['avg_discount'] < -0.50) & df['avg_discount'].notna()
        n_bad_disc = bad_disc.sum()
        df.loc[bad_disc, 'avg_discount'] = np.nan
        print(f'  avg_discount < -0.50 → NaN: {n_bad_disc}')

    print(f'  After trimming: {n0:,} → {len(df):,}')
    return df


def create_estimation_sample(df: pd.DataFrame,
                             winsorize: bool = True) -> pd.DataFrame:
    """Apply Create_Data.do post-magnus steps to create estimation sample.

    Parameters
    ----------
    winsorize : bool
        If True (default), apply DGM-style winsorization after MSc thesis trimming.
        If False, skip winsorization (for sensitivity checks on raw data).
    """

    header('Step 5: Estimation sample (post-magnus processing)')

    # Drop missing go, k, cogs
    n0 = len(df)
    df = df.dropna(subset=['cogs', 'k', 'go']).copy()
    print(f'  Drop missing (go, k, cogs): {n0:,} → {len(df):,}')

    # Drop year < 2005, require lag again
    df = df[df['year'] >= 2005].copy()

    # Panel setup + consecutive year requirement
    df = df.sort_values(['id', 'year'])
    df['lag_year'] = df.groupby('id')['year'].shift(1)
    df['lag_id'] = df.groupby('id')['id'].shift(1)
    df['lgExist'] = ((df['year'] == df['lag_year'] + 1) &
                      (df['id'] == df['lag_id']))
    df = df[df['lgExist'] == True].copy()
    df = df.drop(columns=['lag_year', 'lag_id', 'lgExist'], errors='ignore')
    print(f'  After year>=2005 + consecutive-year: {len(df):,}')

    # Trimming: 1%/99% on sales-cogs ratio
    df['s_g'] = df['rGO'] / df['rCOGS']
    for yr in df['year'].unique():
        mask = df['year'] == yr
        p1 = df.loc[mask, 's_g'].quantile(0.01)
        p99 = df.loc[mask, 's_g'].quantile(0.99)
        df.loc[mask & ((df['s_g'] <= p1) | (df['s_g'] >= p99)), '_trim_sg'] = 1
    df['_trim_sg'] = df['_trim_sg'].fillna(0)
    n_trim_sg = int(df['_trim_sg'].sum())
    df = df[df['_trim_sg'] == 0].copy()
    df = df.drop(columns=['_trim_sg'])
    print(f'  Trim s_g 1/99 pctile: dropped {n_trim_sg}')

    # Trimming: 1%/99% on cost share
    df['costshare'] = df['rCOGS'] / (df['rCOGS'] + df['rK'])
    df = df[df['costshare'].notna() & (df['costshare'] > 0)].copy()
    for yr in df['year'].unique():
        mask = df['year'] == yr
        p1 = df.loc[mask, 'costshare'].quantile(0.01)
        p99 = df.loc[mask, 'costshare'].quantile(0.99)
        df.loc[mask & ((df['costshare'] >= p99) | (df['costshare'] <= p1)),
               '_trim_cs'] = 1
    df['_trim_cs'] = df['_trim_cs'].fillna(0)
    n_trim_cs = int(df['_trim_cs'].sum())
    df = df[df['_trim_cs'] == 0].copy()
    df = df.drop(columns=['_trim_cs', 's_g', 'costshare'])
    print(f'  Trim costshare 1/99 pctile: dropped {n_trim_cs}')

    # DGM-style winsorization (optional)
    if winsorize:
        df = winsorize_and_trim(df)

    print(f'  Final estimation sample: {len(df):,} obs, {df["id"].nunique():,} firms')
    return df


# ========================================================================== #
#  Step 6: Compare with original
# ========================================================================== #

def compare_with_original(new_data: pd.DataFrame) -> None:
    """Compare rebuilt data against original data.dta."""

    header('Step 6: Comparison with original data.dta')

    try:
        orig = pd.read_stata(str(OUT_DIR / 'data.dta'))
    except Exception as e:
        print(f'  Could not load original: {e}')
        return

    print(f'  Original: {len(orig):,} obs, {orig["id"].nunique():,} firms')
    print(f'  Rebuilt:  {len(new_data):,} obs, {new_data["id"].nunique():,} firms')

    # Match on (id, year)
    merged = orig[['id', 'year', 'go', 'k', 'cogs', 'pp_dummy']].merge(
        new_data[['id', 'year', 'go', 'k', 'cogs', 'pp_dummy']],
        on=['id', 'year'], how='outer', suffixes=('_orig', '_new'),
        indicator=True
    )

    n_both = (merged['_merge'] == 'both').sum()
    n_orig_only = (merged['_merge'] == 'left_only').sum()
    n_new_only = (merged['_merge'] == 'right_only').sum()
    print(f'\n  Overlap:')
    print(f'    Both:          {n_both:,}')
    print(f'    Original only: {n_orig_only:,}')
    print(f'    Rebuilt only:  {n_new_only:,}')

    # For matched obs, compare key variables
    both = merged[merged['_merge'] == 'both'].copy()
    if len(both) > 0:
        print(f'\n  Variable differences (matched obs):')
        for var in ['go', 'k', 'cogs']:
            diff = (both[f'{var}_orig'] - both[f'{var}_new']).abs()
            n_diff = (diff > 0.001).sum()
            if n_diff > 0:
                max_diff = diff.max()
                mean_diff = diff[diff > 0.001].mean()
                print(f'    {var}: {n_diff:,} obs differ '
                      f'(max={max_diff:.4f}, mean={mean_diff:.4f})')
            else:
                print(f'    {var}: all identical')

    # Diagnose original-only and rebuilt-only obs
    if n_orig_only > 0:
        orig_only = merged[merged['_merge'] == 'left_only']
        print(f'\n  Original-only obs ({n_orig_only}):')
        print(f'    These exist in original but not in rebuilt data.')
        print(f'    Likely cause: m:m merge created spurious obs from')
        print(f'    cross-product of duplicate (id,year) rows.')
        print(f'    Year range: {orig_only["year"].min()}-{orig_only["year"].max()}')
        print(f'    Unique firms: {orig_only["id"].nunique()}')

    if n_new_only > 0:
        new_only = merged[merged['_merge'] == 'right_only']
        print(f'\n  Rebuilt-only obs ({n_new_only}):')
        print(f'    These exist in rebuilt but not in original.')
        print(f'    Likely cause: deflator merge duplicates in original')
        print(f'    dropped the "wrong" row, removing valid obs.')
        print(f'    Year range: {new_only["year"].min()}-{new_only["year"].max()}')
        print(f'    Unique firms: {new_only["id"].nunique()}')

    # 11 obs with different values — these are the conflict resolutions
    if len(both) > 0:
        diff_mask = (both['go_orig'] - both['go_new']).abs() > 0.001
        if diff_mask.sum() > 0:
            diff_obs = both[diff_mask][['id', 'year',
                                        'go_orig', 'go_new',
                                        'cogs_orig', 'cogs_new']].head(5)
            print(f'\n  Sample of obs with different values (conflict resolutions):')
            print(diff_obs.to_string(index=False))


# ========================================================================== #
#  Main
# ========================================================================== #

def main():
    header('REBUILD DATA FROM RAW CSVs (FIXING m:m MERGE)')

    # Step 1: Load
    print('Step 1: Loading raw CSVs...')
    financial = load_financial()
    ratios = load_ratios()
    firm_info, empl_panel = load_selections()
    deflators = load_deflators()
    tenders = load_tenders()

    # Step 2: Deduplicate
    header('Step 2: Deduplication')
    financial = deduplicate_financial(financial)
    ratios = deduplicate_ratios(ratios)

    # Step 3: Merge
    header('Step 3: Proper 1:1 merge')
    merged = merge_data(financial, ratios, firm_info, empl_panel,
                        deflators, tenders)

    # Step 4: Clean
    magnus = clean_data(merged)

    # Step 5: Estimation sample (winsorized = default)
    est = create_estimation_sample(magnus, winsorize=True)

    # Step 5b: Also save raw (no winsorization) for sensitivity
    est_raw = create_estimation_sample(magnus.copy(), winsorize=False)

    # Step 6: Compare
    compare_with_original(est)

    # Save — replace inf/-inf with NaN for Stata compatibility
    header('Saving')
    magnus_out = OUT_DIR / 'magnus_rebuilt.dta'
    data_out = OUT_DIR / 'data_rebuilt.dta'
    data_raw_out = OUT_DIR / 'data_rebuilt_raw.dta'

    for d, path in [(magnus, magnus_out), (est, data_out),
                     (est_raw, data_raw_out)]:
        # Replace inf with NaN for Stata
        d_save = d.copy()
        for col in d_save.select_dtypes(include=[np.floating]).columns:
            d_save[col] = d_save[col].replace([np.inf, -np.inf], np.nan)
        # Drop columns Stata can't handle (Int64 nullable)
        for col in d_save.columns:
            if d_save[col].dtype == 'Int64':
                d_save[col] = d_save[col].astype('float64')
        d_save.to_stata(str(path), write_index=False, version=118)
        print(f'  Saved: {path}')

    print(f'\n  To use rebuilt data in acf_estimator.py, change data_path to:')
    print(f'    Winsorized (default): {data_out}')
    print(f'    Raw (no winsorization): {data_raw_out}')

    # Summary of new variables
    header('New Variables Summary (estimation sample)')
    new_vars = {
        # DGM-style variables
        'o': 'Log overhead (ln(II-COGS)), DGM "services"',
        'le': 'Log employment (from category panel)',
        'mktshare': 'Market share (rev / nace2-year total)',
        # A. Markov state (firm info set at t-1)
        'pp_cumul_revenue': 'Cumulative PP revenue since panel start',
        'lpp_cumul_revenue': 'Log cumulative PP revenue (0 if never)',
        'pp_stock_3y': 'PP revenue stock (t + t-1 + t-2)',
        'pp_active_3y': 'Active PP in 3y window (dummy)',
        # B. Treatment indicators
        'pp_dummy_L1': 'Won PP contract at t-1 (dummy)',
        'pp_dummy_L2': 'Won PP contract at t-2 (dummy)',
        'pp_ever_2y': 'Won PP in t or t-1 (dummy)',
        'pp_ever_3y': 'Won PP in t, t-1, or t-2 (dummy)',
        'pp_entry_year': 'First year of PP participation',
        'pp_years_since_entry': 'Years since first PP contract',
        'pp_share_L1': 'Lagged PP share of revenue',
        'pp_cumul_contracts': 'Cumulative PP contracts won',
        # Competition measures (from tenders)
        'n_contracts': 'Number of PP contracts won (year t)',
        'avg_bids': 'Avg bidders per contract (competition)',
        'max_bids': 'Max bidders on any contract',
        'single_bid_share': 'Fraction of single-bidder contracts',
        'avg_discount': 'Avg (est_price - bid) / est_price',
        'hhi_revenue': 'HHI of contract values (concentration)',
        # Employment panel
        'empl_mid': 'Employment (midpoint of category)',
        # Firm characteristics
        'foreign': 'Foreign-controlled firm (dummy)',
        'legal_form': 'Legal form (string)',
        'inst_sector': 'Institutional sector (string)',
    }
    for var, desc in new_vars.items():
        if var not in est.columns:
            print(f'  {var:20s}: NOT IN DATA')
            continue
        col = est[var]
        if col.dtype in ['float64', 'int64', 'Int64']:
            nn = col.notna().sum()
            pp_mask = est['pp_dummy'] == 1
            nn_pp = col[pp_mask].notna().sum()
            if nn > 0:
                print(f'  {var:20s}: {desc}')
                print(f'    {"":20s}  N={nn:,} ({nn / len(est):.0%}), '
                      f'PP: {nn_pp:,}/{pp_mask.sum():,} ({nn_pp / max(pp_mask.sum(), 1):.0%}), '
                      f'mean={col.mean():.3f}, med={col.median():.3f}')
            else:
                print(f'  {var:20s}: {desc} — all missing')
        else:
            nn = col.notna().sum()
            print(f'  {var:20s}: {desc} — {nn:,} non-null, '
                  f'{col.nunique()} unique')


if __name__ == '__main__':
    main()
