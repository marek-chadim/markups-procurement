"""Build All-Industry Orbis + Datlab Estimation Panel.

Merges WRDS Orbis financials (10.8M CZ firm-year obs) with Datlab procurement
data (61K firms) to create a multi-industry panel for markup estimation.

Data sources:
  - Orbis Industry Global Financials (all tiers, downloaded Apr 3 2026)
  - Orbis Industry Classifications (NACE Rev.2)
  - Orbis Identifiers (BVDID → ICO crosswalk)
  - Datlab master_tender_analytics (Czech procurement register)
  - OECD/CZSO deflators by NACE 2-digit × year

Author: Marek Chadim (Yale, Tobin Center)
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings('ignore', category=pd.errors.DtypeWarning)

# ========================================================================== #
#  Paths
# ========================================================================== #

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

ORBIS_FINS = INPUT_DIR / 'orbis' / 'orbis_cz_all_tiers.csv'
ORBIS_NACE = INPUT_DIR / 'orbis' / 'orbis_cz_nace.csv'
ORBIS_IDS = INPUT_DIR / 'orbis' / 'orbis_cz_identifiers.csv'
TENDERS_CSV = INPUT_DIR / 'datlab' / 'master_tender_analytics.csv'
DEFLATORS_CSV = INPUT_DIR / 'magnusweb' / 'deflators.csv'

KEEP_VARS = [
    'empl', 'staf', 'opre', 'turn', 'cost', 'mate',
    'fias', 'tfas', 'toas', 'plbt', 'pl', 'depr', 'ebta', 'gros',
]


# ========================================================================== #
#  Step 1: Build Base Orbis Panel
# ========================================================================== #

def build_orbis_base() -> pd.DataFrame:
    """Load Orbis financials, merge NACE + ICO, clean and deduplicate."""
    print('=' * 60)
    print('Step 1: Building base Orbis panel')
    print('=' * 60)

    # --- Load NACE codes ---
    nace = pd.read_csv(ORBIS_NACE, usecols=['bvdid', 'naceccod2', 'nace2_main_section'])
    nace['naceccod2'] = nace['naceccod2'].astype(str)
    nace = nace[~nace['naceccod2'].isin(['nan', 'None', ''])]
    nace['nace4'] = nace['naceccod2'].str[:4]
    nace['nace2'] = pd.to_numeric(nace['naceccod2'].str[:2], errors='coerce')
    nace = nace.dropna(subset=['nace2'])
    nace['nace2'] = nace['nace2'].astype(int)
    nace = nace.drop_duplicates(subset=['bvdid'])
    print(f'  NACE codes: {len(nace):,} firms')

    # --- Load ICO crosswalk ---
    ids = pd.read_csv(ORBIS_IDS, usecols=['bvdid', 'natid_number', 'natid_label'])
    ids = ids[ids['natid_label'] == 'Trade register number']
    ids['ico'] = ids['natid_number'].astype(str).str.strip()
    ids = ids[['bvdid', 'ico']].drop_duplicates(subset=['bvdid'])
    print(f'  ICO crosswalk: {len(ids):,} firms')

    # --- Load financials in chunks (10.8M rows) ---
    print(f'  Loading financials from {ORBIS_FINS.name}...')
    use_cols = ['bvdid', 'closdate_year', 'conscode'] + KEEP_VARS
    # Only request columns that exist in the CSV
    with open(ORBIS_FINS) as f:
        available = f.readline().strip().split(',')
    use_cols = [c for c in use_cols if c in available]
    chunks = pd.read_csv(ORBIS_FINS, chunksize=1_000_000, usecols=use_cols)
    frames = []
    for i, chunk in enumerate(chunks):
        # Filter to 2006-2023
        chunk = chunk[(chunk['closdate_year'] >= 2006) &
                      (chunk['closdate_year'] <= 2023)]
        # Keep unconsolidated only (avoid double-counting)
        if 'conscode' in chunk.columns:
            chunk = chunk[chunk['conscode'].isin(['U1', 'U2', 'LF']) |
                         chunk['conscode'].isna()]
        # Need at least revenue or COGS
        chunk = chunk[chunk['opre'].notna() | chunk['cost'].notna()]
        frames.append(chunk)
        if (i + 1) % 5 == 0:
            print(f'    Processed {(i+1)}M rows...')

    df = pd.concat(frames, ignore_index=True)
    df['year'] = df['closdate_year'].astype(int)
    print(f'  Filtered financials: {len(df):,} rows, {df["bvdid"].nunique():,} firms')

    # --- Handle duplicates (multiple closing dates per year) ---
    # Keep the observation with more data (prefer non-null opre)
    df['has_opre'] = df['opre'].notna().astype(int)
    df = df.sort_values(['bvdid', 'year', 'has_opre'], ascending=[True, True, False])
    before = len(df)
    df = df.drop_duplicates(subset=['bvdid', 'year'], keep='first')
    print(f'  Deduplicated: {before:,} → {len(df):,} ({before - len(df):,} duplicates removed)')

    # --- Merge NACE ---
    df = df.merge(nace[['bvdid', 'nace2', 'nace4']], on='bvdid', how='left')
    n_nace = df['nace2'].notna().sum()
    print(f'  NACE matched: {n_nace:,}/{len(df):,} ({n_nace/len(df):.0%})')

    # --- Merge ICO ---
    df = df.merge(ids, on='bvdid', how='left')
    n_ico = df['ico'].notna().sum()
    print(f'  ICO matched: {n_ico:,}/{len(df):,} ({n_ico/len(df):.0%})')

    # Drop rows without NACE or ICO
    df = df.dropna(subset=['nace2', 'ico'])
    df['nace2'] = df['nace2'].astype(int)
    print(f'  After requiring NACE + ICO: {len(df):,} rows, {df["bvdid"].nunique():,} firms')

    return df


# ========================================================================== #
#  Step 2: Merge Datlab Procurement
# ========================================================================== #

def load_tenders() -> pd.DataFrame:
    """Aggregate tender data to firm-year level (reuses rebuild_data.py logic)."""
    print('\n' + '=' * 60)
    print('Step 2: Loading and aggregating Datlab tenders')
    print('=' * 60)

    df = pd.read_csv(TENDERS_CSV, low_memory=False)
    df = df.rename(columns={'bidder_id': 'ico', 'bid_final_price': 'pp_sales'})

    # Keep Czech firms only (8-digit ICO)
    df['ico'] = df['ico'].astype(str).str.strip()
    df = df[df['ico'].str.len() == 8]

    # Year from signature date or year field
    df['sig_year'] = pd.to_datetime(df['contract_signature_date'], errors='coerce').dt.year
    df['year'] = df['sig_year'].fillna(pd.to_numeric(df['year'], errors='coerce'))
    df = df.dropna(subset=['year'])
    df['year'] = df['year'].astype(int)

    # Fix bids_count == 0
    df.loc[df['bids_count'] == 0, 'bids_count'] = np.nan

    # Single-bidder indicator
    df['single_bid'] = (df['bids_count'] == 1).astype(float)
    df.loc[df['bids_count'].isna(), 'single_bid'] = np.nan

    # Discount
    df['discount'] = np.where(
        df['lot_estimated_price'].notna() & (df['lot_estimated_price'] > 0),
        (df['lot_estimated_price'] - df['pp_sales']) / df['lot_estimated_price'],
        np.nan
    ).clip(-1, 1) if 'lot_estimated_price' in df.columns else np.nan

    # Aggregate to firm-year
    df['has_price'] = df['pp_sales'].notna().astype(int)
    agg = df.groupby(['ico', 'year']).agg(
        pp_sales=('pp_sales', 'sum'),
        n_contracts=('has_price', 'sum'),
        avg_bids=('bids_count', 'mean'),
        max_bids=('bids_count', 'max'),
        single_bid_share=('single_bid', 'mean'),
    ).reset_index()

    print(f'  Tenders: {len(df):,} rows → {len(agg):,} firm-year obs')
    print(f'  Unique firms: {agg["ico"].nunique():,}')
    return agg


def merge_procurement(panel: pd.DataFrame, tenders: pd.DataFrame) -> pd.DataFrame:
    """Merge procurement indicators into the panel."""
    panel = panel.merge(tenders, on=['ico', 'year'], how='left')

    # Treatment variables
    panel['pp_dummy'] = (panel['pp_sales'].notna() & (panel['pp_sales'] > 0)).astype(int)
    panel['pp_sales'] = panel['pp_sales'].fillna(0)
    panel['n_contracts'] = panel['n_contracts'].fillna(0).astype(int)

    # Entry variables
    pp_firms = panel[panel['pp_dummy'] == 1].groupby('ico')['year'].min().reset_index()
    pp_firms.columns = ['ico', 'pp_entry_year']
    panel = panel.merge(pp_firms, on='ico', how='left')
    panel['pp_years_since_entry'] = np.where(
        panel['pp_entry_year'].notna(),
        panel['year'] - panel['pp_entry_year'],
        np.nan
    )

    n_pp = panel['pp_dummy'].sum()
    print(f'  Procurement obs: {n_pp:,}/{len(panel):,} ({n_pp/len(panel):.0%})')
    print(f'  Procurement firms: {panel.loc[panel["pp_dummy"]==1, "ico"].nunique():,}')
    return panel


# ========================================================================== #
#  Step 3: Construct Estimation Variables
# ========================================================================== #

def construct_variables(panel: pd.DataFrame) -> pd.DataFrame:
    """Deflate, take logs, winsorize, construct shares."""
    print('\n' + '=' * 60)
    print('Step 3: Constructing estimation variables')
    print('=' * 60)

    # --- Load deflators ---
    defl = pd.read_csv(DEFLATORS_CSV)
    defl = defl.rename(columns={
        'deflatorPRDP': 'defl_output',
        'deflatorINTP': 'defl_intermed',
        'deflatorGFCP': 'defl_capital',
        'deflatorCPI': 'defl_wages',
        'deflatorGDP': 'defl_gdp',
    })

    # Extend deflators to 2022-2023 using GDP deflator growth
    max_yr = defl['year'].max()
    if max_yr < 2023:
        for ext_yr in range(max_yr + 1, 2024):
            # Use last available year's growth rate
            last = defl[defl['year'] == max_yr].copy()
            prev = defl[defl['year'] == max_yr - 1].copy()
            if len(last) > 0 and len(prev) > 0:
                growth = last.set_index('nace2') / prev.set_index('nace2')
                ext = (last.set_index('nace2') * growth.values).reset_index()
                ext['year'] = ext_yr
                defl = pd.concat([defl, ext], ignore_index=True)
        print(f'  Deflators extended to {defl["year"].max()}')

    # Merge deflators
    panel = panel.merge(defl[['year', 'nace2', 'defl_output', 'defl_intermed',
                               'defl_capital', 'defl_wages', 'defl_gdp']],
                        on=['year', 'nace2'], how='left')

    # For firms without industry deflator, use GDP deflator
    for col in ['defl_output', 'defl_intermed', 'defl_capital', 'defl_wages']:
        panel[col] = panel[col].fillna(panel['defl_gdp'])

    n_defl = panel['defl_output'].notna().sum()
    print(f'  Deflator matched: {n_defl:,}/{len(panel):,} ({n_defl/len(panel):.0%})')

    # --- Deflate ---
    # Orbis reports in thousands of local currency (CZK)
    panel['rGO'] = panel['opre'] / panel['defl_output']
    panel['rCOGS'] = panel['cost'] / panel['defl_intermed']
    panel['rK'] = panel['fias'] / panel['defl_capital']
    panel['rW'] = panel['staf'] / panel['defl_wages']
    panel['rMATE'] = panel['mate'] / panel['defl_intermed']

    # --- Logs ---
    for var, src in [('go', 'rGO'), ('cogs', 'rCOGS'), ('k', 'rK'),
                     ('w', 'rW'), ('le', 'empl'), ('mate_log', 'rMATE')]:
        panel[var] = np.where(panel[src] > 0, np.log(panel[src]), np.nan)

    # --- Winsorize at 2/98 by nace2 ---
    log_vars = ['go', 'cogs', 'k', 'w', 'le']
    for var in log_vars:
        q02 = panel.groupby('nace2')[var].transform(lambda x: x.quantile(0.02))
        q98 = panel.groupby('nace2')[var].transform(lambda x: x.quantile(0.98))
        panel[var] = panel[var].clip(q02, q98)

    # --- Market share ---
    total_rev = panel.groupby(['nace2', 'year'])['rGO'].transform('sum')
    panel['mktshare'] = panel['rGO'] / total_rev

    # --- Require key variables for estimation ---
    has_go = panel['go'].notna()
    has_cogs = panel['cogs'].notna()
    has_k = panel['k'].notna()
    complete = has_go & has_cogs & has_k
    print(f'  Complete (go + cogs + k): {complete.sum():,}/{len(panel):,} ({complete.mean():.0%})')
    print(f'  With le (employees): {panel["le"].notna().sum():,}')
    print(f'  With w (wages): {panel["w"].notna().sum():,}')

    # --- Panel structure: sort and require lag ---
    panel = panel.sort_values(['ico', 'year'])
    panel['has_lag'] = panel.groupby('ico')['year'].diff() == 1
    estimation_sample = panel[complete].copy()

    print(f'\n  Estimation sample: {len(estimation_sample):,} rows, '
          f'{estimation_sample["ico"].nunique():,} firms')

    return panel, estimation_sample


# ========================================================================== #
#  Step 4: Validate Against MagnusWeb
# ========================================================================== #

def validate_magnusweb(panel: pd.DataFrame) -> None:
    """Compare Orbis financials to MagnusWeb for overlapping construction firms."""
    print('\n' + '=' * 60)
    print('Step 4: Validating against MagnusWeb')
    print('=' * 60)

    magnus_path = OUTPUT_DIR / 'data_rebuilt.dta'
    if not magnus_path.exists():
        print('  data_rebuilt.dta not found — skipping MagnusWeb validation')
        return
    magnus = pd.read_stata(magnus_path)
    magnus['ico'] = magnus['id'].astype(int).astype(str)
    magnus['year'] = pd.to_numeric(magnus['year'], errors='coerce').astype(int)

    # Merge on ico + year
    orbis_constr = panel[panel['nace2'].isin([41, 42, 43])].copy()
    merged = orbis_constr.merge(magnus[['ico', 'year', 'go', 'cogs', 'k', 'w']],
                                on=['ico', 'year'], how='inner',
                                suffixes=('_orbis', '_magnus'))

    if len(merged) == 0:
        print('  No overlapping firm-years found!')
        return

    print(f'  Overlapping firm-years: {len(merged):,}')

    for orbis_var, magnus_var, label in [
        ('go_orbis', 'go_magnus', 'GO (log output)'),
        ('cogs_orbis', 'cogs_magnus', 'COGS (log costs)'),
        ('k_orbis', 'k_magnus', 'K (log capital)'),
        ('w_orbis', 'w_magnus', 'W (log wages)'),
    ]:
        valid = merged[[orbis_var, magnus_var]].dropna()
        if len(valid) < 10:
            print(f'  {label}: too few obs ({len(valid)})')
            continue
        corr = valid[orbis_var].corr(valid[magnus_var])
        diff = (valid[orbis_var] - valid[magnus_var]).mean()
        print(f'  {label}: corr={corr:.3f}, mean diff={diff:.3f}, N={len(valid)}')


# ========================================================================== #
#  Step 5: Summary Statistics
# ========================================================================== #

def summary_stats(panel: pd.DataFrame, est: pd.DataFrame) -> pd.DataFrame:
    """Produce summary statistics by industry."""
    print('\n' + '=' * 60)
    print('Step 5: Summary Statistics')
    print('=' * 60)

    rows = []
    for nace in sorted(est['nace2'].unique()):
        sub = est[est['nace2'] == nace]
        rows.append({
            'nace2': int(nace),
            'firms': sub['ico'].nunique(),
            'obs': len(sub),
            'pp_rate': sub['pp_dummy'].mean(),
            'empl_coverage': sub['le'].notna().mean(),
            'mean_go': sub['go'].mean(),
            'mean_cogs': sub['cogs'].mean(),
            'mean_empl': np.exp(sub['le'].dropna()).mean() if sub['le'].notna().any() else np.nan,
        })

    stats = pd.DataFrame(rows).sort_values('obs', ascending=False)
    print(f'\nTop 20 industries by obs:')
    print(stats.head(20).to_string(index=False, float_format='{:.3f}'.format))

    # Totals
    print(f'\nTotal: {est["ico"].nunique():,} firms, {len(est):,} obs')
    print(f'Procurement rate: {est["pp_dummy"].mean():.1%}')
    print(f'Employee coverage: {est["le"].notna().mean():.0%}')

    return stats


# ========================================================================== #
#  Main
# ========================================================================== #

if __name__ == '__main__':
    # Step 1
    panel = build_orbis_base()

    # Step 2
    tenders = load_tenders()
    panel = merge_procurement(panel, tenders)

    # Step 3
    full_panel, est_sample = construct_variables(panel)

    # Step 4
    validate_magnusweb(full_panel)

    # Step 5
    stats = summary_stats(full_panel, est_sample)

    # Save
    print('\n' + '=' * 60)
    print('Saving outputs')
    print('=' * 60)

    # Full panel (all industries)
    out_all = OUTPUT_DIR / 'orbis_panel.dta'
    save_cols = ['ico', 'bvdid', 'year', 'nace2', 'nace4',
                 'empl', 'go', 'cogs', 'k', 'w', 'le', 'mate_log', 'mktshare',
                 'pp_dummy', 'pp_sales', 'n_contracts', 'avg_bids',
                 'single_bid_share', 'pp_entry_year', 'pp_years_since_entry',
                 'rGO', 'rCOGS', 'rK', 'rW']
    save_cols = [c for c in save_cols if c in est_sample.columns]
    est_sample[save_cols].to_stata(out_all, write_index=False, version=118)
    print(f'  {out_all}: {len(est_sample):,} rows')

    # Construction subset
    constr = est_sample[est_sample['nace2'].isin([41, 42, 43])]
    out_constr = OUTPUT_DIR / 'orbis_panel_construction.dta'
    constr[save_cols].to_stata(out_constr, write_index=False, version=118)
    print(f'  {out_constr}: {len(constr):,} rows')

    # Summary stats
    stats.to_csv(OUTPUT_DIR / 'orbis_summary_stats.csv', index=False)
    print(f'  {OUTPUT_DIR / "orbis_summary_stats.csv"}')

    print('\nDone.')
