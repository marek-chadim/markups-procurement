#!/usr/bin/env python3
"""Aggregate 0_raw/external/ into merge-ready CSVs for ABGRS instruments.

Reads:
    0_raw/external/cnb/cnb_fx_*.txt       -> daily FX (pipe-delimited)
    0_raw/external/eurostat/*.json        -> JSON-stat v2.0 time series (geo=CZ)
    0_raw/external/chmi/openmeteo_*.json  -> daily weather, 3 Czech stations

Writes (to 1_data/output/):
    external_panel_annual.csv   -> one row per year (2005-2023), ~20 cols
    external_panel_monthly.csv  -> monthly series (PPI industry only for now)

These CSVs are designed to merge into the firm-year panel on `year` and serve
as additional instruments for the ACF GMM moment set (ABGRS strong exclusion).

Usage:
    /opt/anaconda3/bin/python build_external_panel.py
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore', category=FutureWarning)

SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent           # 1_data
PROJECT_ROOT = MODULE_DIR.parent          # markups-procurement
RAW_DIR = PROJECT_ROOT / '0_raw' / 'external'
OUT_DIR = MODULE_DIR / 'output'
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================ #
#  JSON-stat v2.0 helpers
# ============================================================ #

def _js_strides(size: list[int]) -> list[int]:
    strides = [1] * len(size)
    for i in range(len(size) - 2, -1, -1):
        strides[i] = strides[i + 1] * size[i + 1]
    return strides


def _js_pick(dim_entry: dict, overrides: dict, name: str) -> int:
    cat = dim_entry['category']['index']
    if name in overrides:
        key = overrides[name]
        if key not in cat:
            raise KeyError(
                f"dim '{name}': '{key}' not in {list(cat.keys())}")
        return cat[key]
    if len(cat) == 1:
        return next(iter(cat.values()))
    raise KeyError(
        f"dim '{name}' has {len(cat)} categories {list(cat.keys())[:6]}... "
        f"— must be specified")


def js_extract_series(doc: dict, overrides: dict | None = None) -> dict:
    """Return {time_label: value} for one fixed non-time combination of dims."""
    overrides = overrides or {}
    ids = doc['id']
    size = doc['size']
    dims = doc['dimension']
    values = doc['value']
    strides = _js_strides(size)

    time_pos = ids.index('time')
    idx = [0] * len(ids)
    for i, name in enumerate(ids):
        if i == time_pos:
            continue
        idx[i] = _js_pick(dims[name], overrides, name)

    base = sum(idx[i] * strides[i] for i in range(len(ids)) if i != time_pos)

    time_cat = dims['time']['category']['index']
    time_by_pos = {pos: label for label, pos in time_cat.items()}

    out = {}
    for t in range(size[time_pos]):
        key = str(base + t * strides[time_pos])
        if key in values:
            out[time_by_pos[t]] = values[key]
    return out


# ============================================================ #
#  Source 1: CNB FX (pipe-delimited daily files)
# ============================================================ #

def parse_cnb_fx() -> pd.DataFrame:
    """CNB publishes daily FX with variable column counts — the file header
    can change mid-file when a new currency is added (e.g. CNY on 2005-04-01).
    We parse line by line and track the *current* header to handle this."""
    files = sorted((RAW_DIR / 'cnb').glob('cnb_fx_*.txt'))
    if not files:
        return pd.DataFrame()

    col_name = {
        '1 EUR': 'fx_eur', '1 USD': 'fx_usd',
        '1 CNY': 'fx_cny', '1 CHF': 'fx_chf',
    }
    target_cols = list(col_name.values())

    records = []
    for f in files:
        current_header = None
        for line in f.read_text().splitlines():
            parts = line.split('|')
            if parts[0] == 'Date':
                current_header = parts
                continue
            if current_header is None or len(parts) != len(current_header):
                continue
            try:
                d = pd.to_datetime(parts[0], format='%d.%m.%Y')
            except (ValueError, TypeError):
                continue
            row = {'Date': d}
            for src_key, dst_name in col_name.items():
                if src_key in current_header:
                    i = current_header.index(src_key)
                    try:
                        row[dst_name] = float(parts[i])
                    except ValueError:
                        row[dst_name] = np.nan
                else:
                    row[dst_name] = np.nan
            records.append(row)

    daily = pd.DataFrame(records)
    daily['year'] = daily['Date'].dt.year

    mean = daily.groupby('year')[target_cols].mean()
    std = daily.groupby('year')[target_cols].std()
    std.columns = [c + '_std' for c in std.columns]

    return mean.join(std)


# ============================================================ #
#  Source 2: Eurostat (JSON-stat v2.0)
# ============================================================ #

# (file_stem, overrides, output_col, frequency)
EUROSTAT_SERIES = [
    ('sts_inpp_a_ppi_annual',
     {'indic_bt': 'PRC_PRR', 'nace_r2': 'B-E36', 's_adj': 'NSA', 'unit': 'I21'},
     'ppi_industry', 'annual'),
    ('sts_inpp_m_ppi_monthly',
     {'indic_bt': 'PRC_PRR', 'nace_r2': 'B-E36', 's_adj': 'NSA', 'unit': 'I21'},
     'ppi_industry', 'monthly'),
    ('sts_cc_cs_a_construction_costs',
     {'indic_bt': 'COST', 'cpa2_1': 'CPA_F41001_X_410014',
      's_adj': 'NSA', 'unit': 'I21'},
     'construction_cost_residential', 'annual'),
    ('sts_copi_a_construction_prices',
     {'indic_bt': 'PRC_PRR', 'cpa2_1': 'CPA_F41001_X_410014',
      's_adj': 'NSA', 'unit': 'I21'},
     'construction_prices_new_resid', 'annual'),
    # Note: I15 base gives 2000-2023 full coverage; I21 base only starts 2015.
    ('sts_copr_a_construction_production',
     {'indic_bt': 'PRD', 'nace_r2': 'F', 's_adj': 'NSA', 'unit': 'I15'},
     'construction_production_index', 'annual'),
    ('sts_cobp_a_building_permits',
     {'indic_bt': 'BPRM_SQM', 'cpa2_1': 'CPA_F41001_41002',
      's_adj': 'NSA', 'unit': 'I21'},
     'building_permits_sqm', 'annual'),
    ('nrg_pc_205_electricity',
     {'siec': 'E7000', 'nrg_cons': 'MWH2000-19999', 'unit': 'KWH',
      'tax': 'X_TAX', 'currency': 'EUR'},
     'electricity_price_eur_kwh', 'annual'),
    ('nrg_pc_203_gas',
     {'siec': 'G3000', 'nrg_cons': 'GJ1000-9999', 'unit': 'GJ_GCV',
      'tax': 'X_TAX', 'currency': 'EUR'},
     'gas_price_eur_gj', 'annual'),
    ('irt_st_a_short_rates_annual',
     {'int_rt': 'IRT_M3'},
     'short_rate_3m', 'annual'),
    ('irt_lt_mcby_a_long_rates_annual',
     {},
     'long_rate_mcby', 'annual'),
]


def parse_eurostat() -> tuple[pd.DataFrame, pd.DataFrame]:
    annual_ser: dict = {}
    monthly_ser: dict = {}
    for stem, overrides, col, freq in EUROSTAT_SERIES:
        path = RAW_DIR / 'eurostat' / f'{stem}.json'
        if not path.exists():
            print(f'  WARN: {stem}.json missing')
            continue
        doc = json.loads(path.read_text())
        try:
            series = js_extract_series(doc, overrides)
        except KeyError as e:
            print(f'  WARN: {stem}: {e}')
            continue
        if not series:
            print(f'  WARN: {stem}: zero obs after filter')
            continue
        lo, hi = min(series), max(series)
        print(f'  {col:38s} {freq:7s} {len(series):4d} obs [{lo}..{hi}]')
        if freq == 'monthly':
            monthly_ser[col] = series
            continue
        # Collapse sub-annual (bi-annual energy series) to annual mean
        by_year: dict = {}
        for label, v in series.items():
            yr = str(label)[:4]
            by_year.setdefault(yr, []).append(v)
        annual_ser[col] = {yr: float(np.mean(vs)) for yr, vs in by_year.items()}

    annual_df = pd.DataFrame(annual_ser)
    if not annual_df.empty:
        annual_df.index = annual_df.index.astype(int)
        annual_df.index.name = 'year'
        annual_df = annual_df.sort_index()

    monthly_df = pd.DataFrame(monthly_ser)
    if not monthly_df.empty:
        monthly_df.index.name = 'year_month'
        monthly_df = monthly_df.sort_index()

    return annual_df, monthly_df


# ============================================================ #
#  Source 3: Open-Meteo (CHMI weather proxy)
# ============================================================ #

WEATHER_AGG = {
    'temperature_2m_mean': 'mean',
    'temperature_2m_max':  'mean',
    'temperature_2m_min':  'mean',
    'precipitation_sum':   'sum',
    'snowfall_sum':        'sum',
    'rain_sum':            'sum',
    'wind_speed_10m_max':  'mean',
    'wind_gusts_10m_max':  'mean',
    'sunshine_duration':   'sum',
}


def parse_openmeteo() -> pd.DataFrame:
    files = sorted((RAW_DIR / 'chmi').glob('openmeteo_*.json'))
    if not files:
        return pd.DataFrame()

    frames = []
    for f in files:
        doc = json.loads(f.read_text())
        # filename: openmeteo_<station>_YYYY_YYYY.json
        parts = f.stem.split('_')
        station = '_'.join(parts[1:-2])
        df = pd.DataFrame(doc['daily'])
        df['date'] = pd.to_datetime(df['time'])
        df['year'] = df['date'].dt.year
        df['station'] = station
        frames.append(df)
    big = pd.concat(frames, ignore_index=True)

    per_stn = big.groupby(['station', 'year']).agg(WEATHER_AGG).reset_index()
    national = per_stn.groupby('year').agg({k: 'mean' for k in WEATHER_AGG})
    national.columns = [f'weather_{c}' for c in national.columns]

    big['_frost'] = (big['temperature_2m_min'] < 0).astype(int)
    big['_hot'] = (big['temperature_2m_max'] > 30).astype(int)
    days = big.groupby(['station', 'year'])[['_frost', '_hot']].sum()
    days_nat = days.groupby('year').mean()
    days_nat.columns = ['weather_frost_days', 'weather_hot_days']

    return national.join(days_nat)


# ============================================================ #
#  Main
# ============================================================ #

def main() -> int:
    print('=' * 60)
    print('build_external_panel.py — ABGRS strong-exclusion instruments')
    print('=' * 60)

    print('\n[CNB FX]')
    cnb = parse_cnb_fx()
    print(f'  -> {len(cnb)} years x {len(cnb.columns)} cols')

    print('\n[Eurostat]')
    eur_a, eur_m = parse_eurostat()
    print(f'  -> annual: {len(eur_a)} years x {len(eur_a.columns)} cols')
    print(f'  -> monthly: {len(eur_m)} obs x {len(eur_m.columns)} cols')

    print('\n[Open-Meteo / CHMI proxy]')
    wx = parse_openmeteo()
    print(f'  -> {len(wx)} years x {len(wx.columns)} cols')

    annual = cnb.join(eur_a, how='outer').join(wx, how='outer')
    annual = annual.sort_index().loc[2005:2023]

    out_a = OUT_DIR / 'external_panel_annual.csv'
    annual.to_csv(out_a)
    print(f'\nWrote {out_a.relative_to(PROJECT_ROOT)}')
    print(f'  shape: {annual.shape[0]} rows x {annual.shape[1]} cols')

    if not eur_m.empty:
        out_m = OUT_DIR / 'external_panel_monthly.csv'
        eur_m.to_csv(out_m)
        print(f'Wrote {out_m.relative_to(PROJECT_ROOT)}')
        print(f'  shape: {eur_m.shape[0]} rows x {eur_m.shape[1]} cols')

    print('\nAnnual panel non-null counts per column:')
    for col in annual.columns:
        n = annual[col].notna().sum()
        print(f'  {col:42s} {n:3d} / {len(annual)}')

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
