#!/usr/bin/env python3
"""Fetch external exogenous shifters for ABGRS strong-exclusion instruments.

Sources (all public open-data; ENTSO-E needs a free token):
    1. Eurostat dissemination API — PPI, energy prices, construction costs,
       building permits, interest rates (geo=CZ filter)
    2. Open-Meteo historical archive — daily weather for three Czech stations
       (Prague-Ruzyne, Brno-Turany, Ostrava-Mosnov); serves as CHMI proxy
    3. ENTSO-E Transparency Platform — CZ day-ahead wholesale electricity
       prices (requires ENTSOE_API_TOKEN in env)
    4. CNB — verify existing daily FX file stash

Usage:
    python fetch_external.py                       # fetch missing, skip existing
    python fetch_external.py --force               # re-download everything
    python fetch_external.py --only eurostat       # single source
    python fetch_external.py --only openmeteo,cnb  # subset
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent
EUROSTAT_DIR = ROOT / 'eurostat'
CHMI_DIR = ROOT / 'chmi'
CZSO_DIR = ROOT / 'czso'
ENTSOE_DIR = ROOT / 'entsoe'
CNB_DIR = ROOT / 'cnb'
for _d in (EUROSTAT_DIR, CHMI_DIR, CZSO_DIR, ENTSOE_DIR, CNB_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ============================================================ #
#  Eurostat — dissemination API v1.0
# ============================================================ #

EUROSTAT_BASE = ('https://ec.europa.eu/eurostat/api/dissemination/'
                 'statistics/1.0/data')

# (dataset_id, filename, query-parameter dict)
# Queries are all filtered to geo=CZ; additional filters keep payload under
# the async threshold that triggers the 413 "asynchronous_response" code.
EUROSTAT_DATASETS = [
    # --- Industry producer prices (materials & energy inputs) ---
    ('sts_inpp_a',
     'sts_inpp_a_ppi_annual.json',
     {'geo': 'CZ'}),
    ('sts_inpp_m',
     'sts_inpp_m_ppi_monthly.json',
     {'geo': 'CZ', 'indic_bt': 'PRC_PRR', 'nace_r2': 'B-E36'}),

    # --- Construction-specific producer prices and costs ---
    ('sts_cc_cs_a',
     'sts_cc_cs_a_construction_costs.json',
     {'geo': 'CZ'}),
    ('sts_copi_a',
     'sts_copi_a_construction_prices.json',
     {'geo': 'CZ'}),

    # --- Construction output / demand proxies ---
    ('sts_copr_a',
     'sts_copr_a_construction_production.json',
     {'geo': 'CZ'}),
    ('sts_cobp_a',
     'sts_cobp_a_building_permits.json',
     {'geo': 'CZ'}),

    # --- Energy input costs (non-household) ---
    ('nrg_pc_205',
     'nrg_pc_205_electricity.json',
     {'geo': 'CZ'}),
    ('nrg_pc_203',
     'nrg_pc_203_gas.json',
     {'geo': 'CZ'}),

    # --- Financing / demand shifters ---
    ('irt_st_a',
     'irt_st_a_short_rates_annual.json',
     {'geo': 'CZ'}),
    ('irt_lt_mcby_a',
     'irt_lt_mcby_a_long_rates_annual.json',
     {'geo': 'CZ'}),
]


def _get_eurostat(url: str, params: dict,
                  max_retries: int = 8, wait_s: int = 15) -> dict | None:
    """Eurostat returns HTTP 200 with `error.status=413` while preparing
    a large payload asynchronously. Retry with fixed wait until ready."""
    for attempt in range(1, max_retries + 1):
        try:
            r = requests.get(url, params=params, timeout=90)
        except requests.RequestException as e:
            print(f'    request error (attempt {attempt}): {e}')
            time.sleep(wait_s)
            continue
        if r.status_code != 200:
            print(f'    HTTP {r.status_code} (attempt {attempt})')
            time.sleep(wait_s)
            continue
        try:
            body = r.json()
        except ValueError:
            print(f'    non-JSON body (attempt {attempt})')
            time.sleep(wait_s)
            continue
        err = body.get('error') if isinstance(body, dict) else None
        if err and isinstance(err, list) and err[0].get('status') == 413:
            print(f'    async prep (413), attempt {attempt}, wait {wait_s}s')
            time.sleep(wait_s)
            continue
        if err:
            print(f'    API error: {err}')
            return None
        return body
    print('    gave up after max retries')
    return None


def fetch_eurostat(force: bool = False) -> dict:
    results = {}
    for dataset, filename, params in EUROSTAT_DATASETS:
        out = EUROSTAT_DIR / filename
        if out.exists() and out.stat().st_size > 500 and not force:
            results[dataset] = f'cached ({out.stat().st_size:,} B)'
            print(f'  [eurostat] {dataset}: {results[dataset]}')
            continue
        url = f'{EUROSTAT_BASE}/{dataset}'
        q = {'format': 'JSON', 'lang': 'EN', **params}
        body = _get_eurostat(url, q)
        if body is None:
            results[dataset] = 'failed'
            print(f'  [eurostat] {dataset}: failed')
            continue
        out.write_text(json.dumps(body), encoding='utf-8')
        results[dataset] = f'ok ({out.stat().st_size:,} B)'
        print(f'  [eurostat] {dataset}: {results[dataset]}')
        time.sleep(1.0)
    return results


# ============================================================ #
#  Open-Meteo historical archive — CHMI weather proxy
# ============================================================ #

OPENMETEO_BASE = 'https://archive-api.open-meteo.com/v1/archive'
STATIONS = [
    ('praha_ruzyne',  50.1008, 14.2632),
    ('brno_turany',   49.1513, 16.6994),
    ('ostrava_mosnov', 49.6963, 18.1110),
]
WEATHER_VARS = [
    'temperature_2m_mean', 'temperature_2m_max', 'temperature_2m_min',
    'precipitation_sum', 'snowfall_sum', 'rain_sum',
    'wind_speed_10m_max', 'wind_gusts_10m_max', 'sunshine_duration',
]


def fetch_openmeteo(force: bool = False,
                    start: str = '2005-01-01',
                    end: str = '2023-12-31') -> dict:
    results = {}
    for name, lat, lon in STATIONS:
        out = CHMI_DIR / f'openmeteo_{name}_{start[:4]}_{end[:4]}.json'
        if out.exists() and out.stat().st_size > 1000 and not force:
            results[name] = f'cached ({out.stat().st_size:,} B)'
            print(f'  [open-meteo] {name}: {results[name]}')
            continue
        q = {
            'latitude': lat, 'longitude': lon,
            'start_date': start, 'end_date': end,
            'daily': ','.join(WEATHER_VARS),
            'timezone': 'Europe/Prague',
        }
        try:
            r = requests.get(OPENMETEO_BASE, params=q, timeout=180)
            r.raise_for_status()
            out.write_text(r.text, encoding='utf-8')
            results[name] = f'ok ({out.stat().st_size:,} B)'
            print(f'  [open-meteo] {name}: {results[name]}')
            time.sleep(1.5)
        except requests.RequestException as e:
            print(f'  [open-meteo] {name}: failed ({e})')
            results[name] = f'failed: {e}'
    return results


# ============================================================ #
#  ENTSO-E Transparency Platform — day-ahead CZ prices
# ============================================================ #

ENTSOE_BASE = 'https://web-api.tp.entsoe.eu/api'
CZ_DOMAIN = '10YCZ-CEPS-----N'
ENTSOE_YEARS = range(2008, 2024)  # CZ bidding zone day-ahead from 2008


def fetch_entsoe(force: bool = False) -> dict:
    token = os.environ.get('ENTSOE_API_TOKEN')
    if not token:
        readme = ENTSOE_DIR / 'README.md'
        readme.write_text(
            '# ENTSO-E Transparency Platform — CZ day-ahead prices\n\n'
            'Set `ENTSOE_API_TOKEN` env var and re-run '
            '`python fetch_external.py --only entsoe`.\n\n'
            '## Obtaining a security token\n\n'
            '1. Register at <https://transparency.entsoe.eu/>\n'
            '2. Email transparency@entsoe.eu with subject '
            '"Restful API access" requesting a security token.\n'
            '3. Export the token: `export ENTSOE_API_TOKEN=...`\n\n'
            '## What this script fetches\n\n'
            f'- Day-ahead wholesale electricity prices, CZ bidding zone '
            f'({CZ_DOMAIN})\n'
            f'- Years: {min(ENTSOE_YEARS)}–{max(ENTSOE_YEARS)}\n'
            '- Format: XML per year\n',
            encoding='utf-8',
        )
        print('  [entsoe] no token — wrote README stub')
        return {'entsoe': 'no token; stub written'}
    results = {}
    for year in ENTSOE_YEARS:
        out = ENTSOE_DIR / f'entsoe_dayahead_cz_{year}.xml'
        if out.exists() and out.stat().st_size > 10000 and not force:
            results[str(year)] = f'cached ({out.stat().st_size:,} B)'
            print(f'  [entsoe] {year}: {results[str(year)]}')
            continue
        q = {
            'securityToken': token,
            'documentType': 'A44',
            'in_Domain': CZ_DOMAIN,
            'out_Domain': CZ_DOMAIN,
            'periodStart': f'{year}01010000',
            'periodEnd': f'{year}12312300',
        }
        try:
            r = requests.get(ENTSOE_BASE, params=q, timeout=240)
            r.raise_for_status()
            out.write_text(r.text, encoding='utf-8')
            results[str(year)] = f'ok ({out.stat().st_size:,} B)'
            print(f'  [entsoe] {year}: {results[str(year)]}')
            time.sleep(2.0)
        except requests.RequestException as e:
            print(f'  [entsoe] {year}: failed ({e})')
            results[str(year)] = f'failed: {e}'
    return results


# ============================================================ #
#  CNB — verify existing FX stash (fetcher lives in separate legacy script)
# ============================================================ #

def fetch_cnb(force: bool = False) -> dict:
    fx_files = sorted(CNB_DIR.glob('cnb_fx_*.txt'))
    if not fx_files:
        print('  [cnb] FX files not found')
        return {'fx': 'missing — run the legacy fetch first'}
    years = [int(f.stem.split('_')[-1]) for f in fx_files]
    msg = f'{len(fx_files)} files, {min(years)}-{max(years)}'
    print(f'  [cnb] {msg}')
    return {'fx': msg}


# ============================================================ #
#  Main driver
# ============================================================ #

SOURCES = {
    'eurostat': fetch_eurostat,
    'openmeteo': fetch_openmeteo,
    'entsoe': fetch_entsoe,
    'cnb': fetch_cnb,
}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--force', action='store_true',
                   help='re-download even if file exists')
    p.add_argument('--only', default=None,
                   help='comma-separated source(s): ' + ','.join(SOURCES))
    args = p.parse_args()

    wanted = list(SOURCES) if not args.only else [
        s.strip() for s in args.only.split(',')]
    bad = [s for s in wanted if s not in SOURCES]
    if bad:
        print(f'Unknown source(s): {bad}. Available: {list(SOURCES)}',
              file=sys.stderr)
        return 2

    print('=' * 60)
    print('fetch_external.py — ABGRS strong-exclusion instruments')
    print('=' * 60)

    summary = {}
    for src in wanted:
        print(f'\n[{src}]')
        try:
            summary[src] = SOURCES[src](force=args.force)
        except Exception as e:
            print(f'  [{src}] UNCAUGHT: {e}')
            summary[src] = {'error': str(e)}

    print('\n' + '=' * 60)
    print('Summary')
    print('=' * 60)
    for src, res in summary.items():
        print(f'\n{src}:')
        if isinstance(res, dict):
            for k, v in res.items():
                print(f'  {k}: {v}')
        else:
            print(f'  {res}')

    any_failed = any(
        isinstance(res, dict) and any('failed' in str(v) for v in res.values())
        for res in summary.values()
    )
    return 1 if any_failed else 0


if __name__ == '__main__':
    sys.exit(main())
