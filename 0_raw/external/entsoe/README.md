# ENTSO-E Transparency Platform — CZ day-ahead prices

Set `ENTSOE_API_TOKEN` env var and re-run `python fetch_external.py --only entsoe`.

## Obtaining a security token

1. Register at <https://transparency.entsoe.eu/>
2. Email transparency@entsoe.eu with subject "Restful API access" requesting a security token.
3. Export the token: `export ENTSOE_API_TOKEN=...`

## What this script fetches

- Day-ahead wholesale electricity prices, CZ bidding zone (10YCZ-CEPS-----N)
- Years: 2008–2023
- Format: XML per year
