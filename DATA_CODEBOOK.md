# Data Codebook: Markups and Public Procurement in Czech Construction

**Project:** Chadim (2026), "Markups and Public Procurement"
**Last updated:** April 8, 2026

---

## 1. Raw Data Sources

### 1.1 MagnusWeb (Czech Firm Financials)

**Provider:** Bisnode / Dun & Bradstreet Czech Republic
**Coverage:** Czech firms, CZ-NACE F (Construction), 2005-2021
**Unit:** Firm-year
**Location:** `project/thesis/Code_and_data/magnusweb/`

| File group | Files | Content |
|-|-|-|
| Financial statements | `financial1-6.csv` | Sales, costs, fixed assets, equity |
| Financial ratios | `ratios1-3.csv` | Wage share, labor productivity, contribution margins |
| Firm characteristics | `selections1-5.csv` | NACE codes, employment category, legal form, sector |

**Key raw variables:**

| Variable | Description | Unit |
|-|-|-|
| `ID` | Firm identifier (ICO, Czech trade register number) | 8-digit integer |
| `Year` | Fiscal year | 2005-2023 |
| `Sal - Sales, Outputs` | Total revenue (GO) | CZK |
| `C - Costs` | Cost of goods sold (COGS) | CZK |
| `FA - Fixed assets` | Tangible + intangible fixed assets (K) | CZK |
| `WS - Wages/Sales` | Wage bill / Sales ratio | Share |
| `WVA - Wages/Value Added` | Wage bill / Value Added ratio | Share |
| `LP - Labour Productivity` | Value added per employee | CZK/employee |
| `CM III` | Contribution margin III | CZK |
| NACE code | CZ-NACE Rev.2 industry classification | 2-digit |
| Employment category | Employment size bracket (e.g., "10-19") | Categorical |

**Raw sample:** ~45,000 observations, ~4,400 firms across all financial files.

**Employment categories (midpoints used):**

| Category | Midpoint |
|-|-|
| 1-5 | 3.0 |
| 6-9 | 7.5 |
| 10-19 | 14.5 |
| 20-24 | 22.0 |
| 25-49 | 37.0 |
| 50-99 | 74.5 |
| 100-199 | 149.5 |
| 200-249 | 224.5 |
| 250-499 | 374.5 |
| 500-999 | 749.5 |
| 1000+ | 1000.0 |

### 1.2 Datlab Procurement Register

**Provider:** Datlab s.r.o. (Jiri Skuhrovec, CEO)
**Coverage:** Czech public procurement contracts, all sectors, 2006-2022
**Unit:** Contract-bidder (one row per winning bidder per lot)
**Location:** `project/thesis/Code_and_data/master_tender_analytics_202207251530.csv`
**Pipeline copy:** `markups-procurement/1_data/input/datlab/master_tender_analytics.csv`

| Variable | Description | Unit | Coverage |
|-|-|-|-|
| `title` | Contract title (Czech) | Text | 100% |
| `src` | Data source (profil/vestnik) | Text | 100% |
| `bids_count` | Number of bidders | Integer | ~70% |
| `lot_estimated_price` | Engineer cost estimate (ex-ante) | CZK | 23.4% |
| `bid_final_price` | Realized contract price | CZK | 90.8% |
| `bidder_id` | Winning firm ICO | 8-digit | 100% |
| `bidder_name` | Winning firm name | Text | 100% |
| `year` | Year field | Integer | 100% |
| `contract_signature_date` | Signature date | Date | ~93% |
| `award_date` | Award date | Date | ~60% |

**Raw sample:** 934,722 contracts, 61,475 unique firms.

**Engineer cost estimates:** Prepared by specialized employees of procuring authorities before tender publication. Regression of log(price) on log(estimate) yields R-squared = 0.91 (Baranek & Titl 2024). Available for ~218K contracts (23.4%). Construction has higher coverage (~40%).

### 1.3 Deflators

**Provider:** Czech Statistical Office (CZSO) via Eurostat
**Location:** `project/thesis/Code_and_data/deflators.csv`
**Coverage:** 2,466 rows (NACE2 x year, 1993-2023)

| Deflator | Applied to | Base year |
|-|-|-|
| `deflatorGFCP` | Fixed assets (K) | 2015=1.0 |
| `deflatorINTP` | Intermediate inputs (II) | 2015=1.0 |
| `deflatorPRDP` | Output, COGS | 2015=1.0 |
| `deflatorVALP` | Value added | 2015=1.0 |
| `deflatorGDP` | GDP deflator | 2015=1.0 |
| `deflatorCPI` | Consumer prices | 2015=1.0 |

Deflation is NACE2-year specific: each firm's variables are deflated using its own 2-digit industry deflator for that year.

### 1.4 Orbis Historical (BvD/Moody's)

**Provider:** Bureau van Dijk, accessed via WRDS
**Coverage:** Czech firms, all industries, 2006-2023
**Location:** `project/thesis/replication/outputs/data/orbis_cz_all_tiers.csv` (10.8M obs)
**Status:** Available for extended analysis; not used in baseline pipeline

Key variables: `opre` (operating revenue), `cost` (total costs), `mate` (material costs), `empl` (employees, continuous), `fias` (fixed assets), BVDID, ICO.

---

## 2. Data Processing Pipeline

### 2.1 Step 1: Load and Clean (`rebuild_data.py`)

1. **Concatenate** 6 financial CSVs, 3 ratios CSVs, 5 selections CSVs
2. **Standardize** column names: `ID` -> `id`, `Year` -> `year`, `Sal - Sales, Outputs` -> `sales`, etc.
3. **Employment:** Convert categorical brackets to numeric midpoints (see table above)
4. **Merge** financial + ratios + selections on (`id`, `year`)
5. **Deduplicate:** Keep observation with most non-null financial variables per (`id`, `year`)

### 2.2 Step 2: Procurement Merge

1. **Load** tender data from Datlab CSV
2. **Aggregate** to firm-year: sum `bid_final_price` -> `pp_sales`; count contracts -> `n_contracts`
3. **Compute** firm-year competition: mean `bids_count` -> `avg_bids`; share of single-bid contracts -> `single_bid_share`
4. **Merge** with financial panel on `id` (ICO) and `year`
5. **Create** treatment variables: `pp_dummy = 1{pp_sales > 0}`

### 2.3 Step 3: Variable Construction

**Production function inputs (nominal):**

| Variable | Formula | Description |
|-|-|-|
| `GO` | `sales` | Gross output |
| `W` | `WS * sales` | Wage bill (wages/sales ratio x sales) |
| `II` | `(1 - WVA) * (GO - W/(WVA))` | Intermediate inputs (derived from WVA ratio) |
| `COGS` | `costs` | Cost of goods sold (direct materials) |
| `VA` | `GO - II` | Value added |
| `K` | `assets` | Fixed assets (capital stock) |
| `L` | `VA / LP` | Employment (from labor productivity) |
| `O` | `II - COGS` | Overhead / services (residual intermediates) |

**Key identity:** `GO = W + II = W + COGS + O`. COGS is ~50% of GO (direct materials); II is ~82% of GO (includes overhead, rent, utilities).

**Deflated variables** (prefix `r`):

| Variable | Deflator used |
|-|-|
| `rGO` | `deflatorPRDP` (producer prices) |
| `rVA` | `deflatorVALP` (value added) |
| `rII` | `deflatorINTP` (intermediate) |
| `rW` | `deflatorCPI` (consumer prices) |
| `rK` | `deflatorGFCP` (gross fixed capital) |
| `rCOGS` | `deflatorPRDP` (producer prices) |
| `rO` | `deflatorINTP` (intermediate) |

**Log variables** (lowercase): `go = log(rGO)`, `k = log(rK)`, `cogs = log(rCOGS)`, etc.

**Market structure:**

| Variable | Formula |
|-|-|
| `mktshare` | `rGO / sum(rGO)` within NACE2-year |
| `salessharefirm` | Alias for `mktshare` |

### 2.4 Step 4: Procurement Treatment Variables

| Variable | Definition | Coverage |
|-|-|-|
| `pp_dummy` | 1 if firm won any procurement contract in year t | 100% |
| `pp_share` | Procurement revenue / Total revenue in year t | 100% |
| `pp_sales` | Total procurement revenue in year t (CZK) | 100% (0 for non-PP) |
| `n_contracts` | Number of contracts won in year t | 100% |
| `avg_bids` | Mean number of bidders across firm's contracts | 30% (PP firms only) |
| `single_bid_share` | Share of contracts with single bidder | 30% (PP firms only) |
| `avg_discount` | Mean (estimate - price) / estimate | 20% (PP firms with estimates) |
| `hhi_revenue` | Herfindahl of firm's contract values | 30% |
| `pp_sales_L1`, `L2` | Lagged procurement revenue (t-1, t-2) | Lagged |
| `pp_stock_2y`, `3y` | Cumulative PP revenue over 2/3 years | Derived |
| `pp_active_2y`, `3y` | 1 if PP active in any of last 2/3 years | Derived |
| `pp_ever_2y`, `3y` | 1 if ever in PP within 2/3 year window | Derived |
| `pp_entry_year` | First year of procurement participation | Derived |
| `pp_years_since_entry` | Years since first procurement contract | Derived |
| `pp_cumul_contracts` | Lifetime count of procurement contracts | Derived |
| `pp_share_L1` | Lagged procurement revenue share | Lagged |

### 2.5 Step 5: Winsorization and Sample Restrictions

**Winsorization:** Log variables winsorized at 2nd and 98th percentile within NACE2-year cells.

**Trimming rules:**
- Drop if `O/GO > 0.90` (overhead exceeds 90% of output)
- Drop if `K/GO > 5` (implausible capital intensity)
- Drop if `COGS/GO > 1` (costs exceed output)
- Drop if negative deflated values

**Panel requirement:** Firms must have at least 2 consecutive years of data (for Markov productivity transition in ACF).

**Final sample after all restrictions:**

| Statistic | Value |
|-|-|
| Observations | 9,164 |
| Unique firms | 1,498 |
| Years | 2006-2021 |
| NACE 41 (Buildings) | 4,911 obs (53.6%) |
| NACE 42 (Civil engineering) | 818 obs (8.9%) |
| NACE 43 (Specialized) | 3,435 obs (37.5%) |
| Procurement rate (overall) | 42.2% of firm-years |
| Procurement rate (2006) | 34.1% |
| Procurement rate (2021) | 49.5% |

---

## 3. Estimation Outputs

### 3.1 ACF Production Function Estimates

**Script:** `markups-procurement/2_analysis/source/acf_estimator.py`
**Method:** Ackerberg, Caves & Frazer (2015) with Nelder-Mead GMM
**Output:** `paper_markups.csv` / `paper_markups.dta`

| Column | Description |
|-|-|
| `id` | Firm ICO |
| `year` | Year |
| `markup` | Firm-year markup estimate (mu = theta / alpha) |
| `omega` | Log productivity (Markov residual from 2nd stage) |
| `alphahat` | Expenditure share (COGS / GO) |
| `spec` | Specification identifier |
| `nace2` | 2-digit NACE industry |

**Specifications:**

| Spec | Description | Production function | Markov transition |
|-|-|-|-|
| A | Baseline Cobb-Douglas | `go = beta_k * k + beta_cogs * cogs` | `omega_t = rho * omega_{t-1} + xi_t` |
| B | CD + PP in Markov | Same | `omega_t = rho_0 * omega_{t-1} + rho_1 * pp_dummy_{t-1} + xi_t` |
| C | Translog | Includes k^2, cogs^2, k*cogs | Same as A |
| D | Translog + PP Markov | Translog | Same as B |
| E | CD + PP in PF | `go = ... + beta_pp * pp_dummy` | Same as A |
| OLS | Naive OLS | Same as A | None (no 2nd stage) |

**Markup formula:** `mu_jt = theta_v / alpha_v_jt` where `theta_v` is the estimated output elasticity of the variable input (COGS) and `alpha_v_jt = COGS_jt / GO_jt` is the observed expenditure share.

### 3.2 Rel_Price (Engineer Cost Estimate Ratio)

**Script:** `favoritism_decomposition.py`, Panel F
**Source:** Datlab tender data, `lot_estimated_price` and `bid_final_price`

| Variable | Definition |
|-|-|
| `rel_price` | `bid_final_price / lot_estimated_price` per contract |
| `mean_rel_price` | Firm-year average of `rel_price` (trimmed 0.2-5.0) |

**Coverage:** 520 firms in panel matched with estimate data; 385 procurement firms with both markups and Rel_Price; 19,880 contracts, 3,143 firm-year observations.

**Key statistics:**

| Statistic | Value |
|-|-|
| Mean Rel_Price (panel PP firms) | 0.883 |
| Median Rel_Price | 0.900 |
| SD | 0.209 |
| Mean Rel_Price (B&T full sample) | 0.844 |

### 3.3 Production Function Coefficients

**Output:** `paper_pf_estimates.csv`

Baseline Cobb-Douglas (Spec A) by industry:

| Parameter | NACE 41 | NACE 42 | NACE 43 |
|-|-|-|-|
| theta_cogs | ~0.93 | ~0.90 | ~0.85 |
| theta_k | ~0.10 | ~0.12 | ~0.08 |
| Mean markup | ~2.36 | ~2.59 | ~1.81 |
| Procurement premium | ~14% | ~14% | ~14% |

---

## 4. Variable Codebook (Final Panel)

### 4.1 Identifiers and Structure

| Variable | Type | Description |
|-|-|-|
| `id` | int | Firm identifier (ICO, Czech trade register) |
| `year` | int | Fiscal year (2006-2021) |
| `nace` | str | Full NACE code |
| `nace2` | float | 2-digit NACE (41, 42, 43) |
| `legal_form` | str | Legal form of organization |
| `inst_sector` | str | Institutional sector |
| `foreign` | float | Foreign ownership indicator |

### 4.2 Production Function Variables (Log, Deflated)

| Variable | Mean | SD | Min | Max | Description |
|-|-|-|-|-|-|
| `go` | 18.13 | 1.20 | 15.78 | 23.19 | log(rGO) — log real gross output |
| `k` | 15.94 | 1.64 | 12.06 | 21.04 | log(rK) — log real capital |
| `cogs` | 17.41 | 1.13 | 14.97 | 22.21 | log(rCOGS) — log real cost of goods sold |
| `ii` | 17.92 | 1.25 | 15.26 | 23.10 | log(rII) — log real intermediate inputs |
| `va` | 16.35 | 1.23 | 13.71 | 20.90 | log(rVA) — log real value added |
| `w` | 16.36 | 1.16 | 13.93 | 21.03 | log(rW) — log real wage bill |
| `o` | 17.00 | 1.57 | 13.20 | 22.63 | log(rO) — log real overhead |
| `l` | — | — | — | — | log(rL) — log employment (sparse) |
| `le` | 3.70 | 0.91 | 2.67 | 5.93 | log(empl_mid) — log employment midpoint |

### 4.3 Production Function Variables (Nominal, CZK)

| Variable | Mean | SD | Median | Description |
|-|-|-|-|-|
| `GO` | 238M | 942M | 75M | Gross output |
| `K` | 37M | 147M | 9M | Fixed assets |
| `COGS` | 102M | 375M | 35M | Cost of goods sold |
| `II` | 202M | 840M | 60M | Intermediate inputs |
| `VA` | 37M | 123M | 12M | Value added |
| `W` | 33M | 108M | 12M | Wage bill |

### 4.4 Procurement Variables

| Variable | Mean | SD | Description |
|-|-|-|-|
| `pp_dummy` | 0.42 | 0.49 | =1 if any procurement revenue in year t |
| `pp_share` | 0.15 | 0.29 | Procurement revenue / total revenue |
| `pp_sales` | — | — | Total procurement revenue (CZK) |
| `n_contracts` | — | — | Number of contracts won |
| `avg_bids` | — | — | Mean bidders per contract |
| `single_bid_share` | — | — | Share of single-bidder contracts |
| `avg_discount` | 0.12 | — | Aggregate (estimate-price)/estimate |
| `pp_dummy_L1` | — | — | Lagged pp_dummy (t-1) |
| `pp_share_L1` | — | — | Lagged pp_share (t-1) |
| `pp_stock_2y` | — | — | Cumulative PP revenue (t to t-1) |
| `pp_active_2y` | — | — | =1 if PP active in t or t-1 |
| `pp_entry_year` | — | — | First year of PP participation |
| `pp_years_since_entry` | — | — | Years since first PP contract |

### 4.5 Market Structure

| Variable | Mean | SD | Description |
|-|-|-|-|
| `mktshare` | 0.002 | 0.010 | Firm revenue share within NACE2-year |
| `hhi_revenue` | — | — | Herfindahl of firm's contract values |
| `entry` | — | — | =1 in firm's first observed year |
| `exit` | — | — | =1 in firm's last observed year (before 2021) |
| `empl_mid` | — | — | Employment category midpoint |

---

## 5. Key Data Relationships

### 5.1 Variable Input Choice

**COGS is the variable input** for production function estimation, not II.
- COGS = direct materials (~50% of GO). Well-defined accounting item.
- II = intermediates (~82% of GO). Includes overhead (rent, utilities, services) that does not map to a flexible input in the production function.
- Using II without modeling overhead creates alpha mismatch and economically meaningless markups.
- See `variable_input_decision.md` for full rationale.

### 5.2 Markup Identification

Markup `mu_jt = theta_cogs / alpha_cogs_jt` where:
- `theta_cogs` is estimated from ACF GMM (industry-level parameter)
- `alpha_cogs_jt = COGS_jt / GO_jt` is observed in the data (firm-year)

The **procurement premium** is identified from variation in `alpha` (expenditure shares), NOT from `theta` (which is constant within industry under Cobb-Douglas). This is why the premium is robust across 9 estimation methods despite theta non-identification (GNR 2020 theorem).

### 5.3 Merge Keys

| Merge | Key | Notes |
|-|-|-|
| MagnusWeb financials + ratios | `id`, `year` | Inner join |
| Panel + Datlab tenders | `id` = `bidder_id` (ICO), `year` | Left join (non-PP firms get 0) |
| Panel + Deflators | `nace2`, `year` | Left join |
| Panel + Orbis | `id` = national ICO via BVDID mapping | For extended analysis only |

---

## 6. External Benchmarks

| Source | Key statistic | Our comparable |
|-|-|-|
| Baranek & Titl (2024 JLE) | +6% overpricing from connections | Our Rel_Price test: <8% of premium explained |
| Baranek & Titl (2024 JLE) | Mean Rel_Price = 0.844 (full sample) | Our construction: 0.883 |
| Baranek & Titl (2024 JLE) | N = 34,007 contracts | Our tender data: 934,722 contracts |
| KLMS (2025 AER) | US (1-epsilon) = 0.863 | Czech (1-epsilon) = 0.896 |
| KLMS (2025 AER) | US price markup = 1.16 | Czech price markup = 1.12 |

---

## 7. Known Limitations

1. **Employment is categorical** in MagnusWeb (brackets, not continuous). Midpoints used. This prevents KLMS-style labor supply elasticity (theta) estimation. Orbis or CZSO data needed.
2. **Engineer estimates** (`lot_estimated_price`) have 23.4% coverage overall. Coverage is higher for construction and larger contracts.
3. **Political connection indicator** not yet available. Pending from Vitek Titl (B&T 2024). Connections affect ~2% of contracts / 6% of volume.
4. **Bid-level data** (all bidders' prices) not available in structured form. Only winner observed. Prevents KLMS bid inversion.
5. **Single NACE code per firm.** Multi-product firms classified by primary activity. Some misclassification at 2-digit level.
6. **Unbalanced panel.** Firms enter/exit sample. Entry/exit flags provided but survival correction not applied in baseline.
