"""KLMS-Inspired Double Market Power Analysis for Czech Construction.

Implements reduced-form analogs of Kroft, Luo, Mogstad & Setzler (2025, AER)
using Czech MagnusWeb firm data + Datlab procurement records.

Analyses:
  1. Procurement event study (KLMS Table 1 analog)
  2. Labor supply elasticity θ via DiD/2SLS (KLMS Proposition 3)
  3. Revenue elasticity (1-ε) (KLMS eq. 25)
  4. Labor share β_L (KLMS eq. 24)
  5. Double markdown/markup (KLMS Table 3 analog)
  6. Crowd-out test

References
----------
KLMS (2025): Imperfect Competition and Rents in Labor and Product Markets:
    The Case of the Construction Industry, AER 115(9): 2926-2969.

Author: Marek Chadim (Yale, Tobin Center)
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS
from linearmodels.panel import PanelOLS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# ========================================================================== #
#  Paths
# ========================================================================== #

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
REBUILT_DTA = INPUT_DIR / 'data_rebuilt.dta'
ORBIS_DTA = INPUT_DIR / 'orbis_panel_construction.dta'
TABLE_DIR = OUTPUT_DIR / 'tables'
FIG_DIR = OUTPUT_DIR / 'figures'


def winsorize(s: pd.Series, lb: float = 0.005, ub: float = 0.995) -> pd.Series:
    """Winsorize a series at given quantiles (KLMS convention)."""
    q_lo, q_hi = s.quantile(lb), s.quantile(ub)
    return s.clip(q_lo, q_hi)


# ========================================================================== #
#  Data Loading
# ========================================================================== #

class KLMSAnalysis:
    """KLMS-inspired double market power analysis for Czech construction.

    Parameters
    ----------
    data_path : Path or str
        Path to firm-year panel .dta file.
    use_orbis : bool
        If True, load orbis_panel_construction.dta (continuous employment)
        instead of data_rebuilt.dta (MagnusWeb bracket employment).
    """

    def __init__(self, data_path: Optional[Path] = None, use_orbis: bool = False):
        self.use_orbis = use_orbis
        if data_path:
            self.data_path = Path(data_path)
        elif use_orbis:
            self.data_path = ORBIS_DTA
        else:
            self.data_path = REBUILT_DTA
        self.df = self._load_data()
        self.results: Dict[str, any] = {}

    def _load_data(self) -> pd.DataFrame:
        """Load data and construct KLMS-relevant variables."""
        source = 'Orbis' if self.use_orbis else 'MagnusWeb'
        print(f'Loading {source} data from {self.data_path}')
        df = pd.read_stata(self.data_path)

        # Orbis uses 'ico' as firm ID; MagnusWeb uses 'id'
        if self.use_orbis and 'ico' in df.columns and 'id' not in df.columns:
            df = df.rename(columns={'ico': 'id'})

        # Ensure numeric types
        df['id'] = pd.to_numeric(df['id'], errors='coerce')
        df['year'] = pd.to_numeric(df['year'], errors='coerce')
        df = df.dropna(subset=['id', 'year'])
        df['id'] = df['id'].astype(int)
        df['year'] = df['year'].astype(int)

        # Orbis has continuous 'empl' — compute le if missing
        if self.use_orbis and 'empl' in df.columns:
            if 'le' not in df.columns or df['le'].isna().all():
                df['le'] = np.log(df['empl'].clip(lower=0.5))
            # Fill le from empl where le is NaN but empl exists
            mask = df['le'].isna() & df['empl'].notna() & (df['empl'] > 0)
            df.loc[mask, 'le'] = np.log(df.loc[mask, 'empl'].clip(lower=0.5))

        # --- Construct KLMS variables ---

        # Log average wage per worker (w = log wage bill, le = log employment)
        df['avg_wage'] = df['w'] - df['le']  # NaN where le is missing

        # Revenue and cost shares (in levels, from logs)
        df['W_level'] = np.exp(df['w'])      # real wage bill
        df['GO_level'] = np.exp(df['go'])    # real gross output
        df['COGS_level'] = np.exp(df['cogs'])  # real COGS
        df['s_L'] = df['W_level'] / df['GO_level']       # labor share of revenue
        df['s_M'] = df['COGS_level'] / df['GO_level']    # intermediates share

        # Private-market revenue proxy (GO minus procurement sales)
        df['pp_sales_clean'] = df['pp_sales'].fillna(0).clip(lower=0)
        private_rev = df['GO_level'] - df['pp_sales_clean']
        df['private_rev'] = np.where(private_rev > 0, np.log(private_rev), np.nan)

        # First differences (within firm, sorted by year)
        df = df.sort_values(['id', 'year'])
        for var in ['go', 'cogs', 'k', 'w', 'le', 'avg_wage', 'private_rev']:
            df[f'd_{var}'] = df.groupby('id')[var].diff()

        # Procurement entry indicator (for first-differenced DiD)
        # pp_entry = 1 in the year of first procurement entry
        df['pp_entry'] = ((df['pp_years_since_entry'] == 0) &
                          df['pp_entry_year'].notna()).astype(int)

        # Treatment dummy for DiD: any procurement in current year
        df['D'] = df['pp_dummy'].astype(int)

        # Never-treated indicator
        df['never_treated'] = df.groupby('id')['pp_entry_year'].transform(
            lambda x: x.isna().all()
        ).astype(int)

        # Event time (years relative to first procurement entry)
        df['event_time'] = df['pp_years_since_entry']

        # Year dummies for FE regressions
        df['year_cat'] = pd.Categorical(df['year'])

        print(f'  Panel: {df["id"].nunique():,} firms, {len(df):,} obs, '
              f'{df["year"].min()}-{df["year"].max()}')
        print(f'  Employment coverage: {df["le"].notna().sum():,}/{len(df):,} '
              f'({df["le"].notna().mean():.0%})')
        print(f'  Procurement firms: {df.loc[df["D"]==1, "id"].nunique():,}, '
              f'never-treated: {df.loc[df["never_treated"]==1, "id"].nunique():,}')

        return df

    # ================================================================== #
    #  Analysis 1: Procurement Event Study
    # ================================================================== #

    def event_study(self,
                    outcomes: List[str] = ['le', 'w', 'avg_wage', 'go', 'cogs'],
                    window: Tuple[int, int] = (-3, 4),
                    ref_period: int = -1) -> pd.DataFrame:
        """Staggered DiD event study around first procurement entry.

        Follows KLMS Table 1: firm FE + year FE, event-time dummies
        interacted with treatment, omitting reference period.
        """
        print('\n=== Analysis 1: Procurement Event Study ===')

        df = self.df.copy()
        # Keep firms with valid event time OR never-treated
        mask = df['event_time'].notna() | (df['never_treated'] == 1)
        df = df[mask].copy()

        # For never-treated, set event_time to NaN (they contribute to FEs)
        # Create event-time dummies for treated firms
        e_min, e_max = window
        results_list = []

        for outcome in outcomes:
            sub = df.dropna(subset=[outcome]).copy()
            if len(sub) < 100:
                print(f'  {outcome}: too few obs ({len(sub)}), skipping')
                continue

            # Set panel index
            sub = sub.set_index(['id', 'year'])

            # Create event-time dummies (excluding reference period)
            for e in range(e_min, e_max + 1):
                col = f'evt_{e}'
                sub[col] = ((sub['event_time'] == e)).astype(float)
                # Never-treated always 0
                sub.loc[sub['never_treated'] == 1, col] = 0.0

            evt_cols = [f'evt_{e}' for e in range(e_min, e_max + 1)
                        if e != ref_period]
            ref_col = f'evt_{ref_period}'
            if ref_col in sub.columns:
                sub = sub.drop(columns=[ref_col])

            # Firm FE + year FE regression
            try:
                formula_parts = ' + '.join(evt_cols)
                mod = PanelOLS(sub[outcome], sub[evt_cols],
                               entity_effects=True, time_effects=True)
                res = mod.fit(cov_type='clustered', cluster_entity=True)

                for e in range(e_min, e_max + 1):
                    col = f'evt_{e}'
                    if e == ref_period:
                        coef, se, pval = 0.0, 0.0, 1.0
                    else:
                        coef = res.params.get(col, np.nan)
                        se = res.std_errors.get(col, np.nan)
                        pval = res.pvalues.get(col, np.nan)
                    results_list.append({
                        'outcome': outcome, 'event_time': e,
                        'coef': coef, 'se': se, 'pval': pval,
                        'ci_lo': coef - 1.96 * se, 'ci_hi': coef + 1.96 * se,
                    })
                print(f'  {outcome}: N={res.nobs:,}, '
                      f'post-entry effect = {res.params.get("evt_0", np.nan):.4f} '
                      f'(SE {res.std_errors.get("evt_0", np.nan):.4f})')
            except Exception as exc:
                print(f'  {outcome}: regression failed — {exc}')

        results_df = pd.DataFrame(results_list)
        self.results['event_study'] = results_df
        return results_df

    # ================================================================== #
    #  Analysis 2: Labor Supply Elasticity θ
    # ================================================================== #

    def estimate_theta(self) -> Dict[str, float]:
        """Estimate inverse labor supply elasticity θ via KLMS Proposition 3.

        Three approaches:
        1. Wald/ratio estimator: θ = cov(Δw, D)/cov(Δl, D) (simple)
        2. 2SLS: Δw = θ·Δl + year FEs, instrument Δl with D
        3. Panel levels: w = θ·le + firm FE + year FE, instrument le with D

        Note: Czech employment data comes from categorical brackets, so
        first-differenced approaches have limited power. The levels approach
        (3) exploits cross-sectional variation instead.
        """
        print('\n=== Analysis 2: Labor Supply Elasticity θ ===')

        # --- Approach A: First-differenced Wald ---
        df_fd = self.df.dropna(subset=['d_w', 'd_le']).copy()
        # Remove zero Δl (bracket artifacts)
        df_fd_nonzero = df_fd[df_fd['d_le'].abs() > 0.01]
        print(f'  First-differenced sample: {len(df_fd):,} obs '
              f'(non-zero Δl: {len(df_fd_nonzero):,})')

        dw_1 = df_fd.loc[df_fd['D'] == 1, 'd_w'].mean()
        dw_0 = df_fd.loc[df_fd['D'] == 0, 'd_w'].mean()
        dl_1 = df_fd.loc[df_fd['D'] == 1, 'd_le'].mean()
        dl_0 = df_fd.loc[df_fd['D'] == 0, 'd_le'].mean()

        cov_dw_D = np.cov(df_fd['d_w'], df_fd['D'])[0, 1]
        cov_dl_D = np.cov(df_fd['d_le'], df_fd['D'])[0, 1]
        theta_wald = cov_dw_D / cov_dl_D if abs(cov_dl_D) > 1e-10 else np.nan
        print(f'  Wald (first-diff): θ = {theta_wald:.4f}' if not np.isnan(theta_wald)
              else '  Wald (first-diff): θ = NaN (Δl ≈ 0 — bracket-based employment)')
        print(f'    E[Δw|D=1]={dw_1:.4f}, E[Δw|D=0]={dw_0:.4f}, '
              f'E[Δl|D=1]={dl_1:.4f}, E[Δl|D=0]={dl_0:.4f}')

        # --- Approach B: Panel levels with FE-IV ---
        # w = θ·le + firm FE + year FE, instrument le with pp_dummy
        # This exploits cross-sectional variation in employment brackets
        df_lev = self.df.dropna(subset=['w', 'le']).copy()
        df_lev = df_lev.set_index(['id', 'year'])
        theta_panel = np.nan
        theta_panel_se = np.nan
        first_stage_f = np.nan

        try:
            mod = PanelOLS(df_lev['w'], df_lev[['le']],
                          entity_effects=True, time_effects=True)
            res_ols = mod.fit(cov_type='clustered', cluster_entity=True)
            theta_ols = res_ols.params['le']
            print(f'\n  Panel FE OLS (levels): θ_OLS = {theta_ols:.4f} '
                  f'(SE {res_ols.std_errors["le"]:.4f}), N={res_ols.nobs:,}')
        except Exception as exc:
            theta_ols = np.nan
            print(f'  Panel FE OLS failed: {exc}')

        # IV version: instrument le with pp_dummy (procurement shifts labor demand)
        try:
            from linearmodels.panel import PanelOLS as PanelOLS2
            from linearmodels.iv import IV2SLS as IV2SLS2

            # Use cross-section IV: levels with firm dummies
            df_iv = self.df.dropna(subset=['w', 'le']).copy()
            # Demean within firm (manual FE)
            for v in ['w', 'le']:
                df_iv[f'{v}_dm'] = df_iv[v] - df_iv.groupby('id')[v].transform('mean')
            df_iv['D_dm'] = df_iv['D'] - df_iv.groupby('id')['D'].transform('mean')

            year_dummies = pd.get_dummies(df_iv['year'], prefix='yr',
                                          drop_first=True, dtype=float)
            exog = pd.concat([year_dummies], axis=1)

            iv_mod = IV2SLS(
                dependent=df_iv['w_dm'],
                exog=exog,
                endog=df_iv[['le_dm']],
                instruments=df_iv[['D_dm']]
            )
            iv_res = iv_mod.fit(cov_type='clustered', clusters=df_iv['id'])
            theta_panel = iv_res.params['le_dm']
            theta_panel_se = iv_res.std_errors['le_dm']
            first_stage_f = iv_res.first_stage.diagnostics['le_dm']['f.stat']

            print(f'  Panel FE-IV (levels): θ_IV = {theta_panel:.4f} '
                  f'(SE {theta_panel_se:.4f})')
            print(f'    First-stage F = {first_stage_f:.1f}')
            if abs(theta_panel) > 0.01:
                print(f'    LS elasticity 1/θ = {1/theta_panel:.2f}')
                print(f'    Wage markdown 1/(1+θ) = {1/(1+theta_panel):.3f}')
        except Exception as exc:
            print(f'  Panel FE-IV failed: {exc}')

        # Best estimate — θ must be positive (upward-sloping labor supply)
        theta_best = np.nan
        for candidate, label in [
            (theta_panel, 'panel IV'),
            (theta_ols, 'panel OLS'),
            (theta_wald, 'Wald'),
        ]:
            if not np.isnan(candidate) and candidate > 0:
                theta_best = candidate
                print(f'\n  Best θ estimate: {theta_best:.4f} (from {label})')
                break
        if np.isnan(theta_best):
            print('\n  No valid positive θ estimate — Czech bracket-based '
                  'employment data lacks variation for KLMS-style identification')

        result = {
            'theta_wald': theta_wald,
            'theta_ols': theta_ols if 'theta_ols' in dir() else np.nan,
            'theta_panel_iv': theta_panel,
            'theta_panel_iv_se': theta_panel_se,
            'theta_2sls': theta_panel,  # backward compat
            'theta_2sls_se': theta_panel_se,
            'first_stage_f': first_stage_f,
            'theta_best': theta_best,
            'ls_elasticity': 1 / theta_best if abs(theta_best) > 0.01 else np.nan,
            'wage_markdown': 1 / (1 + theta_best) if not np.isnan(theta_best) else np.nan,
            'n_obs': len(df_lev),
            'n_firms': df_lev.index.get_level_values(0).nunique(),
            'dw_1': dw_1, 'dw_0': dw_0, 'dl_1': dl_1, 'dl_0': dl_0,
        }
        self.results['theta'] = result
        return result

    # ================================================================== #
    #  Analysis 3: Revenue Elasticity (1-ε)
    # ================================================================== #

    def estimate_one_minus_epsilon(self) -> Dict[str, float]:
        """Estimate (1-ε) from KLMS eq. 25: go ~ cogs among D=0 firms.

        Revenue elasticity of output in the private product market.
        """
        print('\n=== Analysis 3: Revenue Elasticity (1-ε) ===')

        df = self.df.copy()

        results = {}

        # --- Baseline: OLS among non-procurement firm-years ---
        for label, sample_mask, fe in [
            ('Pooled OLS, D=0', df['D'] == 0, False),
            ('Pooled OLS, all', pd.Series(True, index=df.index), False),
            ('Firm FE, D=0', df['D'] == 0, True),
        ]:
            sub = df[sample_mask].dropna(subset=['go', 'cogs']).copy()
            if len(sub) < 50:
                continue

            if fe:
                sub = sub.set_index(['id', 'year'])
                mod = PanelOLS(sub['go'], sub[['cogs']],
                               entity_effects=True, time_effects=True)
                res = mod.fit(cov_type='clustered', cluster_entity=True)
            else:
                import statsmodels.api as sm
                X = sm.add_constant(sub['cogs'])
                mod = sm.OLS(sub['go'], X)
                res = mod.fit(cov_type='cluster', cov_kwds={'groups': sub['id']})

            one_minus_eps = res.params['cogs']
            se = res.std_errors['cogs'] if hasattr(res, 'std_errors') else res.bse['cogs']

            eps = 1 - one_minus_eps
            demand_elast = -1 / eps if abs(eps) > 0.01 else np.nan
            markup = 1 / one_minus_eps if one_minus_eps > 0.01 else np.nan

            print(f'  {label}: (1-ε) = {one_minus_eps:.4f} (SE {se:.4f}), '
                  f'N = {res.nobs:,}')
            print(f'    ε = {eps:.4f}, demand elasticity = {demand_elast:.2f}, '
                  f'price markup (1-ε)⁻¹ = {markup:.3f}')

            results[label] = {
                'one_minus_eps': one_minus_eps, 'se': se,
                'epsilon': eps, 'demand_elasticity': demand_elast,
                'price_markup': markup, 'n_obs': int(res.nobs),
            }

        # --- By NACE2 ---
        print('\n  By NACE2 (OLS, D=0):')
        for nace in [41, 42, 43]:
            sub = df[(df['D'] == 0) & (df['nace2'] == nace)].dropna(
                subset=['go', 'cogs']).copy()
            if len(sub) < 30:
                continue
            import statsmodels.api as sm
            X = sm.add_constant(sub['cogs'])
            res = sm.OLS(sub['go'], X).fit(cov_type='cluster',
                                           cov_kwds={'groups': sub['id']})
            ome = res.params['cogs']
            se = res.bse['cogs']
            print(f'    NACE {nace}: (1-ε) = {ome:.4f} (SE {se:.4f}), N = {len(sub)}')
            results[f'NACE {nace}'] = {
                'one_minus_eps': ome, 'se': se,
                'epsilon': 1 - ome, 'n_obs': len(sub),
            }

        self.results['one_minus_eps'] = results
        return results

    # ================================================================== #
    #  Analysis 4: Labor Share β_L
    # ================================================================== #

    def estimate_beta_L(self, theta: float, one_minus_eps: float) -> Dict[str, float]:
        """Estimate β_L from KLMS eq. 24: labor share formula.

        β_L = (1+θ)·s_L / ((1-ε) - s_M)

        Also estimates ρ from intermediates regression and derives β_K.
        """
        print('\n=== Analysis 4: Production Function Parameters ===')

        df = self.df[(self.df['D'] == 0)].dropna(
            subset=['s_L', 's_M', 'go', 'cogs']).copy()

        # Labor share formula (KLMS eq. 24)
        # β_L = (1+θ)·s_L / ((1-ε) - s_M)  — winsorized mean
        inner = (1 + theta) * df['s_L'] / (one_minus_eps - df['s_M'])
        inner_w = winsorize(inner.replace([np.inf, -np.inf], np.nan).dropna())
        beta_L = inner_w.mean()

        # Composite returns ρ from intermediates regression (KLMS eq. 11)
        # x = κ + ρ·ℓ + φ  (among firms with employment data)
        rho = np.nan
        rho_se = np.nan
        df_rho = self.df.dropna(subset=['cogs', 'le']).copy()
        if len(df_rho) > 50:
            import statsmodels.api as sm
            X = sm.add_constant(df_rho['le'])
            res = sm.OLS(df_rho['cogs'], X).fit(
                cov_type='cluster', cov_kwds={'groups': df_rho['id']})
            rho = res.params['le']
            rho_se = res.bse['le']
            print(f'  Composite returns ρ: {rho:.4f} (SE {rho_se:.4f}), N = {len(df_rho)}')

        # Derive β_K = (ρ - β_L)/(1+θ)
        beta_K = (rho - beta_L) / (1 + theta)
        returns_to_scale = beta_L + beta_K

        print(f'  β_L (labor output elasticity): {beta_L:.4f}')
        print(f'  β_K (capital share): {beta_K:.4f}')
        print(f'  Returns to scale (β_L + β_K): {returns_to_scale:.3f}')
        print(f'  [KLMS: β_L = 0.499, β_K = 0.474, RTS = 0.97]')

        result = {
            'beta_L': beta_L, 'beta_K': beta_K,
            'rho': rho, 'rho_se': rho_se,
            'returns_to_scale': returns_to_scale,
            'theta_used': theta, 'one_minus_eps_used': one_minus_eps,
        }
        self.results['beta_L'] = result
        return result

    # ================================================================== #
    #  Analysis 5: Double Markdown/Markup
    # ================================================================== #

    def double_markdown(self, theta: float, epsilon: float) -> pd.DataFrame:
        """Compute KLMS Table 3: double markdown/markup decomposition.

        Wage double markdown:  (1+θ)⁻¹ × (1-ε)
        Price double markup:   (1-ε)⁻¹ × (1+θ)
        """
        print('\n=== Analysis 5: Double Markdown/Markup (KLMS Table 3) ===')

        one_minus_eps = 1 - epsilon

        # Wage side
        markdown = 1 / (1 + theta)              # labor-only markdown
        inv_markup = one_minus_eps               # product-side contribution
        double_markdown = markdown * inv_markup  # combined

        # Price side
        markup = 1 / one_minus_eps               # product-only markup
        inv_markdown = 1 + theta                 # labor-side contribution
        double_markup = markup * inv_markdown    # combined

        rows = []

        # Panel A: wage
        rows.append({
            'panel': 'A. Wage markdown',
            'Markdown (1+θ)⁻¹': markdown,
            'Inverse markup (1-ε)': inv_markup,
            'Double markdown': double_markdown,
        })

        # Panel B: price
        rows.append({
            'panel': 'B. Price markup',
            'Markup (1-ε)⁻¹': markup,
            'Inverse markdown (1+θ)': inv_markdown,
            'Double markup': double_markup,
        })

        result_df = pd.DataFrame(rows).set_index('panel')

        print(f'\n  Panel A — Wage:')
        print(f'    Markdown (1+θ)⁻¹   = {markdown:.3f}')
        print(f'    Inverse markup (1-ε) = {inv_markup:.3f}')
        print(f'    Double markdown      = {double_markdown:.3f}')
        print(f'    → Wages are {(1 - double_markdown)*100:.1f}% below value of MPL')

        print(f'\n  Panel B — Price:')
        print(f'    Markup (1-ε)⁻¹      = {markup:.3f}')
        print(f'    Inverse markdown (1+θ) = {inv_markdown:.3f}')
        print(f'    Double markup        = {double_markup:.3f}')
        print(f'    → Prices are {(double_markup - 1)*100:.1f}% above productivity-adjusted wage')

        # Comparison: what if you ignored one market?
        print(f'\n  Counterfactual:')
        print(f'    If ε=0 (no product power): wage markdown = {markdown:.3f} '
              f'({(1-markdown)*100:.1f}% below VMPL)')
        print(f'    If θ=0 (no labor power):   price markup  = {markup:.3f} '
              f'({(markup-1)*100:.1f}% above MC)')
        print(f'    With both: wage {(1-double_markdown)*100:.1f}% below, '
              f'price {(double_markup-1)*100:.1f}% above')

        # KLMS comparison
        print(f'\n  KLMS (2025) comparison:')
        print(f'    KLMS double markdown = 0.693, ours = {double_markdown:.3f}')
        print(f'    KLMS double markup   = 1.443, ours = {double_markup:.3f}')

        self.results['double_markdown'] = result_df
        return result_df

    # ================================================================== #
    #  Analysis 6: Crowd-Out Test
    # ================================================================== #

    def crowd_out_test(self) -> pd.DataFrame:
        """Test whether procurement crowds out private-market production.

        KLMS find 1+θ > ρ → crowd-out of 27% of private output.
        """
        print('\n=== Analysis 6: Crowd-Out Test ===')

        df = self.df.copy()
        # Compare private revenue around procurement entry
        mask = df['event_time'].notna() | (df['never_treated'] == 1)
        df = df[mask].dropna(subset=['private_rev']).copy()

        if len(df) < 100:
            print('  Insufficient data for crowd-out test')
            return pd.DataFrame()

        # Simple comparison: private revenue for treated vs untreated
        treated = df[df['D'] == 1]
        untreated = df[df['D'] == 0]

        print(f'  Mean log(private rev) — treated:   {treated["private_rev"].mean():.3f} '
              f'(N={len(treated)})')
        print(f'  Mean log(private rev) — untreated: {untreated["private_rev"].mean():.3f} '
              f'(N={len(untreated)})')

        # Event study of private revenue (same design as Analysis 1)
        results_df = self.event_study(
            outcomes=['private_rev'],
            window=(-3, 4), ref_period=-1
        )
        self.results['crowd_out'] = results_df
        return results_df

    # ================================================================== #
    #  Plotting
    # ================================================================== #

    def plot_event_study(self, results_df: pd.DataFrame,
                         save_path: Optional[Path] = None) -> None:
        """Plot event study coefficients with 95% CIs."""
        outcomes = results_df['outcome'].unique()
        n = len(outcomes)
        fig, axes = plt.subplots(1, min(n, 3), figsize=(5 * min(n, 3), 4),
                                 squeeze=False)
        axes = axes.flatten()

        labels = {
            'le': 'Log employment', 'w': 'Log wage bill',
            'avg_wage': 'Log avg wage', 'go': 'Log output',
            'cogs': 'Log COGS', 'private_rev': 'Log private revenue',
        }

        for i, outcome in enumerate(outcomes[:3]):
            ax = axes[i]
            sub = results_df[results_df['outcome'] == outcome].sort_values('event_time')
            ax.fill_between(sub['event_time'], sub['ci_lo'], sub['ci_hi'],
                           alpha=0.2, color='steelblue')
            ax.plot(sub['event_time'], sub['coef'], 'o-', color='steelblue',
                   markersize=4)
            ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
            ax.axvline(-0.5, color='red', linewidth=0.5, linestyle=':')
            ax.set_xlabel('Years since procurement entry')
            ax.set_ylabel('Coefficient')
            ax.set_title(labels.get(outcome, outcome))

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, bbox_inches='tight', dpi=150)
            print(f'  Saved: {save_path}')
        plt.close()

    # ================================================================== #
    #  Summary Table
    # ================================================================== #

    def summary_table(self) -> pd.DataFrame:
        """Compile all structural parameter estimates into a summary table."""
        rows = []

        # θ
        if 'theta' in self.results:
            t = self.results['theta']
            theta_val = t.get('theta_best', np.nan)
            theta_se = t.get('theta_panel_iv_se', np.nan)
            theta_str = f"{theta_val:.3f}" if not np.isnan(theta_val) else 'N/A*'
            se_str = f"{theta_se:.3f}" if not np.isnan(theta_se) else ''
            rows.append({
                'Parameter': 'θ (inv. LS elasticity)',
                'Czech estimate': theta_str,
                'SE': se_str,
                'KLMS (2025)': '0.245',
                'KLMS SE': '0.086',
            })
            ls_val = t.get('ls_elasticity', np.nan)
            rows.append({
                'Parameter': '1/θ (LS elasticity)',
                'Czech estimate': f"{ls_val:.2f}" if not np.isnan(ls_val) else 'N/A*',
                'SE': '',
                'KLMS (2025)': '4.08',
                'KLMS SE': '',
            })
            md_val = t.get('wage_markdown', np.nan)
            rows.append({
                'Parameter': 'Wage markdown 1/(1+θ)',
                'Czech estimate': f"{md_val:.3f}" if not np.isnan(md_val) else 'N/A*',
                'SE': '',
                'KLMS (2025)': '0.803',
                'KLMS SE': '0.055',
            })

        # ε — prefer firm FE estimate
        if 'one_minus_eps' in self.results:
            baseline = (self.results['one_minus_eps'].get('Firm FE, D=0')
                       or self.results['one_minus_eps'].get('Pooled OLS, D=0', {}))
            if baseline:
                rows.append({
                    'Parameter': '1-ε (revenue elasticity)',
                    'Czech estimate': f"{baseline['one_minus_eps']:.3f}",
                    'SE': f"{baseline['se']:.3f}",
                    'KLMS (2025)': '0.863',
                    'KLMS SE': '0.015',
                })
                pm = baseline.get('price_markup', np.nan)
                rows.append({
                    'Parameter': 'Price markup (1-ε)⁻¹',
                    'Czech estimate': f"{pm:.3f}" if not np.isnan(pm) else '',
                    'SE': '',
                    'KLMS (2025)': '1.159',
                    'KLMS SE': '',
                })

        # β_L, β_K
        if 'beta_L' in self.results:
            b = self.results['beta_L']
            rows.append({
                'Parameter': 'β_L (labor elasticity)',
                'Czech estimate': f"{b['beta_L']:.3f}",
                'SE': '',
                'KLMS (2025)': '0.499',
                'KLMS SE': '0.192',
            })
            rows.append({
                'Parameter': 'ρ (composite returns)',
                'Czech estimate': f"{b['rho']:.3f}",
                'SE': f"{b['rho_se']:.3f}",
                'KLMS (2025)': '1.089',
                'KLMS SE': '0.017',
            })

        df = pd.DataFrame(rows)
        self.results['summary'] = df
        return df

    # ================================================================== #
    #  Run All
    # ================================================================== #

    def run_all(self) -> None:
        """Execute all analyses and save outputs."""
        print('=' * 60)
        print('KLMS-Inspired Double Market Power Analysis')
        print('Czech Construction (MagnusWeb + Datlab)')
        print('=' * 60)

        # 1. Event study
        evt = self.event_study(
            outcomes=['le', 'w', 'avg_wage', 'go', 'cogs'])
        if not evt.empty:
            evt.to_csv(TABLE_DIR / 'klms_event_study.csv', index=False)
            self.plot_event_study(evt, FIG_DIR / 'klms_event_study.pdf')

        # 2. Labor supply elasticity
        theta_res = self.estimate_theta()
        theta = theta_res.get('theta_best', np.nan)
        if np.isnan(theta) or abs(theta) < 0.01:
            print('  WARNING: θ estimate unreliable, using KLMS value (0.245) '
                  'for downstream calculations')
            theta = 0.245

        # 3. Revenue elasticity
        # Prefer firm FE estimate (removes unobserved firm heterogeneity)
        eps_res = self.estimate_one_minus_epsilon()
        fe_eps = eps_res.get('Firm FE, D=0', {})
        ols_eps = eps_res.get('Pooled OLS, D=0', {})
        baseline_eps = fe_eps if fe_eps else ols_eps
        one_minus_eps = baseline_eps.get('one_minus_eps', 0.863)
        epsilon = 1 - one_minus_eps
        print(f'\n  Using (1-ε) = {one_minus_eps:.4f} '
              f'({"firm FE" if fe_eps else "pooled OLS"} estimate)')

        # 4. β_L
        self.estimate_beta_L(theta, one_minus_eps)

        # 5. Double markdown/markup
        dm = self.double_markdown(theta, epsilon)
        dm.to_csv(TABLE_DIR / 'klms_double_markdown.csv')

        # 6. Crowd-out
        crowd = self.crowd_out_test()
        if not crowd.empty:
            crowd.to_csv(TABLE_DIR / 'klms_crowd_out.csv', index=False)

        # Summary
        summary = self.summary_table()
        summary.to_csv(TABLE_DIR / 'klms_summary.csv', index=False)
        print('\n' + '=' * 60)
        print('SUMMARY OF STRUCTURAL PARAMETERS')
        print('=' * 60)
        print(summary.to_string(index=False))

        print(f'\nOutputs saved to:')
        print(f'  Tables: {TABLE_DIR}/klms_*.csv')
        print(f'  Figures: {FIG_DIR}/klms_*.pdf')


# ========================================================================== #
#  Main
# ========================================================================== #

def write_theta_latex(results_mw, results_orbis):
    """Write LaTeX table for KLMS theta appendix (replaces hardcoded table)."""
    fn = TABLE_DIR / 'klms_theta_appendix.tex'

    def _fmt(v, fmt='%.3f'):
        return fmt % v if not np.isnan(v) else '---'

    lines = [
        r'\begin{table}[htbp]\centering',
        r'\caption{Labor Supply Elasticity Estimates}\label{tab:klms_theta}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lcccc}',
        r'\toprule',
        r'& $\hat{\theta}$ & $1/\theta$ & Markdown & $N$ \\',
        r'\midrule',
    ]

    # Orbis results
    if results_orbis:
        th = results_orbis.get('theta', {})
        lines.append(r'\multicolumn{5}{l}{\emph{Czech Republic (Orbis, continuous employment)}} \\')

        theta_w = th.get('theta_wald', np.nan)
        lines.append(
            r'\quad Construction (Wald DiD) & '
            f'{_fmt(theta_w)} & {_fmt(1/theta_w, "%.2f") if abs(theta_w) > 0.01 else "---"} & '
            f'{_fmt(1/(1+theta_w))} & {th.get("n_obs", ""):,} \\\\'
        )

        theta_iv = th.get('theta_panel_iv', np.nan)
        lines.append(
            r'\quad Construction (Panel FE-IV) & '
            f'{_fmt(theta_iv)} & {_fmt(1/theta_iv, "%.2f") if abs(theta_iv) > 0.01 else "---"} & '
            f'{_fmt(1/(1+theta_iv))} & {th.get("n_obs", ""):,} \\\\[4pt]'
        )

    # MagnusWeb results
    if results_mw:
        th = results_mw.get('theta', {})
        lines.append(r'\multicolumn{5}{l}{\emph{Czech Republic (MagnusWeb, bracket employment)}} \\')
        theta_w = th.get('theta_wald', np.nan)
        lines.append(
            r'\quad Construction (Wald DiD) & '
            f'{_fmt(theta_w) if not np.isnan(theta_w) else "N/A$^*$"} & '
            f'{_fmt(1/theta_w, "%.2f") if not np.isnan(theta_w) and abs(theta_w) > 0.01 else "---"} & '
            f'{_fmt(1/(1+theta_w)) if not np.isnan(theta_w) else "---"} & '
            f'{th.get("n_obs", ""):,} \\\\[4pt]'
        )

    # KLMS US benchmarks
    lines += [
        r'\multicolumn{5}{l}{\emph{United States (KLMS 2025)}} \\',
        r'\quad Construction (DiD) & 0.245 & 4.08 & 0.803 & --- \\',
        r'\quad Construction (RDD) & 0.286 & 3.50 & 0.777 & --- \\',
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r'\item \emph{Notes:} Wald DiD: $\hat{\theta} = (\bar{\Delta w}_1 - \bar{\Delta w}_0)/(\bar{\Delta \ell}_1 - \bar{\Delta \ell}_0)$. Panel FE-IV instruments log employment with procurement entry (firm and year FE). KLMS estimates from Table~2 of Kroft et al.\ (2025). $^*$MagnusWeb bracket employment produces near-zero $\Delta\ell$, precluding identification.',
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}',
    ]

    with open(fn, 'w') as f:
        f.write('\n'.join(lines))
    print(f'  Saved: {fn}')


if __name__ == '__main__':
    import sys

    # Run MagnusWeb baseline
    print('=' * 60)
    print('  KLMS Analysis — MagnusWeb (baseline)')
    print('=' * 60)
    mw = KLMSAnalysis(use_orbis=False)
    mw.run_all()
    results_mw = mw.results

    # Run Orbis (continuous employment)
    if ORBIS_DTA.exists():
        print('\n' + '=' * 60)
        print('  KLMS Analysis — Orbis (continuous employment)')
        print('=' * 60)
        orb = KLMSAnalysis(use_orbis=True)
        orb.run_all()
        results_orbis = orb.results

        # Save Orbis-specific outputs
        orbis_summary = orb.summary_table()
        orbis_summary.to_csv(TABLE_DIR / 'klms_summary_orbis.csv', index=False)
        print(f'\n  Saved: klms_summary_orbis.csv')

        # Generate appendix LaTeX table
        write_theta_latex(results_mw, results_orbis)
    else:
        print(f'\n  Orbis data not found at {ORBIS_DTA} — skipping')
        results_orbis = None
        write_theta_latex(results_mw, None)
