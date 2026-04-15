"""
Paper Results: Production Function Estimates and Procurement Premium
====================================================================

Produces publication-ready results for the markup estimation paper:
  - Table 1: PF coefficients by NACE (baseline + robustness)
  - Table 2: Procurement premium across specifications
  - Firm-level markups saved for treatment effect estimation

Baseline: CD by nace2, survival + pp_dummy in Markov (CWDL 2015)
"""

import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings('ignore')

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from acf_estimator import (
    ACFEstimator, Formulation, Optimization, CWDLExtensions,
    estimate_by_industry, options
)

options.verbose = False

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

# ================================================================ #
#  Data loading and survival construction
# ================================================================ #

DATA_PATH = str(INPUT_DIR / 'data_rebuilt.dta')
OUT_DIR = str(OUTPUT_DIR)


def construct_survival(df: pd.DataFrame) -> pd.DataFrame:
    """Probit-predicted survival probability (CWDL 2015)."""
    df = df.sort_values(['id', 'year']).copy()
    df['next_yr'] = df.groupby('id')['year'].shift(-1)
    df['survival_1'] = ((df['next_yr'] - df['year']) == 1).astype(float)

    probit_data = df[['survival_1', 'k', 'cogs', 'pp_dummy',
                       'year', 'nace2']].dropna()
    X_cols = ['k', 'cogs', 'pp_dummy']
    X = probit_data[X_cols].copy()
    yr_nace = (probit_data['year'].astype(str) + '_'
               + probit_data['nace2'].astype(str))
    for val in sorted(yr_nace.unique())[1:]:
        X[f'fe_{val}'] = (yr_nace == val).astype(float)
    y = probit_data['survival_1'].values
    model = LogisticRegression(max_iter=1000, solver='lbfgs', penalty=None)
    model.fit(X.values, y)
    phat = model.predict_proba(X.values)[:, 1]
    df.loc[probit_data.index, 'survival'] = phat
    return df


def compute_premium(markups_df, df_orig, treatment_var='pp_dummy'):
    """Compute log markup premium: E[log μ | pp=1] - E[log μ | pp=0]."""
    mu = markups_df[markups_df['markup'] > 0].copy()
    mu['lmu'] = np.log(mu['markup'])
    mu = mu.merge(
        df_orig[['id', 'year', treatment_var]].drop_duplicates(),
        on=['id', 'year'], how='left', suffixes=('', '_r')
    )
    if f'{treatment_var}_r' in mu.columns:
        mu[treatment_var] = mu[f'{treatment_var}_r']

    pp0 = mu[mu[treatment_var] == 0]['lmu']
    pp1 = mu[mu[treatment_var] == 1]['lmu']
    if len(pp0) == 0 or len(pp1) == 0:
        return np.nan, np.nan, 0, 0
    raw = pp1.mean() - pp0.mean()
    se = np.sqrt(pp1.var()/len(pp1) + pp0.var()/len(pp0))
    return raw, se, len(pp0), len(pp1)


def compute_regression_premium(markups_df, df_orig, treatment_var='pp_dummy'):
    """Premium from regression: log(μ) ~ pp_dummy + year×nace2 FE."""
    mu = markups_df[markups_df['markup'] > 0].copy()
    mu['lmu'] = np.log(mu['markup'])
    mu = mu.merge(
        df_orig[['id', 'year', treatment_var, 'nace2', 'k', 'cogs']].drop_duplicates(),
        on=['id', 'year'], how='left', suffixes=('', '_r')
    )
    if f'{treatment_var}_r' in mu.columns:
        mu[treatment_var] = mu[f'{treatment_var}_r']
    if 'nace2_r' in mu.columns:
        mu['nace2'] = mu['nace2_r']

    # Build design matrix: pp_dummy + k + cogs + year×nace2 FE
    mu = mu.dropna(subset=['lmu', treatment_var, 'k', 'cogs'])
    yr_nace = mu['year'].astype(str) + '_' + mu['nace2'].astype(str)
    X = pd.DataFrame({'const': 1.0, 'pp': mu[treatment_var].values,
                       'k': mu['k'].values, 'cogs': mu['cogs'].values},
                      index=mu.index)
    for val in sorted(yr_nace.unique())[1:]:
        X[f'fe_{val}'] = (yr_nace == val).astype(float)

    y = mu['lmu'].values
    Xm = X.values.astype(np.float64)
    bhat = np.linalg.lstsq(Xm, y, rcond=None)[0]
    resid = y - Xm @ bhat
    N = len(y)
    k = Xm.shape[1]
    # HC1 robust SE for pp coefficient (index 1)
    XtX_inv = np.linalg.inv(Xm.T @ Xm)
    meat = Xm.T @ np.diag(resid**2) @ Xm
    V = XtX_inv @ meat @ XtX_inv * N / (N - k)
    se_pp = np.sqrt(V[1, 1])
    r2 = 1 - np.sum(resid**2) / np.sum((y - y.mean())**2)
    return bhat[1], se_pp, r2, N


# ================================================================ #
#  Specification definitions
# ================================================================ #

SPECS = {
    'A': {
        'label': 'Base translog (survival + pp in Markov)',
        'short': 'Translog',
        'formulation': Formulation(spec='tl', overidentify=True,
                                   pp_in_markov=True),
        'extensions': CWDLExtensions(
            survival_correction=True,
            markov_interactions=True,
        ),
    },
    'B': {
        'label': 'No survival correction',
        'short': 'No surv.',
        'formulation': Formulation(spec='cd', pp_in_markov=True),
        'extensions': CWDLExtensions(
            survival_correction=False,
            markov_interactions=True,
        ),
    },
    'C': {
        'label': 'No pp in Markov',
        'short': 'No pp Markov',
        'formulation': Formulation(spec='cd', pp_in_markov=False,
                                   pp_interactions=False),
        'extensions': CWDLExtensions(
            survival_correction=True,
            markov_interactions=False,
        ),
    },
    'D': {
        'label': 'Plain ACF (CD)',
        'short': 'Plain CD',
        'formulation': Formulation(spec='cd', pp_in_markov=False,
                                   pp_interactions=False),
        'extensions': CWDLExtensions(),
    },
    'E': {
        'label': 'Cobb-Douglas (survival + pp)',
        'short': 'CD base',
        'formulation': Formulation(spec='cd', pp_in_markov=True),
        'extensions': CWDLExtensions(
            survival_correction=True,
            markov_interactions=True,
        ),
    },
    # Spec G (per-NACE two-step) was implemented and tested 2026-04-15 but
    # per-NACE N_c is too thin for the de-meaned S_hat inverse (NACE 42 had
    # delta_step = 11.35 — Conlon's "too many moments" regime, gmmnotes §
    # Common Questions). Replaced by pooled Spec H below.
}


# ================================================================ #
#  Main
# ================================================================ #

def main():
    print(f'Loading data from {DATA_PATH}')
    df = pd.read_stata(DATA_PATH)
    print(f'Data: {len(df)} obs, {df["id"].nunique()} firms, '
          f'years {df["year"].min()}-{df["year"].max()}')

    # Survival
    df = construct_survival(df)

    # Treatment variables
    has_pp_ever_3y = 'pp_ever_3y' in df.columns
    if not has_pp_ever_3y:
        df['pp_dummy_L1'] = df.groupby('id')['pp_dummy'].shift(1).fillna(0)
        df['pp_dummy_L2'] = df.groupby('id')['pp_dummy'].shift(2).fillna(0)
        df['pp_ever_3y'] = ((df['pp_dummy'] == 1) |
                            (df['pp_dummy_L1'] == 1) |
                            (df['pp_dummy_L2'] == 1)).astype(int)

    naces = sorted(df['nace2'].unique())
    all_results = []
    all_markups = []

    for spec_key, spec_cfg in SPECS.items():
        print(f'\n{"="*60}')
        print(f'  Spec {spec_key}: {spec_cfg["label"]}')
        print(f'{"="*60}')

        for nace in naces:
            df_n = df[df['nace2'] == nace].copy()
            try:
                est = ACFEstimator(
                    data=df_n,
                    formulation=spec_cfg['formulation'],
                    optimization=Optimization(method='nm+bfgs'),
                    extensions=spec_cfg['extensions'],
                )
                res = est.solve()

                row = {
                    'spec': spec_key,
                    'label': spec_cfg['label'],
                    'short': spec_cfg['short'],
                    'nace2': int(nace),
                    'N': res.n_obs,
                }
                for name, coef, se in zip(res.beta_names, res.betas, res.se):
                    row[f'b_{name}'] = coef
                    row[f'se_{name}'] = se

                md = res.data
                row['markup_mean'] = md['markup'].mean()
                row['markup_sd'] = md['markup'].std()
                row['markup_p10'] = md['markup'].quantile(0.1)
                row['markup_p50'] = md['markup'].quantile(0.5)
                row['markup_p90'] = md['markup'].quantile(0.9)
                row['criterion'] = res.gmm_criterion
                row['r2_first'] = res.first_stage_r2

                # Hansen J test (Conlon gmmnotes lines 141-148; required by
                # .claude/rules/abgrs-transparency.md rule #4 for overid specs).
                # The estimator already populates these in solve(); we just plumb
                # them through to the output table so paper §6.1 / line 733 can
                # cite the Python pipeline numbers instead of relying solely on
                # the Stata-side table_misspec_diagnostics.do path.
                row['hansen_j'] = res.hansen_j if res.hansen_j is not None else np.nan
                row['hansen_j_p'] = (res.hansen_j_pvalue
                                      if res.hansen_j_pvalue is not None else np.nan)
                row['n_overid'] = res.n_overid

                # Markov transition parameters
                if res.markov_coefs is not None:
                    for mn, mc in zip(res.markov_names, res.markov_coefs):
                        row[f'markov_{mn}'] = mc

                all_results.append(row)

                keep_cols = ['id', 'year', 'markup']
                if 'omega' in md.columns:
                    keep_cols.append('omega')
                if 'alphahat' in md.columns:
                    keep_cols.append('alphahat')
                mu = md[keep_cols].copy()
                mu['spec'] = spec_key
                mu['nace2'] = int(nace)
                all_markups.append(mu)

                b_k = row.get('b_k', np.nan)
                b_c = row.get('b_cogs', np.nan)
                se_k = row.get('se_k', np.nan)
                se_c = row.get('se_cogs', np.nan)
                rts = b_k + b_c if np.isfinite(b_k) and np.isfinite(b_c) else np.nan
                print(f'  NACE {int(nace)}: β_k={b_k:.4f}({se_k:.4f}), '
                      f'β_cogs={b_c:.4f}({se_c:.4f}), '
                      f'RTS={rts:.3f}, μ̄={row["markup_mean"]:.3f}, '
                      f'N={res.n_obs}')
                if res.n_overid > 0 and res.hansen_j is not None:
                    print(f'    [Hansen J] χ²({res.n_overid}) = {res.hansen_j:.3f}, '
                          f'p = {res.hansen_j_pvalue:.4f}')
                # Wald test on translog cross-terms: H0: β_k²=β_cogs²=β_kcogs=0
                # (CD is nested inside TL). Rejection → TL preferred;
                # fail to reject → CD is adequate. Gives a clean statistical
                # justification for the NACE-42-CD branch of Spec J below.
                if (spec_cfg['formulation'].spec == 'tl'
                        and all(nm in res.beta_names for nm in ('k2','cogs2','kcogs'))):
                    from scipy.stats import chi2 as _chi2
                    idx = [res.beta_names.index(nm) for nm in ('k2','cogs2','kcogs')]
                    b_cross = res.betas[idx]
                    V_cross = res.vcov[np.ix_(idx, idx)]
                    try:
                        V_inv = np.linalg.inv(V_cross)
                    except np.linalg.LinAlgError:
                        V_inv = np.linalg.pinv(V_cross)
                    wald_cd = float(b_cross @ V_inv @ b_cross)
                    wald_cd_p = float(1 - _chi2.cdf(wald_cd, df=3))
                    row['wald_cd'] = wald_cd
                    row['wald_cd_p'] = wald_cd_p
                    print(f'    [Wald TL vs CD cross-terms] χ²(3) = '
                          f'{wald_cd:.3f}, p = {wald_cd_p:.4f}')
                if res.weighting == 'two_step' and res.delta_step is not None:
                    b1 = res.betas_step1
                    print(f'    [two-step] delta_step (max |Δβ|) = {res.delta_step:.6f}')
                    if b1 is not None and len(b1) >= 3:
                        print(f'    [two-step] step1 β_cogs = {b1[2]:.4f}, '
                              f'step2 β_cogs = {res.betas[2]:.4f}, '
                              f'Δ = {res.betas[2] - b1[2]:+.4f}')

            except Exception as e:
                print(f'  NACE {int(nace)}: FAILED — {e}')

    # ============================================================ #
    #  Spec H — Pooled two-step efficient GMM (Conlon gmmnotes 141-148)
    #
    #  Runs ACF once on the full panel with year × nace2 FE, so the
    #  cluster count N_c ≈ 1,500 is adequate for the de-meaned S_hat
    #  inverse. Per-NACE two-step (prior Spec G) was unstable at thin
    #  cluster counts — see conlon_gmm_bootstrap_eval.md memory note.
    # ============================================================ #
    print(f'\n{"="*60}')
    print('  Spec H: Pooled TL two-step efficient (Conlon gmmnotes 141-148)')
    print(f'{"="*60}')

    form_pool = Formulation(spec='tl', overidentify=True,
                            pp_in_markov=True, year_fe=True,
                            nace2_fe=True)
    ext_pool = CWDLExtensions(
        survival_correction=True,
        markov_interactions=True,
        weighting='two_step',
    )

    try:
        est_pool = ACFEstimator(
            data=df, formulation=form_pool,
            optimization=Optimization(method='nm+bfgs'),
            extensions=ext_pool,
        )
        res_pool = est_pool.solve()

        row_h = {
            'spec': 'H',
            'label': 'Pooled TL two-step efficient',
            'short': 'TL 2-step (pool)',
            'nace2': 0,  # sentinel for pooled (not a NACE code)
            'N': res_pool.n_obs,
        }
        for name, coef, se in zip(res_pool.beta_names,
                                    res_pool.betas, res_pool.se):
            row_h[f'b_{name}'] = coef
            row_h[f'se_{name}'] = se

        md_h = res_pool.data
        row_h['markup_mean'] = md_h['markup'].mean()
        row_h['markup_sd'] = md_h['markup'].std()
        row_h['markup_p10'] = md_h['markup'].quantile(0.1)
        row_h['markup_p50'] = md_h['markup'].quantile(0.5)
        row_h['markup_p90'] = md_h['markup'].quantile(0.9)
        row_h['criterion'] = res_pool.gmm_criterion
        row_h['r2_first'] = res_pool.first_stage_r2
        row_h['hansen_j'] = (res_pool.hansen_j if res_pool.hansen_j is not None
                              else np.nan)
        row_h['hansen_j_p'] = (res_pool.hansen_j_pvalue
                                if res_pool.hansen_j_pvalue is not None
                                else np.nan)
        row_h['n_overid'] = res_pool.n_overid
        if res_pool.markov_coefs is not None:
            for mn, mc in zip(res_pool.markov_names, res_pool.markov_coefs):
                row_h[f'markov_{mn}'] = mc
        all_results.append(row_h)

        keep_cols = ['id', 'year', 'markup']
        if 'omega' in md_h.columns:
            keep_cols.append('omega')
        if 'alphahat' in md_h.columns:
            keep_cols.append('alphahat')
        mu_h = md_h[keep_cols].copy()
        mu_h['spec'] = 'H'
        mu_h['nace2'] = 0  # sentinel — this spec is pooled
        all_markups.append(mu_h)

        print(f'\n  Pooled N = {res_pool.n_obs}, clusters = {res_pool.n_clusters}')
        print(f'  betas (step1 → step2):')
        b1_vec = res_pool.betas_step1
        for i, name in enumerate(res_pool.beta_names):
            b1 = b1_vec[i] if b1_vec is not None else np.nan
            b2 = res_pool.betas[i]
            s2 = res_pool.se[i]
            print(f'    {name:20s}  step1={b1:+.4f}  step2={b2:+.4f} '
                  f'({s2:.4f})  Δ={b2-b1:+.4f}')
        print(f'  delta_step (max |Δβ|) = {res_pool.delta_step:.6f}')
        if res_pool.hansen_j is not None:
            print(f'  [Hansen J] χ²({res_pool.n_overid}) = '
                  f'{res_pool.hansen_j:.3f}, '
                  f'p = {res_pool.hansen_j_pvalue:.4f}')
        print(f'  Pooled markup mean = {md_h["markup"].mean():.4f}, '
              f'sd = {md_h["markup"].std():.4f}')

        # Accept/reject verdict vs plan's 0.05 stop-and-report threshold.
        if res_pool.delta_step is not None and res_pool.delta_step < 0.05:
            print(f'  ✓ delta_step {res_pool.delta_step:.4f} < 0.05 — '
                  f'pooled two-step is STABLE and ready for paper integration')
        elif res_pool.delta_step is not None:
            print(f'  ⚠ delta_step {res_pool.delta_step:.4f} ≥ 0.05 — '
                  f'pooled two-step still unstable; investigate before use')
    except Exception as e:
        print(f'  Spec H pooled two-step FAILED: {e}')

    # ============================================================ #
    #  Spec J — Test-justified hybrid (NACE 42 = CD, NACE 41/43 = TL)
    #
    #  Motivated by a nested-model Wald test on the translog
    #  cross-terms (β_k² = β_cogs² = β_kcogs = 0 under H0: CD is adequate):
    #    NACE 41: Wald χ²(3) = 28.4, p < 0.0001 → TL preferred
    #    NACE 42: Wald χ²(3) =  4.9, p = 0.18  → CD adequate (fail to reject)
    #    NACE 43: Wald χ²(3) = 44.7, p < 0.0001 → TL preferred
    #
    #  Hansen J (efficient W) agrees: NACE 42 passes overid at p = 0.86,
    #  the other two reject at ~2% and ~1%. Spec J uses CD only where
    #  both tests say CD is statistically adequate. The 473-obs / 107-firm
    #  NACE 42 sample is over-parameterized for a 6-coefficient translog.
    # ============================================================ #
    print(f'\n{"="*60}')
    print('  Spec J: Test-justified hybrid (NACE 42=CD, NACE 41/43=TL)')
    print(f'{"="*60}')

    j_form_tl = Formulation(spec='tl', overidentify=True, pp_in_markov=True)
    j_form_cd = Formulation(spec='cd', pp_in_markov=True)
    j_ext = CWDLExtensions(survival_correction=True, markov_interactions=True)

    for nace in naces:
        df_n = df[df['nace2'] == nace].copy()
        if int(nace) == 42:
            form_j = j_form_cd
            branch_tag = 'CD'
        else:
            form_j = j_form_tl
            branch_tag = 'TL'
        try:
            est = ACFEstimator(
                data=df_n,
                formulation=form_j,
                optimization=Optimization(method='nm+bfgs'),
                extensions=j_ext,
            )
            res = est.solve()

            row = {
                'spec': 'J',
                'label': 'Hybrid (NACE 42=CD, NACE 41/43=TL)',
                'short': 'TL/CD hybrid',
                'branch': branch_tag,
                'nace2': int(nace),
                'N': res.n_obs,
            }
            for name, coef, se in zip(res.beta_names, res.betas, res.se):
                row[f'b_{name}'] = coef
                row[f'se_{name}'] = se

            md_j = res.data
            row['markup_mean'] = md_j['markup'].mean()
            row['markup_sd'] = md_j['markup'].std()
            row['markup_p10'] = md_j['markup'].quantile(0.1)
            row['markup_p50'] = md_j['markup'].quantile(0.5)
            row['markup_p90'] = md_j['markup'].quantile(0.9)
            row['criterion'] = res.gmm_criterion
            row['r2_first'] = res.first_stage_r2
            row['hansen_j'] = (res.hansen_j if res.hansen_j is not None
                                else np.nan)
            row['hansen_j_p'] = (res.hansen_j_pvalue
                                  if res.hansen_j_pvalue is not None
                                  else np.nan)
            row['n_overid'] = res.n_overid
            if res.markov_coefs is not None:
                for mn, mc in zip(res.markov_names, res.markov_coefs):
                    row[f'markov_{mn}'] = mc
            all_results.append(row)

            keep_cols = ['id', 'year', 'markup']
            if 'omega' in md_j.columns:
                keep_cols.append('omega')
            if 'alphahat' in md_j.columns:
                keep_cols.append('alphahat')
            mu_j = md_j[keep_cols].copy()
            mu_j['spec'] = 'J'
            mu_j['nace2'] = int(nace)
            all_markups.append(mu_j)

            b_k = row.get('b_k', np.nan)
            b_c = row.get('b_cogs', np.nan)
            se_k = row.get('se_k', np.nan)
            se_c = row.get('se_cogs', np.nan)
            rts = b_k + b_c if np.isfinite(b_k) and np.isfinite(b_c) else np.nan
            print(f'  NACE {int(nace)} [{branch_tag}]: β_k={b_k:.4f}({se_k:.4f}), '
                  f'β_cogs={b_c:.4f}({se_c:.4f}), '
                  f'RTS={rts:.3f}, μ̄={row["markup_mean"]:.3f}, '
                  f'N={res.n_obs}')
            if res.n_overid > 0 and res.hansen_j is not None:
                print(f'    [Hansen J] χ²({res.n_overid}) = {res.hansen_j:.3f}, '
                      f'p = {res.hansen_j_pvalue:.4f}')
        except Exception as e:
            print(f'  NACE {int(nace)} [{branch_tag}]: FAILED — {e}')

    # ============================================================ #
    #  OLS baseline
    # ============================================================ #
    print(f'\n{"="*60}')
    print('  OLS (first stage only)')
    print(f'{"="*60}')
    for nace in naces:
        df_n = df[df['nace2'] == nace].dropna(subset=['go', 'k', 'cogs'])
        X = np.column_stack([np.ones(len(df_n)),
                             df_n['k'].values, df_n['cogs'].values])
        y = df_n['go'].values
        bhat = np.linalg.lstsq(X, y, rcond=None)[0]
        resid = y - X @ bhat
        r2 = 1 - np.sum(resid**2) / np.sum((y - y.mean())**2)
        alpha = np.exp(df_n['cogs'].values) / np.exp(df_n['go'].values)
        mu_ols = bhat[2] / alpha
        row = {
            'spec': 'OLS', 'label': 'OLS', 'short': 'OLS',
            'nace2': int(nace), 'N': len(df_n),
            'b_k': bhat[1], 'b_cogs': bhat[2],
            'markup_mean': mu_ols.mean(), 'markup_sd': mu_ols.std(),
            'markup_p10': np.percentile(mu_ols, 10),
            'markup_p50': np.percentile(mu_ols, 50),
            'markup_p90': np.percentile(mu_ols, 90),
            'r2_first': r2,
        }
        all_results.append(row)

        # OLS markups for premium
        mu_df = pd.DataFrame({
            'id': df_n['id'].values, 'year': df_n['year'].values,
            'markup': mu_ols, 'spec': 'OLS', 'nace2': int(nace)
        })
        all_markups.append(mu_df)

        print(f'  NACE {int(nace)}: β_k={bhat[1]:.4f}, β_cogs={bhat[2]:.4f}, '
              f'RTS={bhat[1]+bhat[2]:.3f}, μ̄={mu_ols.mean():.3f}, '
              f'R²={r2:.3f}, N={len(df_n)}')

    # ============================================================ #
    #  Combine markups and compute premiums
    # ============================================================ #
    mu_all = pd.concat(all_markups, ignore_index=True)
    res_df = pd.DataFrame(all_results)

    print(f'\n{"="*60}')
    print('  PROCUREMENT PREMIUM (log markup difference)')
    print(f'{"="*60}')
    print(f'\n  {"Spec":<20} {"Raw Δ":>8} {"SE":>8}  '
          f'{"Reg Δ":>8} {"SE":>8} {"R²":>6} {"N":>6}')
    print(f'  {"-"*75}')

    premium_rows = []
    # Include Spec H (pooled two-step) and Spec J (test-justified hybrid)
    # in the premium comparison even though they aren't in the SPECS dict
    # — they are standalone blocks above.
    premium_spec_keys = list(SPECS.keys()) + ['H', 'J', 'OLS']
    for spec_key in premium_spec_keys:
        mu_spec = mu_all[mu_all['spec'] == spec_key]
        raw, raw_se, n0, n1 = compute_premium(mu_spec, df, 'pp_dummy')
        reg, reg_se, r2, N = compute_regression_premium(mu_spec, df, 'pp_dummy')

        # Also compute with pp_ever_3y
        raw3, raw3_se, _, _ = compute_premium(mu_spec, df, 'pp_ever_3y')
        reg3, reg3_se, r2_3, _ = compute_regression_premium(mu_spec, df, 'pp_ever_3y')

        if spec_key in SPECS:
            label = SPECS[spec_key]['short']
        elif spec_key == 'H':
            label = 'TL 2-step (pool)'
        elif spec_key == 'J':
            label = 'TL/CD hybrid'
        else:
            label = 'OLS'
        print(f'  {label:<20} {raw:>8.4f} {raw_se:>8.4f}  '
              f'{reg:>8.4f} {reg_se:>8.4f} {r2:>6.3f} {N:>6}')

        premium_rows.append({
            'spec': spec_key, 'label': label,
            'raw_premium': raw, 'raw_se': raw_se,
            'reg_premium': reg, 'reg_se': reg_se, 'reg_r2': r2, 'N': N,
            'raw_premium_3y': raw3, 'raw_se_3y': raw3_se,
            'reg_premium_3y': reg3, 'reg_se_3y': reg3_se,
        })

    prem_df = pd.DataFrame(premium_rows)

    # ============================================================ #
    #  Print publication tables
    # ============================================================ #

    # Table 1: PF Coefficients
    print(f'\n\n{"="*60}')
    print('  TABLE 1: Production Function Estimates')
    print(f'{"="*60}')
    cd_specs = ['A', 'B', 'C', 'D', 'OLS']
    for nace in naces:
        print(f'\n  NACE {int(nace)}:')
        print(f'  {"Spec":<20} {"β_k":>8} {"(SE)":>8} {"β_cogs":>8} '
              f'{"(SE)":>8} {"RTS":>6} {"μ̄":>6} {"N":>6}')
        print(f'  {"-"*75}')
        sub = res_df[(res_df['nace2'] == int(nace)) &
                     (res_df['spec'].isin(cd_specs))]
        for _, r in sub.iterrows():
            bk = r.get('b_k', np.nan)
            bc = r.get('b_cogs', np.nan)
            sek = r.get('se_k', np.nan)
            sec = r.get('se_cogs', np.nan)
            rts = bk + bc if np.isfinite(bk) and np.isfinite(bc) else np.nan
            se_k_str = f'({sek:.4f})' if np.isfinite(sek) else ''
            se_c_str = f'({sec:.4f})' if np.isfinite(sec) else ''
            print(f'  {r["short"]:<20} {bk:>8.4f} {se_k_str:>8} '
                  f'{bc:>8.4f} {se_c_str:>8} '
                  f'{rts:>6.3f} {r["markup_mean"]:>6.3f} {r["N"]:>6}')

    # Markov transition parameters
    markov_cols = [c for c in res_df.columns if c.startswith('markov_')]
    if markov_cols:
        print(f'\n  Productivity Process: omega_t = g(omega_{{t-1}}, controls) + xi_t')
        print(f'  {"Spec":<12} {"NACE":>5}', end='')
        for mc in markov_cols:
            print(f' {mc.replace("markov_",""):>12}', end='')
        print()
        for _, r in res_df[res_df['spec'] == 'A'].iterrows():
            print(f'  {"A":<12} {int(r["nace2"]):>5}', end='')
            for mc in markov_cols:
                v = r.get(mc, np.nan)
                print(f' {v:>12.4f}' if np.isfinite(v) else f' {"":>12}', end='')
            print()

    # Translog
    print(f'\n  Translog:')
    tl_sub = res_df[res_df['spec'] == 'E']
    for _, r in tl_sub.iterrows():
        print(f'  NACE {r["nace2"]}: μ̄={r["markup_mean"]:.3f}, '
              f'p10={r["markup_p10"]:.3f}, p90={r["markup_p90"]:.3f}, N={r["N"]}')

    # Table 2: Premium
    print(f'\n\n{"="*60}')
    print('  TABLE 2: Procurement Premium')
    print(f'{"="*60}')
    print(f'\n  {"Spec":<20} {"pp_dummy":>12} {"(SE)":>8}  '
          f'{"pp_ever_3y":>12} {"(SE)":>8}')
    print(f'  {"-"*65}')
    print(f'  {"":20} {"--- Raw log difference ---":>42}')
    for _, r in prem_df.iterrows():
        print(f'  {r["label"]:<20} {r["raw_premium"]:>12.4f} '
              f'({r["raw_se"]:.4f})  '
              f'{r["raw_premium_3y"]:>12.4f} ({r["raw_se_3y"]:.4f})')
    print(f'\n  {"":20} {"--- Regression (+ k, cogs, yr×nace2 FE) ---":>50}')
    for _, r in prem_df.iterrows():
        print(f'  {r["label"]:<20} {r["reg_premium"]:>12.4f} '
              f'({r["reg_se"]:.4f})  '
              f'{r["reg_premium_3y"]:>12.4f} ({r["reg_se_3y"]:.4f})')

    # ============================================================ #
    #  Save outputs
    # ============================================================ #
    res_df.to_csv(f'{OUT_DIR}/paper_pf_estimates.csv', index=False)
    prem_df.to_csv(f'{OUT_DIR}/paper_premiums.csv', index=False)
    mu_all.to_csv(f'{OUT_DIR}/paper_markups.csv', index=False)
    print(f'\n  Saved: paper_pf_estimates.csv, paper_premiums.csv, paper_markups.csv')

    # Save DTA for downstream scripts (paper_tables.do, dlw_treatment_eval.py)
    import os
    for subdir in ['data', 'temp', 'tables', 'figures']:
        os.makedirs(f'{OUT_DIR}/{subdir}', exist_ok=True)

    # Pivot markups to wide format: markup_A, markup_B, ..., omega_A, alphahat
    mu_wide = mu_all.pivot_table(
        index=['id', 'year', 'nace2'], columns='spec',
        values='markup', aggfunc='first'
    ).reset_index()
    mu_wide.columns = [f'markup_{c}' if c not in ('id', 'year', 'nace2') else c
                        for c in mu_wide.columns]
    # Add omega and alphahat from spec A (baseline)
    if 'omega' in mu_all.columns:
        omega_a = mu_all[mu_all['spec'] == 'A'][['id', 'year', 'omega']].copy()
        omega_a = omega_a.rename(columns={'omega': 'omega_A'})
        mu_wide = mu_wide.merge(omega_a, on=['id', 'year'], how='left')
    if 'alphahat' in mu_all.columns:
        alpha_a = mu_all[mu_all['spec'] == 'A'][['id', 'year', 'alphahat']].copy()
        mu_wide = mu_wide.merge(alpha_a, on=['id', 'year'], how='left')
    # Add pp_dummy and controls from original data (needed by paper_tables.do)
    merge_cols = ['id', 'year', 'pp_dummy', 'k', 'cogs', 'go']
    merge_cols = [c for c in merge_cols if c in df.columns]
    mu_wide = mu_wide.merge(df[merge_cols].drop_duplicates(['id', 'year']),
                            on=['id', 'year'], how='left')

    # Save markups (both locations)
    mu_wide.to_stata(f'{OUT_DIR}/paper_markups.dta', write_index=False)
    mu_wide.to_stata(f'{OUT_DIR}/data/paper_markups.dta', write_index=False)
    for nace in mu_all['nace2'].unique():
        mu_nace = mu_wide[mu_wide['nace2'] == int(nace)]
        mu_nace.to_stata(f'{OUT_DIR}/temp/paper_markups_{int(nace)}.dta', write_index=False)

    # Pivot coefficients to wide format: one row per nace2, columns b_k_A, se_k_A, ...
    # Include Hansen J columns so paper_tables.do can merge them by NACE × spec.
    coef_vars = [c for c in res_df.columns
                 if c.startswith('b_') or c.startswith('se_')
                 or c in ('hansen_j', 'hansen_j_p', 'n_overid')]
    coef_wide = res_df.pivot_table(index='nace2', columns='spec',
                                    values=coef_vars, aggfunc='first')
    coef_wide.columns = [f'{var}_{spec}' for var, spec in coef_wide.columns]
    coef_wide = coef_wide.reset_index()
    # Add N_obs from spec A
    n_obs = res_df[res_df['spec'] == 'A'][['nace2', 'N']].rename(columns={'N': 'N_obs'})
    coef_wide = coef_wide.merge(n_obs, on='nace2', how='left')
    coef_wide.to_stata(f'{OUT_DIR}/paper_coefficients.dta', write_index=False)
    coef_wide.to_stata(f'{OUT_DIR}/data/paper_coefficients.dta', write_index=False)

    print(f'  Saved: paper_markups.dta (wide, {len(mu_wide)} obs), '
          f'paper_coefficients.dta (wide, {len(coef_wide)} rows), per-nace DTAs')


if __name__ == '__main__':
    main()
