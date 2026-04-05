"""
CWDL (2015 AER) Style Robustness Table for Czech Construction

Follows the structure of acf_estimation_boot15.do from
Collard-Wexler & De Loecker (2015, AER) "Reallocation and Technology:
Evidence from the US Steel Industry".

CWDL Specs (adapted to our data):
  1. Base: GO + survival correction + pp_dummy in Markov + interactions
  2. No survival: GO + pp_dummy in Markov, no survival
  3. No pp_dummy: GO + survival, no pp_dummy state variable
  4. Plain: GO only (no survival, no pp_dummy in Markov)
  5. OLS: first-stage coefficients only
  6. Value added: VA output measure + survival + pp_dummy

CWDL (2015) key features in the Mata GMM:
  - OMEGA_lag_pol = (CONST, C, SURV, OMEGA_lag, SURV*OMEGA_lag, C*OMEGA_lag)
  - Z = (const, m_lag, l, k), W = I
  - Survival via probit(exit | observables) -> phat in Markov
  - Tech dummy (C) enters Markov with interaction C*omega_lag

Our analogs:
  - C = pp_dummy (procurement participation, analog to minimill tech)
  - SURV = probit-predicted survival probability
  - m = cogs (cost of goods sold, our variable input)
  - l = not separately identified (absorbed in cogs)
  - k = fixed assets
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression

# Import from our estimator
from acf_estimator import (
    ACFEstimator, Formulation, Optimization, CWDLExtensions,
    estimate_by_industry, options
)

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'


def construct_survival(df: pd.DataFrame) -> pd.DataFrame:
    """Construct survival indicator and probit-predicted survival probability.

    Following CWDL (2015):
      1. survival_1 = 1 if firm observed next year
      2. Probit: survival_1 ~ k + cogs + pp_dummy + year×nace2
      3. phat = predicted probability -> 'survival' column
    """
    df = df.sort_values(['id', 'year']).copy()

    # survival = observed in t+1
    df['next_yr'] = df.groupby('id')['year'].shift(-1)
    df['survival_1'] = ((df['next_yr'] - df['year']) == 1).astype(float)

    # Probit (logistic approximation) on observables
    # Following CWDL: probitvar = "lne l k year#techdummy techdummy#c.k ..."
    # Our analog: k, cogs, pp_dummy, year×nace2 FEs
    probit_data = df[['survival_1', 'k', 'cogs', 'pp_dummy',
                       'year', 'nace2']].dropna()

    # build features
    X_cols = ['k', 'cogs', 'pp_dummy']
    X = probit_data[X_cols].copy()
    # year x nace2 dummies
    yr_nace = (probit_data['year'].astype(str) + '_'
               + probit_data['nace2'].astype(str))
    for val in sorted(yr_nace.unique())[1:]:
        X[f'fe_{val}'] = (yr_nace == val).astype(float)

    y = probit_data['survival_1'].values

    # fit logistic regression (probit approximation)
    model = LogisticRegression(max_iter=1000, solver='lbfgs', penalty=None)
    model.fit(X.values, y)
    phat = model.predict_proba(X.values)[:, 1]

    # merge back
    df.loc[probit_data.index, 'survival'] = phat

    n_exit = (df['survival_1'] == 0).sum()
    n_surv = (df['survival_1'] == 1).sum()
    print(f'  Survival probit: {n_surv} survive, {n_exit} exit/last-obs')
    print(f'  phat range: [{phat.min():.3f}, {phat.max():.3f}], '
          f'mean = {phat.mean():.3f}')

    return df


def run_ols_by_group(df, y_var, x_vars, group_var=None):
    """OLS regression, optionally by group. Returns coefficients."""
    results = []
    groups = [None] if group_var is None else sorted(df[group_var].unique())

    for g in groups:
        subset = df if g is None else df[df[group_var] == g]
        subset = subset.dropna(subset=[y_var] + x_vars)
        if len(subset) < len(x_vars) + 1:
            continue
        X = subset[x_vars].values.astype(np.float64)
        X = np.column_stack([np.ones(len(X)), X])
        y = subset[y_var].values.astype(np.float64)
        bhat = np.linalg.lstsq(X, y, rcond=None)[0]
        names = ['const'] + x_vars
        row = {'group': g if g is not None else 'all', 'N': len(subset)}
        for name, coef in zip(names, bhat):
            row[name] = coef
        row['RTS'] = sum(bhat[1:])  # returns to scale (excl. constant)
        results.append(row)

    return pd.DataFrame(results)


def main():
    # ================================================================ #
    #  Load data
    # ================================================================ #
    data_path = str(INPUT_DIR / 'data_rebuilt.dta')
    print(f'Loading data from {data_path}')
    df = pd.read_stata(data_path)
    print(f'Data: {len(df)} obs, {df["id"].nunique()} firms, '
          f'years {df["year"].min()}-{df["year"].max()}')
    print(f'Industries: {sorted(df["nace2"].unique())}')

    # Construct survival
    print('\n--- Constructing survival ---')
    df = construct_survival(df)

    # ================================================================ #
    #  Specification table (following CWDL acf_estimation_boot15.do)
    # ================================================================ #
    specs = {
        '1_base': {
            'label': 'Base (GO + survival + pp_dummy)',
            'formulation': Formulation(spec='cd', pp_in_markov=True),
            'extensions': CWDLExtensions(
                survival_correction=True,
                markov_interactions=True,
            ),
        },
        '2_no_survival': {
            'label': 'No survival correction',
            'formulation': Formulation(spec='cd', pp_in_markov=True),
            'extensions': CWDLExtensions(
                survival_correction=False,
                markov_interactions=True,
            ),
        },
        '3_no_pp': {
            'label': 'No pp_dummy state variable',
            'formulation': Formulation(spec='cd', pp_in_markov=False,
                                       pp_interactions=False),
            'extensions': CWDLExtensions(
                survival_correction=True,
                markov_interactions=False,
            ),
        },
        '4_plain': {
            'label': 'Plain ACF (no survival, no pp_dummy)',
            'formulation': Formulation(spec='cd', pp_in_markov=False,
                                       pp_interactions=False),
            'extensions': CWDLExtensions(),
        },
        '5_tl_base': {
            'label': 'Translog base (GO + survival + pp_dummy)',
            'formulation': Formulation(spec='tl', pp_in_markov=True),
            'extensions': CWDLExtensions(
                survival_correction=True,
                markov_interactions=True,
            ),
        },
        '6_tl_plain': {
            'label': 'Translog plain',
            'formulation': Formulation(spec='tl', pp_in_markov=False,
                                       pp_interactions=False),
            'extensions': CWDLExtensions(),
        },
    }

    # ================================================================ #
    #  Run all specifications by industry
    # ================================================================ #
    all_results = []
    all_markups = []

    for spec_key, spec_cfg in specs.items():
        print('\n' + '=' * 70)
        print(f'  SPEC: {spec_cfg["label"]}')
        print('=' * 70)

        for nace in sorted(df['nace2'].unique()):
            df_n = df[df['nace2'] == nace].copy()
            print(f'\n--- NACE2 = {nace} ---')
            try:
                est = ACFEstimator(
                    data=df_n,
                    formulation=spec_cfg['formulation'],
                    optimization=Optimization(method='nm+bfgs'),
                    extensions=spec_cfg['extensions'],
                )
                res = est.solve()

                # store results
                row = {
                    'spec': spec_key,
                    'label': spec_cfg['label'],
                    'nace2': nace,
                    'N': res.n_obs,
                }
                for name, coef, se in zip(
                    res.beta_names, res.betas, res.se
                ):
                    row[f'b_{name}'] = coef
                    row[f'se_{name}'] = se
                row['markup_mean'] = res.data['markup'].mean()
                row['markup_sd'] = res.data['markup'].std()
                row['markup_p10'] = res.data['markup'].quantile(0.1)
                row['markup_p50'] = res.data['markup'].quantile(0.5)
                row['markup_p90'] = res.data['markup'].quantile(0.9)
                row['criterion'] = res.gmm_criterion
                all_results.append(row)

                # store markups
                mu = res.data[['id', 'year', 'markup']].copy()
                mu['spec'] = spec_key
                mu['nace2'] = nace
                all_markups.append(mu)

                print(res)

            except Exception as e:
                print(f'  Failed: {e}')
                import traceback
                traceback.print_exc()

    # ================================================================ #
    #  OLS baseline (CWDL Spec 6)
    # ================================================================ #
    print('\n' + '=' * 70)
    print('  OLS (no second stage)')
    print('=' * 70)
    for nace in sorted(df['nace2'].unique()):
        df_n = df[df['nace2'] == nace].dropna(subset=['go', 'k', 'cogs'])
        ols = run_ols_by_group(df_n, 'go', ['k', 'cogs'])
        print(f'\n  NACE2 = {nace} (N={len(df_n)})')
        print(f'    k = {ols.iloc[0]["k"]:.4f}, '
              f'cogs = {ols.iloc[0]["cogs"]:.4f}, '
              f'RTS = {ols.iloc[0]["RTS"]:.4f}')

    # OLS by pp_dummy (CWDL Spec 7: separate coefficients)
    print('\n  --- OLS separate by pp_dummy ---')
    for nace in sorted(df['nace2'].unique()):
        df_n = df[df['nace2'] == nace].dropna(subset=['go', 'k', 'cogs'])
        ols_sep = run_ols_by_group(df_n, 'go', ['k', 'cogs'],
                                   group_var='pp_dummy')
        print(f'\n  NACE2 = {nace}:')
        for _, row in ols_sep.iterrows():
            pp = int(row['group'])
            print(f'    pp={pp}: k={row["k"]:.4f}, cogs={row["cogs"]:.4f}, '
                  f'RTS={row["RTS"]:.4f} (N={int(row["N"])})')

    # ================================================================ #
    #  Procurement premium comparison across specs
    # ================================================================ #
    print('\n' + '=' * 70)
    print('  PROCUREMENT PREMIUM ACROSS SPECIFICATIONS')
    print('=' * 70)

    if all_markups:
        mu_all = pd.concat(all_markups, ignore_index=True)

        for spec_key in specs:
            mu_spec = mu_all[mu_all['spec'] == spec_key].copy()
            mu_spec = mu_spec[mu_spec['markup'] > 0].copy()
            if len(mu_spec) == 0:
                continue
            mu_spec['lmu'] = np.log(mu_spec['markup'])

            # merge pp_dummy
            mu_spec = mu_spec.merge(
                df[['id', 'year', 'pp_dummy']].drop_duplicates(),
                on=['id', 'year'], how='left', suffixes=('', '_r')
            )
            if 'pp_dummy_r' in mu_spec.columns:
                mu_spec['pp_dummy'] = mu_spec['pp_dummy_r']

            pp0 = mu_spec[mu_spec['pp_dummy'] == 0]['lmu']
            pp1 = mu_spec[mu_spec['pp_dummy'] == 1]['lmu']
            if len(pp0) > 0 and len(pp1) > 0:
                raw = pp1.mean() - pp0.mean()
                print(f'  {spec_key}: raw log premium = {raw:.4f} '
                      f'(pp=0: {pp0.mean():.4f}, pp=1: {pp1.mean():.4f})')

    # ================================================================ #
    #  Summary table
    # ================================================================ #
    if all_results:
        print('\n' + '=' * 70)
        print('  SUMMARY TABLE (CWDL-style)')
        print('=' * 70)
        res_df = pd.DataFrame(all_results)

        # CD specs
        cd_specs = [s for s in specs if 'tl' not in s]
        tl_specs = [s for s in specs if 'tl' in s]

        for group_name, group_specs in [('Cobb-Douglas', cd_specs),
                                         ('Translog', tl_specs)]:
            print(f'\n  --- {group_name} ---')
            sub = res_df[res_df['spec'].isin(group_specs)]
            if len(sub) == 0:
                continue

            # print by nace2
            for nace in sorted(sub['nace2'].unique()):
                print(f'\n  NACE2 = {nace}:')
                nace_sub = sub[sub['nace2'] == nace]
                for _, row in nace_sub.iterrows():
                    spec_label = row['label']
                    k = row.get('b_k', np.nan)
                    cogs = row.get('b_cogs', np.nan)
                    se_k = row.get('se_k', np.nan)
                    se_cogs = row.get('se_cogs', np.nan)
                    mu_mean = row['markup_mean']
                    mu_sd = row['markup_sd']
                    rts = k + cogs if np.isfinite(k) and np.isfinite(cogs) else np.nan
                    print(f'    {spec_label}')
                    print(f'      k={k:.4f} ({se_k:.4f}), '
                          f'cogs={cogs:.4f} ({se_cogs:.4f}), '
                          f'RTS={rts:.3f}')
                    print(f'      mu: mean={mu_mean:.3f}, sd={mu_sd:.3f}')

        # save
        out_path = str(OUTPUT_DIR / 'cwdl_robustness_results.csv')
        res_df.to_csv(out_path, index=False)
        print(f'\n  Results saved to {out_path}')


if __name__ == '__main__':
    main()
