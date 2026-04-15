"""ABGRS (2025, QJE) strong-exclusion robustness for ACF markup estimation.

Residualize the inputs (k, cogs, pp_dummy) that drive ACF's moment conditions
with respect to year x nace2 fixed effects BEFORE running the estimator. This
enforces ``strong exclusion'' of the instruments from control variation and
tests whether the procurement markup premium is robust to control
misspecification.
"""
import numpy as np
import pandas as pd
from pathlib import Path

from acf_estimator import (
    ACFEstimator,
    Formulation,
    Optimization,
    CWDLExtensions,
    options,
)

options.verbose = False

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'


def residualize(y, X):
    """Residualize y on X and add the original mean back.

    Re-centering keeps the residualized series on roughly the same scale as
    the raw input, which is important for the ACF estimator's first-stage
    polynomial and GMM numerics.
    """
    Xc = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(Xc, y, rcond=None)
    mean_y = y.mean()
    return (y - Xc @ beta) + mean_y


def build_controls(df):
    """Build year and NACE2 dummy matrix (drop baseline categories)."""
    years = sorted(df['year'].unique())
    naces = sorted(df['nace2'].unique())
    X_year = np.column_stack(
        [(df['year'] == y).astype(float).values for y in years[1:]]
    )
    X_nace = np.column_stack(
        [(df['nace2'] == n).astype(float).values for n in naces[1:]]
    )
    return np.column_stack([X_year, X_nace])


def estimate_premium(markups_df, df_orig):
    """Regression premium: log(mu) on pp_dummy + k + cogs + year x nace2 FE.

    Firm-clustered SEs. Uses ORIGINAL (non-residualized) pp_dummy, k, cogs
    from df_orig so the premium is measured against the researcher's
    observed treatment indicator, not the residualized ghost variable.
    """
    mu = markups_df[markups_df['markup'] > 0].copy()
    mu['lmu'] = np.log(mu['markup'])
    merge_cols = ['id', 'year', 'pp_dummy', 'k', 'cogs', 'nace2']
    orig = df_orig[merge_cols].drop_duplicates(['id', 'year'])
    mu = mu.drop(columns=[c for c in ['pp_dummy', 'k', 'cogs', 'nace2']
                          if c in mu.columns])
    mu = mu.merge(orig, on=['id', 'year'], how='left')
    mu = mu.dropna(subset=['lmu', 'pp_dummy', 'k', 'cogs', 'year', 'nace2'])
    X = build_controls(mu)
    Xmat = np.column_stack([
        np.ones(len(mu)),
        mu['pp_dummy'].values,
        mu['k'].values,
        mu['cogs'].values,
        X,
    ])
    y = mu['lmu'].values
    beta, *_ = np.linalg.lstsq(Xmat, y, rcond=None)
    resid = y - Xmat @ beta
    ids = mu['id'].values
    unique_ids = np.unique(ids)
    XtX_inv = np.linalg.inv(Xmat.T @ Xmat)
    meat = np.zeros((Xmat.shape[1], Xmat.shape[1]))
    for fid in unique_ids:
        mask = ids == fid
        Xi = Xmat[mask]
        ri = resid[mask]
        meat += (Xi.T @ ri[:, None]) @ (ri[None, :] @ Xi)
    G = len(unique_ids)
    V = XtX_inv @ meat @ XtX_inv * G / (G - 1)
    return beta[1], np.sqrt(V[1, 1]), len(mu)


def run_acf(df, tag, form):
    """Run ACF estimator, wrapped in try/except for robustness."""
    print(f'\n=== {tag} ===')
    try:
        est = ACFEstimator(
            data=df,
            formulation=form,
            optimization=Optimization(method='nm+bfgs'),
            extensions=CWDLExtensions(survival_correction=True),
            n_starts=3,
        )
        result = est.solve()
        return result, None
    except Exception as exc:
        print(f'  FAILED: {type(exc).__name__}: {exc}')
        return None, exc


def main():
    df = pd.read_stata(str(INPUT_DIR / 'data.dta'))
    df['year'] = df['year'].astype(int)
    df = (
        df.sort_values(['id', 'year'])
        .reset_index(drop=True)
        .copy()
    )
    df = df.dropna(subset=['go', 'k', 'cogs', 'pp_dummy', 'year', 'nace2']).copy()

    first_stage_controls = ['mktshare'] if 'mktshare' in df.columns else []
    form_base = Formulation(
        spec='cd',
        pp_in_markov=True,
        pp_interactions=True,
        year_fe=True,
        nace2_fe=True,
        first_stage_controls=first_stage_controls,
    )

    # Baseline: raw inputs
    result_base, err_base = run_acf(df, 'Baseline ACF (no residualization)', form_base)
    if result_base is None:
        raise RuntimeError(f'Baseline ACF failed: {err_base}')
    prem_base, se_base, n_base = estimate_premium(result_base.data, df)
    print(f'  Premium: {prem_base:.4f} (SE {se_base:.4f}, N={n_base})')

    # Residualized: residualize k, cogs, pp_dummy on year x nace2 FE
    X_ctrl = build_controls(df)
    df_res = df.copy()
    for v in ['k', 'cogs', 'pp_dummy']:
        df_res[v] = residualize(df[v].values, X_ctrl)

    result_res, err_res = run_acf(
        df_res,
        'ABGRS Strong-Exclusion ACF (residualized vs year x nace2)',
        form_base,
    )

    if result_res is None:
        prem_res, se_res, n_res = np.nan, np.nan, 0
        note = f'ACF failed on residualized data: {type(err_res).__name__}: {err_res}'
        print(f'  {note}')
    else:
        prem_res, se_res, n_res = estimate_premium(result_res.data, df)
        print(f'  Premium: {prem_res:.4f} (SE {se_res:.4f}, N={n_res})')
        note = ''

    # Save results
    results = pd.DataFrame([
        {
            'Specification': 'ACF baseline (raw)',
            'Premium': prem_base,
            'SE': se_base,
            'N': n_base,
        },
        {
            'Specification': 'ACF ABGRS residualized (year x nace2)',
            'Premium': prem_res,
            'SE': se_res,
            'N': n_res,
        },
    ])
    results.to_csv(OUTPUT_DIR / 'abgrs_residualization.csv', index=False)

    # LaTeX table
    def fmt(val, fmt_str='{:.4f}'):
        if pd.isna(val):
            return '--'
        return fmt_str.format(val)

    tex = [
        r'\begin{table}[htbp]\centering',
        r'\caption{ABGRS Strong-Exclusion Robustness: Residualized ACF Estimation}'
        r'\label{tab:abgrs}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lccc}',
        r'\toprule',
        r'Specification & Premium & SE & $N$ \\',
        r'\midrule',
        f'ACF baseline (raw inputs) & {fmt(prem_base)} & ({fmt(se_base)}) '
        f'& {n_base:,} \\\\',
        f'ACF ABGRS residualized (year $\\times$ NACE) & {fmt(prem_res)} '
        f'& ({fmt(se_res)}) & {n_res:,} \\\\',
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r'\item \textit{Notes:} Following Andrews et al.\ \cite{ABGRS2025}, I '
        r'residualize the inputs (k, cogs, pp) with respect to year $\times$ '
        r'NACE fixed effects before estimation, enforcing ``strong exclusion'' '
        r'of the instruments from control-variable variation. The premium '
        r'regression (log markup on pp\_dummy + controls + year $\times$ NACE '
        r'FE, firm-clustered SEs) is applied to both baseline and '
        r'residualized markups.',
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}',
    ]
    with open(OUTPUT_DIR / 'tables' / 'abgrs_residualization.tex', 'w') as f:
        f.write('\n'.join(tex))
    print(
        '\nSaved: abgrs_residualization.csv, tables/abgrs_residualization.tex'
    )
    if not pd.isna(prem_res):
        print(f'\nDifference: {abs(prem_res - prem_base):.4f} log points')
    if note:
        print(f'\n{note}')


if __name__ == '__main__':
    main()
