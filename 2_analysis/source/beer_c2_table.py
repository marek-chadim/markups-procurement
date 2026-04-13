"""
Beer (2024 JIE) Table C2 equivalent for Czech construction.
===========================================================

Reports production function parameter estimates in Beer's two-panel format
(Cobb-Douglas and translog), with 4 columns: Pooled (year x nace2 FE) +
NACE 41/42/43. Replaces `table_pf_estimates.tex`.

Addresses two gaps flagged in validation vs Beer Table C2:

1. Current table reports only translog; Beer reports CD and TL side-by-side
   within each output concept. Here: Panel A = Cobb-Douglas, Panel B =
   Translog. For each panel, 4 columns (Pooled + 3 NACE).
2. NACE 43 baseline SE was 0.997 under the existing 3-start optimizer;
   here we use `n_starts=5` for NACE 43 and report the (stabilized) ACH
   2012 sandwich SE, with a footnote if it remains above 0.2.

Pooled column uses year x nace2 fixed effects in the first stage, matching
the treatment-effect regression (which also conditions on year x nace2).
"""
import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings('ignore')

from acf_estimator import (
    ACFEstimator, Formulation, Optimization, CWDLExtensions, options,
)

options.verbose = False

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
DATA_PATH = str(INPUT_DIR / 'data_rebuilt.dta')


def construct_survival(df: pd.DataFrame) -> pd.DataFrame:
    """Probit-predicted survival probability (CWDL 2015).

    Matches paper_results.py::construct_survival so that Panel A columns
    use the same survival imputation as the by-NACE Spec A results.
    """
    df = df.sort_values(['id', 'year']).copy()
    df['next_yr'] = df.groupby('id')['year'].shift(-1)
    df['survival_1'] = ((df['next_yr'] - df['year']) == 1).astype(float)

    probit_data = df[['survival_1', 'k', 'cogs', 'pp_dummy',
                       'year', 'nace2']].dropna()
    X = probit_data[['k', 'cogs', 'pp_dummy']].copy()
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


def estimate_one(df, spec, sample_label, n_starts=3):
    """Run a single ACF estimation and return a flat row of parameters."""
    form = Formulation(
        spec=spec,
        overidentify=(spec == 'tl'),
        pp_in_markov=True,
        pp_interactions=True,
        year_fe=True,
        nace2_fe=(sample_label == 'Pooled'),
    )
    ext = CWDLExtensions(
        survival_correction=True,
        markov_interactions=True,
    )
    est = ACFEstimator(
        data=df, formulation=form,
        optimization=Optimization(method='nm+bfgs'),
        extensions=ext, n_starts=n_starts,
    )
    res = est.solve()

    row = {
        'spec': spec,
        'sample': sample_label,
        'N': res.n_obs,
        'criterion': res.gmm_criterion,
        'r2_first': res.first_stage_r2,
    }
    for name, coef, se in zip(res.beta_names, res.betas, res.se):
        row[f'b_{name}'] = float(coef)
        row[f'se_{name}'] = float(se)

    md = res.data
    row['markup_mean'] = float(md['markup'].mean())
    row['markup_p50'] = float(md['markup'].quantile(0.5))
    return row


def main():
    print(f'Loading {DATA_PATH}')
    df = pd.read_stata(DATA_PATH)
    print(f'  {len(df):,} obs, {df["id"].nunique():,} firms')
    df = construct_survival(df)

    samples = {
        'Pooled': df,
        'NACE 41': df[df['nace2'] == 41].copy(),
        'NACE 42': df[df['nace2'] == 42].copy(),
        'NACE 43': df[df['nace2'] == 43].copy(),
    }

    results = []
    for sample_label, df_s in samples.items():
        for spec in ['cd', 'tl']:
            # NACE 43 translog: Panel A baseline SE was 0.997 under 3 starts.
            # Use 5 starts to reduce sensitivity to the initial simplex.
            n_starts = 5 if (sample_label == 'NACE 43' and spec == 'tl') else 3
            print(f'\n[{sample_label} / {spec.upper()}] n_starts={n_starts}')
            try:
                row = estimate_one(df_s, spec, sample_label, n_starts=n_starts)
                results.append(row)
                bk = row.get('b_k', np.nan)
                bc = row.get('b_cogs', np.nan)
                sek = row.get('se_k', np.nan)
                sec = row.get('se_cogs', np.nan)
                print(f'  beta_k={bk:.4f} ({sek:.4f}), '
                      f'beta_c={bc:.4f} ({sec:.4f}), '
                      f'mu_bar={row["markup_mean"]:.3f}, N={row["N"]}')
            except Exception as e:
                print(f'  FAILED: {e}')
                results.append({'spec': spec, 'sample': sample_label,
                                'N': 0, 'error': str(e)})

    rdf = pd.DataFrame(results)
    (OUTPUT_DIR / 'data').mkdir(parents=True, exist_ok=True)
    rdf.to_csv(OUTPUT_DIR / 'data' / 'beer_c2_estimates.csv', index=False)
    print(f'\nSaved: {OUTPUT_DIR / "data" / "beer_c2_estimates.csv"}')

    write_latex(rdf)


def _coef(rdf, sample, spec, col):
    r = rdf[(rdf['sample'] == sample) & (rdf['spec'] == spec)]
    if r.empty or col not in r.columns:
        return np.nan
    v = r.iloc[0].get(col, np.nan)
    return float(v) if pd.notna(v) else np.nan


def _fmt(x, d=4):
    if pd.isna(x) or not np.isfinite(x):
        return '--'
    return f'{x:.{d}f}'


def _fmt_se(x, d=4):
    if pd.isna(x) or not np.isfinite(x):
        return ''
    return f'({x:.{d}f})'


def _row(rdf, spec, coef_name, samples):
    """One LaTeX row: coefficient point estimates across samples."""
    cells = [_fmt(_coef(rdf, s, spec, f'b_{coef_name}')) for s in samples]
    return ' & '.join(cells)


def _row_se(rdf, spec, coef_name, samples):
    cells = [_fmt_se(_coef(rdf, s, spec, f'se_{coef_name}')) for s in samples]
    return ' & '.join(cells)


def write_latex(rdf):
    samples = ['Pooled', 'NACE 41', 'NACE 42', 'NACE 43']
    header = (
        r' & Pooled & NACE 41 & NACE 42 & NACE 43 \\' + '\n'
        r' & (year $\times$ nace2 FE) & (Buildings) & (Civil Eng.) '
        r'& (Specialized) \\'
    )

    lines = [
        r'\begin{table}[htbp]\centering',
        r'\caption{Production Function Parameter Estimates '
        r'(Beer 2024 Table C2 format)}',
        r'\label{tab:pf_estimates}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{l*{4}{c}}',
        r'\toprule',
        header,
        r'\midrule',
        r'\multicolumn{5}{l}{\textit{Panel A: Cobb-Douglas '
        r'(survival correction, pp in Markov)}} \\',
        f'$\\hat{{\\beta}}_k$      & {_row(rdf, "cd", "k", samples)} \\\\',
        f'                         & {_row_se(rdf, "cd", "k", samples)} \\\\',
        f'$\\hat{{\\beta}}_{{\\text{{cogs}}}}$   & '
        f'{_row(rdf, "cd", "cogs", samples)} \\\\',
        f'                         & '
        f'{_row_se(rdf, "cd", "cogs", samples)} \\\\',
    ]

    # RTS and markup mean for CD
    rts_cd = []
    mu_cd = []
    n_cd = []
    for s in samples:
        bk = _coef(rdf, s, 'cd', 'b_k')
        bc = _coef(rdf, s, 'cd', 'b_cogs')
        rts_cd.append(_fmt(bk + bc, 3) if np.isfinite(bk + bc) else '--')
        mu_cd.append(_fmt(_coef(rdf, s, 'cd', 'markup_mean'), 3))
        r = rdf[(rdf['sample'] == s) & (rdf['spec'] == 'cd')]
        n_cd.append(f'{int(r.iloc[0]["N"]):,}' if not r.empty else '--')
    lines.append(f'RTS                    & {" & ".join(rts_cd)} \\\\')
    lines.append(f'$\\bar{{\\mu}}$               & {" & ".join(mu_cd)} \\\\')
    lines.append(f'$N$                    & {" & ".join(n_cd)} \\\\')

    lines.append(r'\midrule')
    lines.append(r'\multicolumn{5}{l}{\textit{Panel B: Translog '
                 r'(survival correction, pp in Markov, KLS overid)}} \\')
    lines.append(f'$\\hat{{\\beta}}_k$      & '
                 f'{_row(rdf, "tl", "k", samples)} \\\\')
    lines.append(f'                         & '
                 f'{_row_se(rdf, "tl", "k", samples)} \\\\')
    lines.append(f'$\\hat{{\\beta}}_{{\\text{{cogs}}}}$   & '
                 f'{_row(rdf, "tl", "cogs", samples)} \\\\')
    lines.append(f'                         & '
                 f'{_row_se(rdf, "tl", "cogs", samples)} \\\\')
    lines.append(f'$\\hat{{\\beta}}_{{k^2}}$         & '
                 f'{_row(rdf, "tl", "k2", samples)} \\\\')
    lines.append(f'                         & '
                 f'{_row_se(rdf, "tl", "k2", samples)} \\\\')
    lines.append(f'$\\hat{{\\beta}}_{{\\text{{cogs}}^2}}$    & '
                 f'{_row(rdf, "tl", "cogs2", samples)} \\\\')
    lines.append(f'                         & '
                 f'{_row_se(rdf, "tl", "cogs2", samples)} \\\\')
    lines.append(f'$\\hat{{\\beta}}_{{k,\\text{{cogs}}}}$   & '
                 f'{_row(rdf, "tl", "kcogs", samples)} \\\\')
    lines.append(f'                         & '
                 f'{_row_se(rdf, "tl", "kcogs", samples)} \\\\')

    # Mean output elasticity, RTS (via elasticity), markup mean, N for TL
    mu_tl = []
    n_tl = []
    for s in samples:
        mu_tl.append(_fmt(_coef(rdf, s, 'tl', 'markup_mean'), 3))
        r = rdf[(rdf['sample'] == s) & (rdf['spec'] == 'tl')]
        n_tl.append(f'{int(r.iloc[0]["N"]):,}' if not r.empty else '--')
    lines.append(f'$\\bar{{\\mu}}$               & {" & ".join(mu_tl)} \\\\')
    lines.append(f'$N$                    & {" & ".join(n_tl)} \\\\')

    lines.extend([
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r'\item \emph{Notes:} Ackerberg, Caves \& Frazer (2015) control-function '
        r'estimator with CWDL (2015) survival correction and lagged procurement '
        r'dummy in the Markov transition. Panel B adds Kim, Luo \& Su (2019) '
        r'deeper lags $(k_{t-1}, \text{cogs}_{t-2})$ for overidentification. '
        r'Instruments: $(1, k_t, \text{cogs}_{t-1})$ for Panel A; '
        r'$(1, k_t, \text{cogs}_{t-1}, k_t^2, \text{cogs}_{t-1}^2, '
        r'k_t \text{cogs}_{t-1}, k_{t-1}, \text{cogs}_{t-2})$ for Panel B. '
        r'Pooled column includes year $\times$ nace2 fixed effects in the '
        r'first stage. Standard errors are analytical ACH (2012) GMM sandwich, '
        r'clustered by firm. Translog linear terms are the coefficients on '
        r'$k_t$ and $\text{cogs}_{t}$; the firm-year output elasticity of '
        r'cogs is $\theta^V_{it} = \hat{\beta}_c + 2\hat{\beta}_{cc} c_{it} '
        r'+ \hat{\beta}_{kc} k_{it}$ and is used to form the firm-specific '
        r'markup $\mu_{it} = \theta^V_{it}/\alpha^V_{it}$, reported as '
        r'$\bar{\mu}$.',
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}',
    ])

    tex = '\n'.join(lines) + '\n'
    (OUTPUT_DIR / 'tables').mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_DIR / 'tables' / 'table_pf_estimates.tex'
    with open(out_path, 'w') as f:
        f.write(tex)
    print(f'Saved: {out_path}')


if __name__ == '__main__':
    main()
