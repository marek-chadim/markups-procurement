import numpy as np
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

def compute_r2(y, X):
    """R² from regression of y on X."""
    # Add constant
    Xc = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(Xc, y, rcond=None)
    resid = y - Xc @ beta
    sst = np.sum((y - y.mean())**2)
    if sst < 1e-10:
        return np.nan
    return 1 - np.sum(resid**2) / sst

def main():
    df = pd.read_stata(str(INPUT_DIR / 'data.dta'))
    df = df.sort_values(['id', 'year']).copy()
    df['year'] = df['year'].astype(int)

    # Compute lags
    for v in ['k', 'cogs', 'pp_dummy']:
        df[f'L_{v}'] = df.groupby('id')[v].shift(1)
    df['L2_cogs'] = df.groupby('id')['cogs'].shift(2)
    df['L_k_L_cogs'] = df['L_k'] * df['L_cogs']
    df['L_k_sq'] = df['L_k'] ** 2
    df['L_cogs_sq'] = df['L_cogs'] ** 2

    # Drop missing
    instruments = ['L_k', 'L_cogs', 'L_pp_dummy', 'L2_cogs', 'L_k_L_cogs', 'L_k_sq', 'L_cogs_sq']
    df = df.dropna(subset=instruments + ['year', 'nace2']).copy()

    # Build control matrices
    years = sorted(df['year'].unique())
    naces = sorted(df['nace2'].unique())
    X_year = np.column_stack([(df['year'] == y).astype(float).values for y in years[1:]])
    X_nace = np.column_stack([(df['nace2'] == n).astype(float).values for n in naces[1:]])
    X_both = np.column_stack([X_year, X_nace])

    # Compute R² for each instrument under each control set
    results = []
    for inst in instruments:
        z = df[inst].values
        r2_year = compute_r2(z, X_year)
        r2_nace = compute_r2(z, X_nace)
        r2_both = compute_r2(z, X_both)
        results.append({
            'instrument': inst, 'R2_year_FE': r2_year,
            'R2_nace2_FE': r2_nace, 'R2_year_and_nace2': r2_both
        })

    res_df = pd.DataFrame(results)
    print('\nABGRS Strong Exclusion Diagnostic: Partial R² of Instruments on Controls')
    print('=' * 70)
    print(res_df.to_string(index=False, float_format='%.4f'))
    print()
    print('Interpretation:')
    print('  R² close to 0 → strong exclusion holds (instrument ⊥ controls)')
    print('  R² close to 1 → weak exclusion (instrument explained by controls)')

    # Save CSV
    res_df.to_csv(OUTPUT_DIR / 'strong_exclusion_scores.csv', index=False)

    # Save LaTeX table
    tex = [
        r'\begin{table}[htbp]\centering',
        r'\caption{ABGRS Strong Exclusion Diagnostic: Instrument-Control Dependence}\label{tab:strong_exclusion}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lccc}',
        r'\toprule',
        r'Instrument & $R^2$ (Year FE) & $R^2$ (NACE FE) & $R^2$ (Year + NACE) \\',
        r'\midrule'
    ]
    display_names = {
        'L_k': r'$L.k$', 'L_cogs': r'$L.\text{cogs}$',
        'L_pp_dummy': r'$L.pp_{\text{dummy}}$', 'L2_cogs': r'$L^2.\text{cogs}$',
        'L_k_L_cogs': r'$L.k \cdot L.\text{cogs}$',
        'L_k_sq': r'$(L.k)^2$', 'L_cogs_sq': r'$(L.\text{cogs})^2$'
    }
    for _, r in res_df.iterrows():
        nm = display_names.get(r['instrument'], r['instrument'])
        tex.append(f"{nm} & {r['R2_year_FE']:.4f} & {r['R2_nace2_FE']:.4f} & {r['R2_year_and_nace2']:.4f} \\\\")
    tex += [
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r"\item \textit{Notes:} Each row reports the $R^2$ from regressing an ACF instrument on controls. Following Andrews et al. \cite{ABGRS2025}, instruments satisfying ``strong exclusion'' should have $R^2 \approx 0$ (mean-independent of controls). Instruments with $R^2$ close to 1 are dominated by control-variable variation and fail strong exclusion, making the structural estimate fragile under misspecification of the control function.",
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}'
    ]
    with open(OUTPUT_DIR / 'tables' / 'strong_exclusion_scores.tex', 'w') as f:
        f.write('\n'.join(tex))
    print(f'\nSaved: tables/strong_exclusion_scores.tex')

if __name__ == '__main__':
    main()
