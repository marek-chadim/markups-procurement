"""Borusyak & Hull (2021, 2026) randomization inference for the procurement markup premium.

Fisher permutation test within strata. Two strategies:
  A: permute pp_dummy within year
  B: permute pp_dummy within (year, nace2)

Both preserve stratum-level treated share by construction.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

SEED = 42
K = 1000  # number of permutations


def build_design(df, pp_col):
    """Build design matrix: [const, pp, k, cogs, year x nace2 dummies]."""
    yn = df['year'].astype(str) + '_' + df['nace2'].astype(str)
    unique_yn = sorted(yn.unique())[1:]  # drop one for identification
    Xfe = np.column_stack([(yn == v).astype(float).values for v in unique_yn])
    X = np.column_stack([
        np.ones(len(df)), df[pp_col].values,
        df['k'].values, df['cogs'].values, Xfe
    ])
    return X


def estimate_premium(X, y):
    """Return coefficient on pp (column 1)."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta[1]


def permute_within_strata(df, pp_col, strata_cols, rng):
    """Permute pp_col within strata defined by strata_cols."""
    key = df[strata_cols].astype(str).agg('_'.join, axis=1)
    new_pp = df[pp_col].values.copy()
    for k in key.unique():
        mask = (key == k).values
        new_pp[mask] = rng.permutation(new_pp[mask])
    return new_pp


def main():
    rng = np.random.default_rng(SEED)

    # Load data
    df = pd.read_stata(str(INPUT_DIR / 'data.dta'))
    mk = pd.read_stata(str(OUTPUT_DIR / 'data' / 'paper_markups.dta'))
    if hasattr(mk['year'].iloc[0], 'year'):
        mk['year'] = mk['year'].dt.year
    mk['year'] = mk['year'].astype(int)
    df['year'] = df['year'].astype(int)

    panel = df[['id', 'year', 'pp_dummy', 'k', 'cogs', 'nace2']].merge(
        mk[['id', 'year', 'markup_A']], on=['id', 'year']
    )
    panel = panel[(panel['markup_A'] > 0) & panel['markup_A'].notna()].copy()
    panel['log_mu'] = np.log(panel['markup_A'])
    panel = panel.dropna(subset=['log_mu', 'pp_dummy', 'k', 'cogs']).reset_index(drop=True)

    print(f'Sample: N={len(panel):,}, firms={panel["id"].nunique()}, treated={int(panel["pp_dummy"].sum()):,}')

    y = panel['log_mu'].values

    # Baseline premium
    X_obs = build_design(panel, 'pp_dummy')
    beta_obs = estimate_premium(X_obs, y)
    print(f'\nObserved premium: {beta_obs:.4f}')

    # Strategy A: permute within year
    print(f'\nRunning K={K} permutations (Strategy A: within year)...')
    betas_A = np.zeros(K)
    for k_idx in range(K):
        perm_pp = permute_within_strata(panel, 'pp_dummy', ['year'], rng)
        panel_perm = panel.copy()
        panel_perm['pp_perm'] = perm_pp
        X_perm = build_design(panel_perm, 'pp_perm')
        betas_A[k_idx] = estimate_premium(X_perm, y)
        if (k_idx + 1) % 200 == 0:
            print(f'  {k_idx+1}/{K} done')

    # Strategy B: permute within year x nace2
    print(f'\nRunning K={K} permutations (Strategy B: within year x nace2)...')
    betas_B = np.zeros(K)
    for k_idx in range(K):
        perm_pp = permute_within_strata(panel, 'pp_dummy', ['year', 'nace2'], rng)
        panel_perm = panel.copy()
        panel_perm['pp_perm'] = perm_pp
        X_perm = build_design(panel_perm, 'pp_perm')
        betas_B[k_idx] = estimate_premium(X_perm, y)
        if (k_idx + 1) % 200 == 0:
            print(f'  {k_idx+1}/{K} done')

    # RI p-values (two-sided)
    def ri_pvalue(beta_obs, betas_perm):
        n_ge = np.sum(betas_perm >= beta_obs)
        n_lt = np.sum(betas_perm < beta_obs)
        return 2 * min(n_ge, n_lt) / len(betas_perm)

    def ri_ci(beta_obs, betas_perm, alpha=0.05):
        """95% CI via Fisher-style inversion centered on observed estimate."""
        centered = betas_perm - np.mean(betas_perm)
        lo = beta_obs - np.quantile(centered, 1 - alpha / 2)
        hi = beta_obs - np.quantile(centered, alpha / 2)
        return lo, hi

    p_A = ri_pvalue(beta_obs, betas_A)
    ci_A = ri_ci(beta_obs, betas_A)
    p_B = ri_pvalue(beta_obs, betas_B)
    ci_B = ri_ci(beta_obs, betas_B)

    print(f'\n=== RI Results ===')
    print(f'Observed premium: {beta_obs:.4f}')
    print(f'Strategy A (within year):       p={p_A:.4f}, 95% CI=[{ci_A[0]:.4f}, {ci_A[1]:.4f}]')
    print(f'  Permuted mean: {np.mean(betas_A):.4f}, SD: {np.std(betas_A):.4f}')
    print(f'Strategy B (within year x NACE): p={p_B:.4f}, 95% CI=[{ci_B[0]:.4f}, {ci_B[1]:.4f}]')
    print(f'  Permuted mean: {np.mean(betas_B):.4f}, SD: {np.std(betas_B):.4f}')

    # Save results CSV
    results = pd.DataFrame([
        {'method': 'Observed (OLS baseline)', 'premium': beta_obs, 'ci_lo': np.nan, 'ci_hi': np.nan, 'p_value': np.nan, 'K': np.nan},
        {'method': 'RI within year', 'premium': beta_obs, 'ci_lo': ci_A[0], 'ci_hi': ci_A[1], 'p_value': p_A, 'K': K},
        {'method': 'RI within year x NACE', 'premium': beta_obs, 'ci_lo': ci_B[0], 'ci_hi': ci_B[1], 'p_value': p_B, 'K': K}
    ])
    results.to_csv(OUTPUT_DIR / 'borusyak_hull_ri.csv', index=False)

    # LaTeX table
    def fmt_p(p):
        return '< 0.001' if p < 0.001 else f'{p:.3f}'
    tex = [
        r'\begin{table}[htbp]\centering',
        r'\caption{Borusyak-Hull Randomization Inference: Procurement Premium}\label{tab:bh_ri}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lcccc}',
        r'\toprule',
        r'Inference Method & Premium & 95\% CI & p-value & Perms \\',
        r'\midrule',
        f'OLS baseline (firm-clustered SE) & {beta_obs:.4f} & [0.133, 0.144] & $<$ 0.001 & --- \\\\',
        f'RI within year (coarse) & {beta_obs:.4f} & [{ci_A[0]:.4f}, {ci_A[1]:.4f}] & {fmt_p(p_A)} & {K} \\\\',
        f'RI within year $\\times$ NACE (fine) & {beta_obs:.4f} & [{ci_B[0]:.4f}, {ci_B[1]:.4f}] & {fmt_p(p_B)} & {K} \\\\',
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r'\item \textit{Notes:} Randomization inference following Borusyak and Hull \cite{BorusyakHull2021}. The baseline specification regresses $\log \mu^A$ on $pp_{it}$ with controls ($k$, cogs) and year $\times$ NACE~2-digit fixed effects. Under strategy A, pp\_dummy is permuted within each year; under strategy B, within each (year, NACE) stratum. Both preserve the stratum-level procurement share by construction. CI computed by centering the null distribution on zero and shifting to the observed estimate (Fisher inversion). $K = 1000$ permutations.',
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}'
    ]
    with open(OUTPUT_DIR / 'tables' / 'borusyak_hull_ri.tex', 'w') as f:
        f.write('\n'.join(tex))

    # Histogram figure
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, betas, title in zip(axes, [betas_A, betas_B], ['Within year', 'Within year x NACE']):
        ax.hist(betas, bins=40, color='steelblue', alpha=0.7, edgecolor='white')
        ax.axvline(beta_obs, color='red', linestyle='--', linewidth=2, label=f'Observed: {beta_obs:.3f}')
        ax.axvline(0, color='grey', linestyle=':', linewidth=1)
        ax.set_xlabel(r'Permuted $\hat{\beta}$')
        ax.set_ylabel('Frequency')
        ax.set_title(f'Strategy: {title}')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'figures' / 'borusyak_hull_ri_histogram.pdf', dpi=200, bbox_inches='tight')
    print(f'\nSaved: tables/borusyak_hull_ri.tex, figures/borusyak_hull_ri_histogram.pdf')


if __name__ == '__main__':
    main()
