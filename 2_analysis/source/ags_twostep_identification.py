"""
AGS (2017) Numerical Sensitivity + Andrews (2017) Two-Step Identification-Robust CIs.

Part 1: AGS sensitivity Lambda = -(G'WG)^{-1} G'W
  Andrews, Gentzkow, and Shapiro (2017, QJE): "Measuring the Sensitivity of
  Parameter Estimates to Estimation Moments". Lambda_{k,m} gives the local
  derivative of parameter estimate hat(beta)_k with respect to a violation
  of moment m. Large |Lambda| => load-bearing moment.

Part 2: Two-step identification-robust confidence sets (simplified Andrews 2017
  ECMA "Identification-Robust Subvector Inference"). For each coefficient,
  compare the Wald CS (hat(beta) +/- 1.96 * SE) with an S-statistic CS
  (Anderson-Rubin style) obtained by profiling the GMM criterion over a grid.
  If S-stat CS is contained in Wald CS, identification is strong and the Wald
  CS is reliable; otherwise the S-stat CS is preferred.
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize
from scipy.stats import chi2

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from acf_estimator import (
    ACFEstimator, Formulation, Optimization, CWDLExtensions,
    options as acf_options,
)

INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
acf_options.verbose = False


def _reconstruct_instrument_names(form: Formulation, est: ACFEstimator) -> list:
    """Reconstruct the column labels for est._Z based on formulation.

    Mirrors the z_vars construction inside ACFEstimator._prepare_gmm.
    Falls back to generic 'moment_i' names if the reconstruction length
    disagrees with est._Z.shape[1].
    """
    ai_list = [ai for ai in form.additional_inputs
               if ai in est._est_data.columns]
    if form.spec == 'cd':
        z_vars = ['const', 'k', 'Lcogs'] + [f'L{ai}' for ai in ai_list]
    else:
        z_vars = ['const', 'k', 'Lcogs', 'k2', 'Lcogs2', 'kLcogs']
        if form.overidentify:
            z_vars.append('Lk')
            if 'L2cogs' in est._est_data.columns:
                z_vars.append('L2cogs')
    n_expected = est._Z.shape[1]
    if len(z_vars) != n_expected:
        z_vars = [f'moment_{i}' for i in range(n_expected)]
    return z_vars


def main():
    # Load construction data
    df = pd.read_stata(str(INPUT_DIR / 'data_rebuilt.dta'))
    df = df.dropna(subset=['go', 'k', 'cogs', 'pp_dummy',
                           'year', 'nace2']).copy()

    # Restrict to NACE 41 (buildings) -- largest clean subgroup
    df_nace = df[df['nace2'] == 41].copy()
    print(f'NACE 41 sample: N={len(df_nace)}, firms={df_nace["id"].nunique()}')

    # Estimate ACF spec A (translog with pp_dummy in Markov)
    form = Formulation(spec='tl', overidentify=True,
                       pp_in_markov=True, pp_interactions=True,
                       year_fe=True, nace2_fe=False)
    est = ACFEstimator(
        data=df_nace, formulation=form,
        optimization=Optimization(method='nm+bfgs'),
        extensions=CWDLExtensions(survival_correction=True),
        n_starts=3,
    )
    result = est.solve()
    beta_hat = result.betas
    print(f'\nBaseline coefficients:')
    for nm, b, s in zip(result.beta_names, beta_hat, result.se):
        print(f'  {nm:>12s}: {b:+.4f}  (SE={s:.4f})')

    # === PART 1: AGS SENSITIVITY ===
    G = est._numerical_jacobian(beta_hat)
    W = est._W
    N = est._N
    moment_names = _reconstruct_instrument_names(form, est)

    # Lambda = -(G'WG)^-1 G'W  (shape: K_params x M_moments)
    GWG = G.T @ W @ G
    GWG_inv = np.linalg.inv(GWG)
    Lambda = -GWG_inv @ G.T @ W

    print(f'\nJacobian G shape: {G.shape}')
    print(f'Lambda shape: {Lambda.shape}')
    print(f'Moment names: {moment_names}')
    print(f'Param names:  {est._beta_names}')

    # Report top sensitivities per parameter
    print('\n=== AGS Sensitivity (Andrews, Gentzkow, Shapiro 2017) ===')
    n_top = min(5, len(moment_names))
    for k_idx, name in enumerate(est._beta_names):
        top = np.argsort(-np.abs(Lambda[k_idx]))[:n_top]
        print(f'\n  {name}: top-{n_top} moment sensitivities')
        for m_idx in top:
            print(f'    {moment_names[m_idx]:>12s}: '
                  f'Lambda = {Lambda[k_idx, m_idx]:+.4f}')

    # Save Lambda matrix (moments x params)
    ags_df = pd.DataFrame(Lambda.T, index=moment_names,
                          columns=est._beta_names)
    ags_df.to_csv(OUTPUT_DIR / 'ags_twostep_lambda.csv')

    # LaTeX: focus on b_k and b_cogs
    (OUTPUT_DIR / 'tables').mkdir(parents=True, exist_ok=True)
    k_in = 'k' in est._beta_names
    c_in = 'cogs' in est._beta_names
    tex = [
        r'\begin{table}[htbp]\centering',
        r"\caption{AGS (2017) Numerical Sensitivity: "
        r"$\Lambda = -(G'WG)^{-1}G'W$}"
        r"\label{tab:ags_twostep_sensitivity}",
        r'\begin{threeparttable}',
        r'\begin{tabular}{lcc}',
        r'\toprule',
        r'Moment & $\Lambda_{b_k}$ & $\Lambda_{b_{\text{cogs}}}$ \\',
        r'\midrule',
    ]
    for m_idx, nm in enumerate(moment_names):
        nm_clean = nm.replace('_', r'\_')
        L_k = (Lambda[est._beta_names.index('k'), m_idx]
               if k_in else np.nan)
        L_c = (Lambda[est._beta_names.index('cogs'), m_idx]
               if c_in else np.nan)
        tex.append(f'{nm_clean} & {L_k:+.4f} & {L_c:+.4f} \\\\')
    tex += [
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r"\item \textit{Notes:} Entry $(m, k)$ reports "
        r"$\Lambda_{k,m} = -[(G'WG)^{-1}G'W]_{k,m}$, the local sensitivity "
        r"of the ACF parameter estimate $\hat{\beta}_k$ to a violation of "
        r"moment condition $m$ (Andrews, Gentzkow, and Shapiro 2017). "
        r"Large $|\Lambda|$ indicates load-bearing moments. Estimated on "
        r"NACE 41 (buildings), spec A: translog with Kim, Luo \& Su (2019) "
        r"overidentification lags $(L.k, L^2.\text{cogs})$ and $pp_{it-1}$ in Markov.",
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}',
    ]
    with open(OUTPUT_DIR / 'tables' / 'ags_twostep_sensitivity.tex', 'w') as f:
        f.write('\n'.join(tex))
    print('\nSaved: tables/ags_twostep_sensitivity.tex')

    # === PART 2: TWO-STEP IDENTIFICATION-ROBUST CS (simplified Andrews 2017) ===
    se_hat = result.se
    alpha = 0.05
    z_crit = 1.96
    p_moments = G.shape[0]
    chi2_crit = chi2.ppf(1 - alpha, df=p_moments)
    print(f'\nS-stat critical value: chi2_{{{p_moments}, 0.95}} '
          f'= {chi2_crit:.4f}')

    n_grid = 100
    results_2step = []
    for k_idx, name in enumerate(est._beta_names):
        if name == 'const':
            continue
        b_hat = float(beta_hat[k_idx])
        se_k = float(se_hat[k_idx])
        if not np.isfinite(se_k) or se_k <= 0:
            print(f'  {name}: SE not finite, skipping')
            continue

        wald_lo = b_hat - z_crit * se_k
        wald_hi = b_hat + z_crit * se_k

        grid = np.linspace(b_hat - 5 * se_k, b_hat + 5 * se_k, n_grid)
        criterions = np.full(n_grid, np.nan)

        def make_objective(k_idx_fixed, g_val):
            def obj(b_other):
                betas_try = beta_hat.copy()
                betas_try[k_idx_fixed] = g_val
                j = 0
                for orig_idx in range(len(beta_hat)):
                    if orig_idx == k_idx_fixed:
                        continue
                    betas_try[orig_idx] = b_other[j]
                    j += 1
                g_vec = est._moment_vector(betas_try)
                return N * g_vec @ W @ g_vec
            return obj

        start = np.delete(beta_hat, k_idx)
        for i, g_val in enumerate(grid):
            obj_fn = make_objective(k_idx, float(g_val))
            res = minimize(
                obj_fn, start, method='Nelder-Mead',
                options={'maxiter': 200, 'xatol': 1e-4, 'fatol': 1e-4},
            )
            criterions[i] = res.fun

        in_cs = criterions <= chi2_crit
        if in_cs.any():
            s_lo = float(grid[in_cs].min())
            s_hi = float(grid[in_cs].max())
        else:
            s_lo = s_hi = b_hat

        contains = (s_lo >= wald_lo) and (s_hi <= wald_hi)

        results_2step.append({
            'param': name, 'est': b_hat, 'se': se_k,
            'wald_lo': wald_lo, 'wald_hi': wald_hi,
            's_lo': s_lo, 's_hi': s_hi,
            'robust_contained_in_wald': contains,
        })
        print(f'  {name}: Wald [{wald_lo:.4f}, {wald_hi:.4f}], '
              f'S-stat [{s_lo:.4f}, {s_hi:.4f}], contained={contains}')

    ts_df = pd.DataFrame(results_2step)
    ts_df.to_csv(OUTPUT_DIR / 'ags_twostep_ci.csv', index=False)

    tex = [
        r'\begin{table}[htbp]\centering',
        r'\caption{Two-Step Identification-Robust Confidence Sets '
        r'(Andrews 2017)}\label{tab:ags_twostep_ci}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lccccc}',
        r'\toprule',
        r'Parameter & Estimate & Wald 95\% CS & S-stat 95\% CS & '
        r'Contained? \\',
        r'\midrule',
    ]
    for r in results_2step:
        mark = r'$\checkmark$' if r['robust_contained_in_wald'] else r'$\times$'
        tex.append(
            f'${r["param"]}$ & {r["est"]:.4f} & '
            f'[{r["wald_lo"]:.4f}, {r["wald_hi"]:.4f}] & '
            f'[{r["s_lo"]:.4f}, {r["s_hi"]:.4f}] & {mark} \\\\'
        )
    tex += [
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r"\item \textit{Notes:} Wald CS is the standard 95\% interval "
        r"$\hat{\beta} \pm 1.96 \cdot SE$. S-statistic CS "
        r"(Anderson-Rubin style, Andrews 2017) is the set of $\beta$ values "
        r"where the GMM criterion $N \cdot g(\beta)' W g(\beta)$ does not "
        r"exceed the $\chi^2_{p,0.95}$ critical value; other parameters are "
        r"concentrated out. Containment of the S-statistic CS within the "
        r"Wald CS indicates strong identification; violations flag weak "
        r"identification. NACE 41 estimation.",
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}',
    ]
    with open(OUTPUT_DIR / 'tables' / 'ags_twostep_ci.tex', 'w') as f:
        f.write('\n'.join(tex))
    print('\nSaved: tables/ags_twostep_ci.tex')


if __name__ == '__main__':
    main()
