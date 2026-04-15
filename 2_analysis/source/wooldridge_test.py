"""wooldridge_test.py — Joint-GMM identification test of the ACF first stage.

Implements Wooldridge (2009, Economics Letters 104: 112-114) as a direct
identification test for the paper's ACF first stage. The key difference
from a two-step ACF estimator is that Wooldridge stacks the two moment
equations — the first-stage proxy and the second-stage Markov — into a
single joint GMM system with *shared* parameters across equations.

    Eq (3.3):  y_it = alpha_0 + w_it * beta + x_it * gamma + c_it * lambda + e_it
    Eq (3.4):  y_it = eta_0 + w_it * beta + x_it * gamma + f[g(x_{i,t-1}, m_{i,t-1})] + u_it

Under the random-walk case (G=1, rho=1) that Wooldridge highlights as the
linear closed-form, Eq (3.4) reduces to

    y_it = eta_0 + w_it * beta + x_it * gamma + c_{i,t-1} * lambda + u_it

so the two residuals are linear in (alpha_0, eta_0, beta, gamma, lambda).
Joint GMM is solved in closed form: theta_hat = (X'Z W Z'X)^(-1) X'Z W Z'y
with W updated to the efficient firm-clustered weight matrix at the
step-1 solution.

Wooldridge's identification argument (pp. 113-114): even if Eq (3.3)
alone fails to identify beta (the ACF critique of LP), the joint system
(3.3)+(3.4) identifies beta via the lagged instruments in Eq (3.4).
The overidentification test (Hansen J on the full stacked instrument
set) then gives a clean, interpretable diagnostic of whether the
proxy-variable assumption holds.

This script answers a specific question about the paper: is the ACF
first stage identified in the Wooldridge sense? If Wooldridge's beta_cogs
agrees with ACF's (Spec D and Spec E from paper_results.py) AND Hansen J
does not reject, the first stage is identified. If Wooldridge's beta
differs substantially, the ACF two-step is masking an identification
problem that joint GMM surfaces.

Corrections vs. prodest.ado's WRDG implementation:
  1. Shared lambda across equations. prodest.ado lines 316-321 use
     separate Stata {xd: interactionvars} and {xc: lagInteractionvars}
     parameter blocks — different lambda vectors for cit and ci,t-1.
     Wooldridge Eq (3.10) uses the SAME lambda (because g(.,.) is the
     same function, evaluated at different dates). We impose this.
  2. Cluster-robust Hansen J via the firm-clustered sandwich meat,
     not the iid form that Stata's gmm + estat overid computes.
  3. Per-NACE estimation to match the paper's main pipeline.

Outputs
-------
    output/data/wooldridge_identification.csv
    output/tables/wooldridge_identification.tex

References
----------
Wooldridge, J. M. (2009). "On estimating firm-level production functions
    using proxy variables to control for unobservables", Economics
    Letters 104(3): 112-114.
Ackerberg, Caves, and Frazer (2015). "Identification Properties of
    Recent Production Function Estimators", Econometrica 83(6).
Rovigatti, G. (2017). prodest: Stata package for production function
    estimation, version 1.0.4.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------- #
#  Paths and library imports
# --------------------------------------------------------------------- #
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
DATA_DIR = OUTPUT_DIR / 'data'
TAB_DIR = OUTPUT_DIR / 'tables'
for d in (DATA_DIR, TAB_DIR):
    d.mkdir(parents=True, exist_ok=True)

# Reusable Wooldridge joint-GMM engine lives in lib/ so paper_results.py
# can import it directly; this script is only the stand-alone test driver.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from wooldridge import wooldridge_gmm, prepare_lagged_panel  # noqa: E402


# --------------------------------------------------------------------- #
#  Data preparation — delegated to lib/wooldridge.prepare_lagged_panel
# --------------------------------------------------------------------- #

def load_panel() -> pd.DataFrame:
    """Load the rebuilt firm-year panel with lagged variables.

    Delegates the polynomial+lag construction to the library helper so
    paper_results.py and this test script share the exact same
    preparation path.
    """
    df = pd.read_stata(str(INPUT_DIR / 'data_rebuilt.dta'))
    df['year'] = df['year'].astype(int)
    return prepare_lagged_panel(df)


# --------------------------------------------------------------------- #
#  Wooldridge joint GMM (random-walk case, linear closed-form)
# --------------------------------------------------------------------- #

def wooldridge_gmm(
    df_n: pd.DataFrame,
    pf_order: int = 1,
    proxy_order: int = 2,
    use_controls: bool = False,
) -> dict:
    """Wooldridge (2009) joint-GMM identification test, linear closed form.

    Parameters
    ----------
    df_n : pd.DataFrame
        Single-NACE slice of the rebuilt panel, already sorted by firm-year.
        Must include ``go, k, cogs`` and their ``_lag`` counterparts, plus
        ``id`` for clustering. If ``use_controls=True`` also needs
        ``pp_dummy`` and ``pp_dummy_lag``.
    pf_order : int
        Production-function order. ``1`` = Cobb-Douglas (linear in k, cogs).
        ``2`` = translog (adds k^2, cogs^2, k*cogs to the PF).
    proxy_order : int
        Polynomial order for the proxy function ``c(k, cogs)`` that
        approximates productivity. Must be strictly greater than ``pf_order``
        so the proxy terms do not collide with the PF terms. Supported:
        ``2`` (k^2, cogs^2, k*cogs) and ``3`` (adds k^3, cogs^3, k^2*cogs,
        k*cogs^2).
    use_controls : bool
        If True, include pp_dummy (current in Eq 1, lagged in Eq 2) as
        a shared-coefficient control in both equations. This brings the
        Wooldridge model closer to the paper's Spec E (CD + pp in Markov).

    Returns
    -------
    dict with keys:
        theta, theta_names, se, vcov, hansen_j, hansen_p, n_overid,
        N, n_clusters, beta_cogs, se_beta_cogs, beta_k, se_beta_k.
    """
    if proxy_order <= pf_order:
        raise ValueError(
            f'proxy_order ({proxy_order}) must be > pf_order ({pf_order}); '
            f'otherwise proxy polynomial collides with PF terms.'
        )

    # ---- Select and clean the estimation sample ---------------------- #
    needed = ['go', 'k', 'cogs', 'k_lag', 'cogs_lag', 'id', 'year']
    if pf_order >= 2:
        needed += ['k2', 'cogs2', 'kcogs', 'k2_lag', 'cogs2_lag', 'kcogs_lag']
    if proxy_order >= 2:
        needed += ['k2', 'cogs2', 'kcogs', 'k2_lag', 'cogs2_lag', 'kcogs_lag']
    if proxy_order >= 3:
        needed += [
            'k3', 'cogs3', 'k2cogs', 'kcogs2',
            'k3_lag', 'cogs3_lag', 'k2cogs_lag', 'kcogs2_lag',
        ]
    if use_controls:
        needed += ['pp_dummy', 'pp_dummy_lag']
    needed = sorted(set(needed))
    df = df_n.dropna(subset=needed).copy()
    N = len(df)
    if N < 50:
        raise ValueError(f'too few observations after dropping lags: N={N}')

    y = df['go'].values.astype(np.float64)
    cl = df['id'].values

    # ---- Assemble PF and proxy columns (shared across equations) ----- #
    # PF linear always present:
    pf_cols_cur = ['cogs', 'k']
    pf_cols_lag = ['cogs_lag', 'k_lag']
    pf_names = ['cogs', 'k']
    if pf_order >= 2:
        pf_cols_cur += ['k2', 'cogs2', 'kcogs']
        pf_cols_lag += ['k2_lag', 'cogs2_lag', 'kcogs_lag']
        pf_names += ['k2', 'cogs2', 'kcogs']

    # Proxy polynomial c(k, cogs). Choose terms that are NOT already in
    # the PF regressor list — otherwise the PF coefficient is not
    # separately identified from the proxy coefficient.
    if pf_order == 1 and proxy_order == 2:
        proxy_cols_cur = ['k2', 'cogs2', 'kcogs']
        proxy_cols_lag = ['k2_lag', 'cogs2_lag', 'kcogs_lag']
        proxy_names = ['lam_k2', 'lam_cogs2', 'lam_kcogs']
    elif pf_order == 1 and proxy_order == 3:
        proxy_cols_cur = ['k2', 'cogs2', 'kcogs',
                          'k3', 'cogs3', 'k2cogs', 'kcogs2']
        proxy_cols_lag = ['k2_lag', 'cogs2_lag', 'kcogs_lag',
                          'k3_lag', 'cogs3_lag', 'k2cogs_lag', 'kcogs2_lag']
        proxy_names = ['lam_k2', 'lam_cogs2', 'lam_kcogs',
                       'lam_k3', 'lam_cogs3', 'lam_k2cogs', 'lam_kcogs2']
    elif pf_order == 2 and proxy_order == 3:
        # Translog PF already contains the quadratic cross-terms, so the
        # proxy polynomial adds only the cubic terms (k^3, cogs^3,
        # k^2*cogs, k*cogs^2). Hybrid is exactly the paper's translog.
        proxy_cols_cur = ['k3', 'cogs3', 'k2cogs', 'kcogs2']
        proxy_cols_lag = ['k3_lag', 'cogs3_lag', 'k2cogs_lag', 'kcogs2_lag']
        proxy_names = ['lam_k3', 'lam_cogs3', 'lam_k2cogs', 'lam_kcogs2']
    else:
        raise ValueError(
            f'unsupported (pf_order, proxy_order) = ({pf_order}, {proxy_order})'
        )

    # Optional controls (pp_dummy). Shared coefficient across equations.
    ctrl_cols_cur = []
    ctrl_cols_lag = []
    ctrl_names = []
    if use_controls:
        ctrl_cols_cur = ['pp_dummy']
        ctrl_cols_lag = ['pp_dummy_lag']
        ctrl_names = ['delta_pp']

    # ---- Build the 2N x k regressor matrix X_stack ------------------- #
    # Parameter ordering: [alpha_0, eta_0, (pf...), (proxy...), (ctrl...)]
    # alpha_0 is identified from Eq 1 only -> column is (1 ones ; 0s).
    # eta_0 is identified from Eq 2 only -> column is (0s ; 1 ones).
    # The shared parameters get the same column in both blocks.
    one = np.ones(N, dtype=np.float64)
    zero = np.zeros(N, dtype=np.float64)

    X1_cols = [one, zero]  # alpha_0, eta_0
    X2_cols = [zero, one]

    pf_mat_cur = df[pf_cols_cur].values.astype(np.float64)
    pf_mat_lag = df[pf_cols_lag].values.astype(np.float64)
    for j in range(pf_mat_cur.shape[1]):
        # Shared PF coefficients: current in Eq 1, current (w, x) in Eq 2.
        # Wooldridge Eq (3.4): y = ... + w_it * beta + x_it * gamma + ...
        # so the LINEAR PF terms are CURRENT in BOTH equations.
        X1_cols.append(pf_mat_cur[:, j])
        X2_cols.append(pf_mat_cur[:, j])

    # Shared proxy coefficients: c(k, cogs) in Eq 1, c(k_lag, cogs_lag) in Eq 2.
    proxy_mat_cur = df[proxy_cols_cur].values.astype(np.float64)
    proxy_mat_lag = df[proxy_cols_lag].values.astype(np.float64)
    for j in range(proxy_mat_cur.shape[1]):
        X1_cols.append(proxy_mat_cur[:, j])
        X2_cols.append(proxy_mat_lag[:, j])

    if use_controls:
        ctrl_mat_cur = df[ctrl_cols_cur].values.astype(np.float64)
        ctrl_mat_lag = df[ctrl_cols_lag].values.astype(np.float64)
        for j in range(ctrl_mat_cur.shape[1]):
            X1_cols.append(ctrl_mat_cur[:, j])
            X2_cols.append(ctrl_mat_lag[:, j])

    X1 = np.column_stack(X1_cols)  # N x k
    X2 = np.column_stack(X2_cols)  # N x k
    k = X1.shape[1]
    param_names = (
        ['alpha_0', 'eta_0']
        + [f'beta_{nm}' for nm in pf_names]
        + proxy_names
        + ctrl_names
    )
    assert len(param_names) == k, \
        f'param_names length {len(param_names)} != k={k}'

    X_stack = np.vstack([X1, X2])  # 2N x k
    y_stack = np.concatenate([y, y])  # 2N

    # ---- Instruments per Wooldridge Eq (3.5)-(3.6) ------------------- #
    # Eq 1 instruments: (1, w, x, c^o) where c^o is c without duplicate x.
    #                    In our CD case that's (1, cogs, k, k^2, cogs^2, k*cogs, [pp]).
    # Eq 2 instruments: (1, x, w_lag, c_lag, [pp_lag], q_lag).
    # We do NOT add nonlinear functions q beyond the polynomial c itself
    # (keeps the overidentification manageable and matches Wooldridge's
    # linear-closed-form prescription).
    Z1_cols = [one]  # constant
    Z1_cols.append(df['cogs'].values.astype(np.float64))  # w
    Z1_cols.append(df['k'].values.astype(np.float64))  # x
    for j in range(proxy_mat_cur.shape[1]):
        Z1_cols.append(proxy_mat_cur[:, j])
    if use_controls:
        Z1_cols.append(df['pp_dummy'].values.astype(np.float64))
    Z1 = np.column_stack(Z1_cols)  # N x k_z1

    Z2_cols = [one]
    Z2_cols.append(df['k'].values.astype(np.float64))  # contemporaneous state
    Z2_cols.append(df['cogs_lag'].values.astype(np.float64))  # lagged free
    Z2_cols.append(df['k_lag'].values.astype(np.float64))  # lagged state
    for j in range(proxy_mat_lag.shape[1]):
        Z2_cols.append(proxy_mat_lag[:, j])
    if use_controls:
        Z2_cols.append(df['pp_dummy_lag'].values.astype(np.float64))
    Z2 = np.column_stack(Z2_cols)  # N x k_z2

    k_z1 = Z1.shape[1]
    k_z2 = Z2.shape[1]
    k_z = k_z1 + k_z2
    n_overid = k_z - k
    if n_overid < 0:
        raise ValueError(
            f'underidentified: k_z={k_z} < k={k} '
            f'(pf_order={pf_order}, proxy_order={proxy_order})'
        )

    # Block-diagonal Z_stack with Z1 upper-left, Z2 lower-right.
    Z_stack = np.block([
        [Z1, np.zeros((N, k_z2))],
        [np.zeros((N, k_z1)), Z2],
    ])  # 2N x k_z

    # ---- Step 1: linear GMM with initial W = (Z'Z / 2N)^(-1) --------- #
    ZZ = Z_stack.T @ Z_stack / (2 * N)
    W0 = np.linalg.pinv(ZZ)
    XZ = X_stack.T @ Z_stack  # k x k_z
    Zy = Z_stack.T @ y_stack  # k_z
    ZX = Z_stack.T @ X_stack  # k_z x k
    # theta = (X'Z W Z'X)^(-1) X'Z W Z'y
    A1 = XZ @ W0 @ ZX
    b1 = XZ @ W0 @ Zy
    theta1 = np.linalg.solve(A1, b1)

    # ---- Step 2: efficient W from firm-clustered S at theta1 --------- #
    r1_vec = y - X1 @ theta1
    r2_vec = y - X2 @ theta1
    g_eq1 = Z1 * r1_vec[:, None]  # N x k_z1
    g_eq2 = Z2 * r2_vec[:, None]  # N x k_z2
    g_i = np.concatenate([g_eq1, g_eq2], axis=1)  # N x k_z

    uniq_c = np.unique(cl)
    N_c = int(len(uniq_c))
    S = np.zeros((k_z, k_z), dtype=np.float64)
    for c in uniq_c:
        sel = cl == c
        g_c = g_i[sel].sum(axis=0)
        S += np.outer(g_c, g_c)
    S = S / N * (N_c / (N_c - 1))

    cond_S = np.linalg.cond(S)
    if cond_S > 1e14:
        W_eff = np.linalg.pinv(S)
    else:
        try:
            W_eff = np.linalg.inv(S)
        except np.linalg.LinAlgError:
            W_eff = np.linalg.pinv(S)

    # ---- Step 3: re-solve with efficient W --------------------------- #
    A2 = XZ @ W_eff @ ZX
    b2 = XZ @ W_eff @ Zy
    theta_hat = np.linalg.solve(A2, b2)

    # ---- Sandwich variance (Wooldridge 2002, Ch 14) ----------------- #
    # With the un-normalized moment g_i = Z_i r_i and Jacobian G = Z'X,
    # the standard linear-GMM sandwich is
    #     Var(theta_hat) = (X'Z W Z'X)^(-1) (X'Z W S W Z'X) (X'Z W Z'X)^(-1)
    # with NO 1/N scaling (the N cancels because S is O(1) per-cluster
    # and X'Z grows like O(N)). Under efficient W = S^(-1) this collapses
    # to Var(theta_hat) = (X'Z W Z'X)^(-1).
    #
    # We keep the full sandwich form so the formula is valid under
    # finite-sample drift between W_eff and the true S^(-1).
    try:
        bread = np.linalg.inv(XZ @ W_eff @ ZX)
    except np.linalg.LinAlgError:
        bread = np.linalg.pinv(XZ @ W_eff @ ZX)
    meat = (XZ @ W_eff) @ S @ (W_eff @ ZX)
    V = bread @ meat @ bread
    se = np.sqrt(np.maximum(np.diag(V), 0.0))

    # ---- Hansen J at theta_hat --------------------------------------- #
    r1_final = y - X1 @ theta_hat
    r2_final = y - X2 @ theta_hat
    m1 = Z1.T @ r1_final / N
    m2 = Z2.T @ r2_final / N
    m_final = np.concatenate([m1, m2])
    j_stat = float(N * m_final @ W_eff @ m_final)
    if n_overid > 0:
        p_value = float(1.0 - chi2.cdf(j_stat, df=n_overid))
    else:
        p_value = 1.0

    # ---- Package results --------------------------------------------- #
    idx_beta = param_names.index('beta_cogs')
    idx_gamma = param_names.index('beta_k')
    out = dict(
        theta=theta_hat,
        theta_names=param_names,
        se=se,
        vcov=V,
        hansen_j=j_stat,
        hansen_p=p_value,
        n_overid=n_overid,
        N=N,
        n_clusters=N_c,
        beta_cogs=float(theta_hat[idx_beta]),
        se_beta_cogs=float(se[idx_beta]),
        beta_k=float(theta_hat[idx_gamma]),
        se_beta_k=float(se[idx_gamma]),
    )
    return out


# --------------------------------------------------------------------- #
#  Comparison with ACF results from paper_pf_estimates.csv
# --------------------------------------------------------------------- #

def load_acf_reference() -> pd.DataFrame:
    """Load the ACF estimates from paper_pf_estimates.csv for comparison."""
    path = OUTPUT_DIR / 'paper_pf_estimates.csv'
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def find_acf(ref: pd.DataFrame, spec: str, nace: int) -> dict:
    """Extract (beta_k, beta_cogs) from the ACF CSV for a given spec/NACE."""
    if ref.empty:
        return {}
    sub = ref[(ref['spec'] == spec) & (ref['nace2'] == nace)]
    if sub.empty:
        return {}
    row = sub.iloc[0]
    return dict(
        beta_k=row.get('b_k', np.nan),
        se_beta_k=row.get('se_k', np.nan),
        beta_cogs=row.get('b_cogs', np.nan),
        se_beta_cogs=row.get('se_cogs', np.nan),
        N=row.get('N', np.nan),
    )


# --------------------------------------------------------------------- #
#  Main driver
# --------------------------------------------------------------------- #

def main():
    print('=' * 72)
    print(' Wooldridge (2009 EL) joint-GMM identification test')
    print(' ACF first stage identified iff Wooldridge beta agrees with ACF beta')
    print(' AND Hansen J fails to reject overid')
    print('=' * 72)

    df = load_panel()
    print(f'\nRaw panel: N={len(df)}, firms={df["id"].nunique()}, '
          f'years {df["year"].min()}-{df["year"].max()}')

    ref = load_acf_reference()
    if not ref.empty:
        print(f'ACF reference loaded: {len(ref)} rows from '
              f'paper_pf_estimates.csv')

    naces = sorted(int(v) for v in df['nace2'].dropna().unique())

    # Three variants mapped to paper's existing specs for head-to-head
    # comparison. The Wooldridge beta should closely match the ACF beta
    # if the first-stage identification is sound.
    variants = [
        dict(name='WLD-CD-plain',
             pf_order=1, proxy_order=2, use_controls=False,
             acf_spec='D', acf_spec_label='Plain CD (Spec D)'),
        dict(name='WLD-CD-pp',
             pf_order=1, proxy_order=2, use_controls=True,
             acf_spec='E', acf_spec_label='CD base (Spec E)'),
        dict(name='WLD-TL-pp',
             pf_order=2, proxy_order=3, use_controls=True,
             acf_spec='A', acf_spec_label='Translog (Spec A)'),
    ]

    rows = []
    for v in variants:
        print()
        print('-' * 72)
        print(f'  {v["name"]}  (vs ACF {v["acf_spec_label"]})')
        print(f'  pf_order={v["pf_order"]}  proxy_order={v["proxy_order"]}  '
              f'controls={v["use_controls"]}')
        print('-' * 72)
        for nace in naces:
            df_n = df[df['nace2'] == nace]
            try:
                wld = wooldridge_gmm(
                    df_n,
                    pf_order=v['pf_order'],
                    proxy_order=v['proxy_order'],
                    use_controls=v['use_controls'],
                )
            except Exception as e:
                print(f'  NACE {nace}: FAILED — {type(e).__name__}: {e}')
                continue

            acf = find_acf(ref, v['acf_spec'], nace)
            row = dict(
                variant=v['name'],
                acf_spec=v['acf_spec'],
                nace2=nace,
                N=wld['N'],
                n_clusters=wld['n_clusters'],
                n_overid=wld['n_overid'],
                wld_beta_cogs=wld['beta_cogs'],
                wld_se_beta_cogs=wld['se_beta_cogs'],
                wld_beta_k=wld['beta_k'],
                wld_se_beta_k=wld['se_beta_k'],
                hansen_j=wld['hansen_j'],
                hansen_p=wld['hansen_p'],
            )
            if acf:
                row['acf_beta_cogs'] = acf['beta_cogs']
                row['acf_se_beta_cogs'] = acf['se_beta_cogs']
                row['acf_beta_k'] = acf['beta_k']
                row['acf_se_beta_k'] = acf['se_beta_k']
                row['acf_N'] = acf['N']
                row['delta_beta_cogs'] = wld['beta_cogs'] - acf['beta_cogs']
                row['delta_beta_k'] = wld['beta_k'] - acf['beta_k']
                # Approximate t-stat for the difference: conservative (assumes
                # independent estimators).
                pooled_se = np.sqrt(
                    wld['se_beta_cogs'] ** 2 + (acf['se_beta_cogs'] or 0) ** 2
                )
                if pooled_se > 0:
                    row['t_diff_cogs'] = row['delta_beta_cogs'] / pooled_se
                else:
                    row['t_diff_cogs'] = np.nan
            rows.append(row)

            print(
                f'  NACE {nace}: WLD β_cogs={wld["beta_cogs"]:+.4f} '
                f'({wld["se_beta_cogs"]:.4f})  '
                f'β_k={wld["beta_k"]:+.4f} ({wld["se_beta_k"]:.4f})  '
                f'Hansen J χ²({wld["n_overid"]})={wld["hansen_j"]:.2f} '
                f'p={wld["hansen_p"]:.4f}  '
                f'N={wld["N"]}'
            )
            if acf:
                d = row['delta_beta_cogs']
                t = row.get('t_diff_cogs', np.nan)
                # Economic + statistical thresholds. STRONG AGREE: small
                # economic gap AND non-significant statistical difference.
                # WEAK AGREE: small economic gap but tight SEs inflate the
                # t-stat. DISAGREE: large economic gap.
                abs_d = abs(d)
                abs_t = abs(t)
                if abs_d < 0.03:
                    verdict = '→ STRONG AGREE (Δβ<0.03)'
                elif abs_d < 0.05 and abs_t < 2.0:
                    verdict = '→ AGREE (Δβ<0.05, |t|<2)'
                elif abs_d < 0.05:
                    verdict = '→ WEAK AGREE (Δβ<0.05 but tight SEs)'
                else:
                    verdict = '→ DISAGREE (Δβ≥0.05)'
                print(
                    f'    vs ACF {v["acf_spec"]}: β_cogs={acf["beta_cogs"]:+.4f} '
                    f'({acf["se_beta_cogs"] or 0:.4f})  '
                    f'Δβ_cogs={d:+.4f}  |t_diff|={abs_t:.2f}  '
                    f'{verdict}'
                )

    # ---- Save results ------------------------------------------------ #
    out_df = pd.DataFrame(rows)
    csv_path = DATA_DIR / 'wooldridge_identification.csv'
    out_df.to_csv(csv_path, index=False)
    print(f'\n[save] {csv_path.relative_to(OUTPUT_DIR.parent)}')

    # ---- Write LaTeX table (per variant, per NACE) ------------------- #
    tex_path = TAB_DIR / 'wooldridge_identification.tex'
    with open(tex_path, 'w') as f:
        f.write(r'\begin{table}[htbp]\centering' + '\n')
        f.write(r'\caption{Wooldridge (2009) Joint-GMM Identification Test}'
                + '\n')
        f.write(r'\label{tab:wooldridge_test}' + '\n')
        f.write(r'\begin{threeparttable}' + '\n')
        f.write(r'\begin{tabular}{llcccccccc}' + '\n')
        f.write(r'\toprule' + '\n')
        f.write(
            r'Variant & NACE & $N$ & $\hat\beta_{\text{cogs}}^{\text{WLD}}$'
            r' & SE & $\hat\beta_{\text{cogs}}^{\text{ACF}}$ & SE'
            r' & $\Delta\beta$ & $J$ & $p_J$ \\' + '\n'
        )
        f.write(r'\midrule' + '\n')
        for _, r in out_df.iterrows():
            f.write(
                f'{r["variant"]} & {int(r["nace2"])} & {int(r["N"])} '
                f'& {r["wld_beta_cogs"]:+.4f} & ({r["wld_se_beta_cogs"]:.4f}) '
                f'& {r.get("acf_beta_cogs", float("nan")):+.4f} '
                f'& ({r.get("acf_se_beta_cogs", float("nan")):.4f}) '
                f'& {r.get("delta_beta_cogs", float("nan")):+.4f} '
                f'& {r["hansen_j"]:.2f} & {r["hansen_p"]:.4f} \\\\\n'
            )
        f.write(r'\bottomrule' + '\n')
        f.write(r'\end{tabular}' + '\n')
        f.write(r'\begin{tablenotes}\footnotesize' + '\n')
        f.write(
            r'\item \emph{Notes:} Wooldridge (2009 \emph{Econ.\ Lett.}) joint '
            r'GMM in the random-walk Markov case ($G{=}1,\,\rho{=}1$), linear '
            r'closed form with firm-clustered sandwich SEs and efficient '
            r'weight matrix $\hat W = \hat S^{-1}$ from the clustered moment '
            r'covariance. Parameters $\beta, \gamma, \lambda$ are shared '
            r'across the two residual equations (correcting '
            r'\texttt{prodest.ado} which uses separate parameter blocks). '
            r'Eq 1 instruments: $(1, w_{it}, x_{it}, c_{it}^{o})$; Eq 2 '
            r'instruments: $(1, x_{it}, w_{i,t-1}, c_{i,t-1})$ per '
            r'Wooldridge (3.5)--(3.6). Hansen $J$ is the textbook '
            r'efficient-weighted form $J = N\,\bar m^{\prime} '
            r'\hat W\,\bar m \sim \chi^2(\text{overid})$. \textbf{WLD-CD-plain} '
            r'vs paper Spec D; \textbf{WLD-CD-pp} vs Spec E; '
            r'\textbf{WLD-TL-pp} vs Spec A. The first stage is identified '
            r'in the Wooldridge sense if $|\Delta\beta_{\text{cogs}}|$ is '
            r'small relative to the combined SE and Hansen $J$ does not '
            r'reject overidentification.' + '\n'
        )
        f.write(r'\end{tablenotes}' + '\n')
        f.write(r'\end{threeparttable}' + '\n')
        f.write(r'\end{table}' + '\n')
    print(f'[save] {tex_path.relative_to(OUTPUT_DIR.parent)}')

    # ---- Identification verdict summary ------------------------------ #
    print()
    print('=' * 72)
    print(' Identification verdict (per variant, pooled across NACE)')
    print('=' * 72)
    for v in variants:
        sub = out_df[out_df['variant'] == v['name']]
        if sub.empty:
            print(f'  {v["name"]}: no successful runs')
            continue
        max_abs_delta = sub['delta_beta_cogs'].abs().max()
        max_abs_t = sub['t_diff_cogs'].abs().max()
        j_rejects = int((sub['hansen_p'] < 0.05).sum())
        # Economic-significance-based verdict (|Δβ_cogs| threshold).
        if max_abs_delta < 0.03:
            id_status = 'IDENTIFIED (strong)'
        elif max_abs_delta < 0.05:
            id_status = 'IDENTIFIED (weak)'
        else:
            id_status = 'PARTIAL — largest |Δβ_cogs| exceeds 0.05'
        overid_status = (f'Hansen J rejects in {j_rejects}/{len(sub)} NACEs'
                         if j_rejects > 0 else 'Hansen J passes all NACEs')
        print(
            f'  {v["name"]:<14} max |Δβ_cogs|={max_abs_delta:.4f}, '
            f'max |t_diff|={max_abs_t:.2f}'
        )
        print(f'               {id_status}; {overid_status}')


if __name__ == '__main__':
    main()
