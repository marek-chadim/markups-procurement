"""wooldridge.py — Reusable Wooldridge (2009 EL) joint-GMM engine.

Implements Wooldridge's joint-GMM estimator for the proxy-variable
production-function model. Random-walk Markov case (G=1, rho=1) with
linear closed-form estimation, firm-clustered sandwich SEs, and the
textbook efficient-weighted Hansen J overidentification test.

The key conceptual advantage over ACF two-step (as Wooldridge 2009 EL
argues): joint GMM combines the first-stage proxy and the second-stage
Markov moment conditions with *shared* parameters, letting each equation
contribute identifying information for beta. This is more efficient than
two-step (Wooldridge 2009 §2, p. 113) and avoids the ACF first-stage
identification critique by leveraging Eq 2's lagged instruments.

Corrections vs. Rovigatti's prodest.ado (WRDG option):
  1. Shared lambda across the two equations. prodest uses separate
     Stata parameter blocks ({xd:} and {xc:}) for current- and lagged-
     polynomial terms, effectively estimating DIFFERENT lambda vectors.
     Wooldridge's theory uses the same g(.,.) function at both dates,
     so lambda MUST be shared.
  2. Firm-clustered Hansen J (not iid).
  3. Per-NACE + translog support via (pf_order, proxy_order) parameters.
  4. Cwd-independent unit test: linear GMM closed form, no optimization.

References
----------
Wooldridge, J. M. (2009). On estimating firm-level production functions
    using proxy variables to control for unobservables. Economics Letters
    104: 112-114.
Ackerberg, Caves, Frazer (2015). Identification Properties of Recent
    Production Function Estimators. Econometrica 83(6).
Rovigatti, G. (2017). prodest: Stata package for production function
    estimation.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import chi2


def wooldridge_gmm(
    df_n: pd.DataFrame,
    pf_order: int = 1,
    proxy_order: int = 2,
    use_controls: bool = False,
    compute_markups: bool = False,
) -> dict:
    """Wooldridge (2009) joint-GMM identification test, linear closed form.

    Parameters
    ----------
    df_n : pd.DataFrame
        Single-NACE (or pooled) slice of the rebuilt panel. Must include
        columns ``go, k, cogs, id, year`` plus lagged counterparts
        ``k_lag, cogs_lag``. For ``pf_order>=2`` or ``proxy_order>=2``,
        also needs ``k2, cogs2, kcogs`` and their lags. For
        ``proxy_order>=3``, additionally needs ``k3, cogs3, k2cogs,
        kcogs2`` and their lags. If ``use_controls=True``, needs
        ``pp_dummy, pp_dummy_lag``.
    pf_order : int
        Production-function order. ``1`` = Cobb-Douglas. ``2`` = translog
        (PF contains k^2, cogs^2, k*cogs as coefficients to estimate).
    proxy_order : int
        Polynomial order for the proxy function c(k, cogs) that
        approximates productivity. Must be strictly greater than
        ``pf_order`` so proxy and PF columns do not collide.
    use_controls : bool
        Include pp_dummy (current in Eq 1, lagged in Eq 2) as a shared-
        coefficient control.
    compute_markups : bool
        If True, compute firm-year markups from the estimated theta
        vector and include them in the returned dict as ``markups``
        (a pd.DataFrame with ``id, year, theta, markup``).

    Returns
    -------
    dict with keys theta, theta_names, se, vcov, hansen_j, hansen_p,
    n_overid, N, n_clusters, beta_cogs, se_beta_cogs, beta_k, se_beta_k,
    and optionally markups.
    """
    if proxy_order <= pf_order:
        raise ValueError(
            f'proxy_order ({proxy_order}) must be > pf_order ({pf_order})'
        )

    # ---- Select and clean the estimation sample ---------------------- #
    needed = ['go', 'k', 'cogs', 'k_lag', 'cogs_lag', 'id', 'year']
    if pf_order >= 2 or proxy_order >= 2:
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
    alpha_data = np.exp(df['cogs'].values - df['go'].values)  # cost share

    # ---- Assemble PF and proxy columns (shared across equations) ----- #
    pf_cols_cur = ['cogs', 'k']
    pf_cols_lag = ['cogs_lag', 'k_lag']
    pf_names = ['cogs', 'k']
    if pf_order >= 2:
        pf_cols_cur += ['k2', 'cogs2', 'kcogs']
        pf_cols_lag += ['k2_lag', 'cogs2_lag', 'kcogs_lag']
        pf_names += ['k2', 'cogs2', 'kcogs']

    # Proxy polynomial. Must not duplicate PF columns.
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
        proxy_cols_cur = ['k3', 'cogs3', 'k2cogs', 'kcogs2']
        proxy_cols_lag = ['k3_lag', 'cogs3_lag', 'k2cogs_lag', 'kcogs2_lag']
        proxy_names = ['lam_k3', 'lam_cogs3', 'lam_k2cogs', 'lam_kcogs2']
    else:
        raise ValueError(
            f'unsupported (pf_order, proxy_order) = ({pf_order}, {proxy_order})'
        )

    ctrl_names: list = []
    if use_controls:
        ctrl_names = ['delta_pp']

    # ---- Build X_stack and Z_stack ----------------------------------- #
    one = np.ones(N, dtype=np.float64)
    zero = np.zeros(N, dtype=np.float64)

    X1_cols = [one, zero]
    X2_cols = [zero, one]

    pf_mat_cur = df[pf_cols_cur].values.astype(np.float64)
    pf_mat_lag = df[pf_cols_lag].values.astype(np.float64)
    for j in range(pf_mat_cur.shape[1]):
        X1_cols.append(pf_mat_cur[:, j])
        X2_cols.append(pf_mat_cur[:, j])

    proxy_mat_cur = df[proxy_cols_cur].values.astype(np.float64)
    proxy_mat_lag = df[proxy_cols_lag].values.astype(np.float64)
    for j in range(proxy_mat_cur.shape[1]):
        X1_cols.append(proxy_mat_cur[:, j])
        X2_cols.append(proxy_mat_lag[:, j])

    if use_controls:
        X1_cols.append(df['pp_dummy'].values.astype(np.float64))
        X2_cols.append(df['pp_dummy_lag'].values.astype(np.float64))

    X1 = np.column_stack(X1_cols)
    X2 = np.column_stack(X2_cols)
    k = X1.shape[1]
    param_names = (
        ['alpha_0', 'eta_0']
        + [f'beta_{nm}' for nm in pf_names]
        + proxy_names
        + ctrl_names
    )
    assert len(param_names) == k

    X_stack = np.vstack([X1, X2])
    y_stack = np.concatenate([y, y])

    # ---- Instruments per Wooldridge Eq (3.5)-(3.6) ------------------- #
    Z1_cols = [one, df['cogs'].values.astype(np.float64),
               df['k'].values.astype(np.float64)]
    for j in range(proxy_mat_cur.shape[1]):
        Z1_cols.append(proxy_mat_cur[:, j])
    if use_controls:
        Z1_cols.append(df['pp_dummy'].values.astype(np.float64))
    Z1 = np.column_stack(Z1_cols)

    Z2_cols = [one, df['k'].values.astype(np.float64),
               df['cogs_lag'].values.astype(np.float64),
               df['k_lag'].values.astype(np.float64)]
    for j in range(proxy_mat_lag.shape[1]):
        Z2_cols.append(proxy_mat_lag[:, j])
    if use_controls:
        Z2_cols.append(df['pp_dummy_lag'].values.astype(np.float64))
    Z2 = np.column_stack(Z2_cols)

    k_z1 = Z1.shape[1]
    k_z2 = Z2.shape[1]
    k_z = k_z1 + k_z2
    n_overid = k_z - k
    if n_overid < 0:
        raise ValueError(
            f'underidentified: k_z={k_z} < k={k}'
        )

    Z_stack = np.block([
        [Z1, np.zeros((N, k_z2))],
        [np.zeros((N, k_z1)), Z2],
    ])

    # ---- Step 1: initial weight (Z'Z)^(-1) ---------------------------- #
    ZZ = Z_stack.T @ Z_stack / (2 * N)
    W0 = np.linalg.pinv(ZZ)
    XZ = X_stack.T @ Z_stack
    Zy = Z_stack.T @ y_stack
    ZX = Z_stack.T @ X_stack
    A1 = XZ @ W0 @ ZX
    b1 = XZ @ W0 @ Zy
    theta1 = np.linalg.solve(A1, b1)

    # ---- Step 2: clustered Ŝ at theta1 ------------------------------- #
    r1_vec = y - X1 @ theta1
    r2_vec = y - X2 @ theta1
    g_eq1 = Z1 * r1_vec[:, None]
    g_eq2 = Z2 * r2_vec[:, None]
    g_i = np.concatenate([g_eq1, g_eq2], axis=1)
    uniq_c = np.unique(cl)
    N_c = int(len(uniq_c))
    S = np.zeros((k_z, k_z), dtype=np.float64)
    for c in uniq_c:
        sel = cl == c
        g_c = g_i[sel].sum(axis=0)
        S += np.outer(g_c, g_c)
    S = S / N * (N_c / (N_c - 1))

    if np.linalg.cond(S) > 1e14:
        W_eff = np.linalg.pinv(S)
    else:
        try:
            W_eff = np.linalg.inv(S)
        except np.linalg.LinAlgError:
            W_eff = np.linalg.pinv(S)

    # ---- Step 3: efficient re-solve ---------------------------------- #
    A2 = XZ @ W_eff @ ZX
    b2 = XZ @ W_eff @ Zy
    theta_hat = np.linalg.solve(A2, b2)

    # ---- Sandwich variance (no 1/N scaling; un-normalized moments) --- #
    try:
        bread = np.linalg.inv(XZ @ W_eff @ ZX)
    except np.linalg.LinAlgError:
        bread = np.linalg.pinv(XZ @ W_eff @ ZX)
    meat = (XZ @ W_eff) @ S @ (W_eff @ ZX)
    V = bread @ meat @ bread
    se = np.sqrt(np.maximum(np.diag(V), 0.0))

    # ---- Hansen J ---------------------------------------------------- #
    r1_final = y - X1 @ theta_hat
    r2_final = y - X2 @ theta_hat
    m1 = Z1.T @ r1_final / N
    m2 = Z2.T @ r2_final / N
    m_final = np.concatenate([m1, m2])
    j_stat = float(N * m_final @ W_eff @ m_final)
    p_value = (float(1.0 - chi2.cdf(j_stat, df=n_overid))
               if n_overid > 0 else 1.0)

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
        pf_order=pf_order,
        proxy_order=proxy_order,
        use_controls=use_controls,
    )

    # ---- Firm-year markups from theta and cost share ---------------- #
    if compute_markups:
        # For CD: theta_c = beta_cogs (constant across firms).
        # For TL: theta_c_jt = beta_cogs + 2*beta_cogs2*cogs_jt + beta_kcogs*k_jt
        #         (firm-specific, matches acf_estimator._compute_markups).
        beta_cogs = theta_hat[idx_beta]
        if pf_order == 1:
            theta_c = np.full(N, beta_cogs, dtype=np.float64)
        elif pf_order == 2:
            idx_cogs2 = param_names.index('beta_cogs2')
            idx_kcogs = param_names.index('beta_kcogs')
            beta_cogs2 = theta_hat[idx_cogs2]
            beta_kcogs = theta_hat[idx_kcogs]
            cogs_v = df['cogs'].values.astype(np.float64)
            k_v = df['k'].values.astype(np.float64)
            theta_c = beta_cogs + 2.0 * beta_cogs2 * cogs_v + beta_kcogs * k_v
        else:
            raise NotImplementedError(f'pf_order={pf_order} markup formula')

        markup = theta_c / alpha_data
        out['markups'] = pd.DataFrame({
            'id': df['id'].values,
            'year': df['year'].values,
            'theta_c': theta_c,
            'alphahat': alpha_data,
            'markup': markup,
        })

    return out


def prepare_lagged_panel(df: pd.DataFrame) -> pd.DataFrame:
    """Add polynomial and lagged variables required by wooldridge_gmm.

    Matches the convention used by wooldridge_test.py — produces k2, cogs2,
    kcogs, k3, cogs3, k2cogs, kcogs2 and their one-period lags within firm.
    Input df must be sorted by (id, year).
    """
    df = df.copy()
    if 'k2' not in df.columns:
        df['k2'] = df['k'] ** 2
    if 'cogs2' not in df.columns:
        df['cogs2'] = df['cogs'] ** 2
    if 'kcogs' not in df.columns:
        df['kcogs'] = df['k'] * df['cogs']
    if 'k3' not in df.columns:
        df['k3'] = df['k'] ** 3
    if 'cogs3' not in df.columns:
        df['cogs3'] = df['cogs'] ** 3
    if 'k2cogs' not in df.columns:
        df['k2cogs'] = df['k2'] * df['cogs']
    if 'kcogs2' not in df.columns:
        df['kcogs2'] = df['k'] * df['cogs2']

    lag_vars = [
        'k', 'cogs', 'k2', 'cogs2', 'kcogs', 'k3', 'cogs3',
        'k2cogs', 'kcogs2', 'pp_dummy',
    ]
    df = df.sort_values(['id', 'year']).reset_index(drop=True)
    for col in lag_vars:
        if col in df.columns:
            df[f'{col}_lag'] = df.groupby('id')[col].shift(1)
    return df
