"""ags_sensitivity.py — Python port of the AGS (2017 QJE) sensitivity-matrix
library, translated from the reference Matlab implementation at
`references/replications/sensitivity_replication/.../external/lib/matlab/trans/m/`.

The port covers the three functions we need to deploy in the companion
theory note `abgrs_exposition.tex` for operationalizing Proposition 1's
Lipschitz constant $C_2$ into a concrete observation-level diagnostic:

  get_sensitivity            — Λ = −(J'WJ)^{-1}(J'W)
  build_sigma                — (Σ_tt, Σ_tg, Σ_gg) joint asymptotic vcov
  get_standardized_sensitivity — scale by parameter and moment SEs
  compute_composite_lambda   — chain-rule sensitivity of a scalar causal
                               summary τ*(θ) to each Stage A moment

Reference: Andrews, I., Gentzkow, M., and Shapiro, J.M. (2017).
"Measuring the Sensitivity of Parameter Estimates to Estimation Moments."
QJE 132(4): 1553–1592.

Matlab originals (line counts):
  trans/m/get_sensitivity.m (7)
  trans/m/build_sigma.m (7)
  trans/m/get_standardized_sensitivity.m (10)
  trans/m/get_transform_sensitivity.m (10) — transforms moments, not used here
"""

from __future__ import annotations

import numpy as np


def get_sensitivity(
    jacobian: np.ndarray,
    weight: np.ndarray,
    allow_pinv: bool = True,
) -> np.ndarray:
    r"""AGS (2017) sensitivity matrix for a minimum-distance estimator.

    Given the moment-function Jacobian $J = \partial \bar m(\theta)/\partial\theta'$
    (M rows × P cols, where M = #moments, P = #params) and the GMM weighting
    matrix $W$ (M × M, symmetric positive-definite), the AGS sensitivity is

    .. math::
        \Lambda = -(J' W J)^{-1} (J' W)

    of shape P × M. Entry $\Lambda_{p,m}$ is the derivative of $\hat\theta_p$
    with respect to the $m$-th sample moment, evaluated at $\theta_0$.

    If the estimator is under-identified in the strict Jacobian sense
    (rank $(J'WJ)$ < P), we fall back to the Moore-Penrose pseudoinverse,
    which reports the sensitivity for the identified subspace and zero for
    the null directions. This is the natural fallback when the Jacobian is
    rank-deficient because of a weakly identified parameter — our translog
    curvature terms are the motivating case.

    Matlab reference: `trans/m/get_sensitivity.m` (7 lines).
    """
    if jacobian.ndim != 2 or weight.ndim != 2:
        raise ValueError("jacobian must be M×P and weight must be M×M")
    JWJ = jacobian.T @ weight @ jacobian
    rank = np.linalg.matrix_rank(JWJ)
    if rank < JWJ.shape[0]:
        if not allow_pinv:
            raise np.linalg.LinAlgError(
                "J'WJ is rank-deficient; sensitivity matrix is undefined"
            )
        return -np.linalg.pinv(JWJ) @ (jacobian.T @ weight)
    return -np.linalg.solve(JWJ, jacobian.T @ weight)


def build_sigma(
    Lambda: np.ndarray,
    param_vcov: np.ndarray,
    moment_vcov: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""Joint asymptotic variance-covariance matrix of $(\hat\theta, \bar m_N)$.

    Given the AGS sensitivity $\Lambda$ (P×M), parameter vcov $\Sigma_{tt}$
    (P×P), and moment vcov $\Sigma_{gg}$ (M×M), returns:

    - $\Sigma_{tt}$ — parameter vcov (P×P)
    - $\Sigma_{tg} = \Lambda \Sigma_{gg}$ — parameter-moment cross vcov (P×M)
    - $\Sigma_{gg}$ — moment vcov (M×M)
    - $\Sigma$ — full block matrix (P+M × P+M)

    Matlab reference: `trans/m/build_sigma.m` (7 lines).
    """
    Sigma_tt = param_vcov
    Sigma_tg = Lambda @ moment_vcov
    Sigma_gg = moment_vcov
    top = np.hstack([Sigma_tt, Sigma_tg])
    bot = np.hstack([Sigma_tg.T, Sigma_gg])
    Sigma = np.vstack([top, bot])
    return Sigma, Sigma_tt, Sigma_tg, Sigma_gg


def get_standardized_sensitivity(
    Lambda: np.ndarray,
    se_param: np.ndarray,
    se_mom: np.ndarray,
) -> np.ndarray:
    r"""Standardized AGS sensitivity.

    Each entry $\tilde\Lambda_{p,m} = \Lambda_{p,m} \cdot \sigma_m / \sigma_p$,
    where $\sigma_p$ is the parameter SE and $\sigma_m$ is the moment SE.
    Standardized entries are unit-free and report the expected change in
    $\hat\theta_p$ (in parameter standard deviations) induced by a one-SD
    deviation in the $m$-th moment.

    Matlab reference: `trans/m/get_standardized_sensitivity.m`.
    """
    if se_param.ndim != 1 or se_mom.ndim != 1:
        raise ValueError("se_param and se_mom must be 1-D arrays")
    return Lambda * se_mom[np.newaxis, :] / se_param[:, np.newaxis]


def compute_composite_lambda(
    Lambda: np.ndarray,
    psi: np.ndarray,
) -> np.ndarray:
    r"""AGS sensitivity of a scalar causal summary to each Stage A moment.

    Chain rule: if $\tau^*(\theta) = \psi(\theta)$ is a scalar smooth functional
    of $\theta$ with gradient $\psi = \nabla_\theta\tau^*(\theta_0) \in \mathbb{R}^{1\times P}$,
    then to first order in the sample-moment deviation,

    .. math::
        \hat\tau^* - \tau^*(G) \;\approx\; \psi \cdot (\hat\theta - \theta(G))
                               \;\approx\; (\psi \Lambda)\,(\hat m_N - E[m]),

    so the *composite* sensitivity is $\Lambda^\text{comp} = \psi \Lambda$
    (1×M row vector).

    This function implements the operation the AGS Matlab library calls
    the "transform sensitivity of a functional of parameters", which is
    distinct from `get_transform_sensitivity.m` (which transforms moments,
    not parameters).

    Returns a (1 × M) row vector.
    """
    if psi.ndim == 1:
        psi = psi.reshape(1, -1)
    if psi.shape[1] != Lambda.shape[0]:
        raise ValueError(
            f"psi must be 1×P where P={Lambda.shape[0]}; got psi of shape {psi.shape}"
        )
    return psi @ Lambda


def gmm_sandwich_vcov(
    jacobian: np.ndarray,
    weight: np.ndarray,
    moment_vcov: np.ndarray,
    n: int,
) -> np.ndarray:
    r"""Standard GMM sandwich parameter covariance.

    $\hat V_{\theta} = (1/N) (J'WJ)^{-1} J'W \Omega W J (J'WJ)^{-1}$

    This is the complement to `get_sensitivity`: given the same ingredients,
    it returns the asymptotic parameter vcov matrix in the sandwich form.
    Under identity weighting $W = I$, reduces to $(1/N)\Lambda \Omega \Lambda'$.

    Matlab reference: `trans/m/get_vcov.m`.
    """
    JWJ_inv = np.linalg.inv(jacobian.T @ weight @ jacobian)
    meat = jacobian.T @ weight @ moment_vcov @ weight @ jacobian
    return JWJ_inv @ meat @ JWJ_inv / n
