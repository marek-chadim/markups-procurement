"""dml_premium.py — Candidate 1: Partially Linear DML for the 14% Procurement Premium.

Reference: Chernozhukov et al. (2018) *Double/Debiased Machine Learning for
Treatment and Structural Parameters*, Econometrics Journal 21(1): C1–C68.
Notebook: MetricsMLNotebooks/PM4/python_debiased_ml_for_partially_linear_model_growth.ipynb

Model:
    log(μ_it) = α · pp_dummy_it + g(X_it) + firm_i + ε_it

Estimation:
    1. Within-firm demean Y = log(markup_A), D = pp_dummy, and continuous X
    2. Cross-fit g(X) = E[Y|X] and m(X) = E[D|X] via ML (Lasso/RF/GB)
    3. Orthogonal score: α̂ = Σ(Y - ĝ)(D - m̂) / Σ(D - m̂)²  with firm-clustered SE

Two specifications are reported:
    Spec A — BASE_COVARIATES only (mktshare, foreign + year dummies)
             Headline DML premium; expected ~0.13
    Spec B — + omega_A (ACF productivity control)
             Productivity-decomposed premium; expected ~0.03 matching the paper's
             line 341 claim that "productivity explains 76% of the procurement-
             markup association"

Outputs:
    output/tables/dml_premium.tex      LaTeX table for §6.8.1 of the paper
    output/tables/dml_premium.csv      Raw coefficients
    output/figures/dml_premium.pdf     Forest plot of all estimates
"""

from __future__ import annotations
from pathlib import Path

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from dml_core import (
    load_merged_panel, construct_X, firm_demean, cross_fit,
    plr_orthogonal, make_outcome_estimators, make_treatment_estimators,
    cluster_bootstrap, stars, BASE_COVARIATES, PRODUCTIVITY_COVARIATES,
    TAB_DIR, FIG_DIR, DEFAULT_SEED,
)
from style_markups import apply_markups_style, MARKUPS_BLUE, MARKUPS_RED, MARKUPS_GREY


# ---------------------------------------------------------------------------
# DML PLR estimator with cluster SE
# ---------------------------------------------------------------------------

def dml_plr(df: pd.DataFrame, covariates: list[str], label: str) -> pd.DataFrame:
    """Run DML-PLR across Lasso, RF, GB nuisance estimators.

    Returns a DataFrame with 3 rows, one per ML method, reporting:
        point, stderr (asymptotic), boot_se, ci_lo, ci_hi, t, n
    """
    # Demean outcome and treatment by firm
    y_dm = firm_demean(df["log_mu"], df["id"])
    d_dm = firm_demean(df["pp_dummy"].astype(float), df["id"])

    # Build X matrix — pass custom covariates with firm demeaning
    X = construct_X(df, covariates=covariates, fe_policy="demean")

    # Drop the NACE dummies from X (they are collinear with firm after demeaning)
    nace_cols = [c for c in X.columns if c.startswith("nace_")]
    X = X.drop(columns=nace_cols)

    outcome_models = make_outcome_estimators(seed=DEFAULT_SEED)
    # Treatment nuisance: after demeaning, d_dm is continuous in [-1, 1],
    # so use regressor flavor, not classifier.
    treatment_models = make_treatment_estimators(seed=DEFAULT_SEED, binary=False)

    firm = df["id"].values
    rows = []

    for name in ["Lasso", "RF", "GB"]:
        g_hat = cross_fit(X, y_dm, outcome_models[name])
        m_hat = cross_fit(X, d_dm, treatment_models[name], classifier=False)
        res = plr_orthogonal(y_dm, d_dm, g_hat, m_hat, firm)

        # Cluster bootstrap SE — closure capturing g_hat, m_hat
        def point_fn(idx, y=y_dm, d=d_dm, g=g_hat, m=m_hat):
            sY = y[idx] - g[idx]
            sD = d[idx] - m[idx]
            return float(np.sum(sY * sD) / np.sum(sD ** 2))

        boot = cluster_bootstrap(point_fn, firm, n_rep=200, seed=DEFAULT_SEED,
                                  theta_hat=res["point"])

        rows.append(dict(
            spec=label,
            ml=name,
            point=res["point"],
            asym_se=res["stderr"],
            boot_se=boot["se"],
            ci_lo=res["ci_lo"],
            ci_hi=res["ci_hi"],
            boot_ci_lo=boot["ci_lo"],
            boot_ci_hi=boot["ci_hi"],
            t=res["t"],
            n=res["n"],
        ))

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Output: LaTeX table
# ---------------------------------------------------------------------------

def write_table(results: pd.DataFrame, ols_a: float, ols_b: float) -> None:
    """Write the DML premium LaTeX table for §6.8.1."""
    lines = [
        r"\begin{table}[htbp]\centering",
        r"\caption{Double/Debiased Machine Learning Estimates of the Procurement Markup Premium}",
        r"\label{tab:dml_premium}",
        r"\begin{threeparttable}",
        r"\begin{tabular}{llrrrc}",
        r"\toprule",
        r"Specification & ML nuisance & $\hat\alpha$ & Asym.\ SE & Boot.\ SE & 95\% CI \\",
        r"\midrule",
        rf"\emph{{Spec~A}} (baseline) & OLS (reference) & {ols_a:.4f} & & & \\",
    ]
    for _, r in results[results["spec"] == "Spec A"].iterrows():
        from scipy.stats import norm
        p = 2 * norm.sf(abs(r["t"]))
        lines.append(
            rf"& {r['ml']} & {r['point']:.4f}{stars(p)} & "
            rf"({r['asym_se']:.4f}) & ({r['boot_se']:.4f}) & "
            rf"[{r['ci_lo']:.4f}, {r['ci_hi']:.4f}] \\"
        )
    lines += [
        r"\midrule",
        rf"\emph{{Spec~B}} ($+$ ACF productivity) & OLS (reference) & {ols_b:.4f} & & & \\",
    ]
    for _, r in results[results["spec"] == "Spec B"].iterrows():
        from scipy.stats import norm
        p = 2 * norm.sf(abs(r["t"]))
        lines.append(
            rf"& {r['ml']} & {r['point']:.4f}{stars(p)} & "
            rf"({r['asym_se']:.4f}) & ({r['boot_se']:.4f}) & "
            rf"[{r['ci_lo']:.4f}, {r['ci_hi']:.4f}] \\"
        )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\footnotesize",
        r"\item \emph{Notes:} Each row reports the DML estimate of "
        r"$\alpha$ in $\log(\mu_{it}) = \alpha \cdot pp_{it} + g(X_{it}) + "
        r"\text{firm}_i + \varepsilon_{it}$ using cross-fitted ML nuisance "
        r"functions $g(X) = E[Y|X]$ and $m(X) = E[D|X]$ and the "
        r"Neyman-orthogonal score of Chernozhukov et al.\ "
        r"\cite{ChernozhukovDML2018}. Cross-fitting uses 5 folds with "
        r"seed 42; firm fixed effects are absorbed by within-firm demeaning. "
        r"Spec~A conditions on exogenous controls (market share, foreign "
        r"ownership, year dummies); Spec~B adds the ACF-estimated log "
        r"productivity $\hat\omega_{it}$ as a flexible nonparametric "
        r"confounder. Asymptotic standard errors are clustered by firm "
        r"via the influence-function sandwich; bootstrap standard errors "
        r"are firm-block bootstraps with $B = 200$ replications. "
        r"$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$.",
        r"\end{tablenotes}",
        r"\end{threeparttable}",
        r"\end{table}",
    ]
    path = TAB_DIR / "dml_premium.tex"
    path.write_text("\n".join(lines) + "\n")
    results.to_csv(TAB_DIR / "dml_premium.csv", index=False)
    print(f"[write] {TAB_DIR.name}/dml_premium.tex + .csv")


def plot_forest(results: pd.DataFrame, ols_a: float, ols_b: float) -> None:
    """Forest plot of DML estimates vs OLS baseline."""
    apply_markups_style()
    fig, ax = plt.subplots(figsize=(7.5, 4.5))

    labels = []
    points = []
    errs_lo = []
    errs_hi = []
    colors = []

    # OLS reference lines
    ax.axvline(ols_a, color=MARKUPS_GREY, linestyle="--", linewidth=1.0,
                label=f"OLS Spec A ({ols_a:.4f})")
    ax.axvline(ols_b, color=MARKUPS_GREY, linestyle=":", linewidth=1.0,
                label=f"OLS Spec B ({ols_b:.4f})")

    for spec_label, color in [("Spec A", MARKUPS_BLUE), ("Spec B", MARKUPS_RED)]:
        sub = results[results["spec"] == spec_label]
        for _, r in sub.iterrows():
            labels.append(f"{spec_label}: {r['ml']}")
            points.append(r["point"])
            errs_lo.append(r["point"] - r["ci_lo"])
            errs_hi.append(r["ci_hi"] - r["point"])
            colors.append(color)

    y = np.arange(len(labels))
    ax.errorbar(points, y, xerr=[errs_lo, errs_hi], fmt="o",
                 color="black", ecolor="black", markersize=6,
                 capsize=3, linestyle="none")
    for i, c in enumerate(colors):
        ax.plot(points[i], y[i], "o", color=c, markersize=8, zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"$\hat\alpha$ (DML-PLR premium)")
    ax.set_title("DML Partially Linear Premium: Point Estimates and 95\\% CIs")
    ax.axvline(0, color="black", linewidth=0.4)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    path = FIG_DIR / "dml_premium.pdf"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {FIG_DIR.name}/dml_premium.pdf")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    df = load_merged_panel()

    # Compute OLS baselines for the two specifications
    import statsmodels.formula.api as smf
    ols_a = smf.ols("log_mu ~ pp_dummy + C(year) + C(id)", data=df).fit().params["pp_dummy"]
    ols_b = smf.ols("log_mu ~ pp_dummy + omega_A + C(year) + C(id)", data=df).fit().params["pp_dummy"]

    print(f"[OLS baselines]  Spec A (no productivity): {ols_a:.4f}")
    print(f"[OLS baselines]  Spec B (with omega_A):    {ols_b:.4f}")

    # Spec A: exogenous controls only
    print("\n[DML Spec A] mktshare, foreign, year FE, firm FE (demeaned)")
    res_a = dml_plr(df, BASE_COVARIATES, "Spec A")
    print(res_a[["ml", "point", "asym_se", "ci_lo", "ci_hi", "t"]].to_string(index=False))

    # Spec B: add omega_A
    print("\n[DML Spec B] + omega_A (productivity decomposition)")
    res_b = dml_plr(df, PRODUCTIVITY_COVARIATES, "Spec B")
    print(res_b[["ml", "point", "asym_se", "ci_lo", "ci_hi", "t"]].to_string(index=False))

    results = pd.concat([res_a, res_b], ignore_index=True)
    write_table(results, ols_a, ols_b)
    plot_forest(results, ols_a, ols_b)
    print("\n[done]")


if __name__ == "__main__":
    main()
