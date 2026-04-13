"""dml_sensitivity.py — Candidate 2: DML Sensitivity Analysis.

Reference: Chernozhukov, Cinelli, Hazlett, Kuchibhotla, and Robins (2019),
"Long Story Short: Omitted Variable Bias in Causal Machine Learning", NBER WP.
Notebook: MetricsMLNotebooks/AC1/python-sensitivity-analysis-with-sensmakr-and-debiased-ml.ipynb

For the partially linear model Y = α D + g(X) + U, Chernozhukov et al. (2019)
show that an unobserved confounder with partial R² (R²_Y, R²_D) in the outcome
and treatment equations produces a bias bounded by

    |Bias|² ≤ [R²_Y × R²_D / (1 - R²_D)] × (σ_ε² / σ_D_res²)

where σ_ε² is the variance of the final PLR residual (Y - ĝ - α̂ D - residual)
and σ_D_res² is the variance of the residualized treatment D - m̂(X).

Reported diagnostics:
    1. **Robustness value (RV)**: the minimum R² (equal in both equations)
       at which the CI crosses zero
    2. **Bias landscape**: 2D contour plot of bias as a function of (R²_Y, R²_D)
    3. **Benchmarks**: partial R² of observed covariates (notably omega_A,
       the ACF productivity) as anchor points for the "how much is realistic"
       question

This supplements the existing Oster (2019) δ* = −6.05 in the paper.

Outputs:
    output/tables/dml_sensitivity.tex         LaTeX table for §6.8.2
    output/figures/dml_sensitivity_contour.pdf  2D contour plot
    output/tables/dml_sensitivity.csv         Raw diagnostics
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dml_core import (
    load_merged_panel, construct_X, firm_demean, cross_fit,
    plr_orthogonal, make_outcome_estimators, make_treatment_estimators,
    sensitivity_bias, robustness_value, TAB_DIR, FIG_DIR,
    DEFAULT_SEED, BASE_COVARIATES,
)
from style_markups import apply_markups_style, MARKUPS_BLUE, MARKUPS_RED, MARKUPS_GREY


def compute_observed_benchmark(df: pd.DataFrame, covariate: str,
                                  y_dm: np.ndarray, d_dm: np.ndarray) -> dict:
    """Partial R² of an observed covariate in the Y and D equations.

    Serves as a benchmark for "how much would an unobserved confounder
    of similar strength change the conclusion?"
    """
    z = firm_demean(df[covariate], df["id"])
    # R²_Y: how much does z explain of y | (year FE)?
    # Simple: corr(z, y)² with year trend removed
    year_dummies = pd.get_dummies(df["year"], drop_first=True).astype(float).values
    # Residualize y_dm and d_dm and z by year dummies
    from sklearn.linear_model import LinearRegression
    lm_y = LinearRegression().fit(year_dummies, y_dm)
    lm_d = LinearRegression().fit(year_dummies, d_dm)
    lm_z = LinearRegression().fit(year_dummies, z)
    y_res = y_dm - lm_y.predict(year_dummies)
    d_res = d_dm - lm_d.predict(year_dummies)
    z_res = z - lm_z.predict(year_dummies)

    def partial_r2(a, b):
        ss = np.var(a - b * (np.sum(a * b) / np.sum(b * b)))
        tot = np.var(a)
        return max(0.0, 1.0 - ss / tot) if tot > 0 else 0.0

    r2_y = partial_r2(y_res, z_res)
    r2_d = partial_r2(d_res, z_res)
    return dict(covariate=covariate, r2_y=float(r2_y), r2_d=float(r2_d))


def run_sensitivity(df: pd.DataFrame) -> dict:
    """Compute the DML sensitivity diagnostics."""
    y_dm = firm_demean(df["log_mu"], df["id"])
    d_dm = firm_demean(df["pp_dummy"].astype(float), df["id"])
    X = construct_X(df, covariates=BASE_COVARIATES, fe_policy="demean")
    X = X.drop(columns=[c for c in X.columns if c.startswith("nace_")])
    firm = df["id"].values

    # Use GB as the reference nuisance estimator (most flexible)
    outcome_models = make_outcome_estimators(seed=DEFAULT_SEED)
    treatment_models = make_treatment_estimators(seed=DEFAULT_SEED, binary=False)
    g_hat = cross_fit(X, y_dm, outcome_models["GB"])
    m_hat = cross_fit(X, d_dm, treatment_models["GB"], classifier=False)

    res = plr_orthogonal(y_dm, d_dm, g_hat, m_hat, firm)
    point = res["point"]
    stderr = res["stderr"]
    resY = res["resY"]
    resD = res["resD"]

    print(f"\n[baseline] DML point estimate: α̂ = {point:.4f} (SE {stderr:.4f})")

    # Robustness value — minimum equal R² that flips significance
    rv = robustness_value(point, stderr, resY, resD, q=0.0, alpha=0.05)
    # RV for q = 1 (CI boundary rather than point estimate)
    rv_ci = robustness_value(point, stderr, resY, resD, q=1.0, alpha=0.05)
    print(f"[RV] equal-R² robustness value (flip sign):     {rv:.4f}")
    print(f"[RV] equal-R² robustness value (flip CI):        {rv_ci:.4f}")

    # Observed-covariate benchmarks: compute partial R² of each known
    # confounder candidate (productivity, market share, etc.)
    benchmarks = []
    for cov in ["omega_A", "mktshare", "foreign", "k", "cogs"]:
        if cov in df.columns:
            b = compute_observed_benchmark(df, cov, y_dm, d_dm)
            # Bias from a confounder as strong as this observed covariate
            b["implied_bias"] = sensitivity_bias(point, resY, resD, b["r2_y"], b["r2_d"])
            b["adj_point_lo"] = point - b["implied_bias"]
            b["adj_point_hi"] = point + b["implied_bias"]
            benchmarks.append(b)
            print(f"[benchmark] {cov:12s}: R²_Y={b['r2_y']:.4f}, R²_D={b['r2_d']:.4f}, "
                  f"bias={b['implied_bias']:.4f}, adjusted α ∈ "
                  f"[{b['adj_point_lo']:.4f}, {b['adj_point_hi']:.4f}]")

    return dict(
        point=point,
        stderr=stderr,
        rv=rv,
        rv_ci=rv_ci,
        benchmarks=benchmarks,
        resY=resY,
        resD=resD,
    )


def plot_bias_landscape(result: dict) -> None:
    """Contour plot of absolute bias as a function of (R²_Y, R²_D)."""
    apply_markups_style()
    point = result["point"]
    stderr = result["stderr"]
    resY = result["resY"]
    resD = result["resD"]

    # Build grid
    grid_dc = np.linspace(0.001, 0.30, 120)
    grid_yc = np.linspace(0.001, 0.30, 120)
    DC, YC = np.meshgrid(grid_dc, grid_yc)
    BIAS = np.zeros_like(DC)
    for i in range(DC.shape[0]):
        for j in range(DC.shape[1]):
            BIAS[i, j] = sensitivity_bias(point, resY, resD, YC[i, j], DC[i, j])

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    # Contour levels at fractions of the point estimate
    levels = [abs(point) * k for k in [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]]
    cs = ax.contour(DC, YC, BIAS, levels=levels, colors=MARKUPS_GREY, linewidths=0.8)
    ax.clabel(cs, inline=True, fontsize=8, fmt=lambda v: f"|bias|={v:.3f}")

    # Highlight the "flip sign" contour
    cs_flip = ax.contour(DC, YC, BIAS, levels=[abs(point)],
                          colors=MARKUPS_RED, linewidths=2.0)
    ax.clabel(cs_flip, inline=True, fontsize=9,
              fmt=lambda v: f"|bias|={v:.3f} (flips sign)")

    # Benchmark scatter points
    for b in result["benchmarks"]:
        ax.scatter(b["r2_d"], b["r2_y"], s=60, color=MARKUPS_BLUE,
                    edgecolor="black", linewidth=0.5, zorder=5)
        ax.annotate(b["covariate"], (b["r2_d"], b["r2_y"]),
                     xytext=(5, 5), textcoords="offset points", fontsize=8)

    # RV diagonal line
    rv = result["rv"]
    ax.plot([0, 0.3], [0, 0.3], color="black", linestyle=":", linewidth=0.5,
             alpha=0.5)
    if rv > 0 and rv < 0.3:
        ax.scatter([rv], [rv], s=100, color=MARKUPS_RED, marker="*",
                    edgecolor="black", zorder=6)
        ax.annotate(f"RV = {rv:.3f}", (rv, rv),
                     xytext=(8, -3), textcoords="offset points",
                     fontsize=10, fontweight="bold")

    ax.set_xlabel(r"Partial $R^2$ of confounder in treatment equation ($R^2_D$)")
    ax.set_ylabel(r"Partial $R^2$ of confounder in outcome equation ($R^2_Y$)")
    ax.set_title(f"DML Sensitivity: Bias Landscape (α̂ = {point:.4f})")
    ax.set_xlim(0, 0.30)
    ax.set_ylim(0, 0.30)

    fig.tight_layout()
    path = FIG_DIR / "dml_sensitivity_contour.pdf"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {FIG_DIR.name}/dml_sensitivity_contour.pdf")


def write_table(result: dict) -> None:
    lines = [
        r"\begin{table}[htbp]\centering",
        r"\caption{DML Sensitivity Analysis for the Procurement Markup Premium}",
        r"\label{tab:dml_sensitivity}",
        r"\begin{threeparttable}",
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"& $R^2_Y$ & $R^2_D$ & Implied $|\text{bias}|$ & Adj.\ $\hat\alpha$ range \\",
        r"\midrule",
        rf"\emph{{DML baseline}} & & & & $\hat\alpha = {result['point']:.4f}$ "
        rf"(SE {result['stderr']:.4f}) \\",
        r"\midrule",
        r"\emph{Observed-covariate benchmarks:} & & & & \\",
    ]
    for b in result["benchmarks"]:
        lines.append(
            rf"\quad {b['covariate']} & {b['r2_y']:.4f} & {b['r2_d']:.4f} & "
            rf"{b['implied_bias']:.4f} & "
            rf"[{b['adj_point_lo']:.4f}, {b['adj_point_hi']:.4f}] \\"
        )
    lines += [
        r"\midrule",
        rf"\multicolumn{{4}}{{l}}{{\emph{{Robustness value (RV):}} minimum equal "
        rf"$R^2_Y = R^2_D$ to flip sign of $\hat\alpha$}} & {result['rv']:.4f} \\",
        rf"\multicolumn{{4}}{{l}}{{\emph{{Robustness value (RV$_{{CI}}$):}} "
        rf"minimum equal $R^2_Y = R^2_D$ to cross zero CI boundary}} & {result['rv_ci']:.4f} \\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\footnotesize",
        r"\item \emph{Notes:} DML sensitivity analysis following "
        r"Chernozhukov et al.\ \cite{ChernozhukovSensitivity2019}. The "
        r"\emph{robustness value} (RV) is the smallest equal partial $R^2$ "
        r"in both the outcome and treatment equations at which an unobserved "
        r"confounder would eliminate the baseline estimate's significance. "
        r"Observed-covariate benchmarks report the partial $R^2$ of known "
        r"pre-treatment variables in the outcome equation ($R^2_Y$: variation "
        r"in $\log \mu$ explained conditional on year fixed effects) and "
        r"treatment equation ($R^2_D$: variation in procurement participation), "
        r"with the implied bias from a confounder of that strength. "
        r"An unobserved confounder must exceed the strongest observed "
        r"covariate ($\omega_A$, the ACF-estimated productivity) by the ratio "
        r"$\text{RV} / R^2(\omega_A)$ to threaten the estimate---complementing "
        r"the single-parameter Oster~\cite{Oster2019} $\delta^* = -6.05$ "
        r"already reported in Section~\ref{sec:productivity_decomp}.",
        r"\end{tablenotes}",
        r"\end{threeparttable}",
        r"\end{table}",
    ]
    path = TAB_DIR / "dml_sensitivity.tex"
    path.write_text("\n".join(lines) + "\n")

    # Save raw CSV
    rows = [dict(metric="point", value=result["point"]),
            dict(metric="stderr", value=result["stderr"]),
            dict(metric="rv", value=result["rv"]),
            dict(metric="rv_ci", value=result["rv_ci"])]
    for b in result["benchmarks"]:
        rows.append(dict(metric=f"benchmark_{b['covariate']}_r2y", value=b["r2_y"]))
        rows.append(dict(metric=f"benchmark_{b['covariate']}_r2d", value=b["r2_d"]))
        rows.append(dict(metric=f"benchmark_{b['covariate']}_bias", value=b["implied_bias"]))
    pd.DataFrame(rows).to_csv(TAB_DIR / "dml_sensitivity.csv", index=False)
    print(f"[write] {TAB_DIR.name}/dml_sensitivity.tex + .csv")


def main() -> None:
    df = load_merged_panel()
    result = run_sensitivity(df)
    plot_bias_landscape(result)
    write_table(result)
    print("\n[done]")


if __name__ == "__main__":
    main()
