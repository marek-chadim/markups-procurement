"""dml_cate.py — Candidate 3: DML Conditional Average Treatment Effect by Subgroup.

Reference: Chernozhukov et al. (2018) DML Econometrics Journal 21(1): C1–C68,
with AIPW-style doubly robust signal from MetricsMLNotebooks/T/dml-for-conditional-
average-treatment-effect.irnb cell 13.

Given cross-fitted nuisances (ĝ, m̂) from the full sample, the conditional ATE
for subgroup g is

    α̂_g = Σ_{i ∈ g} (Y_i - ĝ_i)(D_i - m̂_i) / Σ_{i ∈ g} (D_i - m̂_i)²

which is the Neyman-orthogonal score restricted to the subgroup. This is more
efficient than re-running DML on each subset because the nuisance functions
borrow strength from the full sample.

Subgroups reported:
    1. NACE 41 buildings / 42 civil eng. / 43 specialized
    2. Pre-2012 reform / post-2012 reform
    3. NACE × reform era (2 × 3 = 6 cells)

Outputs:
    output/tables/dml_cate_heterogeneity.tex   Supplements §7 Table 24
    output/tables/dml_cate_heterogeneity.csv   Raw subgroup coefficients
"""

from __future__ import annotations

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from dml_core import (
    load_merged_panel, construct_X, firm_demean, cross_fit,
    make_outcome_estimators, make_treatment_estimators, TAB_DIR,
    DEFAULT_SEED, BASE_COVARIATES, cluster_bootstrap, stars,
)


def compute_cate_subgroup(y_dm: np.ndarray, d_dm: np.ndarray,
                            g_hat: np.ndarray, m_hat: np.ndarray,
                            mask: np.ndarray, firm_ids: np.ndarray) -> dict:
    """Orthogonal score restricted to a subgroup, with firm-cluster bootstrap SE."""
    resY = y_dm - g_hat
    resD = d_dm - m_hat
    sub_resY = resY[mask]
    sub_resD = resD[mask]
    sub_firm = firm_ids[mask]

    denom = float(np.sum(sub_resD ** 2))
    if denom <= 0:
        return dict(point=float("nan"), se=float("nan"), n=0)
    point = float(np.sum(sub_resY * sub_resD) / denom)

    # Cluster bootstrap SE within the subgroup
    def point_fn(idx, rY=sub_resY, rD=sub_resD):
        sr = rY[idx]
        sd = rD[idx]
        den = np.sum(sd ** 2)
        return float(np.sum(sr * sd) / den) if den > 0 else float("nan")

    boot = cluster_bootstrap(point_fn, sub_firm, n_rep=200, seed=DEFAULT_SEED,
                              theta_hat=point)
    return dict(
        point=point,
        se=boot["se"],
        ci_lo=boot["ci_lo"],
        ci_hi=boot["ci_hi"],
        n=int(mask.sum()),
        n_firms=int(pd.unique(sub_firm).size),
    )


def run_cate(df: pd.DataFrame, ml_name: str = "GB") -> pd.DataFrame:
    """Compute full-sample nuisance with one ML estimator, then CATE by subgroup."""
    y_dm = firm_demean(df["log_mu"], df["id"])
    d_dm = firm_demean(df["pp_dummy"].astype(float), df["id"])
    X = construct_X(df, covariates=BASE_COVARIATES, fe_policy="demean")
    X = X.drop(columns=[c for c in X.columns if c.startswith("nace_")])
    firm = df["id"].values

    om = make_outcome_estimators(seed=DEFAULT_SEED)[ml_name]
    tm = make_treatment_estimators(seed=DEFAULT_SEED, binary=False)[ml_name]
    g_hat = cross_fit(X, y_dm, om)
    m_hat = cross_fit(X, d_dm, tm, classifier=False)

    rows = []

    # Pooled (sanity check: should match dml_premium C1 Spec A for this ml_name)
    r = compute_cate_subgroup(y_dm, d_dm, g_hat, m_hat,
                                np.ones(len(df), dtype=bool), firm)
    rows.append(dict(panel="Pooled", subgroup="All firms", **r))

    # Panel A: by NACE
    for nace, label in [(41, "NACE 41 buildings"),
                         (42, "NACE 42 civil eng."),
                         (43, "NACE 43 specialized")]:
        mask = (df["nace2"] == nace).values
        r = compute_cate_subgroup(y_dm, d_dm, g_hat, m_hat, mask, firm)
        rows.append(dict(panel="NACE", subgroup=label, **r))

    # Panel B: reform era
    for yr_range, label in [((2006, 2011), "Pre-2012 reform"),
                              ((2012, 2021), "Post-2012 reform")]:
        mask = (df["year"].between(yr_range[0], yr_range[1])).values
        r = compute_cate_subgroup(y_dm, d_dm, g_hat, m_hat, mask, firm)
        rows.append(dict(panel="Reform", subgroup=label, **r))

    # Panel C: NACE × reform era
    for nace, nace_label in [(41, "41"), (42, "42"), (43, "43")]:
        for yr_range, era_label in [((2006, 2011), "pre"), ((2012, 2021), "post")]:
            mask = ((df["nace2"] == nace) &
                    df["year"].between(yr_range[0], yr_range[1])).values
            r = compute_cate_subgroup(y_dm, d_dm, g_hat, m_hat, mask, firm)
            rows.append(dict(panel="NACE × Reform",
                             subgroup=f"NACE {nace_label} {era_label}", **r))

    return pd.DataFrame(rows)


def write_table(results: pd.DataFrame, ml_name: str) -> None:
    lines = [
        r"\begin{table}[htbp]\centering",
        rf"\caption{{DML Conditional Average Treatment Effects by Subgroup ({ml_name} Nuisance)}}",
        r"\label{tab:dml_cate}",
        r"\begin{threeparttable}",
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Subgroup & $\hat\alpha_g$ & Boot.\ SE & 95\% CI & $N$ \\",
        r"\midrule",
    ]
    for panel in ["Pooled", "NACE", "Reform", "NACE × Reform"]:
        panel_rows = results[results["panel"] == panel]
        if len(panel_rows) == 0:
            continue
        lines.append(rf"\emph{{{panel}}} & & & & \\")
        for _, r in panel_rows.iterrows():
            ci = f"[{r['ci_lo']:.4f}, {r['ci_hi']:.4f}]" if not pd.isna(r.get('ci_lo', np.nan)) else ""
            lines.append(
                rf"\quad {r['subgroup']} & {r['point']:.4f} & "
                rf"({r['se']:.4f}) & {ci} & {int(r['n']):,} \\"
            )
        lines.append(r"\midrule")
    lines = lines[:-1] + [
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\footnotesize",
        rf"\item \emph{{Notes:}} Each row reports the DML conditional average "
        r"treatment effect $\hat\alpha_g = \sum_{i \in g}(Y_i - \hat g(X_i)) "
        r"(D_i - \hat m(X_i)) / \sum_{i \in g}(D_i - \hat m(X_i))^2$ restricted "
        r"to the indicated subgroup, using nuisance functions $\hat g$ and "
        rf"$\hat m$ cross-fitted on the \emph{{full}} sample with the {ml_name} "
        r"estimator. Standard errors are firm-clustered bootstrap with $B=200$ "
        r"replications. The pooled row confirms consistency with "
        r"Table~\ref{tab:dml_premium}. Panel A decomposes by NACE 2-digit "
        r"sub-industry; Panel B by the 2012 Act~55/2012 single-bidding ban; "
        r"Panel C crosses the two to reveal NACE-specific reform incidence.",
        r"\end{tablenotes}",
        r"\end{threeparttable}",
        r"\end{table}",
    ]
    path = TAB_DIR / "dml_cate_heterogeneity.tex"
    path.write_text("\n".join(lines) + "\n")
    results.to_csv(TAB_DIR / "dml_cate_heterogeneity.csv", index=False)
    print(f"[write] {TAB_DIR.name}/dml_cate_heterogeneity.tex + .csv")


def main() -> None:
    df = load_merged_panel()
    print("\n[DML-CATE] cross-fitting nuisances on full sample (GB)...")
    results = run_cate(df, ml_name="GB")
    print("\n[DML-CATE] subgroup results:")
    print(results[["panel", "subgroup", "point", "se", "n"]].to_string(index=False))

    write_table(results, "GB")
    print("\n[done]")


if __name__ == "__main__":
    main()
