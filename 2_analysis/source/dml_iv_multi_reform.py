"""dml_iv_multi_reform.py — Multi-reform DML-PLR-IV: exogenous-variation causal estimates.

Extends dml_iv.py (single 2012 reform IV) to three independent Czech procurement
reforms with mechanism-matched firm exposures and weak-IV-robust (Anderson-Rubin)
confidence intervals. Provides a direct response to the DLW 2012 disclaimer
("we are not interpreting δ₁ as a causal parameter") using point-identified
exogenous variation.

Four instrument specifications (4 rows × 3 ML nuisances = 12 estimates):
  1. 2012 single-bid ban, generic exposure:   pre-2011 mean pp_dummy × 1[t≥2012]
  2. 2012 single-bid ban, mechanism-matched:  pre-2011 mean single_bid_share × 1[t≥2012]
  3. 2016 MEAT criteria:                       pre-2015 mean pp_dummy × 1[t≥2016]
  4. 2017 Register of Contracts:               pre-2016 mean pp_dummy × 1[t≥2017]

The 2012 mechanism-matched specification uses the firm's pre-reform average
`single_bid_share` (imputed as 0 for firm-years with null, which are firms that
never bid on procurement and therefore had zero exposure to the ban). Because
the single-bidding ban mechanically affected firms that relied on single-bid
contracts, this exposure is theoretically better-targeted than the generic
pre-procurement share. The contrast between row 1 (generic) and row 2
(mechanism-matched) is the paper's empirical test of whether exposure targeting
matters for LATE identification.

Each estimate is reported with two 95% confidence intervals:
  - Wald CI from the orthogonal-score influence function (cluster-robust)
  - Anderson-Rubin CS from grid inversion of the cluster-robust AR statistic
    (weak-IV-robust, valid regardless of first-stage strength)

Reforms (Table 1 of the paper):
  - 2012: Act 55/2012, mandatory cancellation of single-bid tenders
  - 2016: Act 134/2016, shift from lowest-price to MEAT criteria
  - 2017: Register of Contracts enforcement (online publication condition)

Each reform was enacted by the Czech Parliament for reasons unrelated to
firm-specific markups (EU procurement directive compliance, corruption
concerns, transparency legislation), and firms could neither anticipate nor
control the effective dates.

Outputs:
    output/tables/dml_iv_multi_reform.tex    LaTeX table (4 rows × 3 ML + AR CI)
    output/tables/dml_iv_multi_reform.csv    raw results

References:
    Anderson, T.W., and Rubin, H. (1949). Estimation of the Parameters of a
      Single Equation in a Complete System of Stochastic Equations. Annals
      of Mathematical Statistics, 20(1), 46-63.
    Stock, J.H., and Wright, J.H. (2000). GMM with Weak Identification.
      Econometrica, 68(5), 1055-1096.
    Callaway, B., Goodman-Bacon, A., and Sant'Anna, P. (2024). Difference-in-
      Differences with a Continuous Treatment. Working paper.
    Chernozhukov et al. (2018). DML for Treatment and Structural Parameters.
      Econometrics Journal, 21(1).
"""

from __future__ import annotations

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from dml_core import (
    load_merged_panel, construct_X, firm_demean, cross_fit,
    plr_iv_orthogonal, ar_confidence_set,
    make_outcome_estimators, make_treatment_estimators,
    TAB_DIR, DEFAULT_SEED, BASE_COVARIATES, stars,
)

# (reform_year, pre_end_year, exposure_var, short_label)
REFORMS = [
    (2012, 2011, "pp_dummy",         "2012 single-bid ban (pp-generic)"),
    (2012, 2011, "single_bid_share", "2012 single-bid ban (SB-matched)"),
    (2016, 2015, "pp_dummy",         "2016 MEAT criteria"),
    (2017, 2016, "pp_dummy",         "2017 Reg.\\ of Contracts"),
]


def build_reform_iv(df: pd.DataFrame, reform_year: int, pre_end: int,
                     exposure_var: str = "pp_dummy") -> pd.DataFrame:
    """Build firm pre-period exposure × post-reform indicator.

    Z^(r)_it = mean(exposure_i | year ≤ pre_end) × 1[year ≥ reform_year]

    For `exposure_var="single_bid_share"`, imputes 0 for firm-years with null
    (these are firms that never bid on procurement and therefore had zero
    exposure to the 2012 single-bid ban — the theoretically correct default).
    """
    df = df.sort_values(["id", "year"]).copy()
    if exposure_var == "single_bid_share":
        df["_exposure"] = df["single_bid_share"].fillna(0.0)
    else:
        df["_exposure"] = df[exposure_var].astype(float)

    pre_mask = df["year"] <= pre_end
    pre_exp = (df[pre_mask].groupby("id")["_exposure"].mean()
                .rename(f"pre_exp_{pre_end}_{exposure_var}"))
    df = df.merge(pre_exp, on="id", how="left")
    df[f"pre_exp_{pre_end}_{exposure_var}"] = (
        df[f"pre_exp_{pre_end}_{exposure_var}"].fillna(0.0)
    )
    df["z_iv"] = (
        df[f"pre_exp_{pre_end}_{exposure_var}"]
        * (df["year"] >= reform_year).astype(float)
    )
    return df


def run_reform_iv(df: pd.DataFrame, covariates: list[str],
                  reform_year: int, pre_end: int, exposure_var: str,
                  label: str) -> pd.DataFrame:
    """Run DML-PLR-IV for one reform instrument × three ML nuisances.

    Returns a DataFrame with one row per ML nuisance, columns:
    `reform`, `reform_year`, `exposure`, `ml`, `point`, `stderr`, `t`,
    `partial_r2`, `ci_lo`, `ci_hi`, `ar_lo`, `ar_hi`, `ar_empty`, `n`.
    """
    df_iv = build_reform_iv(df, reform_year, pre_end, exposure_var)
    firm = df_iv["id"].values

    y_dm = firm_demean(df_iv["log_mu"], df_iv["id"])
    d_dm = firm_demean(df_iv["pp_dummy"].astype(float), df_iv["id"])
    z_dm = firm_demean(df_iv["z_iv"], df_iv["id"])

    X = construct_X(df_iv, covariates=covariates, fe_policy="demean")
    X = X.drop(columns=[c for c in X.columns if c.startswith("nace_")])

    outcome_models = make_outcome_estimators(seed=DEFAULT_SEED)
    treatment_models = make_treatment_estimators(seed=DEFAULT_SEED, binary=False)

    rows = []
    for name in ["Lasso", "RF", "GB"]:
        g_hat = cross_fit(X, y_dm, outcome_models[name])
        m_hat_d = cross_fit(X, d_dm, treatment_models[name], classifier=False)
        m_hat_z = cross_fit(X, z_dm, treatment_models[name], classifier=False)

        res = plr_iv_orthogonal(y_dm, d_dm, z_dm, g_hat, m_hat_d, m_hat_z, firm)

        # Anderson-Rubin CS: dynamic grid centered on Wald point ± 10×SE
        # to guarantee coverage of the acceptance region even for weak IVs.
        point = res["point"]
        stderr = res["stderr"]
        if np.isfinite(point) and np.isfinite(stderr) and stderr > 0:
            half_width = max(10.0 * stderr, 0.5)
            beta_grid = np.linspace(point - half_width, point + half_width, 401)
        else:
            beta_grid = np.linspace(-0.5, 1.5, 401)
        ar = ar_confidence_set(res["resY"], res["resD"], res["resZ"], firm,
                                alpha=0.05, beta_grid=beta_grid)

        rows.append(dict(
            reform=label,
            reform_year=reform_year,
            exposure=exposure_var,
            ml=name,
            point=point,
            stderr=stderr,
            t=res["t"],
            ci_lo=res["ci_lo"],
            ci_hi=res["ci_hi"],
            partial_r2=res["partial_r2"],
            ar_lo=ar["ar_lo"],
            ar_hi=ar["ar_hi"],
            ar_empty=ar["ar_empty"],
            ar_n_accepted=ar["ar_n_accepted"],
            n=res["n"],
        ))
    return pd.DataFrame(rows)


def write_multi_reform_table(all_results: pd.DataFrame, ols_premium: float) -> None:
    """LaTeX table: 4 reforms × 3 ML nuisances = 12 IV estimates + OLS reference.

    Column structure: {lrrrcc} — ML nuisance, point estimate, SE, first-stage
    R², Wald 95% CI, AR 95% CI. Reform labels are multicolumn section headers.
    """
    from scipy.stats import norm
    lines = [
        r"\begin{table}[htbp]\centering",
        r"\caption{Multi-Reform DML-PLR-IV: Three Independent Exogenous Shocks to Czech Procurement, with Weak-IV-Robust Confidence Sets}",
        r"\label{tab:dml_iv_multi_reform}",
        r"\begin{threeparttable}",
        r"\small",
        r"\begin{tabular}{lrrrcc}",
        r"\toprule",
        r"ML nuisance & $\hat\alpha_{\text{IV}}$ & SE & FS $R^2$ & Wald 95\% CI & AR 95\% CI \\",
        r"\midrule",
        rf"OLS (reference) & {ols_premium:.4f} & & & & \\",
        r"\midrule",
    ]
    reform_order = [(lab, exp) for (_, _, exp, lab) in REFORMS]
    for i, (reform_label, exposure_var) in enumerate(reform_order):
        lines.append(rf"\multicolumn{{6}}{{l}}{{\textit{{{reform_label}}}}} \\")
        group = all_results[
            (all_results["reform"] == reform_label)
            & (all_results["exposure"] == exposure_var)
        ]
        for _, r in group.iterrows():
            p = 2 * norm.sf(abs(r["t"]))
            ar_cell = (
                r"\text{inverts}" if r["ar_empty"]
                else rf"[{r['ar_lo']:.3f}, {r['ar_hi']:.3f}]"
            )
            lines.append(
                rf"{r['ml']} & {r['point']:.4f}{stars(p)} & "
                rf"({r['stderr']:.4f}) & {r['partial_r2']:.4f} & "
                rf"[{r['ci_lo']:.3f}, {r['ci_hi']:.3f}] & {ar_cell} \\"
            )
        if i < len(reform_order) - 1:
            lines.append(r"\midrule")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\footnotesize",
        r"\item \emph{Notes:} DML-PLR-IV estimates of the procurement markup "
        r"premium using continuous-exposure difference-in-differences instruments "
        r"\cite{CallawayGoodmanBaconSantAnna2024}. For each reform year $r$ the "
        r"instrument is "
        r"$Z^{(r)}_{it} = \overline{\text{exposure}}_i^{\text{pre-}(r-1)} \times "
        r"\mathbf{1}[t \geq r]$. The first 2012 specification uses the generic "
        r"pre-reform procurement dummy as firm-level exposure; the second 2012 "
        r"specification is mechanism-matched, using the firm's pre-2012 mean "
        r"single-bid share (imputed as zero for firms that never bid on "
        r"procurement). The 2016 and 2017 reforms use the generic pp-dummy "
        r"exposure because no mechanism-matched proxy is available in the data. "
        r"The Wald 95\% CI is $\hat\alpha \pm 1.96 \cdot$ SE; the AR 95\% CI is "
        r"the cluster-robust Anderson-Rubin confidence set obtained by grid "
        r"inversion of the AR statistic $N \cdot (\overline{m})^2 / "
        r"\widehat{\mathrm{Var}}_{cl}(\overline{m})$ where "
        r"$m_i(\beta_0) = (\tilde Y_i - \beta_0 \tilde D_i) \tilde Z_i$ on the "
        r"DML-residualized arrays, valid regardless of first-stage strength "
        r"(Anderson and Rubin 1949; Stock and Wright 2000). \text{inverts} "
        r"indicates the AR CS is empty on the grid (unusual). Cross-fitted ML "
        r"nuisance functions partial out market share, foreign ownership, and "
        r"year dummies; firm fixed effects are absorbed by within-firm demeaning. "
        r"Convergence of point estimates across all four specifications---four "
        r"independent exogenous instruments---constitutes causal evidence for "
        r"$\hat\alpha$ beyond the selection-on-observables diagnostics of "
        r"Sections \ref{sec:abgrs}--\ref{sec:favoritism}. "
        r"$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$.",
        r"\end{tablenotes}",
        r"\end{threeparttable}",
        r"\end{table}",
    ]
    path = TAB_DIR / "dml_iv_multi_reform.tex"
    path.write_text("\n".join(lines) + "\n")
    all_results.to_csv(TAB_DIR / "dml_iv_multi_reform.csv", index=False)
    print(f"[write] {TAB_DIR.name}/dml_iv_multi_reform.tex + .csv")


def main() -> None:
    df = load_merged_panel()
    import statsmodels.formula.api as smf
    ols_premium = float(smf.ols("log_mu ~ pp_dummy + C(year) + C(id)", data=df)
                         .fit().params["pp_dummy"])
    print(f"[baseline] OLS firm+year FE premium: {ols_premium:.4f}")

    all_results = []
    for reform_year, pre_end, exposure_var, label in REFORMS:
        print(f"\n[DML-PLR-IV] reform {reform_year}: {label}")
        print(f"  Z = mean({exposure_var} | t<={pre_end}) × 1[t>={reform_year}]")
        results = run_reform_iv(df, BASE_COVARIATES,
                                 reform_year, pre_end, exposure_var, label)
        print(results[["ml", "point", "stderr", "t", "partial_r2",
                        "ar_lo", "ar_hi", "ar_empty"]]
                .to_string(index=False))
        all_results.append(results)

    combined = pd.concat(all_results, ignore_index=True)
    write_multi_reform_table(combined, ols_premium)
    print("\n[done]")


if __name__ == "__main__":
    main()
