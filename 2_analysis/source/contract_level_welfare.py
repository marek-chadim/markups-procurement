"""contract_level_welfare.py — contract-level Rel_Price × firm markup regression.

Companion to fiscal_welfare_tenders.py. Where that script produced aggregate
welfare benchmarks, this one links individual procurement contracts to the
winning firm's ACF-estimated markup and tests whether high-markup firms
systematically bid higher relative to engineer estimates.

Merge strategy:
    Datlab bidder_id (= Czech IČO) → firm id in paper_markups.dta
Both are 8-digit integers after dropping 2-character foreign bidders, so the
join is direct on (id, year) — no crosswalk file required.

Regression (pooled, contract-level observations):
    log(Rel_Price_k) = β · log(markup_A_it) + γ · pp_dummy_it + FE + ε_k

where k indexes contracts, (i, t) indexes the winning firm-year. Fixed effects:
    Spec 1: year
    Spec 2: year + NACE 2-digit
    Spec 3: year + firm (within-firm identification)

Standard errors are clustered by firm. Expected sign: β > 0 — high-markup
firms exhibit higher bid-to-estimate ratios on procurement contracts.

Outputs:
    output/tables/contract_level_welfare.tex    regression table for paper
    output/tables/contract_level_welfare.csv    raw coefficients
    output/data/contract_level_merged.csv       merged dataset for diagnostics
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

# ---- Paths -----------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
TENDERS_CSV = SCRIPT_DIR.parents[1] / "1_data" / "input" / "datlab" / "master_tender_analytics.csv"
MARKUPS_DTA = SCRIPT_DIR.parent / "output" / "data" / "paper_markups.dta"
OUTPUT_DIR = SCRIPT_DIR.parent / "output"
TAB_DIR = OUTPUT_DIR / "tables"
DAT_DIR = OUTPUT_DIR / "data"
for d in (TAB_DIR, DAT_DIR):
    d.mkdir(parents=True, exist_ok=True)

# ---- Parameters ------------------------------------------------------------
REL_PRICE_LOWER = 0.2     # match favoritism_decomposition.py trimming
REL_PRICE_UPPER = 5.0
MARKUP_COL = "markup_A"   # paper's headline ACF markup spec


# ---------------------------------------------------------------------------
# Load and merge
# ---------------------------------------------------------------------------

def load_contracts() -> pd.DataFrame:
    """Read Datlab, compute Rel_Price, cast bidder_id to integer IČO."""
    print(f"[load] {TENDERS_CSV.name}")
    df = pd.read_csv(TENDERS_CSV, low_memory=False)
    print(f"  {len(df):,} contracts loaded")

    # Keep rows with both estimate and final price
    mask = (df["lot_estimated_price"].notna()
            & df["bid_final_price"].notna()
            & (df["lot_estimated_price"] > 0)
            & (df["bid_final_price"] > 0))
    df = df.loc[mask].copy()
    df["rel_price"] = df["bid_final_price"] / df["lot_estimated_price"]
    df = df[(df["rel_price"] >= REL_PRICE_LOWER)
            & (df["rel_price"] <= REL_PRICE_UPPER)].copy()

    # Cast bidder_id → int (drop 2-char foreign codes and non-numerics)
    df["bidder_id"] = df["bidder_id"].astype(str)
    df = df[df["bidder_id"].str.len() > 2].copy()
    df["id"] = pd.to_numeric(df["bidder_id"], errors="coerce")
    df = df.dropna(subset=["id", "year"]).copy()
    df["id"] = df["id"].astype(int)
    df["year"] = df["year"].astype(int)

    print(f"  {len(df):,} contracts with valid Rel_Price and integer IČO")
    return df


def load_markups() -> pd.DataFrame:
    """Read paper_markups.dta, return firm-year panel with markup_A."""
    print(f"[load] {MARKUPS_DTA.name}")
    dfm = pd.read_stata(MARKUPS_DTA)
    cols = ["id", "year", MARKUP_COL, "pp_dummy", "nace2"]
    dfm = dfm[cols].copy()
    dfm["id"] = dfm["id"].astype(int)
    dfm["year"] = dfm["year"].astype(int)
    dfm = dfm.dropna(subset=[MARKUP_COL])
    print(f"  {len(dfm):,} firm-year observations with non-null {MARKUP_COL}")
    return dfm


def merge(contracts: pd.DataFrame, markups: pd.DataFrame) -> pd.DataFrame:
    """Inner join contracts to firm-year markups on (id, year)."""
    merged = contracts.merge(markups, on=["id", "year"], how="inner",
                              validate="many_to_one")
    print(f"[merge] {len(merged):,} contracts matched to "
          f"{merged['id'].nunique():,} firms across "
          f"{merged['year'].nunique()} years")

    # De-duplicate exact-match records (same firm, year, bid, estimate).
    # Datlab contains duplicate rows for a handful of large infrastructure
    # contracts with near-identical titles; these would inflate the
    # concentration statistics and the top-N incidence table.
    pre = len(merged)
    merged = merged.drop_duplicates(
        subset=["id", "year", "bid_final_price", "lot_estimated_price"]
    ).copy()
    print(f"[dedup] dropped {pre - len(merged):,} exact-duplicate contracts "
          f"(now {len(merged):,})")

    # Drop non-construction (keep NACE 41, 42, 43)
    merged = merged[merged["nace2"].isin([41, 42, 43])].copy()
    print(f"  {len(merged):,} after NACE F (41/42/43) filter")

    # Log transforms
    merged["log_rel_price"] = np.log(merged["rel_price"])
    merged["log_markup"] = np.log(merged[MARKUP_COL])

    merged.to_csv(DAT_DIR / "contract_level_merged.csv", index=False)
    print(f"[write] {DAT_DIR.name}/contract_level_merged.csv")
    return merged


# ---------------------------------------------------------------------------
# Regression specs
# ---------------------------------------------------------------------------

def _fit_ols_with_fe(df: pd.DataFrame, fe_cols: list[str], label: str) -> dict:
    """Pooled OLS with dummy-variable FE and firm-clustered SE."""
    X_parts = [df[["log_markup", "pp_dummy"]].astype(float)]
    for fe in fe_cols:
        dummies = pd.get_dummies(df[fe], prefix=fe, drop_first=True).astype(float)
        X_parts.append(dummies)
    X = pd.concat(X_parts, axis=1)
    X = sm.add_constant(X)
    y = df["log_rel_price"].astype(float)

    model = sm.OLS(y, X)
    res = model.fit(cov_type="cluster", cov_kwds={"groups": df["id"]})

    return dict(
        spec=label,
        n=int(res.nobs),
        n_firms=int(df["id"].nunique()),
        beta_markup=float(res.params["log_markup"]),
        se_markup=float(res.bse["log_markup"]),
        t_markup=float(res.tvalues["log_markup"]),
        p_markup=float(res.pvalues["log_markup"]),
        beta_pp=float(res.params["pp_dummy"]),
        se_pp=float(res.bse["pp_dummy"]),
        r2=float(res.rsquared),
    )


def run_regressions(df: pd.DataFrame) -> pd.DataFrame:
    """Three progressively stricter fixed-effects specifications."""
    print("\n[regressions]")
    results = []

    # Spec 1: year FE only
    print("  Spec 1: year FE, clustered SE")
    results.append(_fit_ols_with_fe(df, ["year"], "Year FE"))

    # Spec 2: year + NACE FE
    print("  Spec 2: year + NACE FE, clustered SE")
    results.append(_fit_ols_with_fe(df, ["year", "nace2"], "Year + NACE FE"))

    # Spec 3: year + firm FE (within-firm identification)
    # With 1000+ firms this is expensive; demean log_rel_price by firm instead
    print("  Spec 3: year + firm FE via within demeaning, clustered SE")
    results.append(_fit_firm_fe(df))

    return pd.DataFrame(results)


def _fit_firm_fe(df: pd.DataFrame) -> dict:
    """Firm FE via within-firm demeaning + year dummies + clustered SE."""
    d = df.copy()
    # Demean y and x by firm
    for c in ["log_rel_price", "log_markup", "pp_dummy"]:
        d[f"{c}_dm"] = d[c] - d.groupby("id")[c].transform("mean")

    X_parts = [d[["log_markup_dm", "pp_dummy_dm"]].astype(float)]
    dummies = pd.get_dummies(d["year"], prefix="year", drop_first=True).astype(float)
    X_parts.append(dummies)
    X = pd.concat(X_parts, axis=1)

    y = d["log_rel_price_dm"].astype(float)
    res = sm.OLS(y, X).fit(cov_type="cluster", cov_kwds={"groups": d["id"]})

    return dict(
        spec="Year + Firm FE",
        n=int(res.nobs),
        n_firms=int(d["id"].nunique()),
        beta_markup=float(res.params["log_markup_dm"]),
        se_markup=float(res.bse["log_markup_dm"]),
        t_markup=float(res.tvalues["log_markup_dm"]),
        p_markup=float(res.pvalues["log_markup_dm"]),
        beta_pp=float(res.params["pp_dummy_dm"]),
        se_pp=float(res.bse["pp_dummy_dm"]),
        r2=float(res.rsquared),
    )


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _stars(p: float) -> str:
    if p < 0.01:
        return "***"
    if p < 0.05:
        return "**"
    if p < 0.10:
        return "*"
    return ""


def write_table(res: pd.DataFrame) -> None:
    res.to_csv(TAB_DIR / "contract_level_welfare.csv", index=False)

    lines = [
        r"\begin{table}[htbp]\centering",
        r"\caption{Contract-Level Welfare Regression: Firm Markup and Rel.\ Price}",
        r"\label{tab:contract_welfare}",
        r"\begin{threeparttable}",
        r"\begin{tabular}{lccc}",
        r"\toprule",
        r"& (1) & (2) & (3) \\",
        r"& Year FE & Year $+$ NACE FE & Year $+$ Firm FE \\",
        r"\midrule",
    ]

    def row_cells(values):
        return " & ".join(values)

    beta_cells = [f"{r['beta_markup']:.3f}{_stars(r['p_markup'])}" for _, r in res.iterrows()]
    se_cells = [f"({r['se_markup']:.3f})" for _, r in res.iterrows()]
    pp_cells = [f"{r['beta_pp']:.3f}" for _, r in res.iterrows()]
    pp_se_cells = [f"({r['se_pp']:.3f})" for _, r in res.iterrows()]
    r2_cells = [f"{r['r2']:.4f}" for _, r in res.iterrows()]
    n_cells = [f"{int(r['n']):,}" for _, r in res.iterrows()]
    firms_cells = [f"{int(r['n_firms']):,}" for _, r in res.iterrows()]

    lines += [
        rf"log(markup$_{{\text{{ACF}}}}$) & {row_cells(beta_cells)} \\",
        rf" & {row_cells(se_cells)} \\",
        rf"$\mathbb{{1}}$[$pp_{{it}} = 1$] & {row_cells(pp_cells)} \\",
        rf" & {row_cells(pp_se_cells)} \\",
        r"\midrule",
        rf"Contracts ($N$) & {row_cells(n_cells)} \\",
        rf"Firms & {row_cells(firms_cells)} \\",
        rf"$R^2$ (within for (3)) & {row_cells(r2_cells)} \\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\footnotesize",
        r"\item \emph{Notes:} Dependent variable is $\log(\text{Rel.\ Price}_k) = "
        r"\log(\text{bid\_final\_price}_k / \text{lot\_estimated\_price}_k)$ at the "
        r"contract level, trimmed to Rel.\ Price $\in [0.2, 5]$. Firm-level markups "
        r"are the paper's ACF Specification~A (\texttt{markup\_A} in "
        r"\texttt{paper\_markups.dta}), matched to contracts via Datlab "
        r"\texttt{bidder\_id} $\to$ Czech IČO $\to$ MagnusWeb firm \texttt{id}. "
        r"Sample restricted to CZ-NACE F (construction, divisions 41/42/43). "
        r"Column~(3) absorbs firm fixed effects via within-firm demeaning. "
        r"Standard errors in parentheses are clustered by firm. "
        r"$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$.",
        r"\end{tablenotes}",
        r"\end{threeparttable}",
        r"\end{table}",
    ]
    path = TAB_DIR / "contract_level_welfare.tex"
    path.write_text("\n".join(lines) + "\n")
    print(f"[write] {TAB_DIR.name}/contract_level_welfare.tex")


def _truncate(s: str, n: int = 65) -> str:
    """Truncate a string for table display; strip parenthetical procurer
    and replace Czech quote-marks that break OT1 LaTeX encoding."""
    s = str(s)
    # Drop "(zadavatel: ...)" tail for cleaner titles
    if "(zadavatel:" in s:
        s = s.split("(zadavatel:")[0].strip()
    # Replace Czech low-9 quote (U+201E) and close quote (U+201C/201D)
    # with ASCII double-quote to avoid \quotedblbase encoding errors.
    for bad in ("\u201e", "\u201c", "\u201d", "\u2018", "\u2019"):
        s = s.replace(bad, '"' if bad in ("\u201e", "\u201c", "\u201d") else "'")
    # Drop any leftover double-quotes around the title (noise)
    s = s.strip('"')
    if len(s) > n:
        return s[: n - 3] + "..."
    return s


def heterogeneity_analysis(df: pd.DataFrame) -> None:
    """Within-firm pass-through by NACE sub-industry and by reform era.

    Runs column-3 spec (year + firm FE, clustered SE) on 5 subsets:
        NACE 41 (buildings), NACE 42 (civil engineering), NACE 43 (specialized)
        pre-2012 single-bidding ban, post-2012 single-bidding ban.

    Reports the β on log(markup_A) for each subset, which isolates the
    within-firm pass-through of markup fluctuations to contract-level Rel.
    Price within that cell.
    """
    print("\n[heterogeneity]")
    rows = []

    for nace, label in [(41, "NACE 41 buildings"),
                         (42, "NACE 42 civil eng."),
                         (43, "NACE 43 specialized")]:
        sub = df[df["nace2"] == nace].copy()
        if sub["id"].nunique() < 5:
            print(f"  {label}: skipped (too few firms)")
            continue
        r = _fit_firm_fe(sub)
        r["subset"] = label
        rows.append(r)
        print(f"  {label}: β = {r['beta_markup']:+.3f} "
              f"(SE {r['se_markup']:.3f}, p = {r['p_markup']:.4f}, "
              f"N = {r['n']}, firms = {r['n_firms']})")

    for year_range, label in [((2006, 2011), "Pre-2012 reform"),
                                ((2012, 2021), "Post-2012 reform")]:
        sub = df[(df["year"] >= year_range[0]) & (df["year"] <= year_range[1])].copy()
        if sub["id"].nunique() < 5:
            print(f"  {label}: skipped (too few firms)")
            continue
        r = _fit_firm_fe(sub)
        r["subset"] = label
        rows.append(r)
        print(f"  {label}: β = {r['beta_markup']:+.3f} "
              f"(SE {r['se_markup']:.3f}, p = {r['p_markup']:.4f}, "
              f"N = {r['n']}, firms = {r['n_firms']})")

    het = pd.DataFrame(rows)
    het.to_csv(TAB_DIR / "contract_welfare_heterogeneity.csv", index=False)
    print(f"[write] {TAB_DIR.name}/contract_welfare_heterogeneity.csv")

    # --- LaTeX table ---
    lines = [
        r"\begin{table}[htbp]\centering",
        r"\caption{Within-Firm Pass-Through Heterogeneity: NACE Sub-Industries and the 2012 Reform}",
        r"\label{tab:contract_welfare_het}",
        r"\begin{threeparttable}",
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"& $\hat\beta_{\log \mu}$ & SE & Contracts & Firms \\",
        r"\midrule",
        r"\emph{Panel A: By NACE sub-industry} & & & & \\",
    ]
    for label in ["NACE 41 buildings", "NACE 42 civil eng.", "NACE 43 specialized"]:
        r = het[het["subset"] == label]
        if len(r) == 0:
            continue
        r = r.iloc[0]
        lines.append(
            rf"\quad {label} & {r['beta_markup']:+.3f}{_stars(r['p_markup'])} & "
            rf"({r['se_markup']:.3f}) & {int(r['n']):,} & {int(r['n_firms']):,} \\"
        )
    lines += [
        r"\midrule",
        r"\emph{Panel B: Before vs.\ after the 2012 single-bidding ban} & & & & \\",
    ]
    for label in ["Pre-2012 reform", "Post-2012 reform"]:
        r = het[het["subset"] == label]
        if len(r) == 0:
            continue
        r = r.iloc[0]
        lines.append(
            rf"\quad {label} & {r['beta_markup']:+.3f}{_stars(r['p_markup'])} & "
            rf"({r['se_markup']:.3f}) & {int(r['n']):,} & {int(r['n_firms']):,} \\"
        )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\footnotesize",
        r"\item \emph{Notes:} Each row reports the coefficient on $\log(\mu_{\text{ACF}})$ "
        r"from the Table~\ref{tab:contract_welfare} column-(3) specification "
        r"(dependent variable $\log(\text{Rel.\ Price}_k)$; year fixed effects; "
        r"firm fixed effects via within-firm demeaning; standard errors clustered "
        r"by firm) re-estimated on the indicated subset. Panel A partitions "
        r"the merged sample by NACE 2-digit sub-industry. Panel B partitions "
        r"by the 2012 Act~55/2012 single-bidding ban (pre-2012 = 2006--2011; "
        r"post-2012 = 2012--2021). $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$.",
        r"\end{tablenotes}",
        r"\end{threeparttable}",
        r"\end{table}",
    ]
    path = TAB_DIR / "contract_welfare_heterogeneity.tex"
    path.write_text("\n".join(lines) + "\n")
    print(f"[write] {TAB_DIR.name}/contract_welfare_heterogeneity.tex")


def incidence_analysis(df: pd.DataFrame) -> None:
    """Concentration stats + top-20 incidence table for Appendix B.7."""
    d = df.copy()
    d["transfer_czk"] = d["bid_final_price"] * (d[MARKUP_COL] - 1) / d[MARKUP_COL]
    d["bid_mn"] = d["bid_final_price"] / 1e6
    d["est_mn"] = d["lot_estimated_price"] / 1e6
    d["transfer_mn"] = d["transfer_czk"] / 1e6
    d["overspend_mn"] = d["bid_mn"] - d["est_mn"]

    total_bid_bn = d["bid_mn"].sum() / 1e3
    total_transfer_bn = d["transfer_mn"].sum() / 1e3
    print(f"\n[incidence] merged sample: {len(d):,} contracts, "
          f"{total_bid_bn:,.1f} bn CZK bid, {total_transfer_bn:,.1f} bn CZK transfer")

    # --- Concentration statistics ---
    d_sorted = d.sort_values("transfer_czk", ascending=False).reset_index(drop=True)
    cum_share = d_sorted["transfer_czk"].cumsum() / d_sorted["transfer_czk"].sum()
    conc = {}
    for pct in [0.01, 0.05, 0.10, 0.25, 0.50]:
        k = int(np.ceil(len(d_sorted) * pct))
        conc[pct] = dict(k=k, share=float(cum_share.iloc[k - 1]))

    print("[concentration]")
    for pct, v in conc.items():
        print(f"  Top {pct*100:4.1f}% of contracts ({v['k']:>5,}): "
              f"{v['share']*100:5.1f}% of total transfer")

    pd.DataFrame([
        dict(percentile=pct, n_contracts=v["k"], cumulative_share=v["share"])
        for pct, v in conc.items()
    ]).to_csv(TAB_DIR / "welfare_concentration.csv", index=False)
    print(f"[write] {TAB_DIR.name}/welfare_concentration.csv")

    # --- NACE sub-industry breakdown ---
    nace = d.groupby("nace2").agg(
        n=("bid_mn", "size"),
        bid_bn=("bid_mn", lambda x: float(x.sum()) / 1e3),
        transfer_bn=("transfer_mn", lambda x: float(x.sum()) / 1e3),
        mean_mu=(MARKUP_COL, "mean"),
        med_rp=("rel_price", "median"),
    ).reset_index()
    nace["transfer_share_pct"] = nace["transfer_bn"] / nace["bid_bn"] * 100
    nace.to_csv(TAB_DIR / "welfare_by_nace.csv", index=False)
    print(f"[write] {TAB_DIR.name}/welfare_by_nace.csv")
    print(nace.to_string(index=False))

    # --- Top 20 contracts by implied fiscal transfer ---
    top20 = d_sorted.head(20).copy()
    top20["title_short"] = top20["title"].apply(_truncate)

    lines = [
        r"\begin{table}[htbp]\centering",
        r"\caption{Top 20 Construction Contracts by Implied Fiscal Transfer, 2006--2021}",
        r"\label{tab:incidence_top20}",
        r"\begin{threeparttable}",
        r"\footnotesize",
        r"\begin{tabular}{clllrrr}",
        r"\toprule",
        r"\# & Year & Firm & Project title & Bid (mn) & $\mu_{\text{ACF}}$ & Transfer (mn) \\",
        r"\midrule",
    ]
    for i, r in top20.reset_index(drop=True).iterrows():
        firm = str(r["bidder_name"])[:30].replace("&", r"\&").replace("_", r"\_")
        title = str(r["title_short"]).replace("&", r"\&").replace("_", r"\_").replace("#", r"\#")
        lines.append(
            rf"{i+1} & {int(r['year'])} & {firm} & {title} & "
            rf"{r['bid_mn']:,.0f} & {r[MARKUP_COL]:.2f} & {r['transfer_mn']:,.0f} \\"
        )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\footnotesize",
        r"\item \emph{Notes:} Contracts ranked by implied fiscal transfer, defined "
        r"as $\text{bid}_k \times (\mu_{it} - 1) / \mu_{it}$ where $\mu_{it}$ is "
        r"the winning firm's ACF Specification~A markup in the year of contract "
        r"signature. Sample is the $8{,}860$-contract merged panel described in "
        r"Table~\ref{tab:contract_welfare}, restricted to CZ-NACE~F construction. "
        r"All firm names and contract titles are from the public Datlab "
        r"procurement register. \emph{The ``transfer'' column uses absolute markups "
        r"over marginal cost and therefore includes normal returns to capital; "
        r"the policy-relevant counterfactual is Benchmark~B of "
        r"Table~\ref{tab:fiscal_welfare}, which applies the $14\%$ procurement-specific "
        r"premium.} The top 20 contracts are concentrated in large-scale civil-"
        r"engineering projects (NACE~42) for the Czech railway and motorway "
        r"directorates.",
        r"\end{tablenotes}",
        r"\end{threeparttable}",
        r"\end{table}",
    ]
    path = TAB_DIR / "incidence_top_contracts.tex"
    path.write_text("\n".join(lines) + "\n")
    print(f"[write] {TAB_DIR.name}/incidence_top_contracts.tex")

    # Also write raw top-20 CSV for inspection
    top20[["year", "nace2", "bidder_name", "title", "bid_mn", "est_mn",
            "rel_price", MARKUP_COL, "transfer_mn"]].to_csv(
        TAB_DIR / "incidence_top_contracts.csv", index=False
    )


def main() -> None:
    contracts = load_contracts()
    markups = load_markups()
    merged = merge(contracts, markups)

    print(f"\n[sample] log(Rel_Price) mean = {merged['log_rel_price'].mean():.4f}, "
          f"SD = {merged['log_rel_price'].std():.4f}")
    print(f"[sample] log(markup_A)   mean = {merged['log_markup'].mean():.4f}, "
          f"SD = {merged['log_markup'].std():.4f}")
    print(f"[sample] pp_dummy share        = {merged['pp_dummy'].mean():.4f}")

    res = run_regressions(merged)

    print("\n[results]")
    with pd.option_context("display.float_format", "{:10.4f}".format,
                            "display.width", 120):
        print(res.to_string(index=False))

    write_table(res)
    heterogeneity_analysis(merged)
    incidence_analysis(merged)
    print("\n[done]")


if __name__ == "__main__":
    main()
