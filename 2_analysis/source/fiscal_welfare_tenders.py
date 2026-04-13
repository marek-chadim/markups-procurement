"""fiscal_welfare_tenders.py — §7 Fiscal Welfare Implications.

Computes two Kaldor-Hicks welfare benchmarks for Czech public procurement
2003-2022 using the Datlab master tender register:

  1. Engineer-estimate benchmark (narrow):
     savings = Σ max(0, bid_final_price − lot_estimated_price)

  2. Marginal-cost benchmark (structural):
     savings = (μ − 1) / μ × Σ bid_final_price
     where μ = 1.14 is the paper's firm-level procurement markup premium.

The gap between the two benchmarks = baseline industry markup pre-priced
into engineer estimates; only the engineer-estimate benchmark is recoverable
through procurement-competition reforms alone.

Outputs:
    output/tables/fiscal_welfare_aggregates.tex       LaTeX table for paper
    output/tables/fiscal_welfare_aggregates.csv       raw numbers
    output/tables/fiscal_welfare_timeseries.csv       yearly series for figure
    output/figures/fig_fiscal_welfare_timeseries.pdf  dual-axis timeseries plot

Reference: Hendren and Sprung-Keyser 2020 QJE; Finkelstein and Hendren 2020 JEP.
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from style_markups import (
    apply_markups_style, MARKUPS_BLUE, MARKUPS_RED, MARKUPS_GREY,
)

# ---- Paths -----------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
TENDERS_CSV = SCRIPT_DIR.parents[1] / "1_data" / "input" / "datlab" / "master_tender_analytics.csv"
OUTPUT_DIR = SCRIPT_DIR.parent / "output"
FIG_DIR = OUTPUT_DIR / "figures"
TAB_DIR = OUTPUT_DIR / "tables"
for d in (FIG_DIR, TAB_DIR):
    d.mkdir(parents=True, exist_ok=True)

# ---- Parameters ------------------------------------------------------------
MARKUP_PREMIUM = 0.14      # paper's headline firm-level procurement markup premium
REL_PRICE_LOWER = 0.01     # winsorization lower bound (drop <1% of estimate)
REL_PRICE_UPPER = 10.0     # winsorization upper bound (drop >10× estimate)
BILLIONS = 1e9             # CZK → bn CZK display scale


def load_tenders() -> pd.DataFrame:
    """Load the Datlab master tender register."""
    print(f"[load] {TENDERS_CSV.name}")
    df = pd.read_csv(TENDERS_CSV, low_memory=False)
    print(f"  rows = {len(df):,}; years = {int(df['year'].min())}–{int(df['year'].max())}")
    return df


def build_subset(df: pd.DataFrame) -> pd.DataFrame:
    """Rows with both estimate and final price, winsorized on Rel_Price."""
    mask = (df["lot_estimated_price"].notna()
            & df["bid_final_price"].notna()
            & (df["lot_estimated_price"] > 0)
            & (df["bid_final_price"] > 0))
    sub = df.loc[mask].copy()
    sub["rel_price"] = sub["bid_final_price"] / sub["lot_estimated_price"]
    sub = sub[(sub["rel_price"] > REL_PRICE_LOWER) & (sub["rel_price"] < REL_PRICE_UPPER)].copy()
    print(f"[subset] both fields populated and winsorized: {len(sub):,} rows "
          f"({len(sub)/len(df)*100:.1f}% of universe)")
    return sub


def aggregates(df: pd.DataFrame, sub: pd.DataFrame) -> dict:
    """Compute the six headline welfare aggregates."""
    total_final_all = float(df["bid_final_price"].sum()) / BILLIONS  # full-universe spend
    n_years = int(df["year"].max() - df["year"].min() + 1)

    total_final_sub = float(sub["bid_final_price"].sum()) / BILLIONS
    total_est_sub = float(sub["lot_estimated_price"].sum()) / BILLIONS
    aggregate_ratio_sub = total_final_sub / total_est_sub

    # Engineer-estimate benchmark: overspend on Rel_Price > 1 tenders
    over = sub[sub["rel_price"] > 1.0]
    overspend_bn = float((over["bid_final_price"] - over["lot_estimated_price"]).sum()) / BILLIONS
    overspend_pct_of_subset = overspend_bn / total_final_sub * 100

    # Markup benchmark: structural savings on full universe
    markup_savings_bn = total_final_all * (MARKUP_PREMIUM / (1 + MARKUP_PREMIUM))
    markup_savings_pct = MARKUP_PREMIUM / (1 + MARKUP_PREMIUM) * 100

    gap_pp = markup_savings_pct - overspend_pct_of_subset

    return dict(
        n_years=n_years,
        n_universe=len(df),
        n_subset=len(sub),
        subset_share_pct=len(sub) / len(df) * 100,
        total_final_all_bn=total_final_all,
        annual_avg_bn=total_final_all / n_years,
        total_final_sub_bn=total_final_sub,
        total_est_sub_bn=total_est_sub,
        aggregate_ratio_sub=aggregate_ratio_sub,
        n_over=len(over),
        overspend_bn=overspend_bn,
        overspend_annual_bn=overspend_bn / n_years,
        overspend_pct_of_subset=overspend_pct_of_subset,
        markup_savings_bn=markup_savings_bn,
        markup_savings_annual_bn=markup_savings_bn / n_years,
        markup_savings_pct=markup_savings_pct,
        gap_pp=gap_pp,
    )


def write_aggregates_csv(a: dict) -> None:
    df = pd.DataFrame([{k: v for k, v in a.items()}])
    df.to_csv(TAB_DIR / "fiscal_welfare_aggregates.csv", index=False)
    print(f"[write] {TAB_DIR.name}/fiscal_welfare_aggregates.csv")


def write_aggregates_tex(a: dict) -> None:
    """6-row LaTeX table for the paper, formatted as a simple booktabs block."""
    lines = [
        r"\begin{table}[htbp]\centering",
        r"\caption{Fiscal Welfare Benchmarks: Czech Public Procurement, 2003--2022}",
        r"\label{tab:fiscal_welfare}",
        r"\begin{threeparttable}",
        r"\begin{tabular}{lrr}",
        r"\toprule",
        r"& Aggregate (bn CZK) & Share \\",
        r"\midrule",
        rf"Total procurement spending (universe) & "
        rf"{a['total_final_all_bn']:,.0f} & "
        rf"100.0\% \\",
        rf"\quad Annual average ({a['n_years']} years) & "
        rf"{a['annual_avg_bn']:,.1f} & "
        rf"\\",
        rf"Subset with engineer estimate ({a['n_subset']:,} contracts) & "
        rf"{a['total_final_sub_bn']:,.0f} & "
        rf"{a['subset_share_pct']:.1f}\% of universe \\",
        rf"\quad Aggregate final / engineer estimate & "
        rf"{a['aggregate_ratio_sub']:.4f} & "
        rf"\\",
        r"\midrule",
        rf"\textbf{{Benchmark A}} (engineer estimate): overspend & "
        rf"\textbf{{{a['overspend_bn']:,.0f}}} & "
        rf"\textbf{{{a['overspend_pct_of_subset']:.2f}\%}} of subset \\",
        rf"\quad Annual average & "
        rf"{a['overspend_annual_bn']:,.1f} & \\",
        rf"\textbf{{Benchmark B}} (marginal cost, $\mu = 1.14$) & "
        rf"\textbf{{{a['markup_savings_bn']:,.0f}}} & "
        rf"\textbf{{{a['markup_savings_pct']:.2f}\%}} of universe \\",
        rf"\quad Annual average & "
        rf"{a['markup_savings_annual_bn']:,.1f} & \\",
        r"\midrule",
        rf"Gap (B $-$ A, baseline industry markup) & "
        rf"{a['markup_savings_bn'] - a['overspend_bn'] * (a['total_final_all_bn'] / a['total_final_sub_bn']):,.0f} & "
        rf"{a['gap_pp']:.2f} pp \\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\footnotesize",
        r"\item \emph{Notes:} Benchmark A caps each contract at its engineer "
        r"estimate (\texttt{lot\_estimated\_price} in the Datlab register) and "
        r"sums the excess on the subset of tenders with $\text{bid} > \text{estimate}$. "
        r"Benchmark B prices firms at marginal cost using the paper's 14\% "
        r"firm-level procurement markup premium; savings $= \mu/(1+\mu) \times \Sigma \text{bid}$. "
        r"The gap $= \mu/(1+\mu) - $ subset overspend share in percentage points, "
        r"representing the baseline industry markup already priced into engineer "
        r"estimates, which is not recoverable through procurement-competition reforms.",
        r"\end{tablenotes}",
        r"\end{threeparttable}",
        r"\end{table}",
    ]
    path = TAB_DIR / "fiscal_welfare_aggregates.tex"
    path.write_text("\n".join(lines) + "\n")
    print(f"[write] {TAB_DIR.name}/fiscal_welfare_aggregates.tex")


def yearly_series(sub: pd.DataFrame) -> pd.DataFrame:
    """Year-level aggregate ratio and overspend for the time-series figure."""
    g = sub.groupby("year").agg(
        n=("rel_price", "size"),
        rel_price_med=("rel_price", "median"),
        est_sum=("lot_estimated_price", "sum"),
        final_sum=("bid_final_price", "sum"),
    ).reset_index()
    g["aggregate_ratio"] = g["final_sum"] / g["est_sum"]
    g["overspend_bn"] = np.maximum(0.0, g["final_sum"] - g["est_sum"]) / BILLIONS
    # Use the exact tender-level overspend aggregate, not the netted year-level
    over = sub[sub["rel_price"] > 1.0].copy()
    over["excess"] = over["bid_final_price"] - over["lot_estimated_price"]
    yearly_over = over.groupby("year")["excess"].sum().reset_index()
    yearly_over["excess_bn"] = yearly_over["excess"] / BILLIONS
    g = g.merge(yearly_over[["year", "excess_bn"]], on="year", how="left")
    g["excess_bn"] = g["excess_bn"].fillna(0.0)
    g = g[g["year"].between(2006, 2022)].reset_index(drop=True)
    g.to_csv(TAB_DIR / "fiscal_welfare_timeseries.csv", index=False)
    print(f"[write] {TAB_DIR.name}/fiscal_welfare_timeseries.csv ({len(g)} years)")
    return g


def plot_timeseries(yearly: pd.DataFrame) -> None:
    """Dual-axis: bars of annual overspend (left) + line of median Rel_Price (right)."""
    apply_markups_style()
    fig, ax1 = plt.subplots(figsize=(7.0, 4.2))

    years = yearly["year"].to_numpy()

    # Bars: annual overspend (bn CZK)
    ax1.bar(years, yearly["excess_bn"], color=MARKUPS_RED, alpha=0.75,
            label="Annual overspending (bn CZK)")
    ax1.set_xlabel("Year")
    ax1.set_ylabel("Annual overspending (bn CZK)", color=MARKUPS_RED)
    ax1.tick_params(axis="y", labelcolor=MARKUPS_RED)
    ax1.set_xticks(np.arange(2006, 2023, 2))

    # Line: median Rel_Price
    ax2 = ax1.twinx()
    ax2.plot(years, yearly["rel_price_med"], color=MARKUPS_BLUE, marker="o",
             markersize=4, linewidth=1.5, label="Median Rel. Price (right axis)")
    ax2.axhline(1.0, color=MARKUPS_GREY, linestyle=":", linewidth=0.8)
    ax2.set_ylabel("Median bid / engineer estimate", color=MARKUPS_BLUE)
    ax2.tick_params(axis="y", labelcolor=MARKUPS_BLUE)
    ax2.set_ylim(0.75, 1.05)

    # Reform year annotations
    for ry, label in [(2012, "2012\nsingle-bid ban"), (2016, "2016\nMEAT")]:
        ax1.axvline(ry, color=MARKUPS_GREY, linestyle="--", linewidth=0.8, alpha=0.6)
        ax1.text(ry, ax1.get_ylim()[1] * 0.92, label, ha="center", va="top",
                 fontsize=8, color="#555555")

    fig.tight_layout()
    path = FIG_DIR / "fig_fiscal_welfare_timeseries.pdf"
    fig.savefig(path)
    plt.close(fig)
    print(f"[write] {FIG_DIR.name}/fig_fiscal_welfare_timeseries.pdf")


def main() -> None:
    df = load_tenders()
    sub = build_subset(df)
    agg = aggregates(df, sub)

    print("\n[aggregates]")
    for k, v in agg.items():
        if isinstance(v, float):
            print(f"  {k:32s}: {v:12,.4f}")
        else:
            print(f"  {k:32s}: {v:>12}")

    write_aggregates_csv(agg)
    write_aggregates_tex(agg)
    yearly = yearly_series(sub)
    plot_timeseries(yearly)

    print("\n[done]")


if __name__ == "__main__":
    main()
