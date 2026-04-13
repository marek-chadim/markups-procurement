"""Parse the Stata log file produced by panel_treatment_sdid_honestdid.do and
write a LaTeX table with the Goodman-Bacon (2021) decomposition of the
absorbing-treatment TWFE into four canonical comparison categories.

The bacondecomp Stata package stores per-comparison matrices in e(), but
does NOT expose the 4-category aggregation (Always-treated vs timing,
Never-treated vs timing, Earlier-treated vs later-treated, Later-treated
vs earlier-treated) as a summary return value. The printed display output
is the canonical source for the aggregation.

This script extracts the rows between the BEGIN/END BACONDECOMP TABLE
markers in the Stata log, aggregates the Early_v_Late and Late_v_Early
rows by category, keeps the Always_v_timing and Never_v_timing summary
rows as-is, and emits LaTeX.

Input:  2_analysis/output/panel_treatment_sdid_honestdid.log
Output: 2_analysis/output/tables/panel_bacon.tex
"""

from __future__ import annotations

import re
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR.parent / "output"
LOG_FILE   = OUTPUT_DIR / "panel_treatment_sdid_honestdid.log"
TEX_FILE   = OUTPUT_DIR / "tables" / "panel_bacon.tex"

ROW_RE = re.compile(
    r"^\|\s*(?P<label>\S+)\s*\|\s*(?P<beta>-?\.?\d+\.?\d*)\s+(?P<weight>-?\.?\d+\.?\d*)\s*\|$"
)


def parse_bacon_rows(log_text: str) -> list[tuple[str, float, float]]:
    """Return a list of (label, beta, weight) tuples parsed from the
    bacondecomp display table between the BEGIN/END markers."""
    in_block = False
    rows = []
    for line in log_text.splitlines():
        if "BEGIN BACONDECOMP TABLE" in line:
            in_block = True
            continue
        if "END BACONDECOMP TABLE" in line:
            break
        if not in_block:
            continue
        m = ROW_RE.match(line.rstrip())
        if not m:
            continue
        label = m.group("label").strip()
        beta = float(m.group("beta"))
        weight = float(m.group("weight"))
        rows.append((label, beta, weight))
    return rows


def aggregate(rows: list[tuple[str, float, float]]) -> dict[str, dict[str, float]]:
    """Aggregate detailed rows into the 4 canonical categories.

    Returns a dict keyed by category label with sub-dict containing
    'weight' (total) and 'beta' (weighted-average)."""
    buckets = {
        "Early_v_Late": {"weight": 0.0, "wsum": 0.0, "count": 0},
        "Late_v_Early": {"weight": 0.0, "wsum": 0.0, "count": 0},
        "Always_v_timing": {"weight": 0.0, "wsum": 0.0, "count": 0},
        "Never_v_timing": {"weight": 0.0, "wsum": 0.0, "count": 0},
    }
    for label, beta, weight in rows:
        if label in buckets:
            buckets[label]["weight"] += weight
            buckets[label]["wsum"] += weight * beta
            buckets[label]["count"] += 1
    out = {}
    for cat, b in buckets.items():
        if b["weight"] > 0:
            out[cat] = {
                "weight": b["weight"],
                "beta": b["wsum"] / b["weight"],
                "wsum": b["wsum"],
                "count": b["count"],
            }
    return out


CATEGORY_LABELS = {
    "Always_v_timing": r"Always-treated vs.\ timing groups (``clean'')",
    "Never_v_timing":  r"Never-treated vs.\ timing groups (``clean'')",
    "Early_v_Late":    r"Earlier-treated vs.\ later-treated (``forbidden'')",
    "Late_v_Early":    r"Later-treated vs.\ earlier-treated (``forbidden'')",
}

ORDER = ["Always_v_timing", "Never_v_timing", "Early_v_Late", "Late_v_Early"]


def build_tex(agg: dict[str, dict[str, float]]) -> str:
    total_weight = sum(b["weight"] for b in agg.values())
    total_twfe = sum(b["wsum"] for b in agg.values())
    lines = [
        r"\begin{table}[htbp]\centering",
        r"\caption{Goodman-Bacon Decomposition of the Absorbing-Treatment TWFE}\label{tab:panel_bacon}",
        r"\begin{threeparttable}",
        r"\begin{tabular}{p{6cm}ccc}",
        r"\toprule",
        r"Comparison class & Weight & Avg.\ 2$\times$2 $\hat{\beta}$ & Contribution \\",
        r"\midrule",
    ]
    for key in ORDER:
        if key not in agg:
            continue
        b = agg[key]
        label = CATEGORY_LABELS[key]
        weight = b["weight"]
        beta = b["beta"]
        contribution = b["wsum"]
        lines.append(
            f"{label} & {weight:.4f} & {beta:+.4f} & {contribution:+.4f} \\\\"
        )
    lines += [
        r"\midrule",
        f"Total & {total_weight:.4f} & --- & {total_twfe:+.4f} \\\\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\footnotesize\raggedright",
        r"\item \textit{Notes:} Goodman-Bacon (2021) decomposition of the absorbing-treatment TWFE estimator on a balanced sub-panel of 45 firms (675 firm-years), via the \texttt{bacondecomp} Stata package. The decomposition expresses the TWFE estimate as a weighted average of 2$\times$2 DiD comparisons from four classes: timing cohorts vs.\ the always-treated group, timing cohorts vs.\ the never-treated group, earlier cohorts used as controls for later cohorts, and later cohorts used as controls for earlier cohorts. The last two are ``forbidden'' under heterogeneous effects because they use treated units as implicit controls. The contribution column $w \cdot \hat{\beta}$ sums to the total TWFE.",
        r"\end{tablenotes}",
        r"\end{threeparttable}",
        r"\end{table}",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    if not LOG_FILE.exists():
        raise SystemExit(f"Log file not found: {LOG_FILE}")
    rows = parse_bacon_rows(LOG_FILE.read_text())
    if not rows:
        raise SystemExit("No bacondecomp rows parsed from log")
    agg = aggregate(rows)

    print(f"Parsed {len(rows)} bacondecomp rows")
    print("Aggregated categories:")
    for key in ORDER:
        if key in agg:
            b = agg[key]
            print(
                f"  {key:<20s} count={b['count']:3d}  "
                f"weight={b['weight']:.4f}  beta={b['beta']:+.4f}  "
                f"contribution={b['wsum']:+.4f}"
            )
    total_weight = sum(b["weight"] for b in agg.values())
    total_twfe = sum(b["wsum"] for b in agg.values())
    print(f"  {'TOTAL':<20s} weight={total_weight:.4f}  twfe={total_twfe:+.4f}")

    tex = build_tex(agg)
    TEX_FILE.parent.mkdir(parents=True, exist_ok=True)
    TEX_FILE.write_text(tex)
    print(f"Saved: {TEX_FILE.name}")


if __name__ == "__main__":
    main()
