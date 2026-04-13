"""Parse the Stata log file produced by panel_treatment_sdid_honestdid.do and
write a LaTeX table with the HonestDiD robust confidence intervals.

The honestdid Stata package does not expose its robust CI matrix via r()
return values, so the printed display output is the canonical source.
This script extracts the pipe-delimited table between the BEGIN/END HONESTDID
TABLE markers and emits LaTeX.

Input:  2_analysis/output/panel_treatment_sdid_honestdid.log
Output: 2_analysis/output/tables/panel_honestdid.tex
"""

from __future__ import annotations

import re
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR.parent / "output"
LOG_FILE   = OUTPUT_DIR / "panel_treatment_sdid_honestdid.log"
TEX_FILE   = OUTPUT_DIR / "tables" / "panel_honestdid.tex"


def parse_honestdid_rows(log_text: str) -> list[tuple[str, float, float]]:
    """Return a list of (label, lb, ub) tuples from the honestdid table.

    The first row after the header is (Original); subsequent rows have a
    numeric M value in the first column. Label is "Original" for the first
    row, formatted M otherwise.
    """
    in_block = False
    header_seen = False
    rows = []
    row_re = re.compile(
        r"^\|\s*(?P<m>[^\s|]+)\s*\|\s*(?P<lb>-?\d+\.\d+)\s*\|\s*(?P<ub>-?\d+\.\d+)\s*\|(?P<tail>.*)$"
    )
    for line in log_text.splitlines():
        if "BEGIN HONESTDID TABLE" in line:
            in_block = True
            continue
        if "END HONESTDID TABLE" in line:
            break
        if not in_block:
            continue
        if not header_seen:
            if line.strip().startswith("|    M"):
                header_seen = True
            continue
        m = row_re.match(line)
        if not m:
            continue
        m_str = m.group("m").strip()
        lb = float(m.group("lb"))
        ub = float(m.group("ub"))
        if m_str == ".":
            label = "Original OLS CI"
        else:
            label = f"$\\bar{{M}} = {float(m_str):.2f}$"
        rows.append((label, lb, ub))
    return rows


def build_tex(rows: list[tuple[str, float, float]]) -> str:
    lines = [
        r"\begin{table}[htbp]\centering",
        r"\caption{HonestDiD Sensitivity to Parallel-Trend Violations}\label{tab:panel_honestdid}",
        r"\begin{threeparttable}",
        r"\begin{tabular}{lcc}",
        r"\toprule",
        r"Sensitivity bound & Robust 95\% CI lower & Robust 95\% CI upper \\",
        r"\midrule",
    ]
    for label, lb, ub in rows:
        lines.append(f"{label} & {lb:+.4f} & {ub:+.4f} \\\\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\footnotesize",
        r"\item \textit{Notes:} Honest confidence sets for the average post-treatment event-study coefficient under relative-magnitude ($\Delta^{RM}$) restrictions on parallel-trend violations following Rambachan and Roth (2023), implemented via the \texttt{honestdid} Stata package of C\'{a}ceres and Rambachan (2023). The event-study coefficients that feed \texttt{honestdid} come from a two-way fixed-effects regression of $\log \mu^A_{it}$ on relative-year dummies (window $\tau \in [-5, +5]$; $\tau = -1$ as reference), firm and year fixed effects, clustered at the firm level. $\bar{M}$ is the maximum permitted ratio of post-treatment parallel-trend violation to the largest observed pre-treatment violation; $\bar{M} = 0$ reproduces the original OLS confidence set. The breakdown value $\bar{M}^\star$ is the smallest $\bar{M}$ for which the robust CI includes zero.",
        r"\end{tablenotes}",
        r"\end{threeparttable}",
        r"\end{table}",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    if not LOG_FILE.exists():
        raise SystemExit(f"Log file not found: {LOG_FILE}")
    log_text = LOG_FILE.read_text()
    rows = parse_honestdid_rows(log_text)
    if not rows:
        raise SystemExit("No honestdid rows parsed from log")
    tex = build_tex(rows)
    TEX_FILE.parent.mkdir(parents=True, exist_ok=True)
    TEX_FILE.write_text(tex)
    print(f"Saved: {TEX_FILE.name}")
    print("Parsed rows:")
    for label, lb, ub in rows:
        print(f"  {label}: [{lb:+.4f}, {ub:+.4f}]")


if __name__ == "__main__":
    main()
