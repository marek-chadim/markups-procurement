# Benchmark Synthesis Report — Four-Paper Action List

**Generated:** 2026-04-13 (session_apr13b)
**Scope:** Standalone deliverable — no repository edits from the originating session
**Canonical paper:** `markups-procurement/4_paper/source/markups_procurement.tex` (70 pages)
**Related deliverables:** `/Users/marek/.claude/plans/cryptic-riding-haven.md` (Beer-paper critical evaluation), `jmp_readiness_report.md` (prior 4-agent review, 77→92+)

---

## 1. Context

This report synthesizes a critical evaluation of `markups-procurement` against four benchmark papers the user flagged as potentially implementable. All four are already partially engaged by the paper — the question was whether engagement is load-bearing enough or whether additional analyses should be added.

**Papers evaluated:**

| Benchmark | Full citation | Role in the paper | Current engagement |
|-|-|-|-|
| **Beer paper** | De Loecker & Scott (2025), "Markup Estimation using Production and Demand Data: An Application to the US Brewing Industry," NYU/KU Leuven WP | Closest methodological reference for production-approach markups (same author as DLW 2012) | Orphan bibitem at line 1020 — never `\cite`d in main text |
| **KLMS** | Kroft, Luo, Mogstad & Setzler (2025), "Imperfect Competition and Rents in Labor and Product Markets," AER 115(9) 2926–2969 | Closest substantive reference for construction-industry market power | 6 `\cite{KLMS2025}` anchors (lines 68, 177, 817, 848); full `klms_analysis.py` implementation DISABLED at `make.sh:90` |
| **BMY** | Benkard, Miller & Yurukoglu (2026), "The Rise of Market Power: Comment," WP | Canonical reference for markup-level specification sensitivity (motivates the paper's core thesis) | 2 citations (lines 60, 68) as motivation; BMY decomposition already noted in memory; no code replication |
| **NKP** | Norris Keiller, de Paula & Van Reenen (2024), "Production Function Estimation Using Subjective Expectations Data," NBER WP 32725 | Frontier alternative PF estimator (avoids proxy-variable monotonicity assumption) | **Not cited** — zero engagement currently |

**Key tension the action list addresses:** The paper invokes "near-Leontief technology" at §6.2 (line 476) to explain the ADL null, but Spec A is translog. Simultaneously, it relies on GNR non-identification ("premium identified from observed cost shares") while the Beer paper *chose Leontief* for the same reason and NKP *attacks* the GNR argument on optimization-error grounds. This is the sharpest internal tension and is what BMY-R1 and NKP-R2 are designed to resolve.

---

## 2. Locked decisions (from the scoping conversation)

These decisions were ratified via `AskUserQuestion` in the session that produced this synthesis. Any future implementation should respect them unless the user overrides.

| Decision | Choice | Implementation note |
|-|-|-|
| Scope of this pass | **Standalone deliverable** | No repository edits. Actions below are for a future implementation round. |
| §6.2 Leontief/translog tension | **Option 3b — defend translog in a paragraph** | Keep Spec A = translog. Add a ~200-word defense passage; **do not** promote Leontief to main spec (no code/table regeneration). |
| "Closest papers" triplet in intro | **Expand to four** | Keep KLMS/ADL/DGM; add Beer paper as a fourth explicit contrast. Do not swap DGM out. |

---

## 3. Action list (priority-ranked, all four benchmarks)

All actions are **narrative additions or pipeline re-enablement**. None of them require re-estimating the headline premium (0.138) or regenerating existing tables. New computation is limited to BMY-R1 (one new figure) and BMY-R2 (one new robustness row).

### 3.1 Beer paper (~4 hours, 5 actions)

**Beer-R1 — Cite the Beer paper in 3–5 places** (≈30 min)
- **Problem:** `\bibitem{DeLoeckerScott2025}` at line 1020 is an orphan — never `\cite`d in main text.
- **Fix:**
  - Line ~68 (contribution paragraph): expand "three closest papers" to **four**. Add passage: *"Unlike De Loecker & Scott (2025), who validate production-approach markups against a Miller-Weinberg demand system for US beer, I pursue within-production robustness because bid-level losing data are unavailable in the Czech procurement register."*
  - §5.1 nine-methods discussion (~line 170): cite Beer paper as the benchmark for cross-method markup-level agreement, contrast with the level-sensitivity finding.
  - §6.2 ADL null (~line 476): cite Beer paper's Leontief choice as independent evidence that fixed-proportions technology is the norm in inputs-dominated industries.
- **Files:** `markups-procurement/4_paper/source/markups_procurement.tex`
- **Verification:** `grep -c 'DeLoeckerScott2025' markups_procurement.tex` should return **3+** (currently 1, the bibitem alone).

**Beer-R2 — Reframe §7 pass-through as the procurement analog of Beer Figure 6** (≈1 hour)
- **Problem:** §7's β = 0.81 (within-firm contract-level pass-through) is functionally equivalent to Beer paper's 97% retail pass-through, but the paper does not make this mapping explicit.
- **Fix:** Rewrite the §7 intro paragraph to position the regression as "the procurement analog of Beer paper Figure 6 — both ask whether firm-level markups map into the prices buyers pay." Cite Beer paper at this point.
- **Files:** `markups-procurement/4_paper/source/markups_procurement.tex` §7 intro (around line 577).

**Beer-R3 — Defend translog against the near-Leontief argument (option 3b)** (≈1–2 hours)
- **Problem:** §6.2 argues the ADL null is due to "near-Leontief technology," but Spec A is translog. Referee-visible tension.
- **Fix:** Add a ~200-word passage in §3.4 or §6.2 arguing:
  1. Translog with firm-specific elasticities **nests the Leontief limit** as $(\beta_{vv}, \beta_{kv}) \to 0$.
  2. Under GNR non-identification, the premium projects onto the observed cost-share ratio regardless of which point in the translog parameter space the GMM selects.
  3. The CWDL Leontief robustness (§5.2, line 466) recovers premium = 0.155, within 1.7 pp of the translog 0.138, confirming the cost-share channel is load-bearing.
- **Files:** `markups-procurement/4_paper/source/markups_procurement.tex` §3.4 or §6.2.
- **Risk:** Narrative-only; no code change; no regeneration.

**Beer-R4 — Add a "why not demand-side triangulation" disclaimer** (≈30 min)
- **Problem:** Beer paper's core contribution is production + demand triangulation. Czech procurement data don't permit this. Paper should say so explicitly rather than leave it implicit.
- **Fix:** Add one paragraph (likely in §3 or §5) explaining: (a) Czech procurement records winning bids only — full bid vectors are unavailable pre-2017 and have ~50% coverage after, so a structural auction model is infeasible; (b) public buyers are a single demand side with no product differentiation in the BLP sense; (c) scoring-auction complement is a natural future extension once bid-level data are available.
- **Files:** `markups-procurement/4_paper/source/markups_procurement.tex` §3 or §5.

**Beer-R5 — Promote Czech declining-markup contrast to main text** (≈1 hour)
- **Problem:** Sales-weighted aggregate markup declined 2.39 (2007) → 2.14 (2021) in Czech construction, counter-trend to Beer's US rising markups. This fact currently lives only in memory notes (MEMORY.md Key Findings).
- **Fix:** Add one paragraph in §5 or §7 contrasting the Czech decline with the Beer-paper / DLEU rising-US trend. Reference the 2012 single-bid ban as a plausible accelerator.
- **Files:** `markups-procurement/4_paper/source/markups_procurement.tex` §5 or §7.

---

### 3.2 KLMS (~2 hours, 1 action)

**KLMS-R1 — Re-enable `klms_analysis.py` and add an appendix comparison table** (≈2 hours)
- **Current state:** `klms_analysis.py` (886 lines, 6 full methods: event study, θ via DiD/2SLS, 1−ε, β_L, double markdown, crowd-out). DISABLED in `make.sh:90`. Czech results already computed and documented in `memory/klms_study.md`:
  - (1 − ε) firm FE = **0.896** (KLMS US: 0.863) — price markup 1.116 vs 1.159
  - Double wage markdown = **0.720** (KLMS: 0.693) with calibrated $\theta = 0.245$
  - Double price markup = **1.390** (KLMS: 1.443)
  - Crowd-out: small/null in Czech (KLMS US: 27%)
- **Fix:**
  1. Uncomment `run_python klms_analysis.py ...` at `markups-procurement/2_analysis/make.sh:90`
  2. Add an appendix table "Double Market Power Comparison: Czech Construction vs KLMS US Benchmark" with 4 rows (1−ε, θ, double markdown, double markup) and 2 columns (Czech, KLMS US)
  3. **Footnote the θ = 0.245 calibration**: MagnusWeb employment is categorical, Czech θ is not identifiable from first differences; the calibration matches the KLMS US estimate. Orbis continuous employment exists at Yale WRDS (24MB, `orbis_panel_construction.dta`) and would unlock Czech θ identification in a future pass.
- **Files:**
  - `markups-procurement/2_analysis/make.sh:90` (re-enable)
  - `markups-procurement/4_paper/source/markups_procurement.tex` (new appendix table, likely in the existing methods robustness section)
  - `markups-procurement/4_paper/input/tables/klms_comparison.tex` (new table file, emitted by `klms_analysis.py`)
- **Verification:**
  - `bash 2_analysis/make.sh` should run `klms_analysis.py` without errors
  - `klms_comparison.tex` should exist in `4_paper/input/tables/`
  - Compiling the .tex should pull in the new table with no undefined references

**What NOT to do for KLMS:**
- Do not attempt θ estimation from MagnusWeb brackets (Δℓ ≈ 0; the script already detects this and falls back to calibration)
- Do not port the R/Python replication package at `references/replications/KLMS_AER_repkit/` — it's US IRS-data-specific and your Python reimplementation is already more integrated
- Do not attempt bid-level auction inversion — Datlab has winner-only bids pre-2017 and ~50% coverage after, insufficient for Propositions 3–4

---

### 3.3 BMY (~6 hours, 2 actions)

**BMY-R1 — Rolling-window PF coefficient stability figure** (≈4 hours) **[HIGHEST-ROI NEW ANALYSIS]**
- **Problem:** BMY's central finding is that $\hat\theta^{\text{cogs}}$ drifts substantially across 5-year rolling windows within the same industry. The paper's counter-claim — "levels vary but treatment effect is stable" — is precisely the anti-BMY argument, but currently the paper does **not show the PF-coefficient drift**. This is the most direct operationalization of the paper's thesis.
- **Fix:** Add a new script `markups-procurement/2_analysis/source/bmy_rolling_stability.py` that:
  1. For each of NACE 41/42/43, estimates Spec A (translog) on rolling 3-year windows: 2006–2008, 2009–2011, ..., 2019–2021.
  2. Plots $\hat\beta_v^{\text{linear}}$ with 95% CIs per window.
  3. Overlays the procurement premium estimated on the same window.
  4. Emits `2_analysis/output/figures/bmy_rolling_stability.pdf` and `2_analysis/output/tables/bmy_rolling_stability.tex`.
  5. Add a figure to the paper with caption: "PF Coefficients Vary Across Windows; Procurement Premium Does Not."
- **Implementation notes:**
  - `acf_estimator.py` already supports year-range subsetting via `years` argument; wrap it in a loop over windows
  - Keep window size = 3 (rather than BMY's 5) to maximize the number of windows given Chadim's 16-year panel
  - Use `style_markups.py` for colors/fonts
- **Files:**
  - `markups-procurement/2_analysis/source/bmy_rolling_stability.py` (new, ~150 lines)
  - `markups-procurement/2_analysis/make.sh` (add line)
  - `markups-procurement/4_paper/source/markups_procurement.tex` (new figure reference, likely in §5 robustness)
  - `markups-procurement/2_analysis/output/figures/bmy_rolling_stability.pdf` (generated)
  - `markups-procurement/2_analysis/output/tables/bmy_rolling_stability.tex` (generated)

**BMY-R2 — Extreme-overhead trimming robustness** (≈2 hours)
- **Problem:** BMY's `Analyze_Missing.do` documents that firms with extreme SGA/COGS ratios drive markup-level estimates. Chadim can implement the analog: trim on residual overhead $O_{it} = II_{it} - \text{COGS}_{it}$ relative to $GO_{it}$.
- **Fix:**
  1. Add a robustness cell in `paper_results.py` wrapping a trimmed-sample call to `ACFEstimator.fit()` that drops observations with $O/GO > 0.80$ (approximately the 95th percentile).
  2. Report the premium alongside the baseline in an existing robustness table, or add a new row.
  3. The expected finding — premium stable, markup *levels* move — supports the "levels vs differentials" thesis.
- **Files:**
  - `markups-procurement/2_analysis/source/paper_results.py` (add ~20 lines)
  - `markups-procurement/4_paper/source/markups_procurement.tex` (new robustness row or footnote)
  - `markups-procurement/2_analysis/output/tables/bmy_overhead_trim.tex` (generated)

**What NOT to do for BMY:**
- Do not port the Stata `.do` files literally — they are tuned to Compustat's SGA/XSGA separation which Chadim's data bundles into II
- Do not attempt an SGA-sensitivity test (COGS-only vs COGS+SGA) — Czech accounting does not separate overhead from other intermediates; flag as a data limitation in §5 instead
- Do not add a DLEU-vs-BMY trend-decomposition figure — that's a different paper's question; your in-memory BMY within/between decomposition is sufficient background

---

### 3.4 NKP (~45 min, 2 actions)

**NKP-R1 — Cite NKP as a proxy-free alternative estimator** (≈20 min)
- **Problem:** NKP is not cited anywhere in the paper. This is a frontier PF-estimation paper (NBER WP 32725, July 2024) that directly attacks the OP/LP/ACF monotonicity assumption.
- **Fix:** Add a 2–3 sentence passage in §3.4 or §5.1 (spec-sensitivity discussion). Suggested text:

  > *"An alternative approach by Norris Keiller, de Paula, and Van Reenen (2024) avoids the proxy-variable monotonicity assumption entirely by using firm-level probabilistic expectations of future output and inputs — available in recent UK Office for National Statistics surveys — to recover production function parameters without a parametric Markov for productivity. Their method is robust to optimization error and unobserved input-price variation that would compromise OP/LP/ACF, but requires firm-level subjective expectations data unavailable in the Czech setting. I therefore rely on the Ackerberg–Caves–Frazer estimator with Andrews–Barahona–Gentzkow–Rambachan–Shapiro (2025) strong-exclusion diagnostics as the feasible-data alternative to NKP's ideal-data estimator."*

- **Bibliography:** Add `\bibitem{NorrisKeillerdePaulaVanReenen2024} Norris Keiller, A., Á. de Paula, and J. Van Reenen (2024). "Production Function Estimation Using Subjective Expectations Data," NBER Working Paper 32725.`
- **Files:**
  - `markups-procurement/4_paper/source/markups_procurement.tex` (§3.4 or §5.1 + bibliography)

**NKP-R2 — Acknowledge the optimization-error vulnerability in limitations** (≈25 min)
- **Problem:** Chadim's cost-share identification argument at line 173 implicitly assumes firms' COGS choices satisfy frictionless cost minimization. NKP footnote 4 explicitly attacks this assumption as it applies to GNR. A referee will notice.
- **Fix:** Add one sentence in the limitations section (or §6.2 near the GNR identification discussion at line 173). Suggested text:

  > *"A final caveat concerns the GNR identification argument itself: the claim that the premium is pinned down by the observed cost share $\alpha^{\text{cogs}}_{it}$ presumes that firms' COGS choices satisfy frictionless cost minimization. If Czech construction firms systematically deviate from this benchmark — through bulk cement discounts, VAT-timing inventory adjustments, or working capital constraints at bid-cycle boundaries — then $\alpha^{\text{cogs}}_{it}$ is a noisy proxy for the true flexible input share, and the cross-firm identification of the premium inherits that noise. The estimator of Norris Keiller, de Paula, and Van Reenen (2024) is robust to precisely this failure mode, but requires firm-level expectations data that is unavailable here."*

- **Rationale:** This turns a latent referee objection into a disclosed limitation with a pointer to the frontier solution. It also pairs the paper's ABGRS strong-exclusion defense with an orthogonal frontier defense, which is rhetorically stronger than ABGRS alone.
- **Files:**
  - `markups-procurement/4_paper/source/markups_procurement.tex` §6.2 or limitations section

**What NOT to do for NKP:**
- Do not attempt to implement NKP's estimator — Czech firm-level subjective expectations data do not exist, and the CZSO business tendency survey is aggregate/quarter, not firm-level linkable to MagnusWeb IDs
- Do not promote NKP to "fourth closest paper" in the intro — its substantive focus is UK TFP and jobs growth, not procurement markups; it belongs in §3.4/§5.1 methodological discussion, not the framing triplet
- Do not substitute realized lagged values for expectations — that would be neither NKP nor ACF and would be indefensible

---

## 4. Consolidated priority table

All 10 actions ranked by (impact × feasibility / hours). Expected **total: ~12.75 hours** with zero headline re-estimation.

| Rank | Action | Paper | Type | Hours | Impact |
|-|-|-|-|-|-|
| 1 | **BMY-R1**: rolling-window PF coefficient stability figure | BMY | New script + figure | ~4 | **High** — directly operationalizes "levels vary, treatment stable" thesis |
| 2 | **Beer-R3**: defend translog vs near-Leontief (option 3b) | Beer | Narrative | ~1.5 | **High** — removes sharpest internal tension |
| 3 | **Beer-R1**: cite Beer paper in 3–5 places (expand triplet to four) | Beer | Narrative + bibliography | ~0.5 | **High** — closes orphan-bibitem gap |
| 4 | **NKP-R1 + NKP-R2**: cite NKP + acknowledge optimization-error vulnerability | NKP | Narrative + bibliography | ~0.75 | **Medium-High** — pairs ABGRS with orthogonal frontier defense |
| 5 | **KLMS-R1**: re-enable `klms_analysis.py` + appendix comparison table | KLMS | Pipeline + 1 appendix table | ~2 | **Medium** — elevates KLMS from citation to reproducible benchmark |
| 6 | **Beer-R2**: reframe §7 pass-through as Beer Figure 6 analog | Beer | Narrative | ~1 | **Medium** — converts welfare section into methodologically grounded argument |
| 7 | **BMY-R2**: extreme-overhead trimming robustness | BMY | New robustness row | ~2 | **Medium** — defensive against referee objections |
| 8 | **Beer-R4**: "why not demand-side triangulation" disclaimer | Beer | Narrative | ~0.5 | **Medium** — closes a 1-point blind spot |
| 9 | **Beer-R5**: Czech declining-markup contrast with Beer | Beer | Narrative | ~1 | **Medium** — elevates a memory-only fact to a main-text headline |
| — | Grand total | — | — | **~12.75 h** | — |

---

## 5. Recommended sequencing for a future implementation pass

**Phase 1 — Narrative-only, zero-risk batch (~5 hours)**

Start here. Nothing requires code execution; all changes are in `markups_procurement.tex`. The paper should still compile cleanly (`pdflatex` twice) with the same page count ± 2 pages.

1. **Beer-R1** (citations + expand triplet to four)
2. **Beer-R3** (translog defense passage, option 3b)
3. **Beer-R4** (demand-triangulation disclaimer)
4. **NKP-R1** (frontier alternative citation)
5. **NKP-R2** (optimization-error vulnerability)

**Verification checkpoint after Phase 1:**
- `grep -c 'DeLoeckerScott2025' source/markups_procurement.tex` → **3+**
- `grep -c 'NorrisKeillerdePaulaVanReenen2024' source/markups_procurement.tex` → **2+**
- `pdflatex source/markups_procurement.tex` twice → zero citation warnings, page count 70 ± 2
- Read §6.2 in isolation: the translog/Leontief tension should now be explicitly addressed

**Phase 2 — BMY high-ROI analysis (~4 hours)**

6. **BMY-R1** (rolling-window PF stability figure) — this is the single highest-impact new analysis

**Verification checkpoint after Phase 2:**
- `bmy_rolling_stability.pdf` exists in `2_analysis/output/figures/`
- `bmy_rolling_stability.tex` exists in `2_analysis/output/tables/`
- Figure compiles into the paper
- $\hat\beta_v^{\text{linear}}$ visibly drifts across windows; premium band stays flat (this is the expected finding that *validates* the paper's thesis — if it doesn't happen, investigate before shipping)

**Phase 3 — KLMS re-enablement (~2 hours)**

7. **KLMS-R1** (uncomment `make.sh:90`, add appendix table)

**Verification checkpoint after Phase 3:**
- `bash 2_analysis/make.sh` runs `klms_analysis.py` cleanly
- `klms_comparison.tex` exists in `4_paper/input/tables/`
- Appendix table renders in the paper with matching values to `memory/klms_study.md`

**Phase 4 — Defensive robustness (~2 hours)**

8. **BMY-R2** (extreme-overhead trimming)

**Verification checkpoint after Phase 4:**
- New row or footnote in the spec robustness table shows premium stable after $O/GO > 0.80$ trimming
- Markup levels visibly move (expected)

**Phase 5 — Rhetorical polish (~1.5 hours)**

9. **Beer-R2** (§7 pass-through reframing)
10. **Beer-R5** (Czech declining-markup contrast)

**Final verification:**
- `pdflatex` twice, zero warnings, page count 72 ± 2 (growth from Phase 2's new figure + Phase 1 narrative)
- Re-run `/jmp-readiness` — expect +3 to +5 points on literature engagement sub-score

---

## 6. File paths for implementation

### Files likely to be edited
| File | Sections touched | Actions |
|-|-|-|
| `markups-procurement/4_paper/source/markups_procurement.tex` | §1 intro (~line 68), §3.4 (~line 170), §5.1 (~line 170), §5.2 (~line 466), §6.2 (~line 476), §7 (~line 577), limitations, bibliography | Beer-R1–R5, KLMS-R1 table ref, BMY-R1 figure ref, BMY-R2 row, NKP-R1, NKP-R2 |
| `markups-procurement/2_analysis/make.sh` | Line 90 (uncomment); add new line for `bmy_rolling_stability.py` | KLMS-R1, BMY-R1 |
| `markups-procurement/2_analysis/source/paper_results.py` | New trimmed-sample call | BMY-R2 |
| `markups-procurement/2_analysis/source/acf_estimator.py` | Confirm `years` subsetting argument works for rolling windows | BMY-R1 (read-only, likely no edits) |

### Files likely to be created
| File | Purpose | Action |
|-|-|-|
| `markups-procurement/2_analysis/source/bmy_rolling_stability.py` | Rolling-window PF coefficient estimation by NACE | BMY-R1 |
| `markups-procurement/2_analysis/output/figures/bmy_rolling_stability.pdf` | Auto-generated | BMY-R1 |
| `markups-procurement/2_analysis/output/tables/bmy_rolling_stability.tex` | Auto-generated | BMY-R1 |
| `markups-procurement/2_analysis/output/tables/bmy_overhead_trim.tex` | Auto-generated | BMY-R2 |
| `markups-procurement/4_paper/input/tables/klms_comparison.tex` | Emitted by `klms_analysis.py` | KLMS-R1 |

### Files that should NOT be touched
- `markups-procurement/2_analysis/source/klms_analysis.py` — already implements all 6 KLMS methods; no code edits needed
- Any file under `2_analysis/output/` except the generated outputs above
- `references/replications/KLMS_AER_repkit/` — reference material only
- `references/replications/BMYReplication/` — reference material only
- Headline premium `0.138` — do not re-estimate; all actions are additive

---

## 7. Cross-references

- **Plan file (Beer paper critical evaluation, standalone deliverable):** `/Users/marek/.claude/plans/cryptic-riding-haven.md`
- **Prior JMP readiness report:** `markups-procurement/4_paper/jmp_readiness_report.md`
- **Memory entry for this synthesis:** `memory/session_apr13b_benchmark_synthesis.md`
- **Existing KLMS replication:** `markups-procurement/2_analysis/source/klms_analysis.py` (886 lines)
- **Existing KLMS study notes:** `memory/klms_study.md`
- **Data codebook:** `markups-procurement/DATA_CODEBOOK.md`
- **Benchmark PDFs (cleaned into consistent locations in this session):**
  - `references/pdfs/markup_estimation/Beer_Markups.pdf` — De Loecker & Scott (2025)
  - `references/pdfs/markup_estimation/KLMS_2025_AER.pdf` — KLMS (2025) AER published
  - `references/pdfs/markup_estimation/KLMS_2025_AER_appendix.pdf` — KLMS online appendix (model derivations)
  - `references/pdfs/markup_estimation/BMY_2026_Comment.pdf` — BMY (2026)
  - `references/pdfs/w32725.pdf` — Norris Keiller, de Paula & Van Reenen (2024)
- **Benchmark replication packages:**
  - `references/replications/KLMS_AER_repkit/` — KLMS R/Python replication kit (AER article ID 220402; was at `io/220402-V1/`)
  - `references/replications/BMYReplication/` — BMY Stata replication package (6 .do files, 1,517 lines)

---

## 8. Known limitations of this synthesis

This synthesis is **as of 2026-04-13** and reflects the paper at ~70 pages, Spec A = translog, JMP readiness 77→~92. Reviewers in a future session should:

1. **Re-verify the orphan-bibitem claim** — `grep 'DeLoeckerScott2025' markups_procurement.tex` and count. If it's no longer 1, Beer-R1 has been acted on in the interim.
2. **Re-check whether `klms_analysis.py` is still DISABLED** — `grep 'klms_analysis' 2_analysis/make.sh`. If it's uncommented, KLMS-R1 is done.
3. **Re-check the §6.2 translog-vs-Leontief tension** — read §6.2 in isolation. If the defense paragraph is already there, Beer-R3 is done.
4. **Confirm the Czech (1−ε) = 0.896 and double-markdown = 0.720 numbers** are still current by re-running `klms_analysis.py` against the latest `data_rebuilt.dta`. Memory notes may be stale if the underlying data were rebuilt.
5. **Watch for tangential Orbis staging** — if `orbis_panel_construction.dta` gets staged to `2_analysis/input/` (per `data_acquisition.md`), KLMS-R1 can be extended to include actual Czech θ identification, invalidating the calibration footnote.

---

## 9. Session metadata

- **Session ID (human label):** session_apr13b
- **Files moved during cleanup:**
  - `io/kroft-et-al-2025-imperfect-competition-...pdf` → `references/pdfs/markup_estimation/KLMS_2025_AER.pdf`
  - `io/23258.pdf` → `references/pdfs/markup_estimation/KLMS_2025_AER_appendix.pdf`
  - `io/220402-V1/` → `references/replications/KLMS_AER_repkit/`
  - `io/Beer_Markups.txt` → `references/pdfs/markup_estimation/Beer_Markups.txt` (pdftotext extraction)
- **Files created during this session:**
  - `/Users/marek/.claude/plans/cryptic-riding-haven.md` — Beer-paper critical evaluation plan
  - `markups-procurement/4_paper/benchmark_synthesis_report.md` — this document
  - `memory/session_apr13b_benchmark_synthesis.md` — short memory entry linking here
- **MEMORY.md** updated with a pointer to `session_apr13b_benchmark_synthesis.md`
- **Current paper state:** `markups_procurement.tex` compiles to 70 pages, zero citation warnings, headline premium 0.138, JMP readiness ~92. No edits to the paper from this session — the user explicitly scoped all benchmark evaluations as standalone deliverables.
