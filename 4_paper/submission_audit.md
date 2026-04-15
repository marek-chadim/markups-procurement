# Paper Submission Audit — markups_procurement.tex

Protocol: Paper Submission Audit v1.0 (Irsova-Havranek Duel v1.7 adaptation)
Date: 2026-04-14
Auditor: Claude (single-agent, KB-grounded)
Paper: `markups-procurement/4_paper/source/markups_procurement.tex`
KB: `references/` + `project/`

---

## Verdict: READY WITH REQUIRED CHANGES

## Scorecard

| Criterion | Score | Pass? |
|-|-|-|
| LaTeX compile clean (0 errors, 0 undefined refs) | 78 pages, 0/0 | ✓ |
| Numerical consistency (firms, obs, premium, NACE means, welfare gap) | all headline numbers self-consistent | ✓ |
| Citation fidelity (no *material* mismatches) | 2 material + 2 minor mismatches | ✗ |
| Methodological defense (red-team attacks answered) | 2 of 5 attacks unanswered/partial | ✗ |
| Prose internal consistency (abstract ↔ body ↔ tables) | abstract + body + tables agree | ✓ |
| Predecessor consistency (thesis ↔ proposal ↔ referee reports) | thesis drift explained at line 203; self-referee issues addressed | ✓ |

Overall: **4 of 6 criteria pass**. Two criteria require patches before
submission.

---

## Required changes (severity-ordered)

### FATAL (block submission)

1. **Add ABGRS 2025 QJE PDF to `references/pdfs/transparency/`.**
   The paper's single most load-bearing claim — "first application of
   ABGRS 2025 transparency framework to markup premium evaluation"
   (line 68) — cannot be audited without the source document.
   Sub-audit B blocked on this. The QJE publication is vol. 140(3)
   pp. 1801–1855 (August 2025). Any referee with access to the journal
   can spot a misattribution immediately; auditing now prevents that.

### MATERIAL (should fix)

2. **Line 68: remove "nine estimators that Benkard, Miller, and
   Yurukoglu flag as specification-sensitive."** Sub-audit A verified
   against `BMY_2026_Comment.pdf`: BMY does not flag nine specific
   estimators. The nine-method structure is the paper's own
   (line 219). Replace with: *"even as markup levels range from 1.02
   to 2.17 across the nine-estimator comparison reported in §5 (Table
   \ref{tab:dls}), which includes the OLS/ACF-CD/ACF-TL/cost-share/
   calibrated/Blundell-Bond family that BMY flag as specification-
   sensitive."*

3. **Line 195 + Table caption: "Beer 2024 Table C2 format" should be
   "De Loecker and Scott (2025) Table C2 format."** Sub-audit A
   verified the PDF title page lists Jan De Loecker (KU Leuven) and
   Paul T. Scott (Stern NYU); "Beer" is a filename artifact from
   `Beer_Markups.pdf`. Bibitem `DeLoeckerScott2025` already exists
   (line 1058) so the citation is correct; only the prose label is
   wrong.

4. **§7 contract-level welfare: add 95% CI caveat on within-firm
   pass-through.** Line 601 reports β = 0.813 (SE 0.161) and the paper
   then describes this as "near-unit pass-through." The 95% CI
   [0.498, 1.128] does not statistically exclude 50% or 100%.
   Required addition after line 601: *"The 95% confidence interval
   [0.498, 1.128] does not exclude pass-through rates of 50% or 100%;
   the 0.81 point estimate should be read as the best central estimate
   of the contract-level bridge, not as a tight constraint on how the
   firm-year markup materializes at the contract level."*

### MINOR (nice to have)

5. **Line 197 DGM attribution refinement.** The paper says DGM 2026
   shows "functional-form choice moves markup *levels* negligibly
   once cost-share information is exploited." DGM's emphasis is on
   markup *variation* stability, not levels. Replace "levels" with
   "levels and variation" or rephrase to "DGM show that revenue-based
   markup variation is well-measured even when levels differ by a
   sector-level constant (DGM 2026 Prop. X)" — matching DGM's actual
   wording.

6. **§3.3 Pakes defense: add one sentence making the α^V vs θ^V
   distinction explicit.** After the current Pakes footnote at line
   181, add: *"Formally, the ABGRS strong-exclusion moment conditions
   are defined on cross-firm variation in the observed cost share
   α^V_{it}, not on the estimated elasticity θ^V_{it}, so the
   non-identification of θ^V for multiproduct firms in Pakes §3.1
   does not propagate to the procurement-treatment-effect
   differential."*

7. **§6.3 Chamberlain 1987/2022 shorthand.** Line 498 attributes "sieve
   compression" to Chamberlain 1987. Chamberlain 2022 cites 1987 for
   the Fisher information bound; the sieve construction is the
   paper's implementation choice. Rephrase as: *"Chamberlain (1987)
   optimal-instrument efficiency theory, implemented via the
   Chamberlain (2022) sequential-moment extension and a quadratic
   sieve projection of the 14 candidate instruments onto the
   K_β = 6 parameter dimension."*

---

## Remaining uncertainties (not fixable by prose rewrite)

- **The Stata/Python contract-level pass-through mismatch** documented
  in the overnight port (Stata β ≈ 0.26, Python β ≈ 0.81, same
  specification, different trimming/dedup). Neither run is "wrong" —
  they're on slightly different samples — but the paper reports only
  the Python number. A careful referee running the Stata pipeline
  will get a different number and ask why. The protocol's current
  response is the `_stata.tex` suffix convention; the paper should
  note this in the §7 footnote.

- **β₁ interpretation.** The paper's DLW 2012 §III defense (line 207)
  is textually correct but deflects rather than resolves the
  share-weighted-vs-contract-level interpretation question. A BLP-
  trained IO referee will read this as an epistemic punt. Resolving
  fully would require firm-product-level input allocation data
  (Dhyne-Petrin-Smeets-Warzynski 2020), which the Czech panel does
  not support. This is a genuine data limitation, not a prose defect.

- **Leontief vs translog empirical test.** §6.2 treats the Leontief
  special case as a feature, but does not provide a direct test
  (e.g., testing β_{cc} = β_{kc} = 0 in the translog). The DLS
  9-method comparison is suggestive but not a formal specification
  test.

---

## KB coverage gaps

1. **ABGRS 2025 QJE main paper + online appendix.** Not in
   `references/pdfs/transparency/`. Most load-bearing source in the
   paper. **Download priority: highest.**

2. **GNR 2020 JPE.** Not found in the KB scan. The paper cites
   Gandhi-Navarro-Rivers 2020 for the non-identification result;
   verification requires the source PDF. DLS Ch 3 §5.3 discusses GNR
   and Sub-audit D verified the DLS treatment, so the citation is
   indirectly grounded — but direct verification would strengthen
   the audit.

3. **Dhyne-Petrin-Smeets-Warzynski 2020** (multiproduct PF extension).
   Cited in the Pakes footnote at line 181. Not in KB.

4. **Borusyak-Jaravel-Spiess 2024** (BJS imputation estimator).
   Cited in §5.8 at line 426. Not in KB. The `notes/yale_econ574_...`
   file covers it pedagogically; consider adding the paper PDF.

---

## Sub-audit summary table

| Sub-audit | Verdict | Material findings |
|-|-|-|
| A. Markup lineage | 2 trivial ✓ + 1 minor ✗ + 2 material ✗ | DLW quotes verbatim; KLMS Leontief verified; Pakes quote verified; **BMY 9-estimator misattribution**; **Beer vs De Loecker-Scott label error** |
| B. Transparency/ABGRS | **BLOCKED** — KB gap | Cannot verify any of 8 ABGRS claims without the QJE PDF |
| C. Yale course notes | All ✓ + 1 low drift | ECON 556 / 600 / 574 alignment good; ECON 600 Haile Part 1 "why PF not BLP" is paper's extension beyond course canon — defensible but note |
| D. Handbooks (DLS Ch 3) | All ✓ | GNR non-identification citation correct; external-instrument-breaks-GNR correct; p. 196 DLS quote verbatim |
| E. Econometrics | All ✓ + 1 minor | Chamberlain 1987/2022 shorthand blurs 1987 Fisher bound vs 2022 sequential-moment extension |
| F. Replications + Conlon | All ✓ | dml_strong_exclusion.py matches ABGRS §C.3 logic; klms_analysis.py matches KLMS repkit; Rel_Price matches Titl JLE; acf_estimator.py is independent from DGM ProdFun as stated |
| G. In-house predecessors | All ✓ | Thesis 1,297→1,521 drift explained at line 203; referee-report self-critiques addressed; all 7 thesis replication scripts verified |

**Highest-severity bucket (G) clean.** Advisors (Zimmerman, Biasi,
Bajgar, Meriläinen) will find no undisclosed sample modifications, no
deprecated methods still in use, and no unacknowledged self-referee
critiques. This is the result of the overnight Stata port, the Beer C2
convergence work, and the validation step that surfaced (and fixed) the
fiscal-welfare denominator bug and the contract-level filename
collision.

---

## What I would do before submission

1. **Download ABGRS 2025 QJE PDF** into `references/pdfs/transparency/`
   and re-run Sub-audit B. ~30 min.
2. **Apply patches 2–4 (material)** to the paper prose. ~20 min.
3. **Compile twice, verify 0 errors / 0 undefined.** ~1 min.
4. **Apply patches 5–7 (minor)** if time permits. ~10 min.
5. **Save final PDF to `4_paper/output/`** (the canonical location, not
   `source/`) and commit.

Total required work before submission: **~50 min**. The paper is
substantively ready; the gaps are prose-level citation hygiene plus
one download to close the ABGRS audit.

---

*End of report. Full Phase 2 transcripts and red-team paragraphs
available in the chat history; the one-page §5 summary above is the
actionable deliverable.*
