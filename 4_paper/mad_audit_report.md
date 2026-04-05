## MAD Audit Report: markups_procurement.tex
**Protocol:** Irsova & Havranek (2026) MAD v2.0, adapted for Claude Code
**Date:** 2026-04-04
**Agents:** econometrics-reviewer, io-domain-reviewer, literature-positioning-reviewer, devils-advocate

### Executive Summary

The paper presents a comprehensive empirical study of procurement markups in Czech construction with impressive methodological breadth (ACF, panel FE, DLW, 11 unconfoundedness estimators, SDID, KLMS). The headline 14% procurement premium is, however, **critically sensitive to a single specification choice**: including the procurement dummy in the Markov productivity transition. Without it, the regression-adjusted premium collapses to 0.3% (statistically zero). This specification sensitivity, rather than being a weakness to bury, is the paper's most interesting finding --- but the current framing treats the inflated specification as "baseline" without principled justification. The favoritism decomposition (6% + 8%) is acknowledged calibration but risks being overinterpreted.

---

### Critical Issues (address before submission)

**[CRITICAL #1] The procurement dummy in the Markov transition mechanically inflates the premium from 0.3% to 14%.**
- Source: Agents A, B, D (consensus)
- Evidence: "Including the procurement dummy in the Markov transition...nearly triples the raw procurement premium, from 0.077 (plain ACF) to 0.186" (line 222); "No pp in Markov ... Reg. Premium 0.003 (0.002)" (Table 8 / line 610-611)
- Cross-exam: AGREED by 3/4 agents. The Convergence STRENGTH does not survive this critique because all downstream causal methods (panel FE, DLW, unconfoundedness, SDID) use the same generated markup from the Markov-with-pp specification. The convergence is across Step 2 methods, not across Step 1 specifications.
- **Recommended action:** Reframe the paper around the specification sensitivity as the core contribution. Present results for BOTH Markov specifications as co-equal. Discuss why the plain ACF premium is near zero (potentially because productivity absorbs the treatment effect when pp is excluded from the Markov --- i.e., procurement-induced productivity gains are real but attributed to omega rather than appearing in markups). This is the paper's most intellectually interesting finding.

---

### Major Issues (address in revision)

**[MAJOR #1] The favoritism decomposition (14% = 6% + 8%) is calibration from incomparable quantities.**
- Source: Agents A, C, D (consensus)
- Evidence: "approximately 6 percentage points reflect favoritism" (line 426); "benchmarking against Baranek and Titl (2024) rather than a direct test" (line 498)
- Cross-exam: AGREED. The 6% is contract-price overpricing; the 14% is a production-function markup ratio. These are different economic objects. Subtracting them is not formally justified.
- **Recommended action:** Reframe as "suggestive calibration" explicitly in the text (not just the conclusion). Note that contract-price overpricing and revenue-based markups need not be additively decomposable. The competition heterogeneity (single-bid vs multi-bid) and reform timing tests are stronger evidence of the favoritism channel than the arithmetic decomposition.

**[MAJOR #2] COGS markup levels (1.7--2.35) strain plausibility for construction.**
- Source: Agent B
- Evidence: "NACE 41 (Buildings) 2.14...NACE 42 (Civil Eng.) 2.35" (lines 580-581)
- Cross-exam: CONTESTED by literature reviewer --- DGM (2026) show COGS-based markups are systematically high across all methods. But industry profit margins of 3-8% are hard to reconcile with 114-135% markups over marginal cost.
- **Recommended action:** Add a paragraph reconciling markup levels with profit margins. Key point: markups measure price/marginal cost, not profit margins; fixed costs (land, equipment, overhead) drive the wedge. Cite Hall (2018) on this distinction.

**[MAJOR #3] GNR non-identification implies the premium is identified from alpha (data), not theta (estimation).**
- Source: Agents A, D
- Evidence: "proxy variable methods generically fail to identify the flexible input elasticity for gross output production functions" (line 162); variable input choice flips premium sign (line 220)
- Cross-exam: AGREED that theta is not identified. CONTESTED whether the "common bias" assumption holds: D argues the II sign-flip is direct evidence of differential bias. This is the deepest unresolved methodological issue.
- **Recommended action:** Elevate this from a robustness check to a core discussion point. The premium's identification from alpha (data) rather than theta (estimation) is a feature, not a bug --- but only for COGS. Explicitly state: "the premium is identified from cross-firm variation in expenditure shares, which is observed, not estimated."

**[MAJOR #4] SUTVA violation is acknowledged but unquantified.**
- Source: Agent D
- Evidence: "government projects displace 27% of private market output" (line 261)
- Cross-exam: CONTESTED --- the 27% is from US data (KLMS) and may not apply to Czech construction. But the direction is correct: crowding out depresses control markups, inflating the gap.
- **Recommended action:** Add a bounding exercise. If controls lose X% of output to crowding, what's the implied bias? Even a rough calibration (e.g., 27% displacement x average markup x procurement share) would show whether the bias is 1pp or 5pp.

**[MAJOR #5] Scalar unobservable: Hansen J rejects for NACE 41 (p=0.004).**
- Source: Agent B
- Evidence: "The Hansen J test rejects for NACE 41 (p = 0.004)" (line 378/632)
- Cross-exam: AGREED. NACE 41 is the largest sub-industry (4,131 obs). The rejection suggests the ADL correction is insufficient or the market definition is too coarse.
- **Recommended action:** Report the J-test results prominently and discuss. Consider finer market definitions (NACE 4-digit or region x NACE) as robustness.

**[MAJOR #6] KLMS framework uses non-standard estimand.**
- Source: Agent B
- Evidence: "theta = 0.13 (labor supply elasticity 1/theta = 8.0)" (line 438-439)
- Cross-exam: CONTESTED --- the Wald DiD is standard in the KLMS framework, but the exclusion restriction (procurement shifts labor demand but not supply) is untested.
- **Recommended action:** Discuss the exclusion restriction explicitly. Note that if procurement shifts both supply and demand (e.g., construction workers prefer government project stability), theta is biased.

**[MAJOR #7] DLEU Reply (2025) connection underexploited.**
- Source: Agent C
- Evidence: "output elasticities are remarkably stable over time---confirming that the markup trend is driven by the expenditure share" (line 74)
- **Recommended action:** Add a standalone paragraph connecting DLEU's finding (theta stable, all action in alpha) to this paper's specification sensitivity result. This strengthens the paper's contribution.

**[MAJOR #8] Missing references: Doraszelski & Jaumandreu (2018), Flynn, Huo & Kako (2019).**
- Source: Agent C
- **Recommended action:** Add both to bibliography and cite in the ADL/imperfect competition discussion.

---

### Minor Issues (consider addressing)

**[MINOR #1] Unconfoundedness analysis on 26 firms is underpowered.**
- Source: Agents A, D
- Evidence: "Balanced panel of 26 firms, 312 observations" (line 659/282)
- **Recommended action:** Acknowledge power limitations explicitly. Consider presenting this as supplementary rather than a primary identification strategy.

**[MINOR #2] European markup literature positioning gap.**
- Source: Agent C
- **Recommended action:** Add 1-2 sentences citing Diez et al. (2018) or Cavalleri et al. (2019) on European markup trends for context.

---

### Verified Strengths

**[STRENGTH #1] Causal evidence triangulation is impressive.** (Confirmed by 4/4 agents)
- Four identification strategies with different assumptions yield 11-16% premiums. Near-symmetric entry/exit (11.5% vs 11.2%) is particularly compelling. Oster delta* = -5.9 is very strong.
- Caveat: all methods use the same generated markup, so convergence is across Step 2, not Step 1.

**[STRENGTH #2] Frontier positioning is excellent.** (Confirmed by 3/4 agents)
- Engages KLMS (2025), DGM (2026), Demirer (2025), ADL (2024), BMY (2026) directly. DGM cross-validation is exemplary transparency.

**[STRENGTH #3] The specification sensitivity analysis is the paper's most original contribution.** (Confirmed by 4/4 agents)
- Documenting that the Markov transition choice dominates all other specification dimensions (functional form, estimator, variable input) is genuinely new for treatment evaluation.

---

### Irreducible Disagreements

**[IRREDUCIBLE #1] Is the baseline Markov specification (with pp) justified?**
- Position A (Agents A, B, D): No --- it is circular. Conditioning on pp in the Markov absorbs the very variation the paper claims to measure. The "plain ACF" premium of 0.3% is the honest estimate.
- Position B (paper's position): Yes --- CWDL (2015) include export status in the Markov; pp is analogous. Excluding pp from the Markov forces productivity to absorb procurement effects, biasing omega and the resulting markups.
- **Recommended:** Present both specifications as co-equal and discuss the economic interpretation of each. The disagreement is about whether procurement affects productivity dynamics --- this is an empirical question the paper cannot resolve with current data.

---

### Prioritized Revision Action List

1. **[CRITICAL]** Reframe around specification sensitivity: present both Markov specs as co-equal (Section 3.3-3.4)
2. **[MAJOR]** Reframe favoritism decomposition as "suggestive calibration" (Section 6.1)
3. **[MAJOR]** Add paragraph reconciling markup levels with profit margins (Section 3.1, cite Hall 2018)
4. **[MAJOR]** Elevate GNR / "premium identified from alpha" discussion (Section 3.4 or new subsection)
5. **[MAJOR]** Add SUTVA bounding exercise (Section 4.1, after Oster bounds)
6. **[MAJOR]** Report and discuss Hansen J rejection for NACE 41 (Section 5 or Appendix)
7. **[MAJOR]** Add DLEU Reply connection paragraph (Introduction or Section 3.4)
8. **[MAJOR]** Add Doraszelski & Jaumandreu (2018), Flynn et al. (2019) to bibliography
9. **[MAJOR]** Discuss KLMS exclusion restriction (Section 7.1)
10. **[MINOR]** Acknowledge unconfoundedness power limitations (Section 4.3)
11. **[MINOR]** Add European markup literature context (Introduction)
