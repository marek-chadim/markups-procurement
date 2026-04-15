# `markups-procurement` × ABGRS (2025 QJE): Project–Framework Mapping

**Compiled:** 2026-04-15 • **Paper version:** `markups_procurement.tex` at 91 pp., 110 bibitems, 0 undefined refs • **ABGRS citation:** Andrews, I., N. Barahona, M. Gentzkow, A. Rambachan, and J. M. Shapiro (2025), "Structural Estimation Under Misspecification: Theory and Implications for Practice," *Quarterly Journal of Economics* 140(4), 1801–1855, DOI: 10.1093/qje/qjaf018.

---

## §1. Executive summary

The Czech construction markups-procurement paper deploys the ABGRS (2025 QJE) misspecification framework as the econometric backbone of its causal claim. Four deployment sites anchor the paper's defense of the 14% procurement markup premium against production-function specification skepticism: (i) the abstract (L41) and intro (L74, L76) advertise ABGRS approximate causal consistency as the transparency argument; (ii) §3.4 `sec:sensitivity` (L187) uses strong exclusion as the bridge between the Pakes (2021) multiproduct critique and the stable 12.5–14.3% treatment effect; (iii) §6.3 `sec:abgrs` (L445–455) contains a dedicated strong-exclusion subsection implementing the partial-$R^2$ diagnostic on every ACF instrument plus ABGRS's Appendix C.3 residualization (Row 8) and Appendix C.4 Chamberlain sieve (Row 9) recipes; and (iv) Appendix A `app:ags_detail` (L753+) reports the AGS (2017) load-bearing moment matrix $\hat\Lambda$ that ABGRS implicitly require for reporting transparency. A project rule (`abgrs-transparency.md`) promotes three ABGRS requirements to mandatory reporting standards for every estimation table. The mapping below documents which framework concepts are fully deployed, which are partial, and where gaps exist. The headline finding: **the paper passes the ABGRS checklist on the core theorem (Prop 3 ACC + dynamic strong exclusion Def 12) but has 5 substantive-but-low-effort gaps** that the mapping's §7 flags as recommended paper upgrades.

---

## §2. ABGRS framework recap

### 2.1 Two-layer model and causal summaries

ABGRS separate the *nesting model* (a potential-outcomes framework $Y_i = Y_i(D_i, X_i)$, $D_i = D_i(X_i, Z_i)$ that respects exogeneity and exclusion of $Z$ but is otherwise nonparametric, p. 1808) from the *researcher's model* (a finite-dimensional $Y_i = Y^*(D_i, X_i, \xi_i; \theta)$ parameterized by $\theta = (\alpha, \beta)$, where $\alpha$ governs causal effects of $D$ on $Y$ and $\beta$ governs how controls $X$ shift the residual, p. 1811). A *causal summary* $\tau^*(\hat\theta)$ is any generalized weighted average of partial derivatives $\partial Y_i(d, X_i)/\partial d$ (Definition 1, p. 1814) — for our paper, the relevant causal summary is the procurement premium $\tau = E[\log\mu_i(D_i=1, X_i) - \log\mu_i(D_i=0, X_i)]$.

### 2.2 Propositions 1–3

**Proposition 1** (p. 1818, negative result). Without a bound on the distance $\delta(G)$ between the true DGP and the researcher's model, no oracle estimator can guarantee bounded error $|\tau^*(\hat\theta) - \tau(G)| \leq b$ uniformly over causal summaries and DGPs. Arbitrarily wrong models deliver arbitrarily wrong answers.

**Proposition 2** (p. 1820, decomposition). Causally correct specification holds if and only if the true DGP admits the representation $Y_i(d, x) = Y^{**}(d, x, \xi_i + L_i(x); \alpha_0)$ for some unit-specific $L_i(x)$. This is permissive: the researcher's assumption about how $X$ enters can be wrong, provided the causal effect of $D$ on $Y$ is approximately parameterized correctly.

**Proposition 3** (p. 1831, main positive result — **Approximate Causal Consistency**). "If conditional exogeneity holds, then any estimator satisfying strong exclusion and strong identification is approximately causally consistent. Moreover, even if unconditional exogeneity holds, any estimator that is approximately causally consistent must satisfy strong exclusion." (p. 1831). The second half is the necessity direction: strong exclusion is not optional.

### 2.3 Strong exclusion (Definition 4, p. 1826)

The researcher's estimator satisfies strong exclusion if the corresponding estimand solves a moment equation $E_G[f^*_G(X_i, Z_i) R^*(Y_i, D_i, X_i; \theta)] = 0$ where $f^*_G = (f^E_G, f^I_G)'$ satisfies three requirements:

1. **Mean independence**: $E_G[f^E_G(X_i, Z_i) \mid X_i] = 0$ — the excluded-instrument moment component must be mean-zero conditional on $X_i$, not just the raw instrument $Z$.
2. **Minimal excluded dimension**: at least $\dim(\alpha) = \dim(\theta) - \dim(\beta)$ linearly independent rows of $f^E_G$.
3. **Maximal included dimension** (overidentified case): no more than $\dim(\beta)$ moments depending only on $X_i$ (Online App C.1 Lemma 6 extends this to $\text{rank}(W_G \Xi_G W_G') \leq \dim(\beta)$).

### 2.4 Dynamic strong exclusion (Online App D.3.4, p. 36)

ABGRS explicitly extend the framework to dynamic production-function settings. The researcher's Cobb-Douglas dynamic model is $Y_{i,j} = \beta_0 + \alpha D_{i,j} + \beta_1 K_{i,j} + \nu_{i,j}$, $\nu_{i,j} = \beta_2 \nu_{i,j-1} + \xi_{i,j}$, where the innovation $\xi_{i,j}$ "is realized after the dynamic input is chosen but before the static input is chosen in period $j \geq 1$, and it is therefore independent of $X_{i,j}$" (Online App p. 33) — the ACF (2015) timing assumption. **Definition 12** (Online App p. 36) requires $W_G f^*(X_i, Z_i) = [W_G^E f^*; W_G^I f^*]'$ with $E_G[W_G^E f_j^*(X_{i,j}, Z_{i,j}) \mid X_{i,j}] = 0$ and a rank condition. The critical verbatim quote (Online App p. 36): *"A choice of instruments in the spirit of Blundell and Bond (1998/2000), see also Ackerberg, Caves, and Frazer [2015] Section 4.3.3) is $f_j^*(X_i, Z_i) = (1, K_{i,j}, D_{i,j-1}, K_{i,j-1})'$... These instruments cannot satisfy dynamic strong exclusion because they are fully determined by $X_i$"*. This is the theoretical anchor for the paper's §6.3 claim that Blundell-Bond instruments are inadmissible without an external shifter.

### 2.5 Enforcing strong exclusion — Appendix C recipes

- **C.2 Coarsened mean independence**: residualize against a coarsening $\chi_j(X_i)$ instead of the full $X_i$, trading theoretical tightness for practical feasibility.
- **C.3 Residualization**: form $f^E_G(X_i, Z_i) = f_G(X_i, Z_i) - E_G[f_G(X_i, Z_i) \mid X_i]$, where the conditional expectation is estimated nonparametrically (ABGRS's replication code uses group-mean $Z - E[Z|X]$ within product-count × $X$-value bins) or by flexible ML.
- **C.4 Automated nested optimization**: a two-step procedure that adaptively reweights the moment conditions. ABGRS's own implementation (`estimate_rcnl.m` L134–190 + `rcnl_dvarcov.m` L72–90) uses a cluster-sandwich empirical moment covariance $S = \sum_m IV_m' (\xi_m \xi_m') IV_m$ and adaptive weight $W_\text{new} = \text{pinv}(SE \cdot S \cdot SE')$.
- **C.4.1 Chamberlain (1987) sieve extension**: projects the observation-level Jacobian onto a sieve basis of the instrument pool to construct efficient instruments.

### 2.6 Reporting standards (§6 "Conclusion", p. 1850)

ABGRS close with a four-point checklist for applied practitioners:

1. "Researchers should report their choice of instruments $f^*(X, Z)$ and the form of their estimator."
2. "Report the degree of misspecification" of the researcher's model against the true DGP, measured as RMS discrepancy of partial derivatives.
3. "When a researcher does not have access to any excluded, exogenous variables, we recommend they make explicit that their estimator fails to satisfy strong exclusion" (p. 1850).
4. "Conduct inference using the techniques described in Section C.3" or equivalent residualization-based inference.

---

## §3. Concept-by-concept mapping table

| ABGRS concept | ABGRS ref | Paper location | Deployment | Notes |
|-|-|-|-|-|
| Two-layer nesting model | Def 1, p. 1808 | — | **Absent** | Framework is invoked verbally but never rendered explicitly in our paper. Acceptable — readers meet the nesting model in ABGRS. |
| Researcher's model $Y_i = Y^*(D, X, \xi; \theta)$ | p. 1811 | §3 Methodology + eq. (markup) | **Full** | Translog PF is the $Y^*$; $\theta = (\beta_k, \beta_c, \beta_{cc}, \beta_{kk}, \beta_{kc})$; $\alpha = \hat\beta_{pp}$; $\beta = $ remaining PF coefs. |
| Parameter partition $\theta = (\alpha, \beta)$ | p. 1811 | implicit | **Partial** | Not called out explicitly. §3.4 L187 treats the premium as the causal parameter of interest but does not name the $(\alpha, \beta)$ split. |
| Causal summary $\tau^*(\hat\theta)$ | Def 1, p. 1814 | Abstract L41, Intro L68, §5 L213 | **Full** | Procurement premium $\hat\beta_{pp} = 0.138$ is the causal summary. Treatment-effect differential $\partial\log\mu/\partial pp$ is exactly a generalized weighted average of partial derivatives. |
| Distance from causally correct spec $\delta(G)$ | Def 3, p. 1818 | — | **Absent** | Paper does not compute a numerical $\delta$ value. Instead it relies on the qualitative argument that the CD/translog gap is a lower bound on $\delta$. **Gap #G2 — see §6**. |
| Prop 1 (unbounded $\delta$ → no guarantee) | p. 1818 | §3.4 L187 | **Partial** | Implied by the "approximate causal consistency" language but not stated as a negative result. |
| Prop 2 (permissive decomposition) | p. 1820 | — | **Absent** | Would be the theoretical justification for "the paper's identification is robust to functional-form error in how controls enter". Not cited. **Gap #G3**. |
| Prop 3 ACC (main positive result) | p. 1831 | Abstract L41, Intro L74, §3.4 L187, §6 L429, §6.3 L447, Concl L556 | **Full** | 6 deployment sites. "Approximate causal consistency" appears verbatim in the abstract, intro, §3.4, §6.3, and conclusion. |
| Strong exclusion Def 4 | p. 1826 | §6.3 L447 | **Full** | §6.3 paraphrases Def 4 in one sentence. **Minor sharpening needed** — paraphrase drops the joint strong-identification requirement and elides the three sub-requirements. |
| Mean-independence requirement | p. 1827 | §6.3 L449 (partial $R^2 < 0.10$) + Row 8 (OLS residualization) | **Full** | Operationalized as partial $R^2 < 0.10$ on 5,236 firm-year obs. Lagged $k$ at 0.03, lagged COGS at 0.01–0.02, lagged procurement at 0.08–0.09 — all pass. |
| Minimal excluded dimension | p. 1827 | §6.3, Table `tab:adl_instruments` Row 9 footnote | **Partial** | The $K_\beta = 6$ efficient instrument combinations in Row 9 satisfy $\dim(\alpha) = 6$ (five translog coefs + premium). Not explicitly argued in prose. **Gap #G4**. |
| Maximal included dimension + Lemma 6 | Online App C.1 | — | **Absent** | The Kim–Luo–Su (2019) deeper-lag overidentification $(k_{t-1}, c_{t-2})$ could in principle violate Lemma 6. Not verified in the paper. **Gap #G1**. |
| Dynamic nesting model (Def 11) | Online App D.3.1, p. 32 | §6.3 L447 ("Section IV.D extends the framework") | **Full** | The paper cites the dynamic extension correctly. |
| Dynamic strong exclusion (Def 12) | Online App D.3.4, p. 36 | §6.3 L449 | **Full** | The paper's verbatim claim ("Blundell-Bond instruments cannot satisfy dynamic strong exclusion because they are fully determined by the state-variable control set $X_{i,j}$") *exactly* matches ABGRS Online App p. 36. |
| Cobb-Douglas dynamic example | Online App p. 33, eq. $Y_{i,j} = \beta_0 + \alpha D_{i,j} + \beta_1 K_{i,j} + \nu_{i,j}$ | §3 methodology | **Full** | Our translog is a richer parameterization of the same example; the ACF timing assumption matches ABGRS's "innovation realized after dynamic but before static" (Online App p. 33). |
| External shifter / procurement instrument | §3.4 L187, Intro L74, §6.3 L449 | `pp_{it-1}` in the Markov transition | **Full** | Paper's role for `pp_{it-1}` as the external shifter is exactly the role of $Z$ in ABGRS's Cobb-Douglas example. |
| C.2 Coarsening recipe | Online App C.2 | — | **Absent** | Not needed in our setting — the control set is already coarse (year × NACE-2). Explicit non-adoption; see §8. |
| C.3 Residualization | Online App C.3 | §6.3 L449 (DML cross-fit partial $R^2$) + §6.3 L453 (Row 8 OLS-residualized interactions) | **Full** | Implemented twice: (a) as a cross-fitted ML diagnostic (Lasso, RF, GB) in §6.3 L449, (b) as a 14-instrument OLS-residualized pool in `adl_instrument_comparison.py` Row 8. Both satisfy C.3's spirit. **Minor Gap #G5 — parametric vs nonparametric residualization — see §6**. |
| C.4 Automated nested optimization | Online App C.4 | — | **Absent** | We do not implement the full nested optimization. Explicit non-adoption; see §8. |
| C.4.1 Chamberlain sieve | Online App C.4.1 | §6.3 L455 (Row 9) + `adl_instrument_comparison.py` L277 | **Partial** | Our Row 9 uses pyblp's `optimal_instruments='replace'` flag, which computes efficient instruments via sieve projection onto the Jacobian. ABGRS's own replication uses a cluster-sandwich two-step (`rcnl_dvarcov.m` L90: `W_new = pinv(SE*S*SE')`). The two are related but not identical; the pyblp variant is closer to feasible optimal instruments than to the empirical cluster-sandwich. **Gap #G6**. |
| AGS (2017) Λ scaled sensitivity matrix | §3.2, pp. 1815–1816 + Online App B | Appendix A `app:ags_detail` L753 + §6.4 content | **Full** | $\hat\Lambda$ computed for linear translog coefs $\hat\beta_k, \hat\beta_c$. $|\hat\Lambda_{\beta_c, L.c}| = 2421$ identifies lagged-COGS as the load-bearing moment. **Gap #G7** — $\Lambda$ not computed for the causal summary (premium) itself. |
| J-statistic overidentification | §5 empirical app | `tab:adl_instruments` Rows 6/7 + Row 9 "just-identified" note | **Full** | Hansen $J$ reported for Row 6/7 to flag strong-exclusion violations; Row 9 is just-identified. |
| Local misspecification propagation via delta method | §4, p. 1833 | Appendix A L755 | **Partial** | Informal delta-method propagation: "shift in the procurement premium of under 0.5 percentage points." No formal bound computed. |
| Andrews (2017 ECMA) identification-robust CS | cited as complement | `tab:ags_twostep_ci` L761 | **Full** | S-statistic CS 2–6× wider than Wald for first-stage; second-stage premium regression unaffected. |
| Borusyak-Hull (2023) recentering equivalence | Online App D.1 | — | **Absent** | The equivalence between strong-exclusion residualization and Borusyak-Hull shift-share recentering is not invoked. **Gap #G8** — one-sentence addition would strengthen positioning. |
| Reporting standard #1 (instruments + estimator) | §6 p. 1850 | Throughout | **Full** | Every instrument set is named, every estimator is described. |
| Reporting standard #2 (degree of misspecification) | §6 p. 1850 | §3.4 L187 ("markup levels range 1.02–2.17") | **Partial** | Paper reports the *symptom* of misspecification (level variation) but not a numerical $\delta(G)$ value as ABGRS §6.2 operationalize it. **Gap #G2 (restated)**. |
| Reporting standard #3 (strong-exclusion status) | §6 p. 1850 | §6.3 L449 (explicit partial $R^2 < 0.10$) | **Full** | Paper is explicit that strong exclusion is satisfied. |
| Reporting standard #4 (residualization-based inference) | §6 p. 1850 | §6.3 L449 + L453 | **Full** | Both the partial-$R^2$ diagnostic and the OLS-residualized Row 8 instrument are operational. |
| Project rule `abgrs-transparency.md` | — | `.claude/rules/abgrs-transparency.md` | **Full** | Four reporting requirements promoted to mandatory paper-wide standard: load-bearing moments, partial $R^2 < 0.10$, estimand under misspecification, $J$-statistic for overid specs. |

**Total**: 28 concepts × 4 strength levels. **Full**: 17. **Partial**: 7. **Absent**: 4. **Gaps numbered G1–G8 are elaborated in §6.**

---

## §4. Section-by-section paper annotation

### §0 Abstract (L41)
"…all ACF instruments satisfy the Andrews, Barahona, Gentzkow, Rambachan, and Shapiro (2025) strong-exclusion requirement (partial $R^2 < 0.10$), and the Cobb-Douglas and translog specifications are treated as equivalently valid and yield the same treatment effect."

**Map**: Full deployment of Prop 3 ACC + mean-independence requirement. This is the shop window — any referee skimming sees it first. **Alignment**: Correct but loose. The claim "partial $R^2 < 0.10$ is the strong-exclusion requirement" is a *sufficient condition we impose*, not the ABGRS Def 4 literal requirement (which is mean independence $E[f^E_G \mid X] = 0$, verified in sample via partial $R^2$). The 0.10 threshold is our convention, not ABGRS's; ABGRS do not name a specific cutoff. **Sharpen** by adding "(our convention; ABGRS require the stronger condition $E[f^E_G \mid X_i] = 0$, which partial $R^2 < 0.10$ operationalizes)". Severity: **cosmetic**.

### §1 Introduction L74 + L76
L74 invokes ABGRS as the "transparency argument" that ties together the external-instrument identification (breaking Gandhi-Navarro-Rivers non-identification) with the partial-$R^2$ diagnostic. L76 uses ABGRS strong exclusion as the second leg of the three-feature defense against the BLP-trained IO referee's "why not B&T reduced form?" question.

**Map**: Full deployment. **Alignment**: tight. **Sharpen**: none needed. Both paragraphs are pedagogically accurate.

### §3.4 Specification Sensitivity `sec:sensitivity` L187
"Estimand under misspecification. Andrews, Barahona, Gentzkow, Rambachan, and Shapiro formalize when structural estimates remain reliable despite model misspecification: the estimator achieves 'approximate causal consistency' if the instruments satisfy *strong exclusion*—mean independence from the included controls."

**Map**: Full Prop 3 ACC invocation as the Pakes multiproduct defense. **Alignment**: Two minor issues. First, the paraphrase drops the joint "strong identification" requirement of Prop 3 (though implicit in the word "estimator"). Second, "mean independence from the included controls" is technically the requirement on the *residualized moment* $f^E_G$, not on the raw instrument $Z$. **Sharpen** to "the *residualized moment* $f^E_G(X_i, Z_i)$ is mean-independent of $X_i$" — or, even simpler, cite Def 4 explicitly. Severity: **cosmetic**.

### §6 Robustness intro L429
Cites AGS (2017) Λ in `app:ags_detail` and ABGRS strong-exclusion check in §6.3 as the two diagnostics that "confirm the divergence does not invalidate the causal summary."

**Map**: Full. **Alignment**: tight. The language "causal summary" directly invokes ABGRS Def 1 terminology. No change needed.

### §6.2 Specification Stability `sec:specstab` L437
"Row 6 compresses the premium and inflates the Hansen $J$ because [year-level external shifters] violate the ABGRS mean-independence requirement (§6.3)."

**Map**: Full Def 4 citation applied to instrument pool construction. **Alignment**: tight. Row 6 is our concrete example of a "raw, non-residualized" instrument that fails mean independence — exactly what ABGRS warn against.

### §6.3 Strong Exclusion Under Misspecification `sec:abgrs` L445–L455 (the dedicated subsection)
This is the paper's densest ABGRS deployment. L447 states Prop 3 ACC verbatim. L449 implements the diagnostic with partial $R^2 < 0.10$ across Lasso/RF/GB on 5,236 obs. L449 then invokes Online Appendix D.3.4 dynamic strong exclusion and explains that Blundell-Bond instruments fail — **this claim is now verified verbatim against Online App p. 36**. L453 describes Row 8 OLS-residualized construction. L455 describes Row 9 Chamberlain sieve via pyblp.

**Map**: 8 ABGRS concepts invoked in ~500 words. **Alignment**: Mostly tight, with three fixable issues:
- L449 cites "Online Appendix~C.3" for the flexible residualization equivalence to DML. Correct.
- L449 cites "Online Appendix~D.3.4". Correct (verified this session).
- L455 cites "pyblp two-step procedure (Conlon-Gortmaker 2020)" as the Chamberlain sieve. This is the **Row 9 pyblp-vs-cluster-sandwich gap (#G6)** — the reference implementation in ABGRS's own replication code is not pyblp but a cluster-sandwich two-step GMM. Not wrong, but a more direct citation to ABGRS's implementation would be tighter.

### Appendix A `app:ags_detail` L753–L767
Reports $\hat\Lambda$ for linear translog coefs $\hat\beta_k, \hat\beta_c$. Andrews (2017) identification-robust S-statistic CS. Andrews-Chen-Tecchio √J weight-hacking bound. Armstrong (2025) expanded-model bias-adjusted CI.

**Map**: Full deployment of the ABGRS-adjacent transparency literature. **Alignment**: tight on everything reported, but **Gap #G7**: $\hat\Lambda$ is reported for PF *coefficients*, not for the premium $\hat\beta_{pp}$ itself. ABGRS (Prop 3) frame approximate causal consistency as a statement about the *causal summary* $\tau^*$, not about individual parameters. Chain-rule propagation through $\mu_{it} = (\hat\beta_c + 2\hat\beta_{c^2} c_{it} + \hat\beta_{kc} k_{it})/\alpha_{it}$ would give $\Lambda_{\hat\beta_{pp},\cdot}$ directly. Informal delta-method propagation is already noted at L755 ("shift in the procurement premium of under 0.5 percentage points"); Gap #G7 recommends formalizing this into a $\Lambda$ row for the premium.

---

## §5. Replication-code benchmarks

Our `markups-procurement/2_analysis/source/adl_instrument_comparison.py` Rows 6–9 map onto ABGRS's reference implementation as follows.

**Row 6 (raw interactions, no residualization)**: corresponds to ABGRS's "baseline" in Figures 2–3 of the main paper (pp. 1839–1844) — the estimator that becomes "severely median-biased" under misspecification. Our paper correctly reports that Row 6 fails the ABGRS mean-independence requirement and inflates Hansen $J$; this is by construction, the same pedagogical role Row 6 plays in ABGRS's own exposition.

**Row 8 (OLS-residualized interactions)**: our implementation (`adl_instrument_comparison.py` L192–207) uses `_residualize(series, X)` which computes $Z_\text{new} = Z - X\hat\beta$ via `np.linalg.lstsq` on `ctrl_df = [k, L_cogs, pp_dummy, year FE, NACE FE]`. This is **parametric** residualization. ABGRS's own `residualizer.m` implements **nonparametric** residualization: group markets by product-count bins (10 bins), then by identical $X$-value bins within each, then compute $Z - E[Z \mid X\text{-group}]$ (L8–77). The two are asymptotically equivalent when $X$ enters linearly and exhaustively, which holds in our setting because `ctrl_df` is linear in the first-order translog basis. **Alignment**: conceptually correct; empirically could be sharpened by cross-checking against a nonparametric variant. ABGRS justify both in Appendix C.3 ("the conditional expectation can be estimated by OLS when $X$ enters linearly"). Severity: **cosmetic / no action needed**.

**Row 9 (Chamberlain sieve)**: our implementation triggers pyblp's `optimal_instruments='replace'` (L277, `CHAMBERLAIN_LABELS = {'ext_resid + Chamberlain (optimal)'}`). pyblp's `optimal_instruments` follows Conlon-Gortmaker (2020) pp. 1120–1122: it estimates the Jacobian $\partial\xi/\partial\theta$ at first-step GMM and projects the instrument space onto it, yielding $K_\theta$ efficient combinations. ABGRS's own reference implementation (`estimate_rcnl.m` L134–190 + `rcnl_dvarcov.m` L72–90) is structurally different: it runs first-stage GMM with identity weights → computes $\xi_1$ residuals → builds a cluster-sandwich empirical moment covariance $S = \sum_m IV_m'(\xi_m \xi_m') IV_m$ → forms $W_\text{new} = \text{pinv}(SE \cdot S \cdot SE')$ → re-estimates. Both procedures target optimality but via different routes: pyblp's is feasible optimal instruments (FOI), ABGRS's is cluster-robust empirical efficiency. **Alignment**: our citation of "pyblp two-step procedure" is correct for what we implement, but a reader cross-referencing ABGRS's own code will find a different construction. **Gap #G6** — one-sentence addition in §6.3 L455 footnote acknowledging the two variants would close the loop.

**Comparison not attempted**: ABGRS's empirical illustration is Miller-Weinberg 2017 beer demand, which is proprietary data. We cannot re-run the simulation on Czech data. The benchmarking is conceptual, not numerical.

---

## §6. Gap analysis & audit findings

Eight gaps total, numbered G1–G8, rated by severity.

| # | Gap | Severity | Location | Fix effort |
|-|-|-|-|-|
| **G1** | Overidentification rank condition (Lemma 6, Online App C.1) not verified. Kim-Luo-Su (2019) deeper lags $(k_{t-1}, c_{t-2})$ add moments that *could* violate the maximal-included-dimension requirement if the extra lag structure generates functional dependence on the firm × year cell. | **Substantive** | §4.1 L201 + §6.4 L534 | ~30 min: compute $\text{rank}(W_G \Xi_G W_G')$ with $\Xi_G$ estimated from the first-stage OLS residualization; expected to hold by construction since deeper lags are time-shifted, not $X$-determined. |
| **G2** | Degree-of-misspecification $\delta(G)$ not computed. ABGRS §6 recommend reporting a numerical $\delta$ value; we report only the qualitative symptom (markup levels 1.02–2.17 across specs). | Substantive | §3.4 L187 | ~1 h: compute RMS discrepancy of $\partial\log\mu/\partial c_{it}$ between CD and translog. Result should be small because the premium is robust. |
| **G3** | Proposition 2 (permissive decomposition of $L_i(x)$) not cited. This is the theoretical justification for "the premium is identified from cost-share variation even if we are wrong about how NACE/year FE enter". | Cosmetic | §3.4 L185 | ~10 min: one sentence. |
| **G4** | Minimal-excluded-dimension requirement not explicitly argued. We have $\dim(\alpha) = $ [2 translog linear coefs + 3 curvature + premium = 6] and Row 9 delivers $K_\beta = 6$ efficient combinations. This is exactly $\dim(\alpha)$ — worth stating. | Cosmetic | §6.3 L455 | ~5 min: one clause in the Row 9 description. |
| **G5** | Parametric-vs-nonparametric residualization not flagged. Our Row 8 uses OLS; ABGRS's replication uses group-mean. Asymptotically equivalent in our linear-$X$ setting but worth a footnote. | Cosmetic | §6.3 L453 | ~5 min: one footnote sentence. |
| **G6** | Row 9 pyblp-vs-cluster-sandwich distinction. Our implementation uses pyblp FOI; ABGRS's uses cluster-sandwich empirical efficiency. Both valid, but the citation should acknowledge that the two are distinct procedures. | **Substantive** | §6.3 L455 | ~10 min: one-sentence footnote identifying the two variants and justifying the pyblp choice (pyblp's Jacobian-projection is simpler to implement on our translog first-stage; cluster-sandwich gains matter mainly when the moment variance structure is far from homoskedastic). |
| **G7** | AGS $\hat\Lambda$ computed for PF coefs, not for the premium $\hat\beta_{pp}$. ABGRS (Prop 3) frame ACC as a statement about the causal summary, not about individual parameters. The informal delta-method propagation noted at L755 ("under 0.5 pp shift") is the right object but is not operationalized as a $\Lambda$ row. | **Substantive** | `app:ags_detail` L755–L761 | ~45 min: chain-rule propagation $\Lambda_{\hat\beta_{pp}, m} = \sum_\theta (\partial\hat\beta_{pp}/\partial\theta) \Lambda_{\theta, m}$ through the markup formula. Expected to recover the "<0.5 pp" informal claim as a quantified row. |
| **G8** | Borusyak-Hull (2023) recentering equivalence not invoked. ABGRS Online App D.1 shows the equivalence explicitly. This is a positioning opportunity: our Row 8 residualized interactions *are* a Borusyak-Hull recentering of the external shifters, which ties us to the shift-share IV literature. | Cosmetic | §6.3 L453 + §6.3 L455 | ~10 min: one sentence citing Borusyak-Hull (2023) as an alternative framing. (Note: we have previously deliberately avoided shift-share language because the paper uses continuous-exposure DiD instruments in Appendix B — careful wording needed.) |

**Severity distribution**: 3 Substantive (G1, G2, G6, G7), 5 Cosmetic (G3, G4, G5, G8, G1 partial). None blocking.

**Gap #G6 and Gap #G7 together are the most promotable**: G6 is a citation-precision issue at a single line, G7 is a genuine quantitative upgrade that converts an informal delta-method claim into a named row of the transparency matrix — which is itself a reporting standard ABGRS recommend.

---

## §7. Recommended paper edits (prioritized)

Ordered by leverage (impact per minute of effort). **None of these are executed in this plan; they are recommendations for a follow-up task.**

### Tier 1 — Single-session upgrades (~90 min total)

**R1 (priority)** — §6.3 L449: add a one-clause strong-exclusion sharpening.
- Current: "strong exclusion—mean independence from the included controls."
- Proposed: "strong exclusion—the residualized instrument moment $f^E_G(X_i, Z_i) = f_G(X_i, Z_i) - E_G[f_G(X_i, Z_i) \mid X_i]$ is mean-independent of $X_i$ by construction."
- *Why*: tightens the paraphrase to the ABGRS Definition 4 literal requirement. Severity: cosmetic but pedagogically important for referees who re-read Def 4. Gap #G3, G4.

**R2** — §6.3 L455: add a footnote distinguishing pyblp FOI from cluster-sandwich Chamberlain.
- New footnote: "Our Row 9 triggers pyblp's `optimal_instruments='replace'` flag (Conlon-Gortmaker 2020), which implements feasible optimal instruments via Jacobian projection. ABGRS's own replication implementation (their `estimate_rcnl.m` + `rcnl_dvarcov.m`) uses a cluster-sandwich empirical moment covariance $S = \sum_m IV_m'(\xi_m \xi_m') IV_m$ with $W_\text{new} = \text{pinv}(SE \cdot S \cdot SE')$. Both constructions target C.4.1-style efficiency; the pyblp variant is simpler for our translog first-stage and delivers the same 2.4× SE reduction reported in Row 9."
- *Why*: closes Gap #G6. Prevents any careful reader from flagging the citation as loose. ~10 min.

**R3** — `app:ags_detail` L755: add a computed $\Lambda$ row for the premium.
- New short paragraph after L755: compute $\Lambda_{\hat\beta_{pp}, m} = \sum_\theta (\partial\hat\beta_{pp}/\partial\theta) \Lambda_{\theta, m}$ via chain rule through $\mu_{it} = (\hat\beta_c + 2\hat\beta_{cc} c + \hat\beta_{kc} k)/\alpha_{it}$ and the OLS second-stage $\hat\beta_{pp} = [\sum (\log\mu - \bar{\log\mu})(pp - \bar{pp})]/[\sum(pp - \bar{pp})^2]$. Expected result: $\Lambda_{\hat\beta_{pp}, L.c} \approx 10$–20, dominating all other moments, consistent with the "<0.5 pp delta-method shift" informal claim.
- *Why*: closes Gap #G7 and promotes the informal claim to a formal transparency row. ~45 min if `dml_core.py` can be reused; slightly more if bespoke.

**R4** — §3.4 L185: one sentence citing Prop 2 as the permissive-decomposition justification.
- New sentence: "The $\log\mu$-identification argument is robust to the functional form of the NACE-interacted fixed effects: Andrews et al. (Proposition 2) show that causally correct specification requires only that the $\alpha$ parameters are correct; the $L_i(x)$ component governing how controls shift the residual is permissive."
- *Why*: closes Gap #G3. ~10 min.

**R5** — `\bibitem{ABGRS2025}`: update citation from "forthcoming" to final.
- Current: `\emph{Quarterly Journal of Economics}, forthcoming.`
- Proposed: `\emph{Quarterly Journal of Economics}, 140(4), 1801--1855. DOI: 10.1093/qje/qjaf018.`
- *Why*: citation hygiene. ~1 min.

### Tier 2 — Optional deeper upgrades (~3 h total)

**R6** — Compute a numerical $\delta(G)$ value (Gap #G2). Run the CD-vs-translog RMS discrepancy of $\partial\log\mu/\partial c_{it}$ and report it as one sentence in §3.4. Expected to be small (<0.05) because the premium is specification-robust.

**R7** — Verify Lemma 6 overidentification rank condition (Gap #G1). Compute $\text{rank}(W_G \Xi_G W_G')$ using the first-stage residualization already in place. Not strictly required but would close the last Absent-level gap.

---

## §8. Non-adoptions

Three ABGRS components we *deliberately do not* use, with justifications.

1. **C.2 coarsened mean independence**. Not needed — our control set (year × NACE-2 FE, plus firm fixed effects) is already much coarser than ABGRS's Miller-Weinberg product × month × income set, so full residualization against `ctrl_df` is feasible without exploding variance. C.2 is explicitly a compromise for when full residualization is infeasible; it is not.

2. **C.4 full automated nested optimization**. Not implemented because our GMM setting is low-dimensional (6 parameters in Row 9) and the pyblp FOI variant of C.4.1 captures the efficiency gain without the additional coding complexity of cluster-sandwich reweighting.

3. **Formal numerical $\delta(G)$ reporting per §6 standard**. We report the *symptom* (level variation 1.02–2.17 across specs) rather than a computed $\delta$ value. This is a deliberate choice: in our setting, the level variation is a more interpretable misspecification diagnostic than the raw $\delta$ because it maps directly onto the reader's prior about "how wrong could the PF be." Gap #G2 flags this as an optional upgrade if the user wants literal §6 compliance.

All three non-adoptions are internally consistent with the paper's approach and should be stated explicitly if the mapping is ever externalized to a referee or co-author.

---

## Appendix: Deployment density

- **ABGRS citation count**: 9 `\cite{ABGRS2025}` invocations in `markups_procurement.tex`.
- **Dedicated subsection**: §6.3 `sec:abgrs`, ~500 words.
- **Project rule**: `.claude/rules/abgrs-transparency.md` — 4 mandatory reporting requirements for every estimation table.
- **Deployment sites by section**: Abstract (1), Intro (2), §3.4 Sensitivity (1), §6 Robustness intro (1), §6.2 Specstab (1), §6.3 Strong Exclusion (2), Appendix A ags_detail (0 — indirect via AGS 2017 Λ), Conclusion (1).
- **Concept coverage**: 17 Full / 7 Partial / 4 Absent (28 concepts total).
- **Gap count**: 8 numbered (3 substantive, 5 cosmetic, 0 blocking).

*This mapping is a reference document. It does not modify the paper. Use it as a checklist during final JMP polish and as an anchor for conversations about how the project relates to ABGRS 2025 QJE.*
