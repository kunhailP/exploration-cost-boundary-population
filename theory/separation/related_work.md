# Theorem-by-theorem comparison with the closest prior work

Compiled 2026-09-17 from literature checks. "Read" marks the depth at which the prior result was
inspected:
- **full**: theorem statements read in the paper;
- **abstract**: abstract or snippet only. These must be read in full before submission.

Absence claims mean only "not found in the papers we checked".

## Theorem A — information–exploration inequality (population ATE, no extrapolation structure)

| Prior result | Read | What it already gives | What A adds |
|---|---|---|---|
| Simchi-Levi & Wang (AISTATS 2023 / Mgmt Sci), Thm 1: inf max e·√R = Ω(1) | full | The 1/δ² trade-off between regret and ATE error, K-armed, **no contexts**, orders only | Contexts with heterogeneous operational gap g(h). The functional (E√g)² with an explicit constant. Model scope stated (no extrapolation). Adaptive designs covered via van Trees |
| Wan, Kveton & Song (ICML 2022), Lemma 1, Thm 1 | full | Variance-optimal exploration under a safety (value-loss) constraint, solved numerically | Closed form in g; lower bound over all estimators |
| Douglas, Persson & Provost (arXiv 2605.15108), Props 5–6 | full | Neyman-type optimal logging π ∝ π_t√μ for IPW; theory under several information conditions | They minimize variance with no loss price. A prices exploration in g and bounds all estimators |
| Koppel, Bhatt, Zeng & Ganesh (NeurIPS 2025 Constrained Optimization for ML workshop; listed at constrained-opt-ml.github.io/papers) | **abstract only** (2026-09-17 retry: the workshop page has no PDF, OpenReview returned a browser-verification challenge, no arXiv preprint found) | Abstract: extends inverse-gap-weighting contextual bandits so that ATE estimators converge; "explicit tradeoffs between regret and the quality of ATE estimates"; "optimal adaptive hypothesis testing" alongside sublinear regret | Abstract gives no lower bound, no √g allocation (standard IGW weights ∝ 1/gap), nothing on temperature or boundary effects. **Must be read in a browser (openreview.net/pdf?id=KTn5igMOmS) before submission**; it may overlap with the achievability side of Theorem A |
| Cochran optimum allocation with costs | textbook | n_h ∝ N_h S_h/√c_h | The allocation form of A′ is this rule. **We do not claim the allocation as new** |
| Mou, Ding, Wainwright & Bartlett (arXiv 2301.06240), Thm 1 | full | Local minimax lower bound for off-policy linear functionals over convex classes; implies two-point bounds of our type | A couples the information bound to the exploration loss (the Cauchy–Schwarz step). Mou et al. has no loss |
| Hong, Leung & Li (Econometrics J 2020), Thm 1 | full | Finite strata with drifting propensity: n·a_n → ∞ gives a CLT | Superpopulation, general F, all-estimator lower bound |

## Theorem B — single common temperature

| Prior result | Read | What it already gives | What B adds |
|---|---|---|---|
| Douglas, Persson & Provost (2026), Fig. 6 | full | Simulation: softmax logging worse than top-k; MSE → ∞ as it approaches greedy | A rate: loss ≳ n/log²(nδ²) at fixed precision against a constant optimum. Scope limited to one deterministic common temperature per n; random mixtures of finite temperatures can avoid the rate |
| Cesa-Bianchi, Gentile, Lugosi & Neu (NeurIPS 2017), Thms 1–2 | full | Boltzmann exploration with a common schedule can have linear regret; arm-dependent schedules fix it | Different objective (population-effect precision, not regret). B concerns estimation cost at fixed precision |
| Khan & Tamer (Econometrica 2010) | see note below | Irregular identification with unbounded weights; slower-than-√n rates | See the full-text check below |

## Proposition C and Lemma R

| Prior result | Read | Relation |
|---|---|---|
| Cesa-Bianchi et al. (2017) | full | Arm-dependent learning rates as the fix for Boltzmann exploration. C is the analogous context-dependent fix for estimation, with the explicit schedule τ* ∝ \|Δ\| |
| Kantorovich inequality | textbook | Lemma R is a direct application. The contribution is the decision rule E[g]/(E√g)² vs K(ρ) |

## Theorem D — boundary effect at vanishing exploration loss

| Prior result | Read | What it already gives | What D adds |
|---|---|---|---|
| Narita & Yata, *Algorithm as Experiment* (arXiv 2104.12909 v6), Prop 1, Thm 1, eq. (6) | full | e(1−e)-type weights on an approximate propensity score converge to a boundary estimand (Hausdorff measure). Needs undersmoothing nδ_n² → 0; equivalent to local linear RD | The bandwidth is a **system-chosen score-space temperature**, not an analyst-chosen X-space ball. Exploration loss R_n is computed jointly, giving valid inference with R_n → 0 |
| Eckles, Ignatiadis, Wager & Wu (arXiv 2004.09458 v5; Biometrika), noise-induced randomization | full | A sharp cutoff on a noisy running variable with **known noise density and fixed scale**. The estimand is a weighted average of CATE over a latent variable (eq. 6). Hájek-type estimator: consistency (Thm 1), CLT (Thm 2), bias-aware CIs over a sensitivity class (Cor 2). Weights come from a heuristic QP; no rate theorem. No cost or exploration quantity; no TV or indistinguishability result | Boundary CLT: **partial** overlap (design-based weighting and a CLT near a threshold), but no noise scale → 0, no O(τ) bias or nτ³ condition, and the target is over latent U, not the observed H = 0. Vanishing-loss statement and same-model separation: **different** |
| Li, Morgan & Zaslavsky (arXiv 1404.1785 v3) | full | Overlap weights h = e(1−e) (Sec. 4.1, eq. 10), target "ATO". The weights peak where e = 1/2. Exact balance under logistic PS (Thm 3); minimum asymptotic variance under homoskedasticity (Cor 1). Sec. 3 notes that as the PS model sharpens toward 0/1, the weights shrink and precision falls. No limiting estimand, no boundary or RD, no algorithmic or known propensities, no exploration cost | Cite for the weight, "ATO", and "peaks at 1/2". The sharpening limit to a boundary measure and the loss comparison are **not anticipated** |
| Cattaneo–Titiunik–Yu; Imbens & Wager (REStat 2019); Keele & Titiunik (2015) | abstract/full | Boundary RD estimation and rates | D uses no analyst bandwidth and claims no optimal rate |

## The separation (Corollary S; A with D) — the representative claim

Among the papers checked, we found none stating the following in one model (Corollary S), on the same
exploration-loss scale:
- a boundary effect admits valid inference with vanishing cumulative exploration loss;
- the same data cannot consistently estimate the population ATE across model members sharing that
  boundary effect;
- uniform consistency for the population ATE requires divergent loss. Every component has
close antecedents, listed above. The novelty claim is the comparison on a common cost scale and the
common-temperature consequence, **not** any single ingredient.

## Positioning against generative-AI causal work

| Paper | Read | How we position |
|---|---|---|
| Athey, Imbens & Ji (arXiv 2607.05792), Thm 1 and §9 | full | Randomization-based identification using logged unexposed candidates and probability ratios; logprobs or replays estimate assignment probabilities; trimming at κ. Their estimand averages over co-occurring pairs and is trimmed, so it does not contradict A. Our replay remark credits their §9 |
| Fresh & Shin (arXiv 2607.03597) | full | Taxonomy of conversational estimands (assignment, policy, messages, realized features). We cite it and work at the strategy-assignment level |
| Imai & Nakamura (arXiv 2410.00903); Nakamura & Imai (arXiv 2605.07834) | full/abstract | Text-as-treatment with deterministic decoding. We hold the generation kernel fixed and randomize strategies |
| Fong & Grimmer (AJPS 2023); Egami et al. (Sci Adv 2022) | abstract | Latent-treatment confounding in text. Motivates fixing the kernel |

## Khan & Tamer (2010) full-text check
Read: the July 2009 working paper (34 pp., archived from Tamer's former Northwestern page), not the
published Econometrica PDF, which was paywalled. **Theorem numbers may differ in the published version.**
- Setting: weighted moments θ = E[g] with E‖g‖² = ∞; trimmed estimators; iid data; no designs.
  - Sec. 3: Lewbel binary choice.
  - Sec. 4: ATE under unconfoundedness with 0 < p(x) < 1, and the remark that "identification is lost
    when we remove any region in the support of x".
- Thms 3.1 and 4.1: studentized, rate-adaptive CLTs for trimmed estimators, under a Lindeberg condition
  and negligible trimming bias. The variance scale is v(γ) = E[1{|x| ≤ γ}/(p(1 − p))].
- Sec. 4.2: rates for trimmed IPW with known propensities that go to 0 or 1 in the covariate tails, e.g.
  √(n/log n) for logistic regressors with logit propensity, and n^{−1/4} for normal errors with logistic
  regressors.
- Thm A.1: infinite efficiency bound for binary choice (Chamberlain-type).
- No minimax or van Trees bound, no TV/indistinguishability argument, no n·p_min condition, no
  propensity varying with n, no cost or exploration loss.

Verdicts:
- Theorem A (information–exploration inequality): **different**.
- Corollary S (TV at vanishing loss; boundary inference): **different**. "Thin-set identification" is a
  conceptual precursor and should be cited.
- The recoverability characterization and trimmed-IPW rates: **partial**. Same trimming bias–variance
  structure and v(γ), but with a fixed propensity and no iff.

## Still to read in full before submission
Koppel et al. (2025): needs a browser download from OpenReview.

## Recent JCI papers (full-text check, 2026-09-17)

Access: arXiv full texts, except Esterling–Brady–Schwitzgebel (abstract only). Bibliographic data from
Crossref. Theorem numbers for Ravichandran et al. are from arXiv v1 and may differ in the JCI version.

| Paper | Verdict vs (i) boundary CI at vanishing loss, (ii) same-model TV separation, (iii) information–exploration bound, (iv) common temperature | Use |
|---|---|---|
| Ring & Schomaker, JCI 14(1) 20250007, 2026 | (i) related (changes estimand where support exists; heuristic diagnostic, no design), (ii)–(iv) different | cited: formal counterpart for known sharpening policies |
| Humphreys, JCI 13(1) 20240040, 2025 | related: FE estimand weights strata by p(1−p), same as overlap weights; no adaptive design or bounds | cited in boundary-estimand paragraph |
| Ravichandran, Pashley, Libgober & Dasgupta, JCI 12(1) 20230046, 2024 | (iii) **partial overlap on the allocation formula**: budget-constrained A-optimal allocation N_j ∝ S_j/√C_j (Thm 3(a), arXiv v1); no lower bound over estimators or adaptive designs | cited at the attaining allocation; allocation not claimed as new |
| Li & Imai, JCI 12(1) 20230072, 2024 | related: ex-ante vs ex-post ITR evaluation, variance comparison at equal allocation, no cost for off-rule assignments | cited in exploration paragraph |
| Oganisian, JCI 14(1) 20250058, 2026 | different (Bayesian computation of SATE vs PATE) | not cited |
| Lee, Hubbard & Schuler, JCI 14(1) 20240049, 2026 | different (binarized exposures as modified treatment policies) | not cited |
| Azriel, Krieger & Kapelner, JCI 14(1) 20230053, 2026 | different (minimax randomization at fixed allocation) | not cited |
| Esterling, Brady & Schwitzgebel, JCI 13(1) 20240002, 2025 | different (conceptual; abstract only) | not cited |
