# Boundary side of the separation: Theorem D

> **Alignment with the manuscript.** Under the weak form of (D5), only R_n = O(nτ_n²) holds (e.g.
> g(h) = h² gives R_n ~ 1.803·nτ_n³). R_n ~ nτ_n²κf(0)π²/6 needs the strong form φ′(0) = κ > 0. The
> estimator is set to 0, with interval ℝ, on the event that a weight sum is zero; this event has vanishing
> probability. Lyapunov ratios are O((nτ)⁻¹).
>
> **Notation.** In this file p = p_n(h) = P(A = 1 | H = h). The manuscript (`paper/tex/main_jci.tex`) writes this as e_τ(h), and reserves q_τ(h) = P(A ≠ a*(h) | H = h) = sigmoid(−|h|/τ) for the off-greedy probability. The weights are w_1 = A(1 − e_τ) and w_0 = (1 − A)e_τ; the exploration loss is nE[g q_τ].

Status (2026-09-17, rev. 2): full proof written and self-checked. The audit of commit 6d111fc has been
incorporated: the smoothness remark is corrected, the Step 2 constant fixed, Step 5 notation clarified,
the nτ table corrected, and the coverage column renamed. **Not reviewed by a domain expert.**
The version of this result placed in the same model as the population-side impossibility is
`corollary_S.md`.

This is the one-dimensional score version (Δ(h) = h). A multivariate boundary would follow the
coarea argument of the manuscript's localization appendix. It is not claimed here.

Numerical checks: `simulation/src/separation/run.py boundary` and `run.py coverage`.

---

## Setup

For each n, (H_i, A_i, Y_i), i ≤ n, are iid (a triangular array).

**Assumptions**
- (D1) H has density f with f ≤ f̄. f is continuous at 0 and f(0) > 0.
- (D2) |μ_a| ≤ B, and μ_a is L-Lipschitz on [−u_0, u_0], where μ_a(h) = E[Y | A = a, H = h].
- (D3) σ_a²(h) = Var(Y | A = a, H = h) is continuous at 0 with σ_1²(0) + σ_0²(0) > 0, and E[Y⁴ | A, H] ≤ M.
- (D4) A | H ~ Bern(p_n(H)) with p_n(h) = σ(h/τ_n), known, and τ_n → 0.
  This is a **single common temperature**.
- (D5) Greedy action a*(h) = 1{h > 0}. The operational gap is g(h) = φ(|h|) with φ(u) ≤ κ̄u on
  [0, u_0] and g ≤ Ḡ. Strong form: φ is differentiable at 0 with φ′(0) = κ > 0.

No structure extrapolating c = μ_1 − μ_0 away from 0 is assumed. None is needed on this side.

**Notation**
- ω = p(1 − p) = K(H/τ), where K(u) = eᵘ/(1 + eᵘ)².
- Kernel facts: ∫K = 1, ∫|u|K = 2 log 2, ∫K(u)σ(u)du = 1/2, K(u) ≤ e^{−|u|}, and
  ∫_0^∞ vσ(−v)dv = π²/12.
- m_τ = E[ω].

**Targets**
- β_0 = c(0): the fixed boundary effect.
- β_ov,τ = E[ωc(H)]/m_τ: the temperature-indexed overlap effect.

**Estimator.** With w_1 = A(1 − p) and w_0 = (1 − A)p,
- D̂_a = n⁻¹Σ w_{a,i},
- m̂_a = Σ w_{a,i}Y_i / Σ w_{a,i},
- β̂ = m̂_1 − m̂_0.

**Variance estimator.** ŝ² = n⁻¹Σ ψ̂_i², where ψ̂ = w_1(Y − m̂_1)/D̂_1 − w_0(Y − m̂_0)/D̂_0.
The Wald interval is β̂ ± z ŝ/√n.

---

## Theorem D

Assume (D1)–(D5), τ_n → 0 and nτ_n → ∞.

1. **Asymptotic normality.** (β̂ − β_ov,τ_n) / (s_n/√n) ⇒ N(0, 1), where
   s_n² = (σ_1²(0) + σ_0²(0)) / (2τ_n f(0)) · (1 + o(1)).
2. **Bias.** |β_ov,τ − β_0| ≤ C_1τ_n for n large.
3. **Variance estimation.** ŝ²/s_n² →_p 1.
4. **Validity at the fixed boundary effect.** If also nτ_n³ → 0, then (β̂ − β_0)/(ŝ/√n) ⇒ N(0, 1), so
   the Wald interval for β_0 has asymptotic coverage 1 − η.
5. **Exploration loss.**
   R_n = nE[g(H)σ(−|H|/τ_n)] ≤ nτ_n²κ̄f̄π²/6 + nḠe^{−u_0/τ_n}.
   In the strong form, R_n = nτ_n²κf(0)π²/6 · (1 + o(1)).

**Corollary.** For τ_n = n^{−ζ} with ζ ∈ (1/2, 1), the Wald interval for β_0 is asymptotically valid and
R_n → 0 simultaneously. The three rate conditions are:

| Condition | Needed for | Holds when |
|---|---|---|
| nτ_n → ∞ | effective boundary sample grows | ζ < 1 |
| nτ_n³ → 0 | bias negligible relative to SE | ζ > 1/3 |
| nτ_n² → 0 | cumulative exploration loss vanishes | ζ > 1/2 |

The achievable error is |β̂ − β_0| = O_p((nτ_n)^{−1/2}). Under the strong form of (D5),
R_n ≍ nτ_n², so this equals O_p((nR_n)^{−1/4}). **We state this as achievable, not optimal.** No lower
bound on the boundary side is claimed. The statement is in probability. An RMSE bound would need
separate L² control and a convention for the event Σw_a = 0; neither is given here.

**Remark (smoother effects).** Differentiability of c at 0 is **not** enough for O(τ²) bias.
Counterexample (audit; `run.py audit`): H ~ U(−1, 1) and c(h) = 1 + |h|^{3/2}, which is differentiable
at 0. Its bias is 2.034·τ^{3/2} (bias/τ² = 10.2, 14.4, 20.3, 28.8 at τ = .04, .02, .01, .005).

A sufficient condition for O(τ²) is: c is C^{1,1} near 0 (|c(h) − c(0) − c′(0)h| ≤ M_2h²), and f is
Lipschitz near 0 (constant L_f). Then:
- the linear term c′(0)E[ωH]/m_τ is O(τ²): E[ωH] = τ²∫uK(u)f(τu)du = τ²∫uK(u)(f(τu) − f(0))du
  (because ∫uK = 0), and |f(τu) − f(0)| ≤ L_fτ|u| near 0, so E[ωH] = O(τ³) up to an exponentially small
  tail;
- the remainder is ≤ M_2 f̄ τ³ ∫u²K / m_τ = O(τ²).

Under these conditions part 4 holds with nτ⁵ → 0. The main theorem uses only the O(τ) bound.

---

## Proof

### Step 1: kernel integrals
Let q be bounded and continuous at 0, and let j ≥ 1, r, s ≥ 0. Substitute h = τu. The integrand is
dominated by ‖q‖_∞ f̄ K(u), so dominated convergence gives

  E[K(H/τ)^j (1 − p)^r p^s q(H)] = τ f(0) q(0) ∫K^j σ(−u)^r σ(u)^s du + o(τ).

In particular m_τ = τf(0)(1 + o(1)).

### Step 2: bias
β_ov,τ − β_0 = E[ω(c(H) − c(0))]/m_τ. Split the expectation:
- On |H| ≤ u_0: |c(h) − c(0)| ≤ 2L|h|, and E[ω|H|] ≤ f̄τ²∫|u|K = 2 log 2 · f̄τ².
- On |H| > u_0: ω ≤ e^{−u_0/τ}, and |c(h) − c(0)| ≤ 4B.

So |β_ov,τ − β_0| ≤ (4 log 2 · Lf̄τ² + 4Be^{−u_0/τ})/m_τ = (4 log 2 · Lf̄/f(0))τ(1 + o(1)). ∎

### Step 3: linearization
Let μ̄_a = E[ωμ_a]/m_τ. Since E[w_1 | H] = ω and E[w_1Y | H] = ωμ_1, we have E[w_1(Y − μ̄_1)] = 0, and
likewise for arm 0. Define

  T_a = n⁻¹Σ_i w_{a,i}(Y_i − μ̄_a)/m_τ,  ξ_i = [w_{1,i}(Y_i − μ̄_1) − w_{0,i}(Y_i − μ̄_0)]/m_τ.

Then m̂_a − μ̄_a = (m_τ/D̂_a)T_a, μ̄_1 − μ̄_0 = β_ov,τ, and Eξ = 0.

Because w_1w_0 = 0 and E[A(1 − p)² | H] = ω(1 − p),

  s_n² := Var(ξ) = m_τ⁻² E[ω{(1 − p)(σ_1² + (μ_1 − μ̄_1)²) + p(σ_0² + (μ_0 − μ̄_0)²)}].

- μ̄_a → μ_a(0) by Step 1.
- (μ_a(τu) − μ̄_a)² → 0 pointwise and is bounded by 4B², so those terms are o(τ).
- The σ² terms equal τf(0)σ_a²(0)·½ + o(τ), since ∫K(1 − σ) = ∫Kσ = ½.

Therefore s_n² = (σ_1²(0) + σ_0²(0))/(2τf(0)) · (1 + o(1)).

Numerical check: σ_a² = 1 and f(0) = ½ give SD(β̂) ≈ (2/(nτ))^{1/2}.

### Step 4: central limit theorem
- *Fourth moment.* w_a⁴ ≤ w_a ≤ A(1 − p) or (1 − A)p, and E[(Y − μ̄_a)⁴ | A, H] ≤ C(M, B). So
  E[w_a⁴(Y − μ̄_a)⁴] ≤ Cm_τ, hence Eξ⁴ ≤ C′m_τ⁻³ ≍ τ⁻³.
- *Lyapunov.* n⁻¹Eξ⁴/s_n⁴ ≍ n⁻¹τ⁻³/τ⁻² = (nτ)⁻¹ → 0. Therefore n^{−1/2}Σξ_i/s_n ⇒ N(0, 1).
- *Denominators.* Var(D̂_a)/m_τ² ≤ E[w_a²]/(nm_τ²) ≤ 1/(nm_τ) → 0, so D̂_a/m_τ →_p 1 and
  m_τ/D̂_a − 1 = O_p((nτ)^{−1/2}).
- *Decomposition.* β̂ − β_ov = (T_1 − T_0) + (m_τ/D̂_1 − 1)T_1 − (m_τ/D̂_0 − 1)T_0.
  Here T_1 − T_0 = n⁻¹Σξ_i and T_a = O_p((nτ)^{−1/2}). The remainder is O_p((nτ)⁻¹) = o_p(s_n/√n).

Slutsky gives Part 1. ∎

### Step 5: variance estimator
Write ξ_{a,i} = w_{a,i}(Y_i − μ̄_a)/m_τ, so ξ = ξ_1 − ξ_0.

(i) n⁻¹Σξ_i²/s_n² has mean 1 and variance ≤ Eξ⁴/(ns_n⁴) ≍ (nτ)⁻¹ → 0. So it converges to 1 in
probability.

(ii) For arm 1, ψ̂_{1,i} − ξ_{1,i} = w_1(Y − μ̄_1)(1/D̂_1 − 1/m_τ) + w_1(μ̄_1 − m̂_1)/D̂_1. Hence
- n⁻¹Σ[w_1(Y − μ̄_1)(1/D̂_1 − 1/m_τ)]² = (n⁻¹Σξ_{1,i}²)(m_τ/D̂_1 − 1)² = O_p(s_n²) · O_p((nτ)⁻¹);
- n⁻¹Σ[w_1(μ̄_1 − m̂_1)/D̂_1]² ≤ (m̂_1 − μ̄_1)²/D̂_1, using w_1² ≤ w_1 so that n⁻¹Σw_1² ≤ D̂_1.
  By Step 4, (m̂_1 − μ̄_1)² = O_p((nτ)⁻¹), and D̂_1/m_τ →_p 1 with m_τ = τf(0)(1 + o(1)). Hence
  1/D̂_1 = (τf(0))⁻¹(1 + o_p(1)), and this term is O_p(1/(nτ²)). Relative to s_n² ≍ τ⁻¹ it is
  O_p((nτ)⁻¹).

Arm 0 is identical. By Minkowski, n⁻¹Σ(ψ̂_i − ξ_i)² = o_p(s_n²). Then Cauchy–Schwarz gives
|n⁻¹Σψ̂_i² − n⁻¹Σξ_i²| = o_p(s_n²), which is Part 3. ∎

### Step 6: centering at β_0
(β_ov,τ − β_0)/(s_n/√n) = O(τ) · O((nτ)^{1/2}) = O((nτ³)^{1/2}) → 0. Together with Parts 1 and 3 and
Slutsky, this gives Part 4. ∎

### Step 7: exploration loss
The off-greedy probability is σ(−|h|/τ).
- On |h| ≤ u_0: nE[φ(|H|)σ(−|H|/τ)] ≤ nκ̄f̄ · 2τ²∫_0^∞ vσ(−v)dv = nκ̄f̄τ²π²/6.
- Outside, the contribution is ≤ nḠe^{−u_0/τ}.
- Strong form: on v ≤ u_0/τ we have φ(τv)/τ → κv and f(τv) → f(0), dominated by κ̄vf̄σ(−v). Dominated
  convergence gives the asymptotic equality, and the region v > u_0/τ is the exponentially small term
  above.

For τ_n = n^{−ζ}: nτ² → 0 iff ζ > 1/2, and ne^{−u_0n^ζ} → 0 always. ∎

---

## Monte Carlo evidence

Sources: `simulation/results/separation/{coverage,boundary}.csv`, with seeds and settings in the `.json`
files. The earlier exploratory runs (seed 7 and 11) are summarized in the last table.
DGP: H ~ U(−1, 1), μ_0(h) = h, c(h) = 1 + |h|, σ = 1. The kink makes the bias O(τ).

**Coverage diagnostics, ζ = 0.6 (`run.py coverage`, seed 20260919)**

| n | reps | Kish eff. obs/arm | bias/SD | SD MC (pred.) | SE/SD | cover (SE hat) | cover (MC SD) | MCSE |
|---|---|---|---|---|---|---|---|---|
| 10⁴ | 3000 | 37 | 0.02 | .222 (.224) | .985 | .936 | .945 | .004 |
| 10⁵ | 3000 | 96 | 0.04 | .140 (.141) | 1.000 | .951 | .952 | .004 |
| 10⁶ | 1000 | 244 | 0.05 | .092 (.089) | .963 | .934 | .951 | .007 |

**Boundary runs (`run.py boundary`, seed 20260918)**

| ζ | n | nτ = n^{1−ζ} | reps | RMSE | coverage of β_0 | cum. loss (pred.) |
|---|---|---|---|---|---|---|
| .60 | 10⁴ | 39.8 | 1000 | .215 | .949 | .129 (.130) |
| .60 | 10⁵ | 100.0 | 500 | .142 | .944 | .082 (.082) |
| .60 | 10⁶ | 251.2 | 200 | .082 | .970 | .051 (.052) |
| .75 | 10⁴ | 10 | 1000 | .465 | **.887** | .0081 (.0082) |
| .75 | 10⁵ | 17.8 | 500 | .340 | .924 | .0026 (.0026) |
| .75 | 10⁶ | 31.6 | 200 | .247 | .930 | .0008 (.0008) |

**Earlier exploratory runs (`coverage_diag.py`, not in this runner)**

| n | reps | seed | SE/SD | coverage |
|---|---|---|---|---|
| 10⁴ | 3000 | 7 | .986 | .941 |
| 10⁵ | 3000 | 7 | .986 | .946 |
| 10⁶ | 1000 | 11 | .991 | .948 |

**Reading of the evidence (to be carried into the manuscript)**
- Exploration loss matches the Step 7 constant in every cell.
- Bias relative to SD is ≤ 5%, and skewness and kurtosis are near 0. No evidence of target drift or
  non-normality.
- **Finite-sample under-coverage is real when the effective boundary sample is small.** At ζ = .75 and
  n = 10⁴ (nτ = 10), coverage is .887 over 1,000 replications. It recovers slowly: .924 at nτ = 17.8
  and .930 at nτ = 31.6. Here nτ = n^{1−ζ} is a scale, not the Kish effective sample.
  Asymptotic validity (Part 4) says nothing about this regime, and the manuscript must state that.
- At ζ = .6 and n = 10⁶ the first two seeds gave .948 and .934 (1,000 replications each). A
  5,000-replication run (`results/separation/coverage_n1e6_r5000/`, seed 20260920) gives:
  - coverage .9458 (MCSE .0032, from the observed coverage) with the estimated SE;
  - coverage .9496 using the Monte Carlo SD of the same replications (column renamed from
    `coverage_true_sd` to `coverage_mc_sd`);
  - SE/SD .993, MC SD .0896 vs predicted .0892;
  - bias/SD < .001, skewness −.06, excess kurtosis −.05.

  The remaining shortfall (≈ .004, 1.3 MCSE) matches an SE about 0.7% below the sampling SD. The .934
  run is attributed to Monte Carlo variation.
- SE hat agrees with the asymptotic formula: mean SE ≈ .089 = predicted. The deviation of SE/SD from 1
  at n = 10⁶ comes from the Monte Carlo SD exceeding the formula (.092 vs .089). The relative MCSE of an
  SD estimate from 1,000 replications is about 2.2%.

**Draft wording:**
> With ζ = 0.6 (nτ between 40 and 250), coverage of the fixed boundary effect is between 0.936 and
> 0.970 across runs; the runs with at least 3,000 replications give 0.936–0.951. With ζ = 0.75 and nτ between 10 and 32, the Wald interval under-covers
> (0.887–0.930). At n = 10⁶ (nτ ≈ 250), 5,000 replications give 0.946
> (Monte Carlo SE 0.003); the small shortfall is consistent with the standard error being about 0.7%
> below the sampling SD. Lower values seen in smaller runs (0.91 with 200, 0.934 with 1,000
> replications) did not recur at this scale.

For reference, a bias-free normal estimator whose SE is 0.986 × SD has coverage
2Φ(1.96 × 0.986) − 1 ≈ 0.9467.
