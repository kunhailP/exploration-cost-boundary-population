# Population side of the separation: Theorems A, B, Proposition C, Lemma R

> **Alignment with the manuscript (review revision merged 2026-09-17).** The manuscript states the
> following conditions explicitly; they are implicit here:
> - Proposition 1 (= A′): 0 < E√g < ∞, δ_n > 0, nδ_n² > 4C, q = 1/2 where g = 0, and E[g] < ∞ for
>   the uniform-mixing loss.
> - Theorem 4 (= B): g = φ(|Δ|) and E[g] < ∞. The essential supremum of |Δ| is written **D_max**
>   (G is the generation kernel).
> - Lemma 1 (= R): F(g > 0) = 1 and 0 < a ≤ ĝ/g ≤ b < ∞. Capped designs, including a positive mass of
>   g = 0, are treated in the limit E[1/q] → ∞.
> - "Design variance" σ²E[1/q]/n is called the **sparse-exploration criterion** V_sp. It is an
>   allocation criterion, not the estimator MSE.
>
> **Notation.** In this file p(h) is the **off-greedy** probability P(A ≠ a*(h) | H = h). The manuscript writes it as q(h) and uses e(h) = P(A = 1 | H = h) for the probability of strategy 1. The causal assumptions (consistency, iid units without interference, known sequential assignment) are stated in Section 2 of the manuscript.

Status (2026-09-17, rev. 2): proofs written out in full and self-checked. One audit round (review of
commit 6d111fc) has been incorporated: scope of B, endpoints of C, the "iff" in R, and wording of A.2
and of the extrapolation remark. **Not yet reviewed by a domain expert.**
How Theorem A connects to the boundary model of Theorem D is in `corollary_S.md`.
Numerical checks are in `simulation/src/separation/run.py` (subcommands named below).
Scope statements are part of the results and must be carried into the manuscript verbatim.

---

## 0. Setup

**Data.** For i = 1, …, n: context H_i ~ F on a measurable space 𝓗, iid; action A_i ∈ {0,1}; outcome Y_i.

**Deployer.** A measurable greedy action a*(h) ∈ {0,1} and a measurable, bounded *operational* gap
g(h) = V(a*(h), h) − V(1 − a*(h), h) ≥ 0, where V is the deployer's own objective (e.g. a reward-model
score). g is a design input. It is **not** a functional of the human outcome Y. Write ā(h) = 1 − a*(h).

**Designs.** A design is a sequence of assignment kernels
P(A_i = ā(H_i) | H_i, 𝓕_{i−1}) = P_i,  where 𝓕_{i−1} = σ(H_j, A_j, Y_j : j < i).
It may be adaptive. It is *non-adaptive* if P_i = p(H_i) for a fixed p : 𝓗 → [0, 1/2].
Throughout, the treatment is the **strategy choice** A; the generation kernel producing text under
each strategy (including token-level temperature) is held fixed.

**Outcome model** 𝓜(B, σ): Y | A = a, H = h, 𝓕_{i−1} ~ N(μ_a(h), σ²), where μ_0, μ_1 are *any*
measurable functions with |μ_a| ≤ B. No smoothness, no parametric form, and no restriction linking
μ_a(h) across h. **In particular 𝓜 contains no structure that extrapolates effects learned near the
decision boundary to other contexts.**

**Estimand.** θ(μ) = ∫ (μ_1 − μ_0) dF.

**Cumulative exploration loss.** R_n(μ) = E_μ[ Σ_i g(H_i) 1{A_i = ā(H_i)} ].
For non-adaptive designs R_n = n E[g(H) p(H)], which does not depend on μ.

---

## Theorem A (information–exploration inequality)

Let {S_1, …, S_K} be disjoint measurable sets with f_k = F(S_k) > 0, each contained in {a* = 0} or in
{a* = 1}, and let g_k = inf_{h ∈ S_k} g(h). For **every** design (adaptive or not) and **every**
estimator θ̂,

  sup_{μ∈𝓜} E_μ(θ̂ − θ(μ))² · ( sup_{μ∈𝓜} R_n(μ) + (σ²π²/B²) Σ_k g_k ) ≥ σ² ( Σ_k f_k √g_k )².   (A.1)

**Corollaries.**
- (A.2) *No free lunch.* If F(g ≥ γ) > 0 for some γ > 0, then sup_𝓜 MSE → 0 implies sup_𝓜 R_n → ∞.
- (A.3) *Sharp constant from below.* If sup_𝓜 MSE ≤ δ_n² with δ_n → 0, then
  liminf_n δ_n² sup_𝓜 R_n ≥ σ² (E √g(H))².
- (A.4) (A.1) holds for every R ≥ 0. At R = 0 it gives the finite floor
  sup MSE ≥ B²(Σ f_k √g_k)²/(π² Σ g_k). The display "RMSE ≳ σ E√g / √R" is (A.1) in the regime R ≫ C_K.
  **Neither side contains n.**

### Proof of (A.1)

1. *Submodel.* For ϑ = (ϑ_1, …, ϑ_K) ∈ [−B, B]^K define μ^ϑ by
   - μ^ϑ_{a*(h)}(h) = 0 for all h;
   - on S_k, μ^ϑ_{ā(h)}(h) = s(h)ϑ_k with s(h) = +1 if ā(h) = 1 and s(h) = −1 if ā(h) = 0;
   - μ^ϑ_{ā(h)}(h) = 0 off ∪S_k.

   Since S_k lies on one side of the greedy rule, s is constant on S_k. Then μ^ϑ ∈ 𝓜 and
   θ(μ^ϑ) = Σ_k f_k ϑ_k, because on S_k the effect μ_1 − μ_0 equals +ϑ_k whichever arm is off.
2. *Likelihood.* The joint density of the sample factorizes as
   Π_i dF(H_i) · κ_i(A_i | H_i, 𝓕_{i−1}) · φ_σ(Y_i − μ^ϑ_{A_i}(H_i)).
   The assignment kernels κ_i are known functions of observed data and do not involve ϑ. Hence the
   score for ϑ_k is Σ_i 1{H_i ∈ S_k, A_i = ā(H_i)} s(H_i)(Y_i − s(H_i)ϑ_k)/σ². Its summands are
   martingale differences with respect to (𝓕_i), so cross-products have mean zero and the Fisher
   information matrix is diagonal with I_k(ϑ) = E_ϑ[N_k]/σ², where N_k = Σ_i 1{H_i ∈ S_k, A_i = ā(H_i)}.
3. *Prior.* Let λ = ⊗_k λ_B with λ_B(dx) = B⁻¹ cos²(πx/(2B)) dx on [−B, B]. This density vanishes
   at ±B, and its Fisher information is J = π²/B².
4. *Multivariate van Trees* (Gill & Levit 1995) for the linear functional ψ(ϑ) = Σ f_k ϑ_k with
   N̄_k := ∫ E_ϑ[N_k] λ(dϑ):

   sup_𝓜 MSE ≥ ∫ E_ϑ(θ̂ − ψ(ϑ))² λ(dϑ) ≥ Σ_k f_k² / (N̄_k/σ² + J) = Σ_k f_k σ² / a_k,
   where a_k := (N̄_k + σ²J)/f_k.

   (The middle bound is the van Trees inequality maximized over the direction vector.)
5. *Cauchy–Schwarz.* σ Σ_k f_k √g_k = Σ_k √(f_k g_k a_k) · σ√(f_k/a_k), so
   σ²(Σ_k f_k √g_k)² ≤ (Σ_k f_k g_k a_k)(Σ_k f_k σ²/a_k).
6. *Linking to loss.* Σ_k f_k g_k a_k = Σ_k g_k N̄_k + σ²J Σ_k g_k. Since g ≥ g_k on S_k,
   Σ_k g_k N_k ≤ Σ_i g(H_i) 1{A_i = ā(H_i)}. Therefore
   Σ_k g_k N̄_k ≤ ∫ R_n(μ^ϑ) λ(dϑ) ≤ sup_𝓜 R_n.
7. Combining 4–6 gives (A.1). ∎

### Proof of (A.2)–(A.3)

- (A.2): take K = 1 and S_1 = {g ≥ γ} ∩ {a* = j}, where j is a side with F(S_1) = f_1 > 0, and write
  g_1 = inf_{S_1} g (so g_1 ≥ γ). (A.1) gives
  sup R_n ≥ σ² f_1² g_1 / sup MSE − σ²π² g_1/B² ≥ g_1(σ² f_1²/sup MSE − σ²π²/B²).
  The right side → ∞ as sup MSE → 0.
- (A.3): fix K and take level-set cells S_{k,j} = {g ∈ [γ_k, γ_{k+1})} ∩ {a* = j} with positive mass.
  (A.1) gives δ_n² sup R_n ≥ σ²(Σ f √g_k)² − δ_n² C_K, so liminf δ_n² sup R_n ≥ σ²(Σ f_k √g_k)².
  The lower Lebesgue sums Σ F(S) inf_S √g approach E√g as the level mesh → 0 (monotone convergence;
  no continuity of g is needed). Take the supremum over K. ∎

### Scope remarks (to appear in the paper)

1. **Extrapolation structure removes the separation.** (A.2) can fail in any model where the data near
   the decision boundary identify θ. Examples: a constant effect μ_1 − μ_0 ≡ θ, or a parametric model
   whose parameters, and hence θ, are identified from boundary data. Being parametric is not by itself
   enough; what matters is that identification holds. Under intermediate
   restrictions (e.g. global Lipschitz effects) the perturbations in step 1 must respect the restriction,
   and the constants change; compare Mou, Ding, Wainwright & Bartlett (2023) and Khan, Saveski &
   Ugander (2024). We claim (A.1) only for 𝓜(B, σ).
2. **Designs.** (A.1) covers adaptive and sequential designs, with the worst-case expected loss.
3. **The bound is information-theoretic.** It applies to all estimators, not to IPW only.

---

## Theorem A′ (achievability; non-adaptive)

Assume E√g(H) > 0 and nδ_n² → ∞. Define

  p_c(h) = min{1/2, c/√g(h)}  (with p_c = 1/2 where g = 0),  c_n = E√g / (nδ_n²/(σ² + B²) − 4),

and use the Horvitz–Thompson estimator with the known assignment probabilities. Then
sup_𝓜 MSE ≤ δ_n² and

  R_n ≤ (σ² + B²) (E√g)² / δ_n² · (1 + o(1)).

With C = σ² + B² and nδ_n² ≥ 4C, the uniform design p ≡ C/(nδ_n² − 2C) also has sup_𝓜 MSE ≤ δ_n²,
and R_n = C E[g] / δ_n² (1 + o(1)).
(Correction, 2026-09-17 review of 899e1cc: p = C/(nδ_n²) is not enough. At μ_0 = μ_1 ≡ B the HT MSE
equals C/(np(1 − p)) > δ_n², e.g. 0.010204 > 0.01 for B = σ = 1, n = 10⁴, δ = .1.)

**Proof.**
1. HT is unbiased with variance ≤ n⁻¹(σ² + B²) E[1/e_1 + 1/(1 − e_1)] ≤ n⁻¹(σ² + B²)(E[1/p_c] + 2),
   where e_1 = P(A = 1 | H).
2. Pointwise 1/p_c = max{2, √g/c} ≤ √g/c + 2, so E[1/p_c] ≤ E√g/c + 2.
   The choice of c_n gives E[1/p_c] + 2 ≤ nδ_n²/(σ² + B²), hence variance ≤ δ_n².
3. R_n = n E[g p_c] ≤ n c_n E√g = (E√g)² n / (nδ_n²/(σ² + B²) − 4). ∎

**Sandwich.** σ²(E√g)² ≤ liminf δ_n² R*_n ≤ limsup δ_n² R*_n ≤ (σ² + B²)(E√g)².
- The order 1/δ² and the functional (E√g)² appear in both bounds.
- The constants differ: the lower bound has σ², and the HT upper bound retains B². **The optimal
  constant is not established.**
- With finitely many contexts, a stratified difference-in-means estimator should reduce the constant
  toward σ². This is not proved here.

---

## Theorem B (cost of a single common temperature)

**Design class.** a*(h) = 1{Δ(h) > 0} and p_τ(h) = σ(−|Δ(h)|/τ) with a **single finite τ > 0 common
to all contexts** (τ may depend on n); non-adaptive. The operational gap is g = φ(|Δ|).

**Assumptions.**
- (B1) |Δ(H)| has a density ≥ f_min > 0 on [0, u_0].
- (B2) φ(u) ≥ κu on [0, u_0].
- (B3) G := ess sup |Δ(H)| < ∞, and q_η := max_j P(|Δ(H)| ≥ G − η, a*(H) = j) > 0 for every η ∈ (0, G).

**Claim.** Fix η ∈ (0, G). Let θ̂_n be any estimator with sup_𝓜 MSE ≤ δ_n² under p_{τ_n}, where
δ_n² ≤ q_η²B²/(2π²) and 2nδ_n² > q_ησ². δ_n may be fixed. Then

  R_n ≥ κ f_min c_0 · n · min{ u_0², (G − η)² / log²( 2nδ_n² / (q_η σ²) ) },  c_0 = ∫_0^1 v σ(−v) dv > 0.

Consequently R_n / [(σ² + B²)(E√g)²/δ_n²] ≳ nδ_n² / log²(nδ_n²) → ∞ whenever nδ_n² → ∞.

**Proof.**
1. *Temperature must be large.* Apply the proof of (A.1) with the single cell S = {|Δ| ≥ G − η,
   a* = j}, where j attains q_η. There p_τ ≤ e^{−(G−η)/τ}, so N̄ ≤ n q_η e^{−(G−η)/τ}, and step 4 gives
   δ_n² ≥ q_η² σ² / (n q_η e^{−(G−η)/τ} + σ²π²/B²).
   For n large (δ_n² ≤ q_η² B²/(2π²)) this forces n q_η e^{−(G−η)/τ} ≥ q_η² σ²/(2δ_n²), i.e.
   τ_n ≥ (G − η) / log(2nδ_n²/(q_η σ²)).
2. *Large temperature is expensive near the boundary.*
   R_n = n E[φ(|Δ|) σ(−|Δ|/τ)] ≥ nκ f_min ∫_0^{u_0} u σ(−u/τ) du = nκ f_min τ² ∫_0^{u_0/τ} v σ(−v) dv.
   The right side is ≥ nκ f_min c_0 min{τ, u_0}².
3. Combine 1 and 2. ∎

**Finite-context version.** If |Δ| takes finitely many values, with 0 < Δ_min < Δ_max, φ(Δ_min) > 0 and
masses bounded below, then for fixed δ small enough that step 1 applies (δ² ≤ q²B²/(2π²), where q is
the mass of the Δ_max cell), R_n ≥ c n^{1 − Δ_min/Δ_max}. For large δ a constant estimator can meet
the precision target with no exploration, and the statement is vacuous. The proof is the same:
τ ≥ Δ_max / log(C nδ²), and the Δ_min cell contributes n f φ(Δ_min) e^{−Δ_min/τ}/2.

**Scope.** Theorem B concerns **exactly one deterministic common temperature per sample size n**. It
does **not** cover:
- context-dependent temperatures (Proposition C);
- random mixtures of finite temperatures, even with bounded support and chosen independently of h.

**Counterexample for mixtures (audit, 2026-09-17; `run.py audit`).** Take H ~ U(−1, 1), u = |H|,
g = u, and q(u) = σ(−u) (temperature 1). Let p_n = w_n q + (1 − w_n) σ(−nu), a mixture of temperature
1 with weight w_n and temperature 1/n.
- With w_n = E[1/q]/(nδ²), p_n ≥ w_n q gives E[1/p_n]/n ≤ δ².
- R_n ≤ E[1/q] E[uq]/δ² + O(1/n), a constant in n.
- Numerically, at δ = .1: R_n = 46.362 at n = 10⁴, 10⁵, 10⁶.
- The same construction with w_n = E[1/q]/(nδ_n²/(σ² + B²) − 2) gives HT risk ≤ δ_n² uniformly over 𝓜.

So a small, shrinking weight on one moderate temperature escapes the n/log² n cost, and no infinite
temperature is needed. Such mixtures are a legitimate implementation of polynomial exploration. The
lesson of B is about a single common temperature.

---

## Proposition C (context-dependent temperature implements any design)

1. For Δ(h) ≠ 0 and p(h) ∈ (0, 1/2), the finite temperature τ(h) = |Δ(h)| / log((1 − p(h))/p(h)) gives
   σ(−|Δ(h)|/τ(h)) = p(h).
   - p(h) = 1/2 with Δ(h) ≠ 0 needs the extended value τ(h) = ∞.
   - Where Δ(h) = 0, every τ gives p = 1/2, so only p = 1/2 is implementable there.

   Hence context temperatures in (0, ∞] implement exactly the non-adaptive designs with p ∈ (0, 1/2]
   and p = 1/2 on {Δ = 0}. The capped design p* of Theorem A′ satisfies this, because g = φ(|Δ|) = 0 on
   {Δ = 0} gives p* = 1/2 there.
2. For the design of Theorem A′, p*(h) = c_n/√φ(|Δ(h)|) wherever the cap is slack, and
   τ*(h) = |Δ(h)| / ( log(1/c_n) + ½ log φ(|Δ(h)|) + log(1 − p*(h)) ).
   Since c_n → 0, τ*(h) log(1/c_n) → |Δ(h)| for each h with Δ(h) ≠ 0.
   The limit is pointwise for fixed h with g(h) > 0. It does **not** say that the implementing temperature
   is small near the boundary at finite n: where the cap binds (e.g. 0 < |h| ≤ 4c_n² when g = |h|), the
   implementing temperature is ∞. (The earlier slogan "cold at the boundary, hot where the policy is sure"
   was withdrawn after the 899e1cc review.) ∎

---

## Lemma R (misspecified gaps)

Allocate with ĝ instead of g, normalized to the same design variance (cap slack):
p = c/√ĝ with E[1/p] = T. Suppose ĝ/g ∈ [a, b] on {g > 0}, and set ρ = b/a. Then

  R(ĝ) / R(g) ≤ K(ρ) := (1 + √ρ)² / (4√ρ),

The bound is attained whenever the weighted law w below can be split into two parts of mass exactly
1/2, for example when w is non-atomic.

**Proof.**
- R(ĝ) = n E[g p] = n E[√ĝ] E[g/√ĝ] / T, and R(g) = n (E√g)² / T.
- With s = (ĝ/g)^{1/2} ∈ [√a, √b] and the probability measure w ∝ √g dF,
  R(ĝ)/R(g) = E_w[s] E_w[1/s].
- By Kantorovich's inequality this is ≤ (√a + √b)²/(4√(ab)) = K(b/a).
- If w admits a set of mass exactly 1/2, putting s = √a on it and s = √b elsewhere attains the bound. ∎

**Consequences.**
- A common rescaling ĝ = c·g changes nothing: only the **relative** misspecification ρ = b/a across
  contexts matters.
- **Sufficient condition.** If E[g]/(E√g)² > K(ρ), gap-based allocation beats uniform exploration for
  every admissible ĝ.
- **Necessary and sufficient** only when the bound is attained, e.g. when w ∝ √g dF is non-atomic.
  With atoms it can fail. Example (`run.py audit`): F = (.99, .01), g = (1, 100), ρ = 36. Then
  uniform/optimal = 1.675 < K(36) = 2.042, yet the worst admissible misspecification gives only 1.347,
  so gap-based allocation still beats uniform.
- Example (H ~ U(−1, 1), g = |H|; w is non-atomic): E[g]/(E√g)² = 9/8 = K(4). So with ρ ≥ 4, uniform
  exploration is at least as good in the worst case. The worst case is over allocations that misstate the gap by a factor a
  on {|H| < cut} and b elsewhere. **This conclusion is specific to the example.**

---

## Numerical checks

| Claim | Subcommand | Reference result |
|---|---|---|
| (A.1) close to the asymptotic constant in finite contexts | `run.py finite` | bound formula 51.887 vs **asymptotic constant** σ²(Σf√g)²/δ² = 52.361 (f = (.5, .5), g = (.2, 1), σ = 1, B = 5, δ = .1). Both are closed-form evaluations; no finite-sample MSE is measured |
| Design-variance costs (sparse-exploration criterion E[σ²/p]/n = δ²) | `run.py costs` | closed-form / quadrature, not estimator runs |
| Estimator runs (AIPW, Wald coverage), common random numbers across designs | `run.py designs` | Monte Carlo |
| Audit counterexamples (mixture temperatures, R without atomlessness, smoothness remark of D) | `run.py audit` | quadrature |
| Theorem B rates | `run.py costs` (continuous), `run.py finite` (slope → 1 − Δ_min/Δ_max) | |
| Lemma R attained | `run.py kantorovich` | |
