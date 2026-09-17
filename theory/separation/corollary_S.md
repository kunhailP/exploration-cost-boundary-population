# Corollary S — the separation inside one model

> **Notation.** The manuscript writes e_τ = P(A = 1 | H) and q_τ = P(A ≠ a*(H) | H); the causal assumptions (C1)–(C3) are in its Section 2.

Status (2026-09-17, rev. 2): written in response to the audit of commit 6d111fc (Theorems A and D were
proved in different models). Revised after the audit of 6c096a0: γ is now an explicit essential-infimum
assumption (φ is not assumed monotone), the perturbation assumptions are collected in (S-A), the
temperature exponent is ζ (not a, which indexes arms) and the confidence level is 1 − η. Self-checked; **not reviewed by a domain
expert.**

Theorem A proves its impossibility in the large model 𝓜(B, σ) of all bounded measurable mean
functions. Theorem D assumes Lipschitz means near the boundary. A minimax lower bound over a large
model does not transfer automatically to a smaller one. This corollary restates both sides in the
**same** model, the one Theorem D uses.

---

## Common model 𝓜_D

One-dimensional score Δ(h) = h and greedy action a*(h) = 1{h > 0}.

- (D1) H has density f ≤ f̄, continuous at 0, with f(0) > 0.
- (G) Y | A = a, H = h, past ~ N(μ_a(h), σ²).
- (D2) |μ_a| ≤ B, and μ_a is L-Lipschitz on [−u_0, u_0].
- (D5′) g(h) = φ(|h|) ≤ Ḡ, with φ(u) ≤ κ̄u on [0, u_0]. φ need not be monotone.

Under (G), (D3) holds with σ_a² ≡ σ² and E[Y⁴ | A, H] ≤ C(B, σ). So every μ ∈ 𝓜_D satisfies the
hypotheses of Theorem D.

**Assumption (S-A).** There are u_0 < u_1 < u_2 such that
1. F([u_1, u_2]) > 0;
2. γ := ess inf_{h ∈ [u_1, u_2]} g(h) > 0 (essential infimum with respect to F);
3. B > 0 and L > 0.

Item 2 is an assumption, not a consequence of (D5′). For example, φ(u) = ue^{−10u} is positive but
decreasing on part of its range, so φ(u_1) need not bound g on [u_1, u_2]. In the simulations
g(h) = |h|, and γ = u_1.

Define

  b(h) = sin²(π(h − u_1)/(u_2 − u_1)) · 1{u_1 ≤ h ≤ u_2}.

Then 0 ≤ b ≤ 1, b is ℓ_b-Lipschitz on ℝ with ℓ_b = π/(u_2 − u_1), E[b(H)] > 0 (b > 0 on the open
interval and F([u_1, u_2]) > 0, with F having a density), and b ≡ 0 on (−∞, u_1], so in particular near
the boundary. Set t_max := min{B, L/ℓ_b} > 0. (Since b vanishes on [−u_0, u_0], L/ℓ_b is not needed for
(D2); it is kept so the argument also covers a globally L-Lipschitz model.)

**Perturbation family.** For t ∈ [−t_max, t_max] let μ^t_1 ≡ 0 and μ^t_0 = −t b. Then:
- μ^t ∈ 𝓜_D;
- c^t := μ^t_1 − μ^t_0 = t b;
- β_0(μ^t) = c^t(0) = 0 for every t, and μ^t = μ^0 on (−∞, u_1];
- θ(μ^t) = t E[b(H)].

On supp b ⊂ (0, ∞) the greedy action is 1, so every perturbed observation is an **off-greedy**
observation with H ∈ [u_1, u_2], where g ≥ γ F-almost surely.

---

## Corollary S

Assume (S-A). Consider any design in 𝓜_D (adaptive or not), with cumulative exploration loss R_n(μ). Let
N = Σ_i 1{A_i = 0, H_i ∈ [u_1, u_2]}.

**(S1) Boundary side, pointwise.** Take the non-adaptive common-temperature design τ_n = n^{−ζ},
ζ ∈ (1/2, 1). Then R_n → 0, and for **every** μ ∈ 𝓜_D the Wald interval of Theorem D for β_0(μ) has
asymptotic coverage 1 − η. This is Theorem D and its corollary. The guarantee is pointwise in μ; no
uniformity over 𝓜_D is claimed.

**(S2) Population side, pointwise impossibility at vanishing loss.** For every design and every t,

  TV(P^n_{μ^0}, P^n_{μ^t}) ≤ E_{μ^0}[N] ≤ R_n(μ^0)/γ.

Hence if R_n(μ^0) → 0 (in particular under the design of S1), then for every t_1 ∈ (0, t_max] and
every estimator θ̂_n,

  max_{t ∈ {0, t_1}} P_{μ^t}( |θ̂_n − θ(μ^t)| ≥ t_1E[b]/2 ) ≥ (1 − R_n(μ^0)/γ)/2 → 1/2.

So **no estimator is consistent for θ at both μ^0 and μ^{t_1}**, although both lie in 𝓜_D, coincide
on (−∞, u_1], and have the same boundary effect.

**(S3) Population side, uniform.** For every design and every estimator,

  sup_{μ∈𝓜_D} E_μ(θ̂_n − θ(μ))² ≥ (E[b])² / ( sup_{μ∈𝓜_D} R_n(μ)/(σ²γ) + π²/t_max² ).

In particular, uniform consistency for θ over 𝓜_D requires sup_{𝓜_D} R_n → ∞.

**What is not claimed in 𝓜_D.**
- The constant σ²(E√g)² of (A.3) is a lower bound proved in 𝓜(B, σ) only. Whether it survives in 𝓜_D is left
  open.
- No uniform boundary guarantee over 𝓜_D.
- No lower bound on the boundary side.

---

## Proofs

### (S2)
1. *Likelihoods agree off {N ≥ 1}.* Under any design the sample density factorizes as
   Π_i dF(H_i) κ_i(A_i | H_i, 𝓕_{i−1}) φ_σ(Y_i − μ_{A_i}(H_i)). μ^0 and μ^t differ only in μ_0 on
   [u_1, u_2]. So the two densities coincide at every sample point with N = 0.
2. *TV bound.* Let E = {N = 0}. Integrating equal densities over E gives P_0(E) = P_t(E). For any event
   A, |P_0(A) − P_t(A)| = |P_0(A ∩ E^c) − P_t(A ∩ E^c)| ≤ max{P_0(E^c), P_t(E^c)}. Since
   P_t(E^c) = 1 − P_t(E) = P_0(E^c), this gives TV ≤ P_0(N ≥ 1) ≤ E_0[N].
3. *Link to loss.* Each observation counted in N is off-greedy with H_i ∈ [u_1, u_2]. Since H_i ~ F
   and g ≥ γ F-a.s. on [u_1, u_2], g(H_i) ≥ γ almost surely for those observations. Hence
   γN ≤ Σ_i g(H_i)1{A_i = ā(H_i)} a.s., and taking expectations, E_0[N] ≤ R_n(μ^0)/γ.
4. *Le Cam.* |θ(μ^{t_1}) − θ(μ^0)| = t_1E[b]. Any θ̂ satisfies
   P_0(|θ̂ − θ_0| ≥ t_1E[b]/2) + P_{t_1}(|θ̂ − θ_{t_1}| ≥ t_1E[b]/2) ≥ 1 − TV. ∎

### (S3)
1. *Score and information.* The score of the one-parameter family is
   Σ_i 1{A_i = 0} (−b(H_i))(Y_i + t b(H_i))/σ². Its summands are martingale differences, so
   I(t) = E_t[Σ_i 1{A_i = 0} b(H_i)²]/σ² ≤ E_t[N]/σ² ≤ R_n(μ^t)/(γσ²).
   (This uses b² ≤ 1{H ∈ [u_1, u_2]} and step 3 above.)
2. *Prior.* Put λ(dt) = t_max⁻¹ cos²(πt/(2t_max)) dt on [−t_max, t_max]. Its Fisher information is
   π²/t_max².
3. *van Trees.* For ψ(t) = tE[b],
   sup MSE ≥ ∫E_t(θ̂ − ψ(t))²λ(dt) ≥ (E[b])² / (∫I(t)λ(dt) + π²/t_max²).
   Bounding ∫I dλ ≤ sup R_n/(σ²γ) finishes the proof. ∎

---

## Reading for the manuscript

In one model 𝓜_D, the same logged data from a common temperature τ_n = n^{−ζ}, ζ ∈ (1/2, 1):
- give asymptotically valid inference for the boundary effect at every μ, with R_n → 0;
- cannot distinguish two members of 𝓜_D that share the boundary effect and differ in the population
  effect by t_1E[b].

Consistency for the population effect over 𝓜_D requires divergent exploration loss. These are different
causal targets. Success on the boundary effect is **not** a substitute for the population effect; each
answers a different question.

**Identification caveat.** In a submodel where boundary data identify θ (e.g. constant effects), the
perturbation family is excluded, and S2–S3 do not apply.
