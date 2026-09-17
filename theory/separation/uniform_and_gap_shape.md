# Uniform boundary coverage and the shape of the operational gap

Plain-text companion to Appendix A.2 (Proposition 1) and Appendix C (S1) of `paper/tex/main_jci.tex`.
Notation: ω = e(1−e) = Ω(H/τ), m_τ = E ω, ξ = {w_1(Y − μ̄_1) − w_0(Y − μ̄_0)}/m_τ, s_n² = Var ξ.

## Model 𝓜_Lip

Y | A = a, H = h ~ N(μ_a(h), σ²), |μ_a| ≤ B, μ_a L-Lipschitz on ℝ. F satisfies (D1), g satisfies (D5).

## (S1) uniformly over 𝓜_Lip

All constants C depend only on (B, L, σ, f̄). The law of (H_i, A_i) does not depend on μ, so m_τ,
D̂_a and the zero-weight event do not either.

1. **Variance sandwich.** Since (1−e) + e = 1 and |μ_a − μ̄_a| ≤ 2B:
   σ²/m_τ ≤ s_n² ≤ (σ² + 4B²)/m_τ.
2. **Berry–Esseen.** E|ξ|³ ≤ C m_τ^{−2} (w_a³ ≤ w_a, E w_a = m_τ, Gaussian third moment with mean ≤ 2B).
   Kolmogorov distance of Σξ_i/(√n s_n) to N(0,1) ≤ C E|ξ|³/(√n s_n³) ≤ C/(σ³ √(n m_τ)).
3. **Remainder.** Var T_a ≤ (σ² + 4B²)/(n m_τ); Var(D̂_a/m_τ) ≤ 1/(n m_τ). In units of s_n/√n the
   remainder (m_τ/D̂_a − 1)T_a is O_p((n m_τ)^{−1/2}) uniformly.
4. **Variance estimator.** Step 5 of Theorem 1 uses only E ξ⁴ ≤ C m_τ^{−3}, s_n² ≥ σ²/m_τ and the bounds
   in 3, so ŝ²/s_n² → 1 uniformly in probability.
5. **Bias.** |β_ov − β_0| ≤ 2L E[ω|H|]/m_τ ≤ 4 log 2 · L f̄ τ²/m_τ (no tail term under global
   Lipschitz). In units of s_n/√n: ≤ C τ² √(n/m_τ) = O((nτ³)^{1/2}).
6. **Coverage.** Combine 2–5: the coverage error is at most ε + C/√(n m_τ) + 6ε·max φ for all μ and
   large n; let ε ↓ 0.

## Proposition (shape of the operational gap)

(a) φ(u) ≤ κ̄ u^p on [0, u_0] ⇒ R_n ≤ 2 κ̄ f̄ C_p n τ^{p+1} + n Ḡ e^{−u_0/τ},
    C_p = ∫_0^∞ v^p sig(−v) dv = Γ(p+1)(1 − 2^{−p}) ζ_R(p+1)
    (C_{1/2} ≈ 0.678, C_1 = π²/12, C_2 ≈ 1.803; checked numerically).
    τ_n = n^{−ζ} with max{1/3, 1/(p+1)} < ζ < 1: valid Wald interval and R_n → 0.

(b) φ(u) ≥ λ_0 > 0 on (0, u_0] ⇒ R_n ≥ 2 log 2 · λ_0 f(0) n τ (1 + o(1)) → ∞ whenever nτ → ∞
    (∫_0^∞ sig(−v) dv = log 2, dominated convergence).

With a fixed cost per deviation, boundary inference under a common-temperature design therefore also
incurs divergent loss, and the separation (vanishing vs divergent loss) needs a cost that vanishes at the
boundary. This does not compare the minimal costs of the two targets over all designs.
