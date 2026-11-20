# Cross sections

## Conventions

- `ggxxxsec.m` (γγ → χχ) and `xxggxsec.m` (χχ → γγ) use the full BW propagator (s − mφ²)² + mφ²Γφ². `heavyScalar` = `lightScalar` here.
- `xxggxsec(..., parts)`: `'sm'` = γγ + meson pairs; `'had'` = meson pairs; `'all'` = 'sm' + φφ, A′A′, φ′φ′. Julia `sigma_xxhad`, `:sm`, `:had`.
- Detailed balance: σ(γγ → χχ)/σ(χχ → γγ) = 2(s − 4m²)/s (S_γγ = ½). Prefactors /4608 vs /9216 (SS/PS), /64π vs /128π (SP/PP).
- All families: σ(χχ̄ → γγ) = π s Γ_χχ Γ_γγ / (p_χ² D), σ(γγ → χχ̄) = 8π Γ_χχ Γ_γγ / D, D = (s − mφ²)² + mφ²Γφ², Γ_γγ = g_φγγ² m³/64π.
- heavyScalarPS: Γ_φ = D α² mφ³/(4608π³) + gx² √max(0, mφ²/4 − mχ²)/(4π).
- `ggxxrate` = n_eq²⟨σv⟩/4; `Oh2_FI` adds no further ½.

## Pair channels χχ̄ → φφ / jj / A′A′

- Gated on s > 4m_pair² only.
- Waypoints at 4mφ² and meson-pair thresholds (`svavg.m` `extra_wp`, `thermal_average.jl`).
- Identical particles (φφ, π⁰π⁰, ηη, A′A′) carry ½.
- Scalar and pseudoscalar φφ forms checked against a trace calculation.

## SP/PP pole couplings

- `sppp_mixing.m`: mixing masses and angles θ_j = m²_{jφ}/(mφ² − m_j²). Used by `sppp_couplings.m` and `phi_partial_widths.m`.
- m²_{η′φ} = (√6/3) Σg (B0 F + (L6 + L7) 32 B0² Σm/F).
- `[mj, gcx, ggx, Gam] = sppp_couplings(mphi, mx, gx, g, interaction)`.

| case | poles |
|---|---|
| generic g | {π⁰, φ′, η} |
| g(1) = g(2) = g(3) | {η′, φ′} |

Singlet case: m_{η′φ′} = B0 f (g_u + g_d + g_s)/√6; g_{η′χ} = −gx m_{η′φ′}/(2(mφ² − m_η′²)); g_{η′γγ} = α/F_η′ (m_η′³ form); g_{φ′γγ} = g_{η′γγ} m_{η′φ′}/(2(mφ² − m_η′²)); η′ width = `G_etap` + Γ(η′ → χχ).

SP BW: 1/(m_j² − s) → 1/(m_j² − s − i m_j Γ_j), j ∈ {π⁰, η, φ}.

## Resonances — `thermal_pole_integral.m`

Used by `svavg.m` and `ggxxrate.m`.

| pole width | treatment |
|---|---|
| w_p/s_p < 1e-6 | notch log(s_p) ± 10⁻³; add NWA real(σ(s_p)) w_p² kernel(s_p) π/w_p |
| otherwise | in continuum, waypoints at ±{1, 5, 50} FWHM |

NWA vs brute force: < 0.1% (axialVector π⁰, SP η).

| interaction | poles |
|---|---|
| heavyScalar, lightScalar, `*PS` | {φ} |
| SP, PP | {π⁰, φ′, η} |
| axialVector | {π⁰, η, η′} |

## Freeze-in numerics

- gx bracket: FI root < 1e-15 on resonance (mx < mφ/2); ~1e-26 near axialVector meson poles.
- `compare_limits_g.m` loop-branch floor = `log10(gqvec_fi(j)) - 3`.
- FI rate integrated from T0 to T_RH (0.01 GeV). Without the cutoff the BW σ(γγ → χχ̄) → const at s ≫ mφ² and the integral diverges logarithmically.

## Three-body final states

| variable | channels |
|---|---|
| `LOOP_3PI=1` | χχ → 3π, φ → 3π; axialVector χχ̄ → 3π, A′ → 3π |
| `LOOP_PIPIPHI=1` | χχ → ππφ; axialVector χχ̄ → π⁺π⁻A′ |

- Off by default. Off ⇒ results bit-identical to two-body code.
- Read from the environment of the MATLAB process and parpool workers. Julia: `src/models/threebody.jl`, read in `__init__`; `threebody_refresh!()` re-reads.
- Python uses them only in the CMB bound.
- Tests: `test_threepi.m`; Julia `validate_threebody.jl` (W, σ, Γ to 1e-15; ⟨σv⟩ 4e-9).

### 3π (`threepi_weight.m`)

- W(s) = Σ S ∫ dE₁dE₂ |A|², π⁺π⁻π⁰ (S = 1), 3π⁰ (S = 1/6).
- σ = g² W / (128π³ β |s − mφ² + i mφ Γφ|²), × (s − 4mχ²)/s for SP. Γ(φ → 3π) = W(mφ²)/(64π³ mφ).
- π⁰ + φ pieces: (B₀(g_u − g_d)/F) X/(s − m_π²), X = m₊₋² − m_π² (X → s for 3π⁰). LO mixing B₀F (NLO is 2.58× larger and breaks the cancellation).
- η pole added (η′ for g ∝ [1 1 1]).
- LO χPT.

### ππφ (`pipiphi_weight.m`)

- π⁺π⁻φ (S = 1), π⁰π⁰φ (S = ½); σ = g² W/(128π³β), × (s − 4mχ²)/s for SP. Added in `threepi_sigma.m` under `'had'`.
- σ ∝ g_χ² g_q⁴ (3π: g_χ² g_q²; t-channel meson pairs: g_χ⁴ g_q⁴). Threshold √s = 2m_π + mφ.
- Counted in `'sm'`; production by detailed balance with n_φ = n_φ^eq.

### axialVector 3π (`axial3pi_weight.m`)

- WT: transverse, J = 1, contact term, π⁺π⁻π⁰ only. WL: longitudinal, J = 0, π⁰ path × m_π²/(B₀Δg). No interference.
- σ = λ_a²/(32π³) [β WT/(6|D_A|²) + mχ² WL/(β s² m_A⁴)]; Γ(A′ → 3π) = WT(m_A²)/(192π³ m_A).
- Current scaled by m_{π⁰A′}/F (NLO, `corr4` ≈ 2, 4.1× LO in rate).
- Γ_A (χχ̄ + 3π) regulates the J = 1 A′ pole; `svavg` adds the pole and 3π thresholds except in `'gg'`.
- Tests: `test_threepi.m` checks 7–11.

### axialVector π⁺π⁻A′ (`pipiA_weight.m`)

- V_σ = −2iΔa² [ε*_σ − q_σ(q·ε*)/(s − m_π²)]: contact + π⁰ pole. No π⁰π⁰A′, η, η′ path.
- Added to 3π weights in `axial3pi_sigma.m`. Threshold √s = m_A + 2m_π±.
- Normalisation Δa_eff Δa, Δa_eff = m_{π⁰A′}/F, Δa = a_u − a_d; corr4² ≈ 4× LO.
- σ ∝ λ_a² g_q⁴, with an extra 1/m_A² from the longitudinal A′: σ_J=0 ∝ 1/m_A⁶. Saturates at ~1.5× the rest of ⟨σv⟩ at large g_q, small λ_a (underabundant region).
- Tests: `test_threepi.m` checks 12–13; Julia `pipiA_weight` (`WTa`, `WLa`).
