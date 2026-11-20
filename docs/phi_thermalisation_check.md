# φ–bath conversion: π±γ → π±φ and Primakoff

| channel | families | files |
|---|---|---|
| π±γ → π±φ | S⊗S, P⊗S; mφ < 2m_π | `pigamma_phi_xsec.m`, `pigamma_phi_rate.m`; Julia `src/pigamma_phi.jl` |
| e±γ → e±φ, e⁺e⁻ → γφ | all but axialVector | `primakoff_phi_xsec.m`, `primakoff_phi_rate.m`; Julia `src/primakoff_phi.jl` |

Used in:
- `neff_simple` (MATLAB, Julia `src/neff_simple.jl`), on by default. Off: `opts.pigamma = false`, `opts.primakoff = false` (Julia `pigamma=false`, `primakoff=false`).
- `freeze_coupled` as Γ_sc. Off: `LOOP_PHI_SC=0`.

Both rates are independent of mφ (γγ → φ goes as mφ³).

## π±γ → π±φ

Vertex L ⊃ −μ φ π⁺π⁻ − (μ/2) φ π⁰π⁰, μ = B₀(g_u + g_d). a = s − m², b = u − m²:

    M = e μ [ p·ε/(p·k) − p'·ε/(p'·k) ]
    Σ|M|² = −4e²μ² [ m²(a+b)²/(ab)² + (a+b)/(ab) − m_φ²/(ab) ]
    σ = −(e²μ²/8πa²) { (b₊−b₋)(m²/a² + 1/a) + m²(1/b₋ − 1/b₊) + ((2m² − m_φ²)/a + 1) ln(b₊/b₋) }

- s/s_th − 1 < 10⁻²: Σ|M|² integrated numerically.
- R = 4T/(32π⁴) ∫ ds σ (s − m²)² K₁(√s/T)/√s; Γ = R/n_φ^eq.

Tests (`test_pigamma_phi.m`): σ vs explicit 4-vector integration (10⁻⁹), threshold continuity, g_s independence, (g_u + g_d)² scaling, massless rate.

## Primakoff

g_φγγ² = 64π Γ_γγ/mφ³ (at m_ref = 10⁻¹² GeV if mφ = 0). Γ_P ≈ 3.4×10⁻³ g_φγγ² T³; Γ_P/H(T_ν,dec) ≈ 1.4 (g_q/10⁻⁵)² for ĝ = [2 −1 −1].

Massless e±, X = s − mφ², x = −t ∈ [0, X], D = m_D²:

    Σ|M|²(eγ → eφ) = e² g² x (s² + (X − x)²)/(x + D)²          (scalar = pseudoscalar)
    σ_C = α g² I/(16 s²)
    I = [ y²/2 − (2B + D) y + (C + 2BD) ln y + DC/y ]_{y=D}^{D+X},  B = X + D,  C = s² + B²
      → (αg²/8)[ln(s/D) − 7/4]   (s ≫ D)
    σ_A(e⁺e⁻ → γφ) = α g² (1 − m_φ²/s)³/24
    R = T/(32π⁴) ∫_{m_φ²} ds [8 f_e σ_C + 4 f_e² σ_A] s^{3/2} K₁(√s/T),   Γ = R/n_φ^eq

- X/D < 10⁻²: three-term series (continuous to 1.4×10⁻⁶).
- f_e = z²K₂(z)/2, z = m_e/T; m_D² = (4πα/3) T² f_e.
- ∝ g² exactly; unit-coupling shape cached per (mφ, T grid).

Tests: `test_primakoff_phi.m` (10 checks); Julia `validate_simple` (max |ΔNeff| diff 5.8×10⁻⁵), `validate_fc_sc` (max relerr 1.3×10⁻³).

## `neff_simple`: decoupling and late release

- φ equilibrated before T_ν,dec but decoupled by then: free-streaming share `w_free = weq·(1 − wct)`, ρ_φ^eq(T_ν,dec). ΔNeff = 4/7 if relativistic at recombination; < 0 if φ → γγ before recombination.
- Late release: φ held in the bath by Primakoff leaves at T_late, where τ(T) = ∫_{T_rec}^{T} (Γ/H) d ln T = 1 (log-interpolated). T_late = T_ν,dec if τ < 1 throughout; T_rec if τ ≥ 1 at recombination.
- Released only if ρ_eq(mφ, T_late) > 10⁻⁶ ρ_eq(0, T_late).
- Entropy loop runs to T_late/1.5 with T_late as a grid node; bath gives up `w_bath·gsp`.
- vs instant-release calculation at mφ = 10⁻¹² GeV: ≤ 0.05 (≤ 0.004 near the exclusion threshold).
- Outputs: `info.T_late`, `info.late_release`, `info.Gamma_prim_over_H`.

## Results: heavyScalar, g = g_q·[2 −1 −1], T_RH = 10 MeV, mχ = 10 MeV, g_χ = 1

Thermalisation coupling g_q^th = N^{−1/2}, N = ∫_{2 MeV}^{T_RH} (Γ/H) d ln T.

- g_q^th(πγ) ≈ 2.3×10⁻⁶ for mφ ≲ 10 MeV; 1.6×10⁻⁵ at 270 MeV.
- g_q^th(γγ) ∝ mφ^{−3/2}; equals g_q^th(πγ) at mφ ≈ 6 MeV.
- Window [g_q^th(πγ), g_q^th(γγ)] for mφ < 6 MeV is fully excluded: neutron star + stellar below ~0.5 MeV; SN1987A + K → πφ for 0.5–6 MeV (narrowest margin 1–6 MeV).
- Rate ∝ (g_u + g_d)².

ΔNeff at mφ = 1 eV (with Primakoff):

| g_q | 10⁻⁵ | 3.2×10⁻⁵ | 10⁻⁴ | ≥ 3.2×10⁻⁴ |
|---|---|---|---|---|
| ΔNeff | −0.075 | +0.026 | +0.68 | ≈ +1.3 (T_late ≈ 60–140 keV) |

- Without Primakoff: saturates at +0.17.
- mφ = 0.1 MeV, g_q = 10⁻⁵: −0.66 → −0.79. mφ = 3 MeV: negligible.
- `freeze_coupled`: shifts FIMP loci where gx ≲ 10⁻⁷ and mφ ~ 10 mχ ≲ 10 MeV.
