# axialVector

Spin-1 mediator A′. Arguments: `mphi` = m_A′, `gx` = λ_a, `g` = [a_u, a_d, a_s].

## Couplings — `axial_couplings.m`

`[mj, lam_chi, lam_gg, Gam_j] = axial_couplings(mphi, mx, gx, g)`. Used by `ggxxxsec.m`, `xxggxsec.m`, `make_xxsigma.m`, `ggxxrate.m`.

- γγ ↔ χχ goes through SM-meson poles {π⁰, η, η′}; there is no A′ pole.
- λ_jχ = λ_a m_jA′/m_A′² (m_jA′: meson–A′ mixing masses, `axial_mixing.m`).
- Gam_j = SM total width + Γ(j → χχ̄).
- η′: λ_η′γγ = √(64π³ Γ_η′γγ / m_η′³), Γ_η′γγ = 4.34e-6 GeV, Γ_η′ = 1.88e-4 GeV.
- Γ(A′ → χχ̄) = λ_a² m_A′ β³/12π, in `phi_partial_widths.m`. Used for χχ → A′A′ bookkeeping only, not as a propagator width.

## NLO mixing

`corr4` in `axial_mixing.m`: mixing masses carry 1 + 16 B0 L4 Σm_q/F².

## A′ stability

Below min(2m_χ, 3m_π) every A′ channel is closed (γγ: Landau–Yang; ππ: parity; π⁰γ, 3γ: C; meson mixing ∝ k·ε = 0). Then:
- Γ_total = 0; `freeze_coupled` flags A′ stable.
- Ω_A′h²: `info.Omega_phi_h2`, `model_observables` → `Oh2_mediator`.
- `compare_limits_g` shades Ω_A′h² > 0.12 with `SHOW.mediator_relic` (off by default).

## EFT validity

| boundary | condition | where |
|---|---|---|
| mediatorTooLight | max(Gam_j/mj) ≥ 1 → error `axial_couplings:mediatorTooLight`; point → NaN | `axial_couplings.m`, `axial_gx_ceiling`, `axial_gq_ceiling` |
| perturbative unitarity | λ_a m_χ/m_A′ < √(π/2) (Kahlhoefer et al., JHEP 02 (2016) 016) | `model_observables` (MATLAB, Julia `_eft_flags`): `gx_ceiling = min(axial_gx_ceiling, √(π/2) m_A′/m_χ)` |
| small mixing | sin θ_j = m_jA′/m_A′ < `s_mix` (0.3) | `axial_mixing_floor.m` |

- mediatorTooLight: η′ binds first. At m_χ = 0.05, λ_a = 1e-3 the floor is m_A′ ≈ 0.01 GeV. Only mesons with j → χχ̄ open count, so it is `Inf` for m_χ > m_η′/2.
- Unitarity bound sets `unitarity_ok` and `valid` (GAMBIT treats violation as invalid).
- Quark-side bound a_q < √(π/2) m_A′/m_q is weaker than small mixing.
- Scaling: `gq_mix ∝ m_A′`, `gq_unit ∝ m_A′²`; cross near m_A′ ≈ 0.01 GeV.
- `compare_limits_g.m`: small-mixing boundary is an amber dashed line with upward shade.

## γγ ↔ χχ

- Production (`ggxxxsec.m`, `ggxxrate.m`): NWA notches at π⁰, η, η′; prefactor ∝ m_χ². No separate meson-decay term.
- Depletion (`xxggxsec.m`): γγ (3 poles) + A′A′ (t/u-channel, atan form, `real()`, gated m_χ ≥ m_φ) + meson pairs jj (atanh form, gated m_χ ≥ m_j).
- S_γγ: `xxggxsec.m` /16π, `ggxxxsec.m` /8π.
- Poles integrated by `thermal_pole_integral.m` (`cross_sections.md`).
