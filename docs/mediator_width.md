# Mediator width and numerical guards

## `phi_partial_widths(mphi, mx, gx, g, interaction)`

Single source for Γ(φ → …). Consumers read `.total` (`.gg` for inverse decay).

| interaction | channels |
|---|---|
| SS (`heavyScalar`, `lightScalar`) | γγ, χχ̄ (β³), meson pairs |
| `*PS` | γγ, χχ̄, meson pairs |
| `*SP` | γγ via π⁰/η mixing, χχ̄ (β³) |
| `scalarPP`, `*PP` | γγ via mixing, χχ̄ (β¹) |
| `axialVector` | see `axialvector.md` |

- Γ = S|M̄|² p*/(8π mφ²), p* = √(mφ²/4 − m²), S = ½ for γγ.
- axialVector Γ(A′ → χχ̄) uses 1/12π (polarisation average).
- Total width is the BW regulator in `xxggxsec.m`, `svavg.m`, `ggxxrate.m`, `ggxxxsec.m`. `freeze_coupled.m` and `neff_simple.m` use `G.pipi + G.KK + G.eta`.

## Meson pairs (SS, PS only) — `meson_pairs.m`

From L ⊃ −B0 φ tr(GΦ²), η = η₈. Γ(φ → AB) = S μ² √λ(mφ², mA², mB²) / (16π mφ³).

| pair | μ | S |
|---|---|---|
| π⁺π⁻ | B0(g_u + g_d) | 1 |
| π⁰π⁰ | B0(g_u + g_d) | 1/2 |
| K⁺K⁻ | B0(g_u + g_s) | 1 |
| K⁰K̄⁰ | B0(g_d + g_s) | 1 |
| ηη | B0(g_u + g_d + 4g_s)/3 | 1/2 |
| π⁰η | B0(g_u − g_d)/√3 | 1 |

Fields: `pipi` = π⁺π⁻ + π⁰π⁰, `KK` = K⁺K⁻ + K⁰K̄⁰, `eta` = ηη + π⁰η. `xxggxsec.m` χχ → AB uses the same Sμ².

## mφ = 0 guards

mφ → 0 is supported for scalar types. Guarded sites (keep in any new `besselk`, `1/mphi` or width consumer):

1. `phi_partial_widths.m`: meson-pair widths gated by threshold.
2. `neff_simple.m`, Julia `k1_over_k2`: K1/K2 → 0.
3. `freeze_coupled.m` `Y_eq_mass`: z²K2(z) → 2.
4. `phi_partial_widths.m` SS and SP χχ̄ widths gated on mphi > 2mx (0/0 otherwise hangs `ode15s`).

## Julia `_besselkx1` (`thermal_average.jl`)

- `_besselkx1(z) = z > 1e6 ? √(π/2z) : besselkx(1, z)`, used at all 7 scaled-kernel sites.
- Needed for mx ≲ 9.3e-5 GeV (Amos throws above ~1.07e9).
- Asymptote accurate to ~4e-7 at z = 1e6. Near-threshold nodes must not be zeroed.
- MATLAB `besselk` returns 0; not affected.
