# K → πφ, CMB, INTEGRAL, meson → invisible

| key | quantity | bound | applies to |
|---|---|---|---|
| `kphi` | Γ(K⁺ → π⁺φ)·P_inv × NA62 acceptance A(mφ²) | BR(K⁺ → π⁺ inv) < 1.14e-10 | mφ < m_K± − m_π± ≈ 0.354 GeV |
| `cmb` | ⟨σv⟩ at s = 4mχ² | Slatyer energy injection | P⊗S, P⊗P, A⊗A |
| `integral` | ⟨σv⟩ at s = 4mχ² | INTEGRAL γ-line | all |
| `mesoninv` | Γ(j → χχ̄), j = π⁰, η, η′ | 3.43e-17, 3.14e-11, 3.95e-8 GeV | S⊗P, P⊗P, A⊗A |

Excluded where ratio > 1. Interaction → model map as in `kaon_constraint.md`.

## `kphi`

- P_inv = (Γ_χχ̄ + Γ_vis e^{−L/ℓ})/Γ_tot, ℓ = (p_φ/mφ) c/Γ_tot, L = 100 m, E_φ = 50 GeV.
- Visible = γγ + meson pairs (ππ above 0.270 GeV, S⊗S/P⊗S).
- Widths from `phi_partial_widths.m` / `mediator_widths`. A⊗A: P_inv = 1.
- S⊗P/P⊗P: θ_j from `sppp_mixing` (π⁰, η, η′).
- A⊗A: Γ = G_F² F⁴ V_us² V_ud² (λ_A′K + λ_A′π)² λ^{3/2} / (32π m_A′² m_K³).
- mφ > 2mχ: P_inv = 1, bound independent of gx (horizontal band).
- mφ < 2mχ: excluded g_q is a band. Floor `ratio_crit_coupling` (seed g_q = 1e-12, quarter-decade steps both ways); ceiling `ratio_ceiling_coupling` (to g_q = 1e8; NaN = excluded to top).
- On shell, kaon (K → πχχ̄) = kphi × BR(φ → χχ̄) to 1e-5.

## `cmb`, `integral`

- ⟨σv⟩ from analytic low-velocity forms: `sigmavheavyscalar`, `sigmavPcrossP`, `sigmavScrossP`, `sigmavPcrossS`, `sigmavheavyAxialAxial`. Not from `svavg`.
- ⟨σv⟩ ∝ gx²: evaluate at gx = 1, gx_lim = 1/√ratio.
- P⊗S pole width includes meson pairs (`meson_pairs.m`). S⊗S ⟨σv⟩ has no width term.
- Limit curves inlined (~20 points each).
- `LOOP_3PI` / `LOOP_PIPIPHI` (CMB only): s-wave P⊗P 3π (`threepi_weight`), ππφ (`pipiphi_weight`), axialVector J = 0 3π (`axial3pi_weight` W_L). Energy weights κ_3π = 0.52, κ_ππφ = 0.18. With `LOOP_3PI=1` the S⊗P/P⊗P width includes Γ(φ → 3π) (CMB and INTEGRAL).

## `mesoninv`

- S⊗P, P⊗P: g_jχ = gx θ_j (`sppp_mixing.m`). Γ = g² m_j β/8π (P⊗P), g² m_j β³/8π (S⊗P).
- A⊗A: Γ = λ_jχ² mχ² √(m_j²/4 − mχ²)/π (`axial_couplings.m`). Below the meson-width floor λ_jχ is rescaled from a probe (finite ratio).
- S⊗S, P⊗S: NaN.
- Closed for mχ ≥ m_j/2. `which = 'pi0inv' | 'etainv' | 'etapinv' | 'mesoninv'` (max).
- π⁰ is nonzero for g ∝ [1 1 1] (L7 term ∝ m_u − m_d).
- Γ and freeze-in both ∝ (gx g_q)²: gap to FI locus is gx-independent.

## Files

| layer | file | entry points |
|---|---|---|
| MATLAB (reference) | `matlab/mesonandCMB_ratio.m` | `mesonandCMB_ratio(mx, mphi, gx, g, interaction, which)` |
| MATLAB loglike | `matlab/external_loglikes.m` | `kphi_loglike`, `cmb_loglike`, `integral_loglike`, `mesoninv` |
| Python | `analysis/mesonandCMB.py` | `kphi_ratio`, `cmb_ratio`, `integral_ratio`, `meson_invisible_ratio`, `kind_of` |

Python vs MATLAB: ~4e-13.

## Postprocess

- `postprocess.py`: `kphi()`, `cmb()`, `integral()` (via `_mcmb`), `mesoninv()`.
- Datasets `LogLike_<key>`, `excluded_<key>` (in `excluded_any`), provenance `input::ext_sigma_<key>`.
- Width: `EXT_SIGMA_KPHI`, `EXT_SIGMA_CMB`, `EXT_SIGMA_INTEGRAL`, `EXT_SIGMA_MESONINV` (−1 = hard cut; > 0 = soft width).
- Registered in `constraint_io.py` (`CONSTRAINT_COLS`), `plot_all_constraints.py` (`kphi` #a6761d, `cmb` #17becf, `integral` #7570b3).
- `compare_limits_g.m`: `mesoninv` is the 5th ratio row (`SHOW.mesoninv`, navy).
- Existing sidecars: rerun with `POSTPROC_FORCE=1`.
