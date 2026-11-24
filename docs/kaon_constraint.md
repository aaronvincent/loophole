# K⁺ → π⁺χχ̄

r = Γ(K⁺ → π⁺χχ̄) / (1.14e-10 · Γ_K). Excluded where r > 1. Open for mχ < (m_K± − m_π±)/2 ≈ 0.177 GeV, else r = NaN.

## Models

| branch | interactions | φ–χ | φ–quark |
|---|---|---|---|
| S⊗S | `heavyScalar`, `lightScalar`, `scalar` | scalar | scalar (φ BW) |
| P⊗S | `*PS` | pseudoscalar | scalar (φ BW) |
| S⊗P | `*SP` | scalar | pseudoscalar (π⁰/η poles) |
| P⊗P | `scalarPP`, `*PP` | pseudoscalar | pseudoscalar (π⁰/η poles) |
| A⊗A | `axialVector` | axial | axial (π⁰/η/η′ mixing) |

- S⊗S/P⊗S width: Γ_γγ + Γ_χχ̄ + meson pairs (`meson_pairs.m`; Python `mesonandCMB.meson_pair_width`).
- S⊗P/P⊗P: full width (`phi_partial_widths`), `sppp_mixing` NLO mixing with η′ pole, K → πη′ vertex (√6/3)(2 + 3g₈ˢ/g₈)(m_K² − m_π²), g_jχ = +g_χ θ_j. Vanishes for g ∝ [1 1 1].
- A⊗A: signed λ_jπK; A′ propagator in the interference term; Γ_A′ = relic total width. Inner Eχ integral analytic.
- r is not ∝ gx² on shell; root-find with `ratio_crit_coupling`.

## NA62 acceptance

A(s) = Σ_regions [Φ((s_hi − s)/σ) − Φ((s_lo − s)/σ)], s = m²_miss.
- R1 = [0, 0.010], R2 = [0.026, 0.068] GeV² (JHEP 06 (2021) 093); σ = 0.001 GeV².
- π⁰ peak (s = 0.0182) gets weight ~1e-16; π⁰ → invisible belongs to `mesoninv`.
- `na62_acceptance` in `kaon_invisible_ratio.m`, `mesonandCMB_ratio.m`, `kaon_invisible.py`. Also applied to `kphi`.

## Invisible-BR ceiling

1.14e-10. Set by `_KINV` in `kaon_invisible.py` and the matching constant in MATLAB.

## Python port — `analysis/kaon_invisible.py`

```python
kaon_ratio(mx, mphi, gx, gu, gd, gs, interaction) -> r   # NaN where closed
kind_of(interaction) -> "SxS" | "PxS" | "SxP" | "PxP" | "AxA" | None
```

- Fixed Gauss–Legendre in E_π: `KNODES = 64` (`KAON_KNODES`), per NA62 panel.
- Panels split at acceptance edges ± 4σ and s = 0.017, 0.019, 0.075. (0.017, 0.019) and s > 0.075 skipped (`na62_dead` in MATLAB).
- Narrow poles (φ BW in S⊗S/P⊗S, π⁰ in S⊗P/P⊗P, on-shell A′) excised and added analytically: ∫ dE_π/((s − m²)² + w²) = atan(2 m_K δ/w)/(m_K w).
- Parallel: `KAON_CHUNK_PAR` (10 000) per process, `KAON_NPROC` processes (default `SLURM_CPUS_PER_TASK`).
- Memory: `KAON_CHUNK` (50 000 rows per chunk).

Validated against MATLAB on 500 points.

MATLAB S⊗S/P⊗S on-shell pole: NWA notch (`bw_notch_integral`).

## Postprocess

- `kaon()` → `LogLike_kaon` (`_softregion(log10 r, sigma)`), `excluded_kaon` (uint8, in `excluded_any`).
- `EXT_SIGMA_KAON`: −1 (default) = hard cut, −10⁶ per excluded point; > 0 = soft width in log10 r.
- Registered in `constraint_io.py` (`CONSTRAINT_COLS`), `plot_all_constraints.py` (#66a61e, `K⁺→π⁺χχ`).
- Existing sidecars: rerun with `POSTPROC_FORCE=1`.
