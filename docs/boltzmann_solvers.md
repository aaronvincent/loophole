# Boltzmann solvers

## Which solver

| solver | evolves | φ treated as | use for |
|---|---|---|---|
| `freeze_loop` | Y_χ | virtual (s-channel propagator) | m_φ > 50·T_RH |
| `freeze_coupled` | Y_χ, Y_φ | bath species, own abundance | m_φ ≤ 50·T_RH |
| `Oh2_FI` | — (integral, Y² term dropped) | virtual | freeze-in branch only |
| `freeze_out` | Y_χ from equilibrium (arXiv:1204.3622) | virtual | freeze-out branch only |

- `freeze_loop` and `freeze_coupled` start from Y = 0 at T_RH and include production and depletion, so each covers freeze-in, freeze-out and the transition between them. Use one of these for Ωh².
- `Oh2_FI` and `freeze_out` are single-regime cross-checks.
- Threshold: at the default T_RH = 10 MeV, m_φ = 0.5 GeV. Above it φ is Boltzmann-suppressed throughout and the two solvers agree to ≲ 0.1%. Below it on-shell φ is produced in the bath and must be tracked.
- `freeze_coupled` costs ~5× more per call (0.015–0.03 s vs 0.003–0.006 s).
- Automatic dispatch at m_φ = 50·T_RH: `model_observables` (MATLAB; Julia `observables.jl`), `compare_limits_g.m`. MATLAB `freeze_loop.m` and `freeze_coupled.m` warn when called on the other side of the threshold.
- GAMBIT scans use Julia `model_observables`, so dispatch is automatic.
- In `parfor`, broadcast `solver_name` as a string, not a function handle.

## Common

- λ = ½ √(π/45) M_Pl m_χ ⟨σv⟩. Y counts χ + χ̄ (g_χ = 4); ½ is the Dirac factor. Same in every solver and `relic.jl`.
- Y_eq = (45/4π⁴)(g_χ/g_s) x² K2(x); asymptotic expansion for x > 30.
- Integrated in τ = log x with `ode15s`, explicit Jacobian, RelTol 1e-6, AbsTol 1e-25.
- ⟨σv⟩ tables interpolated in log.
- IC: Y = 0 at x_init = mx/T_RH (T_RH = 0.01 GeV default). First step Δτ = 10⁻³/max(1, x_init) is analytic: Y = Y* tanh(√(a c_f) Δτ), Y* = √(a/c_f), a = c_p Y_eq².
- `production = 'sm'` (default) or `'gg'` (photon-only). Set via function argument, `batch_params.production`, `RUNSET.production`; Julia `production=:gg`.

## `freeze_loop.m`

dY/dx = (1/x²) · correction · (g_s/√g_ρ) · (λ_prod Y_eq² − λ_full Y²)

- λ_prod: `svavg(..., 'sm')` = γγ + meson pairs (ππ, KK; jj for SP/PP/axialVector).
- λ_full: `svavg(...)` = 'sm' + φφ, A′A′, φ′φ′.
- Meson channels are 10⁴–10⁶× γγ; they dominate production for mχ ≳ 50–70 MeV (T_RH = 10 MeV).

## `Oh2_FI`

- Y = ∫ 2R/(sH) (1 + d ln g_s / 3 d ln T) d ln T, R = n_eq²⟨σv⟩/4 (`ggxxrate` + meson term).
- s₀ = 2891.2 cm⁻³, ρ_c/h² = 1.054e-5 GeV cm⁻³.
- FI/loop = 0.99–1.01 at weak coupling.

## Relic roots (`compare_limits_g.m`)

Ωh²(g_q) ∝ g_q² rising, then ∝ 1/g_q². Log grid in g_q, `fzero` at each sign change. Smallest root = FIMP, largest = WIMP, NaN if Ωh² never reaches 0.12.

## `freeze_coupled.m`

ρ = Y_χ/Y_χ,eq:

```
dY_χ/dτ = cp_s Y_χeq² − (cf_s + c_d) Y_χ² + 2 c_xx (Y_φ − Y_φeq ρ²) + c_rev Y_φ²
dY_φ/dτ = −c_SM (Y_φ − Y_φeq) − c_xx (Y_φ − Y_φeq ρ²) − c_mes Y_φ + c_d Y_χ² − c_rev Y_φ²
```

| term | definition |
|---|---|
| c_xx | Γ(φ → χχ̄)⟨K1/K2⟩/H |
| c_SM | [(Γ_γγ + Γ_mesons if 'sm')⟨K1/K2⟩ + Γ_sc]/H |
| Γ_sc | Primakoff (`primakoff_phi_rate`, not axialVector) + π±γ → π±φ (`pigamma_phi_rate`, S⊗S/P⊗S, mφ < 2m_π, 'sm') + A′π ↔ ππ (`apipi_phi_rate`, axialVector, 'sm'); 120 log points, 5 keV – T_RH |
| c_mes | φ → mesons, no inverse; 'gg' only |
| cp_s, cf_s | SM production / depletion minus the on-shell φ part, 2c_xx(Y_φeq/Y_χeq²)·BR(φ → SM); clamped ≥ 0 |
| c_d, c_rev | χχ̄ ↔ φφ (A′A′, φ′φ′) from `sv_darkpair_grid`; c_d = c_rev (Y_φeq/Y_χeq)² |
| c_AA | axialVector: A′A′ ↔ γγ, π⁺π⁻, K⁺K⁻, K⁰K̄⁰; term −c_AA(Y_φ² − Y_φeq²), c_AA = 2λ⟨σv⟩ (`AA_sm_xsec.m`, `sv_AAsm_grid.m`) |

- g_φ = 3 for axialVector.
- Reverse ⟨σv⟩ normalised by φ's K2; only the reverse table is interpolated (exact detailed balance).
- A′π ↔ ππ: axial 3π contact term crossed to 2 → 2, normalised by Δa_eff = m_{π⁰A′}/F. Γ/H ≈ 4×10¹¹ g_q² at mA′ = 0.1, T = 10 MeV. Example (mA′ = 0.1, λ_χ = 0.01): FIMP roots for mχ = 0.07–0.11 move from g_q ~ 10⁻⁶–10⁻⁴ to ~10⁻⁷; mχ = 0.11 WIMP root 0.066 → 0.018.

Env switches: `LOOP_PHI_SC=0` (no Γ_sc), `LOOP_APIPI=0` (no A′π ↔ ππ), `LOOP_AA_SM=0` (no c_AA), `LOOP_FC_FULL2D=1` (no fast-φ elimination), `LOOP_FC_NOEXCH=1` (no dark-exchange shortcut).

### Mediator relic

`[Oh2, info] = freeze_coupled(...)`:
- `info.Omega_phi_h2` = m_φ s₀ Y_φ(x_final)/ρ_c (as if stable)
- `info.phi_stable` (Γ_φ t₀ < 1), `info.overclosed` (stable and > 0.12)
- `model_observables`: `obs.Oh2_mediator`, `obs.mediator_overclosed`
- `compare_limits_g`: `SHOW.mediator_relic`; `.dat` columns `gq_med_lo/hi`, `sg_med_lo/hi`
- axialVector, m_A′ = 0.1, gx = 0.01, on the χ relic locus: Ω_A′h² ≲ 0.01 for mχ < m_A′, ≈ 1 at mχ = 0.11

### Fast φ

Condition: (Γ_φ,tot⟨K1/K2⟩ + Γ_sc)/H > 10³ over the range.
- Y_φ set to its quasi-steady value (closed-form quadratic in Y_φ − Y_φeq); only χ integrated.
- Elastic loops χχ̄ → φ → χχ̄ and χχ̄ → φφ → χχ̄ cancel exactly.
- Coefficients may be negative (χχ̄ → φφ → 4χ); `ann_logY` handles this.

### Dark-exchange shortcut

Condition: no fast-φ path and χχ̄ ↔ φφ > 10⁹ per e-fold (e.g. axialVector, mχ > m_A′/2).
- Y_φ = r Y_χ, r = Y_φeq/Y_χeq; evolve N = Y_χ + Y_φ.
- Switches to the 2×2 system when exchange < 10³ per e-fold.
- Agrees with 2×2 to ≲ 0.1% where both run; 10–500× faster.

### Initial conditions

- Y_φ(T_RH) = Y_φeq (1 − e^{−Γ_SM dτ₀/H}).
- Y_χ: analytic seed (above). Fast-φ path evaluates the saturated seed at ρ = 1. Integrated as W = log Y_χ when seed > Y*/2, linear otherwise.
- ICs are built at x0 = exp(τ0).
- `snap_unit_ratio` (tol 1e-12) sets ⟨σv⟩_gg/⟨σv⟩_full = 1 exactly when equal to round-off (S⊗P, P⊗P).
- Annihilation term: `cf·e^W·expm1(log(cp/cf) − 2(W − W_eq))` (`ann_logY`).

### φφ ↔ χχ̄ effect size

- Light φ in bath (P⊗S, mφ = 1e-9, mχ = 30–55 MeV, g_q ~ 10⁻⁵–10⁻⁴): Ωh² ×4–13.
- φ out of bath, freeze-out branch, mφ = 0.1 freeze-in: ≲ 1%.

### Julia

`relic.jl`, `src/axial_bath.jl`, `freeze_coupled(m; return_info=true)`. Agreement at mA′ = 0.1: σ, ⟨σv⟩ 10 digits; rate 3×10⁻⁹; Ω_χ ≤ 0.25%; Ω_A′ ≤ 1.5%.

Tests: `test_freeze_coupled_phiphi.m`, Julia `validate_fc_sc`.
