# σ_χN — `DMnucxsec.m`

## Interaction names

| prefix | form |
|---|---|
| `heavyScalar*` | heavy-mediator limit (1/mφ⁴ contact) |
| `lightScalar*` | light-mediator closed forms (SS: IR regulator q0 = mx·vrel; PS/SP/PP mφ-independent) |
| `scalar`, `scalarPS`, `scalarSP`, `scalarPP` | full angular-integrated form, any mφ |
| `axialVector` | full |

- `full_sigma_name.m` maps `heavyScalar*`/`lightScalar*` → `scalar*`. `model_observables.m`, `compare_limits*.m`, `dump_m5.m` use it: reported σ_χN is always the full form.
- `scalar`, `scalarPS`, `scalarSP` are σ_χN selectors only, not relic models. `scalarPP` is both.
- Julia `sigma_nucleon` is full-form only.

## Nucleon couplings (`nucleon_couplings.m`)

| quark vertex | interactions | coupling |
|---|---|---|
| scalar | SS, `*PS` | hp = mN Σ_i g_i f_{p,i}/m_i |
| pseudoscalar | `*SP`, `*PP` | h̃_P = mN Σ_i (Δ_i/m_i)(g_i − ⟨g/m⟩/⟨1/m⟩) |
| axial | `axialVector` | a_N = Σ_q a_q Δ_q |

- Δq^p = (0.847, −0.407, −0.035): FLAG 2024, arXiv:2411.04268, N_f = 2+1. Same set for pseudoscalar and axial.
- f_p = (0.02, 0.026, 0.118), arXiv:1202.1292.
- g ∝ [1 1 1] ⇒ h̃_P = 0 ⇒ SP/PP σ_χN = 0.

## Full form

σ = ∫_{−T}^{0} dσ/dt, T = 4μ²vrel², vrel = 220/3e5, propagator 1/(mφ² − t)².

Scalar family: `scalar_angint(x, al, be) = J0 + (al+be)J1 + al·be·J2`, x = T/mφ².

| | al | be | coupling |
|---|---|---|---|
| SS | 4mx²/mφ² | 4mN²/mφ² | hp |
| SP | 4mx²/mφ² | 0 | h̃_P |
| PS | 0 | 4mN²/mφ² | hp |
| PP | 0 | 0 | h̃_P |

Limits: SS → gx²hp²μ²/(π mφ⁴) heavy, ∝ 1/mφ² light. PP ∝ v⁴ heavy, mφ-independent light. PS/SP ∝ v² heavy, ∝ ln(1/mφ²) light.

axialVector: dσ/dΩ = λ_a²a_N² G(s,t) / (32π² mA′⁴ s (mA′² − t)²), G = G0 + G1 t + G2 t². Heavy → 3λ_a²a_N²μ²/(π mA′⁴); light → λ_a²a_N²μ²/(π mA′⁴).

`jmoments(x)` (Julia `_jmoments`): series below x = 0.05 (J0, J1 cancel to O(x³), O(x²)).

## Validation

- MATLAB vs Julia (`observables.jl`): ~4e-16.
- Closed form vs quadrature: ~3e-11 (`verify_dmnucxsec_fullrange.m`).
