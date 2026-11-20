module LoopAbundance

# --- foundation (Milestone 1) ---
include("thermo.jl")        # gstarofT
include("special.jl")       # besselk guards: k2_of, z2k2, k1_over_k2
include("constants.jl")     # physical constants (depends on gstarofT for s0)

# --- model interface + built-in models ---
include("models/interface.jl")
include("models/threebody.jl")   # SP/PP 3π / ππφ, axial 3π / π+π−A′ (off by default)
include("models/couplings.jl")   # SP/PP + axial coupling factories (M5)
include("models/builtins.jl")
include("models/widths.jl")
include("models/cross_sections.jl")

# --- thermal averaging (Milestone 2) ---
include("quadrature.jl")        # gl_nodes_weights, loggrid_edges
include("thermal_average.jl")   # thermal_pole_integral, svavg, svavg_grid, ggxxrate

# --- relic solvers (Milestone 3) ---
include("thermo_interp.jl")     # gstar_interp, _pchip
include("pigamma_phi.jl")       # pi gamma -> pi phi sigma + thermal rate (S⊗S/P⊗S)
include("primakoff_phi.jl")     # e gamma -> e phi, ee -> gamma phi (phi gamma gamma vertex)
include("axial_bath.jl")       # axialVector A'A' -> SM, A'pi -> pi pi (freeze_coupled)
include("relic.jl")             # freeze_loop, freeze_out, Oh2_FI, freeze_coupled

# --- ΔNeff + observables ---
include("neff_simple.jl")    # simple CMB ΔNeff (Boehm-Dolan-McCabe + φ from χ)
include("observables.jl")       # sigma_nucleon, model_observables

# LOOP_3PI / LOOP_PIPIPHI are read at load time, not precompile time.
__init__() = (threebody_refresh!(); nothing)

# thermodynamics
export gstarofT

# model interface
export CrossSectionModel
export sigma_ggxx, sigma_xxgg, sigma_xxall, poles, mediator_widths
export dof_dm, dof_mediator
export mass_dm, mass_mediator, coupling, quark_couplings

# built-in models
export HeavyScalar, LightScalar, HeavyScalarPS, LightScalarPS
export HeavyScalarSP, LightScalarSP, ScalarPP, HeavyScalarPP, LightScalarPP
export AxialVector, AxialMediatorTooLight
export model_from_interaction

# three-body final states (LOOP_3PI / LOOP_PIPIPHI)
export threebody_refresh!, threepi_enabled, pipiphi_enabled
export threepi_weight, pipiphi_weight, threebody_sigma, axial3pi_weight, axial3pi_sigma, pipiA_weight

# axialVector EFT-validity ceilings/floors
export axial_gx_ceiling, axial_eta_ceiling, axial_mixing_floor
export sppp_mixing_angles

# thermal averaging
export svavg, svavg_grid, svavg_grid_pair, svavg_phiphi, ggxxrate, thermal_pole_integral
export sigma_xxphiphi, sigma_xxdark, sv_darkpair_grid

# relic solvers
export freeze_loop, freeze_out, Oh2_FI, freeze_coupled
export gstar_interp

# ΔNeff + observables
export AA_sm_xsec, sv_AAsm_grid, apipi_phi_rate
export neff_simple, pigamma_phi_xsec, pigamma_phi_rate, primakoff_phi_xsec, primakoff_phi_rate
export sigma_nucleon, model_observables

end # module
