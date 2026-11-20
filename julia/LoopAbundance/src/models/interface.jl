# The generic cross-section-model interface.
#
# Everything downstream (thermal averaging, all relic solvers, ΔNeff) talks to
# a model ONLY through the methods below — it never branches on an interaction
# string.  A user supplies their own `struct MyModel <: CrossSectionModel` plus
# these methods and the solvers work unchanged.  The built-in models
# (HeavyScalar, HeavyScalarPS, ...) are concrete instances that reproduce the
# original MATLAB `interaction`-string physics and validate the interface.

"""
    CrossSectionModel

Abstract supertype for a DM model.  A concrete model carries its parameters
(mx, mφ, gx, quark couplings g, plus any precomputed vertex couplings) and
implements the cross sections, resonance poles, widths, and dof below.

Required methods:
  sigma_ggxx(m, s)    γγ→χχ cross section [GeV⁻²]   (freeze-in production)
  sigma_xxgg(m, s)    χχ→γγ, s-channel only         (production partner / 'gg')
  sigma_xxall(m, s)   χχ→all open channels          (depletion)
  sigma_xxhad(m, s)   χχ→meson pairs only           (optional; default 0)
  poles(m)            Vector of (s_p = m_p², w_p = m_p·Γ_p) for the integrator
  mediator_widths(m)  NamedTuple (gg, xx, pipi, KK, eta, total) [GeV]
  dof_dm(m)           DM internal dof g_χ
  dof_mediator(m)     mediator internal dof g_φ

`sigma_*` take a scalar Mandelstam `s` [GeV²] and return σ [GeV⁻²]; broadcast
over an s-grid with `sigma_ggxx.(m, svec)`.
"""
abstract type CrossSectionModel end

# Parameter accessors — concrete models store these fields by convention.
mass_dm(m::CrossSectionModel)        = m.mx
mass_mediator(m::CrossSectionModel)  = m.mphi
coupling(m::CrossSectionModel)       = m.gx
quark_couplings(m::CrossSectionModel) = m.g

# Defaults shared by most models (override per type as needed).
dof_dm(::CrossSectionModel)       = 4.0   # Dirac fermion: 2 spin × (χ + χ̄)
dof_mediator(::CrossSectionModel) = 1.0   # real scalar

# Interface methods — concrete models must implement these.
function sigma_ggxx end
function sigma_xxgg end
function sigma_xxall end
function poles end
function mediator_widths end

# χχ̄→φφ (mediator-pair) cross section [GeV⁻²], the t/u-channel piece isolated
# from sigma_xxall.  Used ONLY by svavg_phiphi, for the dark-sector rate
# (χχ↔φφ) — it is not a detailed-balance partner of
# any production channel.  Default 0 (model has no such channel); concrete
# models with a mediator-pair final state override.  Port of
# make_xxsigma(...,'phiphi').
sigma_xxphiphi(::CrossSectionModel, ::Real) = 0.0

# χχ̄→meson pairs (ππ, KK; SP/PP/axialVector jj) [GeV⁻²], the SM-hadron piece of
# sigma_xxall.  With sigma_xxgg it makes the 'sm' channel: the detailed-balance
# partner of χ production from the equilibrium SM bath.  Default 0.  Port of
# xxggxsec(...,'had').
sigma_xxhad(::CrossSectionModel, ::Real) = 0.0
