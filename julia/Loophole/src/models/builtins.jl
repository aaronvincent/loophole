# Built-in cross-section models — concrete instances of the CrossSectionModel
# interface that reproduce the original MATLAB `interaction`-string physics.
#
# Milestone 1 ships the single-pole scalar models (HeavyScalar / LightScalar /
# HeavyScalarPS / LightScalarPS).  SP/PP and AxialVector follow in M5.

# All four single-pole scalar models share the same fields and a constructor
# that accepts any Real params and any 3-element g (Vector or Tuple).  An inner
# constructor (not an outer one) avoids colliding with Julia's auto-generated
# convert constructor.  heavy / light are the SAME formula (the Breit-Wigner
# propagator covers both the s ≪ mφ² and s ~ mφ² regimes); they are distinct
# types only so the model round-trips to/from the legacy interaction string.
macro defmodel(name)
    quote
        struct $(esc(name)) <: CrossSectionModel
            mx::Float64
            mphi::Float64
            gx::Float64
            g::NTuple{3,Float64}
            $(esc(name))(mx, mphi, gx, g) =
                new(mx, mphi, gx, (Float64(g[1]), Float64(g[2]), Float64(g[3])))
        end
    end
end

@defmodel HeavyScalar
@defmodel LightScalar
@defmodel HeavyScalarPS
@defmodel LightScalarPS

# Unions for shared implementations (scalar-scalar vs scalar-pseudoscalar).
const ScalarSS = Union{HeavyScalar, LightScalar}
const ScalarPS = Union{HeavyScalarPS, LightScalarPS}

# Single Breit-Wigner φ pole for all four; width = φ total width.
poles(m::Union{ScalarSS,ScalarPS}) = [(m.mphi^2, m.mphi * mediator_widths(m).total)]

# =====================================================================
# Multi-pole models (M5): SP/PP family + AxialVector.
# =====================================================================
# SP/PP models store the precomputed meson-pole couplings (mj, gcx, ggx, Gam)
# from _sppp_couplings so the cross sections / widths / poles never recompute.
# heavy/light/scalar variants are distinct types only for round-tripping to the
# legacy interaction string; the physics is identical within SP vs PP.
for (name, isPP) in ((:HeavyScalarSP, false), (:LightScalarSP, false),
                     (:ScalarPP, true), (:HeavyScalarPP, true), (:LightScalarPP, true))
    @eval begin
        struct $name <: CrossSectionModel
            mx::Float64
            mphi::Float64
            gx::Float64
            g::NTuple{3,Float64}
            mj::Vector{Float64}
            gcx::Vector{Float64}
            ggx::Vector{Float64}
            Gam::Vector{Float64}
        end
        function $name(mx, mphi, gx, g)
            gt = (Float64(g[1]), Float64(g[2]), Float64(g[3]))
            mj, gcx, ggx, Gam = _sppp_couplings(Float64(mphi), Float64(mx), Float64(gx), gt, $isPP)
            $name(Float64(mx), Float64(mphi), Float64(gx), gt, mj, gcx, ggx, Gam)
        end
    end
end

const PoleSP   = Union{HeavyScalarSP, LightScalarSP}
const PolePP   = Union{ScalarPP, HeavyScalarPP, LightScalarPP}
const PoleSPPP = Union{PoleSP, PolePP}

_is_pp(::PoleSP) = false
_is_pp(::PolePP) = true

# AxialVector: spin-1 mediator A' (mφ=m_A', gx=λ_a, g=[a_u a_d a_s]).  Stores the
# A'-meson mixing couplings; construction throws AxialMediatorTooLight on the
# ultralight-A' validity floor.  g_φ = 3 (massive vector).
struct AxialVector <: CrossSectionModel
    mx::Float64
    mphi::Float64
    gx::Float64
    g::NTuple{3,Float64}
    mj::Vector{Float64}
    lam_chi::Vector{Float64}
    lam_gg::Vector{Float64}
    Gam_j::Vector{Float64}
end
function AxialVector(mx, mphi, gx, g)
    gt = (Float64(g[1]), Float64(g[2]), Float64(g[3]))
    mj, lam_chi, lam_gg, Gam_j = _axial_couplings(Float64(mphi), Float64(mx), Float64(gx), gt)
    AxialVector(Float64(mx), Float64(mphi), Float64(gx), gt, mj, lam_chi, lam_gg, Gam_j)
end

dof_mediator(::AxialVector) = 3.0   # massive vector A' (3 polarizations)

# Meson s-channel poles (s_p = m_j², w_p = m_j·Γ_j) — matches sv_poles/gg_poles.
poles(m::PoleSPPP)    = [(m.mj[j]^2, m.mj[j] * m.Gam[j]) for j in eachindex(m.mj)]
poles(m::AxialVector) = [(m.mj[j]^2, m.mj[j] * m.Gam_j[j]) for j in eachindex(m.mj)]

"""
    model_from_interaction(interaction, mx, mphi, gx, g) -> CrossSectionModel

Bridge from the legacy MATLAB `interaction` string to a model object, so
`model_observables` and the scan drivers keep the original calling convention.
"""
function model_from_interaction(interaction::AbstractString, mx, mphi, gx, g)
    if     interaction == "heavyScalar";    return HeavyScalar(mx, mphi, gx, g)
    elseif interaction == "lightScalar";    return LightScalar(mx, mphi, gx, g)
    elseif interaction == "heavyScalarPS";  return HeavyScalarPS(mx, mphi, gx, g)
    elseif interaction == "lightScalarPS";  return LightScalarPS(mx, mphi, gx, g)
    elseif interaction == "heavyScalarSP";  return HeavyScalarSP(mx, mphi, gx, g)
    elseif interaction == "lightScalarSP";  return LightScalarSP(mx, mphi, gx, g)
    elseif interaction == "scalarPP";       return ScalarPP(mx, mphi, gx, g)
    elseif interaction == "heavyScalarPP";  return HeavyScalarPP(mx, mphi, gx, g)
    elseif interaction == "lightScalarPP";  return LightScalarPP(mx, mphi, gx, g)
    elseif interaction == "axialVector";    return AxialVector(mx, mphi, gx, g)
    else
        error("model_from_interaction: interaction '$interaction' not recognized")
    end
end
