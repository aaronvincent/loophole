# Precomputed g_s, g_rho, d(ln g_s)/d(ln T) interpolants for fast lookup inside
# an ODE RHS (replaces 4 gstarofT calls per evaluation).  Port of
# gstar_interp.m.  Interpolants are keyed on log(T) (callers pass log(T)).

using DataInterpolations: CubicSpline, PCHIPInterpolation, ExtrapolationType

struct GstarInterp{F}
    gs::F
    grho::F
    dlngs::F
end

# evaluate at log(T)
gs_at(S::GstarInterp, lT)    = S.gs(lT)
grho_at(S::GstarInterp, lT)  = S.grho(lT)
dlngs_at(S::GstarInterp, lT) = S.dlngs(lT)

"""
    gstar_interp(Tlo, Thi; n=400) -> GstarInterp

Cubic-spline interpolants of g_s, g_ρ, and d ln g_s/d ln T over a log-T grid
padded ±0.05 dex beyond [Tlo, Thi] (so the ODE may probe the endpoints).
"""
function gstar_interp(Tlo::Real, Thi::Real; n::Integer=400)
    pad = 0.05
    lT = range(log10(Tlo) - pad, log10(Thi) + pad; length=n)
    Tg = 10.0 .^ lT
    gs = gstarofT.(Tg, 1)
    grho = gstarofT.(Tg, 0)
    dTd = Tg .* 1e-5
    dlngs = (gstarofT.(Tg .+ dTd, 1) .- gstarofT.(Tg .- dTd, 1)) ./ (2 .* dTd) .* (Tg ./ gs)
    lTln = log.(Tg)                       # natural-log abscissa (ascending)
    ext = ExtrapolationType.Constant
    return GstarInterp(
        CubicSpline(gs, lTln; extrapolation=ext),
        CubicSpline(grho, lTln; extrapolation=ext),
        CubicSpline(dlngs, lTln; extrapolation=ext),
    )
end

# pchip on a log-T grid, with constant extrapolation just past the endpoints
# (the relic ODEs map their τ-endpoints to the grid boundaries up to roundoff).
_pchip(t, u) = PCHIPInterpolation(u, t; extrapolation=ExtrapolationType.Constant)
