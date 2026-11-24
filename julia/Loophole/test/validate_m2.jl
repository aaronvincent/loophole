# M2 validation: thermal averages (svavg, svavg_grid, ggxxrate) vs MATLAB.
# Regenerate the reference with test/matlab/dump_m2.m.
using Loophole
using DelimitedFiles
using Printf

const REF = joinpath(@__DIR__, "ref")

# tag → (model-type, mx, mphi, gx, g)
function model_for(tag)
    isps = startswith(tag, "hsps")
    mx = tag == "hs_pp" ? 20.0 : (tag == "hs_nopole" ? 8.0 : 0.1)
    T = isps ? HeavyScalarPS : HeavyScalar
    return T(mx, 10.0, 1.0, [1, 1, 1] / 3)
end

# ⟨σv⟩ below ~1e-20 GeV⁻² is physically zero (Boltzmann-suppressed); MATLAB's
# adaptive integral floors such points to 0 at AbsTol 1e-20.  Treat both-tiny as
# agreeing so the x≫1 underflow tail doesn't masquerade as a relerr.
function relerr(a, b; atol=1e-20)
    abs(a - b) <= atol && return 0.0
    return abs(a - b) / max(abs(a), abs(b), 1e-300)
end

raw = readdlm(joinpath(REF, "m2_ref.csv"), ','; skipstart=1)
n = size(raw, 1)

# group rows by (tag) and cache grid evaluations per (tag, kind)
tags = unique(String.(raw[:, 1]))
# precompute the grid results per tag
Tall = sort(unique(Float64.(raw[:, 3])))
gridcache = Dict{Tuple{String,String},Dict{Float64,Float64}}()
for tag in tags
    m = model_for(tag)
    sg_all = svavg_grid(m, Tall)
    sg_gg  = svavg_grid(m, Tall; channel=:gg)
    gr     = ggxxrate(m, Tall)
    for (kind, vals) in (("svgrid_all", sg_all), ("svgrid_gg", sg_gg), ("ggrate", gr))
        d = Dict{Float64,Float64}()
        for (i, T) in enumerate(Tall); d[T] = vals[i]; end
        gridcache[(tag, kind)] = d
    end
end

emax = Dict(k => 0.0 for k in ("svgrid_all", "svgrid_gg", "sv_all", "sv_gg", "ggrate"))
worst = Dict(k => ("", 0.0) for k in keys(emax))
for r in 1:n
    tag = String(raw[r, 1]); kind = String(raw[r, 2]); T = raw[r, 3]; ref = raw[r, 4]
    m = model_for(tag)
    jl = if kind == "svgrid_all"; gridcache[(tag, "svgrid_all")][T]
    elseif kind == "svgrid_gg";   gridcache[(tag, "svgrid_gg")][T]
    elseif kind == "ggrate";      gridcache[(tag, "ggrate")][T]
    elseif kind == "sv_all";      svavg(m, T)
    elseif kind == "sv_gg";       svavg(m, T; channel=:gg)
    else; error("unknown kind $kind"); end
    e = relerr(jl, ref)
    if e > get(emax, kind, 0.0); emax[kind] = e; worst[kind] = (tag, T); end
end

println("=== M2 thermal averages: max relerr per kind ===")
allmax = 0.0
for kind in ("svgrid_all", "svgrid_gg", "sv_all", "sv_gg", "ggrate")
    e = get(emax, kind, 0.0)
    global allmax = max(allmax, e)
    @printf("  %-11s relerr = %.3e  (worst: %s T=%.3g)\n", kind, e, worst[kind][1], worst[kind][2])
end
@printf("\nOVERALL max relerr = %.3e  -> %s\n", allmax,
        allmax < 1e-4 ? "PASS (<1e-4)" : (allmax < 1e-2 ? "OK (<1%)" : "FAIL"))
