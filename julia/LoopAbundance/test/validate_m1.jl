# M1 validation: compare Julia cross sections + widths against the MATLAB
# reference dumps in test/ref/ (regenerate with test/matlab/dump_m1.m).
using LoopAbundance
using DelimitedFiles
using Printf

const REF = joinpath(@__DIR__, "ref")

model_for(tag) = startswith(tag, "hsps") ?
    HeavyScalarPS : HeavyScalar

# tag → (mx, mphi, gx, g)
function params(tag)
    pp = endswith(tag, "_pp")
    mx = pp ? 20.0 : 0.1
    return (mx, 10.0, 1.0, [1, 1, 1] / 3)
end

relerr(a, b) = (a == 0 && b == 0) ? 0.0 : abs(a - b) / max(abs(a), abs(b), 1e-300)

# --- widths ---
wraw = readdlm(joinpath(REF, "m1_widths.csv"), ','; skipstart=1)
println("=== widths ===")
wmax = 0.0
for r in 1:size(wraw, 1)
    tag = String(wraw[r, 1])
    mx, mphi, gx, g = params(tag)
    m = model_for(tag)(mx, mphi, gx, g)
    w = mediator_widths(m)
    ref = (gg=wraw[r,2], xx=wraw[r,3], pipi=wraw[r,4], KK=wraw[r,5], total=wraw[r,6])
    for f in (:gg, :xx, :pipi, :KK, :total)
        e = relerr(getfield(w, f), ref[f])
        global wmax = max(wmax, e)
    end
    @printf("  %-8s total: jl=%.6e ml=%.6e  relerr=%.2e\n", tag, w.total, ref.total, relerr(w.total, ref.total))
end
@printf("widths max relerr = %.3e\n\n", wmax)

# --- cross sections ---
raw = readdlm(joinpath(REF, "m1_ref.csv"), ','; skipstart=1)
println("=== cross sections ===")
emax = Dict("ggxx"=>0.0, "xxgg"=>0.0, "xxall"=>0.0)
worst = Dict("ggxx"=>("",0.0), "xxgg"=>("",0.0), "xxall"=>("",0.0))
for r in 1:size(raw, 1)
    tag = String(raw[r, 1]); s = raw[r, 2]
    mx, mphi, gx, g = params(tag)
    m = model_for(tag)(mx, mphi, gx, g)
    jl = (ggxx=sigma_ggxx(m, s), xxgg=sigma_xxgg(m, s), xxall=sigma_xxall(m, s))
    ref = (ggxx=raw[r,3], xxgg=raw[r,4], xxall=raw[r,5])
    for f in (:ggxx, :xxgg, :xxall)
        e = relerr(getfield(jl, f), getfield(ref, f))
        k = String(f)
        if e > emax[k]; emax[k] = e; worst[k] = (tag, s); end
    end
end
for k in ("ggxx", "xxgg", "xxall")
    @printf("  %-6s max relerr = %.3e  (worst: %s s=%.4g)\n", k, emax[k], worst[k][1], worst[k][2])
end

allmax = max(wmax, maximum(values(emax)))
@printf("\nOVERALL max relerr = %.3e  -> %s\n", allmax, allmax < 1e-6 ? "PASS (<1e-6)" : (allmax < 1e-2 ? "OK (<1%)" : "FAIL"))
