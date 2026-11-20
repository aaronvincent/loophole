# M3 validation: relic solvers (freeze_loop/out, Oh2_FI, freeze_coupled) vs
# MATLAB.  Regenerate the reference with test/matlab/dump_m3.m.
using LoopAbundance
using DelimitedFiles
using Printf

const REF = joinpath(@__DIR__, "ref")

raw = readdlm(joinpath(REF, "m3_ref.csv"), ','; skipstart=1)

# both-NaN agree; otherwise relative error (Oh2 spans many decades)
function relerr(a, b)
    (isnan(a) && isnan(b)) && return 0.0
    (isnan(a) || isnan(b)) && return Inf
    return abs(a - b) / max(abs(a), abs(b), 1e-300)
end

emax = Dict{String,Float64}()
worst = Dict{String,String}()
for r in 1:size(raw, 1)
    solver = String(raw[r, 1]); ic = String(raw[r, 2])
    mx, mphi, gx, T_RH, ref = raw[r, 3], raw[r, 4], raw[r, 5], raw[r, 6], raw[r, 7]
    m = model_from_interaction(ic, mx, mphi, gx, [1, 1, 1] / 3)
    jl = if solver == "freeze_loop";       freeze_loop(m; T_RH=T_RH)
    elseif solver == "freeze_out";          freeze_out(m)
    elseif solver == "Oh2_FI";              Oh2_FI(m, T_RH)
    elseif solver == "freeze_coupled";      freeze_coupled(m; T_RH=T_RH)
    else; error("unknown solver $solver"); end
    e = relerr(jl, ref)
    if e > get(emax, solver, 0.0)
        emax[solver] = e
        worst[solver] = @sprintf("%s mx=%.3g mphi=%.3g gx=%.3g  jl=%.6g ml=%.6g", ic, mx, mphi, gx, jl, ref)
    end
end

println("=== M3 relic solvers: max relerr per solver ===")
allmax = 0.0
for solver in ("freeze_loop", "freeze_out", "Oh2_FI", "freeze_coupled")
    e = get(emax, solver, 0.0)
    global allmax = max(allmax, e)
    @printf("  %-15s relerr = %.3e\n      worst: %s\n", solver, e, get(worst, solver, "—"))
end
@printf("\nOVERALL max relerr = %.3e  -> %s\n", allmax,
        allmax < 1e-2 ? "PASS (<1%)" : (allmax < 5e-2 ? "OK (<5%)" : "FAIL"))
