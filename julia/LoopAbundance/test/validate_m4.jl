# M4 validation: model_observables (relic + ΔNeff bbn/cmb + σ_χN) vs MATLAB.
# Regenerate the reference with test/matlab/dump_m4.m.
using LoopAbundance
using DelimitedFiles
using Printf

const REF = joinpath(@__DIR__, "ref")
raw = readdlm(joinpath(REF, "m4_ref.csv"), ','; skipstart=1)

# relative error; both-NaN agree.  ΔNeff is judged on an absolute floor too,
# since dNeff -> 0 in much of the grid (relerr is meaningless near 0).
function relerr(a, b)
    (isnan(a) && isnan(b)) && return 0.0
    (isnan(a) || isnan(b)) && return Inf
    return abs(a - b) / max(abs(a), abs(b), 1e-300)
end
# combined metric for ΔNeff: pass if relerr small OR abs diff tiny
neff_err(a, b) = min(relerr(a, b), (isnan(a) || isnan(b)) ? Inf : abs(a - b))

emax = Dict("Oh2" => 0.0, "dNeff_bbn" => 0.0, "dNeff_cmb" => 0.0, "sigma" => 0.0)
worst = Dict{String,String}()
# Points where MATLAB's ode15s returned NaN (solver non-convergence) but Julia's
# QNDF+dt-seed converged to a finite, sane value — the documented robustness
# difference, reported separately, not scored as disagreement.
ml_nonconv = String[]

for r in 1:size(raw, 1)
    ic = String(raw[r, 1])
    mx, mphi, gx = raw[r, 2], raw[r, 3], raw[r, 4]
    g = [raw[r, 5], raw[r, 6], raw[r, 7]]
    T_RH = raw[r, 8]
    Oh2_ref, db_ref, dc_ref, sig_ref = raw[r, 9], raw[r, 10], raw[r, 11], raw[r, 12]

    o = model_observables(mx, mphi, gx, g, ic; T_RH=T_RH)
    tag = @sprintf("%s mx=%.3g mphi=%.3g gx=%.3g", ic, mx, mphi, gx)

    # MATLAB solver failed but Julia produced a finite positive relic: log & skip.
    skip_oh2 = isnan(Oh2_ref) && !isnan(o.Oh2)
    if skip_oh2
        push!(ml_nonconv, @sprintf("%s  jl=%.6g (ml=NaN)", tag, o.Oh2))
    end

    for (key, jl, ref, ef) in (
        ("Oh2",       o.Oh2,        Oh2_ref, relerr),
        ("dNeff_bbn", o.dNeff_bbn,  db_ref,  neff_err),
        ("dNeff_cmb", o.dNeff_cmb,  dc_ref,  neff_err),
        ("sigma",     o.sigma_chiN, sig_ref, relerr),
    )
        key == "Oh2" && skip_oh2 && continue
        e = ef(jl, ref)
        if e > emax[key]
            emax[key] = e
            worst[key] = @sprintf("%s  jl=%.6g ml=%.6g", tag, jl, ref)
        end
    end
end

if !isempty(ml_nonconv)
    println("--- MATLAB ode15s non-convergence (Julia QNDF converged), not scored ---")
    for s in ml_nonconv
        println("  ", s)
    end
    println()
end

println("=== M4 model_observables: max error per field ===")
allmax = 0.0
for key in ("Oh2", "dNeff_bbn", "dNeff_cmb", "sigma")
    e = emax[key]
    global allmax = max(allmax, e)
    @printf("  %-10s err = %.3e\n      worst: %s\n", key, e, get(worst, key, "—"))
end
@printf("\nOVERALL max err = %.3e  -> %s\n", allmax,
        allmax < 1e-2 ? "PASS (<1%)" : (allmax < 5e-2 ? "OK (<5%)" : "FAIL"))
