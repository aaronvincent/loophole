# M5 validation: multi-pole models (SP/PP family + axialVector).  Cross sections
# (sigma_ggxx / sigma_xxgg / sigma_xxall), mediator widths, meson poles, and
# sigma_nucleon vs MATLAB.  Regenerate the reference with test/matlab/dump_m5.m.
using Loophole
using DelimitedFiles
using Printf

const REF = joinpath(@__DIR__, "ref")
raw = readdlm(joinpath(REF, "m5_ref.csv"), ','; skipstart=1)

function relerr(a, b)
    (isnan(a) && isnan(b)) && return 0.0
    (isnan(a) || isnan(b)) && return Inf
    # floor tiny magnitudes (sub-threshold / near-zero cross sections) — MATLAB
    # and Julia both round these to ~0; relerr there is meaningless.
    (abs(a) < 1e-30 && abs(b) < 1e-30) && return 0.0
    return abs(a - b) / max(abs(a), abs(b), 1e-300)
end

# evaluate the Julia quantity matching a reference row
function jl_value(m, q, arg)
    if     q == "ggxx";   return sigma_ggxx(m, arg)
    elseif q == "xxgg";   return sigma_xxgg(m, arg)
    elseif q == "xxall";  return sigma_xxall(m, arg)
    elseif q == "wgg";    return mediator_widths(m).gg
    elseif q == "wxx";    return mediator_widths(m).xx
    elseif q == "wtot";   return mediator_widths(m).total
    elseif q == "pole_s"; return poles(m)[Int(arg)][1]
    elseif q == "pole_w"; return poles(m)[Int(arg)][2]
    elseif q == "signuc"; return sigma_nucleon(m)
    elseif q == "Oh2";    return model_observables(m; T_RH=0.01).Oh2
    elseif q == "dNbbn";  return model_observables(m; T_RH=0.01).dNeff_bbn
    elseif q == "dNcmb";  return model_observables(m; T_RH=0.01).dNeff_cmb
    else error("unknown quantity $q")
    end
end

# ΔNeff is near zero across this grid; judge on abs OR rel (like validate_m4).
neff_err(a, b) = min(relerr(a, b), (isnan(a) || isnan(b)) ? Inf : abs(a - b))

emax = Dict{String,Float64}()
worst = Dict{String,String}()
# cache one model per (interaction, gtag, mx, mphi, gx) — construction is cheap
# but this also lets a construction error surface clearly.
for r in 1:size(raw, 1)
    ic = String(raw[r, 1]); gt = String(raw[r, 2])
    mx, mphi, gx = raw[r, 3], raw[r, 4], raw[r, 5]
    q = String(raw[r, 6]); arg = raw[r, 7]; ref = raw[r, 8]
    g = gt == "deg" ? [1, 1, 1] / 3 : [0.5, 0.3, 0.2]
    m = model_from_interaction(ic, mx, mphi, gx, g)
    jl = jl_value(m, q, arg)
    e = (q in ("dNbbn", "dNcmb")) ? neff_err(jl, ref) : relerr(jl, ref)
    key = "$ic/$q"
    if e > get(emax, key, 0.0)
        emax[key] = e
        worst[key] = @sprintf("%s %s mx=%.3g mphi=%.3g gx=%.3g arg=%.4g  jl=%.6g ml=%.6g",
                              ic, gt, mx, mphi, gx, arg, jl, ref)
    end
end

println("=== M5 multi-pole models: max relerr per (interaction/quantity) ===")
allmax = 0.0
for key in sort(collect(keys(emax)))
    e = emax[key]
    global allmax = max(allmax, e)
    e > 1e-6 && @printf("  %-22s relerr = %.3e   %s\n", key, e, get(worst, key, ""))
end
@printf("\n(quantities with relerr < 1e-6 omitted)\n")
@printf("OVERALL max relerr = %.3e  -> %s\n", allmax,
        allmax < 1e-2 ? "PASS (<1%)" : (allmax < 5e-2 ? "OK (<5%)" : "FAIL"))

# --- axialVector validity floor: ultralight A' must throw (matches MATLAB's
#     axial_couplings:mediatorTooLight error) ---
function floor_throws()
    try
        AxialVector(0.05, 1e-10, 1e-3, [0.5, 0.3, 0.2])
        return false
    catch e
        return e isa AxialMediatorTooLight
    end
end
@printf("axialVector validity floor throws AxialMediatorTooLight: %s\n",
        floor_throws() ? "PASS" : "FAIL")
