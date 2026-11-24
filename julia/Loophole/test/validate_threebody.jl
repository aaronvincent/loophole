# SP/PP three-body final states (χχ̄→3π, φ→3π, χχ̄→ππφ) and the axialVector
# χχ̄→3π, A′→3π, χχ̄→π+π−A′ vs MATLAB, toggles on, plus the toggles-off checks.  Regenerate the reference with
#   LOOP_3PI=1 LOOP_PIPIPHI=1 matlab -batch "dump_threebody"   (test/matlab/)
using Loophole
using DelimitedFiles
using Printf

ENV["LOOP_3PI"] = "1"; ENV["LOOP_PIPIPHI"] = "1"
threebody_refresh!()

const REF = joinpath(@__DIR__, "ref")
raw = readdlm(joinpath(REF, "threebody_ref.csv"), ','; skipstart=1)

function relerr(a, b)
    (isnan(a) && isnan(b)) && return 0.0
    (isnan(a) || isnan(b)) && return Inf
    (a == 0 && b == 0) && return 0.0
    return abs(a - b) / max(abs(a), abs(b), 1e-300)
end

function jl_value(m, q, arg)
    isSP = m isa Union{HeavyScalarSP,LightScalarSP}
    if     q == "w3pi";   return mediator_widths(m).pipi
    elseif q == "wtot";   return mediator_widths(m).total
    elseif q == "W3pi";   return threepi_weight(arg, m.g)
    elseif q == "Wppf";   return pipiphi_weight(arg, m.mphi, m.g)
    elseif q == "WT";     return axial3pi_weight(arg, m.g)[1]
    elseif q == "WL";     return axial3pi_weight(arg, m.g)[2]
    elseif q == "WTa";    return pipiA_weight(arg, m.mphi, m.g)[1]
    elseif q == "WLa";    return pipiA_weight(arg, m.mphi, m.g)[2]
    elseif q == "sig3" && m isa AxialVector
        return axial3pi_sigma(arg, m.mx, m.mphi, m.gx, mediator_widths(m).total, m.g)
    elseif q == "sig3";   return threebody_sigma(arg, m.mx, m.mphi, m.gx, m.Gam[2], m.g, isSP)
    elseif q == "xxall";  return sigma_xxall(m, arg)
    elseif q == "xxhad";  return Loophole.sigma_xxhad(m, arg)
    elseif q == "sv_all"; return svavg(m, arg)
    elseif q == "sv_sm";  return svavg(m, arg; channel=:sm)
    else error("unknown quantity $q")
    end
end

# fixed-order quadratures, identical in both languages: machine-level agreement.
# svavg is adaptive (different rules in the two languages): looser.
tol(q) = startswith(q, "sv_") ? 1e-4 : 1e-8

emax = Dict{String,Float64}(); worst = Dict{String,String}()
nonzero = Dict{String,Int}()
for r in 1:size(raw, 1)
    ic = String(raw[r, 1]); gt = String(raw[r, 2])
    mx, mphi, gx = raw[r, 3], raw[r, 4], raw[r, 5]
    q = String(raw[r, 6]); arg = raw[r, 7]; ref = raw[r, 8]
    g = gt == "deg" ? [1, 1, 1] / 3 : gt == "ax" ? 0.01 * [2, -1, -1] : [0.5, 0.3, 0.2]
    m = model_from_interaction(ic, mx, mphi, gx, g)
    jl = jl_value(m, q, arg)
    e = relerr(jl, ref)
    ref != 0 && (nonzero[q] = get(nonzero, q, 0) + 1)
    key = "$ic/$q"
    if e >= get(emax, key, -1.0)
        emax[key] = e
        worst[key] = @sprintf("%s mx=%.3g mphi=%.3g gx=%.3g arg=%.4g jl=%.10g ml=%.10g",
                              gt, mx, mphi, gx, arg, jl, ref)
    end
end

ok = true
for key in sort(collect(keys(emax)))
    q = split(key, "/")[2]
    pass = emax[key] <= tol(q)
    global ok &= pass
    @printf("  %-5s %-24s max relerr = %.3e   %s\n", pass ? "PASS" : "FAIL", key, emax[key], worst[key])
end
println("  nonzero reference rows per quantity: ", nonzero)

# toggles off: every three-body piece is exactly zero
ENV["LOOP_3PI"] = "0"; ENV["LOOP_PIPIPHI"] = "0"
threebody_refresh!()
for ic in ("heavyScalarSP", "heavyScalarPP"), g in ([1, 1, 1] / 3, [0.5, 0.3, 0.2])
    m = model_from_interaction(ic, 0.2, 1.0, 1.0, g)
    off = threepi_weight(2.0, g) == 0 && pipiphi_weight(2.0, 0.02, g) == 0 &&
          mediator_widths(m).pipi == 0 && Loophole._sppp_threebody(m, 2.0) == 0 &&
          Loophole._model_extra_thresh(m) == ()
    global ok &= off
    @printf("  %-5s toggles off: %s g=%s all three-body pieces zero\n", off ? "PASS" : "FAIL", ic, g)
end
let g = 0.01 * [2, -1, -1], m = model_from_interaction("axialVector", 0.3, 1.0, 0.01, g)
    off = axial3pi_weight(2.0, g) == (0.0, 0.0) && pipiA_weight(2.0, 0.1, g) == (0.0, 0.0) &&
          mediator_widths(m).pipi == 0 &&
          Loophole._model_extra_thresh(m) == () &&
          length(Loophole._sv_poles(m, :all)) == length(poles(m))
    global ok &= off
    @printf("  %-5s toggles off: axialVector all 3pi, pi pi A' pieces zero, no A' pole\n", off ? "PASS" : "FAIL")
end
println("validate_threebody: ", ok ? "PASS" : "FAIL")
