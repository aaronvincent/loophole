# neff_simple (simple CMB ΔNeff) vs MATLAB.  Regenerate the reference with
# test/matlab/dump_simple.m.
using Loophole
using DelimitedFiles
using Printf

const REF_SIMPLE = joinpath(@__DIR__, "ref", "simple_ref.csv")
raw = readdlm(REF_SIMPLE, ','; skipstart=1)

# ΔNeff is judged on an absolute scale (it passes through 0); weights likewise.
absdiff(a, b) = (isnan(a) && isnan(b)) ? 0.0 : ((isnan(a) || isnan(b)) ? Inf : abs(a - b))

keys_ = ("dNeff", "dNeff_heat", "dNeff_phi_rad", "w_phi", "w_chi")
emax  = Dict(k => 0.0 for k in keys_)
worst = Dict{String,String}()
for r in 1:size(raw, 1)
    ic = String(raw[r, 1])
    mx, mphi, gx = raw[r, 2], raw[r, 3], raw[r, 4]
    g = [raw[r, 5], raw[r, 6], raw[r, 7]]
    T_RH = raw[r, 8]
    ref = raw[r, 9:13]
    jl = try
        m = model_from_interaction(ic, mx, mphi, gx, g)
        d, I = neff_simple(m; T_RH=T_RH)
        [d, I.dNeff_heat, I.dNeff_phi_rad, I.w_phi, I.w_chi]
    catch
        fill(NaN, 5)
    end
    tag = @sprintf("%s mx=%.3g mphi=%.3g gx=%.3g g1=%.3g", ic, mx, mphi, gx, g[1])
    for (k, key) in enumerate(keys_)
        e = absdiff(jl[k], ref[k])
        if e > emax[key]
            emax[key] = e
            worst[key] = @sprintf("%s  jl=%.6g ml=%.6g", tag, jl[k], ref[k])
        end
    end
end

println("=== validate_simple: neff_simple Julia vs MATLAB ($(size(raw,1)) points) ===")
const TOL_SIMPLE = 2e-3          # absolute, on ΔNeff and on the weights
npass = 0
for key in keys_
    ok = emax[key] <= TOL_SIMPLE
    global npass += ok
    @printf("  %-14s max |diff| = %.3e  %s\n", key, emax[key], ok ? "PASS" : "FAIL")
    haskey(worst, key) && println("      worst: ", worst[key])
end
println("$(npass)/$(length(keys_)) passed")
println("validate_simple: ", npass == length(keys_) ? "PASS" : "FAIL")
