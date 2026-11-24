# freeze_coupled with φ's 2→2 conversions (Primakoff, πγ→πφ) vs MATLAB, on and
# off (ENV LOOP_PHI_SC=0).  Regenerate the reference with test/matlab/dump_fc_sc.m.
using Loophole
using DelimitedFiles
using Printf

raw = readdlm(joinpath(@__DIR__, "ref", "fc_sc_ref.csv"), ','; skipstart=1)
relerr(a, b) = (isnan(a) && isnan(b)) ? 0.0 : ((isnan(a) || isnan(b)) ? Inf :
               abs(a - b) / max(abs(a), abs(b), 1e-300))
emax = Dict("on" => 0.0, "off" => 0.0);  worst = Dict{String,String}()
for r in 1:size(raw, 1)
    ic = String(raw[r, 1]);  mx, mphi, gx = raw[r, 2], raw[r, 3], raw[r, 4]
    g = [raw[r, 5], raw[r, 6], raw[r, 7]];  T_RH = raw[r, 8]
    m = model_from_interaction(ic, mx, mphi, gx, g)
    for (key, env, col) in (("on", "", 9), ("off", "0", 10))
        ENV["LOOP_PHI_SC"] = env
        jl = freeze_coupled(m; T_RH=T_RH)
        e = relerr(jl, raw[r, col])
        if e > emax[key]
            emax[key] = e
            worst[key] = @sprintf("%s mx=%.3g mphi=%.3g gx=%.3g g1=%.3g  jl=%.6g ml=%.6g",
                                  ic, mx, mphi, gx, g[1], jl, raw[r, col])
        end
    end
end
ENV["LOOP_PHI_SC"] = ""
println("=== validate_fc_sc: freeze_coupled ± φ conversions, Julia vs MATLAB ($(size(raw,1)) points) ===")
ok = true
for key in ("on", "off")
    @printf("  %-4s max relerr = %.3e\n      worst: %s\n", key, emax[key], get(worst, key, "—"))
    global ok &= emax[key] < 5e-2
end
println("validate_fc_sc: ", ok ? "PASS (<5%)" : "FAIL")
