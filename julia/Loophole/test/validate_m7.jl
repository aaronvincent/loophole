# M7 validation: the physics merged from `main` into the Julia port — the
# pseudoscalar nucleon coupling h̃_P for SP/PP (DMnucxsec), the axialVector
# EFT-validity ceilings/floors (axial_gx_ceiling / axial_eta_ceiling /
# axial_mixing_floor) and the isolated χχ̄→A′A′ piece (sigma_xxphiphi / svavg_phiphi).
#
# Reference: test/ref/m7_ref.csv, regenerate with
#   matlab -batch "run('test/matlab/dump_m7.m')"
using Loophole
using DelimitedFiles
using Printf

const REF = joinpath(@__DIR__, "ref")

raw = readdlm(joinpath(REF, "m7_ref.csv"), ',', skipstart=1)
const R = Dict{String,Float64}(String(raw[i, 1]) => Float64(raw[i, 2]) for i in 1:size(raw, 1))
ref(k) = haskey(R, k) ? R[k] : error("m7_ref.csv missing key '$k' — regenerate with dump_m7.m")

relerr(a, b) = abs(b) < 1e-300 ? abs(a) : abs(a / b - 1)
npass = 0; nfail = 0
function check(name, got, key; rtol=1e-6)
    r = ref(key)
    e = relerr(got, r)
    ok = e < rtol
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-28s got=%.10g ref=%.10g relerr=%.2e  %s\n", name, got, r, e,
            ok ? "OK" : "FAIL")
end
function check_bool(name, ok)
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-28s %s\n", name, ok ? "OK" : "FAIL")
end

println("=== M7: h̃_P nucleon σ_N (SP/PP), g=[2 -1 -1]/3, mx=0.05, mφ=0.1, gx=1 ===")
g = [2.0, -1.0, -1.0] / 3
sN(it) = sigma_nucleon(model_from_interaction(it, 0.05, 0.1, 1.0, g))   # GeV⁻²
for it in ("heavyScalarPP", "lightScalarPP", "heavyScalarSP", "lightScalarSP", "scalarPP")
    check(it, sN(it), "signuc_$it")
end
# degenerate g → η′-singlet projection zeroes h̃_P → σ_N = 0
gd = [1.0, 1.0, 1.0] / 3
check("heavyScalarPP (deg g → 0)",
      sigma_nucleon(model_from_interaction("heavyScalarPP", 0.05, 0.1, 1.0, gd)),
      "signuc_deg_heavyScalarPP"; rtol=1e-30)

println("=== M7: axialVector EFT ceilings/floors (mφ=mx=0.05, g=[1 1 1]) ===")
check("axial_gx_ceiling", axial_gx_ceiling(0.05, 0.05, [1.0, 1, 1]), "axial_gx_ceiling")
em, mm, jb = axial_mixing_floor(0.05, 0.05, [1.0, 1, 1], 0.3)
check("axial_mixing_floor eta_mix", em, "axial_mixing_floor_eta_mix")
check("axial_mixing_floor mxmix",   mm, "axial_mixing_floor_mxmix")
check_bool("axial_mixing_floor jbind", jb == Int(ref("axial_mixing_floor_jbind")))
check("axial_eta_ceiling", axial_eta_ceiling(0.05, 0.05, 1e-3, [1.0, 1, 1]), "axial_eta_ceiling")

println("=== M7: χχ̄→A′A′ thermal average (svavg_phiphi axialVector) ===")
ma = AxialVector(0.05, 0.02, 1e-3, (1.0, 1.0, 1.0))
check("svavg_phiphi(T=mx/20)", svavg_phiphi(ma, 0.05 / 20), "svavg_phiphi"; rtol=1e-4)

@printf("\nM7 OVERALL: %d passed, %d failed -> %s\n", npass, nfail,
        nfail == 0 ? "PASS" : "FAIL")
nfail == 0 || error("M7 validation failed")
