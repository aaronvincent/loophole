# Precompile workload for a PackageCompiler sysimage (HPC startup optimization).
#
# Exercises the solver / ΔNeff / observable code paths for representative
# interactions so the resulting sysimage has them compiled — cutting the per-rank
# worker startup from ~seconds to ~1s and removing first-call JIT mid-scan.
#
# Build the sysimage with gambit/build_sysimage.sh, then run with
#   export LOOP_JULIA_SYSIMAGE=$PWD/gambit/loopabundance_sys.so

using Loophole

# (interaction, mx, mphi, gx, g) — valid points that don't trip the axial floor.
const _PTS = [
    ("heavyScalar",   0.1,  10.0,  1.0,    [1.0, 1, 1] ./ 3),   # freeze_loop
    ("heavyScalarPS", 0.1,  10.0,  1.0,    [1.0, 1, 1] ./ 3),   # single-pole PS
    ("heavyScalarSP", 0.05, 0.8,   0.1,    [2.0, -1, -1] ./ 3), # multi-pole SP
    ("scalarPP",      0.05, 0.8,   0.1,    [1.0, 1, 1] ./ 3),   # multi-pole PP (η′ singlet)
    ("heavyScalar",   0.005, 1e-4, 1.0,    [1.0, 1, 1] ./ 3),   # freeze_coupled (mφ ≤ 10·T_RH)
    ("axialVector",   0.05, 1e-9,  1e-10,  2e-10 .* [1.0, 1, 1]),# multi-pole axial
]

for (it, mx, mphi, gx, g) in _PTS
    try
        m = model_from_interaction(it, mx, mphi, gx, g)
        model_observables(m)
    catch
        # warmup only — ignore any point that errors (e.g. an EFT-floor throw)
    end
end
