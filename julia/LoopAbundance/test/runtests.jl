# Single entry point for the MATLAB->Julia port validation suites.
#
#   julia --project=. test/runtests.jl              # all suites
#   julia --project=. test/runtests.jl m1 m5        # a subset
#
# Each suite is a standalone script (test/validate_<tag>.jl) comparing the Julia
# port against a MATLAB reference in test/ref/.  They are run as SUBPROCESSES:
# they use top-level `global` state and m7 calls error() on failure, so importing
# them into one session would let the first failure abort the rest.
#
# Exit code is 1 if any suite fails -- the individual scripts (except m7) only
# print "FAIL" and exit 0, so this runner is what makes them usable from CI.
#
# Reference provenance:
#   m1 m2 m3 m4 m5 m7 simple fc_sc   test/matlab/dump_<tag>.m -> test/ref/<tag>_ref.csv
#   threebody               dump_threebody.m, run with LOOP_3PI=1 LOOP_PIPIPHI=1
# Every suite regenerates the same way: matlab -batch "dump_<tag>" from test/matlab/.
using Printf

const TESTDIR = @__DIR__

# tag => description
const SUITES = [
    "m1"  => "cross sections + widths (SS/PS)",
    "m2"  => "thermal averages (svavg, ggxxrate)",
    "m3"  => "relic solvers (loop/out/FI/coupled)",
    "m4"  => "model_observables (relic + dNeff + sigma)",
    "m5"  => "multi-pole models (SP/PP + axialVector)",
    "m7"  => "h_P, axial EFT validity",
    "simple" => "neff_simple (simple CMB dNeff)",
    "fc_sc" => "freeze_coupled phi conversions (Primakoff + pi gamma), on/off",
    "threebody" => "SP/PP 3pi + pi pi phi, axialVector 3pi (LOOP_3PI/LOOP_PIPIPHI on)",
]

# A suite reports its verdict as "-> PASS/OK/FAIL" (m1-m7) or one leading
# "PASS"/"FAIL" per check.  Collect every
# verdict; worst one wins.
const VERDICT = r"->\s*(PASS|OK|FAIL)|validate_\w+:\s*(PASS|FAIL)|^\s+(PASS|FAIL)\s"m
const METRIC  = r"(?:max relerr|max err|max loglike relerr)[^=]*=\s*([0-9.eE+-]+)"

rank(v) = v == "FAIL" ? 2 : v == "OK" ? 1 : 0
label(r) = r == 2 ? "FAIL" : r == 1 ? "OK" : "PASS"

function run_suite(tag)
    path = joinpath(TESTDIR, "validate_$tag.jl")
    isfile(path) || return (missing, "no such suite", 0.0, String[])
    buf = IOBuffer()
    t0 = time()
    # Don't let a nonzero exit (m7) throw; we classify from the output.
    p = run(pipeline(ignorestatus(`$(Base.julia_cmd()) --project=$(dirname(TESTDIR)) $path`),
                     stdout=buf, stderr=buf))
    dt = time() - t0
    out = String(take!(buf))

    verdicts = String[]
    for m in eachmatch(VERDICT, out)
        push!(verdicts, something(m.captures...))
    end
    metrics = [m.captures[1] for m in eachmatch(METRIC, out)]
    # No verdict line at all => the script died before reporting.
    worst = isempty(verdicts) ? 2 : maximum(rank.(verdicts))
    if isempty(verdicts) && p.exitcode != 0
        return (worst, "crashed (exit $(p.exitcode))", dt, split(out, '\n'))
    end
    return (worst, isempty(metrics) ? "" : "max " * metrics[argmax(parse.(Float64, metrics))],
            dt, split(out, '\n'))
end

tags = isempty(ARGS) ? first.(SUITES) : ARGS
info = Dict(SUITES)

println("\nRunning $(length(tags)) validation suite(s)\n")
results = Dict{String,Any}()
for tag in tags
    print("  $(rpad(tag, 5)) ... ")
    flush(stdout)
    r, metric, dt, lines = run_suite(tag)
    results[tag] = (r, metric, dt, lines)
    @printf("%-5s  %6.1fs  %s\n", r === missing ? "?" : label(r), dt, metric)
end

println("\n" * "="^74)
@printf("%-5s %-40s %-6s %s\n", "SUITE", "COVERS", "STATUS", "NOTE")
println("="^74)
nfail = 0
for tag in tags
    r, metric, dt, _ = results[tag]
    desc = get(info, tag, "?")
    st = r === missing ? "?" : label(r)
    r === 2 && (global nfail += 1)
    @printf("%-5s %-40s %-6s %s\n", tag, desc, st, metric)
end
println("="^74)

if nfail > 0
    println("\nFailing suites -- last 12 lines each:\n")
    for tag in tags
        results[tag][1] === 2 || continue
        println("--- $tag " * "-"^(70 - length(tag)))
        for l in filter(!isempty, results[tag][4])[max(1, end - 11):end]
            println("  ", l)
        end
        println()
    end
end

@printf("%d suite(s) run, %d failing\n", length(tags), nfail)
exit(nfail == 0 ? 0 : 1)
