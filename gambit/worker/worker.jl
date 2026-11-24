#!/usr/bin/env julia
# Persistent Loophole worker for the GAMBIT-light C++ user library.
#
# Loads Loophole once, then serves one observable request per input line over
# a trivial whitespace-delimited stdin/stdout protocol (no JSON dependency).  One
# worker is spawned per MPI rank by the C++ side (cpp/julia_worker.cpp); the few-
# second startup is paid once per scan.
#
# Run standalone:
#   julia --project=julia/Loophole gambit/worker/worker.jl
# then type, e.g.:
#   axialVector 0.05 1e-9 1e-10 2e-10 2e-10 2e-10 0.01 0.3
#
# Protocol
#   request : <interaction> <mx> <mphi> <gx> <gu> <gd> <gs> [<T_RH> [<s_mix>]]
#   reply   : OK <Oh2> <dNeff_bbn> <dNeff_cmb> <sigma_chiN> <eft_valid> <sin_mix> <gx_ceiling> <solver>
#           | ERR <message-with-underscores>
#   control : "quit" / "exit" / EOF  -> the worker exits.
# T_RH defaults to 0.01, s_mix to 0.3.  eft_valid is 1.0/0.0; solver is 0
# (freeze_loop) or 1 (freeze_coupled).  NaN/Inf are printed verbatim and parsed
# on the C++ side (which treats non-finite Oh2 as an invalid point).

using Loophole
using Printf
using Logging

const NOUT = 8   # Oh2 dNeff_bbn dNeff_cmb sigma_chiN eft_valid sin_mix gx_ceiling solver

function handle(line::AbstractString)
    toks = split(strip(line))
    isempty(toks) && return nothing
    (toks[1] in ("quit", "exit")) && return :quit
    if length(toks) < 7
        return "ERR too_few_fields_(need_interaction_mx_mphi_gx_gu_gd_gs)"
    end
    try
        interaction = String(toks[1])
        mx, mphi, gx = parse.(Float64, toks[2:4])
        g  = parse.(Float64, toks[5:7])
        T_RH = length(toks) >= 8 ? parse(Float64, toks[8]) : 0.01
        s_mix = length(toks) >= 9 ? parse(Float64, toks[9]) : 0.3

        m = model_from_interaction(interaction, mx, mphi, gx, g)
        o = model_observables(m; T_RH=T_RH, s_mix=s_mix)

        eft_valid  = o.eft.valid ? 1.0 : 0.0
        solvercode = o.solver === :freeze_loop ? 0.0 : 1.0
        return @sprintf("OK %.10e %.10e %.10e %.10e %.1f %.10e %.10e %.1f",
                        o.Oh2, o.dNeff_bbn, o.dNeff_cmb, o.sigma_chiN,
                        eft_valid, o.eft.sin_mix, o.eft.gx_ceiling, solvercode)
    catch e
        msg = replace(sprint(showerror, e), r"\s+" => "_")
        return "ERR " * msg
    end
end

function main()
    # Drop @warn-level noise (e.g. SciMLBase's "dt forced below epsilon, Aborting"
    # for parameter points the ODE can't integrate — those are handled via the
    # retcode → NaN → invalid_point, so the warning is benign but floods the scan
    # log).  This is pure log routing: it does NOT change any numerics.  (Do NOT
    # pass verbose=false to solve() instead — in this SciML version that flips the
    # abort retcode and turns finite results into NaN.)  Exceptions still surface
    # (caught in handle → ERR) since they are not log messages.
    global_logger(SimpleLogger(stderr, Logging.Error))

    # One-line startup banner to stderr (keeps the stdout protocol clean) so a
    # CPU/sysimage mismatch is visible in the rank log instead of only surfacing
    # as a hard "Unable to find compatible target in cached code image" abort.
    # Sys.CPU_NAME is the microarchitecture the runtime resolved on THIS node; if
    # the sysimage lacks a compatible target the process would already have died
    # before reaching here, so seeing this line means the image loaded cleanly.
    let img = try unsafe_string(Base.JLOptions().image_file) catch; "" end
        println(stderr, "WORKER node_cpu=", Sys.CPU_NAME,
                " sysimage=", isempty(img) ? "<none>" : basename(img),
                " threads=", Threads.nthreads(),
                " 3pi=", Int(threepi_enabled()), " pipiphi=", Int(pipiphi_enabled()))
        flush(stderr)
    end

    # Touch the solvers once so the first real request isn't paying full JIT cost
    # mid-scan (best-effort; ignore any error).
    try
        model_observables(model_from_interaction("heavyScalar", 0.1, 1.0, 1.0, [1.0,1,1]/3))
    catch
    end
    println("READY"); flush(stdout)
    for line in eachline(stdin)
        r = handle(line)
        r === nothing && continue
        r === :quit && break
        println(r); flush(stdout)
    end
end

main()
