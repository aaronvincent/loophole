// Manager for the persistent LoopAbundance Julia worker process.
//
// One worker per OS process (i.e. one per MPI rank): the GAMBIT-light user lib
// spawns `julia --project=<proj> worker.jl` on first use, then exchanges one
// request/reply line per likelihood point over pipes.  Keeps the validated Julia
// untouched and avoids embedding the Julia runtime inside GAMBIT.
//
// Configuration via environment variables (set before launching gambit_light;
// see README.md and the example YAMLs):
//   LOOP_JULIA_BIN       julia executable           (default: "julia" on PATH)
//   LOOP_JULIA_PROJECT   path to julia/LoopAbundance (required)
//   LOOP_WORKER_JL       path to worker/worker.jl    (required)
//   LOOP_JULIA_SYSIMAGE  optional PackageCompiler sysimage (-J…) to cut startup
//                        (recommended at HPC scale; see gambit/build_sysimage.sh)
// The worker is launched single-threaded by default (JULIA_NUM_THREADS=1,
// overridable) so N MPI ranks don't oversubscribe a node.
#pragma once

#include <string>
#include <map>
#include <mutex>
#include <cstdio>
#include <sys/types.h>

namespace loop_like {

struct Observables {
    bool ok = false;            // false => worker ERR / dead / unparsable
    std::string error;
    double Oh2 = 0.0, dNeff_bbn = 0.0, dNeff_cmb = 0.0, sigma_chiN = 0.0;
    double eft_valid = 1.0, sin_mix = 0.0, gx_ceiling = 0.0, solver = 0.0;
};

class JuliaWorker {
public:
    // Process-wide singleton (lazy, thread-safe per C++11 static-init rules).
    static JuliaWorker& instance();

    // Evaluate one model point.  Cached on the request string so the several
    // UserLogLikes registered for the same point only solve once.
    Observables eval(const std::string& interaction,
                     double mx, double mphi, double gx,
                     double gu, double gd, double gs,
                     double T_RH, double s_mix);

    ~JuliaWorker();

private:
    JuliaWorker();
    JuliaWorker(const JuliaWorker&) = delete;
    JuliaWorker& operator=(const JuliaWorker&) = delete;

    bool spawn();                                   // fork/exec the worker; wait for READY
    void shutdown();                                // "quit" + reap the child
    Observables query(const std::string& request);  // raw IO (assumes mu_ held)
    static Observables parse(const std::string& reply);

    std::string julia_bin_, project_, worker_jl_, sysimage_;
    pid_t pid_ = -1;
    FILE* to_child_ = nullptr;     // child stdin  (we write requests)
    FILE* from_child_ = nullptr;   // child stdout (we read replies)
    bool alive_ = false;
    std::mutex mu_;

    std::map<std::string, Observables> cache_;
    static constexpr size_t CACHE_MAX = 256;
};

}  // namespace loop_like
