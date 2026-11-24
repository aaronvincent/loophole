#include "julia_worker.hpp"

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cerrno>
#include <csignal>
#include <stdexcept>
#include <vector>
#include <unistd.h>
#include <sys/wait.h>
#include <spawn.h>

// Access to the process environment for posix_spawn's envp.  On macOS a shared
// library must reach `environ` through _NSGetEnviron() (a direct reference does
// not link); elsewhere the global symbol is fine.
#if defined(__APPLE__)
#  include <crt_externs.h>
#  define LOOP_ENVIRON (*_NSGetEnviron())
#else
extern char** environ;
#  define LOOP_ENVIRON environ
#endif

namespace loop_like {

namespace {
    std::string env_or(const char* name, const std::string& dflt) {
        const char* v = std::getenv(name);
        return (v && *v) ? std::string(v) : dflt;
    }
    // Read one '\n'-terminated line from a FILE*; returns false at EOF/error.
    bool read_line(FILE* f, std::string& out) {
        out.clear();
        char* buf = nullptr;
        size_t cap = 0;
        ssize_t n = ::getline(&buf, &cap, f);
        if (n < 0) { std::free(buf); return false; }
        while (n > 0 && (buf[n-1] == '\n' || buf[n-1] == '\r')) --n;
        out.assign(buf, n);
        std::free(buf);
        return true;
    }
}  // namespace

JuliaWorker& JuliaWorker::instance() {
    static JuliaWorker w;       // lazy, thread-safe init (C++11)
    return w;
}

JuliaWorker::JuliaWorker() {
    julia_bin_ = env_or("LOOP_JULIA_BIN", "julia");
    project_   = env_or("LOOP_JULIA_PROJECT", "");
    worker_jl_ = env_or("LOOP_WORKER_JL", "");
    sysimage_  = env_or("LOOP_JULIA_SYSIMAGE", "");   // optional PackageCompiler sysimage
    // Writing to a dead worker must return EPIPE, not kill the GAMBIT process.
    std::signal(SIGPIPE, SIG_IGN);
}

JuliaWorker::~JuliaWorker() { shutdown(); }

bool JuliaWorker::spawn() {
    if (project_.empty() || worker_jl_.empty()) {
        std::fprintf(stderr,
            "[Loophole] LOOP_JULIA_PROJECT and LOOP_WORKER_JL must be set "
            "(see README.md).\n");
        return false;
    }

    int in_pipe[2], out_pipe[2];     // parent->child stdin, child->parent stdout
    if (::pipe(in_pipe) != 0 || ::pipe(out_pipe) != 0) {
        std::fprintf(stderr, "[Loophole] pipe() failed\n");
        return false;
    }

    // Build argv: julia [--sysimage=…] --startup-file=no --project=… worker.jl
    std::string proj = "--project=" + project_;
    std::vector<std::string> args = {julia_bin_};
    if (!sysimage_.empty()) args.push_back("--sysimage=" + sysimage_);
    args.push_back("--startup-file=no");
    args.push_back(proj);
    args.push_back(worker_jl_);
    std::vector<char*> argv;
    for (auto& a : args) argv.push_back(const_cast<char*>(a.c_str()));
    argv.push_back(nullptr);

    // Build envp = the parent environment + JULIA_NUM_THREADS=1 (added only if
    // unset, so it stays overridable) — defaulting to a single-threaded worker so
    // N MPI ranks don't each spawn a thread pool and oversubscribe the node.  We
    // construct envp explicitly rather than calling setenv() so we never mutate
    // the (multithreaded) parent's environment.
    std::vector<std::string> envs;
    bool have_threads = false;
    for (char** e = LOOP_ENVIRON; e && *e; ++e) {
        envs.emplace_back(*e);
        if (std::strncmp(*e, "JULIA_NUM_THREADS=", 18) == 0) have_threads = true;
    }
    if (!have_threads) envs.emplace_back("JULIA_NUM_THREADS=1");
    std::vector<char*> envp;
    for (auto& s : envs) envp.push_back(const_cast<char*>(s.c_str()));
    envp.push_back(nullptr);

    // Spawn via posix_spawn rather than fork()/exec(): GAMBIT calls the
    // likelihood from inside an OpenMP parallel region, and fork() in a
    // multithreaded process followed by malloc (here, the argv/exec setup) before
    // exec is undefined and corrupts the heap.  posix_spawn is async-safe for
    // multithreaded callers and does the dup2/close + exec atomically.  The child
    // gets a fresh image, so any copied MPI/network state is discarded before use.
    posix_spawn_file_actions_t fa;
    ::posix_spawn_file_actions_init(&fa);
    ::posix_spawn_file_actions_adddup2(&fa, in_pipe[0],  STDIN_FILENO);
    ::posix_spawn_file_actions_adddup2(&fa, out_pipe[1], STDOUT_FILENO);
    ::posix_spawn_file_actions_addclose(&fa, in_pipe[0]);
    ::posix_spawn_file_actions_addclose(&fa, in_pipe[1]);
    ::posix_spawn_file_actions_addclose(&fa, out_pipe[0]);
    ::posix_spawn_file_actions_addclose(&fa, out_pipe[1]);

    pid_t pid = -1;
    int rc = ::posix_spawnp(&pid, julia_bin_.c_str(), &fa, nullptr,
                            argv.data(), envp.data());
    ::posix_spawn_file_actions_destroy(&fa);

    // --- parent: keep our ends, close the child ends ---
    ::close(in_pipe[0]);
    ::close(out_pipe[1]);

    if (rc != 0) {
        std::fprintf(stderr, "[Loophole] posix_spawnp(%s) failed: %s\n",
                     julia_bin_.c_str(), std::strerror(rc));
        ::close(in_pipe[1]);
        ::close(out_pipe[0]);
        return false;
    }

    to_child_   = ::fdopen(in_pipe[1], "w");
    from_child_ = ::fdopen(out_pipe[0], "r");
    pid_ = pid;
    if (!to_child_ || !from_child_) {
        std::fprintf(stderr, "[Loophole] fdopen() failed\n");
        shutdown();
        return false;
    }

    // Wait for the worker's READY handshake (it prints precompile noise to
    // stderr, so the first stdout line is READY).
    std::string line;
    while (read_line(from_child_, line)) {
        if (line == "READY") { alive_ = true; return true; }
    }
    std::fprintf(stderr, "[Loophole] worker exited before READY\n");
    shutdown();
    return false;
}

void JuliaWorker::shutdown() {
    if (to_child_) {
        std::fprintf(to_child_, "quit\n");
        std::fflush(to_child_);
        std::fclose(to_child_);
        to_child_ = nullptr;
    }
    if (from_child_) { std::fclose(from_child_); from_child_ = nullptr; }
    if (pid_ > 0) { int st = 0; ::waitpid(pid_, &st, 0); pid_ = -1; }
    alive_ = false;
}

Observables JuliaWorker::parse(const std::string& reply) {
    Observables o;
    if (reply.rfind("OK ", 0) == 0) {
        // OK Oh2 dNeff_bbn dNeff_cmb sigma_chiN eft_valid sin_mix gx_ceiling solver
        double v[8];
        int n = std::sscanf(reply.c_str() + 3, "%lf %lf %lf %lf %lf %lf %lf %lf",
                            &v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &v[6], &v[7]);
        if (n == 8) {
            o.ok = true;
            o.Oh2 = v[0]; o.dNeff_bbn = v[1]; o.dNeff_cmb = v[2]; o.sigma_chiN = v[3];
            o.eft_valid = v[4]; o.sin_mix = v[5]; o.gx_ceiling = v[6]; o.solver = v[7];
        } else {
            o.error = "unparsable OK reply: " + reply;
        }
    } else if (reply.rfind("ERR", 0) == 0) {
        o.error = (reply.size() > 4) ? reply.substr(4) : "worker error";
    } else {
        o.error = "unexpected reply: " + reply;
    }
    return o;
}

Observables JuliaWorker::query(const std::string& request) {
    if (!alive_ && !spawn()) {
        Observables o; o.error = "worker spawn failed"; return o;
    }
    // Send request; read one reply line.  On IO failure, respawn once and retry.
    for (int attempt = 0; attempt < 2; ++attempt) {
        if (std::fprintf(to_child_, "%s\n", request.c_str()) < 0 ||
            std::fflush(to_child_) != 0) {
            shutdown();
            if (attempt == 0 && spawn()) continue;
            Observables o; o.error = "worker write failed"; return o;
        }
        std::string reply;
        if (!read_line(from_child_, reply)) {
            shutdown();
            if (attempt == 0 && spawn()) continue;
            Observables o; o.error = "worker died (no reply)"; return o;
        }
        return parse(reply);
    }
    Observables o; o.error = "worker unrecoverable"; return o;
}

Observables JuliaWorker::eval(const std::string& interaction,
                              double mx, double mphi, double gx,
                              double gu, double gd, double gs,
                              double T_RH, double s_mix) {
    char buf[512];
    std::snprintf(buf, sizeof buf, "%s %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g",
                  interaction.c_str(), mx, mphi, gx, gu, gd, gs, T_RH, s_mix);
    std::string request(buf);

    std::lock_guard<std::mutex> lock(mu_);
    auto it = cache_.find(request);
    if (it != cache_.end()) return it->second;

    Observables o = query(request);
    if (cache_.size() >= CACHE_MAX) cache_.clear();   // simple bounded cache
    cache_.emplace(request, o);
    return o;
}

}  // namespace loop_like
