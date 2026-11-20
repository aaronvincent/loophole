// Standalone verification driver — mimics how GAMBIT-light loads and calls the
// user library, so the integration can be checked WITHOUT building the GAMBIT
// executable.  It dlopen's loop_loglike.so, registers one loglike by name,
// calls it on a point, and prints the loglike + observables.  It also provides
// the gambit_light_invalid_point/error/warning symbols that GAMBIT normally
// supplies (resolved at load via dynamic_lookup).
//
// Build:
//   g++ -std=c++17 gambit/test/loop_driver.cpp \
//       -I $GAMBIT_LIGHT_DIR/gambit_light_interface/include -ldl \
//       -o /tmp/loop_driver
//
// Run (export the LOOP_JULIA_* env first, see README.md):
//   /tmp/loop_driver gambit/cpp/loop_loglike.so loop_heavyScalar \
//       <mx> <mphi> <gx> <gu> <gd> <gs> <relic_mode> <relic_mu> <relic_sigma>
//
// Examples:
//   ... loop_heavyScalar 0.1 10 1e-13 0.3333 0.3333 0.3333 0 0.12 0.012   -> loglike -50 (Gaussian, Oh2~0)
//   ... loop_axialVector 0.05 1e-9 1e-10 2e-10 2e-10 2e-10 1 0.12 0.012    -> Oh2=0.092
//   ... loop_axialVector 0.001 1e-2 1.0 1 1 1 1 0.12 0.012                 -> invalid_point (EFT floor)
#include "gambit_light_interface.h"
#include <dlfcn.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>

extern "C" void gambit_light_invalid_point(const char* m) {
    throw std::runtime_error(std::string("invalid_point: ") + m);
}
extern "C" void gambit_light_error(const char* m)   { std::cerr << "[error] "   << m << "\n"; }
extern "C" void gambit_light_warning(const char* m) { std::cerr << "[warning] " << m << "\n"; }

static t_loglike_fcn_cpp g_fn = nullptr;
extern "C" int capture_cb(const char* /*name*/, void* fptr) {
    g_fn = reinterpret_cast<t_loglike_fcn_cpp>(fptr);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0]
                  << " <loop_loglike.so> <func_name> [mx mphi gx gu gd gs relic_mode relic_mu relic_sigma]\n";
        return 2;
    }
    const std::string so   = argv[1];
    const std::string func = argv[2];   // e.g. loop_heavyScalar
    std::vector<std::string> names = {"mx","mphi","gx","gu","gd","gs","relic_mode","relic_mu","relic_sigma"};
    std::vector<double> vals;
    for (int i = 3; i < 3 + (int)names.size() && i < argc; ++i) vals.push_back(std::stod(argv[i]));

    void* h = dlopen(so.c_str(), RTLD_NOW | RTLD_GLOBAL);
    if (!h) { std::cerr << "dlopen: " << dlerror() << "\n"; return 1; }
    std::string sym = "gambit_light_register_loglike_" + func;
    auto regfn = reinterpret_cast<t_fcn_loglike>(dlsym(h, sym.c_str()));
    if (!regfn) { std::cerr << "dlsym " << sym << ": " << dlerror() << "\n"; return 1; }
    regfn(func.c_str(), capture_cb);
    if (!g_fn) { std::cerr << "registration did not capture a function\n"; return 1; }

    std::map<std::string,double> output;
    try {
        double ll = g_fn(names, vals, output);
        std::cout << "loglike = " << ll << "\n";
    } catch (const std::exception& e) {
        std::cout << "loglike = (invalid) " << e.what() << "\n";
    }
    for (auto& kv : output) std::cout << "  " << kv.first << " = " << kv.second << "\n";
    return 0;
}
