// GAMBIT-light user library for the Loophole dark-matter relic / Delta-Neff
// model.  One registered log-likelihood per interaction (selected via `func_name`
// in the YAML); each computes the observables through a persistent Julia worker
// and returns the v1 likelihood (relic Omega h^2 + Delta Neff).
//
// Inputs (named, from UserModel): mx, mphi, gx, gu, gd, gs.
// Optional tunables (pass as `fixed_value` params; defaults below):
//   T_RH=0.01, s_mix=0.3,
//   relic_mode=1 (0 Gaussian / 1 upper-limit), relic_mu=0.12, relic_sigma=0.012,
//   BBN  Delta-Neff: neff_bbn_mode=0 (two-sided Gaussian), neff_bbn_mu=-0.14,
//                    neff_bbn_sigma=0.21  (arXiv:2401.15054),
//   CMB  Delta-Neff: neff_cmb_mode=0 (two-sided Gaussian), neff_cmb_mu=-0.054,
//                    neff_cmb_sigma=0.17  (Planck N_eff=2.99+/-0.17 vs SM 3.044).
// Outputs (added to the output map, list under `output:` to print):
//   Oh2, dNeff_bbn, dNeff_cmb, sigma_chiN, eft_valid, sin_mix, gx_ceiling, solver.

#include "gambit_light_interface.h"
#include "julia_worker.hpp"
#include "likelihoods.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Look up a named input; throw a clear runtime_error if a required one is absent.
double require(const std::map<std::string,double>& in, const char* key) {
    auto it = in.find(key);
    if (it == in.end())
        throw std::runtime_error(std::string("missing required input '") + key + "'");
    return it->second;
}
double opt(const std::map<std::string,double>& in, const char* key, double dflt) {
    auto it = in.find(key);
    return (it == in.end()) ? dflt : it->second;
}

// Shared core: evaluate one interaction at one point.
double loop_core(const std::string& interaction,
                 const std::vector<std::string>& names,
                 const std::vector<double>& vals,
                 std::map<std::string,double>& output)
{
    std::map<std::string,double> in;
    for (size_t i = 0; i < names.size(); ++i) in[names[i]] = vals[i];

    const double mx   = require(in, "mx");
    const double mphi = require(in, "mphi");
    const double gx   = require(in, "gx");
    const double gu   = require(in, "gu");
    const double gd   = require(in, "gd");
    const double gs   = require(in, "gs");

    const double T_RH  = opt(in, "T_RH", 0.01);
    const double s_mix = opt(in, "s_mix", 0.3);

    loop_like::Observables obs = loop_like::JuliaWorker::instance().eval(
        interaction, mx, mphi, gx, gu, gd, gs, T_RH, s_mix);

    if (!obs.ok) {
        gambit_light_invalid_point(obs.error.c_str());
        return -std::numeric_limits<double>::infinity();
    }

    // Always expose the observables (for printing and any further YAML likelihoods).
    output["Oh2"]        = obs.Oh2;
    output["dNeff_bbn"]  = obs.dNeff_bbn;
    output["dNeff_cmb"]  = obs.dNeff_cmb;
    output["sigma_chiN"] = obs.sigma_chiN;
    output["eft_valid"]  = obs.eft_valid;
    output["sin_mix"]    = obs.sin_mix;
    output["gx_ceiling"] = obs.gx_ceiling;
    output["solver"]     = obs.solver;

    // Points outside the EFT's domain of validity are unphysical:
    //   axialVector  -- above the unitarity/meson-mixing ceiling, or sin(theta) >= s_mix;
    //   SP / PP      -- max_j |theta_j| >= s_mix, i.e. the pseudoscalar mediator's
    //                   mixing with pi0/eta/eta' is no longer a small-angle expansion
    //                   (every coupling in that sector is first order in theta).
    // S(x)S and P(x)S go through the chiral loop factor D and are always valid.
    if (obs.eft_valid == 0.0) {
        gambit_light_invalid_point("EFT-invalid point (unitarity ceiling / large meson mixing)");
        return -std::numeric_limits<double>::infinity();
    }
    // A non-finite relic means the solver did not converge here.
    if (!std::isfinite(obs.Oh2)) {
        gambit_light_invalid_point("relic solver returned non-finite Omega h^2");
        return -std::numeric_limits<double>::infinity();
    }

    const int    relic_mode  = static_cast<int>(opt(in, "relic_mode", 1));
    const double relic_mu    = opt(in, "relic_mu", 0.12);
    const double relic_sigma = opt(in, "relic_sigma", 0.012);
    // Both epochs use a two-sided Gaussian on the measured Delta Neff:
    // BBN -0.14 +/- 0.21 (arXiv:2401.15054); CMB -0.054 +/- 0.17 (Planck
    // N_eff = 2.99 +/- 0.17 vs SM 3.044).  Two-sided is required because the
    // model's CMB Delta Neff is signed and goes negative (photon heating from
    // late phi decay) -- a one-sided upper cap would never penalize it.
    const int    bbn_mode    = static_cast<int>(opt(in, "neff_bbn_mode", 0));
    const double bbn_mu      = opt(in, "neff_bbn_mu", -0.14);
    const double bbn_sigma   = opt(in, "neff_bbn_sigma", 0.21);
    const int    cmb_mode    = static_cast<int>(opt(in, "neff_cmb_mode", 0));
    const double cmb_mu      = opt(in, "neff_cmb_mu", -0.054);
    const double cmb_sigma   = opt(in, "neff_cmb_sigma", 0.17);

    double loglike = loop_like::relic_loglike(obs.Oh2, relic_mode, relic_mu, relic_sigma);
    loglike += loop_like::neff_loglike(obs.dNeff_bbn, bbn_mode, bbn_mu, bbn_sigma);
    loglike += loop_like::neff_loglike(obs.dNeff_cmb, cmb_mode, cmb_mu, cmb_sigma);
    return loglike;
}

}  // namespace

// --- one registered log-likelihood per interaction (pick via YAML func_name) ---
#define LOOP_LOGLIKE(NAME)                                                        \
    double loop_##NAME(const std::vector<std::string>& names,                     \
                       const std::vector<double>& vals,                           \
                       std::map<std::string,double>& output)                      \
    { return loop_core(#NAME, names, vals, output); }                             \
    GAMBIT_LIGHT_REGISTER_LOGLIKE(loop_##NAME)

LOOP_LOGLIKE(heavyScalar)
LOOP_LOGLIKE(lightScalar)
LOOP_LOGLIKE(heavyScalarPS)
LOOP_LOGLIKE(lightScalarPS)
LOOP_LOGLIKE(heavyScalarSP)
LOOP_LOGLIKE(lightScalarSP)
LOOP_LOGLIKE(scalarPP)
LOOP_LOGLIKE(heavyScalarPP)
LOOP_LOGLIKE(lightScalarPP)
LOOP_LOGLIKE(axialVector)
