// Likelihood terms for the Loophole GAMBIT-light user library.
// v1: relic Omega h^2 (Gaussian or upper-limit) + dark-sector Delta Neff.
// All thresholds / sigmas / mode are passed in (set as fixed_value params in the
// GAMBIT-light YAML), so they live in config rather than the binary.
#pragma once

namespace loop_like {

// Relic-abundance log-likelihood.
//   mode 0 -> Gaussian:    -0.5*((Oh2-mu)/sigma)^2   (chi = all of DM)
//   mode 1 -> upper limit:  0 if Oh2 <= mu, else -0.5*((Oh2-mu)/sigma)^2
//             (no overproduction; chi may be a sub-component)
double relic_loglike(double Oh2, int mode, double mu, double sigma);

// Delta-Neff log-likelihood (mirrors relic_loglike's mode flag).
//   mode 0 -> two-sided Gaussian:  -0.5*((dNeff-mu)/sigma)^2
//             (mu = measured Delta Neff central value, e.g. BBN -0.14 +/- 0.21
//              from arXiv:2401.15054)
//   mode 1 -> one-sided cap:        0 if dNeff <= mu, else -0.5*((dNeff-mu)/sigma)^2
//             (mu = turn-on point of the half-Gaussian, a loose upper bound)
// Call once per epoch (BBN, CMB) and sum.
double neff_loglike(double dNeff, int mode, double mu, double sigma);

}  // namespace loop_like
