#include "likelihoods.hpp"
#include <cmath>
#include <limits>

namespace loop_like {

double relic_loglike(double Oh2, int mode, double mu, double sigma)
{
    // Non-finite Oh2 (solver failure) is handled as an invalid point upstream;
    // guard defensively so a stray NaN never propagates into the scanner.
    if (!std::isfinite(Oh2))
        return -std::numeric_limits<double>::infinity();
    if (sigma <= 0.0)
        return 0.0;

    const double d = (Oh2 - mu) / sigma;
    if (mode == 1 && Oh2 <= mu)   // upper limit: no penalty below the target
        return 0.0;
    return -0.5 * d * d;          // Gaussian (mode 0), or the >mu tail (mode 1)
}

double neff_loglike(double dNeff, int mode, double mu, double sigma)
{
    if (!std::isfinite(dNeff))    // missing/failed Neff -> no contribution
        return 0.0;
    if (sigma <= 0.0)
        return 0.0;
    if (mode == 1 && dNeff <= mu) // one-sided cap: no penalty below the bound
        return 0.0;
    const double d = (dNeff - mu) / sigma;
    return -0.5 * d * d;          // two-sided Gaussian (mode 0), or the >mu tail (mode 1)
}

}  // namespace loop_like
