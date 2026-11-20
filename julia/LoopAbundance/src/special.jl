# Special-function helpers shared across solvers.
#
# Wraps SpecialFunctions.besselk / besselkx and centralizes the asymptotic /
# limiting forms that the MATLAB code applies in several places to avoid
# besselk underflow (x > 30) and the massless-mediator z → 0 limit
# (0·Inf = NaN).  Keep these limits exactly — four besselk/(1/mphi)
# sites depend on them (docs/mediator_width.md).

using SpecialFunctions: besselk, besselkx

"""
    k2_of(x) -> K₂(x)

Modified Bessel K₂ with the same x > 30 asymptotic switchover used in
freeze_loop.m's `Y_eq_of` (error < 1e-5 already at x = 30).
"""
function k2_of(x::Real)
    if x > 30
        return exp(-x + 0.5 * log(pi / (2x))) * (1 + 15 / (8x))
    else
        return besselk(2, x)
    end
end

"""
    z2k2(z) -> z² K₂(z)

The combination entering the massive-species equilibrium yield Y_eq.  Limits:
z → 0 (relativistic) → 2; z > 30 → scaled asymptotic to avoid underflow.
Mirrors the guards in freeze_coupled.m `Y_eq_mass`.
"""
function z2k2(z::Real)
    if z < 1e-3
        return 2.0                      # relativistic limit z²K₂ → 2
    elseif z > 30
        return z^2 * exp(-z + 0.5 * log(pi / (2z))) * (1 + 15 / (8z))
    else
        return z^2 * besselk(2, z)
    end
end

"""
    k1_over_k2(z) -> K₁(z)/K₂(z)

Thermal-average decay factor.  z → 0 limit is 0 (a massless mediator gives
0, not NaN).  Uses scaled Bessels so
the ratio is stable at large z.
"""
function k1_over_k2(z::Real)
    if z == 0
        return 0.0
    elseif z > 1e6
        # Amos argument reduction throws / loses all accuracy past ~1e8 even for
        # the scaled Bessels (MATLAB silently returns the ratio → 1).  Use the
        # large-z asymptotic K_ν ~ √(π/2z)e^{-z}(1+(4ν²-1)/8z): error ~1/z² here.
        return (1 + 3 / (8z)) / (1 + 15 / (8z))
    elseif z > 30
        # K₁/K₂ → 1 - 3/(2z) + ... ; use scaled Bessels for stability
        return besselkx(1, z) / besselkx(2, z)
    else
        return besselk(1, z) / besselk(2, z)
    end
end
