# Effective relativistic degrees of freedom g_*(T).
#
# Port of gstarofT_vec.m (function name `gstarofT`).  Five-tanh interpolation
# fit from arXiv:0910.1066.  `sflag == 0` returns g_rho (energy density dof);
# any nonzero `sflag` returns g_s (entropy dof).  T in GeV.

"""
    gstarofT(T, sflag) -> g_*

Effective relativistic dof.  `sflag == 0` → g_ρ (energy), nonzero → g_s (entropy).
Scalar `T`; broadcast with `gstarofT.(Tvec, sflag)` for a grid.
"""
function gstarofT(T::Real, sflag::Integer)
    lT = log(T)
    if sflag == 0
        a0 = 1.21
        a1 = (0.572, 0.330, 0.579, 0.138, 0.108)
        a2 = (-8.77, -2.95, -1.80, -0.162, 3.76)
        a3 = (0.682, 1.01, 0.165, 0.934, 0.869)
    else
        a0 = 1.36
        a1 = (0.498, 0.327, 0.579, 0.140, 0.109)
        a2 = (-8.74, -2.89, -1.79, -0.102, 3.82)
        a3 = (0.693, 1.01, 0.155, 0.963, 0.907)
    end
    acc = a0
    @inbounds for k in 1:5
        acc += a1[k] * (1 + tanh((lT - a2[k]) / a3[k]))
    end
    return exp(acc)
end

# Allow Bool / Float sflag transparently (MATLAB passed 0 / 1).
gstarofT(T::Real, sflag) = gstarofT(T, Int(sflag != 0))
