# φ ↔ photon conversion on the e± plasma through the φγγ vertex: Primakoff
# e±γ → e±φ (Debye-screened) and e⁺e⁻ → γφ.  Port of matlab/primakoff_phi_xsec.m
# and matlab/primakoff_phi_rate.m (formulae there and in
# docs/phi_thermalisation_check.md).  Used by neff_simple and freeze_coupled as
# part of φ's SM link.

using QuadGK: quadgk
using SpecialFunctions: besselk
# K₁/K₂ (scaled) via thermal_average.jl's Bessels.jl wrappers: SpecialFunctions'
# besselkx(ν::Int, x) takes the complex AMOS path (~340 ns vs 20 ns).

const _PR_ALPHA = 1 / 137.036
const _PR_ME = 0.51099895e-3

"""
    primakoff_phi_xsec(s, mphi, D) -> (σ_C, σ_A)  [GeV⁻²]

σ(e γ → e φ) and σ(e⁺e⁻ → γ φ) at g_φγγ = 1 GeV⁻¹, massless electrons,
Debye mass² `D`.  σ_C = α I/(16 s²) (closed form for X/D ≥ 10⁻², X = s − m_φ²;
three-term series in X/D below); σ_A = α (1 − m_φ²/s)³/24.
"""
function primakoff_phi_xsec(s::Real, mphi::Real, D::Real)
    s > 0 || return 0.0, 0.0
    X = max(s - mphi^2, 0.0)
    I = 0.0
    if X > 0
        if X / D >= 1e-2
            B = X + D;  C = s^2 + B^2
            I = (X^2 + 2D * X) / 2 - (2B + D) * X + (C + 2B * D) * log1p(X / D) -
                C * X / (D + X)
        else
            m0 = s^2 * X^2 / 2 + X^4 / 12
            m1 = s^2 * X^3 / 3 + X^5 / 30
            m2 = s^2 * X^4 / 4 + X^6 / 60
            I = (m0 - 2m1 / D + 3m2 / D^2) / D^2
        end
    end
    sC = _PR_ALPHA * I / (16 * s^2)
    sA = s > mphi^2 ? _PR_ALPHA * (1 - mphi^2 / s)^3 / 24 : 0.0
    return sC, sA
end

const _PR_CACHE = Dict{Tuple{Float64,Float64},NTuple{2,Float64}}()
const _PR_LOCK = ReentrantLock()

function _primakoff_unit(T::Float64, mphi::Float64)
    key = (T, mphi)
    v = lock(_PR_LOCK) do
        get(_PR_CACHE, key, (NaN, NaN))
    end
    isnan(v[1]) || return v
    Gc = 0.0;  Ga = 0.0
    if T >= _PR_ME / 60
        ze = _PR_ME / T
        fe = ze < 1e-3 ? 1.0 : ze^2 * besselk(2, ze) / 2
        D  = (4pi * _PR_ALPHA / 3) * T^2 * fe
        xth = mphi
        w  = x -> 2x^4 * _besselkx1(x / T) * exp(-(x - xth) / T)
        pts = (xth, (xth .+ T .* (0.01, 0.1, 1.0, 5.0))..., xth + 80T)
        Ic = first(quadgk(x -> w(x) * primakoff_phi_xsec(x^2, mphi, D)[1], pts...;
                          rtol=1e-8, atol=0.0))
        Ia = first(quadgk(x -> w(x) * primakoff_phi_xsec(x^2, mphi, D)[2], pts...;
                          rtol=1e-8, atol=0.0))
        z = mphi / T
        n_scaled = z < 1e-4 ? T^3 / pi^2 : mphi^2 * T * _besselkx2(z) / (2pi^2)
        Gc = 8fe    * T / (32pi^4) * Ic / n_scaled
        Ga = 4fe^2  * T / (32pi^4) * Ia / n_scaled
    end
    lock(_PR_LOCK) do
        length(_PR_CACHE) > 200_000 && empty!(_PR_CACHE)
        _PR_CACHE[key] = (Gc, Ga)
    end
    return Gc, Ga
end

"""
    primakoff_phi_rate(T, mphi, g2) -> Γ [GeV]

Per-φ rate `R/n_φ^eq` of e±γ ↔ e±φ plus e⁺e⁻ ↔ γφ, g2 = g_φγγ² [GeV⁻²]
(Γ_γγ = g_φγγ² m_φ³/(64π)).  Electron mass only through f_e = z²K₂(z)/2 on the
scatterer densities and the Debye mass; zero for T < m_e/60.
"""
function primakoff_phi_rate(T::Real, mphi::Real, g2::Real)
    g2 > 0 || return 0.0
    Gc, Ga = _primakoff_unit(Float64(T), Float64(mphi))
    return g2 * (Gc + Ga)
end
