# π±γ → π±φ for a scalar φ–quark coupling (S⊗S / P⊗S): cross section and the
# per-φ thermal rate.  Port of matlab/pigamma_phi_xsec.m and
# matlab/pigamma_phi_rate.m (formulae: docs/phi_thermalisation_check.md).
# Used by neff_simple as part of φ's SM link.

using QuadGK: quadgk
# K₁/K₂ (scaled) via thermal_average.jl's Bessels.jl wrappers (fast real path; see
# primakoff_phi.jl).

const _PG_ALPHA = 1 / 137.036          # as the MATLAB original (not ALPHA_EM = 1/137)
const _PG_THR_SWITCH = 1e-2

function _pg_bounds(s, m, mphi)
    a = s - m^2
    L = sqrt(max((s - (m + mphi)^2) * (s - (m - mphi)^2), 0.0))
    c = (s + m^2) * (s + mphi^2 - m^2)
    bp = mphi^2 - (c - a * L) / (2s)
    bm = mphi^2 - (c + a * L) / (2s)
    return a, bp, bm
end

"""
    pigamma_phi_xsec(s, mphi, g) -> σ(π⁺γ→π⁺φ) [GeV⁻²]

LO χPT, `μ = B0 (g_u + g_d)`, photon-polarisation averaged.  Closed form above
s/s_th − 1 = 10⁻², direct quadrature of Σ|M|² over b below (the closed form
cancels catastrophically near threshold for light φ).
"""
function pigamma_phi_xsec(s::Real, mphi::Real, g)
    m = MPI_PM
    mu = B0 * (g[1] + g[2])
    pref = 4pi * _PG_ALPHA * mu^2 / (8pi)
    sth = (m + mphi)^2
    s > sth || return 0.0
    a, bp, bm = _pg_bounds(s, m, mphi)
    if (s - sth) / sth >= _PG_THR_SWITCH
        I = (bp - bm) * (m^2 / a^2 + 1 / a) + m^2 * (1 / bm - 1 / bp) +
            ((2m^2 - mphi^2) / a + 1) * log(bp / bm)
        return -pref * I / a^2
    end
    bm < bp || return 0.0
    f = b -> -(m^2 * ((a + b) / (a * b))^2 + (a + b) / (a * b) - mphi^2 / (a * b))
    return pref * first(quadgk(f, bm, bp; rtol=1e-10, atol=0.0)) / a^2
end

# Unit-coupling (g_u + g_d = 1) rate, cached per (mphi, T): the rate is exactly
# ∝ (g_u + g_d)² and neff_simple calls it on a fixed T grid.
const _PG_CACHE = Dict{Tuple{Float64,Float64},Float64}()
const _PG_LOCK = ReentrantLock()

function _pigamma_phi_rate_unit(T::Float64, mphi::Float64)
    key = (T, mphi)
    v = lock(_PG_LOCK) do
        get(_PG_CACHE, key, NaN)
    end
    isnan(v) || return v
    m = MPI_PM
    G = 0.0
    if T >= m / 60
        xth = m + mphi
        gu = (1.0, 0.0, 0.0)
        f = x -> pigamma_phi_xsec(x^2, mphi, gu) * (x^2 - m^2)^2 *
                 _besselkx1(x / T) * exp(-(x - xth) / T) * 2
        pts = (xth, (xth .+ T .* (0.01, 0.1, 1.0, 5.0))..., xth + 60T)
        I = first(quadgk(f, pts...; rtol=1e-7, atol=0.0))
        R_scaled = 4T / (32pi^4) * I
        z = mphi / T
        if z < 1e-4
            G = R_scaled * exp(-xth / T) / (T^3 / pi^2)
        else
            G = R_scaled * exp(-(xth - mphi) / T) / (mphi^2 * T * _besselkx2(z) / (2pi^2))
        end
    end
    lock(_PG_LOCK) do
        length(_PG_CACHE) > 200_000 && empty!(_PG_CACHE)
        _PG_CACHE[key] = G
    end
    return G
end

"""
    pigamma_phi_rate(T, mphi, g) -> Γ [GeV]

Per-φ rate of π±γ ↔ π±φ, `R / n_φ^eq` (Maxwell–Boltzmann, Gondolo–Gelmini for
unequal masses, g_π g_γ = 4).  Zero for T < m_π/60.
"""
pigamma_phi_rate(T::Real, mphi::Real, g) =
    (g[1] + g[2])^2 * _pigamma_phi_rate_unit(Float64(T), Float64(mphi))
