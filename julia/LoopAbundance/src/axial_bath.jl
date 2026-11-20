# axialVector A′ ↔ bath links for freeze_coupled (docs/notes/AA_to_gammagamma.tex):
#   AA_sm_xsec      A′A′ → γγ (charged π/K loop) and A′A′ → π⁺π⁻, K⁺K⁻, K⁰K̄⁰ (tree)
#   sv_AAsm_grid    log ⟨σv⟩ of the two, Gondolo–Gelmini, on freeze_coupled's T grid
#   apipi_phi_rate  per-A′ rate of A′π ↔ ππ (contact of the 3π current, crossed)
# Ports of AA_sm_xsec.m, sv_AAsm_grid.m, apipi_phi_rate.m.

# Charged-scalar loop function A0(τ) = (f(τ) − τ)/τ², A0(0) = 1/3; complex above τ = 1.
function _A0_scalar(tau::Real)
    tau < 1e-3 && return complex(1 / 3 + 8tau / 45 + 4tau^2 / 35)
    if tau <= 1
        f = asin(sqrt(tau))^2
        return complex((f - tau) / tau^2)
    end
    r = sqrt(1 - 1 / tau)
    f = -0.25 * (log((1 + r) / (1 - r)) - im * pi)^2
    return (f - tau) / tau^2
end

"""
    AA_sm_xsec(s, mA, g) -> (σ_γγ, σ_MM)   [GeV⁻²]

A′A′ → γγ through the charged-meson loop and A′A′ → meson pairs at tree level,
from L_aa = −Σ_{i<j}(a_i − a_j)² A′·A′ M_ij M_ji.  P(s) = 2 + (s − 2mA²)²/(4mA⁴);
σ_γγ = |c|² s P/(9π β_A), σ_MM = Δa⁴ P β_M/(36π s β_A).  Port of AA_sm_xsec.m.
"""
function AA_sm_xsec(s::Real, mA::Real, g)
    (mA > 0 && s > 4mA^2) || return 0.0, 0.0
    alpha = 1 / 137.036
    bA = sqrt(1 - 4mA^2 / s)
    P = 2 + (s - 2mA^2)^2 / (4mA^4)
    Dpi = g[1] - g[2]; DK = g[1] - g[3]; DK0 = g[2] - g[3]
    c = alpha / (16pi) * (Dpi^2 * _A0_scalar(s / (4MPI_PM^2)) / MPI_PM^2 +
                          DK^2 * _A0_scalar(s / (4MK^2)) / MK^2)
    sgg = abs2(c) * s * P / (9pi * bA)
    smm = 0.0
    for (D, m) in ((Dpi, MPI_PM), (DK, MK), (DK0, MK0))
        (D != 0 && s > 4m^2) || continue
        smm += D^4 * P * sqrt(1 - 4m^2 / s) / (36pi * s * bA)
    end
    return sgg, smm
end

"""
    sv_AAsm_grid(mA, g, Tgrid) -> Vector  (natural log ⟨σv⟩, −Inf where zero)

Gondolo–Gelmini ⟨σv⟩ of A′A′ → SM, ∫σ(s−4mA²)√s K₁(√s/T) ds/(8mA⁴ T K₂²), piecewise
between the thresholds 2mA, 2m_π±, 2m_K±, 2m_K0 (√s = b + uT, u on 400 log nodes,
trapezoid).  Port of sv_AAsm_grid.m.
"""
function sv_AAsm_grid(mA::Real, g, Tgrid::AbstractVector)
    nT = length(Tgrid)
    lsv = fill(-Inf, nT)
    mA > 0 || return lsv
    u = vcat(1e-10, 10.0 .^ range(-9, log10(80); length=399))
    mth = [2MPI_PM, 2MK, 2MK0]
    b = vcat(2mA, sort(mth[mth .> 2mA]))
    lz2k2(z) = z > 30 ? 2log(z) - z + 0.5log(pi / (2z)) + log1p(15 / (8z)) :
               z < 1e-3 ? log(2.0) : log(z^2 * besselk(2, z))
    F = similar(u)
    for (k, T) in enumerate(Tgrid)
        I = 0.0
        for kb in eachindex(b)
            b_hi = kb < length(b) ? b[kb+1] : Inf
            @inbounds for (i, ui) in enumerate(u)
                rs = b[kb] + ui * T; s = rs^2
                if rs >= b_hi
                    F[i] = 0.0
                else
                    sgg, smm = AA_sm_xsec(s, mA, g)
                    F[i] = (sgg + smm) * (s - 4mA^2) * rs * besselkx(1, rs / T) * exp(-ui) * 2rs * T
                end
            end
            It = 0.0
            @inbounds for i in 1:length(u)-1
                It += 0.5 * (F[i] + F[i+1]) * (u[i+1] - u[i])
            end
            I += It * exp(-(b[kb] - 2mA) / T)
        end
        (I > 0 && isfinite(I)) || continue
        lsv[k] = log(I) - 2mA / T - log(8) - 5log(T) - 2lz2k2(mA / T)
    end
    return lsv
end

# --- A′π → ππ ----------------------------------------------------------------
_lamk(s, a, b) = (s - (sqrt(a) + sqrt(b))^2) * (s - (sqrt(a) - sqrt(b))^2)

# σ at Δa_eff = 1: A′π⁰ → π⁺π⁻ (isotropic) and A′π± → π±π⁰ (angle-averaged)
function _apipi_sig1(s, mA)
    li = _lamk(s, mA^2, MPI0^2); lf = _lamk(s, MPI_PM^2, MPI_PM^2)
    (li > 0 && lf > 0) || return 0.0
    M2 = 4 / FPI^2 * li / (4mA^2)
    return M2 / 3 * sqrt(lf) / (16pi * s * sqrt(li))
end
function _apipi_sig2(s, mA)
    li = _lamk(s, mA^2, MPI_PM^2); lf = _lamk(s, MPI_PM^2, MPI0^2)
    (li > 0 && lf > 0) || return 0.0
    rs = sqrt(s)
    Ek = (s + mA^2 - MPI_PM^2) / (2rs); E0 = (s + MPI0^2 - MPI_PM^2) / (2rs)
    pi2 = li / (4s); pf2 = lf / (4s)
    M2 = 4 / FPI^2 * (((Ek * E0)^2 + pi2 * pf2 / 3) / mA^2 - MPI0^2)
    return M2 / 3 * sqrt(pf2) / (16pi * s * sqrt(pi2))
end

const _APIPI_CACHE = Dict{Tuple{Float64,Float64},Float64}()
const _APIPI_LOCK = ReentrantLock()

function _apipi_rate_unit(T::Float64, mA::Float64)
    key = (T, mA)
    v = lock(_APIPI_LOCK) do
        get(_APIPI_CACHE, key, NaN)
    end
    isnan(v) || return v
    G = 0.0
    if T >= MPI_PM / 60 && mA > 0
        th1 = max(mA + MPI0, 2MPI_PM); th2 = max(mA + MPI_PM, MPI_PM + MPI0)
        xth = min(th1, th2)
        function gg(sig, th, m2)
            f = x -> sig(x^2, mA) * _lamk(x^2, mA^2, m2^2) * besselkx(1, x / T) *
                     exp(-(x - th) / T) * 2
            pts = (th * (1 + 1e-12), (th .+ T .* (0.01, 0.1, 1.0, 5.0))..., th + 60T)
            return first(quadgk(f, pts...; rtol=1e-7, atol=0.0))
        end
        I1 = gg(_apipi_sig1, th1, MPI0); I2 = gg(_apipi_sig2, th2, MPI_PM)
        R_scaled = 3T / (32pi^4) * (I1 * exp(-(th1 - xth) / T) + 2I2 * exp(-(th2 - xth) / T))
        n_scaled = 3mA^2 * T * besselkx(2, mA / T) / (2pi^2)
        G = R_scaled * exp(-(xth - mA) / T) / n_scaled
    end
    lock(_APIPI_LOCK) do
        length(_APIPI_CACHE) > 200_000 && empty!(_APIPI_CACHE)
        _APIPI_CACHE[key] = G
    end
    return G
end

"""
    apipi_phi_rate(T, mA, g) -> Γ [GeV]

Per-A′ rate of A′π⁰ ↔ π⁺π⁻ and A′π± ↔ π±π⁰ (contact term of the axial 3π current;
the π⁰ pole is ∝ k·ε = 0 on shell), normalised by Δa_eff = m_{π⁰A′}/F as Γ(A′→3π).
`R / n_A′^eq`, Maxwell–Boltzmann, g_A′ = 3.  Zero for T < m_π/60.  Port of apipi_phi_rate.m.
"""
function apipi_phi_rate(T::Real, mA::Real, g)
    m_pi0Ap, _, _ = _axial_mixings((Float64(g[1]), Float64(g[2]), Float64(g[3])))
    return (m_pi0Ap / FPI)^2 * _apipi_rate_unit(Float64(T), Float64(mA))
end
