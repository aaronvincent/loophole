# Three-particle final states of the pseudoscalar φ–SM vertex (SP/PP): χχ̄→3π,
# φ→3π and χχ̄→ππφ at LO χPT; and of the axial-vector A′–SM vertex: χχ̄→3π,
# A′→3π, χχ̄→π+π−A′.  Ports of threepi_weight.m, pipiphi_weight.m,
# dalitz_moments.m, threepi_sigma.m, axial3pi_weight.m, pipiA_weight.m,
# axial3pi_sigma.m.  Physics:
# docs/cross_sections.md.
#
# Off by default.  LOOP_3PI=1 (3π, both χχ̄→3π and φ/A′→3π) and LOOP_PIPIPHI=1
# (ππφ, π+π−A′) in the process environment turn them on.  The flags are read when the
# module is loaded (`__init__`, so at runtime even from a sysimage), not at
# precompile time; `threebody_refresh!()` re-reads them.

const _THREEPI_ON = Ref(false)
const _PIPIPHI_ON = Ref(false)

"""
    threebody_refresh!()

Re-read `LOOP_3PI` / `LOOP_PIPIPHI` from `ENV` ("1" turns a channel on).
Called at module load; call it again after changing `ENV` in a session.
"""
function threebody_refresh!()
    _THREEPI_ON[] = get(ENV, "LOOP_3PI", "") == "1"
    _PIPIPHI_ON[] = get(ENV, "LOOP_PIPIPHI", "") == "1"
    return (threepi = _THREEPI_ON[], pipiphi = _PIPIPHI_ON[])
end
threepi_enabled() = _THREEPI_ON[]
pipiphi_enabled() = _PIPIPHI_ON[]

_sppp_degenerate(g) = g[1] == g[2] && g[2] == g[3]

# I_n(s) = ∫dE Y(E) X(E)^n over the energy E of the m1 particle; the other two
# have mass m2.  Y = sqrt((E²−m1²)(1 − 4m2²/w)), w = s + m1² − 2E√s,
# X = w − mπ0².  Y has square-root zeros at both ends, so t = (1 − cos θ)/2 and
# the midpoint rule in θ converge spectrally.
function _dalitz_moments(s::Real, m1, m2)
    n = 32
    rs = sqrt(s)
    Elo = m1
    Ehi = (s + m1^2 - 4m2^2) / (2rs)
    I0 = 0.0; I1 = 0.0; I2 = 0.0
    @inbounds for k in 1:n
        th = pi * (k - 0.5) / n
        t  = (1 - cos(th)) / 2
        dt = sin(th) / 2 * pi / n
        E  = Elo + (Ehi - Elo) * t
        w  = s + m1^2 - 2E * rs
        Y  = sqrt(max((E^2 - m1^2) * (1 - 4m2^2 / w), 0.0))
        X  = w - MPI0^2
        wt = (Ehi - Elo) * dt * Y
        I0 += wt; I1 += wt * X; I2 += wt * X^2
    end
    return I0, I1, I2
end

# W(s) = Σ_ch S ∫dE1 dE2 |A|² over π+π−π0 (S = 1) and 3π0 (S = 1/6), mediator
# propagator stripped.  σ = gx² W /(128π³ β |s − mφ² + i mφ Γφ|²) (× (s−4mx²)/s
# for SP); Γ(φ→3π) = W(mφ²)/(64π³ mφ).  LO mixing on purpose (see
# threepi_weight.m); η pole generically, η′ for g ∝ [1 1 1].
function threepi_weight(s::Real, g)
    threepi_enabled() || return 0.0
    f  = FPI
    dm = M_DN - M_UP
    D  = g[1] - g[2]
    if _sppp_degenerate(g)
        mj = M_ETAP; Gj = G_ETAP
        aj = B0 * f * sqrt(2/3) * sum(g)
        cj = -sqrt(6) / 3 * B0 * dm / f^2       # 3! × the coefficient in L
        bj = -sqrt(6) / 9 * B0 * dm / f^2
    else
        mj = M_ETA; Gj = G_ETA
        aj = B0 * f * (g[1] + g[2] - 2g[3]) / sqrt(3)
        cj = -B0 * dm / (sqrt(3) * f^2)
        bj = -B0 * dm / (3 * sqrt(3) * f^2)
    end
    d2 = mj^2 - MPI0^2
    pole = aj / (s - mj^2 + im * mj * Gj)
    W = 0.0
    if s > 9MPI0^2                                   # π0π0π0, constant on the Dalitz plot
        A = B0 * D / f * s / (s - MPI0^2) + pole * cj
        W += abs2(A) * _dalitz_moments(s, MPI0, MPI0)[1] / 6
    end
    if s > (2MPI_PM + MPI0)^2                        # π+π−π0: A = P + Q X
        P = pole * bj * (1 - s / d2)
        Q = B0 * D / f / (s - MPI0^2) + pole * bj * 3 / d2
        I0, I1, I2 = _dalitz_moments(s, MPI0, MPI_PM)
        W += abs2(P) * I0 + 2real(P * conj(Q)) * I1 + abs2(Q) * I2
    end
    return W
end

# W(s) over π+π−φ (S = 1) and π0π0φ (S = 1/2), every propagator kept, so
# σ = gx² W/(128π³ β) (× (s−4mx²)/s for SP).  A = Σ_j U_φj V_j/(s − m_j²),
# O(g²), a function of m_ππ² only.  Production through detailed balance assumes
# φ is in the bath, while the :sm rate drops χχ̄→φφ and the coupled solver's Y_φ
# has no scattering source.  See docs/cross_sections.md, "Freeze-in consistency".
function pipiphi_weight(s::Real, mphi, g)
    pipiphi_enabled() || return 0.0
    f  = FPI
    b  = B0 * f
    p2 = MPI0^2
    D  = g[1] - g[2]
    Sg = g[1] + g[2]
    deg = _sppp_degenerate(g)
    if deg
        t0   = -b * sqrt(2/3) * sum(g) / (mphi^2 - M_ETAP^2)
        Rm   = 1 / (s - M_ETAP^2 + im * M_ETAP * G_ETAP)
        Nphi = -2t0 * (sqrt(6) * b * Sg - p2 * t0) / (3f^2)
        Nj   =   t0 * (sqrt(6) * b * Sg - 2p2 * t0) / (3f^2)
    else
        tp = -b * D / (mphi^2 - p2)
        te = -b * (g[1] + g[2] - 2g[3]) / sqrt(3) / (mphi^2 - M_ETA^2)
        Re = 1 / (s - M_ETA^2 + im * M_ETA * G_ETA)
    end
    W = 0.0
    if s > (2MPI_PM + mphi)^2                        # π+π−φ: A = P + Q X
        Rf = 1 / (s - mphi^2)
        if deg
            P = Nphi * Rf + Nj * Rm
            Q = 0.0
        else
            Rp = 1 / (s - p2)
            a  = mphi^2 + p2 + s - 3p2
            P  = Rf * (-2sqrt(3) * b * Sg * te - 2b * D * tp + p2 * te^2 - a * tp^2) / (3f^2) +
                 Rp * tp * (b * D + a * tp) / (3f^2) +
                 Re * te * (sqrt(3) * b * Sg - p2 * te) / (3f^2)
            Q  = tp^2 / f^2 * (Rf - Rp)
        end
        I0, I1, I2 = _dalitz_moments(s, mphi, MPI_PM)
        W += abs2(P) * I0 + 2real(P * conj(Q)) * I1 + abs2(Q) * I2
    end
    if s > (2MPI0 + mphi)^2                          # π0π0φ, constant on the Dalitz plot
        Rf = 1 / (s - mphi^2)
        if deg
            A = Nphi * Rf + Nj * Rm
        else
            A = -Rf * (2sqrt(3) * b * Sg * te + 6b * D * tp - p2 * te^2 - 3p2 * tp^2) / (3f^2) +
                1 / (s - p2) * tp * (b * D - p2 * tp) / f^2 +
                Re * te * (sqrt(3) * b * Sg - p2 * te) / (3f^2)
        end
        W += abs2(A) * _dalitz_moments(s, mphi, MPI0)[1] / 2
    end
    return W
end

# Γ(φ→3π); 0 below threshold (so mφ = 0 gives 0, not 0/0) or with LOOP_3PI off.
_phi_3pi_width(mphi, g) =
    (threepi_enabled() && mphi > 3MPI0) ? threepi_weight(mphi^2, g) / (64pi^3 * mphi) : 0.0

# χχ̄ → 3π + ππφ.  Gphi is the mediator's total width; isSP selects the scalar χ
# vertex's extra (s − 4mx²)/s.
function threebody_sigma(s::Real, mx, mphi, gx, Gphi, g, isSP::Bool)
    W3 = threepi_weight(s, g)
    sig = (W3 / ((s - mphi^2)^2 + mphi^2 * Gphi^2) + pipiphi_weight(s, mphi, g)) *
          gx^2 / (128pi^3) / sqrt(1 - 4mx^2 / s)
    return isSP ? sig * (s - 4mx^2) / s : sig
end

# Axial vector: Dalitz weights (WT, WL) of the transverse (J = 1, π+π−π0 only)
# and longitudinal (J = 0) parts of V^μ = ⟨3π|A^μ|0⟩, times Δa² with
# Δa = m_{π⁰A′}/F, the A′–π⁰ coupling the γγ channel uses.
#   σ = λ_a²/(32π³) [β WT/(6|D_A|²) + mx² WL/(β s² m_A′⁴)]
#   Γ(A′→3π) = WT(m_A′²)/(192π³ m_A′)
# q·V is the pseudoscalar π⁰ path × mπ²/(B0 Δg) (PCAC); the isospin-breaking
# η→3π path is omitted.  Port of axial3pi_weight.m.
function axial3pi_weight(s::Real, g)
    threepi_enabled() || return 0.0, 0.0
    f   = FPI
    Da2 = (_axial_mixings(g)[1] / f)^2
    cL  = Da2 * MPI0^4 / f^2
    WT = 0.0; WL = 0.0
    if s > 9MPI0^2                                   # π0π0π0: J = 0 only
        WL += cL * s^2 / (s - MPI0^2)^2 * _dalitz_moments(s, MPI0, MPI0)[1] / 6
    end
    if s > (2MPI_PM + MPI0)^2                        # π+π−π0; |p0|² quadratic in X
        I0, I1, I2 = _dalitz_moments(s, MPI0, MPI_PM)
        p02 = (s^2 * I0 - 2s * I1 + I2) / (4s) - MPI0^2 * I0
        WT += 4Da2 / f^2 * p02
        WL += cL * I2 / (s - MPI0^2)^2
    end
    return WT, WL
end

# Axial vector: Dalitz weights (WT, WL) of V_σ = ⟨π+π−A′|A_σ|0⟩, summed over
# the A′ polarizations:  V_σ = −2i Δa Δa0 [ε*_σ − q_σ (q·ε*)/(s − mπ0²)], the
# A′A′π+π− contact plus the π0 pole, Δa0 = a_u − a_d (LO).  L4 rescales a tree
# amplitude once, so Δa = m_{π⁰A′}/F appears once, as in the 3π current.
# Depends on E_A alone; no π0π0A′, η or η′ path.  Port of pipiA_weight.m.
function pipiA_weight(s::Real, mA, g)
    pipiphi_enabled() || return 0.0, 0.0
    s > (2MPI_PM + mA)^2 || return 0.0, 0.0
    c4 = 4 * (_axial_mixings(g)[1] / FPI)^2 * (g[1] - g[2])^2
    I0, I1, I2 = _dalitz_moments(s, mA, MPI_PM)
    c = s + mA^2 - MPI0^2                            # E_A = (c − X)/(2√s)
    K = (c^2 * I0 - 2c * I1 + I2) / (4s) - mA^2 * I0  # ∫ Y |k_A|²
    return c4 * (3I0 + K / mA^2), c4 * MPI0^4 * s * K / (mA^2 * (s - MPI0^2)^2)
end

# χχ̄ → 3π + π+π−A′ for the axial vector; GA is Γ(A′) total.  J = 1 and J = 0
# do not interfere.  Port of axial3pi_sigma.m.
function axial3pi_sigma(s::Real, mx, mA, gx, GA, g)
    WT, WL = axial3pi_weight(s, g) .+ pipiA_weight(s, mA, g)
    beta = sqrt(1 - 4mx^2 / s)
    return gx^2 / (32pi^3) * (beta * WT / (6 * ((s - mA^2)^2 + mA^2 * GA^2)) +
                              mx^2 * WL / (beta * s^2 * mA^4))
end

# Γ(A′→3π): only the transverse current decays an on-shell A′.
_axial_3pi_width(mA, g) =
    (threepi_enabled() && mA > 3MPI0) ? axial3pi_weight(mA^2, g)[1] / (192pi^3 * mA) : 0.0

# Thresholds of the enabled channels, as ⟨σv⟩ integration waypoints.
function _threebody_thresholds(mphi)
    t3 = threepi_enabled() ? (9MPI0^2, (2MPI_PM + MPI0)^2) : ()
    tf = pipiphi_enabled() ? ((2MPI0 + mphi)^2, (2MPI_PM + mphi)^2) : ()
    return (t3..., tf...)
end
