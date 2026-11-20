# Pole-coupling factories for the multi-pole models (SP/PP and AxialVector).
# Ports of sppp_couplings.m, axial_couplings.m, and phi_partial_widths.m's
# phi_gg_coupling helper.  Single source of truth for the meson-pole masses,
# couplings, and Breit-Wigner widths; the model structs precompute these once
# at construction and the cross sections / widths / poles read the stored arrays.

# Two-meson final states of a scalar φ coupled to quark masses.  Lowest-order
# χPT with M_q → M_q + gφ gives L ⊃ −B0 φ tr(G Φ²) (η = η₈, no η–η′ mixing):
# vertex amplitude μ and phase-space factor S (½ for identical pairs), so
# Γ(φ→AB) = S μ² √λ(mφ²,mA²,mB²)/(16π mφ³).  Port of meson_pairs.m.
function meson_pairs(g)
    gu, gd, gs = g[1], g[2], g[3]
    return ((mA = MPI_PM, mB = MPI_PM, mu = B0 * (gu + gd),             S = 1.0, family = :pipi),
            (mA = MPI0,   mB = MPI0,   mu = B0 * (gu + gd),             S = 0.5, family = :pipi),
            (mA = MK,     mB = MK,     mu = B0 * (gu + gs),             S = 1.0, family = :KK),
            (mA = MK0,    mB = MK0,    mu = B0 * (gd + gs),             S = 1.0, family = :KK),
            (mA = M_ETA,  mB = M_ETA,  mu = B0 * (gu + gd + 4gs) / 3,   S = 0.5, family = :eta),
            (mA = MPI0,   mB = M_ETA,  mu = B0 * (gu - gd) / sqrt(3),   S = 1.0, family = :eta))
end

# √λ(s, mA², mB²)/√s: twice the CM momentum (= √(s − 4m²) for equal masses)
_two_p(s, mA, mB) = sqrt(max(0.0, (s - (mA + mB)^2) * (s - (mA - mB)^2))) / sqrt(s)


# NLO pseudoscalar meson–mediator mass mixings m²_{π⁰φ̃}, m²_{ηφ̃} (port of
# sppp_couplings.m).  The η expression is 1/√3 times the π⁰ bracket, term by
# term: the chiral trace Tr[M_g Φ] = (g_u−g_d)π⁰/√2 + (g_u+g_d−2g_s)η₈/√6 gives
# the η₈ coefficient a relative √2/√6 = 1/√3 — the same algebra behind the
# √3F/3 in m_{ηA′} (_axial_couplings below).
function _sppp_mixings(g)
    summ      = M_UP + M_DN + M_S
    K         = B0 * FPI + 32 * B0^2 * L6 / FPI * summ
    sg        = g[1] + g[2] + g[3]
    m2_pi0phi = (g[1] - g[2]) * K + 32 * B0^2 * L7 / FPI * sg * (M_UP - M_DN)
    m2_etaphi = ((g[1] + g[2] - 2g[3]) * K +
                 32 * B0^2 * L7 / FPI * sg * (M_UP + M_DN - 2 * M_S)) / sqrt(3)
    # η′ (SU(3) singlet).  √6/3 = √(2/3) is the η₀ coefficient of Tr[g_q Π].
    # Replaced B0*FPI/sqrt(6) (LO only, and √6/3 → 1/√6 — together 4.3× low) on
    # 2026-08-12, when the draft appendix carrying that form was commented out.
    m2_etapphi = sqrt(6) / 3 * sg * (B0 * FPI + (L6 + L7) * 32 * B0^2 / FPI * summ)
    return m2_pi0phi, m2_etaphi, m2_etapphi
end

"""
    sppp_mixing_angles(mphi, g) -> (theta, theta_max, mj)

Small-mixing angles θ_j = m²_{jφ}/(mφ² − m_j²) for j = π⁰, η, η′, and
`theta_max = max|θ_j|`.  The 2026-08-12 draft states the SP/PP small-mixing
limit as |m²_{jφ}| ≪ |mφ² − m_j²| for all three poles at once, so `theta_max`
is the expansion parameter the whole sector rests on; `model_observables`
rejects points where it exceeds `s_mix`.  The η′ is included even though the
generic pole set is {π⁰, φ′, η}: the condition is on the validity of the
diagonalization, not on which poles are kept.  Mirrors sppp_mixing.m.
"""
function sppp_mixing_angles(mphi, g)
    m2 = _sppp_mixings(g)
    mj = (MPI0, M_ETA, M_ETAP)
    theta = ntuple(i -> m2[i] / (mphi^2 - mj[i]^2), 3)
    return theta, maximum(abs, theta), mj
end

# Effective φ′→γγ coupling for SP/PP (port of phi_partial_widths.m's
# phi_gg_coupling).  Generic: φ′ mixes with π⁰/η (∝ g_u−g_d, g_u+g_d−2g_s).  For
# g_u=g_d=g_s both vanish and it is induced via the η′ singlet (∝ g_u+g_d+g_s).
function _phi_gg_coupling(mphi, g)
    if g[1] == g[2] && g[2] == g[3]
        F_etap = sqrt(ALPHA_EM^2 * M_ETAP^3 / (64pi^3 * G_ETAP_GG))
        m_mix  = _sppp_mixings(g)[3]
        # ALPHA_EM/(pi*F_etap): F_etap is a decay constant, so it takes the same
        # 1/pi as the gpigamma convention below.  Fixed 2026-08-05.
        return (ALPHA_EM / (pi * F_etap)) * m_mix / (mphi^2 - M_ETAP^2)
    else
        gpigamma  = ALPHA_EM / pi / FPI
        getagamma = ALPHA_EM / sqrt(3) / pi / FPI
        m2_pi0phi, m2_etaphi, _ = _sppp_mixings(g)
        # Both terms same sign: g_{φ′γγ} = -Σ_j g_{jγγ} θ_j.  The π⁰ term used to
        # carry the opposite sign; fixed 2026-08-10.
        return -gpigamma  * m2_pi0phi / (mphi^2 - MPI0^2) -
                getagamma * m2_etaphi / (mphi^2 - M_ETA^2)
    end
end

# Mediator (gg, xx, total) widths for the SP/PP family.  isPP selects the χχ̄
# vertex velocity factor: β¹ (PP, pseudoscalar χ vertex) vs β³ (SP, scalar χ
# vertex).  φ→ππ is parity forbidden; `pipi` carries Γ(φ→3π), which is 0
# unless LOOP_3PI=1 (models/threebody.jl).
function _sppp_widths(mphi, mx, gx, g, isPP::Bool)
    gphigg = _phi_gg_coupling(mphi, g)
    gg = gphigg^2 * mphi^3 / (64pi)
    if isPP
        xx = gx^2 / (4pi) * sqrt(max(0.0, mphi^2 / 4 - mx^2))
    else
        xx = mphi > 2mx ? gx^2 * max(0.0, mphi^2 - 4mx^2)^1.5 / (8pi * mphi^2) : 0.0
    end
    pipi = _phi_3pi_width(mphi, g)
    return (gg = gg, xx = xx, pipi = pipi, total = gg + xx + pipi)
end

"""
    _sppp_couplings(mphi, mx, gx, g, isPP) -> (mj, gcx, ggx, Gam)

s-channel pole couplings/widths for SP/PP: the generic {π⁰, φ′, η} set, or the
SU(3)-singlet {η′, φ′} set when `g[1]==g[2]==g[3]` exactly (where the π⁰/η mixing
and the induced φ′γγ coupling all vanish).  Port of sppp_couplings.m.
"""
function _sppp_couplings(mphi, mx, gx, g, isPP::Bool)
    Gphi = _sppp_widths(mphi, mx, gx, g, isPP).total

    if g[1] == g[2] && g[2] == g[3]
        # degenerate g: η′ singlet + mediator φ′
        F_etap      = sqrt(ALPHA_EM^2 * M_ETAP^3 / (64pi^3 * G_ETAP_GG))
        g_etapgamma = ALPHA_EM / (pi * F_etap)   # 1/pi: as gpigamma, see sppp_couplings.m
        m_mix       = _sppp_mixings(g)[3]
        denom       = (mphi^2 - M_ETAP^2)

        g_etapx        = -gx * m_mix / denom
        gphiprimex     = gx
        gphiprimegamma = g_etapgamma * m_mix / denom      # induced φ′→γγ

        Gam_etapxx = isPP ?
            g_etapx^2 * sqrt(max(0.0, M_ETAP^2 / 4 - mx^2)) / (4pi) :      # β¹
            g_etapx^2 * max(0.0, M_ETAP^2 - 4mx^2)^1.5 / (8pi * M_ETAP^2)  # β³
        Gam_etap = G_ETAP + Gam_etapxx

        mj  = [M_ETAP, mphi]
        gcx = [g_etapx, gphiprimex]
        ggx = [g_etapgamma, gphiprimegamma]
        Gam = [Gam_etap, Gphi]
    else
        # generic g: π⁰, mediator φ′, η
        gpigamma       = ALPHA_EM / pi / FPI
        getagamma      = ALPHA_EM / sqrt(3) / pi / FPI
        m2_pi0phi, m2_etaphi, _ = _sppp_mixings(g)
        # θ_j = m²_{jφ}/(mφ²-m_j²) for BOTH poles, so g_{jχ} = +gx θ_j with the
        # same sign for π⁰ and η, and g_{φ′γγ} = -Σ_j g_{jγγ} θ_j.  Until
        # 2026-08-10 the π⁰ block carried the opposite sign, so π⁰ and η
        # interfered the wrong way round.
        gpix           = gx * m2_pi0phi / (mphi^2 - MPI0^2)
        gphiprimex     = gx
        gphiprimegamma = -gpigamma  * m2_pi0phi / (mphi^2 - MPI0^2) -
                          getagamma * m2_etaphi / (mphi^2 - M_ETA^2)
        getax          = gx * m2_etaphi / (mphi^2 - M_ETA^2)

        mj  = [MPI0, mphi, M_ETA]
        gcx = [gpix, gphiprimex, getax]
        ggx = [gpigamma, gphiprimegamma, getagamma]
        Gam = [GPI, Gphi, G_ETA]
    end
    return mj, gcx, ggx, Gam
end

# AxialVector validity-floor error (port of axial_couplings:mediatorTooLight):
# λ_jχ ∝ 1/m_A′² blows the meson→χχ̄ width past its mass at small m_A′, where the
# mixing-EFT / Breit-Wigner picture is meaningless.  Thrown at construction; the
# M6 scan catches it and marks the point NaN.
struct AxialMediatorTooLight <: Exception
    msg::String
end
Base.showerror(io::IO, e::AxialMediatorTooLight) = print(io, "axial_couplings:mediatorTooLight: ", e.msg)

# A′–meson mixing masses (m_{π⁰A′}, m_{ηA′}, m_{η′A′}) [GeV].  Shared by
# _axial_couplings and the axial 3π current (threebody.jl).  Uses the corr4 L4 `+` sign (see
# docs/axialvector.md).  Port of axial_mixing.m.
function _axial_mixings(g)
    summ  = M_UP + M_DN + M_S
    corr4 = 1 + 16 * B0 * L4 / FPI^2 * summ            # "+" branch (sign ambiguity)
    au, ad, as = g[1], g[2], g[3]
    m_pi0Ap  = FPI            * ((au - ad)        * corr4 + 16 * B0 * L5 * (au * M_UP - ad * M_DN)              / FPI^2)
    m_etaAp  = sqrt(3) * FPI / 3 * ((au + ad - 2as)  * corr4 + 16 * B0 * L5 * (au * M_UP + ad * M_DN - 2as * M_S) / FPI^2)
    m_etapAp = sqrt(6) * FPI / 3 * ((au + ad + as)   * corr4 + 16 * B0 * L5 * (au * M_UP + ad * M_DN + as * M_S)  / FPI^2)
    return m_pi0Ap, m_etaAp, m_etapAp
end

"""
    _axial_couplings(mphi, mx, gx, g) -> (mj, lam_chi, lam_gg, Gam_j)

A′–pseudoscalar (π⁰, η, η′) mixing couplings and BW widths for AxialVector.
`mphi=m_A′`, `gx=λ_a`, `g=[a_u a_d a_s]`.  λ_jχ = λ_a·m_jA′/m_A′².  Meson BW
widths = SM total + Γ(j→χχ̄).  Throws `AxialMediatorTooLight` once any meson
width reaches its mass (the EFT validity floor).  Uses the corr4 L4 `+` sign
(unresolved ambiguity — see docs/axialvector.md).  Port of axial_couplings.m.
"""
function _axial_couplings(mphi, mx, gx, g)
    F_etap     = sqrt(ALPHA_EM^2 * M_ETAP^3 / (64pi^3 * G_ETAP_GG))
    # Data-driven F_eta as of the 2026-08-11 draft Eq. (lambdagamma): the LO
    # octet value ALPHA_EM/(sqrt(3) pi FPI) gives Gamma(eta->gg) = 0.17 keV
    # against a measured 0.52 keV.  1.73x in amplitude.  The draft's
    # pseudoscalar section still carries the LO form, so _sppp_couplings below
    # is deliberately NOT changed to match.
    F_eta      = sqrt(ALPHA_EM^2 * M_ETA^3 / (64pi^3 * G_ETA_GG))
    lam_pigg   = ALPHA_EM / (pi * FPI)
    lam_etagg  = ALPHA_EM / (pi * F_eta)
    lam_etapgg = ALPHA_EM / (pi * F_etap)   # 1/pi: all three share the Gamma=lam^2 m^3/(64pi) convention

    m_pi0Ap, m_etaAp, m_etapAp = _axial_mixings(g)

    mj      = [MPI0, M_ETA, M_ETAP]
    lam_gg  = [lam_pigg, lam_etagg, lam_etapgg]
    lam_chi = (gx / mphi^2) .* [m_pi0Ap, m_etaAp, m_etapAp]

    Gam_SM = [GPI, G_ETA, G_ETAP]
    Gam_xx = lam_chi .^ 2 .* mx^2 .* sqrt.(max.(0.0, mj .^ 2 ./ 4 .- mx^2)) ./ pi
    Gam_j  = Gam_SM .+ Gam_xx

    rmax, jmax = findmax(Gam_j ./ mj)
    if rmax >= 1
        jn = ("pi0", "eta", "etaprime")
        throw(AxialMediatorTooLight(
            "$(jn[jmax]) BW width exceeds its mass (Gamma/m = $(round(rmax, sigdigits=3))) " *
            "at m_A'=$mphi GeV, mx=$mx, lambda_a=$gx; mixing-EFT/Breit-Wigner invalid here."))
    end
    return mj, lam_chi, lam_gg, Gam_j
end

# Probe `_axial_couplings` at a λ_a small enough to clear the mediatorTooLight
# floor (λ_jχ is linear in λ_a and the mixing masses are λ_a-independent, so any
# sub-floor probe recovers them exactly).  Returns the probe value + the arrays,
# or `nothing` if the floor is unreachable even at vanishing λ_a.
function _axial_probe(mphi, mx, g)
    gx_probe = 1.0
    for _ in 1:200
        try
            mj, lam_chi, lam_gg, Gam_j = _axial_couplings(mphi, mx, gx_probe, g)
            return (gx_probe, mj, lam_chi, lam_gg, Gam_j)
        catch e
            e isa AxialMediatorTooLight || rethrow(e)
            gx_probe *= 1e-6                       # step well below the floor
        end
    end
    return nothing
end

"""
    axial_gx_ceiling(mphi, mx, g) -> gx_max

Largest axialVector coupling gx (=λ_a) for which the meson-mixing EFT is
self-consistent: solves max_j Γ_j/m_j = 1 for gx (the companion of
`_axial_couplings`'s mediatorTooLight throw).  `Inf` if all χχ̄ channels are
closed (mⱼ ≤ 2mx), `NaN` if the floor is unreachable.  Port of axial_gx_ceiling.m.
"""
function axial_gx_ceiling(mphi, mx, g)
    pr = _axial_probe(mphi, mx, g)
    pr === nothing && return NaN
    gx_probe, mj, lam_chi, _, Gam_j = pr
    m_jAp  = abs.(lam_chi) .* mphi^2 ./ gx_probe          # |λ_jχ| = gx |m_jA′|/mφ²
    kin    = mx^2 .* sqrt.(max.(0.0, mj .^ 2 ./ 4 .- mx^2)) ./ pi
    Gam_xx = lam_chi .^ 2 .* kin
    Gam_SM = Gam_j .- Gam_xx
    headroom = max.(0.0, mj .- Gam_SM)                    # mⱼ − Γ_SM,j (~mⱼ)
    gx_max_j = fill(Inf, length(mj))
    @inbounds for j in eachindex(mj)
        kin[j] > 0 && (gx_max_j[j] = (mphi^2 / m_jAp[j]) * sqrt(headroom[j] / kin[j]))
    end
    return minimum(gx_max_j)
end

"""
    axial_mixing_floor(mphi, mx, g, smax=0.3) -> (eta_mix, mxmix, jbind)

Small-mixing-breakdown boundary: η at which the largest A′–meson mixing angle
sin θ_j = η|m_jA′|/m_A′ reaches `smax`.  `mxmix = max_j|m_jA′|` (so a gx-scan can
form sin θ = mxmix/mφ).  Overlay diagnostic only — never gates a scan.  Port of
axial_mixing_floor.m.
"""
function axial_mixing_floor(mphi, mx, g, smax=0.3)
    pr = _axial_probe(mphi, mx, g)
    pr === nothing && return (NaN, NaN, NaN)
    gx_probe, _, lam_chi, _, _ = pr
    m_jAp = abs.(lam_chi) .* mphi^2 ./ gx_probe
    mxmix, jbind = findmax(m_jAp)
    eta_mix = mxmix > 0 ? smax * mphi / mxmix : Inf
    return (eta_mix, mxmix, jbind)
end

"""
    axial_eta_ceiling(mphi, mx, gx, g) -> eta_max

EFT-validity ceiling on the η-scan normalization (companion of
`axial_gx_ceiling` for the η-scan): largest η with `_axial_couplings(mphi, mx,
gx, η·g)` valid, bisected on log10 η.  `Inf` if η=1e12 is still valid, `NaN` if
even η=1e-30 is invalid.  Port of axial_eta_ceiling.m.
"""
function axial_eta_ceiling(mphi, mx, gx, g)
    valid(le) = _axial_eta_valid(mphi, mx, gx, 10.0^le, g)
    lo = -30.0; hi = 12.0
    valid(hi) && return Inf
    valid(lo) || return NaN
    for _ in 1:60
        mid = 0.5 * (lo + hi)
        valid(mid) ? (lo = mid) : (hi = mid)
    end
    return 10.0^lo
end

function _axial_eta_valid(mphi, mx, gx, eta, g)
    try
        _axial_couplings(mphi, mx, gx, eta .* g)
        return true
    catch e
        e isa AxialMediatorTooLight && return false
        rethrow(e)
    end
end
