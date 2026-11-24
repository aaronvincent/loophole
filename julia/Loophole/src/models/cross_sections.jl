# Cross sections for the single-pole scalar models.
#
# Ports of the heavyScalar/lightScalar and …PS branches of ggxxxsec.m (γγ→χχ),
# make_xxsigma.m ('gg' piece, χχ→γγ), and xxggxsec.m (χχ→all).  Each takes a
# scalar Mandelstam s [GeV²]; broadcast for a grid.  Heavy/light share formulas.
#
# Detailed-balance / S_γγ convention (see docs/cross_sections.md):
#   every γγ→χχ production cross section carries a factor 2 relative to its
#   χχ→γγ partner, i.e. sigma_ggxx/sigma_xxgg = 2(s−4mχ²)/s.  The 2 is
#   identical-photon statistics S_γγ = 1/2, which applies in the annihilation
#   direction only.  Concretely: /4608 vs /9216 π³ (SS, PS), /64π vs /128π
#   (SP, PP), /8π vs /16π (AxialVector).

_prop2(m::CrossSectionModel, s) = (s - m.mphi^2)^2 + m.mphi^2 * mediator_widths(m).total^2

# =====================================================================
# Scalar-scalar (heavyScalar / lightScalar)
# =====================================================================

# γγ → χχ  (freeze-in production)
function sigma_ggxx(m::ScalarSS, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    D = dscript(m.g)
    return m.gx^2 * ALPHA_EM^2 * D / (4608 * pi^3) *
           s * (s - 4 * m.mx^2) / _prop2(m, s) * sqrt((s - 4 * m.mx^2) / s)
end

# χχ → γγ  (s-channel only; production partner / 'gg' channel)
function sigma_xxgg(m::ScalarSS, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    D = dscript(m.g)
    # /9216 (not /4608): the identical-photon S_gg = 1/2 lives here.  This is
    # now exactly pi s Gamma_xx Gamma_gg / (p_chi^2 D).
    return m.gx^2 * ALPHA_EM^2 * D / (9216 * pi^3) *
           s * (s - 4 * m.mx^2) / _prop2(m, s) * sqrt(s / max(s - 4 * m.mx^2, eps()))
end

# χχ → all  (depletion): γγ + φφ (t-channel) + ππ + KK
function sigma_xxall(m::ScalarSS, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    mx, mphi, gx, g = m.mx, m.mphi, m.gx, m.g
    D = dscript(g)
    p2 = _prop2(m, s)
    sig = 0.0
    # φφ (t-channel), open for s > 4mφ².
    # s > 4mφ² is the only gate: the physical threshold, with s spanning
    # [4mx², …) in a thermal bath, so the channel is reachable for mx < mφ at a
    # cost of exp(-2(mφ-mx)/T).  thermal_average.jl carries 4mφ² as a waypoint
    # so the adaptive rule resolves the cusp the range then crosses.
    if s > 4 * mphi^2
        sig += _sig_pp_scalar(s, mx, mphi, gx)
    end
    # γγ (s-channel BW)
    sig += gx^2 * ALPHA_EM^2 * D / (9216 * pi^3) *
           s * (s - 4 * mx^2) / p2 * sqrt(s / (s - 4 * mx^2))
    return sig + sigma_xxhad(m, s)
end

# χχ → meson pairs (meson_pairs): S μ² times the π⁺π⁻ form, with
# √(s − 4m²) → √λ(s, mA², mB²)/√s
function sigma_xxhad(m::ScalarSS, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    mx, gx = m.mx, m.gx
    p2 = _prop2(m, s)
    sig = 0.0
    for P in meson_pairs(m.g)
        s > (P.mA + P.mB)^2 || continue
        sig += P.S * P.mu^2 * gx^2 / (32pi) * _two_p(s, P.mA, P.mB) * sqrt(s - 4 * mx^2) / (s * p2)
    end
    return sig
end

# χχ̄ → φφ isolated piece (dark-sector rate via svavg_phiphi).  Opens at s > 4mφ² — matches make_xxsigma(...,'phiphi').
# Callers needing mx ≥ mφ impose it themselves.
function sigma_xxphiphi(m::ScalarSS, s::Real)
    s > 4 * m.mphi^2 || return 0.0
    return _sig_pp_scalar(s, m.mx, m.mphi, m.gx)
end

# χχ̄ → φφ, scalar φ-χ-χ vertex (closed form, xxggxsec.m / make_xxsigma.m).
# Above threshold (4mφ²−s) < 0, so sqrt/atan go complex; the physical result is
# real(...) — matches MATLAB, which let the algebra go complex and took real().
function _sig_pp_scalar(s::Real, mx, mphi, gx)
    kin   = sqrt((s - 4mx^2) * (s - 4mphi^2))
    A     = -32mx^4 + 16mx^2 * (s - mphi^2) + 6mphi^4 - 4mphi^2 * s + s^2
    atarg = sqrt(complex((s - 4mx^2) * (4mphi^2 - s))) / (2mphi^2 - s)
    t1    = 2 * (mx^2 * (s - 4mphi^2) + mphi^4) * A * atan(atarg) /
            ((s - 4mx^2)^(3/2) * (2mphi^2 - s) * sqrt(complex(4mphi^2 - s)))
    t2    = (16mx^4 + 2mx^2 * (s - 8mphi^2) + 3mphi^4) / (4mx^2 - s)
    return real(gx^4 * kin * (t1 + t2) / (32pi * s * (mx^2 * (s - 4mphi^2) + mphi^4)))
end

# =====================================================================
# Scalar-pseudoscalar (heavyScalarPS / lightScalarPS)
# =====================================================================

function sigma_ggxx(m::ScalarPS, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    D = dscript(m.g)
    # /4608 (not /9216): identical-photon factor 2 on the production side, same
    # convention as ScalarSS — see the header note.
    return m.gx^2 * ALPHA_EM^2 * D / (4608 * pi^3) *
           s^2 / _prop2(m, s) * sqrt((s - 4 * m.mx^2) / s)
end

function sigma_xxgg(m::ScalarPS, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    D = dscript(m.g)
    # /9216 (not /4608): identical-photon S_gg = 1/2, as in ScalarSS above.
    return m.gx^2 * ALPHA_EM^2 * D / (9216 * pi^3) *
           s^2 / _prop2(m, s) * sqrt(s / max(s - 4 * m.mx^2, eps()))
end

function sigma_xxall(m::ScalarPS, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    mx, mphi, gx, g = m.mx, m.mphi, m.gx, m.g
    D = dscript(g)
    p2 = _prop2(m, s)
    sig = 0.0
    if s > 4 * mphi^2
        sig += _sig_pp_ps(s, mx, mphi, gx)
    end
    sig += gx^2 * ALPHA_EM^2 * D / (9216 * pi^3) * s^2 / p2 * sqrt(s / (s - 4mx^2))
    return sig + sigma_xxhad(m, s)
end

function sigma_xxhad(m::ScalarPS, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    mx, gx = m.mx, m.gx
    p2 = _prop2(m, s)
    sig = 0.0
    for P in meson_pairs(m.g)
        s > (P.mA + P.mB)^2 || continue
        sig += P.S * P.mu^2 * gx^2 / (32pi) * _two_p(s, P.mA, P.mB) / (sqrt(s - 4mx^2) * p2)
    end
    return sig
end

function sigma_xxphiphi(m::ScalarPS, s::Real)
    s > 4 * m.mphi^2 || return 0.0
    return _sig_pp_ps(s, m.mx, m.mphi, m.gx)
end

# χχ̄ → φφ (or → jj), pseudoscalar χ vertex.  atanh argument can leave (−1,1),
# so do it in complex and take real().
#
# REPLACED 2026-08-11 together with xxggxsec.m.  The previous closed form was
# negative and ~1/v_rel in the non-relativistic regime; the true result is
# positive and ~v_rel (a pseudoscalar χ vertex has no s-wave at all).  This form
# matches a direct tree-level trace calculation to all digits and already
# carries the identical-particle S = 1/2.
function _sig_pp_ps(s::Real, mx, mphi, gx)
    kin = sqrt((s - 4mx^2) * (s - 4mphi^2))
    t1  = -sqrt((s - 4mphi^2) / (s - 4mx^2)) *
          (2mx^2 * (s - 4mphi^2) + 3mphi^4) /
          (2 * (mx^2 * (s - 4mphi^2) + mphi^4))
    t2  = (6mphi^4 - 4mphi^2 * s + s^2) / ((s - 4mx^2) * (s - 2mphi^2)) *
          atanh(complex(kin / (s - 2mphi^2)))
    return real(gx^4 * (t1 + t2) / (16pi * s))
end

# =====================================================================
# SP/PP family — pseudoscalar mediator via meson poles (sppp_couplings).
# γγ↔χχ runs through the {π⁰,φ′,η} (or singlet {η′,φ′}) s-channel poles, with
# complex BW propagators summed in `_sppp_sumprop`.  Ports of the SP/PP branches
# of ggxxxsec.m, xxggxsec.m, make_xxsigma.m.
# =====================================================================

# Σ_j g_cχ·g_γγ / (m_j² − s − i m_j Γ_j)
function _sppp_sumprop(m::PoleSPPP, s::Real)
    sp = zero(ComplexF64)
    @inbounds for j in eachindex(m.mj)
        sp += m.gcx[j] * m.ggx[j] / (m.mj[j]^2 - s - im * m.mj[j] * m.Gam[j])
    end
    return sp
end

# χχ̄ → 3π + ππφ (models/threebody.jl), part of both 'had' and 'all'; exactly
# 0.0 with LOOP_3PI and LOOP_PIPIPHI off.  m.Gam[2] is the mediator's width in
# both pole orderings.
function _sppp_threebody(m::PoleSPPP, s::Real)
    (threepi_enabled() || pipiphi_enabled()) || return 0.0
    return threebody_sigma(s, m.mx, m.mphi, m.gx, m.Gam[2], m.g, !_is_pp(m))
end

# --- SP (scalar χ vertex, β³) ---
function sigma_ggxx(m::PoleSP, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    # /64pi (not /128pi): identical-photon factor 2 on the production side.
    return sqrt(s) * (s - 4m.mx^2)^1.5 / (64pi) * abs2(_sppp_sumprop(m, s))
end
function sigma_xxgg(m::PoleSP, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    return s^1.5 * sqrt(s - 4m.mx^2) / (128pi) * abs2(_sppp_sumprop(m, s))
end
function sigma_xxall(m::PoleSP, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    sig = sigma_xxgg(m, s)                                   # γγ piece
    @inbounds for j in eachindex(m.mj)
        if s > 4 * m.mj[j]^2
            sig += _sig_pp_scalar(s, m.mx, m.mj[j], m.gcx[j])  # χχ̄→jj (t-channel)
        end
    end
    return sig + _sppp_threebody(m, s)
end
# meson pairs only: j = 2 is φ' (couplings.jl ordering), a dark final state
function sigma_xxhad(m::PoleSP, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    sig = 0.0
    @inbounds for j in eachindex(m.mj)
        if j != 2 && s > 4 * m.mj[j]^2
            sig += _sig_pp_scalar(s, m.mx, m.mj[j], m.gcx[j])
        end
    end
    return sig + _sppp_threebody(m, s)
end

# --- PP (pseudoscalar χ vertex, β¹) ---
function sigma_ggxx(m::PolePP, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    # /64pi (not /128pi): identical-photon factor 2 on the production side.
    return s^2 / (64pi) * abs2(_sppp_sumprop(m, s)) * sqrt((s - 4m.mx^2) / s)
end
function sigma_xxgg(m::PolePP, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    return s^2 / (128pi) * abs2(_sppp_sumprop(m, s)) * sqrt(s / (s - 4m.mx^2))
end
function sigma_xxall(m::PolePP, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    sig = sigma_xxgg(m, s)                                   # γγ piece
    @inbounds for j in eachindex(m.mj)
        if s > 4 * m.mj[j]^2
            sig += _sig_pp_ps(s, m.mx, m.mj[j], m.gcx[j])      # χχ̄→jj (t-channel)
        end
    end
    return sig + _sppp_threebody(m, s)
end
function sigma_xxhad(m::PolePP, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    sig = 0.0
    @inbounds for j in eachindex(m.mj)
        if j != 2 && s > 4 * m.mj[j]^2
            sig += _sig_pp_ps(s, m.mx, m.mj[j], m.gcx[j])
        end
    end
    return sig + _sppp_threebody(m, s)
end

"""
    sigma_xxdark(m, s) -> σ(χχ̄ → dark mediator pair) [GeV⁻²]

`sigma_xxall − (sigma_xxgg + sigma_xxhad)`, built directly: φφ (S⊗S, P⊗S), A′A′
(axialVector), the φ′φ′ term (j = 2) of the SP/PP jj sum.  All t-channel.  Feeds
freeze_coupled's χχ̄ ↔ φφ rates (`sv_darkpair_grid`).  Port of xxggxsec(...,'dark').
"""
sigma_xxdark(m::CrossSectionModel, s::Real) = sigma_xxphiphi(m, s)
function sigma_xxdark(m::PoleSP, s::Real)
    (s <= 4 * m.mx^2 || s <= 4 * m.mj[2]^2) && return 0.0
    return _sig_pp_scalar(s, m.mx, m.mj[2], m.gcx[2])
end
function sigma_xxdark(m::PolePP, s::Real)
    (s <= 4 * m.mx^2 || s <= 4 * m.mj[2]^2) && return 0.0
    return _sig_pp_ps(s, m.mx, m.mj[2], m.gcx[2])
end

# =====================================================================
# AxialVector — spin-1 A' via {π⁰,η,η′} SM-meson mixing poles.
# γγ↔χχ runs through the meson poles (NOT an A' pole); depletion adds A'A'
# (t/u-channel χ exchange) and meson pairs.  Ports of the axialVector branches.
# =====================================================================

# Σ_j λ_jχ·λ_jγγ / (s − m_j² + i m_j Γ_j)  (axial sign convention; |·|² is the
# same as the SP/PP convention but kept verbatim to match MATLAB).
function _axial_sumprop(m::AxialVector, s::Real)
    sp = zero(ComplexF64)
    @inbounds for j in 1:3
        sp += m.lam_chi[j] * m.lam_gg[j] / (s - m.mj[j]^2 + im * m.mj[j] * m.Gam_j[j])
    end
    return sp
end

function sigma_ggxx(m::AxialVector, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    return m.mx^2 * s^2 / (16pi) * abs2(_axial_sumprop(m, s)) * sqrt((s - 4m.mx^2) / s)
end
function sigma_xxgg(m::AxialVector, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    return m.mx^2 * s^2 / (32pi) * abs2(_axial_sumprop(m, s)) * sqrt(s / (s - 4m.mx^2))
end
function sigma_xxall(m::AxialVector, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    mx, mphi, gx = m.mx, m.mphi, m.gx
    sig = sigma_xxgg(m, s)                                   # γγ piece
    # A'A' (t/u-channel χ exchange)
    if s > 4 * mphi^2
        sig += _sig_axial_pp(s, mx, mphi, gx)
    end
    return sig + sigma_xxhad(m, s)                           # meson pairs jj
end
function sigma_xxhad(m::AxialVector, s::Real)
    s <= 4 * m.mx^2 && return 0.0
    sig = 0.0
    @inbounds for j in 1:3
        if s > 4 * m.mj[j]^2
            sig += _sig_axial_mesonpair(s, m.mx, m.mj[j], m.lam_chi[j])
        end
    end
    # χχ̄ → 3π, π+π−A′ (threebody.jl); exactly 0.0 with both toggles off
    if threepi_enabled() || pipiphi_enabled()
        sig += axial3pi_sigma(s, m.mx, m.mphi, m.gx, mediator_widths(m).total, m.g)
    end
    return sig
end

# χχ̄ → A'A' isolated piece (dark-sector rate via svavg_phiphi).  Same A'A' closed form as sigma_xxall, opening
# at s > 4mφ² like every other pair channel.  Depends
# only on gx (=λ_a) and the masses, NOT on the meson mixing, so it never trips
# the mediatorTooLight floor.  Port of make_xxsigma(...,'phiphi') (sig_pp_axial).
function sigma_xxphiphi(m::AxialVector, s::Real)
    s > 4 * m.mphi^2 || return 0.0
    return _sig_axial_pp(s, m.mx, m.mphi, m.gx)
end

# χχ̄ → A'A' (t/u-channel χ exchange).  Above threshold the inner sqrt goes
# imaginary (atan→complex); MATLAB let it go complex and took real().
function _sig_axial_pp(s::Real, mx, mphi, gx)
    sqneg = sqrt(complex((4mphi^2 - s) * (s - 4mx^2)))      # imaginary for s>4mφ²
    atrm  = atan(sqneg / (2mphi^2 - s))
    base  = mphi^4 - 4mphi^2 * mx^2 + mx^2 * s
    poly  = 16mphi^8 - 16mphi^6 * (4mx^2 + s) + 2mphi^4 * s * (8mx^2 + 5s) -
            4mphi^2 * s * (-16mx^4 - 2mx^2 * s + s^2) - 16mx^4 * s^2 + s^4
    num   = -8mphi^8 + 8mphi^6 * (2mx^2 + s) +
            mphi^4 * (64mx^4 - 28mx^2 * s - 3s^2) +
            4mphi^2 * mx^2 * s * (3s - 16mx^2) +
            2base * poly * atrm / ((2mphi^2 - s) * sqneg) -
            2mx^2 * s^2 * (s - 4mx^2)
    return real(gx^4 * sqrt((s - 4mphi^2) / (s - 4mx^2)) * num / (32pi * mphi^4 * s * base))
end

# χχ̄ → jj (j = π⁰/η/η′), axial-vector meson-pair closed form.  atanh arg can
# leave (−1,1); done complex, real() taken (matches MATLAB).
function _sig_axial_mesonpair(s::Real, mx, mj, lam)
    kin = sqrt((s - 4mj^2) * (s - 4mx^2))
    den = mj^4 + mx^2 * (s - 4mj^2)
    t1  = kin * (mx^2 * (4mj^4 + 8mj^2 * s - 2s^2) - 3mj^4 * s) / den
    t2  = 2 * (s * (6mj^4 - 4mj^2 * s + s^2) + 4mx^2 * (2mj^4 - 4mj^2 * s + s^2)) *
          atanh(complex(kin / (2mj^2 - s))) / (2mj^2 - s)
    return real(lam^4 * (t1 + t2) / (32pi * s))
end
