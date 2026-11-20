# Mediator partial widths — port of phi_partial_widths.m.
#
# Returns a NamedTuple (gg, xx, pipi, KK, eta, nn, ll, ee, total) [GeV]: pipi = π⁺π⁻ + π⁰π⁰,
# KK = K⁺K⁻ + K⁰K̄⁰, eta = ηη + π⁰η.  The total width regulates the Breit-Wigner
# s-channel propagator in the cross sections.

# Two-meson widths of the scalar φ–SM vertex (S⊗S, P⊗S), from meson_pairs:
# Γ = S μ² √λ(mφ², mA², mB²)/(16π mφ³), gated on threshold (mφ → 0 stays finite).
function _meson_widths(mphi, g)
    w = Dict(:pipi => 0.0, :KK => 0.0, :eta => 0.0)
    for P in meson_pairs(g)
        mphi > P.mA + P.mB || continue
        lam = (mphi^2 - (P.mA + P.mB)^2) * (mphi^2 - (P.mA - P.mB)^2)
        w[P.family] += P.S * P.mu^2 * sqrt(lam) / (16pi * mphi^3)
    end
    return w[:pipi], w[:KK], w[:eta]
end

# Scalar-scalar (heavyScalar / lightScalar): γγ + χχ̄ (β³) + ππ + KK.
function mediator_widths(m::ScalarSS)
    mphi, mx, gx, g = m.mphi, m.mx, m.gx, m.g
    D = dscript(g)
    gg = ALPHA_EM^2 * D * mphi^3 / (4608 * pi^3)
    # β³ scalar-vertex χχ̄ width, same form as the SP branch (identical φχχ̄
    # vertex); 1/mφ² gated on threshold (closed channel → 0, not Inf·0 = NaN).
    xx = mphi > 2mx ? gx^2 * (mphi^2 - 4mx^2)^1.5 / (8pi * mphi^2) : 0.0
    pipi, KK, eta = _meson_widths(mphi, g)
    return (gg = gg, xx = xx, pipi = pipi, KK = KK, eta = eta, total = gg + xx + pipi + KK + eta)
end

# Scalar-pseudoscalar coupling (…PS): γγ (/4608) + χχ̄ + meson pairs.
function mediator_widths(m::ScalarPS)
    mphi, mx, gx, g = m.mphi, m.mx, m.gx, m.g
    D = dscript(g)
    gg = ALPHA_EM^2 * D * mphi^3 / (4608 * pi^3)
    xx = gx^2 / (4pi) * sqrt(max(0.0, mphi^2 / 4 - mx^2))
    pipi, KK, eta = _meson_widths(mphi, g)
    return (gg = gg, xx = xx, pipi = pipi, KK = KK, eta = eta, total = gg + xx + pipi + KK + eta)
end

# SP/PP family: γγ via π⁰/η (or η′) mixing + χχ̄ (β³ for SP, β¹ for PP).  Shares _sppp_widths with the coupling factory.
function mediator_widths(m::PoleSPPP)
    w = _sppp_widths(m.mphi, m.mx, m.gx, m.g, _is_pp(m))
    return (gg = w.gg, xx = w.xx, pipi = w.pipi, KK = 0.0, eta = 0.0, total = w.total)
end

# AxialVector A': A'→χχ̄ (β³).  No γγ or two-meson width (γγ↔χχ runs through
# the meson poles, not an A' pole).
# λ_a²m_A′β³/12π: the standard polarization-AVERAGED axial-vector width.  The
# draft and docs/notes/axialvector.txt carry /4π, missing the 1/3 over the three
# polarizations of a massive vector.
function mediator_widths(m::AxialVector)
    mphi, mx, gx = m.mphi, m.mx, m.gx
    xx = mphi > 2mx ? gx^2 * max(0.0, mphi^2 - 4mx^2)^1.5 / (12pi * mphi^2) : 0.0
    return (gg = 0.0, xx = xx, pipi = 0.0, KK = 0.0, eta = 0.0, total = xx)
end

