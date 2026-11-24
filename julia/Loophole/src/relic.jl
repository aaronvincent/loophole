# Relic-abundance solvers.  Ports of freeze_loop.m, freeze_out.m, Oh2_FI.m and
# freeze_coupled.m.  All take a CrossSectionModel (never an interaction string).
#
# MATLAB ode15s → OrdinaryDiffEqBDF.FBDF (variable-order BDF, the closest analog)
# with an analytic Jacobian; integration in τ = log(x) for the loop solvers to
# reduce stiffness.  Relic conversion uses the freeze_loop.m hard-coded
# s0 = 2891.2 cm⁻³ and ρ_crit/h² = 1.054e-5 GeV cm⁻³.

using OrdinaryDiffEqBDF: FBDF, QNDF
using SciMLBase: ODEProblem, ODEFunction, solve, successful_retcode, ContinuousCallback, terminate!, ReturnCode
using SpecialFunctions: besselk

const _SQRT_PI45 = sqrt(pi / 45)
# Y counts χ + χ̄ (g_χ = 4) and ⟨σv⟩ is the χχ̄ cross section, so
# ṅ = −(⟨σv⟩/2)(n² − n_eq²): every λ carries the Dirac factor 1/2.
const _LAM_PREF = 0.5 * _SQRT_PI45

_reached(sol, tend) = successful_retcode(sol) && sol.t[end] >= tend - 1e-6

# ⟨σv⟩(T) interpolant in log⟨σv⟩: the meson thresholds make the production
# ⟨σv⟩ ∝ e^(−2m_π/T), which plain pchip on the T grid misses by several %.
# Both channels use it so ⟨σv⟩_prod/⟨σv⟩_full stays exact where they coincide.
function _logpchip(x, y)
    itp = _pchip(x, log.(max.(y, floatmin(Float64))))
    return lt -> exp(itp(lt))
end
_omega_h2(mx, Y) = mx * S0_CM3 * Y / RHO_CRIT_H2

# =====================================================================
# freeze_loop — self-consistent FIMP↔WIMP scalar Boltzmann ODE
# =====================================================================
"""
    freeze_loop(m; T_RH=0.01, production=:sm) -> Ω h²

Self-consistent relic abundance from dY/dx = (1/x²)·corr·(g_s/√g_ρ)·(λ_prod·Y_eq²
− λ_full·Y²), integrated in τ = log(x) from T_RH through freeze-out.  NaN if the
integrator fails to reach x_final.  Port of freeze_loop.m.
"""
function freeze_loop(m::CrossSectionModel; T_RH::Real=0.01, ngrid::Int=200,
                     production::Symbol=:sm)
    mx = mass_dm(m)
    g_chi = dof_dm(m)

    x_init = mx / T_RH
    x_final = max(1000.0, 10 * x_init)
    T_min = mx / x_final

    Tgrid = 10.0 .^ range(log10(T_min), log10(T_RH); length=ngrid)
    logTgrid = log.(Tgrid)
    sv_all, sv_sm = svavg_grid_pair(m, Tgrid; production)
    svFull = _logpchip(logTgrid, sv_all)
    svProd = _logpchip(logTgrid, sv_sm)
    S = gstar_interp(T_min, T_RH)

    Yeq_of(x) = (45 / (4pi^4)) * (g_chi / S.gs(log(mx / x))) * x^2 * k2_of(x)

    function bc(x)
        T = mx / x; lT = log(T)
        pref = (1 + S.dlngs(lT) / 3) * (S.gs(lT) / sqrt(S.grho(lT))) / x
        cp = _LAM_PREF * MPL * mx * max(0.0, svProd(lT)) * pref
        cf = _LAM_PREF * MPL * mx * max(0.0, svFull(lT)) * pref
        return cp, cf, Yeq_of(x)
    end

    function f(u, p, t)
        cp, cf, Yeq = bc(exp(t))
        return cp * Yeq^2 - cf * u^2
    end
    function jac(u, p, t)
        _, cf, _ = bc(exp(t))
        return -2 * cf * u
    end

    # Y = 0 at T_RH, first step taken analytically (see freeze_loop.m): with
    # constant coefficients dY/dτ = a − cf·Y² gives Y = Y*·tanh(√(a·cf)·Δτ),
    # continuous in the coupling, so Ωh² is too.
    tau0 = log(x_init) + 1e-3 / max(1.0, x_init)
    dtau0 = tau0 - log(x_init)
    cp0, cf0, Yeq0 = bc(exp(tau0))
    a0 = cp0 * Yeq0^2
    Y0 = cf0 > 0 ? sqrt(a0 / cf0) * tanh(sqrt(a0 * cf0) * dtau0) : a0 * dtau0

    tspan = (tau0, log(x_final))
    # Y is unconstrained here and the annihilation term −cf·u² is even in u, so
    # there is no restoring force once a step undershoots below zero: −cf·u²
    # keeps driving Y down and the solve runs away to Ωh² ~ −1e5.  MATLAB's
    # ode15s does not undershoot at these points, so it is an FBDF robustness
    # issue, not a shared one -- freeze_loop.m needs nothing here.
    #
    # The fix is to reject the negative steps, NOT to tighten reltol.  Tightening
    # to 1e-8 does remove the negatives, but it makes FBDF fail to reach t_end at
    # a subset of points, which returns NaN -- fatal in a scan, because Diver
    # never falls back to a random vector for a seeded population slot (see
    # init.f90: the initial_guesses copy sits inside the retry loop), so one NaN
    # seed burns max_initialisation_attempts and hard-quits the run.  That bit
    # the mphi=10 scans, where mphi > 50*T_RH sends every point through here.
    sol = solve(ODEProblem(ODEFunction(f; jac=jac), Y0, tspan), FBDF();
                reltol=1e-4, abstol=1e-18, maxiters=10_000_000,
                isoutofdomain=(u, p, t) -> u < zero(u))
    _reached(sol, tspan[2]) || return NaN
    return _omega_h2(mx, sol.u[end])
end

# =====================================================================
# freeze_out — thermal WIMP via Steigman-Dasgupta-Beacom log-Y ODE
# =====================================================================
"""
    freeze_out(m) -> Ω h²

Equilibrium-IC freeze-out: dW/dx = (λ/x²)·corr·(g_s/√g_ρ)·(e^{2W_eq−W} − e^W),
W = ln Y, x = mx/T from 1 to 1000.  Port of freeze_out.m.
"""
function freeze_out(m::CrossSectionModel)
    mx = mass_dm(m)
    g_chi = dof_dm(m)

    Tgrid = 10.0 .^ range(log10(mx / 1000), log10(mx); length=200)
    svInterp = _pchip(log.(Tgrid), svavg_grid(m, Tgrid))
    S = gstar_interp(mx / 1000, mx)

    W0 = log(0.145 * g_chi / S.gs(log(mx))) + 1.5 * log(1.0) - 1.0

    function rhs(W, p, x)
        T = mx / x; lT = log(T)
        corr = 1 + S.dlngs(lT) / 3
        λ = _LAM_PREF * MPL * mx * max(0.0, svInterp(lT))
        W_eq = log(0.145 * g_chi / S.gs(lT)) + 1.5 * log(x) - x
        return (λ / x^2) * corr * (S.gs(lT) / sqrt(S.grho(lT))) * (exp(2W_eq - W) - exp(W))
    end

    sol = solve(ODEProblem(rhs, W0, (1.0, 1000.0)), FBDF();
                reltol=1e-10, abstol=1e-12, maxiters=10_000_000)
    _reached(sol, 1000.0) || return NaN
    return _omega_h2(mx, exp(sol.u[end]))
end

# =====================================================================
# Oh2_FI — linearized freeze-in integral
# =====================================================================
"""
    Oh2_FI(m, T_RH; production=:sm) -> Ω h²

Linearized freeze-in abundance: freeze_loop's equation without the Y² term.
With R = n_eq²⟨σv⟩/4 (`ggxxrate`'s normalisation, and that of the meson-pair
term from detailed balance) the Dirac-normalised production is dn/dt = 2R, so
Y = ∫ 2R/(sH)·(1 + dln g_s/3 dln T) dln T from max(T0, min(mx/200, T_RH/10³)) to
T_RH, and Ωh² = `_omega_h2(mx, Y)`.  The rate is interpolated in log.
`production=:gg` keeps photons only.  Port of Oh2_FI.m.
"""
function Oh2_FI(m::CrossSectionModel, T_RH::Real; production::Symbol=:sm)
    mx = mass_dm(m)
    # Below mx/200 the rate is < e⁻⁴⁰⁰ of its peak.
    T_lo = max(T0, min(mx / 200, T_RH / 1e3))
    Tgrid = 10.0 .^ range(log10(T_lo), log10(T_RH); length=300)
    R = ggxxrate(m, Tgrid)
    # ⟨σv⟩_had on a coarse grid above mx/60, interpolated in log.
    T_lo_h = max(T_lo, min(mx / 60, T_RH / 100))
    Th = 10.0 .^ range(log10(T_lo_h), log10(T_RH); length=60)
    sv_h = production === :sm ? svavg_grid(m, Th; channel=:had) : zeros(length(Th))
    if any(>(0), sv_h)
        lsv = _pchip(log.(Th), log.(max.(sv_h, floatmin(Float64))))
        for (k, T) in enumerate(Tgrid)
            T >= T_lo_h || continue
            neq = 4 / (2pi^2) * mx^2 * T * besselk(2, mx / T)
            R[k] += neq^2 * exp(lsv(log(T))) / 4
        end
    end
    lR = _pchip(log.(Tgrid), log.(max.(R, floatmin(Float64))))
    S = gstar_interp(T_lo, T_RH)
    dYdlnT = lT -> begin
        T = exp(lT)
        s = 2pi^2 / 45 * S.gs(lT) * T^3
        H = sqrt(4pi^3 / 45 * S.grho(lT)) * T^2 / MPL
        2 * exp(lR(lT)) / (s * H) * (1 + S.dlngs(lT) / 3)
    end
    Y = first(quadgk(dYdlnT, log(T_lo), log(T_RH); rtol=1e-8))
    return _omega_h2(mx, Y)
end

# =====================================================================
# freeze_coupled — coupled (Y_χ, Y_φ) Boltzmann ODE (4 branches)
# =====================================================================

# K₁/K₂ with freeze_coupled.m's exact asymptotics (distinct from special.jl's).
function _k1k2_fc(z)
    if z > 30
        return (1 - 3 / (8z)) / (1 + 15 / (8z))
    elseif z < 1e-3
        return z / 2
    else
        return besselk(1, z) / besselk(2, z)
    end
end

# Annihilation part of dW/dτ in log-Y form, cancellation-free.  The direct
# `cp*exp(2Weq - W) - cf*exp(W)` subtracts two terms that are individually ~1e31
# (cf ~ 1e32 once Γ/H ~ 1e29 and ⟨σv⟩ is BW-saturated); with no dark channel
# open every χχ channel is an SM one, so cp == cf identically and at W = Weq the RHS is pure
# roundoff (~1e15/e-fold, a few ulp), which the stiff solver cannot integrate --
# it aborts on step 1 and the caller sees NaN.  Factoring out cf·e^W turns the
# difference into one expm1, exact at cp == cf, W == Weq.  See the long note at
# freeze_coupled.m:`ann_logY`.
@inline function _ann_logY(cp, cf, W, Weq)
    # freeze_coupled's eliminated fast-φ coefficients can have either sign
    # (χχ̄ → φφ → 4χ makes the net Y² term negative): then nothing cancels and the
    # direct form is used.
    (cf > 0 && cp >= 0) ? cf * exp(W) * expm1(log(_snap_unit_ratio(cp / cf)) - 2 * (W - Weq)) :
                          cp * exp(2Weq - W) - cf * exp(W)
end

# ⟨σv⟩_prod/⟨σv⟩_full snapped to exactly 1 when the two agree to round-off.
# Whenever no dark channel is open they are the same quantity computed by
# different code paths (svavg_grid with/without the :sm filter, then separate pchip fits), so they
# differ in the last ulp -- and cf*exp(W) ~ 1e31 turns that into rhs ~ 1e15.
# 1e-12 is ~1e4 ulp yet ~1e11 below any real channel difference (S⊗S ~ 0.85).
@inline _snap_unit_ratio(r) = abs(r - 1) <= 1e-12 ? one(r) : r

"""
    freeze_coupled(m; T_RH=0.01, production=:sm) -> Ω h²

Coupled (Y_χ, Y_φ) Boltzmann solver for mφ ≲ T_RH (real φ a bath species).
Four runtime branches (adiabatic-φ × log-Y/lin-Y IC); adiabatic shortcut locks
Y_φ to Y_φ_eq when Γ_φ·K1/K2/H > 10³.  Port of freeze_coupled.m.

φ's 2→2 conversions on the bath, Γ_sc(T), join c_dec_eq (not c_dec_xx):
e±γ ↔ e±φ and e⁺e⁻ ↔ γφ (`primakoff_phi_rate`, non-axial models) and
π±γ ↔ π±φ (`pigamma_phi_rate`, S⊗S/P⊗S, mφ < 2m_π, production = :sm).  They
are tabulated on a fixed log-T grid (5 keV..T_RH) and enter the Y_φ IC test
and the adiabatic test (minimum over the range).  `ENV["LOOP_PHI_SC"] = "0"`
turns them off.  For axialVector, A′π ↔ ππ (`apipi_phi_rate`, production = :sm)
joins them (`ENV["LOOP_APIPI"] = "0"` alone off), and A′A′ ↔ γγ, meson pairs enters
Y_φ as −c_AA (Y_φ² − Y_φeq²) (`sv_AAsm_grid`; `ENV["LOOP_AA_SM"] = "0"` off).

`return_info = true` returns `(Ω h², info)`, info a NamedTuple with the mediator
left at x_final: `Y_phi_end`, `Omega_phi_h2`, `Gamma_phi`, `phi_stable`
(Γ t₀ < 1), `overclosed` (stable, Ω_φ h² > 0.12) and `tau_phi` [s].  As
freeze_coupled.m's second output.
"""
function freeze_coupled(m::CrossSectionModel; T_RH::Real=0.01, debug::Bool=false, ngrid::Int=200,
                        production::Symbol=:sm, return_info::Bool=false)
    Om, info = _freeze_coupled(m, T_RH, debug, ngrid, production)
    return return_info ? (Om, info) : Om
end

function _freeze_coupled(m::CrossSectionModel, T_RH::Real, debug::Bool, ngrid::Int,
                         production::Symbol)
    mx = mass_dm(m); mphi = mass_mediator(m)
    g_chi = dof_dm(m); g_phi = dof_mediator(m)

    x_init = mx / T_RH
    x_final = max(1000.0, 10 * x_init)
    T_min = mx / x_final

    Wd = mediator_widths(m)
    Gamma_gg = Wd.gg; Gamma_xx = Wd.xx
    Gamma_meson = Wd.pipi + Wd.KK + Wd.eta
    Gtot = Gamma_gg + Gamma_xx + Gamma_meson
    # mediator left at x_final (freeze_coupled.m phi_info)
    function phi_info(Yp)
        Ye = max(Yp, 0.0); Om = _omega_h2(mphi, Ye); st = Gtot * 6.61e41 < 1
        return (; Y_phi_end=Ye, Omega_phi_h2=Om, Gamma_phi=Gtot, phi_stable=st,
                overclosed=st && Om > 0.12, tau_phi=6.582e-25 / Gtot)
    end
    info_nan = (; Y_phi_end=NaN, Omega_phi_h2=NaN, Gamma_phi=NaN, phi_stable=false, overclosed=false,
                tau_phi=NaN)
    mes_eq = production === :sm          # mesons in the bath: φ ↔ ππ/KK inverse decay alive

    # φ ↔ bath 2→2 conversions (Primakoff, πγ→πφ) on a fixed, mx-independent grid.
    Tsc = 10.0 .^ range(log10(5e-6), log10(T_RH); length=120)
    Gsc = zeros(length(Tsc))
    use_sc = get(ENV, "LOOP_PHI_SC", "") != "0"
    if use_sc && !(m isa AxialVector)
        g2 = _nsimple_gphigg2(m, Wd)
        Gsc .+= primakoff_phi_rate.(Tsc, mphi, g2)
    end
    if use_sc && mes_eq && mphi < 2MPI_PM && (m isa ScalarSS || m isa ScalarPS)
        gq = quark_couplings(m)
        Gsc .+= pigamma_phi_rate.(Tsc, mphi, Ref(gq))
    end
    # axialVector: A′π → ππ, the A′'s tree-level link to the pion bath (apipi_phi_rate)
    if use_sc && mes_eq && m isa AxialVector && get(ENV, "LOOP_APIPI", "") != "0"
        Gsc .+= apipi_phi_rate.(Tsc, mphi, Ref(m.g))
    end
    has_sc = any(>(0), Gsc)
    Gsc_itp = has_sc ? _logpchip(log.(Tsc), Gsc) : (lt -> 0.0)
    Gamma_sc(T) = (!has_sc || T < Tsc[1]) ? 0.0 :
                  (v = Gsc_itp(log(min(T, Tsc[end]))); v <= 2floatmin(Float64) ? 0.0 : v)

    Tgrid = 10.0 .^ range(log10(T_min), log10(T_RH); length=ngrid)
    logTgrid = log.(Tgrid)
    # Depletion into SM final states (:sm) and production (production channel);
    # the dark pair's depletion is the same forward rate that feeds Y_φ (below),
    # so χχ̄ ↔ φφ conserves χ + φ exactly.
    sv_prod = svavg_grid(m, Tgrid; channel=production)
    sv_smd = production === :sm ? sv_prod : svavg_grid(m, Tgrid; channel=:sm)
    svFull = _logpchip(logTgrid, sv_smd)
    svProd = _logpchip(logTgrid, sv_prod)
    # χχ̄ ↔ φφ (φφ / A'A' / φ'φ'): only the reverse table; the forward rate is
    # derived from it in logs, cd = crev (Y_φeq/Y_χeq)², so detailed balance holds
    # to round-off (freeze_coupled.m).
    _, lsvR = sv_darkpair_grid(m, Tgrid)
    has_dark = any(isfinite, lsvR)
    lsvR = [isfinite(v) ? v : log(floatmin(Float64)) - 50 for v in lsvR]
    svDarkR = has_dark ? _pchip(logTgrid, lsvR) : (lt -> -Inf)
    # axialVector: A′A′ → γγ (charged-meson loop) and meson pairs (sv_AAsm_grid)
    has_AA = m isa AxialVector && get(ENV, "LOOP_AA_SM", "") != "0"
    if has_AA
        lsvA = sv_AAsm_grid(mphi, m.g, Tgrid)
        has_AA = any(isfinite, lsvA)
        lsvA = [isfinite(v) ? v : log(floatmin(Float64)) - 50 for v in lsvA]
    end
    svAA = has_AA ? _pchip(logTgrid, lsvA) : (lt -> -Inf)
    S = gstar_interp(T_min, T_RH)

    hubble(T) = T^2 / MPL * sqrt(8 * pi^3 * gstarofT(T, 0) / 90)
    Yeq_mass(x, mm, gi) = (45 / (4pi^4)) * (gi / S.gs(log(mx / x))) * z2k2(mm * x / mx)
    function logYeq_mass(x, mm, gi)
        z = mm * x / mx
        lz2k2 = z > 30 ? 2log(z) - z + 0.5log(pi / (2z)) + log1p(15 / (8z)) :
                z < 1e-3 ? log(2.0) : log(z^2 * besselk(2, z))
        return log(45 / (4pi^4)) + log(gi / S.gs(log(mx / x))) + lz2k2
    end
    # exp(l)·Yⁿ for Y > 0, in logs, capped; 0 otherwise
    ex(l, Y, n) = Y > 0 ? exp(min(l + n * log(Y), 700.0)) : 0.0

    adiabatic_phi = (Gtot * _k1k2_fc(mphi / T_RH) / hubble(T_RH)) > 1e3
    if has_sc && !adiabatic_phi
        # without Γ_sc the T_RH value is the minimum; with it, scan the range
        Tchk = 10.0 .^ range(log10(T_min), log10(T_RH); length=100)
        adiabatic_phi = minimum(T -> (Gtot * _k1k2_fc(mphi / T) + Gamma_sc(T)) / hubble(T), Tchk) > 1e3
    end
    get(ENV, "LOOP_FC_FULL2D", "") == "1" && (adiabatic_phi = false)   # tests

    # Every coefficient at x (freeze_coupled.m boltz_coef).  cp_s/cf_s: production
    # and SM depletion with the on-shell φ removed (φ is explicit); lp = log(Y_φeq/Y_χeq²).
    function bc(x)
        T = mx / x; lT = log(T)
        g_rho = S.grho(lT); g_s = S.gs(lT)
        corr = 1 + S.dlngs(lT) / 3
        lam = _LAM_PREF * MPL * mx * corr * (g_s / sqrt(g_rho)) / x
        cp = lam * max(0.0, svProd(lT)); cf = lam * max(0.0, svFull(lT))
        H_val = T^2 / MPL * sqrt(8 * pi^3 * g_rho / 90)
        K1K2 = _k1k2_fc(mphi / T)
        c_xx = corr * Gamma_xx * K1K2 / H_val
        c_SM = corr * (Gamma_gg + mes_eq * Gamma_meson) * K1K2 / H_val
        has_sc && (c_SM += corr * Gamma_sc(T) / H_val)
        c_mes = corr * (!mes_eq) * Gamma_meson * K1K2 / H_val
        lYc = logYeq_mass(x, mx, g_chi); lYp = logYeq_mass(x, mphi, g_phi)
        lp = lYp - 2lYc
        if has_dark
            lsr = svDarkR(lT)
            crev = lam * exp(min(lsr, 700.0))                    # φφ → χχ̄, on Y_φ²
            cd = lam * exp(min(lsr + 2(lYp - lYc), 700.0))       # χχ̄ → φφ, on Y_χ²
        else
            crev = 0.0; cd = 0.0
        end
        # A′A′ ↔ SM on (Y_φ² − Y_φeq²): 2λ = s/H, no Dirac ½
        c_AA = has_AA ? 2lam * exp(min(svAA(lT), 700.0)) : 0.0
        if c_xx > 0 && Gtot > 0
            Ron = exp(min(log(2c_xx) + lp, 700.0))
            R_p = Ron * (Gamma_gg + mes_eq * Gamma_meson) / Gtot
            R_f = Ron * (Gamma_gg + Gamma_meson) / Gtot
        else
            R_p = 0.0; R_f = 0.0
        end
        cp_s = max(cp - R_p, 0.0)
        cf_s = cp_s + max((cf - cp) - (R_f - R_p), 0.0)
        return (; cp_s, cf_s, cd, crev, c_xx, c_SM, c_mes, c_AA, lp, lYc, lYp,
                Yeqc = exp(lYc), Yeqp = exp(lYp))
    end

    # χ equation as cpE·Y_χeq² − cfE·Y_χ² at fixed Y_φ; P_phi = d(production)/dY_φ
    function chi_coefs(C, Yphi)
        lYc2 = -2C.lYc
        cpE = C.cp_s + (C.c_xx > 0 ? ex(log(2C.c_xx) + lYc2, Yphi, 1) : 0.0) +
                       (C.crev > 0 ? ex(log(C.crev) + lYc2, Yphi, 2) : 0.0)
        cfE = C.cf_s + C.cd + (C.c_xx > 0 ? exp(min(log(2C.c_xx) + C.lp, 700.0)) : 0.0)
        P_phi = 2C.c_xx + 2C.crev * max(Yphi, 0.0)
        # log of the production in Y units, as a log-sum: never divided by
        # Y_χeq², which underflows late while φ lingers (cpE above is capped there)
        l1 = C.cp_s > 0 ? log(C.cp_s) + 2C.lYc : -Inf
        l2 = (Yphi > 0 && C.c_xx > 0) ? log(2C.c_xx) + log(Yphi) : -Inf
        l3 = (Yphi > 0 && C.crev > 0) ? log(C.crev) + 2log(Yphi) : -Inf
        lm = max(l1, l2, l3)
        lP = isfinite(lm) ? lm + log(exp(l1 - lm) + exp(l2 - lm) + exp(l3 - lm)) : -Inf
        return cpE, cfE, P_phi, lP
    end
    # _ann_logY with the production as lP = log(cp·Y_χeq²) (freeze_coupled.m ann_logP)
    function ann_logP(lP, cf, W, Weq)
        if cf > 0 && isfinite(lP)
            lr = lP - 2Weq - log(cf)
            abs(lr) <= 1e-12 && (lr = 0.0)
            arg = lr - 2(W - Weq)
            return arg > 30 ? exp(lP - W) - cf * exp(W) : cf * exp(W) * expm1(arg)
        end
        return isfinite(lP) ? exp(lP - W) : -cf * exp(W)
    end
    # fast φ eliminated at its exact quasi-steady value: (Tφ, dTφ/dY_χ, δ), where
    # Tφ is everything χ exchanges with φ, from φ's balance (freeze_coupled.m
    # chi_exch_adia): vanishes at ρ = 1 through ρ² − 1 = expm1, nothing linearised
    function chi_exch_adia(C, Ychi)
        Yp = C.Yeqp
        q = C.crev + C.c_AA                        # c_AA (2Y_φeq δ + δ²) from A′A′ ↔ SM
        D1 = C.c_SM + C.c_xx + C.c_mes + 2q * Yp
        A = C.c_xx * Yp + C.crev * Yp^2
        # A(ρ²−1), crev Y_φeq²(ρ²−1): expm1 near equilibrium, A ρ² − A in logs far
        # from it (ρ² overflows once χ has frozen out)
        Crho2 = Ychi > 0 ? C.cd * Ychi^2 : 0.0
        Arho2 = Crho2 + (C.c_xx > 0 ? ex(log(C.c_xx) + C.lp, Ychi, 2) : 0.0)
        lr = Ychi > 0 ? 2(log(Ychi) - C.lYc) : -Inf
        if abs(lr) <= 1
            Ar = A * expm1(lr); Cr = C.crev * Yp^2 * expm1(lr)
        else
            Ar = Arho2 - A; Cr = Crho2 - C.crev * Yp^2
        end
        Bp = Ar - C.c_mes * Yp                     # δ from crev δ² + D′δ − B′ = 0
        if D1 > 0
            if q > 0
                sq = sqrt(max(D1^2 + 4q * Bp, 0.0))
                delta = 2Bp / (D1 + sq)
                fp = sq > 0 ? (2(D1 + sq) - 4q * Bp / sq) / (D1 + sq)^2 : 1 / D1
            else
                delta = Bp / D1; fp = 1 / D1
            end
        else
            delta = 0.0; fp = 0.0
        end
        if abs(delta) <= 0.5Yp
            # φ near its equilibrium (whatever χ does): δ form, 0 at ρ = 1
            G = C.c_SM + C.c_mes + C.crev * Yp + 2C.c_AA * Yp
            q2 = C.crev + 2C.c_AA
            T = -2G * delta - 2C.c_mes * Yp + Cr - q2 * delta^2
            dT = Ychi > 0 ? (-2G - 2q2 * delta) * fp * (2Arho2 / Ychi) + 2Crho2 / Ychi : 0.0
        else
            # φ largely empty: Y_φ itself, root of crev Y² + D Y − Q = 0 (no subtraction)
            D = C.c_SM + C.c_xx + C.c_mes
            Q = C.c_SM * Yp + Arho2 + C.c_AA * Yp^2
            sq = sqrt(D^2 + 4q * Q)
            Yph = D + sq > 0 ? 2Q / (D + sq) : Yp
            delta = Yph - Yp
            T = -2C.c_SM * delta - 2C.c_mes * Yph + Crho2 - C.crev * Yph^2 - 2C.c_AA * (Yph^2 - Yp^2)
            dT = (Ychi > 0 && sq > 0) ?
                 (-2C.c_SM - 2C.c_mes - 2(C.crev + 2C.c_AA) * Yph) * (2Arho2 / Ychi) / sq + 2Crho2 / Ychi : 0.0
        end
        return T, dT, delta
    end
    phi_rhs(C, Ychi, Yphi) = -C.c_SM * (Yphi - C.Yeqp) - C.c_xx * (Yphi - ex(C.lp, Ychi, 2)) -
                             C.c_mes * Yphi + C.cd * Ychi^2 - C.crev * max(Yphi, 0.0)^2 -
                             C.c_AA * (max(Yphi, 0.0)^2 - C.Yeqp^2)
    J22(C, Yphi) = -(C.c_SM + C.c_xx + C.c_mes) - 2(C.crev + C.c_AA) * max(Yphi, 0.0)

    # Y_χ = 0 at T_RH, first step analytic (freeze_coupled.m).  Every IC is built
    # at x0 = exp(tau0), the x the RHS reconstructs at t0.
    tau0 = log(x_init) + 1e-3 / max(1.0, x_init)
    dtau0 = tau0 - log(x_init)
    x0 = exp(tau0)
    # Y_φ IC: like χ, nothing before T_RH -- the same first step over dtau0,
    # Y_φeq (1 − exp(−Γ_SM dtau0/H)); continuous in the couplings (freeze_coupled.m)
    Gsm_RH = (Gamma_gg + mes_eq * Gamma_meson) * _k1k2_fc(mphi / T_RH) + Gamma_sc(T_RH)
    C0 = bc(x0)
    # A′A′ → SM run backwards fills the A′ at c_AA Y_φeq per e-fold (linearised)
    Y_phi_0 = Yeq_mass(x0, mphi, g_phi) * -expm1(-(Gsm_RH / hubble(T_RH) + C0.c_AA * C0.Yeqp) * dtau0)
    Yeq_chi0 = C0.Yeqc
    if adiabatic_phi
        cp0 = C0.cp_s; cf0 = C0.cf_s                   # φ exchange vanishes at ρ = 1
        a_lin = max(cp0 * Yeq_chi0^2 + chi_exch_adia(C0, 0.0)[1], 0.0)
        # fast dark exchange (χχ̄ → φφ → 4χ for mφ > 2mx) saturates any seed within
        # dtau0: start saturated (freeze_coupled.m)
        dark_sat = a_lin > 0 && abs(chi_exch_adia(C0, a_lin * dtau0)[2]) * dtau0 > 1
    else
        dark_sat = false
        _, cf0, _, lP0 = chi_coefs(C0, Y_phi_0)
        a_lin = exp(lP0)
        cp0 = exp(min(lP0 - 2C0.lYc, 700.0))          # seed only (ratio, tanh argument)
    end
    if cf0 > 0
        ratio = _snap_unit_ratio(cp0 / cf0)
        seed_f = tanh(sqrt(cp0 * cf0) * Yeq_chi0 * dtau0)
        # A seed within 1e-10 of Y* IS the steady state (QNDF cannot start from
        # the ~1e-9 residual otherwise).
        1 - seed_f < 1e-10 && (seed_f = 1.0)
        use_logY = seed_f > 0.5
        W_chi_0 = C0.lYc + 0.5 * log(min(1.0, ratio)) + log(seed_f)
        Y_chi_0 = sqrt(a_lin / cf0) * tanh(sqrt(a_lin * cf0) * dtau0)
    else
        use_logY = false
        Y_chi_0 = a_lin * dtau0
    end
    # fast φ whose exchange alone carries χ past half its equilibrium within dtau0
    # (or the dark 2→4 runaway): saturated seed (freeze_coupled.m)
    if adiabatic_phi && (dark_sat || a_lin * dtau0 > 0.5Yeq_chi0)
        use_logY = true
        W_chi_0 = C0.lYc + (cf0 > 0 ? 0.5 * log(min(1.0, _snap_unit_ratio(cp0 / cf0))) : 0.0)
    end
    Y0vec = use_logY ? [W_chi_0, Y_phi_0] : [Y_chi_0, Y_phi_0]

    # Dark-exchange shortcut (freeze_coupled.m): φ's own rates slow, χχ̄ ↔ φφ far
    # faster than H (axialVector with mx > mφ/2 has Γ(A′) = 0 and ~1e20 H).  χ and φ
    # sit in relative chemical equilibrium, Y_φ = r Y_χ, r = Y_φeq/Y_χeq; only
    # N = Y_χ + Y_φ evolves slowly.  Stage 1 integrates log N until the exchange
    # falls below 1e3 per e-fold, stage 2 continues the 2×2 system from the split.
    exch_fast = has_dark && !adiabatic_phi && get(ENV, "LOOP_FC_NOEXCH", "") != "1" &&
                C0.cd * Yeq_chi0 > 1e9 && C0.crev * C0.Yeqp > 1e9
    if exch_fast
        function split_exch(C, WN)
            lrr = C.lYp - C.lYc
            l1r = lrr > 30 ? lrr + log1p(exp(-lrr)) : log1p(exp(lrr))
            return exp(WN - l1r), exp(WN + lrr - l1r), l1r
        end
        rhsN = (u, p, t) -> begin
            C = bc(exp(t)); Yc, Yp_, l1r = split_exch(C, u)
            X = C.c_xx * (Yp_ - ex(C.lp, Yc, 2)) - C.c_SM * (Yp_ - C.Yeqp) - C.c_mes * Yp_ -
                C.c_AA * (Yp_^2 - C.Yeqp^2)
            _ann_logY(C.cp_s * exp(-2l1r), C.cf_s * exp(-2l1r), u, C.lYc + l1r) + X * exp(-u)
        end
        jacN = (u, p, t) -> begin
            C = bc(exp(t)); Yc, Yp_, l1r = split_exch(C, u)
            # φ's own links, d(X e^-u)/du = (dX/du − X) e^-u: the SM link can be fast
            # (Γ_sc ~ 1e2 H) without φ being eliminated, and is stiff then
            Xc = ex(C.lp, Yc, 2)
            X  = C.c_xx * (Yp_ - Xc) - C.c_SM * (Yp_ - C.Yeqp) - C.c_mes * Yp_ -
                 C.c_AA * (Yp_^2 - C.Yeqp^2)
            dX = C.c_xx * (Yp_ - 2Xc) - (C.c_SM + C.c_mes) * Yp_ - 2C.c_AA * Yp_^2
            -C.cp_s * exp(2C.lYc - u) - C.cf_s * exp(u - 2l1r) + (dX - X) * exp(-u)
        end
        # N's own first step, around its equilibrium Y_χeq (1+r)
        _, _, l1r0 = split_exch(C0, 0.0)
        if C0.cf_s > 0
            fN = tanh(sqrt(C0.cp_s * C0.cf_s) * Yeq_chi0 * dtau0); 1 - fN < 1e-10 && (fN = 1.0)
            WN0 = C0.lYc + l1r0 + 0.5log(min(1.0, _snap_unit_ratio(C0.cp_s / C0.cf_s))) +
                  log(max(fN, floatmin(Float64)))
        else
            WN0 = log(max(C0.cp_s * Yeq_chi0^2 * dtau0 + Y_phi_0, floatmin(Float64)))
        end
        condN = (u, t, integ) -> begin
            C = bc(exp(t)); Yc, Yp_, _ = split_exch(C, u)
            log(max(min(C.cd * Yc, C.crev * Yp_), floatmin(Float64))) - log(1e3)
        end
        cbN = ContinuousCallback(condN, nothing, integ -> terminate!(integ))   # on down-crossing
        solN = solve(ODEProblem(ODEFunction(rhsN; jac=jacN), WN0, (tau0, log(x_final))), QNDF();
                     reltol=1e-6, abstol=1e-10, maxiters=1_000_000, callback=cbN)
        Yc, Yp_, _ = split_exch(bc(exp(solN.t[end])), solN.u[end])
        if solN.retcode != ReturnCode.Terminated
            _reached(solN, log(x_final)) || return (NaN, info_nan)
            return (_omega_h2(mx, Yc), phi_info(Yp_))      # fast to the end
        end
        rhs2 = (u, p, t) -> begin
            C = bc(exp(t)); _, cfE, _, lP = chi_coefs(C, u[2])
            [ann_logP(lP, cfE, u[1], C.lYc), phi_rhs(C, exp(u[1]), u[2])]
        end
        jac2 = (u, p, t) -> begin
            C = bc(exp(t)); _, cfE, P_phi, lP = chi_coefs(C, u[2]); W = u[1]
            [(-exp(lP - W) - cfE * exp(W))  (P_phi * exp(-W));
             (2C.c_xx * ex(C.lp, exp(W), 2) + 2C.cd * exp(2W))  J22(C, u[2])]
        end
        sol2 = solve(ODEProblem(ODEFunction(rhs2; jac=jac2), [log(Yc), Yp_], (solN.t[end], log(x_final))),
                     QNDF(); reltol=1e-6, abstol=1e-10, maxiters=1_000_000)
        _reached(sol2, log(x_final)) || return (NaN, info_nan)
        return (_omega_h2(mx, exp(sol2.u[end][1])), phi_info(sol2.u[end][2]))
    end

    tspan = (tau0, log(x_final))

    # assemble (rhs, jac, u0, reltol, abstol, extract) for the active branch
    local rhsf, jacf, u0, rtol, atol, extract
    if adiabatic_phi && use_logY
        rhsf = (u, p, t) -> begin
            C = bc(exp(t)); T, _, _ = chi_exch_adia(C, exp(u))
            _ann_logY(C.cp_s, C.cf_s, u, C.lYc) + T * exp(-u)
        end
        jacf = (u, p, t) -> begin
            C = bc(exp(t)); T, dT, _ = chi_exch_adia(C, exp(u))
            -C.cp_s * exp(2C.lYc - u) - C.cf_s * exp(u) + dT - T * exp(-u)
        end
        # rtol 1e-8: the exact-minus-linear remainder varies fast near T_RH and is not
        # in the Jacobian; QNDF at 1e-6 under-resolves it (×0.46 at heavyScalar
        # mx = 1 MeV, mφ = 10 MeV, gx = 1), ode15s at 1e-6 does not.
        u0 = Y0vec[1]; rtol = 1e-8; atol = 1e-10; extract = s -> exp(s.u[end])
    elseif adiabatic_phi
        rhsf = (u, p, t) -> begin
            C = bc(exp(t))
            C.cp_s * C.Yeqc^2 - C.cf_s * u^2 + chi_exch_adia(C, u)[1]
        end
        jacf = (u, p, t) -> begin
            C = bc(exp(t)); _, dT, _ = chi_exch_adia(C, u)
            -2C.cf_s * u + dT
        end
        u0 = Y0vec[1]; rtol = 1e-4; atol = 1e-18; extract = s -> s.u[end]
    elseif use_logY
        rhsf = (u, p, t) -> begin
            C = bc(exp(t)); _, cfE, _, lP = chi_coefs(C, u[2])
            [ann_logP(lP, cfE, u[1], C.lYc), phi_rhs(C, exp(u[1]), u[2])]
        end
        jacf = (u, p, t) -> begin
            C = bc(exp(t)); _, cfE, P_phi, lP = chi_coefs(C, u[2]); W = u[1]
            [(-exp(lP - W) - cfE * exp(W))  (P_phi * exp(-W));
             (2C.c_xx * ex(C.lp, exp(W), 2) + 2C.cd * exp(2W))  J22(C, u[2])]
        end
        u0 = Y0vec; rtol = 1e-6; atol = 1e-10; extract = s -> exp(s.u[end][1])
    else
        rhsf = (u, p, t) -> begin
            C = bc(exp(t)); _, cfE, _, lP = chi_coefs(C, u[2])
            [exp(lP) - cfE * u[1]^2, phi_rhs(C, u[1], u[2])]
        end
        jacf = (u, p, t) -> begin
            C = bc(exp(t)); _, cfE, P_phi, _ = chi_coefs(C, u[2])
            [(-2cfE * u[1])  P_phi;
             (2C.c_xx * ex(C.lp, u[1], 1) + 2C.cd * u[1])  J22(C, u[2])]
        end
        u0 = Y0vec; rtol = 1e-4; atol = 1e-18; extract = s -> s.u[end][1]
    end

    prob = ODEProblem(ODEFunction(rhsf; jac=jacf), u0, tspan)
    # maxiters 1e6: a closed, stiff dark sector (axialVector, Γ(A′) = 0, gx ~ 1) is not
    # integrable here and QNDF took minutes to report it at 1e7; converging points
    # take ≲ 1e4 steps.
    sol = solve(prob, QNDF(); reltol=rtol, abstol=atol, maxiters=1_000_000)
    # QNDF (the ode15s/NDF analog) can still abort on a stiff non-steady IC
    # transient; retry with a tiny seed dt (lets it ramp like ode15s does).
    if !_reached(sol, tspan[2])
        sol = solve(prob, QNDF(); reltol=rtol, abstol=atol, maxiters=1_000_000,
                    dt=1e-12, force_dtmin=true)
        # A forced-dt run that "reaches the end" in a couple of steps has not
        # integrated anything: it returns the IC.  Report failure instead.
        length(sol.t) < 5 && return (NaN, info_nan)
    end
    if debug
        @info "freeze_coupled" adiabatic_phi use_logY retcode=sol.retcode tend=sol.t[end] tgoal=tspan[2] u0 seed_z=(cf0 > 0 ? sqrt(cp0 * cf0) * Yeq_chi0 * dtau0 : NaN) rhs0=rhsf(u0, nothing, tau0) nsteps=length(sol.t)
    end
    _reached(sol, tspan[2]) || return (NaN, info_nan)
    Ychi_end = extract(sol)
    Yphi_end = adiabatic_phi ? bc(x_final).Yeqp + chi_exch_adia(bc(x_final), Ychi_end)[3] :
                               sol.u[end][2]
    return (_omega_h2(mx, Ychi_end), phi_info(Yphi_end))
end
