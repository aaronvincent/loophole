# Simple CMB ΔNeff: Boehm-Dolan-McCabe photon heating from species in contact
# with the plasma around neutrino decoupling, plus the φ that χ makes and leaves
# as free-streaming radiation.  Port of matlab/neff_simple.m.  Keep the two in step; test/validate_simple.jl checks
# them against each other.

using QuadGK: quadgk

"""
    neff_simple(m; T_RH=0.01, T_nu_dec=2e-3, Neff_SM=3.046, species=:phi_chi,
                T_BBN=1e-3, pigamma=true, primakoff=true) -> (dNeff, info)

Signed CMB-epoch ΔNeff.  Each species enters with a contact weight
`w = [1 - exp(-∫_{T_ν}^{T_RH} Γ/H dlnT)]·[1 - exp(-∫_{T_ν/100}^{T_ν} Γ/H dlnT)]`
(equilibrated by T_ν,dec, times still in contact after it) in

    N_eff = Neff_SM [1 + (7/22) Σ w_i (g_i/2) F_i(m_i/T_ν,dec)]^(-4/3),

φ through γγ↔φ and γγ↔φφ, χ through γγ↔χχ̄ or through φ.  φ that χ makes (χχ̄→φφ, χχ̄→φ)
relaxes toward (1-w_φ)ρ_φ^eq at a rate capped by χ's own photon link, and
free-streams; one entropy balance from T_ν,dec to recombination then replaces
(11/4)(1+X).  `species=:phi` drops χ.  CMB epoch only.

S⊗S / P⊗S with m_φ < 2m_π add π±γ ↔ π±φ to φ's link (`pigamma=false` drops
it).  A φ that equilibrates before T_ν,dec but loses contact by then has
decoupled with the neutrinos: its share `w_free = weq·(1 − wct)` joins the
free-streaming φ at T_ν,dec as ρ_φ^eq (ΔNeff = 4/7 if massless).
Non-axial models add e±γ ↔ e±φ and e⁺e⁻ ↔ γφ through the φγγ vertex
(`primakoff=false` drops them).  The φ still in contact after T_ν,dec stays in
the photon bath until T_late (optical depth from recombination up = 1), where
a still-relativistic φ is released as free radiation.
"""
function neff_simple(m::CrossSectionModel; T_RH::Real=0.01, T_nu_dec::Real=2e-3,
                     Neff_SM::Real=3.046, species::Symbol=:phi_chi, T_BBN::Real=1e-3,
                     pigamma::Bool=true, primakoff::Bool=true)
    species in (:phi, :phi_chi) || throw(ArgumentError("species must be :phi or :phi_chi"))
    mphi = mass_mediator(m); mx = mass_dm(m)
    g_chi = dof_dm(m); g_phi = dof_mediator(m)
    T_nu  = T_nu_dec
    zeta3 = 1.2020569
    Hof  = T -> 1.66 * sqrt(gstarofT(T, 0)) * T^2 / MPL
    ngam = T -> (2zeta3 / pi^2) * T^3

    # Rates on one grid, T_nu/100 .. T_RH, split at T_nu.
    Tb = _logspace(T_nu / 100, T_nu, 30)
    Ta = _logspace(T_nu, T_RH, 30)
    Tall = vcat(Tb, Ta[2:end]);  Hall = Hof.(Tall)
    nb = length(Tb);  ib = 1:nb;  ia = nb:length(Tall)
    weq = G -> 1 - exp(-_trapz(log.(Tall[ia]), G[ia] ./ Hall[ia]))   # equilibrated by T_nu
    wct = G -> 1 - exp(-_trapz(log.(Tall[ib]), G[ib] ./ Hall[ib]))   # contact after it
    wt  = G -> weq(G) * wct(G)

    # --- φ: inverse decay γγ ↔ φ, and γγ ↔ φφ ---
    Gw   = mediator_widths(m)
    Ginv = T -> Gw.gg * k1_over_k2(mphi / T)
    sigpp = _nsimple_sigpp(m, Gw)
    Gpp  = T -> _nsimple_ggphiphi(T, mphi, sigpp, g_phi)
    has_pg = pigamma && mphi < 2MPI_PM && (m isa ScalarSS || m isa ScalarPS)
    gq   = quark_couplings(m)
    Gpgf = T -> has_pg ? pigamma_phi_rate(T, mphi, gq) : 0.0
    g2_pr = (primakoff && !(m isa AxialVector)) ? _nsimple_gphigg2(m, Gw) : 0.0
    Gprf = T -> primakoff_phi_rate(T, mphi, g2_pr)
    Gid  = Ginv.(Tall);  Gp = Gpp.(Tall);  Gpg = Gpgf.(Tall);  Gpr = Gprf.(Tall)
    Gi   = Gid .+ Gp .+ Gpg .+ Gpr                   # φ's whole SM link
    r_inv = Gid[nb] / Hall[nb]
    r_pp  = Gp[nb] / Hall[nb]
    r_pr  = Gpr[nb] / Hall[nb]
    r_pg  = Gpg[nb] / Hall[nb]
    phi_coupled = Gi[nb] >= Hall[nb]
    w_phi  = wt(Gi)
    w_free = weq(Gi) * (1 - wct(Gi))                 # equilibrated, then decoupled before T_ν

    # Late decoupling of the bath-coupled φ: T_late where the optical depth below
    # T_ν, counted up from recombination, reaches 1 (continuous in the coupling).
    T_rec = 2.6e-10
    Tlo = _logspace(T_rec, T_nu, 60)
    rlo = (Ginv.(Tlo) .+ Gpp.(Tlo) .+ Gpgf.(Tlo) .+ Gprf.(Tlo)) ./ Hof.(Tlo)
    # Segments use the logarithmic mean of Γ/H (exact for a power law): the
    # trapezoid overshoots where the rate dies as e^{-m_e/T}, putting T_late 3-10% low.
    rm = map(rlo[1:end-1], rlo[2:end]) do r1, r2
        ra = 0.5 * (r1 + r2)
        (r1 > 0) != (r2 > 0) && return 0.0
        (r1 > 0 && r2 > 0 && abs(r2 - r1) > 1e-6 * ra) ? (r2 - r1) / log(r2 / r1) : ra
    end
    tau = vcat(0.0, cumsum(rm .* diff(log.(Tlo))))
    k1 = findfirst(>=(1), tau)
    T_late = if k1 === nothing
        T_nu
    elseif k1 == 1
        T_rec
    else
        exp(log(Tlo[k1-1]) + (1 - tau[k1-1]) / (tau[k1] - tau[k1-1]) * log(Tlo[k1] / Tlo[k1-1]))
    end
    late_rel = w_phi > 1e-10 &&
               _nsimple_rho(mphi, T_late, g_phi) > 1e-6 * _nsimple_rho(0.0, T_late, g_phi)

    # γγ→χχ̄ on both grids in ONE ggxxrate call.  Most of its cost is σ(s) on the
    # shared s-nodes, and each T's value is independent of the others in the call
    # (both grids end at T_RH, which sets s_max), so this is bit-identical to two
    # calls at ~2/3 the cost.
    Tg  = _logspace(T_BBN, T_RH, 40)
    Tq  = species == :phi_chi ? vcat(Tg, Tall) : Tg
    Rq  = max.(0.0, ggxxrate(m, Tq)) ./ ngam.(Tq)

    # "ever coupled" over [T_BBN, T_RH]: portal max(Γ₂₂, Γ_inv) vs H.
    G22 = Rq[1:length(Tg)]
    phi_ever_coupled = maximum(max.(G22, Ginv.(Tg) .+ Gpp.(Tg) .+ Gpgf.(Tg) .+ Gprf.(Tg)) ./ Hof.(Tg)) >= 1

    # --- χ ---
    r_22 = NaN; r_cp = NaN; w_chi = 0.0; w_chi_eq = 0.0; chi_coupled = false
    Gc = zeros(length(Tall))
    if species == :phi_chi
        Gc   = Rq[length(Tg)+1:end]
        r_22 = Gc[nb] / Hall[nb]
        if any(Gi .> Gc)                             # the link through φ can matter
            # per-χ rate: a χ annihilates on χ̄ only, n_χ̄ = n_eq/2 (g_χ = 4 counts both)
            Gcp = [0.5 * _nsimple_n(mx, T, g_chi, 1.0) * svavg_phiphi(m, T) for T in Tall]
            Gcp[.!isfinite.(Gcp)] .= 0.0             # min() would otherwise drop a NaN
            r_cp = Gcp[nb] / Hall[nb]
            Gc   = max.(Gc, min.(Gcp, Gi))
        end
        w_chi    = wt(Gc)
        w_chi_eq = weq(Gc)
        chi_coupled = Gc[nb] >= Hall[nb]
    end

    # --- heating formula (Boehm-Dolan-McCabe) ---
    y_phi = mphi / T_nu;  F_phi = _nsimple_F(y_phi, -1.0)
    y_chi = mx / T_nu;    F_chi = _nsimple_F(y_chi, 1.0)
    X = (7 / 22) * ((g_phi / 2) * F_phi * w_phi + (g_chi / 2) * F_chi * w_chi)
    dNeff_heat = Neff_SM * ((1 + X)^(-4 / 3) - 1)

    # --- φ made by χ: free-streaming dark radiation ---
    Q = 0.0; P_free = 1.0; dNeff_phichi = 0.0; Tphi_Tnu = 0.0
    TgTnu3 = (11 / 4) * (1 + X)
    chi_chain = species == :phi_chi && w_chi_eq * (1 - w_phi) > 1e-10
    if chi_chain || w_free > 1e-10 || late_rel
        T_end = min(T_nu / 100, mx / 50)
        late_rel && (T_end = min(T_end, T_late / 1.5))
        Tc  = _logspace(T_RH, T_end, 60)             # descending
        if late_rel && T_late < T_nu                 # release exactly at T_late, not a grid step below
            Tc = sort!(vcat(Tc[abs.(log.(Tc ./ T_late)) .> 1e-9], T_late); rev=true)
        end
        rq  = [_nsimple_rho(mphi, T, g_phi) for T in Tc]
        if chi_chain
            svE = [svavg_phiphi(m, T; moment=:energy) for T in Tc]
            svE[.!isfinite.(svE)] .= 0.0
            nx  = [_nsimple_n(mx,   T, g_chi,  1.0) for T in Tc]
            np  = [_nsimple_n(mphi, T, g_phi, -1.0) for T in Tc]
            # reaction density n_χ n_χ̄ ⟨σvE⟩ = nx²⟨σvE⟩/4: nx counts χ + χ̄
            CE  = nx .^ 2 ./ 4 .* max.(svE, 0.0) .+ Gw.xx * mphi .* np
            kc  = CE ./ Hof.(Tc) ./ rq               # relaxation rate per ln T
            kc[.!(rq .> 0)] .= 0.0                   # no φ reservoir left
            Gcl = exp.(_interp_linx(log.(Tall), log.(max.(Gc, floatmin(Float64))), log.(Tc)))
            kc  = min.(kc, Gcl ./ Hof.(Tc))          # chain bottleneck
        else
            kc  = zeros(length(Tc))
        end
        Q_free = w_free * _nsimple_rho(mphi, T_nu, g_phi) / T_nu^4
        # gs of the photon-contact bath
        gsb0 = [2 + 3.5 * _nsimple_F(0.511e-3 / T, 1.0) + w_chi * 3.5 * _nsimple_F(mx / T, 1.0)
                for T in Tc]
        gsp  = [(7 / 8) * g_phi * _nsimple_F(mphi / T, -1.0) for T in Tc]   # φ's share × w_bath
        sphi = q -> g_phi * (30q / (pi^2 * g_phi))^(3 / 4)
        S_i = NaN; Tnu_T = 1.0; w_bath = w_phi
        for i in 2:length(Tc)
            T = Tc[i]
            if T < T_nu && isnan(S_i)
                Q   = max(Q, Q_free)
                S_i = (11 / 2) * (1 + X) + sphi(Q)
            end
            if !isnan(S_i)
                Tnu_T = ((gsb0[i] + w_bath * gsp[i]) / max(S_i - sphi(Q), 1e-3))^(1 / 3)   # T_ν / T_γ
            end
            if late_rel && w_bath > 0 && T <= T_late && !isnan(S_i)
                # bath φ leaves as free radiation.  The bath gives up exactly the
                # entropy it carried for φ (w_bath gsp), not sphi(Q)'s thermal-shape
                # estimate of it, so S_i is rebased to keep S_i - sphi(Q) = bath.
                Qn  = Q + w_bath * rq[i] / (Tnu_T * T)^4
                S_i = S_i - sphi(Q) + sphi(Qn) - w_bath * gsp[i] / Tnu_T^3
                Q   = Qn
                w_bath = 0.0
            end
            Qeq = (1 - w_bath) * rq[i] / (Tnu_T * T)^4
            dlt = log(Tc[i-1] / T)
            Q   = Qeq + (Q - Qeq) * exp(-0.5 * (kc[i-1] + kc[i]) * dlt)
        end
        if isnan(S_i)
            Q = max(Q, Q_free); S_i = (11 / 2) * (1 + X) + sphi(Q)
        end
        # φ that decays (or is re-absorbed) into photons before recombination
        # returns its entropy: drop it from the free-streaming budget.
        Tl  = _logspace(T_rec, T_nu, 60)
        Gsm = Gw.gg + Gw.pipi + Gw.KK + Gw.eta
        P_free = exp(-_trapz(log.(Tl), Gsm .* k1_over_k2.(mphi ./ Tl) ./ Hof.(Tl)))
        Qf = P_free * Q
        TgTnu3   = (S_i - sphi(Qf)) / 2              # (T_γ/T_ν)³ at recombination
        Tphi_Tnu = (30Qf / (pi^2 * g_phi))^(1 / 4)
        y_rec = mphi / (Tphi_Tnu * T_rec / TgTnu3^(1 / 3))
        rho1nu_ref = (7 / 8) * (pi^2 / 15) * (4 / 11)^(4 / 3)
        if Qf > 0
            dNeff_phichi = Qf * _nsimple_radfrac(y_rec) / TgTnu3^(4 / 3) / rho1nu_ref
        end
    end
    Neff  = Neff_SM * ((4 / 11) * TgTnu3)^(-4 / 3) + dNeff_phichi
    dNeff = Neff - Neff_SM

    info = (phi_coupled = phi_coupled, w_phi = w_phi, chi_coupled = chi_coupled,
            w_chi = w_chi, w_chi_eq = w_chi_eq, phi_ever_coupled = phi_ever_coupled,
            Gamma_inv_over_H = r_inv, Gamma_pp_over_H = r_pp,
            Gamma_pg_over_H = r_pg, Gamma_prim_over_H = r_pr, w_phi_free = w_free,
            Gamma22_over_H = r_22, Gamma_chiphi_over_H = r_cp,
            y_phi = y_phi, y_chi = y_chi, F_phi = F_phi, F_chi = F_chi, X = X,
            dNeff_heat = dNeff_heat, dNeff_phi_rad = dNeff_phichi,
            Tphi_over_Tnu = Tphi_Tnu, P_free = P_free, T_late = T_late,
            late_release = late_rel, Neff = Neff)
    return dNeff, info
end

# Entropy of one species per massless Weyl-pair fermion (sgn +1 fermion, −1 boson).
function _nsimple_F(y, sgn)
    y > 100 && return 0.0
    f = u -> (4u^2 - y^2) * sqrt(max(0.0, u^2 - y^2)) / (exp(u) + sgn)
    return 30 / (7pi^4) * first(quadgk(f, y, y + 60; rtol=1e-8, atol=0.0))
end

# g_φγγ² from Γ_gg = g_φγγ² mφ³/(64π); at mφ = 0 from a tiny mφ (g_φγγ stays
# finite), rebuilding the model with the type's (mx, mphi, gx, g) constructor.
function _nsimple_gphigg2(m::CrossSectionModel, Gw)
    mphi = mass_mediator(m);  mref = max(mphi, 1e-12)
    if mphi < mref
        Gw = mediator_widths(typeof(m).name.wrapper(mass_dm(m), mref, coupling(m),
                                                    quark_couplings(m)))
    end
    return 64pi * Gw.gg / mref^3
end

# σ(γγ→φφ) as a function of s, or `nothing` if it vanishes.
# S⊗S / P⊗S: the charged-meson loop, boxes plus seagull triangle, through the
# low-energy theorem: σ = α²Y² s β/(36864π³),
# Y = Σ_j μ_j²/m_j⁴, μ_π = B0(g_u+g_d), μ_K = B0(g_u+g_s).  Photon exchange is
# O(α) smaller and dropped there.  Valid for √s ≪ 2m_π.
# Other models: t/u-channel photon exchange with g_φγγ² from Γ_gg.
function _nsimple_sigpp(m::CrossSectionModel, Gw)
    mphi = mass_mediator(m)
    if m isa Union{ScalarSS, ScalarPS}
        g = quark_couplings(m)
        Y = B0^2 * ((g[1] + g[2])^2 / MPI_PM^4 + (g[1] + g[3])^2 / MK^4)
        Y > 0 || return nothing
        return s -> _nsimple_sig_ggphiphi_box(s, mphi, Y)
    end
    g2 = _nsimple_gphigg2(m, Gw)
    g2 > 0 || return nothing
    return s -> _nsimple_sig_ggphiphi(s, mphi, g2)
end

# Per-φ rate of γγ ↔ φφ: 2R/n_φ^eq, R = T/(16π⁴) ∫ds σ s^{3/2} K₁(√s/T)
# (ggxxrate's normalisation); with s = (xT)², R = T⁶/(8π⁴) ∫_{x₀} σ x⁴ K₁(x) dx.
function _nsimple_ggphiphi(T, mphi, sig, g_phi)
    sig === nothing && return 0.0
    x0 = 2mphi / T
    x0 > 200 && return 0.0
    f = x -> sig((x * T)^2) * x^4 * _besselkx1(x) * exp(-(x - x0))   # Bessels.jl path (thermal_average.jl)
    I = first(quadgk(f, x0, x0 + 80; rtol=1e-6, atol=0.0))
    n = _nsimple_n(mphi, T, g_phi, -1.0)
    return n > 0 ? 2 * T^6 / (8pi^4) * exp(-x0) * I / n : 0.0
end

# σ(γγ→φφ), charged-meson loop (low-energy theorem, includes the ½ for φφ).
function _nsimple_sig_ggphiphi_box(s, m, Y)
    b = sqrt(max(0.0, 1 - 4m^2 / s))
    return ALPHA_EM^2 * Y^2 / (36864pi^3) * s * b
end

# σ(γγ→φφ) by photon exchange; r = m²/s.
function _nsimple_sig_ggphiphi(s, m, g2)
    r = m^2 / s
    r >= 0.25 && return 0.0
    b = sqrt(1 - 4r)
    L = r == 0 ? 0.0 : log((1 - 2r + b) / (2r))     # r² L → 0
    return g2^2 / (1024pi) * s * (b * (5 - 4r + 6r^2) - 8r^2 * (1 - 3r) * (1 - r) / (1 - 2r) * L)
end

function _nsimple_n(mm, T, g, sgn)
    y = mm / T
    y > 100 && return 0.0
    f = u -> u * sqrt(max(0.0, u^2 - y^2)) / (exp(u) + sgn)
    return g / (2pi^2) * T^3 * first(quadgk(f, y, y + 60; rtol=1e-6, atol=0.0))
end

# Equilibrium boson energy density.
function _nsimple_rho(mm, T, g)
    y = mm / T
    y > 100 && return 0.0
    f = u -> u^2 * sqrt(max(0.0, u^2 - y^2)) / (exp(u) - 1)
    return g / (2pi^2) * T^4 * first(quadgk(f, y, y + 60; rtol=1e-6, atol=0.0))
end

# Kinetic energy of a free-streaming, relativistically decoupled boson,
# f = 1/(e^{p/T}-1), over its massless value.  y = m/T.
function _nsimple_radfrac(y)
    y == 0 && return 1.0
    isfinite(y) || return 0.0
    f = u -> u^4 / (sqrt(u^2 + y^2) + y) / expm1(u)   # (E - m) u² f, no cancellation
    return first(quadgk(f, 0.0, 60.0; rtol=1e-6, atol=0.0)) / (pi^4 / 15)
end

_logspace(a, b, n) = 10.0 .^ range(log10(a), log10(b); length=n)

function _trapz(x, y)
    s = 0.0
    @inbounds for i in 2:length(x)
        s += 0.5 * (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    return s
end

# Linear interpolation with linear extrapolation (MATLAB interp1 'linear','extrap');
# x ascending.
function _interp_linx(x, y, xq)
    n = length(x)
    map(xq) do q
        j = clamp(searchsortedlast(x, q), 1, n - 1)
        y[j] + (y[j+1] - y[j]) * (q - x[j]) / (x[j+1] - x[j])
    end
end
