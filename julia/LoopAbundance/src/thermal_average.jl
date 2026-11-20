# Thermally averaged rates: the shared resonance-aware integrator
# (thermal_pole_integral), the Gondolo-Gelmini ⟨σv⟩ for χχ̄→γγ/all (svavg /
# svavg_grid), and the freeze-in γγ→χχ rate (ggxxrate).  Ports of
# thermal_pole_integral.m, svavg.m, svavg_grid.m, ggxxrate.m.
#
# Everything dispatches through the model interface (cross sections + poles),
# never on an interaction string.

using QuadGK: quadgk
import Bessels

# Bessel K₁/K₂ on the hot path (≈44k K₁ evaluations per relic solve) come from
# Bessels.jl, not SpecialFunctions: the latter routes a *real* argument through
# the complex Fortran AMOS kernel (133 ns/call vs 5.5 ns), which alone was ~72%
# of a scan point.  Agreement with AMOS is 3e-15 relative, and Bessels returns
# the correct asymptotic instead of throwing at large z, so the *accuracy* guard
# the AMOS path needed (scaled z>1e6 → √(π/2z), itself good only to 4e-7) is
# gone: besselk1x stays accurate past z = 1e12.
#
# The two branches that remain are not accuracy guards:
#   z > 700  unscaled K₁ is below realmin there, so 0 is exact to underflow, and
#            skipping the call is worth ~10% of a scan point — 87.6% of
#            ggxxrate's nodes land in this range.
#   z ≤ 0    Bessels throws a DomainError at 0 where SpecialFunctions returns
#            Inf.  A massless mediator (mphi = 0 → a pole at z = 0) must stay a
#            NaN-producing point, not an exception that kills the scan point.
#
# NB: only `besselk1`/`besselk1x` are collision-free names — `Bessels.besselkx`
# must stay qualified, since special.jl imports SpecialFunctions.besselkx into
# this same module.
_besselk1(z)  = z > 700 ? 0.0 : (z > 0 ? Bessels.besselk1(z) : Inf)
_besselk0(z)  = z > 700 ? 0.0 : (z > 0 ? Bessels.besselk0(z) : Inf)
_besselkx1(z) = z > 0 ? Bessels.besselk1x(z) : Inf
_besselkx0(z) = z > 0 ? Bessels.besselk0x(z) : Inf
_besselkx2(z) = Bessels.besselkx(2, z)

# ---------------------------------------------------------------------------
# MOMENT of the thermal s-integral.  `:number` is the zeroth moment, the usual
# reaction density / <sigma v> the relic solvers consume.  `:energy` is the
# first moment: each reaction is weighted by the lab-frame total energy of the
# incoming pair, whose mean at fixed s is
#
#     <E_tot>_s = sqrt(s) K0(sqrt(s)/T) / K1(sqrt(s)/T) + 2T
#
# (the total-energy distribution at fixed s goes as sqrt(E^2-s) exp(-E/T) for
# E >= sqrt(s)).  So the ONLY change is the Bessel kernel,
#
#     K1(sqrt(s)/T)  ->  sqrt(s) K0(sqrt(s)/T) + 2 T K1(sqrt(s)/T),
#
# with the same integrand, the same poles and the same narrow-width residues --
# hence a parameter rather than a second copy of the pole bookkeeping.  The
# Gondolo-Gelmini measure that produces the K1 kernel is unchanged by the
# initial-state mass, so the substitution is identical in svavg, svavg_phiphi
# and ggxxrate.
#
# `_bkern` uses UNSCALED Bessels (ggxxrate's convention); `_bkernx` uses SCALED
# ones (svavg's, where the exp((2mx - sqrt(s))/T) factor outside compensates and
# the K0/K1 scalings are the same exp(+x), so the combination is consistent).
@inline _bkern(sq, T, moment::Symbol) = moment === :number ?
    _besselk1(sq / T) : (sq * _besselk0(sq / T) + 2T * _besselk1(sq / T))
@inline _bkernx(sq, T, moment::Symbol) = moment === :number ?
    _besselkx1(sq / T) : (sq * _besselkx0(sq / T) + 2T * _besselkx1(sq / T))
_check_moment(moment::Symbol) = moment in (:number, :energy) ||
    throw(ArgumentError("moment must be :number or :energy"))

# Meson-pair / φφ / A'A' thresholds that kink the 'all' cross section; passed
# as continuum waypoints so the adaptive rule resolves the cusps (svavg.m:63-64).
const _SV_EXTRA_THRESH = (4 * MPI_PM^2, 4 * MPI0^2, 4 * M_ETA^2, 4 * MK^2, 4 * MK0^2,
                          (MPI0 + M_ETA)^2, 4 * M_ETAP^2)
# Model-specific additions: the three-body thresholds when enabled.
_model_extra_thresh(m::CrossSectionModel) = ()
_model_extra_thresh(m::PoleSPPP) = _threebody_thresholds(m.mphi)
_model_extra_thresh(m::AxialVector) =
    ((threepi_enabled() ? (9MPI0^2, (2MPI_PM + MPI0)^2) : ())...,
     (pipiphi_enabled() ? ((2MPI_PM + m.mphi)^2,) : ())...)

# Poles of the ⟨σv⟩ integrand for a channel: the model's, plus the A′ pole of
# the axial J = 1 χχ̄→3π piece in the channels that carry it (svavg.m sv_poles).
_sv_poles(m::CrossSectionModel, channel::Symbol) = poles(m)
function _sv_poles(m::AxialVector, channel::Symbol)
    ps = poles(m)
    if channel !== :gg && threepi_enabled() && m.mphi > 3MPI0
        ps = vcat(ps, [(m.mphi^2, m.mphi * mediator_widths(m).total)])
    end
    return ps
end

"""
    thermal_pole_integral(logintegrand, sig, kernel_at, s_lo, s_max, pole_s, pole_w;
                          RelTol=1e-6, AbsTol=1e-20, logdelta=1e-3,
                          narrow_tol=1e-5, extra_wp=())

Continuum + narrow-resonance integral in log(s).  Narrow poles (w_p/s_p <
narrow_tol) are notched out of the continuum and added analytically via the NWA
residue `real(sig(s_p))·w_p²·kernel_at(s_p)·π/w_p`; resolvable poles are left in
and straddled with a ±{1,5,50}-FWHM waypoint stack.  Port of
thermal_pole_integral.m.
"""
function thermal_pole_integral(logintegrand, sig, kernel_at, s_lo, s_max,
                               pole_s, pole_w;
                               RelTol=1e-6, AbsTol=1e-20, logdelta=1e-3,
                               narrow_tol=1e-5, extra_wp=())
    L_lo = log(s_lo)
    L_hi = log(s_max)
    I = 0.0
    L_hi <= L_lo && return I

    # --- classify poles ---
    notch = Tuple{Float64,Float64}[]            # (lo, hi) log-s windows
    wp_log = Float64[]
    for e in extra_wp
        (e > s_lo && e < s_max) && push!(wp_log, log(e))
    end
    # Graded breakpoints just above the lower limit, mirroring loggrid_edges, so
    # the adaptive rule resolves the thermal threshold boundary layer (width
    # ~ T/mx in s).  MATLAB's `integral` resolves this via its own subdivision;
    # QuadGK needs the hint, else a sharp spike at L_lo is sampled as ≈ 0.
    let span = L_hi - L_lo
        for d in exp10.(range(-4, log10(span); length=14))
            (0 < d < span) && push!(wp_log, L_lo + d)
        end
    end
    nwa_idx = Int[]
    for p in eachindex(pole_s)
        sp = pole_s[p]; wp = pole_w[p]
        Lp = log(sp)
        (sp > s_lo && Lp > L_lo && Lp < L_hi) || continue
        if wp > 0 && (wp / sp) < narrow_tol
            push!(notch, (Lp - logdelta, Lp + logdelta))
            push!(nwa_idx, p)
        else
            fw = max(wp, 1e-12 * sp)
            for c in (sp .+ (-50, -5, -1, 0, 1, 5, 50) .* fw)
                if c > s_lo
                    Lc = log(c)
                    (Lc > L_lo && Lc < L_hi) && push!(wp_log, Lc)
                end
            end
        end
    end

    # --- continuum on [L_lo, L_hi] minus the notched windows ---
    bounds = Float64[L_lo, L_hi]
    for (nl, nh) in notch
        push!(bounds, max(L_lo, nl)); push!(bounds, min(L_hi, nh))
    end
    bounds = sort(unique(bounds))
    for ii in 1:length(bounds)-1
        a = bounds[ii]; b = bounds[ii+1]
        b <= a && continue
        mid = 0.5 * (a + b)
        in_notch = any(nl < mid < nh for (nl, nh) in notch)
        in_notch && continue
        wseg = sort(filter(x -> a < x < b, wp_log))
        I += isempty(wseg) ?
             first(quadgk(logintegrand, a, b; rtol=RelTol, atol=AbsTol)) :
             first(quadgk(logintegrand, a, wseg..., b; rtol=RelTol, atol=AbsTol))
    end

    # --- narrow-width pole contributions ---
    for q in nwa_idx
        sp = pole_s[q]; wp = pole_w[q]
        I += real(sig(sp)) * wp^2 * kernel_at(sp) * pi / wp
    end
    return I
end

# split a model's poles into parallel (s_p, w_p) arrays
function _pole_arrays(m::CrossSectionModel, channel::Symbol)
    ps = _sv_poles(m, channel)
    return [p[1] for p in ps], [p[2] for p in ps]
end

# select the χχ̄ cross section for a channel
function _xx_sigma(m, channel)
    channel === :gg  && return s -> sigma_xxgg(m, s)
    channel === :had && return s -> sigma_xxhad(m, s)
    channel === :sm  && return s -> sigma_xxgg(m, s) + sigma_xxhad(m, s)
    channel === :all && return s -> sigma_xxall(m, s)
    throw(ArgumentError("channel must be :all, :sm, :gg or :had"))
end

"""
    svavg(m, T; channel=:all, moment=:number) -> ⟨σv⟩ [GeV⁻²]

Gondolo-Gelmini thermally averaged χχ̄ annihilation cross section at temperature
T.  `channel=:sm` is every SM final state, γγ + meson pairs (production side of
the relic solvers, by detailed balance with the equilibrium SM bath); `:gg` and
`:had` are its two pieces; `:all` is every open channel (depletion).  `moment=:energy`
returns `⟨σv E_tot⟩` in GeV⁻¹ instead — see `_bkernx`.  Port of svavg.m.
"""
function svavg(m::CrossSectionModel, T::Real; channel::Symbol=:all,
               moment::Symbol=:number)
    _check_moment(moment)
    mx = mass_dm(m); mphi = mass_mediator(m)
    sig = _xx_sigma(m, channel)
    K2scaled = _besselkx2(mx / T)
    integrand = s -> sig(s) * (s - 4mx^2) * sqrt(s) *
                     _bkernx(sqrt(s), T, moment) * exp((2mx - sqrt(s)) / T)
    logintegrand = ls -> integrand(exp(ls)) * exp(ls)

    s_lo = 4mx^2
    s_max = max(100 * T^2, 100 * mphi^2, 100 * mx^2, 1e2)

    pole_s, pole_w = _pole_arrays(m, channel)
    extra_wp = (_SV_EXTRA_THRESH..., 4 * mphi^2, _model_extra_thresh(m)...)
    kernel_at = sp -> (sp - 4mx^2) * sqrt(sp) *
                      _bkernx(sqrt(sp), T, moment) * exp((2mx - sqrt(sp)) / T)
    I = thermal_pole_integral(logintegrand, sig, kernel_at, s_lo, s_max,
                              pole_s, pole_w; RelTol=1e-5, AbsTol=1e-20,
                              narrow_tol=1e-5, extra_wp=extra_wp)
    return I / (8mx^4 * T * K2scaled^2)
end

# A pole goes to the narrow-width approximation when w_p/s_p = Γ_p/m_p is below
# this.  Tighter than thermal_pole_integral's 1e-5 default: the notch that comes
# with the NWA is a fixed ±logdelta window in log s, so it double-counts the
# Breit-Wigner tails lying outside it, at the ~2/(π·logdelta·m_p/Γ_p) level —
# 6e-4 at this cut, 1.5e-3 at the adaptive path's.  A fixed grid can afford to
# resolve everything wider instead of paying that.
const _NWA_CUT = 1e-6

# Classify a model's poles for the fixed-grid fast paths.  A narrow pole is
# notched out of the continuum (a ±logdelta window in log s) and returned in
# `nwa_s`/`nwa_w` for the analytic residue; a resolvable one stays in the
# continuum and gets a ladder of panel edges around it.  Poles outside
# (s_lo, s_max) are dropped.
function _classify_poles(ps, s_lo, s_max; logdelta=1e-3)
    L_extra = Float64[]
    notches = Tuple{Float64,Float64}[]
    nwa_s = Float64[]; nwa_w = Float64[]
    for (sp, wp) in ps
        (sp > s_lo && sp < s_max) || continue
        if wp > 0 && wp / sp < _NWA_CUT
            Lp = log(sp)
            push!(notches, (Lp - logdelta, Lp + logdelta))
            push!(L_extra, Lp - logdelta, Lp + logdelta)
            push!(nwa_s, sp); push!(nwa_w, wp)
        else
            # Half-decade ladder of ±offsets, from the width out to s_p itself.
            # A fixed ±100·w_p stack (what the single-pole path used) truncates
            # the Lorentzian tail, which carries ~2/(π·100) of the line — a 6e-3
            # error on a pole that dominates ⟨σv⟩.  The floor at _NWA_CUT·s_p
            # both bounds the ladder at ~12 rungs (anything narrower was notched
            # instead) and keeps a ZERO-width pole — Γ_total = 0 with every
            # channel closed — from looping forever on d = 0.
            d = max(wp, _NWA_CUT * sp)
            while d < sp
                for c in (sp - d, sp + d)
                    (s_lo < c < s_max) && push!(L_extra, log(c))
                end
                d *= sqrt(10)
            end
        end
    end
    return L_extra, notches, nwa_s, nwa_w
end

# Graded log-s panels (loggrid_edges) with fixed `npt`-point Gauss-Legendre
# nodes in each, skipping every panel that falls inside a notch window.  The
# membership test is on the panel MIDPOINT, so a threshold edge landing inside a
# notch — or two notches overlapping — still drops the right panels.
function _gl_panels(L_lo, L_hi, L_extra, notches; npt=24)
    edges = loggrid_edges(L_lo, L_hi, L_extra)
    xb, wb = gl_nodes_weights(npt)
    ls = Float64[]; wl = Float64[]
    for j in 1:length(edges)-1
        a = edges[j]; b = edges[j+1]
        b <= a && continue
        mid = (a + b) / 2
        any(nl < mid < nh for (nl, nh) in notches) && continue
        half = (b - a) / 2
        append!(ls, mid .+ half .* xb)
        append!(wl, half .* wb)
    end
    return exp.(ls), wl
end

# Shared ⟨σv⟩ node grid (channel-independent): graded log-s panels with
# meson/threshold/BW breakpoints and the notches of any narrow poles, then fixed
# 24-point Gauss-Legendre nodes per panel.  Returns the integration nodes
# (s_nodes, sq), panel weights (wl), and the NWA bookkeeping reused by the per-T
# loop.  Both svavg_grid and svavg_grid_pair build the integrand from these.
function _sv_pole_nodes(m::CrossSectionModel, Tgrid::AbstractVector)
    mx = mass_dm(m); mphi = mass_mediator(m)

    s_lo = 4mx^2
    Tmax = maximum(Tgrid)
    s_max = max(100 * Tmax^2, 100 * mphi^2, 100 * mx^2, 1e2)

    # Every threshold that kinks σ_all, i.e. the adaptive path's waypoint set.
    # The η/η′ pairs matter only for the SP/PP families, but costing two extra
    # panels everywhere is cheaper than branching on which channels a model has.
    bp = filter(x -> s_lo < x < s_max, (_SV_EXTRA_THRESH..., 4 * mphi^2, _model_extra_thresh(m)...))
    L_extra = collect(Float64, log.(bp))

    # :all is the superset: a :gg production channel loses only the A′ notch
    # window, where σ_γγ has no pole.
    Lp, notches, nwa_s, nwa_w = _classify_poles(_sv_poles(m, :all), s_lo, s_max)
    append!(L_extra, Lp)

    s_nodes, wl = _gl_panels(log(s_lo), log(s_max), L_extra, notches)
    return (s_nodes = s_nodes, sq = sqrt.(s_nodes), wl = wl, mx = mx,
            nwa_s = nwa_s, nwa_w = nwa_w, nwa_sq = sqrt.(nwa_s))
end

# The σ-independent half of the ⟨σv⟩ integrand on the shared nodes.
_sv_base(nd, sig) =
    nd.wl .* nd.s_nodes .* real.(sig.(nd.s_nodes)) .* (nd.s_nodes .- 4nd.mx^2) .* nd.sq

# Smooth (non-resonant) factor of the NWA residue at each narrow pole, at T.
@inline function _sv_nwa_smooth(nd, T, q)
    sqp = nd.nwa_sq[q]; mx = nd.mx
    return (nd.nwa_s[q] - 4mx^2) * sqp * _besselkx1(sqp / T) * exp((2mx - sqp) / T)
end

# Single-channel ⟨σv⟩(T) given the shared node grid `nd` and per-channel σ `sig`.
function _sv_from_nodes(nd, sig, Tgrid::AbstractVector)
    mx = nd.mx
    base = _sv_base(nd, sig)
    sig_pole = [real(sig(sp)) for sp in nd.nwa_s]
    sv = similar(Tgrid, Float64)
    for (k, T) in enumerate(Tgrid)
        kern = _besselkx1.(nd.sq ./ T) .* exp.((2mx .- nd.sq) ./ T)
        I = sum(base .* kern)
        for q in eachindex(nd.nwa_s)
            I += sig_pole[q] * nd.nwa_w[q] * pi * _sv_nwa_smooth(nd, T, q)
        end
        K2scaled = _besselkx2(mx / T)
        sv[k] = max(0.0, I / (8mx^4 * T * K2scaled^2))
    end
    return sv
end

"""
    svavg_grid(m, Tgrid; channel=:all) -> Vector

Vectorized ⟨σv⟩(T) over a whole grid, on the fixed Gauss-Legendre fast path:
σ(s) is evaluated once on a shared node grid and reused for every T.  Port of
svavg_grid.m.
"""
function svavg_grid(m::CrossSectionModel, Tgrid::AbstractVector; channel::Symbol=:all)
    nd = _sv_pole_nodes(m, Tgrid)
    return _sv_from_nodes(nd, _xx_sigma(m, channel), Tgrid)
end

"""
    svavg_grid_pair(m, Tgrid; production=:sm) -> (sv_all, sv_prod)

Both ⟨σv⟩ channels (`:all` depletion and the production channel, `:sm` or the
photon-only `:gg`) on one grid, sharing
the node grid **and** the per-T Bessel/exp kernel (the dominant cost) between
them and computing only the two cross sections separately — roughly halving the
work vs two `svavg_grid` calls.  Used by the relic solvers, which always need
both channels.
"""
function svavg_grid_pair(m::CrossSectionModel, Tgrid::AbstractVector; production::Symbol=:sm)
    nd = _sv_pole_nodes(m, Tgrid)
    mx = nd.mx
    sig_all = _xx_sigma(m, :all); sig_sm = _xx_sigma(m, production)
    base_all = _sv_base(nd, sig_all)
    base_sm  = _sv_base(nd, sig_sm)
    sigp_all = [real(sig_all(sp)) for sp in nd.nwa_s]
    sigp_sm  = [real(sig_sm(sp))  for sp in nd.nwa_s]
    sv_all = similar(Tgrid, Float64); sv_sm = similar(Tgrid, Float64)
    for (k, T) in enumerate(Tgrid)
        kern = _besselkx1.(nd.sq ./ T) .* exp.((2mx .- nd.sq) ./ T)   # shared
        I_all = sum(base_all .* kern); I_sm = sum(base_sm .* kern)
        for q in eachindex(nd.nwa_s)
            resid = nd.nwa_w[q] * pi * _sv_nwa_smooth(nd, T, q)
            I_all += sigp_all[q] * resid; I_sm += sigp_sm[q] * resid
        end
        K2scaled = _besselkx2(mx / T)
        denom = 8mx^4 * T * K2scaled^2
        sv_all[k] = max(0.0, I_all / denom); sv_sm[k] = max(0.0, I_sm / denom)
    end
    return sv_all, sv_sm
end

"""
    svavg_phiphi(m, T; moment=:number) -> ⟨σv⟩_{χχ̄→φφ} [GeV⁻²]

Thermally averaged χχ̄→φφ cross section (the smooth t-channel mediator-pair
piece, no s-channel pole).  Used by neff_simple for χ's link through φ and the
φ that χ makes.  Returns 0 when the channel is closed (mx < mφ) or the
Bessel average underflows.  Same Gondolo-Gelmini normalization (scaled Bessels)
as `svavg`.  `moment=:energy` returns `⟨σv E_tot⟩` in GeV⁻¹.  Port of svavg_phiphi.m.
"""
function svavg_phiphi(m::CrossSectionModel, T::Real; moment::Symbol=:number)
    _check_moment(moment)
    mx = mass_dm(m); mphi = mass_mediator(m)
    K2scaled = _besselkx2(mx / T)
    (isfinite(K2scaled) && K2scaled > 0) || return 0.0

    sig = s -> sigma_xxphiphi(m, s)
    integrand = s -> sig(s) * (s - 4mx^2) * sqrt(s) *
                     _bkernx(sqrt(s), T, moment) * exp((2mx - sqrt(s)) / T)
    logintegrand = ls -> integrand(exp(ls)) * exp(ls)

    s_lo = max(4mx^2, 4mphi^2)
    s_max = max(100 * T^2 + s_lo, 100 * mx^2, 1e2)
    s_max <= s_lo && return 0.0

    I = first(quadgk(logintegrand, log(s_lo), log(s_max); rtol=1e-4, atol=1e-30))
    return I / (8mx^4 * T * K2scaled^2)
end

"""
    sv_darkpair_grid(m, Tgrid) -> (lsv_fwd, lsv_rev)

Natural logs of ⟨σv⟩ for χχ̄ → dark pair (forward, Gondolo-Gelmini, = svavg) and
for its reverse, normalised by the mediator's K2:
⟨σv⟩_rev = (g_χ/g_φ)² I / (8 mφ⁴ T K₂(mφ/T)²), I = ∫ σ (s−4mx²)√s K₁(√s/T) ds over
s ≥ 4max(mx,mφ)², so ⟨σv⟩_rev Y_φeq² = ⟨σv⟩_fwd Y_χeq² exactly.  Each direction
stays O(1) where it is exothermic, at any m/T (a difference of ⟨σv⟩_all and
⟨σv⟩_sm cannot).  −Inf where zero.  √s = √s_lo + uT on 400 log-spaced u in
[1e-10, 80], trapezoid.  Port of sv_darkpair_grid.m.
"""
function sv_darkpair_grid(m::CrossSectionModel, Tgrid::AbstractVector)
    mx = mass_dm(m); mphi = mass_mediator(m)
    g_chi = dof_dm(m); g_phi = dof_mediator(m)
    nT = length(Tgrid)
    lf = fill(-Inf, nT); lr = fill(-Inf, nT)
    rs_lo = 2max(mx, mphi)
    u = vcat(1e-10, 10.0 .^ range(-9, log10(80); length=399))
    lz2k2(z) = z > 30 ? 2log(z) - z + 0.5log(pi / (2z)) + log1p(15 / (8z)) :
               z < 1e-3 ? log(2.0) : log(z^2 * besselk(2, z))
    for (k, T) in enumerate(Tgrid)
        F = similar(u)
        @inbounds for (i, ui) in enumerate(u)
            rs = rs_lo + ui * T; s = rs^2
            F[i] = sigma_xxdark(m, s) * (s - 4mx^2) * rs * _besselkx1(rs / T) * exp(-ui) * 2rs * T
        end
        It = 0.0
        @inbounds for i in 1:length(u)-1
            It += 0.5 * (F[i] + F[i+1]) * (u[i+1] - u[i])
        end
        (It > 0 && isfinite(It)) || continue
        lI = log(It) - rs_lo / T
        lf[k] = lI - log(8) - 5log(T) - 2lz2k2(mx / T)
        lr[k] = lI - log(8) - 5log(T) - 2lz2k2(mphi / T) + 2log(g_chi / g_phi)
    end
    return lf, lr
end

"""
    ggxxrate(m, Tvec; moment=:number) -> Vector

Thermally averaged γγ→χχ production rate Γ(T) (FI normalization), = T/(16π⁴)·∫
σ_γγ→χχ·s^{3/2}·K₁(√s/T) ds, in GeV⁴.  `moment=:energy` returns the energy
injection rate in GeV⁵ instead — see `_bkern`.  Fixed Gauss-Legendre nodes
shared across the whole T grid.  Port of ggxxrate.m / ggxxrate_fast.m.
"""
function ggxxrate(m::CrossSectionModel, Tvec::AbstractVector;
                  moment::Symbol=:number)
    _check_moment(moment)
    mx = mass_dm(m); mphi = mass_mediator(m)
    sig = s -> sigma_ggxx(m, s)

    s_lo = 4mx^2
    Tmax = maximum(Tvec)
    # 100mx² (as in svavg) keeps the window ≳ 8·max(T, mx) above threshold; without
    # it a heavy χ at low T (mx > 500 GeV at T_RH = 0.01) put s_lo above s_max.
    s_max = max(100 * Tmax^2, 100 * mphi^2, 100 * mx^2, 1e6)

    L_extra, notches, nwa_s, nwa_w = _classify_poles(poles(m), s_lo, s_max)

    s_nodes, wl = _gl_panels(log(s_lo), log(s_max), L_extra, notches)
    sq = sqrt.(s_nodes)
    base = wl .* real.(sig.(s_nodes)) .* s_nodes .^ 2.5             # N

    # Narrow-width residue: the same delta-function limit with the kernel at the
    # pole.  s^(3/2) → m_p³ is factored out, so only the Bessel combination
    # differs between the moments.
    nwa_sq = sqrt.(nwa_s)
    nwa_num = [real(sig(nwa_s[q])) * nwa_w[q] * pi * nwa_sq[q]^3
               for q in eachindex(nwa_s)]

    G = similar(Tvec, Float64)
    for (k, T) in enumerate(Tvec)
        KK = _bkern.(sq, T, moment)                   # unscaled (matches svavg's scalar path)
        I = sum(base .* KK)
        for q in eachindex(nwa_s)
            I += nwa_num[q] * _bkern(nwa_sq[q], T, moment)
        end
        G[k] = T / 16 / pi^4 * I
    end
    return G
end
