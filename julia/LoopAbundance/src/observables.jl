# DM-nucleon cross section and the single-point observable evaluator.
# Ports of DMnucxsec.m and model_observables.m.

# --- nucleon inputs (DMnucxsec.m) --------------------------------------
const _DU = 0.842     # proton axial charges Δq^p
const _DD = -0.427
const _DS = -0.085
const _FPU = 0.02     # proton scalar couplings (constants.jl has the same;
const _FPD = 0.026    # kept local to mirror DMnucxsec.m's own values)
const _FPS = 0.118
const _MUP = 0.0022   # quark masses [GeV]
const _MDN = 0.00469
const _MS  = 0.095
const _VREL = 220 / 3e5

# scalar nucleon form factor h_p, axial nucleon coupling a_N, reduced mass μ.
# (_hp, _htildeP below, and _a_N here are the single source for the nucleon
# couplings — mirror of nucleon_couplings.m.)
_hp(g) = MN * (g[1] * _FPU / _MUP + g[2] * _FPD / _MDN + g[3] * _FPS / _MS)
# Scalar NEUTRON coupling: isospin swaps the u<->d sigma-term form factors, each
# 1/m_i kept with its own quark.  So this is NOT _hp([g[2], g[1], g[3]]).
_hn(g) = MN * (g[1] * _FPD / _MUP + g[2] * _FPU / _MDN + g[3] * _FPS / _MS)
_a_N(g) = g[1] * _DU + g[2] * _DD + g[3] * _DS
_mu_red(mx) = MN * mx / (MN + mx)

# Pseudoscalar nucleon coupling h̃_P (used by the SP/PP interactions in place of
# the scalar h_p): proton spin content Δn_i weighted by 1/m_i with the
# SU(3)-singlet (η′) piece projected out, (g_i − ⟨g/m⟩/⟨1/m⟩).  Port of
# DMnucxsec.m's htildeP.  Note these spin components differ from the axial
# charges _DU/_DD/_DS used by the axialVector σ_N.
const _DNU = 0.84      # proton quark spin components Δn_i
const _DND = -0.44
const _DNS = -0.03
function _htildeP(g)
    invm_sum = 1 / _MUP + 1 / _MDN + 1 / _MS
    g_proj   = (g[1] / _MUP + g[2] / _MDN + g[3] / _MS) / invm_sum   # singlet projection
    return MN * (_DNU / _MUP * (g[1] - g_proj) +
                 _DND / _MDN * (g[2] - g_proj) +
                 _DNS / _MS  * (g[3] - g_proj))
end

# Pseudoscalar NEUTRON coupling h̃_N: isospin swaps the u and d SPIN CONTENTS
# only — each 1/m_i stays with its own quark, and g_proj (a property of the
# coupling direction, not of the nucleon) is unchanged.  So h̃_N is NOT
# _htildeP([g2, g1, g3]).  Used by the pseudoscalar 5th-force bound, which
# constrains g_P² and |g_P g_N| separately.  Mirror of nucleon_couplings.m.
function _htildeN(g)
    invm_sum = 1 / _MUP + 1 / _MDN + 1 / _MS
    g_proj   = (g[1] / _MUP + g[2] / _MDN + g[3] / _MS) / invm_sum
    return MN * (_DND / _MUP * (g[1] - g_proj) +
                 _DNU / _MDN * (g[2] - g_proj) +
                 _DNS / _MS  * (g[3] - g_proj))
end

# =====================================================================
# sigma_nucleon — DM-nucleon cross section [GeV⁻²].  Always the FULL angular-
# integrated form (valid across the whole mediator-mass range); the formula
# differs only by vertex structure (SS/PS/SP/PP/axial).  There is NO heavy/light
# split here — HeavyScalar* and LightScalar* give the same full σ.  The MATLAB
# DMnucxsec.m additionally offers heavy- and light-mediator approximations selected by
# name ('heavyScalar*'/'lightScalar*'); the canonical/scan σ is the full form, here.
# =====================================================================
# Scalar/pseudoscalar χ-coupling × scalar/pseudoscalar nucleon-coupling, valid across
# the FULL mediator-mass range (no heavy/light propagator split). Integrate
#   dσ/dΩ = gx² hN²/(64π²) (a−t)(b−t)/(s (mφ²−t)²)
# over solid angle at reference _VREL (t = −q², p_cm = μ·v, t ∈ [−T,0], T = 4μ²v²):
#   σ = gx² hN²/(64π μ²v² s) · mφ² · J(x,α,β),  x=T/mφ², α=a/mφ², β=b/mφ².
# Each vertex sets one numerator factor — a pseudoscalar χ zeroes α (the 4mx² term),
# a pseudoscalar nucleon zeroes β (the 4mN² term) and swaps hp→htildeP:
#   SS:  α=4mx²/mφ², β=4mN²/mφ², hp      → μ²/mφ⁴ heavy; 1/mφ² light
#   PS (ps-χ × scalar-N):  α=0, β=4mN²/mφ², hp      → ~v² (p-wave) heavy; ~ln(1/mφ²) light
#   PP (ps × ps):          α=0, β=0,         htildeP → ~v⁴ heavy; mφ-independent light
# Port of DMnucxsec.m.
function _sigma_scalar_core(m, al, be, hN)
    mu = _mu_red(m.mx)
    s  = (m.mx + MN)^2
    x  = 4 * mu^2 * _VREL^2 / m.mphi^2
    I  = m.mphi^2 * _scalar_angint(x, al, be)
    return m.gx^2 * hN^2 / (64pi * mu^2 * _VREL^2 * s) * I
end
_sigma_scalar(m) = _sigma_scalar_core(m, 4*m.mx^2/m.mphi^2, 4*MN^2/m.mphi^2, _hp(m.g))
sigma_nucleon(m::HeavyScalar) = _sigma_scalar(m)
sigma_nucleon(m::LightScalar) = _sigma_scalar(m)

# Dimensionless angular integral J(x,α,β) = x + (α+β−2)ln(1+x) + (1−α−β+αβ)x/(1+x)
# = J0 + (α+β)J1 + αβ J2, with J0,J1,J2 the cancellation-safe t-moments from _jmoments.
_scalar_angint(x, al, be) = (j = _jmoments(x); j[1] + (al + be) * j[2] + al * be * j[3])

# Cancellation-safe moments of 1/(mφ²−t)² over t∈[−T,0] (x=T/mφ²): J0=x−2ln(1+x)+x/(1+x)
# (~x³/3), J1=ln(1+x)−x/(1+x) (~x²/2), J2=x/(1+x). J0,J1 cancel to O(x³)/O(x²), so use
# their convergent series below x=0.05 to keep full precision at heavy mφ (x~1e-9);
# J2 is cancellation-free. Shared by the scalar family and axialVector. See DMnucxsec.m.
function _jmoments(x)
    J2 = x / (1 + x)
    if x < 0.05
        J0 = sum((-1)^(n+1) * (n-2)/n * x^n for n in 3:12)   # x³/3 − x⁴/2 + …
        J1 = sum((-1)^n     * (n-1)/n * x^n for n in 2:12)   # x²/2 − 2x³/3 + …
    else
        L  = log1p(x)
        J0 = x - 2L + J2
        J1 = L - J2
    end
    return (J0, J1, J2)
end
# PS: pseudoscalar χ × scalar nucleon → α=0, scalar nucleon coupling hp.
sigma_nucleon(m::HeavyScalarPS) = _sigma_scalar_core(m, 0.0, 4*MN^2/m.mphi^2, _hp(m.g))
sigma_nucleon(m::LightScalarPS) = _sigma_scalar_core(m, 0.0, 4*MN^2/m.mphi^2, _hp(m.g))
# SP: scalar χ × pseudoscalar nucleon → β=0, pseudoscalar nucleon coupling htildeP.
sigma_nucleon(m::HeavyScalarSP) = _sigma_scalar_core(m, 4*m.mx^2/m.mphi^2, 0.0, _htildeP(m.g))
sigma_nucleon(m::LightScalarSP) = _sigma_scalar_core(m, 4*m.mx^2/m.mphi^2, 0.0, _htildeP(m.g))
# PP: pseudoscalar χ × pseudoscalar nucleon → α=β=0, pseudoscalar nucleon coupling htildeP.
sigma_nucleon(m::Union{HeavyScalarPP,ScalarPP}) = _sigma_scalar_core(m, 0.0, 0.0, _htildeP(m.g))
sigma_nucleon(m::LightScalarPP)                 = _sigma_scalar_core(m, 0.0, 0.0, _htildeP(m.g))
# AxialVector: spin-dependent, full mediator-mass range.  a_N = Σ_q a_q Δ_q^p.
# Integrate dσ/dΩ = λ_a² a_N²/(32π² m_A'⁴ s (m_A'²−t)²) · G(s,t) over solid angle:
#   σ = λ_a² a_N²/(32π μ²v² m_A'⁴ s) · I,  I = G0/m_A'²·J2 − G1·J1 + G2·m_A'²·J0,
# G(s,t)=G0+G1 t+G2 t² (mn=MN). Heavy (m_A'²≫T) → 3λ_a²a_N²μ²/(π m_A'⁴) (G0→24mx²mN²);
# light (m_A'→0) → λ_a²a_N²μ²/(π m_A'⁴) — smooth crossover across the mediator-mass range.
function sigma_nucleon(m::AxialVector)
    a_N = _a_N(m.g)
    mu  = _mu_red(m.mx)
    s   = (m.mx + MN)^2
    x   = 4 * mu^2 * _VREL^2 / m.mphi^2
    J0, J1, J2 = _jmoments(x)
    G0 = m.mphi^4 * (2*MN^4 + 4*MN^2*(5*m.mx^2 - s) + 2*m.mx^4 - 4*m.mx^2*s + 2*s^2)
    G1 = -16*MN^2*m.mx^2*m.mphi^2 + m.mphi^4 * (2*s - 4*MN^2 - 4*m.mx^2)
    G2 = 8*MN^2*m.mx^2 + m.mphi^4
    I  = G0 * J2 / m.mphi^2 - G1 * J1 + G2 * m.mphi^2 * J0
    return m.gx^2 * a_N^2 / (32pi * mu^2 * _VREL^2 * m.mphi^4 * s) * I
end

# =====================================================================
# model_observables — single-point evaluator (relic, ΔNeff bbn/cmb, σ_χN).
# =====================================================================
"""
    model_observables(m; T_RH=0.01, production=:sm) -> NamedTuple

One-stop observable map for a model point.  Returns a NamedTuple:
  Oh2         relic abundance (authoritative loop solver: freeze_loop if
              mφ > 50·T_RH, else freeze_coupled)
  dNeff_bbn   NaN (neff_simple has no BBN epoch; the scan counts NaN as no
              contribution)
  dNeff_cmb   signed ΔNeff at the CMB epoch
  sigma_chiN  DM-nucleon cross section [cm²]
  Oh2_mediator         Ω h² of a mediator stable on the age of the Universe
              (freeze_coupled's info; in practice the axialVector A′ below
              min(2mχ, 3m_π)); NaN from freeze_loop
  mediator_overclosed  stable mediator with Oh2_mediator > 0.12 (not folded
              into Oh2 or eft.valid)
  solver      :freeze_loop or :freeze_coupled
  neff        the neff_simple info NamedTuple (`nothing` if it errored)
  errors      Dict{Symbol,String} of caught failures (empty if all succeeded)

Each observable is computed in its own try/catch; a failure degrades that field
to NaN rather than aborting the point.  Port of model_observables.m.
"""
function model_observables(m::CrossSectionModel; T_RH::Real=0.01, s_mix::Real=0.3,
                           production::Symbol=:sm)
    mphi = mass_mediator(m)
    errors = Dict{Symbol,String}()

    # --- relic abundance (loop solver, same mφ-vs-T_RH dispatch as compare_*) ---
    solver = mphi > 50 * T_RH ? :freeze_loop : :freeze_coupled
    Oh2, Oh2_mediator, mediator_overclosed = try
        if solver === :freeze_loop
            freeze_loop(m; T_RH=T_RH, production), NaN, false
        else
            O, fi = freeze_coupled(m; T_RH=T_RH, production, return_info=true)
            O, fi.Omega_phi_h2, fi.overclosed
        end
    catch e
        errors[:Oh2] = sprint(showerror, e)
        NaN, NaN, false
    end

    # --- ΔNeff (neff_simple): CMB epoch only ---
    dNeff_bbn, dNeff_cmb, neff = try
        d, ninfo = neff_simple(m; T_RH=T_RH)
        NaN, d, ninfo
    catch e
        errors[:neff] = sprint(showerror, e)
        NaN, NaN, nothing
    end

    # --- DM-nucleon cross section [cm²] ---
    sigma_chiN = try
        sigma_nucleon(m) * HBARC^2
    catch e
        errors[:sigma] = sprint(showerror, e)
        NaN
    end

    # --- EFT-validity flags (axialVector-only; trivially valid otherwise) ---
    eft = _eft_flags(m, s_mix, errors)

    return (Oh2 = Oh2, dNeff_bbn = dNeff_bbn, dNeff_cmb = dNeff_cmb,
            sigma_chiN = sigma_chiN, Oh2_mediator = Oh2_mediator,
            mediator_overclosed = mediator_overclosed, solver = solver, neff = neff, eft = eft,
            errors = errors)
end

# EFT-validity flags (mirrors model_observables.m's obs.eft and compare_limits_g.m).
# axialVector has two independent boundaries:
#   (1) unitarity ceiling — gx ≤ min(axial_gx_ceiling, √(π/2) mφ/mx) (gx-dependent).
#       axial_gx_ceiling keeps every meson BW width below its mass; it exists only
#       for mesons with j→χχ open (mj > 2mx), so it is Inf for mx > m_η′/2.  The
#       second term is perturbative unitarity of χχ→χχ through the longitudinal
#       A′, λ_a mx/m_A′ < √(π/2) (Kahlhoefer, Schmidt-Hoberg, Schwetz & Vogl,
#       JHEP 02 (2016) 016; same λ_a as Γ(A′→χχ̄) = λ_a² m_A′ β³/12π).  It bounds
#       the 1/m_A′⁴ of the J = 0 χχ̄→3π channel where the meson ceiling has gone.
#   (2) small-mixing breakdown — sin θ = max_j|m_jA′|/mφ < s_mix (gx-INdependent).
# PoleSP/PolePP rest on the same kind of expansion (the pseudoscalar mediator's
# mixing with π⁰/η/η′), so they get (2) with sin_mix = max_j|θ_j| from
# `sppp_mixing_angles`; there is no unitarity ceiling on that side.
# ScalarSS/ScalarPS go through the chiral loop factor 𝒟, not a mixing angle, so
# they really are trivially valid.
function _eft_flags(m::CrossSectionModel, s_mix, errors)
    valid = (valid = true, unitarity_ok = true, gx_ceiling = Inf,
             small_mixing_ok = true, sin_mix = NaN, s_mix = float(s_mix))
    if m isa PoleSP || m isa PolePP
        try
            _, theta_max, _ = sppp_mixing_angles(m.mphi, m.g)
            ok = theta_max < s_mix
            return (valid = ok, unitarity_ok = true, gx_ceiling = Inf,
                    small_mixing_ok = ok, sin_mix = theta_max, s_mix = float(s_mix))
        catch e
            errors[:eft] = sprint(showerror, e)
            return (valid = false, unitarity_ok = true, gx_ceiling = Inf,
                    small_mixing_ok = false, sin_mix = NaN, s_mix = float(s_mix))
        end
    end
    m isa AxialVector || return valid
    try
        gx_ceil = min(axial_gx_ceiling(m.mphi, m.mx, m.g), sqrt(pi / 2) * m.mphi / m.mx)
        _, mxmix, _ = axial_mixing_floor(m.mphi, m.mx, m.g, s_mix)
        sin_mix = mxmix / m.mphi
        unitarity_ok = m.gx <= gx_ceil             # NaN ceiling → false (invalid)
        small_mixing_ok = sin_mix < s_mix
        return (valid = unitarity_ok && small_mixing_ok, unitarity_ok = unitarity_ok,
                gx_ceiling = gx_ceil, small_mixing_ok = small_mixing_ok,
                sin_mix = sin_mix, s_mix = float(s_mix))
    catch e
        # _axial_couplings hard-errors when the mediator is too light (the
        # unitarity boundary itself): the point is EFT-invalid, not unknown.
        errors[:eft] = sprint(showerror, e)
        return (valid = false, unitarity_ok = false, gx_ceiling = Inf,
                small_mixing_ok = true, sin_mix = NaN, s_mix = float(s_mix))
    end
end

# Convenience: legacy MATLAB calling convention (interaction string).
model_observables(mx, mphi, gx, g, interaction::AbstractString; T_RH::Real=0.01) =
    model_observables(model_from_interaction(interaction, mx, mphi, gx, g); T_RH=T_RH)
