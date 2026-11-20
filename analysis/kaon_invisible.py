#!/usr/bin/env python3
r"""Vectorized rare-kaon invisible-decay constraint: K+ -> pi+ chi chi.

r = Gamma(K+ -> pi+ chi chi) / [K+ -> pi+ + (invisible) BR limit];  r > 1 is
excluded.  Open only for mx < (m_K+- - m_pi+-)/2 ~ 0.177 GeV (r = NaN when the
channel is closed or the interaction has no meson-decay implementation).

This is a NumPy-vectorized port of matlab/kaon_invisible_ratio.m, built to run over the many-million-row scan
tables in postprocess.py (per-point adaptive quadrature would take hours).  The
physics constants and the invisible-BR limit (1.14e-10) are kept byte-for-byte
identical to those modules, so the SMOOTH-region results agree to quadrature
tolerance.

Method
------
Each model reduces to a 1D integral over the pion energy E_pi in
[m_pi+-, E_pi^max].  We use a single fixed Gauss-Legendre rule (KNODES points),
mapped per point onto that point's interval, evaluated as one (npoints, KNODES)
array -- so millions of points integrate in a handful of vectorized passes.

Razor-thin, in-range resonances are unresolvable by any fixed grid, so each is
excised over a window [pole +/- delta] and added back analytically,
    int_{lo}^{hi} dE_pi / ((s - m^2)^2 + w^2)
        = (atan(u_lo/w) - atan(u_hi/w)) / (2 m_K w),   u = s(E_pi) - m^2,
freezing the smooth numerator at the pole.  delta is tied to the GL node spacing
so the outer panels stay well resolved for ANY width, narrow or broad.  Every
pole that can go on shell is cored: the s-channel phi BW (SxS/PxS), the pi0 AND
the mediator (SxP/PxP, the latter on shell for 2 mx < mphi < mK - mpi+), and the
pi0 in AxA.

MATLAB (matlab/kaon_invisible_ratio.m) is the reference.  Verified 2026-08-04
across all five interactions over mx in {0.02, 0.05, 0.1} x mphi in
{0.05, 0.2, 1.0}.  Agreement is ~1e-6 or better wherever the mediator pole is
well inside the s range or out of it.  It degrades only as that pole approaches
the UPPER endpoint E_pi^max (i.e. mphi -> 2 mx from above), where the BW pole and
the vanishing phase-space factor Y collide: ~1% at mphi = 1.25 * (2 mx), and up
to 27% at mphi = 2 mx exactly.  Moving 2.5% off that point in mphi restores 1e-6.
Both schemes are straining on a genuinely endpoint-singular configuration there;
neither is obviously right, and a scan point landing on it is measure-zero.

Faithfulness caveats (inherited / necessary):
  * The SxP/PxP amplitude's phi term and AxA's A' propagator carry no width in
    the source (bare 1/(s - m^2)); adaptive quadrature there returns an
    integrator-noise value near the pole (e.g. a spurious r ~ 1e9 even at
    gx = 1e-6, where the true rate is negligible).  We add a small physical width
    so the vectorized evaluator is finite and deterministic.  Smooth regions are
    unaffected; the r>1 mask matches the source at 99%+ of a validation grid, the
    exceptions being exactly these near-pole points where the regularized (and
    more physical) result is r < 1 rather than the source's pole-noise r >> 1.
  * AxA keeps only the real part of the amplitude (the source interference term
    is imaginary); AxA is flagged unvalidated upstream and is treated the same
    way here.
"""
import os

import numpy as np

from mesonandCMB import meson_pair_width, _sppp_theta, _relic_phi_widths

# --- self-contained constants (identical to kaon_invisible_ratio.m;
# verified 2026-08-04, MATLAB is the reference) ---------------------------------
_MK    = 493.68e-3
_MPI   = 139.57e-3           # charged pion
_MPI0  = 134.98e-3           # neutral pion
_META  = 547.86e-3
_METAP = 958e-3
_F     = 92.4e-3
_B0    = 0.13957 ** 2 / (0.0022 + 0.00469)   # mpi+^2/(mu+md) ~ 2.8273 GeV; was 2.6
_L4    = 1.87e-3
_L5    = 1.22e-3
_MU    = 2.2e-3
_MD    = 4.7e-3
_MS    = 95e-3
_E     = 0.3
_GF    = 1.17e-5
_VUS   = 0.222
_VUD   = 0.974
_G8    = 3.62
_G8S   = -0.38 * 3.62
_GPI   = 7.8e-9
_GETA  = 1.31e-6
_GETAP = 0.196e-3
_HBAR  = (197e-3 * 1e-13) / 3e10             # GeV*s (hbarc/c)
_KINV  = (1.0 / 1.238e-8) * 1.14e-10         # K+ invisible-decay rate limit [s^-1]

MX_THRESHOLD = (_MK - _MPI) / 2.0            # ~0.177 GeV; channel closed above
# NA62 acceptance (as kaon_invisible_ratio.m): every rate
# carries A(s), s = m_miss^2: R1 [0, 0.010] U R2 [0.026, 0.068] GeV^2 (K+ -> pi+ nu
# nu signal regions), edges smeared by the resolution 0.001 GeV^2.  The pi0 peak
# (s = 0.0182) is 8 sigma into the gap; a pole on an edge counts half.
NA62_SIG = 1e-3
# panel edges in s for the fixed GL rule: each region edge +- 4 sigma, plus the
# 7-sigma points 0.017 / 0.019 (gap) and 0.075 (above R2), where A < 1e-12.  The
# panels between those -- (0.017, 0.019), holding the pi0 pole, and s > 0.075 --
# are skipped, as na62_dead in kaon_invisible_ratio.m.
_NA62_PANELS = (0.006, 0.014, 0.017, 0.019, 0.022, 0.030, 0.064, 0.072, 0.075)
_NA62_DEAD = ((0.017, 0.019), (0.075, np.inf))


def na62_acceptance(s):
    from scipy.special import erfc
    Phi = lambda x: 0.5 * erfc(-x / np.sqrt(2.0))
    sg = NA62_SIG
    return (Phi((0.010 - s) / sg) - Phi((0.0 - s) / sg)
            + Phi((0.068 - s) / sg) - Phi((0.026 - s) / sg))

_KIND = {
    "heavyScalar": "SxS", "lightScalar": "SxS", "scalar": "SxS",
    "heavyScalarPS": "PxS", "lightScalarPS": "PxS", "scalarPS": "PxS",
    "heavyScalarSP": "SxP", "lightScalarSP": "SxP", "scalarSP": "SxP",
    "heavyScalarPP": "PxP", "lightScalarPP": "PxP", "scalarPP": "PxP",
    "axialVector": "AxA",
}

KNODES = int(os.environ.get("KAON_KNODES", 64))   # Gauss-Legendre order per E_pi panel (was 256 on one panel; 64 per NA62 panel agrees with MATLAB to 4e-6)
_CORE  = 4.0          # core half-window = _CORE * (interval / KNODES)
_x_gl, _w_gl = np.polynomial.legendre.leggauss(KNODES)   # nodes/weights on [-1,1]


def kind_of(interaction):
    """Model branch (SxS/PxS/SxP/PxP/AxA) for an interaction string, or None."""
    return _KIND.get(interaction)


def _epimax(mx):
    return (_MK ** 2 + _MPI ** 2 - 4.0 * mx ** 2) / (2.0 * _MK)


def _na62_panels(mx):
    """Per-point E_pi panels [(a, b), ...] covering [mpi, epimax], split at
    _NA62_PANELS so each acceptance edge sits mid-panel; empty ones have b == a.
    The dead panels (_NA62_DEAD, A < 1e-12) are left out."""
    E = lambda sv: (_MK ** 2 + _MPI ** 2 - sv) / (2.0 * _MK)
    emax = _epimax(mx)
    sb = (0.0,) + _NA62_PANELS + (np.inf,)
    out = []
    for s_lo, s_hi in zip(sb[:-1], sb[1:]):
        if (s_lo, s_hi) in _NA62_DEAD:
            continue
        a = np.full_like(mx, max(_MPI, E(s_hi)) if np.isfinite(s_hi) else _MPI)
        b = np.minimum(emax, E(s_lo))
        out.append((np.minimum(a, b), np.maximum(b, np.minimum(a, b))))
    return out


def _pole_epi(mphi):
    """E_pi of the on-shell mediator, s = mphi^2 (outside [mpi, epimax] when off shell)."""
    return (_MK ** 2 + _MPI ** 2 - mphi ** 2) / (2.0 * _MK)


def _default_panels(mx):
    return [(a, b, True) for a, b in _na62_panels(mx)]


def _accw(f):
    """f(Epi) x NA62 acceptance at s(Epi)."""
    return lambda Epi: f(Epi) * na62_acceptance(_MK ** 2 + _MPI ** 2 - 2.0 * _MK * Epi)


def _Yf(Epi, mx):
    num = (Epi ** 2 - _MPI ** 2) * (_MK ** 2 + _MPI ** 2 - 2.0 * Epi * _MK - 4.0 * mx ** 2)
    den = (_MK ** 2 + _MPI ** 2 - 2.0 * Epi * _MK)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.sqrt(np.clip(num / den, 0.0, None))


def _gl_panel(f, a, b):
    """Integrate f(Epi) over per-point interval [a, b] (arrays, shape (n,)) with
    the shared GL rule.  f maps an (n, K) node array to (n, K).  Zero-width
    panels (a == b) return 0 automatically."""
    a = a[:, None]
    b = b[:, None]
    Epi = 0.5 * (b - a) * _x_gl[None, :] + 0.5 * (b + a)
    return 0.5 * (b[:, 0] - a[:, 0]) * (f(Epi) * _w_gl[None, :]).sum(axis=1)


def _integrate_with_core(integrand, a, b, pole, wbw, num_core):
    """int_a^b integrand dEpi, excising a small window around `pole` (per point)
    and adding the analytic Breit-Wigner core there.

      integrand(Epi)  full (regularized) integrand, (n,K)->(n,K)
      a, b            interval ends, shape (n,)
      pole            resonance E_pi location, shape (n,); NaN/out-of-range -> no core
      wbw             Breit-Wigner half-width w = m*Gamma (s units), shape (n,)
      num_core(Epi)   smooth numerator N with local shape N/((s-m^2)^2 + w^2);
                      the core adds N(pole) * atan(2 mK delta / w)/(mK w)
    """
    n = a.shape[0]
    in_rng = np.isfinite(pole) & (pole > a) & (pole < b) & (wbw > 0.0)
    delta = _CORE * (b - a) / KNODES
    lo = np.where(in_rng, np.maximum(pole - delta, a), a)
    hi = np.where(in_rng, np.minimum(pole + delta, b), a)   # a==a -> panel B full below
    # panel A = [a, lo]; panel B = [hi, b].  For non-core points lo=a, hi=a ->
    # A empty, B = [a, b] = full interval.
    IA = _gl_panel(integrand, a, lo)
    IB = _gl_panel(integrand, hi, b)
    core = np.zeros(n)
    if in_rng.any():
        pr = np.where(in_rng, pole, a)                       # safe eval location
        Ncore = num_core(pr[:, None])[:, 0]
        with np.errstate(divide="ignore", invalid="ignore"):
            core_val = Ncore * np.arctan(2.0 * _MK * delta / wbw) / (_MK * wbw)
        core = np.where(in_rng, core_val, 0.0)
    return IA + IB + core


def _integrate_with_cores(integrand, a, b, specs):
    """int_a^b integrand dEpi with one or two Breit-Wigner cores excised and added
    back analytically.  SxP/PxP passes two (a pi0 pole and, for
    2 mx < mphi < mK - mpi+, an on-shell mediator pole); AxA passes one (pi0).

      specs   list of (pole, wbw, num_core), each as in _integrate_with_core

    Over an excised window [lo, hi] the exact Lorentzian integral is
      (atan(u_lo/w) - atan(u_hi/w)) / (2 mK w),   u = s(Epi) - m_j^2
                                                    = 2 mK (pole - Epi),
    which reduces to the symmetric atan(2 mK delta/w)/(mK w) of
    _integrate_with_core and stays exact when a window is clipped by a or b.
    Overlapping windows (only when mphi is within ~1 MeV of m_pi0) are merged so
    the continuum is excised exactly once; both cores are still added, since each
    is its own resonance term of |amp|^2 and the NWA drops the interference
    either way.  Same construction as the notch loop in kaon_invisible_ratio.m,
    so the two ports resolve the same poles.
    """
    n = a.shape[0]
    delta = _CORE * (b - a) / KNODES
    los, his, cores = [], [], []
    for pole, wbw, num_core in specs:
        in_rng = np.isfinite(pole) & (pole > a) & (pole < b) & (wbw > 0.0)
        # out-of-range -> degenerate window parked at b, so it sorts last and
        # contributes only empty panels
        lo = np.where(in_rng, np.maximum(pole - delta, a), b)
        hi = np.where(in_rng, np.minimum(pole + delta, b), b)
        safe_pole = np.where(in_rng, pole, a)
        with np.errstate(divide="ignore", invalid="ignore"):
            u_lo = 2.0 * _MK * (pole - lo)
            u_hi = 2.0 * _MK * (pole - hi)
            cv = (num_core(safe_pole[:, None])[:, 0]
                  * (np.arctan(u_lo / wbw) - np.arctan(u_hi / wbw))
                  / (2.0 * _MK * wbw))
            # Subtract-and-add: the analytic core is only the j-th resonance
            # term, N_j(pole) * int L_j.  Replacing the WHOLE integrand over the
            # window by it would discard the other poles' continuum there, which
            # is what left this branch several % off MATLAB/Julia near the range
            # ends.  Integrate the remainder integrand - N_j(Epi) L_j(Epi) over
            # the window as well; the sharp part is gone, so the GL rule handles
            # what is left.
            def _rem(Epi, _nc=num_core, _p=safe_pole[:, None], _w=wbw[:, None]):
                u = 2.0 * _MK * (_p - Epi)
                return integrand(Epi) - _nc(Epi) / (u * u + _w * _w)
            rem = _gl_panel(_rem, np.where(in_rng, lo, a), np.where(in_rng, hi, a))
        los.append(lo)
        his.append(hi)
        cores.append(np.where(in_rng, cv + rem, 0.0))

    if len(specs) == 1:
        return (_gl_panel(integrand, a, los[0]) + _gl_panel(integrand, his[0], b)
                + cores[0])

    (lo1, lo2), (hi1, hi2) = los, his
    swap = lo2 < lo1
    lo1, lo2 = np.where(swap, lo2, lo1), np.where(swap, lo1, lo2)
    hi1, hi2 = np.where(swap, hi2, hi1), np.where(swap, hi1, hi2)
    overlap = lo2 <= hi1
    hi1 = np.where(overlap, np.maximum(hi1, hi2), hi1)       # merge into window 1
    lo2 = np.where(overlap, hi1, lo2)                        # window 2 -> degenerate
    hi2 = np.where(overlap, hi1, hi2)

    return (_gl_panel(integrand, a, lo1) + _gl_panel(integrand, hi1, lo2)
            + _gl_panel(integrand, hi2, b) + cores[0] + cores[1])


def _ratio_sxs_pxs(mx, mphi, gx, Cu, Cd, Cs, kind, panels=None):
    """SxS / PxS: scalar phi-quark vertex, s-channel phi Breit-Wigner.

    panels: [(a, b, cores), ...] E_pi panels to integrate (default: every live
    NA62 panel with the pole core on); the ratio is linear in their sum, which is
    how _ratio_chunk splits the work (see there)."""
    if panels is None:
        panels = _default_panels(mx)
    delta = -_GF * _F ** 2 * _VUS * _VUD * _G8 / np.sqrt(2.0)
    Dnew = 0.5 * (_B0 * (Cu + Cd) / _MPI ** 2 + _B0 * (Cu + Cs) / _MK ** 2) ** 2
    # gamma gamma + phi -> pi pi, K K, eta eta, pi0 eta (open above 2 m_pi0);
    # chi chibar added below
    Gbase = ((Dnew * _E ** 4 * mphi ** 3) / (73728.0 * np.pi ** 5)
             + meson_pair_width(mphi, Cu, Cd, Cs))
    open_chi = mphi ** 2 > 4.0 * mx ** 2
    with np.errstate(invalid="ignore"):
        if kind == "SxS":
            gterm = gx ** 2 * np.power(np.abs(mphi ** 2 - 4.0 * mx ** 2), 1.5) / (8.0 * np.pi * mphi ** 2)
        else:
            gterm = gx ** 2 * np.sqrt(np.abs(mphi ** 2 - 4.0 * mx ** 2)) / (8.0 * np.pi)
    Gphi = Gbase + np.where(open_chi, gterm, 0.0)
    wbw = mphi * Gphi                                        # BW half-width (s units)

    def _num(Epi):                                          # numerator N(Epi)*Yf
        base = (_MK ** 2 + _MPI ** 2 - 2.0 * _MK * Epi)
        if kind == "SxS":
            base = base - 4.0 * mx[:, None] ** 2
        return base * _Yf(Epi, mx[:, None])

    def _integrand(Epi):
        s = _MK ** 2 + _MPI ** 2 - 2.0 * _MK * Epi
        den = (s - mphi[:, None] ** 2) ** 2 + (mphi[:, None] ** 2) * (Gphi[:, None] ** 2)
        return _num(Epi) / den

    Epi_res = _pole_epi(mphi)
    Fint = sum(_integrate_with_core(_accw(_integrand), a, b, Epi_res, wbw, _accw(_num))
               if cores else _gl_panel(_accw(_integrand), a, b)
               for a, b, cores in panels)

    gpiK = delta * _B0 * ((Cs + Cu) * _MPI - (Cu + Cd) * _MK) / (_MK - _MPI)
    Gamma = gx ** 2 * gpiK ** 2 * Fint / (32.0 * np.pi ** 3 * _MK) / _HBAR
    return Gamma / _KINV


def _ratio_sxp_pxp(mx, mphi, gx, Cu, Cd, Cs, kind, panels=None):
    """SxP / PxP: pseudoscalar phi-quark vertex, pi0/eta/eta' mixing poles.

    theta_j from _sppp_theta (sppp_mixing.m: NLO, with eta'), g_{j chi} = +gx theta_j
    for every j, K -> pi j vertices from eq. Klagrangian, mediator total width
    (chi chibar + gamma gamma).  Until 2026-11-17: LO B0 F mixings, no eta', a
    spurious minus on the pi0 chi coupling, and a gamma-gamma-only phi width.  As
    kaon_invisible_ratio.m.  panels: as _ratio_sxs_pxs."""
    if panels is None:
        panels = _default_panels(mx)
    mphi2 = mphi ** 2
    th = _sppp_theta(mphi, Cu, Cd, Cs)
    gpichi = gx * th[0]
    getachi = gx * th[1]
    getapchi = gx * th[2]
    gphichi = gx
    prefac = (_GF * _F * _G8 * _VUS * _VUD) ** 2 / (64.0 * np.pi ** 3 * _MK)
    ps_sub = 4.0 * mx ** 2 if kind == "SxP" else np.zeros_like(mx)

    vis, xx = _relic_phi_widths(kind, mx, mphi, gx, Cu, Cd, Cs)
    Gphi_reg = vis + xx
    wphi = (mphi * Gphi_reg)[:, None]
    c_th = th[:, :, None]
    c_getapchi = getapchi[:, None]
    c_gpichi = gpichi[:, None]
    c_gphichi = gphichi[:, None]
    c_getachi = getachi[:, None]
    c_mphi2 = mphi2[:, None]
    c_ps = ps_sub[:, None]
    c_mx = mx[:, None]

    def _gpipik(Epi):
        return 2.0 * _MK * Epi - _MK ** 2

    def _getapik(Epi):
        return (_MK ** 2 - 6.0 * _MK * Epi + 2.0 * _MPI ** 2) / np.sqrt(3.0)

    def _getappik(Epi):
        return np.sqrt(6.0) / 3.0 * (2.0 + 3.0 * _G8S / _G8) * (_MK ** 2 - _MPI ** 2) + 0.0 * Epi

    def _gphipik(Epi):
        return -c_th[0] * _gpipik(Epi) - c_th[1] * _getapik(Epi) - c_th[2] * _getappik(Epi)

    def _PSf(Epi):
        return _MK ** 2 + _MPI ** 2 - 2.0 * _MK * Epi - c_ps

    def _amp2(Epi):
        s = _MK ** 2 + _MPI ** 2 - 2.0 * _MK * Epi
        amp = (c_gpichi * _gpipik(Epi)) / (s - _MPI0 ** 2 + 1j * _MPI0 * _GPI) \
            + (c_gphichi * _gphipik(Epi)) / (s - c_mphi2 + 1j * wphi) \
            + (c_getachi * _getapik(Epi)) / (s - _META ** 2 + 1j * _META * _GETA) \
            + (c_getapchi * _getappik(Epi)) / (s - _METAP ** 2 + 1j * _METAP * _GETAP)
        return np.abs(amp) ** 2

    def _integrand(Epi):
        return prefac * _PSf(Epi) * _Yf(Epi, c_mx) * _amp2(Epi)

    def _num_core(Epi):        # diagonal pi0 numerator: prefac*PSf*Yf*|gpichi*gpipik|^2
        nj = c_gpichi * _gpipik(Epi)
        return prefac * _PSf(Epi) * _Yf(Epi, c_mx) * np.abs(nj) ** 2

    def _num_core_phi(Epi):    # diagonal mediator numerator
        nj = c_gphichi * _gphipik(Epi)
        return prefac * _PSf(Epi) * _Yf(Epi, c_mx) * np.abs(nj) ** 2

    # Two poles can be on shell in range (the eta pole s=meta^2 > s_max never is):
    # the pi0, and the MEDIATOR for 2 mx < mphi < mK - mpi+ ~ 0.354 GeV.  The
    # mediator one is regularized by wphi above but is far narrower than the node
    # spacing (Gamma_phi/mphi ~ 1e-5), so leaving it to the fixed GL rule loses
    # most of the resonance -- up to 22% low at mphi = 0.35 against MATLAB/Julia.
    # The pi0 pole (s = mpi0^2) lies in the skipped NA62 gap panel, so only the
    # mediator is cored; _num_core is kept for reference.
    Epi_phi = _pole_epi(mphi)
    Icont = sum(_integrate_with_cores(_accw(_integrand), a, b, [
        (Epi_phi, mphi * Gphi_reg, _accw(_num_core_phi))])
        if cores else _gl_panel(_accw(_integrand), a, b)
        for a, b, cores in panels)
    return Icont / _HBAR / _KINV


def _ratio_axa(mx, mphi, gx, Cu, Cd, Cs, panels=None):
    """AxA: axial-vector A', pi0/eta/eta' A'-mixing.  Draft eq. (kpichichiA),
    rechecked 2026-11-17: signed lambda_jpiK (common i dropped), the A' propagator
    in the interference term, s - mA^2 = mK^2 + mpi^2 - 2 Epi mK in the |A'|^2
    term, lambda_A'K^2 mA^4 mK^4, Gamma_A' = chi chibar width (no width before).
    The E_chi integral is analytic (quadratic); the pi0 pole and the on-shell A'
    (2 mx < mA < mK - mpi) are cored.  As kaon_invisible_ratio.m.
    panels: as _ratio_sxs_pxs."""
    if panels is None:
        panels = _default_panels(mx)
    mmed = mphi
    mdm = mx
    mmed2 = mmed ** 2

    mpi0A = _F * ((Cu - Cd) * (1.0 + (16.0 * _B0 * _L4) / _F ** 2 * (_MU + _MD + _MS))
                  + (16.0 * _B0 * _L5 * (Cu * _MU - Cd * _MD)) / _F ** 2)
    metaA = (np.sqrt(3.0) * _F) / 3.0 * ((Cu + Cd - 2.0 * Cs) * (1.0 + (16.0 * _B0 * _L4) / _F ** 2 * (_MU + _MD + _MS))
                                         + 1.0 / _F ** 2 * 16.0 * _B0 * _L5 * (Cu * _MU + Cd * _MD - 2.0 * Cs * _MS))
    meta0A = (np.sqrt(6.0) * _F) / 3.0 * ((Cu + Cd + Cs) * (1.0 + (16.0 * _B0 * _L4) / _F ** 2 * (_MU + _MD + _MS))
                                          + 1.0 / _F ** 2 * 16.0 * _B0 * _L5 * (Cu * _MU + Cd * _MD + Cs * _MS))
    lak = _G8 * (Cs + Cu) + _G8S * (Cu + Cd + Cs)
    lap = _G8 * (Cd + Cu) + _G8S * (Cu + Cd + Cs)
    GA = np.where(mmed > 2.0 * mdm,
                  gx ** 2 * np.clip(mmed2 - 4.0 * mdm ** 2, 0.0, None) ** 1.5 / (12.0 * np.pi * mmed2), 0.0)

    Cu1, Cd1, Cs1 = Cu[:, None], Cd[:, None], Cs[:, None]
    mmed21, mdm1, gx1, GA1 = mmed2[:, None], mdm[:, None], gx[:, None], GA[:, None]
    lak1, lap1 = lak[:, None], lap[:, None]
    lapk, lamk = lak1 + lap1, lak1 - lap1
    lj = [(gx * mpi0A / mmed2)[:, None], (gx * metaA / mmed2)[:, None], (gx * meta0A / mmed2)[:, None]]
    mA_j = [mpi0A[:, None], metaA[:, None], meta0A[:, None]]
    mj = (_MPI0, _META, _METAP)
    Gj = (_GPI, _GETA, _GETAP)

    def _pikbulk(Epi):
        return (_G8 * ((Cd1 + Cu1) * (Epi * _MK - _MPI ** 2) - (Cs1 + Cu1) * (Epi * _MK - _MK ** 2))
                + _G8S * (_MK ** 2 - _MPI ** 2) * (Cu1 + Cd1 + Cs1))

    def _lpk(Epi):                      # lambda_jpiK / i, real and signed
        bulk = _pikbulk(Epi)
        return (_G8 / _F ** 3 * (_MK ** 2 - 2.0 * Epi * _MK) + 2.0 * mA_j[0] / (_F ** 2 * mmed21) * bulk,
                2.0 * np.sqrt(3.0) * _G8 / (3.0 * _F ** 3) * (_MK ** 2 - _MPI ** 2) + 2.0 * mA_j[1] / (_F ** 2 * mmed21) * bulk,
                np.sqrt(6.0) / (3.0 * _F ** 3) * (2.0 * _G8 + 3.0 * _G8S) * (_MK ** 2 - _MPI ** 2)
                + 2.0 * mA_j[2] / (_F ** 2 * mmed21) * bulk)

    def _s(Epi):
        return _MK ** 2 + _MPI ** 2 - 2.0 * Epi * _MK

    def _DA(Epi):
        return (_s(Epi) - mmed21) ** 2 + mmed21 * GA1 ** 2

    pref = (_GF ** 2 * _F ** 8) / (128.0 * np.pi ** 3 * _MK) * _VUS ** 2 * _VUD ** 2 / _HBAR
    c3 = -8.0 * gx1 ** 2 / (_F ** 4 * mmed21 ** 2)

    def _moments(Epi):
        with np.errstate(invalid="ignore", divide="ignore"):
            d = 0.5 * np.sqrt(np.clip((Epi ** 2 - _MPI ** 2) * (_s(Epi) - 4.0 * mdm1 ** 2) / _s(Epi), 0.0, None))
        c = 0.5 * (_MK - Epi)
        return 2.0 * d, 2.0 * c * d, 2.0 * d * (3.0 * c ** 2 + d ** 2) / 3.0

    def _N3(Epi):                       # |A'|^2 term times DA (smooth across the pole)
        I0, I1, I2 = _moments(Epi)
        Q2 = 4.0 * mmed21 ** 2 * _MK ** 2 * lapk ** 2
        Q1 = 4.0 * mmed21 ** 2 * _MK ** 2 * (Epi - _MK) * lapk ** 2
        Q0 = (8.0 * Epi ** 3 * _MK ** 3 * mdm1 ** 2 * lamk ** 2
              + 4.0 * Epi ** 2 * _MK ** 2 * mdm1 ** 2 * lamk * (2.0 * mmed21 * lamk + _MK ** 2 * (lap1 - 5.0 * lak1) - _MPI ** 2 * (lak1 - 5.0 * lap1))
              - 2.0 * Epi * _MK * (mmed21 ** 2 * (_MK ** 2 * lapk ** 2 + 4.0 * lak1 * lap1 * mdm1 ** 2)
                                   + 8.0 * mmed21 * mdm1 ** 2 * lamk * (lak1 * _MK ** 2 - lap1 * _MPI ** 2)
                                   - 4.0 * mdm1 ** 2 * (lak1 * _MK ** 4 * (2.0 * lak1 - lap1) + _MK ** 2 * _MPI ** 2 * (lak1 ** 2 - 4.0 * lak1 * lap1 + lap1 ** 2) + lap1 * _MPI ** 4 * (2.0 * lap1 - lak1)))
              + mmed21 ** 2 * (_MK ** 4 * lapk ** 2 + _MK ** 2 * (_MPI ** 2 * lapk ** 2 - 4.0 * lak1 ** 2 * mdm1 ** 2) - 4.0 * lap1 ** 2 * _MPI ** 2 * mdm1 ** 2)
              + 8.0 * mmed21 * mdm1 ** 2 * (lak1 * _MK ** 2 - lap1 * _MPI ** 2) ** 2
              - 4.0 * mdm1 ** 2 * (_MK ** 2 + _MPI ** 2) * (lak1 * _MK ** 2 - lap1 * _MPI ** 2) ** 2)
        return pref * c3 * (Q0 * I0 + Q1 * I1 + Q2 * I2)

    def _outer(Epi):
        s = _s(Epi)
        L = _lpk(Epi)
        A = 0.0
        S2 = 0.0
        for k in range(3):
            num = lj[k] * L[k]
            A = A + num / (s - mj[k] ** 2 + 1j * mj[k] * Gj[k])
            S2 = S2 + num * ((mmed21 - s) * (mj[k] ** 2 - s) + GA1 * Gj[k] * mmed21 ** 0.5 * mj[k]) \
                / ((mj[k] ** 2 - s) ** 2 + Gj[k] ** 2 * mj[k] ** 2)
        T12 = (8.0 * mdm1 ** 2 * s * np.abs(A) ** 2
               - 32.0 * gx1 * mdm1 ** 2 * (mmed21 - s) * (_MK * Epi * lamk - _MK ** 2 * lak1 + _MPI ** 2 * lap1) * S2
               / (_F ** 2 * mmed21 * _DA(Epi)))
        I0 = _moments(Epi)[0]
        # zero phase space (I0 = 0) with the A' propagator also vanishing there
        # (mA = 2 mx, Gamma_A' = 0: pole on the endpoint) is 0/0 -> take 0.
        with np.errstate(invalid="ignore", divide="ignore"):
            v = pref * T12 * I0 + _N3(Epi) / _DA(Epi)
        return np.where(I0 > 0.0, v, 0.0)

    def _num_core_pi0(Epi):
        """Diagonal pi0 term of |A|^2 (the only piece resonant at s = mpi0^2)."""
        N = lj[0] * _lpk(Epi)[0]
        return pref * _moments(Epi)[0] * 8.0 * mdm1 ** 2 * _s(Epi) * N ** 2

    # pi0 pole: in the skipped gap panel (see _ratio_sxp_pxp); A' pole cored.
    Epi_A = _pole_epi(mmed)
    Gamma = sum(_integrate_with_cores(_accw(_outer), a, b, [(Epi_A, mmed * GA, _accw(_N3))])
                if cores else _gl_panel(_accw(_outer), a, b)
                for a, b, cores in panels)
    return Gamma / _KINV


CHUNK = int(os.environ.get("KAON_CHUNK", 50_000))   # points per pass; AxA builds
# ~40 live (CHUNK x KNODES) temporaries, so peak ~ CHUNK*KNODES*8*40 B
# (~4 GB at the default).  Lower KAON_CHUNK on a memory-tight box.


def _ratio_chunk(mx, mphi, gx, gu, gd, gs, kind):
    """The ratio is linear in the sum over E_pi panels, so integrate panel by
    panel on just the points it is non-empty for, and run the pole-core machinery
    (3 GL passes) only on those whose mediator pole lies in it; the rest take one
    GL pass.  (Until 2026-11-18 every point ran every panel with both the pi0 and
    the mediator core, ~35 passes: 10x slower than before the NA62 acceptance.)"""
    if kind in ("SxS", "PxS"):
        f = lambda j, p: _ratio_sxs_pxs(mx[j], mphi[j], gx[j], gu[j], gd[j], gs[j], kind, p)
    elif kind in ("SxP", "PxP"):
        f = lambda j, p: _ratio_sxp_pxp(mx[j], mphi[j], gx[j], gu[j], gd[j], gs[j], kind, p)
    else:
        f = lambda j, p: _ratio_axa(mx[j], mphi[j], gx[j], gu[j], gd[j], gs[j], p)
    r = np.zeros(mx.shape)
    pole = _pole_epi(mphi)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for a, b in _na62_panels(mx):
            live = b > a
            inp = live & (pole > a) & (pole < b)
            for sel, cores in ((live & ~inp, False), (inp, True)):
                j = np.flatnonzero(sel)
                if j.size:
                    r[j] += f(j, [(a[j], b[j], cores)])
    return r


def _ratio_job(args):
    return _ratio_chunk(*args)


def _nproc():
    """Worker processes for kaon_ratio: KAON_NPROC, else the SLURM allocation."""
    v = os.environ.get("KAON_NPROC") or os.environ.get("SLURM_CPUS_PER_TASK") or "1"
    try:
        return max(1, int(v))
    except ValueError:
        return 1


def kaon_ratio(mx, mphi, gx, gu, gd, gs, interaction):
    """r = Gamma(K+->pi+ chi chi)/limit for each point (array).  r>1 excluded;
    NaN where the channel is closed (mx >= 0.177 GeV) or the model has no
    meson-decay implementation.  All inputs broadcast to a common shape.

    Evaluated only on open-channel points, in chunks of CHUNK to bound the
    (npoints x KNODES) working arrays.  With more than one process (_nproc:
    KAON_NPROC, else SLURM_CPUS_PER_TASK) the chunks go to a fork pool, each
    worker taking chunks of KAON_CHUNK_PAR (default 10 000, ~0.8 GB) points."""
    kind = kind_of(interaction)
    mx, mphi, gx, gu, gd, gs = np.broadcast_arrays(
        np.asarray(mx, float), np.asarray(mphi, float), np.asarray(gx, float),
        np.asarray(gu, float), np.asarray(gd, float), np.asarray(gs, float))
    shape = mx.shape
    if kind is None:
        return np.full(shape, np.nan)
    mx = mx.ravel(); mphi = mphi.ravel(); gx = gx.ravel()
    gu = gu.ravel(); gd = gd.ravel(); gs = gs.ravel()
    r = np.full(mx.shape, np.nan)
    idx = np.flatnonzero(mx < MX_THRESHOLD)     # open-channel points only
    nproc = _nproc()
    if nproc > 1 and idx.size > 2 * int(os.environ.get("KAON_CHUNK_PAR", 10_000)):
        import multiprocessing as mp
        step = int(os.environ.get("KAON_CHUNK_PAR", 10_000))
        parts = [idx[s:s + step] for s in range(0, idx.size, step)]
        jobs = ((mx[j], mphi[j], gx[j], gu[j], gd[j], gs[j], kind) for j in parts)
        with mp.get_context("fork").Pool(nproc) as pool:
            for j, rj in zip(parts, pool.imap(_ratio_job, jobs, chunksize=1)):
                r[j] = rj
        return r.reshape(shape)
    for s in range(0, idx.size, CHUNK):
        j = idx[s:s + CHUNK]
        r[j] = _ratio_chunk(mx[j], mphi[j], gx[j], gu[j], gd[j], gs[j], kind)
    return r.reshape(shape)
