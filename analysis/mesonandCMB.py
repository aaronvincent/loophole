#!/usr/bin/env python3
r"""Vectorized meson-decay and indirect-detection limits beyond K+ -> pi+ chi chi.

Three exclusion ratios (signal / bound; r > 1 excluded), each a NumPy-vectorized
port of matlab/mesonandCMB_ratio.m:

  kphi_ratio      Gamma(K+ -> pi+ phi) * P_inv / [K+->pi+ + invisible limit].
                  On-shell mediator production; open only for mphi < m_K+- -
                  m_pi+- ~ 0.354 GeV (NaN above).  P_inv is the chance the phi
                  looks invisible (kphi_invisible_fraction); gx enters only there.
  cmb_ratio       <sigma v>_low-v / Slatyer CMB energy-injection bound.  ~gx^2;
                  NaN for SxS/SxP (velocity-suppressed -> no CMB bound).
  integral_ratio  <sigma v>_low-v / INTEGRAL telescope bound.  ~gx^2; all models.
  meson_invisible_ratio
                  Gamma(j -> chi chibar) / Gamma(j -> invisible) limit for
                  j = pi0, eta, eta' (or the largest of the three).  The meson
                  couples to chi by mixing with the mediator: SxP/PxP/AxA only,
                  NaN for SxS/PxS.  ~gx^2.

Unlike kaon_invisible.py these are closed forms (no quadrature), so the port is a
direct vectorized expression -- no fixed-grid rule, no chunking.  Constants and
the K+ invisible-BR limit (1.14e-10) are byte-for-byte identical to the MATLAB /
Julia ports, and the <sigma v> is the reference's analytic low-velocity
(s = 4 m_chi^2) expression, NOT a v->0 limit of the relic svavg.
"""
import os

import numpy as np

# --- self-contained constants (identical to mesonandCMB_ratio.m;
# verified 2026-08-04, MATLAB is the reference) --------------------------------
_E     = 0.3
_ALPHA = 1.0 / 137.0
_B0    = 0.13957 ** 2 / (0.0022 + 0.00469)   # mpi+^2/(mu+md) ~ 2.8273 GeV; was 2.6
_F     = 92.4e-3
_MPI0  = 134.98e-3
_META  = 547.86e-3
_METAP = 958e-3                              # eta'
_MPI   = 139.57e-3                           # charged pion
_MK    = 493.68e-3                           # charged kaon
_MP    = 0.938
_GF    = 1.17e-5
_VUS   = 0.222
_VUD   = 0.974
_G8    = 3.62
_G8S   = -0.38 * 3.62
_L4    = 1.87e-3
_L5    = 1.22e-3
_L6    = 1.46e-3
_L7    = -0.39e-3
_MU    = 2.2e-3
_MD    = 4.7e-3
_MS    = 95e-3
_HBARC = 197e-3 * 1e-13                      # GeV*cm
_CC    = 3e10                                # cm/s
_HBAR  = _HBARC / _CC                        # GeV*s
_GAMMAPI0 = _HBAR / 8.52e-17
_GETA     = 1.31e-6
_GETA0    = 0.188e-3
_KINV     = (1.0 / 1.238e-8) * 1.14e-10      # K+ invisible-decay rate limit [s^-1]

KPHI_THRESHOLD = _MK - _MPI                  # ~0.354 GeV; K->pi phi closed above

_KIND = {
    "heavyScalar": "SxS", "lightScalar": "SxS", "scalar": "SxS",
    "heavyScalarPS": "PxS", "lightScalarPS": "PxS", "scalarPS": "PxS",
    "heavyScalarSP": "SxP", "lightScalarSP": "SxP", "scalarSP": "SxP",
    "heavyScalarPP": "PxP", "lightScalarPP": "PxP", "scalarPP": "PxP",
    "axialVector": "AxA",
}

# --- inlined digitized limit curves (log10 mx/GeV vs log10 <sigma v>/cm^3 s^-1)
_SLATYER_X = np.array([
    -5.299492385786800, -5.0152284263959400, -4.770123277737490,
    -4.501812907904280, -4.171138506163890, -3.908629441624370,
    -3.537345902828140, -3.09934735315446, -2.973168963016680,
    -2.742567077592460, -2.4829586656997800, -2.265409717186370,
    -2.029006526468460, -1.789702683103700, -1.361856417693980,
    -0.9905728788977540, -0.45250181290790600, -0.014503263234229500,
    0.2871646120377070, 0.6932559825960830])
_SLATYER_Y = np.array([
    -32.62642947903430, -32.27573062261750, -31.99618805590850,
    -31.726810673443500, -31.38627700127070, -31.058449809402800,
    -30.49682337992380, -29.73189326556540, -29.513341804320200,
    -29.337992376111800, -29.363405336721700, -29.320203303684900,
    -29.19059720457430, -29.007623888183000, -28.606099110546400,
    -28.21982210927570, -27.57941550190600, -26.97712833545110,
    -26.601016518424400, -26.260482846251600])
_INTEGRAL_X = np.array([
    -1.5702280912364900, -1.312124849939980, -1.3085234093637500,
    -1.0036014405762300, -0.9987995198079230, -0.701080432172869,
    -0.6974789915966390, -0.22448979591836700, -0.21848739495798300,
    0.2557022809123650]) - 3.0
_INTEGRAL_Y = np.array([
    -32.47761194029850, -31.970149253731300, -32.471641791044800,
    -31.86268656716420, -31.907462686567200, -31.31641791044780,
    -31.01492537313430, -30.074626865671600, -30.555223880597000,
    -29.608955223880600])


def kind_of(interaction):
    """Model branch (SxS/PxS/SxP/PxP/AxA) for an interaction string, or None."""
    return _KIND.get(interaction)


def _curve(mx, xp, yp):
    # log10 <sigma v> bound, NaN outside the digitized data range.
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.log10(np.asarray(mx, dtype=float))
    return 10.0 ** np.interp(lg, xp, yp, left=np.nan, right=np.nan)


def _pp_couplings(Cu, Cd, Cs):
    # phi-meson^2 couplings for the pseudoscalar phi-quark vertex (L6/L7-corrected).
    mpi0phi2 = ((Cu - Cd) * (_B0 * _F + 32 * _B0**2 * _L6 / _F * (_MU + _MD + _MS))
                + 32 * _B0**2 * _L7 / _F * (Cu + Cd + Cs) * (_MU - _MD))
    # eta = 1/sqrt(3) x the pi0 bracket (chiral trace; see sppp_couplings.m).
    metaphi2 = ((Cu + Cd - 2 * Cs) * (_B0 * _F + 32 * _B0**2 * _L6 / _F * (_MU + _MD + _MS))
                + 32 * _B0**2 * _L7 / _F * (Cu + Cd + Cs) * (_MU + _MD - 2 * _MS)) / np.sqrt(3)
    return mpi0phi2, metaphi2


# phi -> two-meson pairs of the scalar phi-SM vertex (SxS/PxS): port of
# matlab/meson_pairs.m (masses, B0 and the coefficients mu as there).
#   (mA, mB, S, mu(Cu, Cd, Cs)/B0)
_B0MP = 0.13957 ** 2 / (0.0022 + 0.00469)
_MESON_PAIRS = (
    (0.13957,  0.13957,  1.0, lambda u, d, s: u + d),                 # pi+ pi-
    (0.13498,  0.13498,  0.5, lambda u, d, s: u + d),                 # pi0 pi0
    (0.493677, 0.493677, 1.0, lambda u, d, s: u + s),                 # K+ K-
    (0.497611, 0.497611, 1.0, lambda u, d, s: d + s),                 # K0 K0bar
    (0.547862, 0.547862, 0.5, lambda u, d, s: (u + d + 4 * s) / 3),   # eta eta
    (0.13498,  0.547862, 1.0, lambda u, d, s: (u - d) / np.sqrt(3)),  # pi0 eta
)


def meson_pair_width(mphi, Cu, Cd, Cs):
    """Gamma(phi -> pi pi, K K, eta eta, pi0 eta) [GeV], scalar phi-SM vertex:
    sum of S mu^2 sqrt(lambda(mphi^2, mA^2, mB^2)) / (16 pi mphi^3), each gated on
    its threshold (so mphi = 0 gives 0).  As phi_partial_widths.m / widths.jl."""
    mphi = np.asarray(mphi, dtype=float)
    G = np.zeros(np.broadcast(mphi, Cu, Cd, Cs).shape)
    with np.errstate(invalid="ignore", divide="ignore"):
        for mA, mB, S, c in _MESON_PAIRS:
            lam = (mphi**2 - (mA + mB)**2) * (mphi**2 - (mA - mB)**2)
            w = S * (_B0MP * c(Cu, Cd, Cs))**2 * np.sqrt(np.clip(lam, 0.0, None)) / (16 * np.pi * mphi**3)
            G = G + np.where(mphi > mA + mB, w, 0.0)
    return G


def _phi_gg_xx_width(kind, mx, mphi, gx, Cu, Cd, Cs):
    """Gamma(phi -> gamma gamma) + Gamma(phi -> chi chibar), scalar phi-SM vertex:
    SxS beta^3, PxS beta^1 chi chibar width."""
    Dnew = 0.5*(_B0 * (Cu + Cd) / _MPI**2 + _B0 * (Cu + Cs) / _MK**2)**2
    Gbase = (Dnew * _E**4 * mphi**3) / (73728 * np.pi**5)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.clip(mphi**2 - 4 * mx**2, 0.0, None)
        if kind == "SxS":
            gxx = gx**2 * d**1.5 / (8 * np.pi * mphi**2)
        else:
            gxx = gx**2 * np.sqrt(d) / (8 * np.pi)
    return Gbase + np.where(mphi**2 > 4 * mx**2, gxx, 0.0)


# Mediator widths as the relic code has them (phi_partial_widths.m /
# sppp_mixing.m / widths.jl), with ITS constants -- not the kphi/CMB set above.
_RALPHA, _RB0, _RF = 1.0 / 137, 0.13957**2 / (0.0022 + 0.00469), 0.0924
_RMPI0, _RMPIPM, _RMK, _RMETA, _RMETAP = 0.13498, 0.13957, 0.493677, 0.547862, 0.95778
_RMU, _RMD, _RMS = 0.0022, 0.00469, 0.095


def _relic_phi_widths(kind, mx, mphi, gx, Cu, Cd, Cs):
    """(Gamma_visible, Gamma_chichibar) [GeV] of phi: visible = gamma gamma + meson
    pairs.  SxS/PxS/SxP/PxP; as phi_partial_widths.m."""
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.clip(mphi**2 - 4 * mx**2, 0.0, None)
        if kind in ("SxS", "PxS"):
            D = 0.5 * (_RB0 * (Cu + Cd) / _RMPIPM**2 + _RB0 * (Cu + Cs) / _RMK**2)**2
            gg = _RALPHA**2 * D * mphi**3 / (4608 * np.pi**3)
            vis = gg + meson_pair_width(mphi, Cu, Cd, Cs)
        else:
            sg = Cu + Cd + Cs
            K = _RB0 * _RF + 32 * _RB0**2 * _L6 / _RF * (_RMU + _RMD + _RMS)
            m2pi = (Cu - Cd) * K + 32 * _RB0**2 * _L7 / _RF * sg * (_RMU - _RMD)
            m2eta = ((Cu + Cd - 2 * Cs) * K
                     + 32 * _RB0**2 * _L7 / _RF * sg * (_RMU + _RMD - 2 * _RMS)) / np.sqrt(3)
            g_gen = (-_RALPHA / np.pi / _RF * m2pi / (mphi**2 - _RMPI0**2)
                     - _RALPHA / np.sqrt(3) / np.pi / _RF * m2eta / (mphi**2 - _RMETA**2))
            # SU(3)-singlet direction Cu=Cd=Cs: induced via eta' (phi_gg_coupling)
            F_etap = np.sqrt(_RALPHA**2 * _RMETAP**3 / (64 * np.pi**3 * 4.34e-6))
            m2etap = (np.sqrt(6) / 3 * sg
                      * (_RB0 * _RF + (_L6 + _L7) * 32 * _RB0**2 / _RF * (_RMU + _RMD + _RMS)))
            g_sing = _RALPHA / (np.pi * F_etap) * m2etap / (mphi**2 - _RMETAP**2)
            gphigg = np.where((Cu == Cd) & (Cd == Cs), g_sing, g_gen)
            vis = gphigg**2 * mphi**3 / (64 * np.pi)
        if kind in ("SxS", "SxP"):
            xx = np.where(mphi > 2 * mx, gx**2 * d**1.5 / (8 * np.pi * mphi**2), 0.0)
        else:
            xx = gx**2 / (4 * np.pi) * np.sqrt(np.clip(mphi**2 / 4 - mx**2, 0.0, None))
    return vis, xx


def kphi_invisible_fraction(kind, mx, mphi, gx, Cu, Cd, Cs):
    """Chance an on-shell phi from K+ -> pi+ phi looks invisible at NA62:
    P_inv = 1 - BR_vis (1 - exp(-L/ell)), ell = (p_phi/mphi) c / Gamma_tot -- it
    decays to chi chibar, or is still undecayed after L = 100 m (K decay point to
    the LKr) at E_phi = 50 GeV (75 GeV K+ less a ~25 GeV pi+).  AxA: 1 (A' decays
    only to chi chibar).  As mesonandCMB_ratio.m."""
    if kind == "AxA":
        return np.ones(np.broadcast(mx, mphi, gx, Cu).shape)
    vis, xx = _relic_phi_widths(kind, mx, mphi, gx, Cu, Cd, Cs)
    tot = vis + xx
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        ell = np.sqrt(50.0**2 - mphi**2) / mphi * 1.973269804e-16 / tot
        # (xx + vis e^{-L/ell})/tot, not 1 - BR_vis (1 - e^{-L/ell}): the latter
        # cancels to exactly 0 once BR(chi chibar) < 1e-16.
        P = (xx + vis * np.exp(-100.0 / ell)) / tot
    return np.where((tot > 0) & (vis > 0), P, 1.0)


def kphi_ratio(mx, mphi, gx, gu, gd, gs, interaction):
    """Gamma(K+ -> pi+ phi) * P_inv / (K+->pi+ + invisible limit); r > 1
    excluded.  NaN where mphi >= m_K+- - m_pi+- (channel closed)."""
    kind = kind_of(interaction)
    mx, mphi, gx, gu, gd, gs = np.broadcast_arrays(
        *(np.asarray(a, dtype=float) for a in (mx, mphi, gx, gu, gd, gs)))
    if kind is None:
        return np.full(mx.shape, np.nan)
    Cu, Cd, Cs = gu, gd, gs
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        kin = np.sqrt(np.clip((_MK**2 + _MPI**2 - mphi**2)**2 - 4 * _MK**2 * _MPI**2, 0.0, None))
        if kind in ("SxS", "PxS"):
            delta = -_GF * _F**2 * _VUS * _VUD * _G8 / np.sqrt(2)
            cpik = delta * _B0 * ((Cs + Cu) * _MPI - (Cu + Cd) * _MK) / (_MK - _MPI)
            Gamma = cpik**2 * kin / (16 * np.pi * _MK**3) / _HBAR
        elif kind in ("SxP", "PxP"):
            # theta_j from _sppp_theta (sppp_mixing.m: NLO, pi0, eta AND eta'); the
            # eta' vertex (sqrt(6)/3)(2 + 3 g8s/g8)(mK^2 - mpi^2) is from eq.
            # Klagrangian, missing until 2026-11-17.  As mesonandCMB_ratio.m.
            th = _sppp_theta(np.atleast_1d(mphi), np.atleast_1d(Cu), np.atleast_1d(Cd),
                             np.atleast_1d(Cs)).reshape((3,) + mphi.shape)
            Epi = (_MK**2 + _MPI**2 - mphi**2) / (2 * _MK)
            gphipik = (-(2 * _MK * Epi - _MK**2) * th[0]
                       - (_MK**2 - 6 * _MK * Epi + 2 * _MPI**2) / np.sqrt(3) * th[1]
                       - np.sqrt(6) / 3 * (2 + 3 * _G8S / _G8) * (_MK**2 - _MPI**2) * th[2])
            # No Epi factor: the draft's bracket [(mK^2+mpi^2-mphi^2)/(2mK)] is the
            # EVALUATION ARGUMENT Epi* of gphipik (used above), not a multiplicative
            # one -- as a factor the width is dimensionally wrong (GeV^2), and the
            # NWA limit of the 3-body rate through the phi pole reproduces this form
            # with no Epi.  Removed from mesonandCMB_ratio.m/.jl 2026-07-22; this
            # port kept it until 2026-08-04 (~4x at small mphi).
            Gamma = ((_GF * _F * _G8 * _VUS * _VUD)**2 / (32 * np.pi * _MK**3) * gphipik**2
                     * kin / _HBAR)
        else:  # AxA
            # Covariant polarisation sum: (lak + lap)^2 lambda/(4 mA^2); the rest-frame
            # bracket used until 2026-11-17 was x1.69-1.81 low.  As mesonandCMB_ratio.m.
            lak = _G8 * (Cs + Cu) + _G8S * (Cu + Cd + Cs)
            lap = _G8 * (Cd + Cu) + _G8S * (Cu + Cd + Cs)
            Gamma = (_GF**2 * _F**4 * _VUS**2 * _VUD**2 * (lak + lap)**2 * kin**3
                     / (32 * np.pi * mphi**2 * _MK**3) / _HBAR)
        r = Gamma * kphi_invisible_fraction(kind, mx, mphi, gx, Cu, Cd, Cs) / _KINV
    # x NA62 acceptance at m_miss^2 = mphi^2 (kaon_invisible.na62_acceptance; 2026-11-17)
    from kaon_invisible import na62_acceptance
    r = r * na62_acceptance(mphi**2)
    return np.where(mphi >= KPHI_THRESHOLD, np.nan, r)


def _sigmav(kind, mx, mphi, gx, Cu, Cd, Cs, cmb=False):
    # Low-velocity <sigma v> [cm^3/s] used by the CMB / INTEGRAL bounds.  Analytic
    # threshold (s = 4 m_chi^2) expressions from mesonandCMBlimits.py.
    # cmb=True adds the s-wave three-body annihilation (_threebody_cmb_sv), each
    # channel weighted by its deposited-energy fraction relative to gamma gamma.
    sapprox = 4 * mx**2
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        if kind == "SxS":
            Dnew = 0.5*(_B0 * (Cu + Cd) / _MPI**2 + _B0 * (Cu + Cs) / _MK**2)**2
            jfac = (2.42e17) / (1.30e23)
            # 4pi: the source expression is dsigma/dOmega.
            # _HBARC**2 * _CC (ONE factor of c), as in every other branch: the
            # bracket is [GeV^-2], so GeV^-2 -> cm^2 needs (hbar c)^2 and
            # cm^2 -> cm^3/s a single c.  The source (mesonandCMBlimits.py:418)
            # had `*c * hbarc^2 * c`; the second c has nothing to convert and made
            # this branch 3e10 too large.  Cf. SxP below -- same p-wave structure,
            # same jfac, one _CC.  Fixed 2026-08-04.
            return (4 * np.pi * (gx**2 * _E**4 * Dnew) / (589824 * np.pi**6)
                    * (8 * mx**4 * jfac) / (4 * mx**2 - mphi**2)**2 * _HBARC**2 * _CC)
        if kind == "PxS":
            Dnew = 0.5*(_B0 * (Cu + Cd) / _MPI**2 + _B0 * (Cu + Cs) / _MK**2)**2
            # gamma gamma + chi chibar + meson pairs
            GphiPS = (_phi_gg_xx_width(kind, mx, mphi, gx, Cu, Cd, Cs)
                      + meson_pair_width(mphi, Cu, Cd, Cs))
            return (4 * np.pi * gx**2 * _ALPHA**2 * Dnew / (36864 * np.pi**4) * (4 * mx**2)**2
                    / ((4 * mx**2 - mphi**2)**2 + (mphi * GphiPS)**2) * 2 * _HBARC**2 * _CC)
        if kind in ("SxP", "PxP"):
            isPP = kind == "PxP"

            def _mediator_width(gchi, ggam):
                Ggg = ggam**2 * mphi**3 / (64 * np.pi)
                if isPP:
                    Gxx = gchi**2 * np.sqrt(np.clip(mphi**2 / 4 - mx**2, 0.0, None)) / (4 * np.pi)
                else:
                    Gxx = gchi**2 * np.clip(mphi**2 - 4 * mx**2, 0.0, None)**1.5 / (8 * np.pi * mphi**2)
                return np.where(mphi**2 >= 4 * mx**2, Gxx + Ggg, Ggg) + G3pi

            # phi -> 3pi, as in the relic widths (phi_partial_widths.m); 0 when off
            G3pi = _phi_3pi_width(mphi, Cu, Cd, Cs)

            # --- generic pole set {pi0, eta, phi'} ---------------------------
            mpi0phi2, metaphi2 = _pp_couplings(Cu, Cd, Cs)
            gphichi = gx
            gpichi = gx * mpi0phi2 / (mphi**2 - _MPI0**2)
            getachi = gx * metaphi2 / (mphi**2 - _META**2)
            ggammapi = _E**2 / (4 * np.pi**2 * _F)
            ggammaeta = _E**2 / (4 * np.sqrt(3) * np.pi**2 * _F)
            ggammaphi = (-ggammapi * mpi0phi2 / (mphi**2 - _MPI0**2)
                         - ggammaeta * metaphi2 / (mphi**2 - _META**2))
            Gphi = _mediator_width(gphichi, ggammaphi)
            amp_gen = (gpichi * ggammapi / (sapprox - _MPI0**2 + 1j * _MPI0 * _GAMMAPI0)
                       + getachi * ggammaeta / (sapprox - _META**2 + 1j * _META * _GETA)
                       + gphichi * ggammaphi / (sapprox - mphi**2 + 1j * mphi * Gphi))

            # --- SU(3)-singlet direction Cu=Cd=Cs: {eta', phi'} --------------
            # There the pi0/eta mixings and the induced phi'gamma coupling all
            # vanish, so amp_gen == 0 and chi would decouple from photons.  The
            # surviving channel is the singlet eta';
            # mirrors sppp_couplings.m so the relic and indirect pipelines agree.
            _GETAP_GG, _GETAP_TOT = 4.34e-6, 1.88e-4
            F_etap = np.sqrt(_ALPHA**2 * _METAP**3 / (64 * np.pi**3 * _GETAP_GG))
            g_etapgamma = _ALPHA / (np.pi * F_etap)   # 1/pi: see sppp_couplings.m
            # eta' singlet mixing mass: sqrt(6)/3 * sum(g) * (B0 F + (L6+L7) NLO), the draft body's form (4.3x the old B0*F/sqrt(6), whose appendix was commented out 2026-08-12)
            m_mix = (np.sqrt(6) / 3 * (Cu + Cd + Cs)
                     * (_B0 * _F + (_L6 + _L7) * 32 * _B0**2 / _F * (_MU + _MD + _MS)))
            denom = (mphi**2 - _METAP**2)
            g_etapx = -gx * m_mix / denom
            gphiprimegamma = g_etapgamma * m_mix / denom
            if isPP:
                Gam_etapxx = g_etapx**2 * np.sqrt(np.clip(_METAP**2 / 4 - mx**2, 0.0, None)) / (4 * np.pi)
            else:
                Gam_etapxx = g_etapx**2 * np.clip(_METAP**2 - 4 * mx**2, 0.0, None)**1.5 / (8 * np.pi * _METAP**2)
            Gphi_s = _mediator_width(gx, gphiprimegamma)
            amp_sing = (g_etapx * g_etapgamma
                        / (sapprox - _METAP**2 + 1j * _METAP * (_GETAP_TOT + Gam_etapxx))
                        + gx * gphiprimegamma / (sapprox - mphi**2 + 1j * mphi * Gphi_s))

            amp = np.where((Cu == Cd) & (Cd == Cs), amp_sing, amp_gen)
            # 4pi: the source expression is dsigma/dOmega.
            if isPP:
                sv = 4 * np.pi * sapprox**2 / (512 * np.pi**2) * np.abs(amp)**2 * 2 * _HBARC**2 * _CC
                if cmb:
                    Gtot = np.where((Cu == Cd) & (Cd == Cs), Gphi_s, Gphi)
                    sv = sv + _threebody_cmb_sv("PxP", mx, mphi, gx, Cu, Cd, Cs, Gtot)
                return sv
            jfac = (2.42e17) / (1.30e23)
            return (4 * np.pi * sapprox**1.5 * mx * jfac / (512 * np.pi**2)
                    * np.abs(amp)**2 * _HBARC**2 * _CC)
        # AxA
        Getapgg = (2.22 / 100) * 0.196e-3
        Fetap = np.sqrt((_ALPHA**2 * _METAP**3) / (64 * np.pi**3 * Getapgg))
        lampi0gamma = _E**2 / (4 * np.pi**2 * _F)
        # data-driven F_eta (2026-08-11 draft Eq. lambdagamma); cf. axial_couplings.m
        F_eta = np.sqrt(_ALPHA**2 * _META**3 / (64 * np.pi**3 * 5.16e-7))
        lametagamma = _E**2 / (4 * np.pi**2 * F_eta)
        lameta0gamma = _E**2 / (4 * np.pi**2 * Fetap)   # 4*pi**2, as pi0/eta above
        mpi0Ap = _F * ((Cu - Cd) * (1 + (16 * _B0 * _L4) / _F**2 * (_MU + _MD + _MS))
                       + (16 * _B0 * _L5 * (Cu * _MU - Cd * _MD)) / _F**2)
        metaAp = (np.sqrt(3) * _F) / 3 * ((Cu + Cd - 2 * Cs) * (1 + (16 * _B0 * _L4) / _F**2 * (_MU + _MD + _MS))
                                          + 1 / _F**2 * 16 * _B0 * _L5 * (Cu * _MU + Cd * _MD - 2 * Cs * _MS))
        meta0Ap = (np.sqrt(6) * _F) / 3 * ((Cu + Cd + Cs) * (1 + (16 * _B0 * _L4) / _F**2 * (_MU + _MD + _MS))
                                           + 1 / _F**2 * 16 * _B0 * _L5 * (Cu * _MU + Cd * _MD + Cs * _MS))
        lampichi = gx * mpi0Ap / mphi**2
        lametachi = gx * metaAp / mphi**2
        lameta0chi = gx * meta0Ap / mphi**2
        amp = (lampichi * lampi0gamma / (sapprox - _MPI0**2 + 1j * _MPI0 * _GAMMAPI0)
               + lametachi * lametagamma / (sapprox - _META**2 + 1j * _META * _GETA)
               + lameta0chi * lameta0gamma / (sapprox - _METAP**2 + 1j * _METAP * _GETA0))
        # 4pi: the source expression is dsigma/dOmega.
        sv = 4 * np.pi * mx**2 * sapprox**2 / (128 * np.pi**2) * np.abs(amp)**2 * 2 * _HBARC**2 * _CC
        if cmb:
            sv = sv + _threebody_cmb_sv("AxA", mx, mphi, gx, Cu, Cd, Cs, None)
        return sv


# ---------------------------------------------------------------------------
# Three-particle annihilation (LOOP_3PI / LOOP_PIPIPHI, read from the
# environment like the relic code).  Ports of matlab/threepi_weight.m,
# pipiphi_weight.m, axial3pi_weight.m and dalitz_moments.m, with the relic
# pipeline's constants (as in the meson-invisible section below), so the
# channels are the ones the relic abundance used.
#
# CMB only, and only the s-wave pieces: P(x)P chi chi -> 3pi and pi pi phi, and
# the axial J = 0 chi chi -> 3pi.  S(x)P carries (s - 4mx^2)/s and the axial
# J = 1 part carries beta^2, so both vanish at s = 4 mx^2; the INTEGRAL curve is
# a gamma-line bound, which a pion continuum does not feed.
#
# The Slatyer curve bounds gamma gamma, which deposits all 2mx.  Each channel is
# weighted by the fraction of its energy that ends up in photons and e+-,
# assuming f_eff(e) ~ f_eff(gamma) at these energies:
#   pi0 -> gamma gamma: all of it.  pi+- -> mu nu keeps (mpi^2 + mmu^2)/(2 mpi^2)
#   = 0.787 of the energy, mu -> e nu nu passes 7/20 of that: 0.275.
#   3pi:     (1 + 2 x 0.275)/3 = 0.52, the pi+pi-pi0 value, used for 3pi0 too.
#   pi pi phi: the pions' 2/3 share at the charged value, 2/3 x 0.275 = 0.18; the
#            phi's energy is not counted (it may decay after recombination).
# Both are lower bounds on the deposited energy, so the bound is conservative.
# ---------------------------------------------------------------------------
KAPPA_3PI = 0.52
KAPPA_PIPIPHI = 0.18
_T_MPI0, _T_MPIPM = 0.13498, 0.13957
_T_META, _T_METAP = 0.547862, 0.95778
_T_GETA, _T_GETAP = 1.31e-6, 1.88e-4
_T_L4, _T_L5 = 1.87e-3, 1.22e-3


def _threepi_on():
    return os.environ.get("LOOP_3PI", "") == "1"


def _pipiphi_on():
    return os.environ.get("LOOP_PIPIPHI", "") == "1"


def _dalitz_moments(s, m1, m2):
    """I_n(s) = int dE Y(E) X(E)^n over the m1 energy; port of dalitz_moments.m.
    s and m1 broadcast; the other two particles have mass m2."""
    n = 32
    th = np.pi * (np.arange(1, n + 1) - 0.5) / n
    t = (1 - np.cos(th)) / 2
    dt = np.sin(th) / 2 * np.pi / n
    s = np.asarray(s, dtype=float)[..., None]
    m1 = np.asarray(m1, dtype=float)[..., None]
    rs = np.sqrt(s)
    Elo = m1
    Ehi = (s + m1**2 - 4 * m2**2) / (2 * rs)
    E = Elo + (Ehi - Elo) * t
    w = s + m1**2 - 2 * E * rs
    Y = np.sqrt(np.clip((E**2 - m1**2) * (1 - 4 * m2**2 / w), 0.0, None))
    X = w - _T_MPI0**2
    wt = (Ehi - Elo) * dt * Y
    return wt.sum(-1), (wt * X).sum(-1), (wt * X**2).sum(-1)


def _threepi_weight(s, Cu, Cd, Cs):
    """W(s) over pi+pi-pi0 and 3pi0, mediator propagator stripped; port of
    threepi_weight.m (eta pole generically, eta' for Cu = Cd = Cs)."""
    B0, f, dm = _RB0, _RF, _RMD - _RMU
    mpi0, mpm = _T_MPI0, _T_MPIPM
    D = Cu - Cd
    deg = (Cu == Cd) & (Cd == Cs)
    mj = np.where(deg, _T_METAP, _T_META)
    Gj = np.where(deg, _T_GETAP, _T_GETA)
    aj = np.where(deg, B0 * f * np.sqrt(2 / 3) * (Cu + Cd + Cs),
                  B0 * f * (Cu + Cd - 2 * Cs) / np.sqrt(3))
    cj = np.where(deg, -np.sqrt(6) / 3 * B0 * dm / f**2, -B0 * dm / (np.sqrt(3) * f**2))
    bj = np.where(deg, -np.sqrt(6) / 9 * B0 * dm / f**2, -B0 * dm / (3 * np.sqrt(3) * f**2))
    d2 = mj**2 - mpi0**2
    pole = aj / (s - mj**2 + 1j * mj * Gj)
    A = B0 * D / f * s / (s - mpi0**2) + pole * cj
    W = np.where(s > 9 * mpi0**2, np.abs(A)**2 * _dalitz_moments(s, mpi0, mpi0)[0] / 6, 0.0)
    P = pole * bj * (1 - s / d2)
    Q = B0 * D / f / (s - mpi0**2) + pole * bj * 3 / d2
    I0, I1, I2 = _dalitz_moments(s, mpi0, mpm)
    Wpm = np.abs(P)**2 * I0 + 2 * np.real(P * np.conj(Q)) * I1 + np.abs(Q)**2 * I2
    return W + np.where(s > (2 * mpm + mpi0)**2, Wpm, 0.0)


def _pipiphi_weight(s, mphi, Cu, Cd, Cs):
    """W(s) over pi+pi-phi and pi0pi0phi, every propagator kept; port of
    pipiphi_weight.m."""
    B0, f = _RB0, _RF
    b, p2 = B0 * f, _T_MPI0**2
    mpm, metap, meta = _T_MPIPM, _T_METAP, _T_META
    D, Sg = Cu - Cd, Cu + Cd
    deg = (Cu == Cd) & (Cd == Cs)
    Rf = 1 / (s - mphi**2)
    Rp = 1 / (s - p2)
    # degenerate g ~ [1 1 1]: eta' and phi poles
    t0 = -b * np.sqrt(2 / 3) * (Cu + Cd + Cs) / (mphi**2 - metap**2)
    Rm = 1 / (s - metap**2 + 1j * metap * _T_GETAP)
    Nphi = -2 * t0 * (np.sqrt(6) * b * Sg - p2 * t0) / (3 * f**2)
    Nj = t0 * (np.sqrt(6) * b * Sg - 2 * p2 * t0) / (3 * f**2)
    # generic: pi0, eta and phi poles
    tp = -b * D / (mphi**2 - p2)
    te = -b * (Cu + Cd - 2 * Cs) / np.sqrt(3) / (mphi**2 - meta**2)
    Re = 1 / (s - meta**2 + 1j * meta * _T_GETA)
    a = mphi**2 + p2 + s - 3 * p2
    P_gen = (Rf * (-2 * np.sqrt(3) * b * Sg * te - 2 * b * D * tp + p2 * te**2 - a * tp**2) / (3 * f**2)
             + Rp * tp * (b * D + a * tp) / (3 * f**2)
             + Re * te * (np.sqrt(3) * b * Sg - p2 * te) / (3 * f**2))
    Q_gen = tp**2 / f**2 * (Rf - Rp)
    A_deg = Nphi * Rf + Nj * Rm
    P = np.where(deg, A_deg, P_gen)
    Q = np.where(deg, 0.0, Q_gen)
    I0, I1, I2 = _dalitz_moments(s, mphi, mpm)
    Wpm = np.abs(P)**2 * I0 + 2 * np.real(P * np.conj(Q)) * I1 + np.abs(Q)**2 * I2
    A00 = np.where(deg, A_deg,
                   -Rf * (2 * np.sqrt(3) * b * Sg * te + 6 * b * D * tp - p2 * te**2 - 3 * p2 * tp**2) / (3 * f**2)
                   + Rp * tp * (b * D - p2 * tp) / f**2
                   + Re * te * (np.sqrt(3) * b * Sg - p2 * te) / (3 * f**2))
    W00 = np.abs(A00)**2 * _dalitz_moments(s, mphi, _T_MPI0)[0] / 2
    return (np.where(s > (2 * mpm + mphi)**2, Wpm, 0.0)
            + np.where(s > (2 * _T_MPI0 + mphi)**2, W00, 0.0))


def _axial3pi_wl(s, Cu, Cd, Cs):
    """Longitudinal (J = 0) Dalitz weight W_L of the axial 3pi current; port of
    axial3pi_weight.m with axial_mixing.m's m_{pi0 A'}."""
    f, mpi0, mpm = _RF, _T_MPI0, _T_MPIPM
    corr4 = 1 + 16 * _RB0 * _T_L4 / f**2 * (_RMU + _RMD + _RMS)
    m_pi0A = f * ((Cu - Cd) * corr4 + 16 * _RB0 * _T_L5 * (Cu * _RMU - Cd * _RMD) / f**2)
    cL = (m_pi0A / f)**2 * mpi0**4 / f**2
    WL = np.where(s > 9 * mpi0**2,
                  cL * s**2 / (s - mpi0**2)**2 * _dalitz_moments(s, mpi0, mpi0)[0] / 6, 0.0)
    I2 = _dalitz_moments(s, mpi0, mpm)[2]
    return WL + np.where(s > (2 * mpm + mpi0)**2, cL * I2 / (s - mpi0**2)**2, 0.0)


def _phi_3pi_width(mphi, Cu, Cd, Cs):
    """Gamma(phi -> 3pi) [GeV] = W(mphi^2)/(64 pi^3 mphi); 0 below 3 mpi0 or when off."""
    if not _threepi_on():
        return np.zeros(np.shape(mphi))
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        W = _threepi_weight(mphi**2, Cu, Cd, Cs)
        return np.where(mphi > 3 * _T_MPI0, W / (64 * np.pi**3 * mphi), 0.0)


def _threebody_cmb_sv(kind, mx, mphi, gx, Cu, Cd, Cs, Gphi):
    """Energy-weighted s-wave <sigma v> [cm^3/s] of the three-body channels at
    s = 4 mx^2 (v_rel = 2 beta), to add to the gamma gamma <sigma v> for the CMB.
    Gphi: the mediator's total width (PxP; unused for AxA)."""
    s0 = 4 * mx**2
    sv = np.zeros(np.shape(mx))
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        if kind == "PxP":
            if _threepi_on():
                D2 = (s0 - mphi**2)**2 + mphi**2 * Gphi**2
                sv = sv + KAPPA_3PI * gx**2 * _threepi_weight(s0, Cu, Cd, Cs) / (64 * np.pi**3 * D2)
            if _pipiphi_on():
                sv = sv + KAPPA_PIPIPHI * gx**2 * _pipiphi_weight(s0, mphi, Cu, Cd, Cs) / (64 * np.pi**3)
        elif kind == "AxA" and _threepi_on():
            sv = sv + KAPPA_3PI * gx**2 * mx**2 * _axial3pi_wl(s0, Cu, Cd, Cs) / (16 * np.pi**3 * s0**2 * mphi**4)
    # NaN -> 0: the sum feeds the gamma gamma <sigma v>, which a failed
    # three-body evaluation must not erase.
    return np.nan_to_num(sv, nan=0.0) * _HBARC**2 * _CC


def _indirect_ratio(mx, mphi, gx, gu, gd, gs, interaction, which):
    kind = kind_of(interaction)
    mx, mphi, gx, gu, gd, gs = np.broadcast_arrays(
        *(np.asarray(a, dtype=float) for a in (mx, mphi, gx, gu, gd, gs)))
    if kind is None or (which == "cmb" and kind in ("SxS", "SxP")):
        return np.full(mx.shape, np.nan)
    sv = _sigmav(kind, mx, mphi, gx, gu, gd, gs, cmb=(which == "cmb"))
    lim = _curve(mx, _SLATYER_X, _SLATYER_Y) if which == "cmb" \
        else _curve(mx, _INTEGRAL_X, _INTEGRAL_Y)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return sv / lim


def cmb_ratio(mx, mphi, gx, gu, gd, gs, interaction):
    """<sigma v> / Slatyer CMB bound; r > 1 excluded.  NaN for SxS/SxP."""
    return _indirect_ratio(mx, mphi, gx, gu, gd, gs, interaction, "cmb")


def integral_ratio(mx, mphi, gx, gu, gd, gs, interaction):
    """<sigma v> / INTEGRAL bound; r > 1 excluded.  NaN outside the data range."""
    return _indirect_ratio(mx, mphi, gx, gu, gd, gs, interaction, "integral")


# ---------------------------------------------------------------------------
# Neutral meson -> invisible.  Port of mesonandCMB_ratio.m's
# meson_invisible_widths, which takes its couplings from the relic pipeline's
# sppp_mixing.m (SxP/PxP) and axial_couplings.m (AxA); their constants are used
# here, NOT the reference-script constants above, so all three ports agree.
# ---------------------------------------------------------------------------
MESONINV_LIMITS = (3.43e-17, 3.14e-11, 3.95e-8)    # Gamma(j -> invisible) [GeV]
_MJ = np.array([0.13498, 0.547862, 0.95778])       # pi0, eta, eta'
_RB0 = 0.13957**2 / (0.0022 + 0.00469)
_RF = 0.0924
_RMU, _RMD, _RMS = 0.0022, 0.00469, 0.095


def _sppp_theta(mphi, Cu, Cd, Cs):
    # sppp_mixing.m: theta_j = m2_{j phi} / (mphi^2 - m_j^2), shape (3, ...)
    summ = _RMU + _RMD + _RMS
    K = _RB0 * _RF + 32 * _RB0**2 * _L6 / _RF * summ
    sg = Cu + Cd + Cs
    m2 = np.stack([
        (Cu - Cd) * K + 32 * _RB0**2 * _L7 / _RF * sg * (_RMU - _RMD),
        ((Cu + Cd - 2 * Cs) * K + 32 * _RB0**2 * _L7 / _RF * sg * (_RMU + _RMD - 2 * _RMS)) / np.sqrt(3),
        np.sqrt(6) / 3 * sg * (_RB0 * _RF + (_L6 + _L7) * 32 * _RB0**2 / _RF * summ)])
    return m2 / (mphi**2 - _MJ[:, None]**2)


def _axial_lam_chi(mphi, gx, Cu, Cd, Cs):
    # axial_couplings.m: lambda_{j chi} = lambda_a m_{jA'} / m_A'^2, shape (3, ...)
    corr4 = 1 + 16 * _RB0 * _L4 / _RF**2 * (_RMU + _RMD + _RMS)
    m = np.stack([
        _RF * ((Cu - Cd) * corr4 + 16 * _RB0 * _L5 * (Cu * _RMU - Cd * _RMD) / _RF**2),
        np.sqrt(3) * _RF / 3 * ((Cu + Cd - 2 * Cs) * corr4
                                + 16 * _RB0 * _L5 * (Cu * _RMU + Cd * _RMD - 2 * Cs * _RMS) / _RF**2),
        np.sqrt(6) * _RF / 3 * ((Cu + Cd + Cs) * corr4
                                + 16 * _RB0 * _L5 * (Cu * _RMU + Cd * _RMD + Cs * _RMS) / _RF**2)])
    return gx / mphi**2 * m


def meson_invisible_widths(mx, mphi, gx, gu, gd, gs, interaction):
    """Gamma(j -> chi chibar) [GeV], shape (3, N) for j = pi0, eta, eta'.  NaN where
    the channel is closed (mx >= m_j/2) or the model has no meson-chi coupling."""
    kind = kind_of(interaction)
    mx, mphi, gx, gu, gd, gs = (a.ravel() for a in np.broadcast_arrays(
        *(np.asarray(a, dtype=float) for a in (mx, mphi, gx, gu, gd, gs))))
    mj = _MJ[:, None]
    b2 = 1 - 4 * mx**2 / mj**2                        # beta^2; open where > 0
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        if kind in ("SxP", "PxP"):
            # g_{j chi} = gx theta_j; pseudoscalar (PxP, beta^1) / scalar (SxP, beta^3)
            p = 0.5 if kind == "PxP" else 1.5
            G = (gx * _sppp_theta(mphi, gu, gd, gs))**2 * mj * np.maximum(b2, 0) ** p / (8 * np.pi)
        elif kind == "AxA":
            lam = _axial_lam_chi(mphi, gx, gu, gd, gs)
            # No validity floor here: the MATLAB/Julia ports probe below
            # axial_couplings' meson-width floor and rescale (lambda ~ gx), which
            # is exactly this unclipped expression.
            G = lam**2 * mx**2 * np.sqrt(np.maximum(0.0, mj**2 / 4 - mx**2)) / np.pi
        else:
            return np.full((3, mx.size), np.nan)
    return np.where(b2 > 0, G, np.nan)


def meson_invisible_ratio(mx, mphi, gx, gu, gd, gs, interaction, which="max"):
    """Gamma(j -> chi chibar) / limit_j; r > 1 excluded.  which = 'pi0' | 'eta' |
    'etap' | 'max' (largest of the three; NaN only if all three are NaN)."""
    shape = np.broadcast(*(np.asarray(a) for a in (mx, mphi, gx, gu, gd, gs))).shape
    r = meson_invisible_widths(mx, mphi, gx, gu, gd, gs, interaction) \
        / np.array(MESONINV_LIMITS)[:, None]
    if which == "max":
        with np.errstate(invalid="ignore"):
            out = np.where(np.isnan(r).all(axis=0), np.nan,
                           np.nanmax(np.where(np.isnan(r), -np.inf, r), axis=0))
    else:
        out = r[{"pi0": 0, "eta": 1, "etap": 2}[which]]
    return out.reshape(shape)
