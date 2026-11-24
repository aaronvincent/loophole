#!/usr/bin/env python3
"""Build a combined Neff + external-constraint SIDECAR for a Loophole scan.

A relic-only scan (neff_*_sigma = 0) is mapped by the sampler, then ALL offline
likelihoods are applied here -- WITHOUT re-running and WITHOUT touching the raw
scan file.

Applied likelihoods:
  relic      recomputed from output::Oh2 (identical to the raw scan's term)
  neff_bbn   two-sided Gaussian on dNeff_bbn  (BBN epoch; the worker writes NaN,
             which counts as no contribution)
  neff_cmb   two-sided Gaussian on dNeff_cmb  (CMB epoch)
  sn1987a    SN1987A cooling -- closed exclusion band in (log10 mphi, log10 g_N)
  neutronstar  NS cooling (Damiano) -- one-sided upper bound on scalar g_N (SS/PS)
  stellar      Stellar cooling, TRGB + horizontal branch (Hardy & Lasenby) --
               one-sided upper bound on alpha_phiNN = g^2/(4 pi), scalar
               coupling, g_bar = (hp+hn)/2 (SS/PS)
  fifthforce   5th-force bound: g_A^2 vs range (axialVector), or the
               pseudoscalar g_P^2 / |g_P g_N| spin-spin curves (SP/PP)
  bullet       Bullet-cluster DM self-interaction bound sigma_self/m < bound
  kaon         Rare K+ -> pi+ chi chi vs the K+->pi+invisible limit (mx<0.177 GeV)
  kphi         K+ -> pi+ phi (on-shell mediator) vs the same limit (mphi<0.354 GeV),
               weighted by the chance the phi looks invisible at NA62
  cmb          DM <sigma v> vs the Slatyer CMB bound (~gx^2; not SxS/SxP)
  integral     DM <sigma v> vs the INTEGRAL telescope bound (~gx^2; light DM)
  mesoninv     pi0 / eta / eta' -> chi chibar vs the invisible-decay limits, via
               the meson-mediator mixing (SxP/PxP/AxA; ~gx^2)

The raw scan <tag>.hdf5 is opened READ-ONLY.  The sidecar <tag>_postprocessed.hdf5
is self-contained: all columns the plot scripts need, contiguous (no GAMBIT
chunks=(100,)) so reads are fast.  Atomic write (temp + rename) -- crash-safe.

Sidecar columns:
    LogLike            = relic + BBN + CMB + SN1987A + NS + stellar + 5th-force
    LogLike_isvalid    (copied from the raw scan)
    LogLike_orig       the original relic-only LogLike
    LogLike_relic / LogLike_neff_bbn / LogLike_neff_cmb
    LogLike_sn1987a / LogLike_neutronstar / LogLike_stellar / LogLike_fifthforce
    LogLike_bullet / LogLike_kaon / LogLike_kphi / LogLike_cmb / LogLike_integral
    LogLike_mesoninv
    excluded_sn1987a / excluded_neutronstar / excluded_stellar /
    excluded_fifthforce / excluded_bullet / excluded_kaon / excluded_kphi /
    excluded_cmb / excluded_integral / excluded_mesoninv / excluded_any
                       crisp boolean masks (uint8)
    <scan axes>(+_isvalid), output::Oh2/dNeff_bbn/dNeff_cmb/sigma_chiN  copied
    input::relic_* , input::neff_bbn/cmb_* , input::ext_*  (applied params)

Usage:
    python postprocess.py IN [group]
        IN     a raw scan file, a directory, or a glob
        group  HDF5 group (default: data)

    Batch: skips files lacking '/<group>/LogLike' and (unless POSTPROC_FORCE=1)
    raw files that already have a <tag>_postprocessed.hdf5 sidecar.

Env overrides:
    NEFF_BBN_MU    (default -0.14)    NEFF_BBN_SIGMA  (0.21)
    NEFF_CMB_MU    (default -0.054)   NEFF_CMB_SIGMA  (0.17)
    EXT_SIGMA_SN, EXT_SIGMA_NS, EXT_SIGMA_STELLAR (Hardy-Lasenby, log10 alpha),
    EXT_SIGMA_AXIAL, EXT_SIGMA_FF_PSEUDO (5th force), EXT_SIGMA_BULLET,
    EXT_SIGMA_KAON, EXT_SIGMA_KPHI, EXT_SIGMA_CMB, EXT_SIGMA_INTEGRAL,
    EXT_SIGMA_MESONINV
                   all default -1: a HARD cut (soft widths are in log10(r)
                   around r=1, or log10 distance for the region constraints)

    Every EXT_SIGMA_* knob takes three kinds of value:
        < 0   HARD cut: -EXT_HARD_PENALTY (default 1e6) on every excluded point,
              so LogLike agrees exactly with the excluded_* flag  (the default).
              Right for a sharp experimental limit (e.g. a measured branching
              ratio), where a factor-of-2 "1 sigma" tolerance has no meaning.
        > 0   soft one-sided Gaussian of that width, in decades
        = 0   constraint OFF -- no penalty at all (NOT a sharp cut)
        The applied value is recorded per row in input::ext_sigma_*, so a
        sidecar says which mode produced it.
        e.g.  EXT_SIGMA_KAON=0.3 python postprocess.py ...   (soft kaon bound)
    EXT_INTERACTION  override the auto-detected interaction name
    POSTPROC_FORCE   set to 1 to rebuild an existing sidecar
"""

import glob as globmod
import os
import re
import sys
import time

import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kaon_invisible
import mesonandCMB

# Big HDF5 chunk cache: GAMBIT's chunks=(100,) make raw column reads slow on a
# network FS without it.
_RDCC = dict(rdcc_nbytes=512 * 1024**2, rdcc_nslots=2_000_003)

_DATA_DIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))

_AXES = ["input::mx", "input::mphi", "input::gx", "input::gq",
         "input::gu", "input::gd", "input::gs", "output::sigma_chiN"]
_OBS  = ["output::Oh2", "output::dNeff_bbn", "output::dNeff_cmb"]

# Suffixes that are sidecars (never treated as raw inputs).
_SIDECAR_SUFFIXES = ("_postprocessed.hdf5",)


# ---------------------------------------------------------------------------
# Likelihood pieces (mirrors cpp/likelihoods.cpp)
# ---------------------------------------------------------------------------

def relic_loglike(Oh2, mode, mu, sigma):
    if sigma <= 0.0:
        return np.where(np.isfinite(Oh2), 0.0, -np.inf)
    d = (Oh2 - mu) / sigma
    out = -0.5 * d * d
    if mode == 1:
        out = np.where(Oh2 <= mu, 0.0, out)
    return np.where(np.isfinite(Oh2), out, -np.inf)


def neff_loglike(dNeff, mode, mu, sigma):
    if sigma <= 0.0:
        return np.zeros_like(np.asarray(dNeff, dtype=float))
    d = (dNeff - mu) / sigma
    out = -0.5 * d * d
    if mode == 1:
        out = np.where(dNeff > mu, out, 0.0)
    return np.where(np.isfinite(dNeff), out, 0.0)


# Penalty applied by a HARD constraint (sigma < 0) to every excluded point.
# Finite, not -inf: every downstream mask uses np.isfinite(LogLike), so -inf
# would DROP the point rather than exclude it.  1e6 sits far below any real
# LogLike spread in these scans (~1e5), and several stacked hard cuts still
# cannot overflow float64.
HARD_PENALTY = float(os.environ.get("EXT_HARD_PENALTY", 1e6))
# Every EXT_SIGMA_* knob's default: a hard cut.  < 0 hard, > 0 soft width, 0 off.
EXT_SIGMA_DEFAULT = -1.0


def _softregion(depth, sigma):
    """One-sided penalty on `depth` = how far INSIDE the excluded region a point
    is, in decades of log10(ratio) (or distance to the contour, for the region
    constraints).  Three modes, selected by the sign of `sigma`:

        sigma > 0   SOFT: -0.5*(depth/sigma)^2, a one-sided Gaussian of width
                    sigma decades.  depth<=0 (allowed) -> 0.
        sigma == 0  OFF: no penalty at all (used by relic-only scans).
        sigma < 0   HARD: -HARD_PENALTY everywhere inside, 0 outside.  The
                    likelihood then agrees exactly with the excluded_* flag,
                    which is what you want for a sharp experimental limit.

    NOTE the sigma==0 convention: it disables the constraint, it does NOT make
    it sharp.  Use a negative sigma for a hard cut.
    """
    depth = np.asarray(depth, float)
    out = np.zeros_like(depth)
    if sigma == 0:
        return out
    m = np.isfinite(depth) & (depth > 0)
    if sigma < 0:
        out[m] = -HARD_PENALTY
    else:
        out[m] = -0.5 * (depth[m] / sigma) ** 2
    return out


def _smode(sigma):
    """Human-readable mode for a sigma knob, for the run log."""
    return "HARD" if sigma < 0 else ("off" if sigma == 0 else f"soft {sigma:g} dex")


def const_input(g, name, default=None):
    key = f"input::{name}"
    if key in g:
        arr = np.asarray(g[key])
        # Fixed-value params are constant on VALID rows, but row 0 can be an
        # unwritten fill/invalid row (all zeros) for some scanners -- notably
        # paraprof, whose orchestrator rank leaves the first slot empty.  Reading
        # g[key][0] there silently returned relic_sigma=0, which disabled the
        # relic Gaussian (relic_loglike short-circuits to 0 for sigma<=0) and
        # made the whole box look on-shell.  Pick the first VALID finite row.
        if "LogLike_isvalid" in g:
            valid = np.asarray(g["LogLike_isvalid"]).astype(bool)
            vv = arr[valid]
            vv = vv[np.isfinite(vv)]
            if vv.size:
                return float(vv[0])
        return float(arr[0])
    if default is None:
        raise KeyError(f"missing '{key}' and no default given")
    print(f"  note: '{key}' absent -> using default {default}")
    return float(default)


# ---------------------------------------------------------------------------
# Nucleon couplings (mirrors matlab/nucleon_couplings.m)
# ---------------------------------------------------------------------------

MN  = 0.938
FPU, FPD, FPS   = 0.02, 0.026, 0.118
MUP, MDN, MS    = 0.0022, 0.00469, 0.095
DU,  DD,  DS    = 0.847, -0.407, -0.035   # FLAG 2024 N_f=2+1, one set (as nucleon_couplings.m)
DNU, DND, DNS   = DU, DD, DS
HBARC_M         = 1.97e-16                    # GeV*m
HBARC           = 1.97e-14                    # GeV*cm
GEV_TO_G        = 1.783e-24                   # grams per GeV (for sigma/m in cm^2/g)


def _hp(gu, gd, gs):
    return MN * (gu * FPU / MUP + gd * FPD / MDN + gs * FPS / MS)


def _hn(gu, gd, gs):
    """Scalar NEUTRON coupling: isospin swaps the u<->d sigma-term form factors,
    each 1/m_i kept with its own quark.  So this is NOT _hp(gd, gu, gs)."""
    return MN * (gu * FPD / MUP + gd * FPU / MDN + gs * FPS / MS)


def _htildeP(gu, gd, gs):
    invm  = 1 / MUP + 1 / MDN + 1 / MS
    gproj = (gu / MUP + gd / MDN + gs / MS) / invm
    return MN * (DNU / MUP * (gu - gproj)
               + DND / MDN * (gd - gproj)
               + DNS / MS  * (gs - gproj))


def _htildeN(gu, gd, gs):
    """Pseudoscalar NEUTRON coupling: isospin swaps the u and d SPIN CONTENTS
    only -- each 1/m_i stays with its own quark, and gproj (a property of the
    coupling direction, not of the nucleon) is unchanged.  So this is NOT
    _htildeP(gd, gu, gs)."""
    invm  = 1 / MUP + 1 / MDN + 1 / MS
    gproj = (gu / MUP + gd / MDN + gs / MS) / invm
    return MN * (DND / MUP * (gu - gproj)
               + DNU / MDN * (gd - gproj)
               + DNS / MS  * (gs - gproj))


def _a_N(gu, gd, gs):
    return gu * DU + gd * DD + gs * DS


def nucleon_kind(interaction):
    if interaction in ("heavyScalar", "lightScalar", "scalar",
                       "heavyScalarPS", "lightScalarPS", "scalarPS"):
        return "scalar"
    if interaction in ("heavyScalarSP", "lightScalarSP", "scalarSP",
                       "heavyScalarPP", "lightScalarPP", "scalarPP"):
        return "pseudo"
    if interaction == "axialVector":
        return "axial"
    return ""


def dm_vertex_kind(interaction):
    """phi-chi-chi (DM-side) vertex, from the FIRST suffix letter -- what the
    self-interaction (bullet) cross section depends on.  Distinct from
    nucleon_kind, which groups by the phi-SM vertex."""
    if interaction in ("heavyScalar", "lightScalar", "scalar",
                       "heavyScalarSP", "lightScalarSP", "scalarSP"):
        return "scalar"
    if interaction in ("heavyScalarPS", "lightScalarPS", "scalarPS",
                       "heavyScalarPP", "lightScalarPP", "scalarPP"):
        return "pseudo"
    if interaction == "axialVector":
        return "axial"
    return ""


# ---------------------------------------------------------------------------
# Geometry helpers (vectorized; mirrors external_loglikes.m)
# ---------------------------------------------------------------------------

def _seg_dist(x, y, ax, ay, bx, by):
    abx, aby = bx - ax, by - ay
    denom = abx * abx + aby * aby
    t = np.zeros_like(x) if denom == 0 else \
        np.clip(((x - ax) * abx + (y - ay) * aby) / denom, 0.0, 1.0)
    return np.hypot(x - ax - t * abx, y - ay - t * aby)


def _dist_to_polyline(x, y, px, py):
    d = np.full(np.asarray(x).shape, np.inf)
    for i in range(len(px) - 1):
        d = np.minimum(d, _seg_dist(x, y, px[i], py[i], px[i+1], py[i+1]))
    return d


def _inpolygon(x, y, px, py):
    x, y = np.asarray(x, float), np.asarray(y, float)
    inside = np.zeros(x.shape, dtype=bool)
    n, j = len(px), len(px) - 1
    for i in range(n):
        xi, yi, xj, yj = px[i], py[i], px[j], py[j]
        with np.errstate(divide="ignore", invalid="ignore"):
            cross = (yi > y) != (yj > y)
            xint  = (xj - xi) * (y - yi) / (yj - yi) + xi
            inside ^= cross & (x < xint)
        j = i
    return inside


def _region_penalty(x, y, poly_x, poly_y, edge_x, edge_y, sigma):
    inside = _inpolygon(x, y, poly_x, poly_y)
    depth  = np.zeros_like(x)
    if inside.any():
        depth[inside] = _dist_to_polyline(x[inside], y[inside], edge_x, edge_y)
    return _softregion(depth, sigma), inside


def _above_curve_penalty(x, y, cx0, cy0, sigma, extend_left=False):
    """One-sided 'excluded above the curve' region.  extend_left continues the
    bound flat at cy0[0] for all x below the data (cx0 sorted ascending), as
    SN1987A does; that flat continuation also sets the soft-edge depth."""
    y_top = np.max(cy0) + 50.0
    if extend_left:
        x_left = -30.0                       # far below any scan floor
        poly_x = np.concatenate([cx0, [cx0[-1], x_left, x_left]])
        poly_y = np.concatenate([cy0, [y_top,   y_top,  cy0[0]]])
        edge_x = np.concatenate([[x_left], cx0])
        edge_y = np.concatenate([[cy0[0]], cy0])
    else:
        poly_x = np.concatenate([cx0, [cx0[-1], cx0[0]]])
        poly_y = np.concatenate([cy0, [y_top,   y_top]])
        edge_x, edge_y = cx0, cy0
    return _region_penalty(x, y, poly_x, poly_y, edge_x, edge_y, sigma)


# ---------------------------------------------------------------------------
# Curve loader
# ---------------------------------------------------------------------------

_CURVE_CACHE = {}


def load_curve(fpath):
    if fpath not in _CURVE_CACHE:
        _CURVE_CACHE[fpath] = np.loadtxt(fpath)
    return _CURVE_CACHE[fpath]


# ---------------------------------------------------------------------------
# Per-constraint evaluators
# ---------------------------------------------------------------------------

def sn1987a(interaction, mphi, gu, gd, gs, sigma, data_dir):
    kind = nucleon_kind(interaction)
    n = mphi.shape
    ll, excl = np.zeros(n), np.zeros(n, dtype=bool)
    if kind == "scalar":
        gN, fname = _hp(gu, gd, gs), "SN1987AScalarlogmphiloggN.dat"
    elif kind == "pseudo":
        gN, fname = _htildeP(gu, gd, gs), "SN1987AlogmphiloggN_pseudo.dat"
    else:
        return ll, excl, False

    with np.errstate(divide="ignore"):
        log_gN = np.log10(np.abs(gN))
    x = np.log10(mphi)
    curve = load_curve(os.path.join(data_dir, fname))
    cx0, cy0 = curve[:, 0], curve[:, 1]

    X_LEFT = -30.0
    poly_x = np.concatenate([cx0, [X_LEFT, X_LEFT]])
    poly_y = np.concatenate([cy0, [cy0[-1], cy0[0]]])
    edge_x = np.concatenate([cx0, [cx0[0]]])
    edge_y = np.concatenate([cy0, [cy0[0]]])

    fin = np.isfinite(log_gN)
    ll_f, ins_f = _region_penalty(x[fin], log_gN[fin], poly_x, poly_y,
                                  edge_x, edge_y, sigma)
    ins_f &= (x[fin] <= np.max(cx0))
    ll_f   = np.where(ins_f, ll_f, 0.0)
    ll[fin], excl[fin] = ll_f, ins_f
    return ll, excl, True


def neutronstar(interaction, mphi, gu, gd, gs, sigma, data_dir):
    n = mphi.shape
    ll, excl = np.zeros(n), np.zeros(n, dtype=bool)
    if nucleon_kind(interaction) != "scalar":
        return ll, excl, False
    with np.errstate(divide="ignore"):
        log_gN = np.log10(np.abs(_hn(gu, gd, gs)))   # neutron-dominated core
    curve = load_curve(os.path.join(data_dir, "NeutronStarLoggNLogmphi.dat"))
    cx0 = curve[:, 0] - 9.0
    cy0 = curve[:, 1]
    fin = np.isfinite(log_gN)
    ll_f, ins_f = _above_curve_penalty(np.log10(mphi[fin]), log_gN[fin],
                                       cx0, cy0, sigma)
    ll[fin], excl[fin] = ll_f, ins_f
    return ll, excl, True


def stellar(interaction, mphi, gu, gd, gs, sigma, data_dir):
    """Combined stellar-cooling bound (TRGB + horizontal branch) from Hardy &
    Lasenby, quoted as alpha_phiNN = g^2/(4 pi) vs mphi assuming equal
    phi-p / phi-n coupling.  First-pass fudge: plug in the average
    g_bar = (hp + hn)/2.  Data: log10(mphi/eV), log10(alpha_phiNN); excluded
    ABOVE (extra stellar-core energy loss), extended flat to the left in mphi
    below the data (as SN1987A).  SS/PS only (scalar quark vertex)."""
    n = mphi.shape
    ll, excl = np.zeros(n), np.zeros(n, dtype=bool)
    if nucleon_kind(interaction) != "scalar":
        return ll, excl, False
    g_bar = 0.5 * (_hp(gu, gd, gs) + _hn(gu, gd, gs))
    with np.errstate(divide="ignore"):
        log_alpha = np.log10(g_bar ** 2 / (4.0 * np.pi))
    curve = load_curve(os.path.join(data_dir, "alphaphiNNvsmphieV.dat"))
    cx0 = curve[:, 0] - 9.0                   # log10(mphi/eV) -> log10(mphi/GeV)
    cy0 = curve[:, 1]
    fin = np.isfinite(log_alpha)
    ll_f, ins_f = _above_curve_penalty(np.log10(mphi[fin]), log_alpha[fin],
                                       cx0, cy0, sigma, extend_left=True)
    ll[fin], excl[fin] = ll_f, ins_f
    return ll, excl, True


def fifthforce(interaction, mphi, gu, gd, gs, sigma, data_dir,
               sigma_pseudo=0.3):
    """Long-range 5th force -- gx-independent (a nucleon-nucleon force), so it
    depends only on (mphi, g).  Excluded ABOVE the envelope in every case.
    axialVector uses the g_A^2-vs-range envelope; SP/PP use the two pseudoscalar
    spin-spin curves (see _fifthforce_pseudo)."""
    n = mphi.shape
    ll, excl = np.zeros(n), np.zeros(n, dtype=bool)
    kind = nucleon_kind(interaction)
    if kind == "pseudo":
        return _fifthforce_pseudo(mphi, gu, gd, gs, sigma_pseudo, data_dir)
    if kind != "axial":
        return ll, excl, False
    with np.errstate(divide="ignore"):
        log_gA2     = np.log10(_a_N(gu, gd, gs) ** 2)
        log_lambda  = np.log10(HBARC_M / mphi)
    curve = load_curve(os.path.join(data_dir, "axialConstr.dat"))
    cx0, cy0 = curve[:, 0], curve[:, 1]
    fin = np.isfinite(log_gA2)
    ll_f, ins_f = _above_curve_penalty(log_lambda[fin], log_gA2[fin],
                                       cx0, cy0, sigma)
    ll[fin], excl[fin] = ll_f, ins_f
    return ll, excl, True


def _fifthforce_pseudo(mphi, gu, gd, gs, sigma, data_dir):
    """Pseudoscalar (SP/PP) 5th force: two independent spin-spin curves in
    x = log10(mphi/GeV), y = log10 of the constrained product --
      PProtonProtonFifthForce.dat  -> g_P^2      (proton-proton)
      PProtonNeutronFifthForce.dat -> |g_P g_N|  (proton-neutron)
    with g_P = htildeP, g_N = htildeN.  Both files are digitized in
    log10(mphi/eV), hence the -9.  Independent experiments, so the penalties
    ADD and a point is excluded if EITHER curve excludes it.  htildeP =
    htildeN = 0 on the singlet direction g ~ [1 1 1], where both are vacuous."""
    n = mphi.shape
    ll, excl = np.zeros(n), np.zeros(n, dtype=bool)
    with np.errstate(divide="ignore"):
        gP = _htildeP(gu, gd, gs)
        gN = _htildeN(gu, gd, gs)
        prods = ((np.log10(gP ** 2),          "PProtonProtonFifthForce.dat"),
                 (np.log10(np.abs(gP * gN)),  "PProtonNeutronFifthForce.dat"))
    x = np.log10(mphi)
    for y, fname in prods:
        curve = load_curve(os.path.join(data_dir, fname))
        cx0, cy0 = curve[:, 0] - 9, curve[:, 1]   # col1 = log10(mphi/eV)
        fin = np.isfinite(y)
        if not fin.any():
            continue
        ll_f, inside_f = _above_curve_penalty(x[fin], y[fin], cx0, cy0, sigma)
        ll[fin]  += ll_f
        excl[fin] |= inside_f
    return ll, excl, True


def _transfer_L(x):
    """log(1+x) - x/(1+x); series sum_{n>=2} (-1)^n (n-1)/n x^n below x = 0.1."""
    x = np.asarray(x, dtype=float)
    xs = np.minimum(x, 0.1)
    ser = sum((-1) ** n * (n - 1) / n * xs ** n for n in range(2, 25))
    with np.errstate(divide="ignore", invalid="ignore"):
        full = np.log1p(x) - x / (1.0 + x)
    return np.where(x < 0.1, ser, full)


def _transfer_P(x):
    """[6(x+1)log(1+x) + (x^2-3x-6)x] / (x^2 (x+1)); series below x = 0.1."""
    x = np.asarray(x, dtype=float)
    xs = np.minimum(x, 0.1)
    ser = 6.0 * sum((-1) ** n * xs ** (n - 2) / (n * (n - 1)) for n in range(4, 27)) / (xs + 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        full = (6.0 * (x + 1.0) * np.log1p(x) + (x ** 2 - 3.0 * x - 6.0) * x) / (x ** 2 * (x + 1.0))
    return np.where(x < 0.1, ser, full)


def _sigma_self_over_m(interaction, mx, mphi, gx, v):
    """sigma_self/m_chi [cm^2/g], Born transfer cross section at relative velocity
    v (beta=v/c).  R = mx v / mphi, x = R^2.  Depends on (mx, mphi, gx) and the
    DM-side vertex.  Mirrors external_loglikes.m sigma_self_over_m.  NaN where not impl.
      scalar (SS/SP): gx^4/(2pi mx^2 v^4) L(x)
      pseudo (PS/PP): gx^4/(32 mx^2) P(x)
      axial  (AA):    6 gx^4/(pi mx^2 v^4) L(x)"""
    kind = dm_vertex_kind(interaction)
    x = (mx * v / mphi) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        if kind == "scalar":
            sigma_T = gx ** 4 / (2.0 * np.pi * mx ** 2 * v ** 4) * _transfer_L(x)
        elif kind == "pseudo":
            sigma_T = gx ** 4 / (32.0 * mx ** 2) * _transfer_P(x)
        elif kind == "axial":
            sigma_T = 6.0 * gx ** 4 / (np.pi * mx ** 2 * v ** 4) * _transfer_L(x)
        else:
            return np.full(np.shape(mx), np.nan)
    return sigma_T * HBARC ** 2 / (mx * GEV_TO_G)


def bullet(interaction, mx, mphi, gx, sigma, v, bound):
    """Bullet-cluster DM self-interaction bound: sigma_self/m_chi < bound [cm^2/g],
    excluded ABOVE.  Returns (loglike, excluded_mask, applicable)."""
    n = np.shape(mx)
    ll, excl = np.zeros(n), np.zeros(n, dtype=bool)
    if dm_vertex_kind(interaction) == "":
        return ll, excl, False
    som = _sigma_self_over_m(interaction, mx, mphi, gx, v)
    fin = np.isfinite(som) & (som > 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        depth = np.log10(som[fin]) - np.log10(bound)   # >0 over the bound
    ll[fin] = _softregion(depth, sigma)
    excl[fin] = depth > 0.0
    return ll, excl, True


def kaon(interaction, mx, mphi, gx, gu, gd, gs, sigma):
    """Rare-kaon invisible-decay bound, K+ -> pi+ chi chi: excluded where the
    ratio r = Gamma(K+->pi+ chi chi)/[K->pi+invisible limit] > 1.  Open only for
    mx < 0.177 GeV (r = NaN otherwise -> not excluded).  Soft penalty in
    log10(r).  Returns (loglike, excluded_mask, applicable).  See
    kaon_invisible.py for the vectorized ratio and its validation."""
    n = np.shape(mx)
    ll, excl = np.zeros(n), np.zeros(n, dtype=bool)
    if kaon_invisible.kind_of(interaction) is None:
        return ll, excl, False
    r = kaon_invisible.kaon_ratio(mx, mphi, gx, gu, gd, gs, interaction)
    fin = np.isfinite(r) & (r > 0.0)           # r<=0 / NaN -> not excluded
    with np.errstate(divide="ignore", invalid="ignore"):
        depth = np.log10(r[fin])               # >0 where r>1 (excluded)
    ll[fin] = _softregion(depth, sigma)
    excl[fin] = r[fin] > 1.0
    return ll, excl, True


def kaon_reuse_check(prev_path, raw_path, group_name, interaction, sigma,
                     npts=None, valid=None):
    """Can the kaon columns of an existing sidecar stand in for recomputing the
    K+ -> pi+ chi chi ratio?  Returns (ok, why).

    The ratio depends only on the raw inputs, the interaction and
    kaon_invisible.py, and excluded_kaon (r > 1) does not depend on sigma.  So
    the previous sidecar is reusable when it is newer than both the raw file and
    kaon_invisible.py and describes the same rows.  A soft target additionally
    needs a soft source: its penalty rescales as (sigma_old/sigma_new)^2, while
    a hard or off source has lost the depth.  Pass npts/valid for the full row
    check; without them (a preview) only the row count is compared.
    POSTPROC_REUSE_KAON=0 disables reuse; =force skips only the kaon_invisible.py
    age test, for an edit known not to touch this interaction's branch."""
    reuse = os.environ.get("POSTPROC_REUSE_KAON", "1")
    if reuse in ("0", "", "false", "False"):
        return False, "POSTPROC_REUSE_KAON=0"
    if kaon_invisible.kind_of(interaction) is None:
        return False, "no kaon implementation for this interaction"
    if not os.path.isfile(prev_path):
        return False, "no previous sidecar"
    t_prev = os.path.getmtime(prev_path)
    if t_prev <= os.path.getmtime(raw_path):
        return False, "previous sidecar older than the raw scan"
    if reuse != "force" and t_prev <= os.path.getmtime(kaon_invisible.__file__):
        return False, "previous sidecar older than kaon_invisible.py"
    try:
        with h5py.File(prev_path, "r") as h:
            pg = h[group_name]
            need = ("LogLike", "LogLike_isvalid", "LogLike_kaon", "excluded_kaon",
                    "input::ext_sigma_kaon")
            missing = [k for k in need if k not in pg]
            if missing:
                return False, f"previous sidecar lacks {', '.join(missing)}"
            if pg["LogLike"].attrs.get("interaction") != interaction:
                return False, "previous sidecar is a different interaction"
            n_prev = pg["excluded_kaon"].shape[0]
            if npts is None:
                with h5py.File(raw_path, "r") as r:
                    npts = r[group_name]["LogLike"].shape[0]
            if n_prev != npts:
                return False, f"previous sidecar has {n_prev} rows, raw has {npts}"
            s_old = float(pg["input::ext_sigma_kaon"][0])
            if sigma > 0 and not s_old > 0:
                return False, f"soft target needs a soft source (previous sigma {s_old:g})"
            if valid is not None and not np.array_equal(
                    np.asarray(pg["LogLike_isvalid"]).astype(bool), np.asarray(valid).astype(bool)):
                return False, "previous sidecar has a different valid mask"
    except OSError as e:
        return False, f"previous sidecar unreadable ({e})"
    return True, f"previous sigma {s_old:g}"


def kaon_from_sidecar(prev_path, group_name, sigma, vb):
    """The kaon (loglike, excluded_mask, applicable) that kaon() would return,
    rebuilt from a sidecar that kaon_reuse_check() accepted, on the rows vb.

    Hard target: -HARD_PENALTY on excluded_kaon, which is exactly kaon()'s set
    (depth = log10 r > 0 <=> r > 1).  Off: zeros.  Soft: the old penalty times
    (sigma_old/sigma_new)^2."""
    with h5py.File(prev_path, "r") as h:
        pg = h[group_name]
        excl = np.asarray(pg["excluded_kaon"])[vb].astype(bool)
        s_old = float(pg["input::ext_sigma_kaon"][0])
        if sigma < 0:
            ll = np.where(excl, -HARD_PENALTY, 0.0)
        elif sigma == 0:
            ll = np.zeros(excl.shape)
        else:
            ll = np.asarray(pg["LogLike_kaon"], float)[vb]
            if sigma != s_old:
                ll = ll * (s_old / sigma) ** 2
    return ll, excl, True


def _mcmb(ratio_fn, interaction, mx, mphi, gx, gu, gd, gs, sigma):
    """Shared evaluator for a mesonandCMB ratio (kphi/cmb/integral): excluded
    where ratio > 1, soft penalty in log10(ratio).  ratio_fn returns NaN where
    the constraint is closed / not applicable.  Returns (loglike, mask, applicable)."""
    n = np.shape(mx)
    ll, excl = np.zeros(n), np.zeros(n, dtype=bool)
    if mesonandCMB.kind_of(interaction) is None:
        return ll, excl, False
    r = ratio_fn(mx, mphi, gx, gu, gd, gs, interaction)
    fin = np.isfinite(r) & (r > 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        depth = np.log10(r[fin])
    ll[fin] = _softregion(depth, sigma)
    excl[fin] = r[fin] > 1.0
    # 'applicable' is False if the ratio is NaN everywhere (e.g. cmb for SxS/SxP,
    # or K->pi phi with a heavy mediator) so the summary print stays informative.
    return ll, excl, bool(fin.any())


def kphi(interaction, mx, mphi, gx, gu, gd, gs, sigma):
    """K+ -> pi+ phi (on-shell mediator) vs the K+->pi+ + invisible limit; open
    only for mphi < ~0.354 GeV; weighted by the chance the phi looks invisible
    (decays to chi chibar, or escapes NA62), which is where gx enters.  See
    mesonandCMB.py."""
    return _mcmb(mesonandCMB.kphi_ratio, interaction, mx, mphi, gx, gu, gd, gs, sigma)


def cmb(interaction, mx, mphi, gx, gu, gd, gs, sigma):
    """DM annihilation <sigma v> vs the Slatyer CMB bound; ~gx^2; not applicable
    to SxS/SxP (velocity-suppressed).  See mesonandCMB.py."""
    return _mcmb(mesonandCMB.cmb_ratio, interaction, mx, mphi, gx, gu, gd, gs, sigma)


def integral(interaction, mx, mphi, gx, gu, gd, gs, sigma):
    """DM annihilation <sigma v> vs the INTEGRAL telescope bound; ~gx^2; all
    models, light DM only.  See mesonandCMB.py."""
    return _mcmb(mesonandCMB.integral_ratio, interaction, mx, mphi, gx, gu, gd, gs, sigma)


def mesoninv(interaction, mx, mphi, gx, gu, gd, gs, sigma):
    """pi0 / eta / eta' -> chi chibar vs the invisible-decay limits (the largest
    of the three ratios); SxP/PxP/AxA only, open for mx < m_j/2; ~gx^2.  See
    mesonandCMB.py."""
    return _mcmb(mesonandCMB.meson_invisible_ratio, interaction, mx, mphi, gx, gu, gd, gs, sigma)


# ---------------------------------------------------------------------------
# Interaction auto-detection
# ---------------------------------------------------------------------------

def detect_interaction(h5, fallback_path):
    try:
        y = h5["metadata"]["YAML"][0]
        if isinstance(y, bytes):
            y = y.decode()
        m = re.search(r"func_name:\s*loop_(\w+)", y)
        if m:
            return m.group(1)
    except Exception:
        pass
    return os.path.basename(fallback_path).split("_")[0]


# ---------------------------------------------------------------------------
# Core: process one raw file
# ---------------------------------------------------------------------------

def sidecar_path(raw_path):
    return os.path.splitext(raw_path)[0] + "_postprocessed.hdf5"


def _reuse_source(raw_path):
    """The sidecar kaon_reuse_check() may reuse -- the one about to be replaced."""
    return sidecar_path(raw_path)


def process_one(raw_path, group_name):
    """Build <tag>_postprocessed.hdf5 for one raw scan.  Raw is read-only."""
    side = sidecar_path(raw_path)
    tmp  = side + ".tmp"

    nb_mu  = float(os.environ.get("NEFF_BBN_MU",    -0.14))
    nb_sig = float(os.environ.get("NEFF_BBN_SIGMA",   0.21))
    nc_mu  = float(os.environ.get("NEFF_CMB_MU",   -0.054))
    nc_sig = float(os.environ.get("NEFF_CMB_SIGMA",   0.17))
    # External constraints: < 0 hard cut (default), > 0 soft width, 0 off.
    sig_sn = float(os.environ.get("EXT_SIGMA_SN",     EXT_SIGMA_DEFAULT))
    sig_ns = float(os.environ.get("EXT_SIGMA_NS",     EXT_SIGMA_DEFAULT))
    sig_st = float(os.environ.get("EXT_SIGMA_STELLAR", EXT_SIGMA_DEFAULT))  # soft: log10(alpha) width
    sig_ax = float(os.environ.get("EXT_SIGMA_AXIAL",  EXT_SIGMA_DEFAULT))
    sig_fp = float(os.environ.get("EXT_SIGMA_FF_PSEUDO", EXT_SIGMA_DEFAULT))
    sig_bc = float(os.environ.get("EXT_SIGMA_BULLET", EXT_SIGMA_DEFAULT))
    v_bc   = float(os.environ.get("EXT_V_BULLET",     0.01))   # beta = v/c
    bound_bc = float(os.environ.get("EXT_BOUND_BULLET", 1.0))  # cm^2/g
    sig_ka = float(os.environ.get("EXT_SIGMA_KAON",   EXT_SIGMA_DEFAULT))    # soft: width in log10(r)
    sig_kp = float(os.environ.get("EXT_SIGMA_KPHI",     EXT_SIGMA_DEFAULT))
    sig_cm = float(os.environ.get("EXT_SIGMA_CMB",      EXT_SIGMA_DEFAULT))
    sig_in = float(os.environ.get("EXT_SIGMA_INTEGRAL", EXT_SIGMA_DEFAULT))
    sig_mi = float(os.environ.get("EXT_SIGMA_MESONINV", EXT_SIGMA_DEFAULT))

    print(f"Reading (read-only): {raw_path}  (group '/{group_name}')")
    t_read = time.time()
    with h5py.File(raw_path, "r", **_RDCC) as src:
        g = src[group_name]

        interaction = (os.environ.get("EXT_INTERACTION")
                       or detect_interaction(src, raw_path))
        kind = nucleon_kind(interaction)
        print(f"interaction : {interaction}  (nucleon coupling: {kind or 'none'})"
              + ("" if os.environ.get("EXT_INTERACTION") else "  [auto-detected]"))

        Oh2   = np.asarray(g["output::Oh2"])
        dbbn  = np.asarray(g["output::dNeff_bbn"])
        dcmb  = np.asarray(g["output::dNeff_cmb"])
        mphi  = np.asarray(g["input::mphi"])
        mx    = np.asarray(g["input::mx"])
        gx    = np.asarray(g["input::gx"])
        gu    = np.asarray(g["input::gu"])
        gd    = np.asarray(g["input::gd"])
        gs    = np.asarray(g["input::gs"])
        LL    = np.asarray(g["LogLike"])
        valid = np.asarray(g["LogLike_isvalid"])
        vb    = valid.astype(bool)
        npts  = LL.shape[0]

        relic_mode  = int(round(const_input(g, "relic_mode")))
        stored_relic_mu    = const_input(g, "relic_mu")
        stored_relic_sigma = const_input(g, "relic_sigma")
        # The APPLIED (physical) relic cut can be decoupled from the scan's target
        # width.  Widened-sigma scans deliberately broaden the relic Gaussian so
        # Diver's DE covers the ridge boundaries densely; the true Planck cut
        # (0.12 +/- 0.012) is reapplied here via POSTPROC_RELIC_SIGMA/_MU.  Default
        # = the scan's stored value (no decoupling), so unwidened runs are
        # unaffected.
        relic_mu    = float(os.environ.get("POSTPROC_RELIC_MU",    stored_relic_mu))
        relic_sigma = float(os.environ.get("POSTPROC_RELIC_SIGMA", stored_relic_sigma))
        if relic_mu != stored_relic_mu or relic_sigma != stored_relic_sigma:
            print(f"  relic cut OVERRIDE: applying mu={relic_mu}, sigma={relic_sigma} "
                  f"(scan stored mu={stored_relic_mu}, sigma={stored_relic_sigma})")

        # Stored Neff params (relic-only scans have sigma=0 -> disabled).
        sb_mode = int(round(const_input(g, "neff_bbn_mode",   0)))
        sb_mu   = const_input(g, "neff_bbn_mu",  -0.14)
        sb_sig  = const_input(g, "neff_bbn_sigma", 0.0)
        sc_mode = int(round(const_input(g, "neff_cmb_mode",   0)))
        sc_mu   = const_input(g, "neff_cmb_mu",  -0.054)
        sc_sig  = const_input(g, "neff_cmb_sigma", 0.0)

        relic_LL = relic_loglike(Oh2, relic_mode, relic_mu, relic_sigma)
        print(f"Valid points: {vb.sum():,} of {npts:,}  (read {time.time()-t_read:.1f}s)")

        # Self-check: relic + stored_Neff should equal the raw LogLike exactly.
        # This validates the RAW chain, so it must use the scan's STORED relic
        # params -- not any POSTPROC_RELIC_* override applied above.
        relic_LL_stored = (relic_LL
                           if (relic_mu == stored_relic_mu
                               and relic_sigma == stored_relic_sigma)
                           else relic_loglike(Oh2, relic_mode,
                                              stored_relic_mu, stored_relic_sigma))
        recon = (relic_LL_stored
                 + neff_loglike(dbbn, sb_mode, sb_mu, sb_sig)
                 + neff_loglike(dcmb, sc_mode, sc_mu, sc_sig))
        fin   = vb & np.isfinite(LL) & np.isfinite(recon)
        resid = np.nanmax(np.abs(recon[fin] - LL[fin])) if fin.any() else np.nan
        print(f"Reconstruction check  max|relic + stored-Neff - raw LogLike| = {resid:.3g}")
        if not (resid < 1e-6):
            print("  WARNING: raw LogLike != relic + neff(stored params); the raw")
            print("  chain may carry an extra likelihood term -- inspect first.")

        # --- Neff ---
        print(f"Applying Neff: BBN two-sided mu={nb_mu} sigma={nb_sig}; "
              f"CMB two-sided mu={nc_mu} sigma={nc_sig}")
        ll_bbn = neff_loglike(dbbn, 0, nb_mu, nb_sig)
        ll_cmb = neff_loglike(dcmb, 0, nc_mu, nc_sig)
        print(f"  BBN penalized (<-0.5): {np.mean(ll_bbn[vb] < -0.5):.3%}; "
              f"CMB penalized (<-0.5): {np.mean(ll_cmb[vb] < -0.5):.3%}")

        # --- external constraints (valid points only: invalid rows carry sentinel
        #     inputs that would produce NaN warnings and are masked out anyway) ---
        ll_sn = np.zeros(npts); ex_sn = np.zeros(npts, bool)
        ll_ns = np.zeros(npts); ex_ns = np.zeros(npts, bool)
        ll_st = np.zeros(npts); ex_st = np.zeros(npts, bool)
        ll_ff = np.zeros(npts); ex_ff = np.zeros(npts, bool)
        ll_bc = np.zeros(npts); ex_bc = np.zeros(npts, bool)
        ll_ka = np.zeros(npts); ex_ka = np.zeros(npts, bool)
        ll_kp = np.zeros(npts); ex_kp = np.zeros(npts, bool)
        ll_cm = np.zeros(npts); ex_cm = np.zeros(npts, bool)
        ll_in = np.zeros(npts); ex_in = np.zeros(npts, bool)
        ll_mi = np.zeros(npts); ex_mi = np.zeros(npts, bool)
        if vb.any():
            a = (interaction, mphi[vb], gu[vb], gd[vb], gs[vb])
            s, e, ap_sn = sn1987a(*a,     sig_sn, _DATA_DIR)
            ll_sn[vb], ex_sn[vb] = s, e
            s, e, ap_ns = neutronstar(*a, sig_ns, _DATA_DIR)
            ll_ns[vb], ex_ns[vb] = s, e
            s, e, ap_st = stellar(*a, sig_st, _DATA_DIR)
            ll_st[vb], ex_st[vb] = s, e
            s, e, ap_ff = fifthforce(*a,  sig_ax, _DATA_DIR, sig_fp)
            ll_ff[vb], ex_ff[vb] = s, e
            s, e, ap_bc = bullet(interaction, mx[vb], mphi[vb], gx[vb],
                                 sig_bc, v_bc, bound_bc)
            ll_bc[vb], ex_bc[vb] = s, e
            # The kaon ratio is most of the runtime and ignores every likelihood
            # knob, so take it from the previous sidecar when that is still valid.
            prev = _reuse_source(raw_path)
            reuse, why = kaon_reuse_check(prev, raw_path, group_name, interaction,
                                          sig_ka, npts=npts, valid=valid)
            t_ka = time.time()
            if reuse:
                s, e, ap_ka = kaon_from_sidecar(prev, group_name, sig_ka, vb)
                print(f"kaon        : REUSED from {os.path.basename(prev)} ({why})")
            else:
                s, e, ap_ka = kaon(interaction, mx[vb], mphi[vb], gx[vb],
                                   gu[vb], gd[vb], gs[vb], sig_ka)
                print(f"kaon        : computed ({why}; {time.time()-t_ka:.0f}s)")
            ll_ka[vb], ex_ka[vb] = s, e
            ma = (interaction, mx[vb], mphi[vb], gx[vb], gu[vb], gd[vb], gs[vb])
            s, e, ap_kp = kphi(*ma, sig_kp)
            ll_kp[vb], ex_kp[vb] = s, e
            s, e, ap_cm = cmb(*ma, sig_cm)
            ll_cm[vb], ex_cm[vb] = s, e
            s, e, ap_in = integral(*ma, sig_in)
            ll_in[vb], ex_in[vb] = s, e
            s, e, ap_mi = mesoninv(*ma, sig_mi)
            ll_mi[vb], ex_mi[vb] = s, e
        else:
            ap_sn = ap_ns = ap_st = ap_ff = ap_bc = ap_ka = False
            ap_kp = ap_cm = ap_in = ap_mi = False

        ex_any = (ex_sn | ex_ns | ex_st | ex_ff | ex_bc | ex_ka | ex_kp
                  | ex_cm | ex_in | ex_mi)
        print(f"Applicable: sn1987a={ap_sn} ({_smode(sig_sn)}), "
              f"neutronstar={ap_ns} ({_smode(sig_ns)}), "
              f"stellar={ap_st} ({_smode(sig_st)}), "
              f"fifthforce={ap_ff} (axial {_smode(sig_ax)}, pseudo {_smode(sig_fp)}), "
              f"bullet={ap_bc} ({_smode(sig_bc)}, v={v_bc}, bound={bound_bc}), "
              f"kaon={ap_ka} ({_smode(sig_ka)}), kphi={ap_kp} ({_smode(sig_kp)}), "
              f"cmb={ap_cm} ({_smode(sig_cm)}), integral={ap_in} ({_smode(sig_in)}), "
              f"mesoninv={ap_mi} ({_smode(sig_mi)})")
        _hard = [n for n, sg in (("sn1987a", sig_sn), ("neutronstar", sig_ns),
                                 ("stellar", sig_st),
                                 ("fifthforce", sig_ax), ("fifthforce_pseudo", sig_fp),
                                 ("bullet", sig_bc), ("kaon", sig_ka),
                                 ("kphi", sig_kp), ("cmb", sig_cm),
                                 ("integral", sig_in), ("mesoninv", sig_mi)) if sg < 0]
        if _hard:
            print(f"  HARD cut (-{HARD_PENALTY:g} per excluded point): {', '.join(_hard)}")
        if vb.any():
            print(f"  excluded: sn={ex_sn[vb].mean():.3%}, "
                  f"ns={ex_ns[vb].mean():.3%}, stellar={ex_st[vb].mean():.3%}, "
                  f"ff={ex_ff[vb].mean():.3%}, "
                  f"bullet={ex_bc[vb].mean():.3%}, kaon={ex_ka[vb].mean():.3%}, "
                  f"kphi={ex_kp[vb].mean():.3%}, cmb={ex_cm[vb].mean():.3%}, "
                  f"integral={ex_in[vb].mean():.3%}, mesoninv={ex_mi[vb].mean():.3%}, "
                  f"any={ex_any[vb].mean():.3%}")

        # Combined LogLike.
        new_total = (relic_LL + ll_bbn + ll_cmb + ll_sn + ll_ns + ll_st + ll_ff
                     + ll_bc + ll_ka + ll_kp + ll_cm + ll_in + ll_mi)
        LL_new = LL.copy()
        LL_new[vb] = new_total[vb]

        # --- write sidecar atomically -----------------------------------------
        if os.path.exists(tmp):
            os.remove(tmp)
        print(f"Writing sidecar: {side}")
        t_write = time.time()
        try:
            with h5py.File(tmp, "w") as dst:
                dg = dst.create_group(group_name)
                # scan axes + _isvalid companions
                for nm in _AXES:
                    if nm in g:
                        dg.create_dataset(nm, data=g[nm][...])
                        ivn = nm + "_isvalid"
                        dg.create_dataset(ivn,
                                          data=g[ivn][...] if ivn in g else valid)
                # observables
                for nm in _OBS:
                    if nm in g:
                        dg.create_dataset(nm, data=g[nm][...])
                # likelihood columns
                dg.create_dataset("LogLike",              data=LL_new)
                dg.create_dataset("LogLike_isvalid",      data=valid)
                dg.create_dataset("LogLike_orig",         data=LL)
                dg.create_dataset("LogLike_relic",        data=relic_LL)
                dg.create_dataset("LogLike_neff_bbn",     data=ll_bbn)
                dg.create_dataset("LogLike_neff_cmb",     data=ll_cmb)
                dg.create_dataset("LogLike_sn1987a",      data=ll_sn)
                dg.create_dataset("LogLike_neutronstar",  data=ll_ns)
                dg.create_dataset("LogLike_stellar",      data=ll_st)
                dg.create_dataset("LogLike_fifthforce",   data=ll_ff)
                dg.create_dataset("LogLike_bullet",       data=ll_bc)
                dg.create_dataset("LogLike_kaon",         data=ll_ka)
                dg.create_dataset("LogLike_kphi",         data=ll_kp)
                dg.create_dataset("LogLike_cmb",          data=ll_cm)
                dg.create_dataset("LogLike_integral",     data=ll_in)
                dg.create_dataset("LogLike_mesoninv",     data=ll_mi)
                # exclusion masks (uint8)
                dg.create_dataset("excluded_sn1987a",     data=ex_sn.astype(np.uint8))
                dg.create_dataset("excluded_neutronstar", data=ex_ns.astype(np.uint8))
                dg.create_dataset("excluded_stellar",     data=ex_st.astype(np.uint8))
                dg.create_dataset("excluded_fifthforce",  data=ex_ff.astype(np.uint8))
                dg.create_dataset("excluded_bullet",      data=ex_bc.astype(np.uint8))
                dg.create_dataset("excluded_kaon",        data=ex_ka.astype(np.uint8))
                dg.create_dataset("excluded_kphi",        data=ex_kp.astype(np.uint8))
                dg.create_dataset("excluded_cmb",         data=ex_cm.astype(np.uint8))
                dg.create_dataset("excluded_integral",    data=ex_in.astype(np.uint8))
                dg.create_dataset("excluded_mesoninv",    data=ex_mi.astype(np.uint8))
                dg.create_dataset("excluded_any",         data=ex_any.astype(np.uint8))
                dg["LogLike"].attrs["interaction"] = interaction
                # constant params (length-1, plot reads [0])
                params = {
                    "relic_mode": float(relic_mode),
                    "relic_mu": relic_mu, "relic_sigma": relic_sigma,
                    "neff_bbn_mode": 0.0, "neff_bbn_mu": nb_mu,
                    "neff_bbn_sigma": nb_sig,
                    "neff_cmb_mode": 0.0, "neff_cmb_mu": nc_mu,
                    "neff_cmb_sigma": nc_sig,
                    "ext_sigma_sn": sig_sn, "ext_sigma_ns": sig_ns,
                    "ext_sigma_stellar": sig_st,
                    "ext_sigma_axial": sig_ax,
                    "ext_sigma_ff_pseudo": sig_fp,
                    "ext_sigma_bullet": sig_bc, "ext_v_bullet": v_bc,
                    "ext_bound_bullet": bound_bc,
                    "ext_sigma_kaon": sig_ka,
                    "ext_sigma_kphi": sig_kp, "ext_sigma_cmb": sig_cm,
                    "ext_sigma_integral": sig_in,
                    "ext_sigma_mesoninv": sig_mi,
                }
                for k, v in params.items():
                    dg.create_dataset(f"input::{k}", data=np.array([v], np.float64))
            os.replace(tmp, side)
        except BaseException:
            if os.path.exists(tmp):
                os.remove(tmp)
            raise

    print(f"done ({time.time()-t_write:.1f}s).  "
          f"Sidecar carries relic+Neff+SN1987A+NS+stellar+5th-force; raw scan untouched.")
    return {"sidecar": side, "npts": npts, "valid": int(vb.sum())}


# ---------------------------------------------------------------------------
# Batch / CLI
# ---------------------------------------------------------------------------

def resolve_inputs(in_arg):
    """Expand IN to sorted list of raw scan files, skipping all sidecars."""
    if os.path.isdir(in_arg):
        files = globmod.glob(os.path.join(in_arg, "*.hdf5"))
    elif any(c in in_arg for c in "*?[") and not os.path.exists(in_arg):
        files = globmod.glob(in_arg)
    else:
        return [in_arg]
    return sorted(f for f in files
                  if not any(f.endswith(s) for s in _SIDECAR_SUFFIXES))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    in_arg     = sys.argv[1]
    group_name = sys.argv[2] if len(sys.argv) > 2 else "data"
    force      = os.environ.get("POSTPROC_FORCE", "0") not in ("0", "", "false", "False")

    inputs = resolve_inputs(in_arg)
    if not inputs:
        sys.exit(f"no scan files matched '{in_arg}'")

    is_glob = any(c in in_arg for c in "*?[") and not os.path.exists(in_arg)
    if len(inputs) == 1 and not os.path.isdir(in_arg) and not is_glob:
        process_one(inputs[0], group_name)
        return

    print(f"Batch: {len(inputs)} raw file(s)\n")
    done, skipped, failed = [], [], []
    for raw in inputs:
        try:
            with h5py.File(raw, "r") as h:
                if group_name not in h or "LogLike" not in h[group_name]:
                    print(f"-- skip {raw}: no '/{group_name}/LogLike'\n")
                    skipped.append(raw)
                    continue
            if (not force) and os.path.exists(sidecar_path(raw)):
                print(f"-- skip {raw}: sidecar exists (POSTPROC_FORCE=1 to rebuild)\n")
                skipped.append(raw)
                continue
            process_one(raw, group_name)
            done.append(raw)
        except Exception as e:
            print(f"!! FAILED {raw}: {e}\n")
            failed.append(raw)
        else:
            print()

    print("=" * 60)
    print(f"Batch summary: {len(done)} sidecars built, {len(skipped)} skipped, "
          f"{len(failed)} failed (of {len(inputs)}).")
    if failed:
        for f in failed:
            print(f"  FAILED: {f}")
        sys.exit(1)


if __name__ == "__main__":
    main()
