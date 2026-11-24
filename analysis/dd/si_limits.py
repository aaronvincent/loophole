#!/usr/bin/env python3
"""
Spin-independent DM-nucleon direct-detection limits, as an overlay for the
Loophole (m_chi, sigma_chiN) panels.

Data: Ciaran O'Hare's DirectDetectionPlots compilation, fetched into ./limits/
by ``fetch_si_limits.sh``.  Every file is two columns, ``m_chi [GeV]`` and
``sigma_p^SI [cm^2]``; the ``AllLimits-*`` / ``*_Combined`` envelopes pad the
masses they do not cover with a placeholder sigma = 1 cm^2, which is dropped
here (:data:`PLACEHOLDER`).

Two kinds of curve, because they bound sigma in different directions:

``upper``   an ordinary exclusion: everything ABOVE the curve is excluded, over
            the mass interval the experiment covers.  All the underground
            detectors, and the Milky-Way-satellite bound.
``band``    a closed contour: the excluded region is the polygon INTERIOR.  XQC
            is the only one here -- its upper edge is the overburden ceiling,
            above which DM never reaches the detector.

The overlay is built in log10(m), log10(sigma) space, which is what the
plot_all_constraints panels plot in.

CAVEATS (read before believing a figure made with this):

  * Published SI limits assume a CONTACT interaction (heavy mediator) and the
    standard halo.  This scan runs m_phi from 1e-9 to 10 GeV, and at sub-GeV
    m_chi the momentum transfer is q ~ mu*v ~ 0.1 MeV .. 0.1 GeV.  Wherever
    m_phi < q the published bound would have to be re-derived with the light-
    mediator form factor; the curve drawn here is then only indicative.
  * The scan's ``output::sigma_chiN`` is the full angular-integrated per-nucleon
    cross-section from DMnucxsec (full propagator, no contact limit), so it is
    the same OBSERVABLE as sigma_p^SI but not, in general, in the same
    kinematic regime -- see the point above.
  * MW satellites is a structure-formation bound, not direct detection; it is
    the only curve here that covers the whole m_chi range of this model space,
    so it is carried separately and drawn in its own style.
"""

import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LIMITS_DIR = os.environ.get("DD_LIMITS_DIR", os.path.join(HERE, "limits"))

# The AllLimits/_Combined envelopes pad uncovered masses with sigma = 1 cm^2.
PLACEHOLDER = 1.0

# name -> (filename, kind, pretty label, reference)
#
# The reference string is transcribed from the "#" header of the upstream file,
# NOT looked up independently -- several of these curves are digitised from a
# figure, and a couple carry more than one paper (a combined limit).  Check the
# header before putting one in a bibliography:  grep -h '^#' limits/<file>.txt
CATALOG = {
    "XENON1T-Migdal":     ("XENON1T-Migdal.txt",     "upper", "XENON1T (Migdal)",   "1907.12771"),
    "CRESST":             ("CRESST.txt",             "upper", "CRESST-III",         "1904.00498"),
    "DarkSide-50_2022":   ("DarkSide-50_2022.txt",   "upper", "DarkSide-50 (2022)", "2207.11966"),
    "DarkSide":           ("DarkSide.txt",           "upper", "DarkSide-50",        "1802.06994"),
    "CDMSLite":           ("CDMSLite.txt",           "upper", "CDMSlite",           "1808.09098; 0912.3592"),
    "CDEX":               ("CDEX.txt",               "upper", "CDEX-10",            "PRL 120, 241301"),
    "NEWS-G":             ("NEWS-G.txt",             "upper", "NEWS-G",             "1706.04934"),
    "PandaX_LowMass":     ("PandaX_LowMass.txt",     "upper", "PandaX (low mass)",  "2212.10067; 2207.04883"),
    "PandaX_2024_Combined": ("PandaX_2024_Combined.txt", "upper", "PandaX-4T",      "2408.00664; 2212.10067; 2207.04883"),
    "XENONnT_2024_Combined": ("XENONnT_2024_Combined.txt", "upper", "XENONnT",      "2409.17868; 2303.14729"),
    "XENONnT_2025_Combined": ("XENONnT_2025_Combined.txt", "upper", "XENON (2025 envelope)",
                              "none given upstream; header mislabels the column sigma_n_SD"),
    "LZ_2024":            ("LZ_2024.txt",            "upper", "LZ (2024)",          "2410.17036"),
    "XENON1T":            ("XENON1T.txt",            "upper", "XENON1T",            "1805.12562; 2012.02846"),
    "DEAP-3600":          ("DEAP-3600.txt",          "upper", "DEAP-3600",          "1902.04048"),
    "PICO60":             ("PICO60.txt",             "upper", "PICO-60",            "1702.07666; PRD 100, 022001"),
    "EDELWEISS":          ("EDELWEISS.txt",          "upper", "EDELWEISS",          "1607.03367"),
    "XQC":                ("XQC.txt",                "band",  "XQC (overburden)",   "PRD 76, 042007, rescaled by 2010.07240"),
    "MWSatellites":       ("MWSatellites.txt",       "upper", "MW satellites",
                           "1904.10000 (Nadler et al., ApJL 878, 32)"),
    "AllLimits-2022":     ("AllLimits-2022.txt",     "upper", "best SI limit (2022)", "compilation; no header"),
    "AllLimits-2024":     ("AllLimits-2024.txt",     "upper", "best SI limit (2024)", "compilation; no header"),
}

# Default envelope members: every underground detector in the catalogue.  The
# envelope is their pointwise minimum, so adding an experiment can only ever
# strengthen it -- and the ones that matter for m_chi < 0.2 GeV are
# XENON1T-Migdal (from 0.104 GeV) and CRESST-III (from 0.171 GeV).
DEFAULT_ENVELOPE = [
    "XENON1T-Migdal", "CRESST", "DarkSide-50_2022", "DarkSide", "CDMSLite",
    "CDEX", "NEWS-G", "PandaX_LowMass", "PandaX_2024_Combined",
    "XENONnT_2024_Combined", "XENONnT_2025_Combined", "LZ_2024", "XENON1T",
    "DEAP-3600", "PICO60", "EDELWEISS",
]

# Bounds that are NOT underground direct detection but live in the same file
# tree and are the only thing covering m_chi below 0.1 GeV.  Drawn separately.
COSMO_ASTRO = ["MWSatellites"]
BANDS = ["XQC"]


def available():
    """Catalogue entries whose file is present in ``LIMITS_DIR``."""
    return [k for k, v in CATALOG.items()
            if os.path.exists(os.path.join(LIMITS_DIR, v[0]))]


def load(name):
    """``(m [GeV], sigma [cm^2])`` for one catalogue entry, placeholders dropped.

    Rows are returned in file order -- for a ``band`` entry that order IS the
    polygon, so it must not be sorted.
    """
    fn, kind, _label, _ref = CATALOG[name]
    path = os.path.join(LIMITS_DIR, fn)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found -- run ./fetch_si_limits.sh (or fetch_si_limits.sh {name})")
    d = np.loadtxt(path)
    if d.ndim != 2 or d.shape[1] < 2:
        raise ValueError(f"{path}: expected 2 columns, got shape {d.shape}")
    m, s = d[:, 0].astype(float), d[:, 1].astype(float)
    keep = np.isfinite(m) & np.isfinite(s) & (m > 0) & (s > 0) & (s < PLACEHOLDER)
    if kind == "band":
        # A polygon vertex list: dropping interior rows would tear the contour,
        # so only reject outright-invalid points (there are none in practice).
        keep = np.isfinite(m) & np.isfinite(s) & (m > 0) & (s > 0)
    return m[keep], s[keep]


def curve_log(name):
    """``(log10 m, log10 sigma)`` for one entry, ``upper`` curves sorted in m."""
    m, s = load(name)
    if CATALOG[name][1] == "upper":
        o = np.argsort(m)
        m, s = m[o], s[o]
    return np.log10(m), np.log10(s)


def envelope_log(names=None, n=2000, mlim=None):
    """Pointwise best (lowest) limit over ``names``, in log10-log10.

    Returns ``(lm, ls)`` on a uniform log-mass grid; ``ls`` is NaN at masses no
    member covers, so plotting it leaves an honest gap instead of extrapolating
    a bound where none exists.  ``mlim`` is an optional ``(log10 m_lo, log10
    m_hi)`` clip -- pass the panel's x range to avoid gridding decades that are
    then thrown away.
    """
    names = list(names if names is not None else DEFAULT_ENVELOPE)
    names = [n_ for n_ in names if os.path.exists(os.path.join(LIMITS_DIR, CATALOG[n_][0]))]
    if not names:
        raise FileNotFoundError("no envelope members present -- run ./fetch_si_limits.sh")
    curves = [curve_log(n_) for n_ in names]
    lo = min(c[0].min() for c in curves)
    hi = max(c[0].max() for c in curves)
    if mlim is not None:
        lo, hi = max(lo, mlim[0]), min(hi, mlim[1])
    if not (hi > lo):
        return np.array([]), np.array([])
    lm = np.linspace(lo, hi, n)
    best = np.full(n, np.inf)
    for cm, cs in curves:
        # np.interp clamps outside the data range, which would smear an
        # experiment's endpoint across the whole grid -- mask to its support.
        inside = (lm >= cm.min()) & (lm <= cm.max())
        if not inside.any():
            continue
        best[inside] = np.minimum(best[inside], np.interp(lm[inside], cm, cs))
    best[~np.isfinite(best)] = np.nan
    return lm, best


def covered_span(names=None):
    """``(log10 m_lo, log10 m_hi)`` union actually covered by ``names``."""
    names = list(names if names is not None else DEFAULT_ENVELOPE)
    names = [n_ for n_ in names if os.path.exists(os.path.join(LIMITS_DIR, CATALOG[n_][0]))]
    curves = [curve_log(n_) for n_ in names]
    return (min(c[0].min() for c in curves), max(c[0].max() for c in curves))


class PowerOfTenFormatter:
    """``10^n`` tick labels, as a picklable stand-in for the closure that
    ``plot_all_constraints.style_log_axis`` installs.

    A figure cannot be pickled while a local lambda sits on one of its axes, and
    caching the finished panel is what lets the overlay be restyled without
    re-reading a multi-GB sidecar.  It lives here, in a module that is imported
    rather than run, because a class defined in ``__main__`` does not survive
    the ``runpy`` hand-off either.  The labels are byte-identical.
    """

    def __call__(self, v, _pos=None):
        return r"$10^{%d}$" % int(round(v))


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
# Achromatic on purpose.  plot_all_constraints reserves blue for "allowed" and
# gives each constraint family a hue (yellow = stars, teal = lab, vermillion =
# early universe, wine = DM astrophysics); a sixth hue would not survive the
# alpha stacking.  Direct detection is also a different KIND of statement here:
# it is not in the scan likelihood, it is an external curve laid over the
# result.  Black line + open hatch reads on top of every fill and in greyscale.
DD_COLOR = os.environ.get("DD_COLOR", "#000000")
ASTRO_COLOR = os.environ.get("DD_ASTRO_COLOR", "#404040")


def draw_upper(ax, lm, ls, ylim, color=DD_COLOR, label=None, hatch="xx",
               lw=2.2, ls_style="-", alpha=0.75, zorder=5.0, fill="band",
               band_frac=0.055, ceiling=None):
    """Draw an ``upper`` exclusion: the curve, plus the excluded side of it.

    ``fill`` picks how much of the excluded side is drawn:

    ``"band"``  (default) hatch a strip of ``band_frac`` x the y span just above
                the curve -- the usual "hatched side of the line" convention.
                These panels are ~60 decades tall and already carry eight
                translucent constraint fills, so hatching all the way up buries
                them for no extra information.
    ``"full"``  hatch everything from the curve to ``ylim[1]``.
    ``"none"``  the line only.

    ``ceiling`` (log10 cm^2), with ``fill="band"``, replaces the fixed strip by
    a hatch from the curve up to that sigma.  An underground detector is blind
    above some cross-section -- the DM thermalises in the rock and atmosphere
    before reaching it -- so an SI limit is really a band, not a half-plane.
    The true ceiling is experiment- and mass-dependent; a single value is an
    approximation, which is why no edge line is drawn along it.

    NaN gaps in ``ls`` (masses no experiment covers) are drawn as separate
    pieces, so the overlay stops where the data stop.
    """
    if fill is True:
        fill = "full"
    elif fill is False:
        fill = "none"
    handles = []
    # A curve that runs off the TOP of the panel has nothing to say inside it:
    # keep drawing it and the band clamp below would invert fill_between and
    # hatch the whole panel.  One that runs off the BOTTOM excludes everything
    # visible, so clamp it to the bottom edge instead of dropping it.
    ls = np.maximum(ls, ylim[0])
    good = np.isfinite(ls) & (ls < ylim[1])
    if not good.any():
        return handles
    # Contiguous runs of covered mass.
    pieces, cur = [], []
    for i, g in enumerate(good):
        if g:
            cur.append(i)
        elif cur:
            pieces.append(np.array(cur)); cur = []
    if cur:
        pieces.append(np.array(cur))
    for k, idx in enumerate(pieces):
        if idx.size < 2:
            continue
        x, y = lm[idx], ls[idx]
        if fill != "none":
            if fill == "full":
                top = np.full_like(y, ylim[1])
            elif ceiling is not None:
                top = np.maximum(np.full_like(y, min(ceiling, ylim[1])), y)
            else:
                top = np.minimum(y + band_frac * (ylim[1] - ylim[0]), ylim[1])
            ax.fill_between(x, y, top, facecolor="none", edgecolor=color,
                            hatch=hatch, linewidth=0.0, alpha=alpha, zorder=zorder)
        (line,) = ax.plot(x, y, color=color, lw=lw, ls=ls_style, zorder=zorder + 0.2,
                          solid_capstyle="round",
                          label=label if (k == 0 and label) else None)
        if k == 0 and label:
            handles.append(line)
    return handles


def draw_band(ax, name, color=DD_COLOR, label=None, hatch="..", lw=1.6,
              alpha=0.75, zorder=5.0):
    """Draw a ``band`` entry (XQC): the closed polygon, hatched interior.

    The polygon is drawn whole and clipped by the axes; its vertex order is the
    contour, so :func:`load` must not sort it."""
    m, s = load(name)
    x, y = np.log10(m), np.log10(s)
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    if (x.min() > xlim[1] or x.max() < xlim[0] or
            y.min() > ylim[1] or y.max() < ylim[0]):
        return []                       # entirely off the panel: no line, no legend row
    ax.fill(x, y, facecolor="none", edgecolor=color, hatch=hatch, linewidth=0.0,
            alpha=alpha, zorder=zorder)
    (line,) = ax.plot(np.r_[x, x[:1]], np.r_[y, y[:1]], color=color, lw=lw,
                      ls="--", zorder=zorder + 0.2, label=label)
    return [line] if label else []
