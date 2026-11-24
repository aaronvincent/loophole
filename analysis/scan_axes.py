"""Shared axis registry + auto-windowing for the Loophole scan plotters.

Single source of truth (used by both ``plot_2D_profiles.py`` and
``plot_allowed_regions.py``) for:

  * which scan quantities are plottable and how each maps to a plot axis
    (``log`` = log10 for positive quantities; ``symlog`` = signed log10 for the
    quark couplings, whose ``double_log_flat`` prior changes sign);
  * pretty axis labels;
  * **data-driven axis bounds** (so a frame is never clipped to a hardcoded
    window that does not match the scan -- the eta-scan cut-off bug);
  * which input parameters *actually varied* in a given scan, with
    perfectly-correlated (``same_as``-tied) axes collapsed -- so eta / gu / gd /
    gs planes appear only for the scans that move them.

Each scan therefore gets an axis for every varied model parameter automatically.
"""

import numpy as np

# Linear threshold for the signed-log transform (below the gu/gd/gs prior gap of
# 1e-11, so the whole sampled range is in the logarithmic part of the map).
SYMLOG_LIN = 1.0e-12

# key -> (hdf5 dataset, transform mode, pretty axis label)
KNOWN = {
    "mx":         ("input::mx",          "log",    r"$\log_{10}(m_\chi\,/\,\mathrm{GeV})$"),
    "mphi":       ("input::mphi",        "log",    r"$\log_{10}(m_\phi\,/\,\mathrm{GeV})$"),
    "gx":         ("input::gx",          "log",    r"$\log_{10} g_\chi$"),
    "eta":        ("input::eta",         "log",    r"$\log_{10}\eta$"),
    "gu":         ("input::gu",          "symlog", r"signed $\log_{10} g_u$"),
    "gd":         ("input::gd",          "symlog", r"signed $\log_{10} g_d$"),
    "gs":         ("input::gs",          "symlog", r"signed $\log_{10} g_s$"),
    "sigma_chiN": ("output::sigma_chiN", "log",    r"$\log_{10}(\sigma_{\chi N}\,/\,\mathrm{cm}^2)$"),
}

# Preference order for de-duplicating perfectly-correlated axes (a master like
# eta wins over the gu/gd/gs it ties together) and for ordering default planes.
ORDER = ["mx", "mphi", "gx", "eta", "gu", "gd", "gs", "sigma_chiN"]

MASSES = ("mx", "mphi")
COUPLINGS = ("gx", "eta", "gu", "gd", "gs")
QUARK_COUPLINGS = ("gu", "gd", "gs")

# Derived axes: quantities NOT stored as their own dataset, computed on the fly
# from several scanned columns.  Used for the gquark scans, where the quark
# couplings gu/gd/gs are varied independently (no `input::eta`): we want a plot
# axis for the *magnitude* of the coupling vector, eta = |(g_u, g_d, g_s)|.
# key -> (source axis keys, transform mode, pretty axis label)
DERIVED = {
    "gmag": (("gu", "gd", "gs"), "log",
             r"$\log_{10}\eta,\ \ \eta\equiv\sqrt{g_u^2+g_d^2+g_s^2}$"),
}


def derived_available(avail):
    """Derived axes whose every source axis is present in ``avail``."""
    return {k: v for k, v in DERIVED.items()
            if all(s in avail for s in v[0])}


def compute_derived(key, data):
    """Compute a derived quantity array from already-read source datasets."""
    srcs = DERIVED[key][0]
    if key == "gmag":
        # Computed on the FULL (unfiltered) arrays, which still carry the hdf5
        # fill sentinels on invalid rows; squaring those overflows to inf.  Those
        # rows are always dropped by the caller's `valid` mask, so suppress the
        # (harmless) overflow warning rather than let it spam the logs.
        with np.errstate(over="ignore", invalid="ignore"):
            return np.sqrt(sum(np.asarray(data[s], float) ** 2 for s in srcs))
    raise KeyError(f"unknown derived axis: {key}")


def transform(vals, mode):
    """Map a raw quantity to its plot axis. ``symlog`` is monotonic in ``vals``
    and smooth through zero (signed decades from ``SYMLOG_LIN``)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        if mode == "symlog":
            return np.sign(vals) * np.log10(1.0 + np.abs(vals) / SYMLOG_LIN)
        return np.log10(vals)


def available(group_keys):
    """KNOWN entries whose dataset is present in the file, in registry order."""
    keys = set(group_keys)
    return {k: KNOWN[k] for k in ORDER if KNOWN[k][0] in keys}


def auto_bounds(tvals, pad_frac=0.04, round_to=0.5):
    """Data-driven ``[lo, hi]`` in transformed space: a robust percentile span,
    padded slightly and rounded outward to a clean grid.  Replaces the hardcoded
    windows that clipped scans whose data lived elsewhere.

    The lower edge is then clamped to the true data minimum so a hard prior wall
    (e.g. the ``mx`` scan floor at 1e-6) sits flush with the frame instead of
    leaving an empty margin below the scan region.  Because the rounded bound is
    derived from a robust 0.05-percentile, this only ever pulls the edge *in*
    toward real points -- genuine low outliers stay cropped."""
    v = np.asarray(tvals)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return (-1.0, 1.0)
    vmin = float(v.min())
    lo, hi = np.percentile(v, 0.05), np.percentile(v, 99.95)
    if hi <= lo:
        lo, hi = vmin, float(v.max())
    if hi <= lo:
        hi = lo + 1.0
    span = hi - lo
    lo -= pad_frac * span
    hi += pad_frac * span
    lo = np.floor(lo / round_to) * round_to
    hi = np.ceil(hi / round_to) * round_to
    lo = max(lo, vmin)
    return (float(lo), float(hi))


def _is_varied(raw):
    v = raw[np.isfinite(raw)]
    return v.size > 0 and np.unique(v).size > 1


def varied_axes(raw_by_key):
    """Input axes that actually move, with perfectly-correlated ones collapsed
    (keep the registry-preferred member -- e.g. eta over the gu/gd/gs it ties).

    ``sigma_chiN`` is excluded: it is a derived output, handled separately.
    Correlation is computed on the *raw* sampled values, so ``same_as``-tied
    parameters (exactly proportional) give |corr| = 1 while independently
    log-sampled parameters give ~0, regardless of the likelihood surface.
    """
    cand = [k for k in ORDER
            if k in raw_by_key and k != "sigma_chiN" and _is_varied(raw_by_key[k])]
    kept = []
    for k in cand:
        rk = raw_by_key[k]
        dup = False
        for j in kept:
            rj = raw_by_key[j]
            m = np.isfinite(rk) & np.isfinite(rj)
            if m.sum() > 2 and np.std(rk[m]) > 0 and np.std(rj[m]) > 0:
                # A near-constant column can still make corrcoef divide by ~0;
                # the |corr| test below simply won't fire, so ignore the warning.
                with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
                    r = abs(np.corrcoef(rk[m], rj[m])[0, 1])
                if r > 0.99999:
                    dup = True
                    break
        if not dup:
            kept.append(k)
    return kept


def default_pairs(varied, have_sigma):
    """Sensible ``(x, y)`` planes from the varied axes: each coupling vs each
    mass, the quark couplings against each other, the two masses against each
    other, and sigma_chiN vs each mass."""
    masses = [k for k in MASSES if k in varied]
    couplings = [k for k in COUPLINGS if k in varied]
    pairs = []
    for m in masses:
        for c in couplings:
            pairs.append((m, c))
    # Quark-coupling planes.  The gquark scans move gu/gd/gs independently, and
    # the pi0/eta-mixing physics lives in their differences (g_u - g_d, ...), so
    # the mutual planes carry structure no mass plane shows.  In an eta-tied scan
    # they collapse to a single axis and none of these pairs appear.
    quarks = [k for k in QUARK_COUPLINGS if k in varied]
    for i, a in enumerate(quarks):
        for b in quarks[i + 1:]:
            pairs.append((a, b))
    if "mx" in masses and "mphi" in masses:
        pairs.append(("mx", "mphi"))
    if have_sigma:
        for m in masses:
            pairs.append((m, "sigma_chiN"))
    return pairs
