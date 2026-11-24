#!/usr/bin/env python3
"""
Shared loader for the Loophole constraint sidecar.

A relic-only scan is mapped by Diver, then postprocess.py applies all offline
likelihoods in one pass into a single self-contained sidecar ``<base>_postprocessed.hdf5`` carrying

    LogLike_relic, LogLike_neff_bbn, LogLike_neff_cmb,
    LogLike_sn1987a, LogLike_neutronstar, LogLike_stellar, LogLike_fifthforce,
    LogLike (= their sum), LogLike_isvalid, the scan axes + output::sigma_chiN,
    output::Oh2/dNeff_bbn/dNeff_cmb, and length-1 input:: parameter columns.

This module locates that sidecar from the raw scan or the sidecar path and
returns the relic + per-constraint
likelihoods, axes, valid mask and interaction.  All plotters go through
``load`` so they cannot drift in how they read the columns.
"""

import os
import re

import numpy as np
import h5py

import scan_axes

# Canonical constraint key -> sidecar dataset name.  Order = overlay order.
CONSTRAINT_COLS = {
    "neff_bbn":    "LogLike_neff_bbn",
    "neff_cmb":    "LogLike_neff_cmb",
    "sn1987a":     "LogLike_sn1987a",
    "neutronstar": "LogLike_neutronstar",
    "stellar":     "LogLike_stellar",
    "fifthforce":  "LogLike_fifthforce",
    "bullet":      "LogLike_bullet",
    "kaon":        "LogLike_kaon",
    "kphi":        "LogLike_kphi",
    "cmb":         "LogLike_cmb",
    "integral":    "LogLike_integral",
    "mesoninv":    "LogLike_mesoninv",
}

# Constraints the plots show as ONE region: key -> members, folded into the
# first member's overlay slot.  The sidecar keeps every member's column; summing
# the hard-cut LogLikes leaves total_LL unchanged and penalizes a point iff any
# member does, so the merged region is the union of the members' exclusions.
MERGED = {
    "meson": ("kaon", "kphi", "mesoninv"),   # K->pi chi chi, K->pi phi, pi0/eta/eta'->inv
}


def _merge(comp):
    """Fold each MERGED group into a single entry, keeping overlay order."""
    for key, members in MERGED.items():
        merged = sum(comp[m] for m in members)
        out = {}
        for name, c in comp.items():
            if name == members[0]:
                out[key] = merged
            elif name not in members:
                out[name] = c
        comp = out
    return comp


def _strip_base(in_arg):
    """Run base path: drop a trailing _postprocessed + .hdf5."""
    return re.sub(r"(_postprocessed)?\.hdf5$", "", in_arg)


def locate(in_arg):
    """Find the constraint sidecar for a run.

    Returns ``(base, kind, paths)`` where kind is 'postprocessed' or None
    (nothing found), and ``paths`` is {'pp': path} or {}.
    """
    base = _strip_base(in_arg)
    pp = base + "_postprocessed.hdf5"
    if os.path.exists(pp):
        return base, "postprocessed", {"pp": pp}
    return base, None, {}


class Constraints:
    """Container for the loaded sidecar arrays (all full-length, unfiltered)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    @property
    def active(self):
        """Constraint keys that actually penalize some valid point."""
        return [n for n, c in self.comp.items()
                if np.any(c[self.valid] < 0.0)]


def load(in_arg, group_name="data"):
    """Load relic + per-constraint likelihoods, axes, valid mask, interaction.

    Returns a :class:`Constraints` with attributes:

        avail        {axis_key: (dataset, mode, label)} present in the sidecar
        axis_raw     {axis_key: ndarray}  full-length (unfiltered)
        valid        bool ndarray         validity & finiteness of relic+total
        relic_LL     ndarray              LogLike_relic
        comp         {name: ndarray}      per-constraint LogLike (0 if absent),
                                          MERGED groups folded into one key
        total_LL     ndarray              relic + sum(comp)
        interaction  str
        source       str (for logging)

    Raises FileNotFoundError if no sidecar is found.
    """
    base, kind, paths = locate(in_arg)
    if kind is None:
        raise FileNotFoundError(
            f"no constraint sidecar for base {base!r}: expected "
            f"{base}_postprocessed.hdf5.  Build it with postprocess.py.")

    comp = {}
    with h5py.File(paths["pp"], "r") as h:
        g = h[group_name]
        relic_LL = g["LogLike_relic"][...]
        valid = g["LogLike_isvalid"][...].astype(bool)
        for name, dset in CONSTRAINT_COLS.items():
            comp[name] = g[dset][...] if dset in g else np.zeros_like(relic_LL)
        avail = scan_axes.available(g.keys())
        axis_raw = {k: g[avail[k][0]][...] for k in avail}
        interaction = g["LogLike"].attrs.get(
            "interaction", os.path.basename(base).split("_")[0])
    source = paths["pp"]

    comp = _merge(comp)
    total_LL = relic_LL + sum(comp.values())
    valid = valid & np.isfinite(relic_LL) & np.isfinite(total_LL)

    return Constraints(base=base, kind=kind, avail=avail, axis_raw=axis_raw,
                       valid=valid, relic_LL=relic_LL, comp=comp,
                       total_LL=total_LL, interaction=interaction, source=source)
