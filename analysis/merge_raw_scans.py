#!/usr/bin/env python3
"""Merge raw GAMBIT scan HDF5 files for the SAME physics config into one file.

The relic plots are coverage / profile-likelihood projections, so the plotted
fill is just the UNION of evaluated points.  Concatenating the per-point columns
of several scanners (Diver + t-walk + paraprof) for the same model therefore
gives a strictly denser plot -- each point keeps its own LogLike / *_isvalid, so
there is no double-counting.  Feed the merged file to postprocess.py exactly like
a single run.

Why this is a tool and not a one-liner -- the guards, which matter because scans
fail and configs differ:

  * SKIP empty / unreadable / LogLike-less inputs (failed or interrupted scans)
    with a warning, instead of aborting the whole merge.
  * REFUSE to mix different configs: the set of VARIED (scanned) parameters must
    match across inputs, so you can't accidentally merge an mphi-scanned (full)
    run with an mphi-fixed (mphi10) run.
  * Concatenate the INTERSECTION of per-point columns; warn on anything dropped.

Usage:
    merge_raw_scans.py OUT.hdf5 IN1.hdf5 IN2.hdf5 [...] [--group data]

The output basename should keep the interaction as its first "_"-token (e.g.
axialVector_relic_merged.hdf5) so postprocess.py's interaction auto-detect works.
"""
import argparse
import sys

import numpy as np
import h5py


# Scanned params clipped to their prior range at merge time.  Out-of-box rows are
# boundary-overshoot evaluations from paraprof's volume-sampling stretch move,
# which proposes theta[j] + z*(theta[k]-theta[j]) with no clip/reflect and
# evaluates it before the accept/reject test (master.py:437-438), so GAMBIT prints
# the out-of-box point.  Diver reflects into the cube and t-walk skips out-of-cube
# proposals, so those scanners never leak.
#
# The bounds below are only a FALLBACK: the real ones are read per file from its
# own metadata (see _prior_ranges).  They have to be, because the ranges differ
# per config -- gx floors are 1e-13 or 1e-16, mphi is [1e-9,10] or [1e-3,10] or
# fixed 10.0 -- so any single project-wide constant necessarily over-clips some
# configs and under-clips others.  These values are the widest range each param
# takes across the yamls, i.e. the most permissive safe default.
PRIOR_RANGES = {
    "mx":   (1.0e-6, 1.0),
    "gx":   (1.0e-16, 1.0),
    "mphi": (1.0e-9, 10.0),
    "gu":   (-1.0, 1.0),
    "gd":   (-1.0, 1.0),
    "gs":   (-1.0, 1.0),
}


def _meta_str(d):
    """Decode a length-1 metadata string dataset to str."""
    v = d[...]
    if getattr(v, "shape", ()):
        v = v[0]
    return v.decode() if isinstance(v, bytes) else str(v)


def _parse_bounds(s):
    """'[1.0e-9, 10.0]' -> (1e-9, 10.0).  Also handles double_log_flat_join's
    4-entry 'ranges: [-1.0, -1e-11, 1e-11, 1.0]' -> (-1.0, 1.0): the outer pair is
    the box, the inner pair is just the flat bridge across zero.  None if
    unparseable."""
    try:
        vals = [float(x) for x in s.strip().strip("[]").split(",")]
    except (ValueError, AttributeError):
        return None
    if len(vals) < 2:
        return None
    return (min(vals), max(vals))


def _prior_ranges(h):
    """Prior bounds {param: (lo, hi)} for one raw scan file, read from its own
    metadata (Parameters::UserModel::pN::{name,range,ranges}), falling back to the
    PRIOR_RANGES constants for anything the metadata doesn't declare.

    fixed_value params are skipped: they aren't scanned dims, so no scanner can
    overshoot them, and comparing a printed double against a re-parsed decimal
    string only risks false drops.
    """
    out = dict(PRIOR_RANGES)
    meta = h.get("metadata")
    if meta is None:
        return out
    params = {}
    for k in meta.keys():
        bits = k.split("::")
        # Parameters::UserModel::p3::range -> params["p3"]["range"]
        if len(bits) == 4 and bits[0] == "Parameters":
            try:
                params.setdefault(bits[2], {})[bits[3]] = _meta_str(meta[k])
            except Exception:
                continue
    for p in params.values():
        name = p.get("name")
        if name not in PRIOR_RANGES or "fixed_value" in p:
            continue
        b = _parse_bounds(p.get("ranges") or p.get("range") or "")
        if b:
            out[name] = b
    return out


def _scan_inputs(g, nrows, valid, ranges):
    """ONE pass over the input:: columns, returning (varied, keep, {param: n_dropped}).

    Reads each column exactly once because these reads dominate the merge
    wallclock: GAMBIT's printer stores datasets with chunks=(100,), i.e. 800-byte
    chunks, so a full column read is ~nrows/100 separate round trips to a network
    FS.  Deriving the varied-param signature and the prior-range keep mask from
    the same read (rather than a column pass each) cuts that traffic ~3x.

      * varied: params that actually vary on valid rows -- the scanned dims, used
        as the config signature so we don't merge unlike scans.
      * keep:   drop rows whose scanned input sits outside its declared prior
        range (boundary-overshoot artifacts).  Non-finite inputs are kept
        (can't judge).  Empty `ranges` (--no-clip) keeps everything.
    """
    varied = set()
    keep = np.ones(nrows, bool)
    dropped = {}
    for k in g.keys():
        if not k.startswith("input::") or k.endswith("_isvalid"):
            continue
        p = k.split("::", 1)[1]
        v = np.asarray(g[k])
        a = v[valid]
        if a.size and np.unique(a[np.isfinite(a)]).size > 3:
            varied.add(p)
        if p in ranges:
            lo, hi = ranges[p]
            fin = np.isfinite(v)
            # ABSOLUTE tolerances, so the box widens outward at both ends and an
            # exact boundary sample survives.  A relative lo*(1-1e-9) would move a
            # NEGATIVE floor *inward* and drop the legitimate gu = -1.0 edge.
            bad = fin & ((v < lo - abs(lo) * 1e-9) | (v > hi + abs(hi) * 1e-9))
            n = int(bad.sum())
            if n:
                dropped[p] = n
            keep &= ~bad
    return varied, keep, dropped


def _per_point_cols(g, nrows):
    """Datasets under the group that are per-point (first axis == nrows)."""
    cols = set()
    for k in g.keys():
        d = g[k]
        if isinstance(d, h5py.Dataset) and d.ndim >= 1 and d.shape[0] == nrows:
            cols.add(k)
    return cols


def inspect(path, group, clip=True):
    """Return per-file metadata for a healthy input, or None if unusable.

    Carries the prior-range keep MASK (a bool array, ~1 byte/row) rather than the
    input columns themselves, so merge() never has to re-read them.
    """
    try:
        with h5py.File(path, "r") as h:
            if group not in h:
                print(f"  SKIP {path}: no group '{group}'")
                return None
            g = h[group]
            if "LogLike" not in g or "LogLike_isvalid" not in g:
                print(f"  SKIP {path}: no LogLike/LogLike_isvalid (failed scan?)")
                return None
            nrows = g["LogLike"].shape[0]
            if nrows == 0:
                print(f"  SKIP {path}: 0 rows (empty / never flushed)")
                return None
            valid = np.asarray(g["LogLike_isvalid"]).astype(bool)
            cols = _per_point_cols(g, nrows)
            ranges = _prior_ranges(h) if clip else {}
            varied, keep, dropped = _scan_inputs(g, nrows, valid, ranges)
            print(f"  OK   {path}: rows={nrows:,} valid={int(valid.sum()):,} "
                  f"cols={len(cols)} varied={sorted(varied)}")
            return dict(path=path, nrows=nrows, cols=cols, varied=varied,
                        keep=keep, dropped=dropped)
    except Exception as e:  # corrupt / unreadable -> skip, don't abort
        print(f"  SKIP {path}: {type(e).__name__}: {str(e)[:70]}")
        return None


def merge(out_path, in_paths, group="data", dry_run=False, clip=True):
    print(f"Inspecting {len(in_paths)} input(s):")
    good = [m for m in (inspect(p, group, clip=clip) for p in in_paths) if m is not None]
    if not good:
        sys.exit("No usable inputs -- nothing to merge.")

    ref = good[0]
    accepted = [ref]
    for m in good[1:]:
        if m["varied"] != ref["varied"]:
            print(f"  REFUSE {m['path']}: varied params {sorted(m['varied'])} "
                  f"!= reference {sorted(ref['varied'])} (different config)")
            continue
        accepted.append(m)

    cols = set.intersection(*(m["cols"] for m in accepted))
    missing = set.union(*(m["cols"] for m in accepted)) - cols
    if missing:
        print(f"  note: dropping {len(missing)} column(s) not present in all "
              f"accepted inputs: {sorted(missing)[:8]}{' ...' if len(missing)>8 else ''}")

    # Masks came from inspect()'s single pass -- see _scan_inputs on why the input
    # columns must not be read a second time here.
    masks = [m["keep"] for m in accepted]
    for m in accepted:
        nout = int(m["keep"].sum())
        if nout != m["nrows"]:
            per = ", ".join(f"{p}:{n:,}" for p, n in sorted(m["dropped"].items()))
            print(f"  clip {m['path']}: dropped {m['nrows'] - nout:,} "
                  f"of {m['nrows']:,} rows outside prior range [{per}]")

    files = [h5py.File(m["path"], "r") for m in accepted]
    try:
        total = int(sum(int(msk.sum()) for msk in masks))
        print(f"Merging {len(accepted)} run(s) -> {total:,} rows, {len(cols)} columns")
        if dry_run:
            print(f"[dry-run] would write {out_path} "
                  f"({', '.join(m['path'] for m in accepted)})")
            return

        with h5py.File(out_path, "w") as ho:
            og = ho.create_group(group)
            for c in sorted(cols):
                parts = [np.asarray(f[group][c])[msk]
                         for f, msk in zip(files, masks)]
                arr = np.concatenate(parts, axis=0)
                og.create_dataset(c, data=arr)
                # carry through the interaction attr etc. from the reference
                for ak, av in files[0][group][c].attrs.items():
                    og[c].attrs[ak] = av
            # copy the metadata group from the reference (YAML, interaction, ...)
            if "metadata" in files[0]:
                files[0].copy("metadata", ho)
    finally:
        for f in files:
            f.close()
    print(f"Wrote {out_path}  ({total:,} rows).  postprocess it like a single run.")


def _selftest():
    """Cheap correctness check: build two tiny matching-schema files, merge,
    assert the row count adds up and values are preserved in order."""
    import os, tempfile
    d = tempfile.mkdtemp()

    # _parse_bounds: 2-entry range, and double_log_flat_join's 4-entry ranges
    # (outer pair is the box; the inner pair is only the zero-bridge).
    assert _parse_bounds("[1.0e-9, 10.0]") == (1e-9, 10.0)
    assert _parse_bounds("[-1.0, -1.0e-11, 1.0e-11, 1.0]") == (-1.0, 1.0)
    assert _parse_bounds("garbage") is None

    paths = []
    for i, n in enumerate((6, 5)):
        p = os.path.join(d, f"scalarPP_relic_part{i}.hdf5")
        with h5py.File(p, "w") as h:
            g = h.create_group("data")
            mx = np.linspace(0.1, 1.0, n)                  # varied, all in range
            mphi = np.logspace(-8, 0, n)                   # varied, all in range
            gu = np.linspace(-1.0, 1.0, n)                 # varied; hits BOTH edges
            if i == 1:
                mx[0] = 9.9e-7                             # mx overshoot   -> drop
                mphi[1] = 6.0e8                            # mphi overshoot -> drop
                gu[2] = -1.1e30                            # gu overshoot   -> drop
            g["LogLike"] = np.full(n, -float(i))
            g["LogLike_isvalid"] = np.ones(n, "int32")
            g["input::mx"] = mx
            g["input::gx"] = np.logspace(-15, -2, n)       # varied; below the 1e-13
            g["input::mphi"] = mphi                        #   fallback floor, but the
            g["input::gu"] = gu                            #   metadata says 1e-16
            g["output::Oh2"] = np.full(n, 0.12) + i
            # Per-config priors, exactly as GAMBIT's hdf5 printer records them.
            m = h.create_group("metadata")
            for pid, (nm, key, val) in enumerate((
                    ("mx",   "range",  "[1.0e-6, 1.0]"),
                    ("mphi", "range",  "[1.0e-9, 10.0]"),
                    ("gx",   "range",  "[1.0e-16, 1.0]"),
                    ("gu",   "ranges", "[-1.0, -1.0e-16, 1.0e-16, 1.0]")), start=1):
                m[f"Parameters::UserModel::p{pid}::name"] = np.array([nm.encode()])
                m[f"Parameters::UserModel::p{pid}::{key}"] = np.array([val.encode()])
        paths.append(p)

    # metadata beats the fallback constants, per param
    with h5py.File(paths[0], "r") as h:
        r = _prior_ranges(h)
    assert r["gx"] == (1e-16, 1.0), r["gx"]        # from metadata, not the fallback
    assert r["gu"] == (-1.0, 1.0), r["gu"]         # 4-entry double_log_flat_join
    assert r["mphi"] == (1e-9, 10.0), r["mphi"]
    assert r["gd"] == PRIOR_RANGES["gd"], r["gd"]  # undeclared -> fallback

    out = os.path.join(d, "scalarPP_relic_merged.hdf5")
    merge(out, paths)
    with h5py.File(out, "r") as h:
        g = h["data"]
        # 6 + 5 rows, minus the three overshoot rows (mx, mphi, gu) in part1
        assert g["LogLike"].shape[0] == 8, g["LogLike"].shape
        assert g["input::mx"][...].min() >= 1e-6, g["input::mx"][...].min()
        assert g["input::mphi"][...].max() <= 10.0, g["input::mphi"][...].max()
        assert g["input::gu"][...].min() >= -1.0, g["input::gu"][...].min()
        # exact +-1.0 gu edges survive (absolute, outward tolerance)
        assert np.isclose(np.abs(g["input::gu"][...]).max(), 1.0)
        # gx down to 1e-15 is legal HERE (metadata floor 1e-16); the old hardcoded
        # 1e-13 constant would have wrongly eaten it
        assert g["input::gx"][...].min() < 1e-13, g["input::gx"][...].min()
        assert np.isclose(g["output::Oh2"][:6], 0.12).all()
        assert np.isclose(g["output::Oh2"][6:], 1.12).all()

    # --no-clip must reproduce the pre-fix row count (keeps the overshoot rows).
    out_nc = os.path.join(d, "scalarPP_relic_merged_noclip.hdf5")
    merge(out_nc, paths, clip=False)
    with h5py.File(out_nc, "r") as h:
        assert h["data"]["LogLike"].shape[0] == 11, h["data"]["LogLike"].shape
    print("SELFTEST OK")


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        _selftest()
        sys.exit(0)
    ap = argparse.ArgumentParser(description="Merge same-config raw GAMBIT scans.")
    ap.add_argument("out")
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--group", default="data")
    ap.add_argument("--dry-run", action="store_true",
                    help="inspect + apply guards + report the plan; write nothing")
    ap.add_argument("--no-clip", action="store_true",
                    help="keep rows outside the prior range (reproduce old behaviour)")
    a = ap.parse_args()
    merge(a.out, a.inputs, a.group, dry_run=a.dry_run, clip=not a.no_clip)
