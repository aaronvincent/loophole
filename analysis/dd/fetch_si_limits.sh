#!/bin/bash
# ---------------------------------------------------------------------------
# Fetch the spin-independent DM-nucleon limit curves from Ciaran O'Hare's
# DirectDetectionPlots repository into ./limits/.
#
#   https://github.com/cajohare/DirectDetectionPlots/tree/main/data/limits/Nucleon/SI
#
# Each file is two whitespace-separated columns, "# " comments:
#     m_chi [GeV]     sigma_p^SI [cm^2]
#
# These are DOWNLOADED REFERENCE CURVES for a plot overlay; they live inside
# this prototype folder and nothing in the repo's own data/ runs/ oldruns/ trees
# is touched.  Re-running is idempotent (curl overwrites the cached copy).
#
# Usage:
#   ./fetch_si_limits.sh            # the curated set used by the prototype
#   ./fetch_si_limits.sh --all      # every .txt in the upstream SI folder
#   ./fetch_si_limits.sh NAME ...   # named files, with or without .txt
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/limits"
RAW="https://raw.githubusercontent.com/cajohare/DirectDetectionPlots/main/data/limits/Nucleon/SI"
API="https://api.github.com/repos/cajohare/DirectDetectionPlots/contents/data/limits/Nucleon/SI"

# The curated set: everything that reaches below ~10 GeV (this model space tops
# out at m_chi ~ 0.17 GeV), plus the heavy-mass workhorses for context, plus the
# two bounds that cover the whole sub-GeV range (XQC overburden, MW satellites).
CURATED=(
  XENON1T-Migdal CRESST DarkSide-50_2022 DarkSide CDMSLite CDEX NEWS-G
  PandaX_LowMass PandaX_2024_Combined XENONnT_2024_Combined XENONnT_2025_Combined
  LZ_2024 XENON1T DEAP-3600 PICO60 EDELWEISS
  XQC MWSatellites
  AllLimits-2022 AllLimits-2024
)

mkdir -p "$OUT"

if [[ "${1:-}" == "--all" ]]; then
  mapfile -t FILES < <(curl -fsS "$API" |
    python3 -c 'import json,sys; [print(f["name"]) for f in json.load(sys.stdin) if f["name"].endswith(".txt")]')
elif [[ $# -gt 0 ]]; then
  FILES=("$@")
else
  FILES=("${CURATED[@]}")
fi

n_ok=0
for f in "${FILES[@]}"; do
  f="${f%.txt}.txt"
  if curl -fsS -m 60 -o "$OUT/$f" "$RAW/$f"; then
    printf '  %-34s %s\n' "$f" "$(wc -l < "$OUT/$f") lines"
    n_ok=$((n_ok + 1))
  else
    echo "  FAILED: $f" >&2
    rm -f "$OUT/$f"
  fi
done
echo "fetched $n_ok file(s) into $OUT"
