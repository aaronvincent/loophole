#!/usr/bin/env bash
# Build a PackageCompiler sysimage of LoopAbundance for the GAMBIT worker.
# Optional: it cuts each worker's startup from tens of seconds to ~1 s.
#
# Usage (from anywhere):
#   gambit/build_sysimage.sh
#   LOOP_SYSIMG=/path/to/out.so JULIA=/path/to/julia gambit/build_sysimage.sh
#
# The image is multiversioned (cpu_target below) so it also runs on cluster
# nodes whose CPU differs from the build host's.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/.." && pwd)"
JULIA="${JULIA:-julia}"
export LOOP_SYSIMG="${LOOP_SYSIMG:-$here/loopabundance_sys.so}"
echo "target: $LOOP_SYSIMG"

cd "$repo"
"$JULIA" --project=julia/LoopAbundance -e 'import Pkg; Pkg.instantiate()'
"$JULIA" -e '
  import Pkg; Pkg.activate(; temp=true); Pkg.add("PackageCompiler")
  using PackageCompiler
  create_sysimage(["LoopAbundance"];
    project                   = "julia/LoopAbundance",
    sysimage_path             = ENV["LOOP_SYSIMG"],
    precompile_execution_file = "gambit/worker/warmup.jl",
    cpu_target                = "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1);skylake-avx512,base(1)")'

echo "Done: $(ls -lh "$LOOP_SYSIMG")"
