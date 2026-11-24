#!/usr/bin/env bash
# Build the Loophole GAMBIT-light user library (loop_loglike.so).
#
# The user library is a standalone shared object (decoupled from GAMBIT's own
# build): GAMBIT-light dlopen's it via the `user_lib` path in the YAML.  It needs
# only the gambit_light_interface header.
#
# Usage:
#   ./build.sh                 # uses the defaults below
#   GAMBIT_LIGHT_DIR=/path ./build.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$here/.." && pwd)"

GAMBIT_LIGHT_DIR="${GAMBIT_LIGHT_DIR:-$HOME/software/gambit_light}"
inc="$GAMBIT_LIGHT_DIR/gambit_light_interface/include"
if [[ ! -f "$inc/gambit_light_interface.h" ]]; then
  echo "ERROR: gambit_light_interface.h not found under $inc" >&2
  echo "       Set GAMBIT_LIGHT_DIR to your gambit_light checkout." >&2
  exit 1
fi

# IMPORTANT: build with the SAME C++ compiler GAMBIT-light was built with.  The
# user lib receives std::string / std::vector / std::map by reference across the
# dlopen boundary, and libstdc++ (GNU g++) and libc++ (Apple clang) lay those out
# differently — a mismatch corrupts the arguments and segfaults on the first
# point.  So default CXX to GAMBIT's own compiler (from its CMakeCache), then a
# GNU g++, and only fall back to plain `g++` (Apple clang on macOS) as a last
# resort.  Override by exporting CXX.
if [[ -z "${CXX:-}" ]]; then
  cache="$GAMBIT_LIGHT_DIR/build/CMakeCache.txt"
  if [[ -f "$cache" ]]; then
    CXX="$(awk -F= '/^CMAKE_CXX_COMPILER:/{print $2; exit}' "$cache")"
  fi
fi
if [[ -z "${CXX:-}" ]] || ! command -v "$CXX" >/dev/null 2>&1; then
  for c in g++-15 g++-14 g++-13 g++-12 g++-11 g++; do
    if command -v "$c" >/dev/null 2>&1; then CXX="$c"; break; fi
  done
fi
echo "Using CXX=$CXX"
if [[ "$(uname)" == "Darwin" ]] && "$CXX" --version 2>/dev/null | grep -qi "Apple clang"; then
  echo "WARNING: $CXX is Apple clang (libc++).  GAMBIT-light on this machine is" >&2
  echo "         built with GNU g++ (libstdc++); the std::string/vector ABI will" >&2
  echo "         mismatch and segfault.  Install gcc (brew install gcc) or set" >&2
  echo "         CXX to the g++ from GAMBIT's CMakeCache." >&2
fi

out="$here/cpp/loop_loglike.so"

extra=()
if [[ "$(uname)" == "Darwin" ]]; then
  # Allow the interface symbols (gambit_light_invalid_point/error/warning,
  # resolved by GAMBIT at load time) to be undefined in the standalone .so.  Use
  # the -Wl, form so it works through both the clang and GNU g++ drivers.
  extra+=(-Wl,-undefined,dynamic_lookup)
fi

echo "Building $out"
"$CXX" -std=c++17 -O2 -fPIC -shared \
  "$here/cpp/loop_loglike.cpp" \
  "$here/cpp/julia_worker.cpp" \
  "$here/cpp/likelihoods.cpp" \
  -I "$inc" \
  "${extra[@]}" \
  -o "$out"

# Post-build ABI check: the .so and the gambit binary must share a libstdc++
# (GNU g++) major version, or the std::string/vector ABI mismatches and gambit
# segfaults on the first point.  GCC stamps its version into the ELF .comment
# section; read it from both and compare.  Linux-only — readelf/.comment are
# absent on the macOS dev build, where this is a no-op (soft skip).
gcc_major() {
  command -v readelf >/dev/null 2>&1 || return 0
  readelf -p .comment "$1" 2>/dev/null \
    | grep -oiE 'GCC: \(.*\) [0-9]+' | grep -oE '[0-9]+$' | head -1
}
so_gcc="$(gcc_major "$out")"
gambit_bin="$GAMBIT_LIGHT_DIR/gambit"
gambit_gcc="$([[ -x "$gambit_bin" ]] && gcc_major "$gambit_bin")"
if [[ -n "$so_gcc" && -n "$gambit_gcc" ]]; then
  if [[ "$so_gcc" == "$gambit_gcc" ]]; then
    echo "ABI check: loop_loglike.so and gambit both built with gcc $so_gcc — OK."
  else
    echo "ERROR: loop_loglike.so built with gcc $so_gcc but gambit is gcc" >&2
    echo "       $gambit_gcc — libstdc++ ABI mismatch, will segfault on point 1." >&2
    echo "       Rebuild with the matching g++ (export CXX=...) or check that" >&2
    echo "       $GAMBIT_LIGHT_DIR/build/CMakeCache.txt points at gambit's compiler." >&2
    exit 1
  fi
else
  # Either .comment was stripped (common on a release gambit) or readelf is
  # unavailable (macOS): can't compare, so warn rather than fail.
  echo "ABI check: could not read gcc version from both binaries (.comment" >&2
  echo "           stripped or no readelf) — skipping the compiler match check." >&2
fi

cat <<EOF
Built: $out

Before running gambit_light, export the worker settings (or set them in your
shell profile / job script):

  export LOOP_JULIA_BIN=$(command -v julia || echo julia)
  export LOOP_JULIA_PROJECT=$repo_root/julia/Loophole
  export LOOP_WORKER_JL=$here/worker/worker.jl

Then set user_lib in one of gambit/yaml/*.yaml to $out and run gambit on it.
EOF
