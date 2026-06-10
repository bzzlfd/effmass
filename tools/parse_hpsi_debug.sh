#!/bin/bash
#
# parse_hpsi_debug.sh — Extract per-band H|psi> expectation values from test output.
#
# ── Compilation ──────────────────────────────────────────────────────────
# The HPSI_DEBUG blocks in hamiltonian_callable.cpp are guarded by a
# preprocessor flag.  Enable it at CMake configure time:
#
#     cmake -B build -G Ninja -DHPSI_DEBUG=ON
#     cmake --build build
#
# To revert to normal (no debug output):
#
#     cmake -B build -G Ninja -DHPSI_DEBUG=OFF
#     cmake --build build
#
# ── Usage ────────────────────────────────────────────────────────────────
# Run directly — builds table from ctest and saves to a temp file:
#
#     bash tools/parse_hpsi_debug.sh
#
# Pipe ctest output into the script (useful for custom test runs):
#
#     ctest --test-dir build -R test-hpsi-eigen --output-on-failure 2>/dev/null \
#         | bash tools/parse_hpsi_debug.sh
#
# In either case the table is printed to stdout and also saved to a
# temporary file whose path is printed to stderr.
#
# ── Output columns ───────────────────────────────────────────────────────
#   ikpt iband        <T>      <V_loc>       <V_NL>      E_file      E_calc           DE
#
#   <T>        = ⟨ψ|T|ψ⟩          kinetic energy expectation value
#   <V_loc>    = ⟨ψ|V_loc|ψ⟩      local pseudopotential expectation value
#   <V_NL>     = ⟨ψ|V_NL|ψ⟩       nonlocal pseudopotential expectation value
#   E_file     = reference eigenvalue from OUT.EIGEN
#   E_calc     = ⟨ψ|H|ψ⟩ / ⟨ψ|ψ⟩  Rayleigh quotient from the H|psi> code
#   DE         = E_file − E_calc   discrepancy (zero if implementation is correct)
# ──────────────────────────────────────────────────────────────────────────
set -euo pipefail 2>/dev/null || set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
AWKFILE="$SCRIPT_DIR/parse_hpsi_debug.awk"
TMPFILE=$(mktemp /tmp/hpsi_debug_XXXXXXXX.txt)

if [ -p /dev/stdin ]; then
    # stdin is piped — read from pipe
    awk -f "$AWKFILE" | tee "$TMPFILE"
else
    # stdin is terminal, /dev/null, or file — run test ourselves
    (ctest --test-dir "$ROOT/build" -R test-hpsi-eigen --output-on-failure 2>/dev/null || true) \
        | awk -f "$AWKFILE" \
        | tee "$TMPFILE"
fi

echo "→ saved to ${TMPFILE}" >&2
