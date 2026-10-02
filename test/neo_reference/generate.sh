#!/bin/bash
# Regenerate the Fortran NEO reference data used by test/runtests_neo_native.jl.
#
# Not run in CI. Needs a built gacode tree with its environment sourced
# (GACODE_ROOT, GACODE_PLATFORM) and the neo_dump driver from
# utilities/serial_neo (built here if missing).
#
# Usage:  bash test/neo_reference/generate.sh [case ...]
#         (default: every directory here that holds an input.neo, plus the
#          gacode regression cases listed in REG_CASES, copied in if absent)
#
# Per case directory the following are kept:
#   input.neo              the input (hand-written or copied from gacode)
#   out.neo.transport      per-species dke results (e16.8)
#   out.neo.transport_gv   gyroviscous fluxes (e16.8)
#   out.neo.transport_flux GB-normalised fluxes incl. the tgyro block
#   out.neo.prec           check_sum (what gacode's regression test compares)
#   out.neo.dump           full-precision setup arrays (neo_dump.f90; small*, reg12 only)
#   out.neo.f              full-precision solution vector g (small*, reg12 only)
#   out.neo.dump_fcoll     fcoll tables (only kept for the `small` case)
#   out.neo.diagnostic_coll{test,field}  mono-basis matrices if WRITE_CMOMENTS_FLAG=1
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PKG=$(cd "$HERE/../.." && pwd)
: "${GACODE_ROOT:?set GACODE_ROOT (and GACODE_PLATFORM) and source the gacode env first}"
: "${GACODE_PLATFORM:?set GACODE_PLATFORM}"

REG_CASES="reg01 reg02 reg03 reg04 reg05 reg06 reg07 reg08 reg09 reg10 reg11 reg12 reg13 reg14 reg15"

# NEO's UMFPACK/BLAS are OpenMP-threaded; a few threads are plenty for these sizes
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

DUMP="$PKG/utilities/serial_neo/neo_dump"
if [ ! -x "$DUMP" ]; then
    make -C "$PKG/utilities/serial_neo" neo_dump
fi

# copy regression inputs in if absent
for c in $REG_CASES; do
    if [ ! -f "$HERE/$c/input.neo" ]; then
        mkdir -p "$HERE/$c"
        cp "$GACODE_ROOT/neo/tools/input/$c/input.neo" "$HERE/$c/input.neo"
    fi
done

if [ $# -gt 0 ]; then
    CASES="$*"
else
    CASES=$(cd "$HERE" && for d in */; do [ -f "$d/input.neo" ] && echo "${d%/}"; done)
fi

for c in $CASES; do
    dir="$HERE/$c"
    echo "== $c"
    (
        cd "$dir"
        rm -f out.neo.* input.neo.gen
        python "$GACODE_ROOT/neo/bin/neo_parse.py" > parse.log 2>&1
        : > out.neo.run
        "$DUMP" > dump.log 2>&1
        rm -f input.neo.gen parse.log dump.log
        # trim what the tests do not read
        rm -f out.neo.vel out.neo.vel_fourier out.neo.phi out.neo.grid \
              out.neo.theory out.neo.theory_nclass out.neo.species \
              out.neo.diagnostic_geo out.neo.diagnostic_geo2 out.neo.diagnostic_coll \
              out.neo.rotation out.neo.diagnostic_rot out.neo.transport_exp \
              out.neo.gxi out.neo.gxi_t out.neo.gxi_x out.neo.equil out.neo.run
        if [ "$c" != "small" ]; then rm -f out.neo.dump_fcoll; fi
        # the full-precision dumps are large (~24 bytes per value); keep them for the
        # CI cases and one full-size case, the others are checked end to end
        case "$c" in small*|reg12) ;; *) rm -f out.neo.dump out.neo.f ;; esac
        cat out.neo.prec; echo
    )
done
