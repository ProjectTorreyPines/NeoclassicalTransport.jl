# Build the serial (no-MPI) NEO binary out-of-tree against the stock gacode
# neo_lib.a, reusing the compiler/flags of YOUR gacode platform (the same
# make.inc the regular neo build uses), so it works on any system where
# gacode's neo builds.
#
# Prerequisites: a built gacode tree (neo compiled) and its environment
# sourced, so GACODE_ROOT and GACODE_PLATFORM are set.
#
# Usage:
#   make          # -> ./neo_serial
#   make clean
#
# Then point run_neo at it, e.g. in Julia:
#   ENV["NEO_EXECUTABLE"] = joinpath(<this dir>, "neo_serial")

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

FFLAGS = ${FOPT} ${FOMP}

# same support libraries the stock neo executable links (gacode neo/Makefile)
EXTRA_LIBS = \
	${GACODE_ROOT}/f2py/expro/expro_lib.a \
	${GACODE_ROOT}/f2py/geo/geo_lib.a \
	${GACODE_ROOT}/shared/math/math_lib.a \
	${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a \
	${GACODE_ROOT}/shared/nclass/nclass_lib.a

NEO_LIB = ${GACODE_ROOT}/neo/src/neo_lib.a

all: neo_serial

neo_serial: neo_serial.f90 ${NEO_LIB}
	$(FC) $(FMATH) $(FFLAGS) -c neo_serial.f90
	$(FC) $(FFLAGS) -o neo_serial neo_serial.o ${NEO_LIB} $(EXTRA_LIBS) $(LMATH) $(NN_LIB)
	rm -f neo_serial.o

${NEO_LIB}:
	$(error $(NEO_LIB) not found - build gacode neo first)

clean:
	rm -f neo_serial neo_serial.o
