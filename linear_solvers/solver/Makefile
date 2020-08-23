# -*- Makefile -*-
# 
# For use with GNU make.

#-----------------------------------------------------------------------------
# Makefile for building the CMISS solver library, which contains iterative
# linear solvers, and interfaces to the LAPACK, Harwell, Umfpack and SuperLU
# direct solvers.
#
# Original by Stuart Norris 2002-07-29
# Changes:
#-----------------------------------------------------------------------------

BUILDING_SOLVER = true

include ../Makefile.inc

#-----------------------------------------------------------------------------
# End of initialisation: set targets and rules.
#-----------------------------------------------------------------------------

LIB = $(LIBDIR)/libsolver$(SUFFIX).a
EXE = soltest-$(EXESUF)
OBJ = cg.o dense.o jacobi.o entry.o ic0.o ic1.o ilu0.o ilu1.o ildl0.o ildl1.o\
      misc.o matrix.o test.o umfpack_c.o umfpack4.o superLU.o\
      superLU_c$(MP_SUFF).o

LIBLIST = -L$(LIBDIR) -lsolver$(SUFFIX) -lxblas$(SUFFIX) -lumfpack-4.0$(SERIAL)\
          -lsuperlu$(SUFFIX) $(SCILIB)
TESTOBJ = main.o cmlib.o malloc.o fd.o timer.o

$(LIB): $(OBJ)
	$(AR) $(AFLAGS) $(LIB) $(OBJ)

$(EXE): $(TESTOBJ) $(LIB)
	$(LD) -o $@ $(LDFLAGS) $(TESTOBJ) $(LIBLIST)

all: $(EXE)

exe: $(EXE)

# Remove objects so that they aren't reused for a different build.
.INTERMEDIATE: $(OBJ)

# Dependencies from INCLUDE statements
$(OBJ): cbdi02.cmn cbfe01.cmn mach00.inc solver.inc

vpath %.cmn $(CMISS_SRC)
%.cmn:
	$(error $@ not found in $(CMISS_SRC))
vpath %.inc $(CMISS_SRC)
%.inc:
	$(error $@ not found in $(CMISS_SRC))

#-----------------------------------------------------------------------------
# Make directories and clean up
#-----------------------------------------------------------------------------
clean:
	rm -rf $(CLEAN_OBJS)

clobber:
	rm -rf $(CLEAN_OBJS) $(LIB)
