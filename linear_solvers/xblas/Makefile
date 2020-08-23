# -*- Makefile -*-
# 
# For use with GNU make.

#-----------------------------------------------------------------------------
# Makefile for building the eXtended BLAS library, which contains routines
# for initialising arrays, for finding max/min values of vectors, integer
# equivalents of BLAS routines, and parallel verions of level 1 BLAS routines.
#
# Original by Stuart Norris 2002-07-29
# Changes:
#-----------------------------------------------------------------------------

BUILDING_XBLAS = true

include ../Makefile.inc

#-----------------------------------------------------------------------------
# End of initialisation: set targets and rules.
#-----------------------------------------------------------------------------

LIB = $(LIBDIR)/libxblas$(SUFFIX).a
OBJ = icopy.o \
	pdamax.o pdasum.o pdaxpy.o pdcopy.o \
	pddot.o pdload.o pdnrm2.o pdscal.o
# OBJ = camax.o camin.o cload.o damax.o damin.o dload.o dmaxmn.o iamax.o iamin.o\
#       iasum.o icopy.o iiamax.o iload.o imaxmn.o iswap.o lcopy.o lload.o lswap.o\
#       pcamax.o pcamin.o pcload.o pdamax.o pdamin.o pdasum.o pdaxpy.o pdcopy.o\
#       pddot.o pdload.o pdmaxmn.o pdnrm2.o pdscal.o pdswap.o piamax.o piamin.o\
#       piasum.o picopy.o piiamax.o piload.o pimaxmn.o piswap.o plcopy.o plload.o\
#       plswap.o psamax.o psamin.o psload.o psmaxmn.o pzamax.o pzamin.o pzload.o\
#       samax.o samin.o sload.o smaxmn.o zamax.o zamin.o zload.o

$(LIB): $(OBJ)
	$(AR) $(AFLAGS) $(LIB) $(OBJ)

# Remove objects so that they aren't reused for a different build.
.INTERMEDIATE: $(OBJ)

#-----------------------------------------------------------------------------
# Make directories and clean up
#-----------------------------------------------------------------------------
clean:
	rm -rf $(CLEAN_OBJS)

clobber:
	rm -rf $(CLEAN_OBJS) $(LIB)
