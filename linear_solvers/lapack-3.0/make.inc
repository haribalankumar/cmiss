# -*- Makefile -*-

####################################################################
#  LAPACK make include file.                                       #
#  LAPACK, Version 3.0                                             #
#  June 30, 1999                                                  #
####################################################################

#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
#BLASLIB      = ../../blas$(PLAT).a
#BLASLIB     = -xlic_lib=sunperf_mt
#BLASLIB     = -xlic_lib=sunperf
#LAPACKLIB    = lapack$(PLAT).a
#TMGLIB       = tmglib$(PLAT).a
#EIGSRCLIB    = eigsrc$(PLAT).a
#LINSRCLIB    = linsrc$(PLAT).a

#-----------------------------------------------------------------------------
# Makefile for building the SuperLU 2.0 solver library.
#
# Original by Stuart Norris 2002-08-25
# Changes:
#-----------------------------------------------------------------------------

BUILDING_LAPACK = true

#-----------------------------------------------------------------------------
# End of initialisation: set targets and rules.
#-----------------------------------------------------------------------------

#
#  The name of the libraries to be created/linked to
#
TMGLIB       = libtmg$(SUFFIX).a
LAPACKLIB    = $(LIBDIR)/liblapack$(SUFFIX).a

#
#  The archiver and the flag(s) to use when building archive (library)
#  Convert between the rest of the worlds (and make's) default names
#  and those used in the SuperLU library.
#
ARCH         = $(AR)
ARCHFLAGS    = $(AFLAGS)
FORTRAN	     = $(FC)
LOADER       = $(CC)
LOADOPTS     =
