# -*- Makefile -*-

############################################################################
#
#  Program:         SuperLU_MT
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   August 15, 1997
#
#  Modified:	    September 1, 1999 version 1.0
#
############################################################################

#-----------------------------------------------------------------------------
# Makefile for building the SuperLU 2.0 solver library.
#
# Original by Stuart Norris 2002-08-25
# Changes:
#-----------------------------------------------------------------------------

BUILDING_SUPERLU_MP = true

#-----------------------------------------------------------------------------
# End of initialisation: set targets and rules.
#-----------------------------------------------------------------------------

#
#  The name of the libraries to be created/linked to
#
TMGLIB       = libtmg$(SUFFIX).a
SUPERLULIB   = $(LIBDIR)/libsuperlu$(SUFFIX).a

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
