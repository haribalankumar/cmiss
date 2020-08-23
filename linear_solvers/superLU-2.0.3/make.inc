# -*- Makefile -*-

############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################

#-----------------------------------------------------------------------------
# Makefile for building the SuperLU 2.0 solver library.
#
# Original by Stuart Norris 2002-08-25
# Changes:
#-----------------------------------------------------------------------------

BUILDING_SUPERLU = true

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
