# ============================================================================
#
#  Makefile:  Dendy's Symmetric Black Box Multigrid (2D and 3D)
#
#  Author:  J. David Moulton
#  Date:    May 4, 1999.
#
#  ******* CMISS Changes Follow ********
#  * Author:  Travis Austin            *
#  * Date:    January 14, 2005.        *
#  *************************************
#
#
#
#
#
#
#  ****************************************
#
#  This is a makefile for Joel Dendy's (Symmetric, Standard Coarsening)
#  Black Box multigrid. To build and install the library make sure you
#  have done the following,
# 
#    - set the appropriate environment variables 
#      (see below and/or INSTALL)
#    - customized the compiler flags as necessary
#      (see INSTALL)
#
#  Then at the prompt type the commands,
#
#    make
#
#  The library ( libboxmg$(SUFFIX).a ) will be in ../lib/$(PLATFORM).
#
#  Environment Variables:
#  ----------------------
#
#  BOXMGdist   - location of BOXMG distribution
#  BOXMGlibp   - location to build the BOXMG library
#            
#  OS          - Operating system (e.g., Linux, SunOS )
#  CPU         - CPU (e.g., i686, sun4u)
#
# ============================================================================

SHELL=/bin/sh

BUILDING_BOXMG=true

ifndef BOXMGdist
  BOXMGdist = .
endif

ifndef EXTRASdist
  EXTRASdist = $(BOXMGdist)/extras
endif

ifndef SOLVdist
  SOLVdist = ../
endif

# ============================================================================

# -------------------------
#  Linear Solvers Macros:
# -------------------------

include $(SOLVdist)/Makefile.inc

# ------------------------------
#  Addition BOXMG Global Macros:
# ------------------------------

include $(BOXMGdist)/make/global

# -------------------------
#  MAKE Definition:
# -------------------------

ifeq ($(OS),AIX)
  MAKE = gmake
  MFLAGS = 
else
  MAKE = make
  MFLAGS = --no-print-directory
endif

# -------------------------
#  Library Build Path:
# -------------------------

BOXMGlibp = $(LIBDIR)

BMG_LIB = $(BOXMGlibp)/$(BOXMGlibn)

# ======================
#  Tools: 
# ======================

MV = mv 
RM = rm
CP = cp
CD = cd

MKDIR = mkdir

CHMOD = chmod
CHGRP = chgrp

# ========================================================================= 
# -------------------------
#   TARGETS:
# -------------------------

makewithflags:
	@echo ""
	@echo "Building the boxmg library .... "
	@( $(MAKE) $(MFLAGS) $(BMG_LIB) )
	@echo ""
	@echo "Building the boxmg-extras library .... "
	@( $(MAKE) $(MFLAGS) lib-extras )
	@echo ""


# -------------------------
#  BoxMG library
# -------------------------

$(BMG_LIB):
	@( $(CD) $(SRC2D); $(MAKE) $(MFLAGS) )
	@( $(CD) $(SRC3D); $(MAKE) $(MFLAGS) )
	@( $(CD) $(SRCUT); $(MAKE) $(MFLAGS) )

clean: 
	@( $(CD) $(SRC2D); $(MAKE) $(MFLAGS) clean )
	@( $(CD) $(SRC3D); $(MAKE) $(MFLAGS) clean )
	@( $(CD) $(SRCUT); $(MAKE) $(MFLAGS) clean )

# -------------------------
#  Extras Library
# -------------------------

lib-extras:
	@( $(CD) $(EXTRASdist); $(MAKE) $(MFLAGS) ) 

lib-extras-clean: 
	@( $(CD) $(EXTRASdist); $(MAKE) $(MFLAGS) clean )

# ===========================================================================

