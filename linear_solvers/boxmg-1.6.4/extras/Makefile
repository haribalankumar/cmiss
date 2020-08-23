# ============================================================================
#
#  Makefile:  Extras for Dendy's Black Box Multigrid 
#
#  Author:  J. David Moulton
#  Date:    Someday, 1993.	
#
#  ****************************************
#
#  This is a makefile for the extras that are needed for Joel Dendy's 
#  (Symmetric with Standard Coarsening) Black Box multigrid. Specifically,
#  pieces of the well known packages LAPACK (with self-contained BLAS),
#  are included.  These codes are all freely available from netlib
#
#             http://www.netlib.org/
#
#  The versions here have minor changes to improve portability (in our
#  opinion).  These changes should be benign, but no warrenty is made
#  implied (see usual legal stuff in LICENSE).  These routines are 
#  provided only for your convenience, optimized versions may be available
#  on your system.  
#
#  ******************************************
#
#  To build the extras library make sure you have done the following,
# 
#    - set the appropriate environment variables 
#      (see top-level INSTALL file of boxmg)
#    - customized the compiler flags as necessary
#      (see top-level INSTALL file of boxmg)
#
#  Then at the prompt type the command,
#
#    make
#
#  At the moment there are no tests for this portion of the library
#  If you have concerns there may be appropriate tests available 
#  from netlib.  
#
#  If you plan on using a vendor optimized version of LAPACK, you
#  may wish to exclude this from the boxmg-extras library (i.e., to 
#  avoid mixups in the linking phase).  This is easily accomplished
#  by passing arguments to make,
#
#    make LAPACK=no
#
#  Similarly, if you don't want this random number generator, use
#  RANDOM=no.  By default all components are included in the library.
#
#  Environment Variables:
#  ----------------------
#
#  BOXMGdist    - location of boxmg distribution
#                 ( default: ../boxmg)
#
#  EXTRASdist   - location of boxmg-extras distribution
#                 ( default: $(BOXMGdist)/extras )
#  EXTRASlibp   - location to build the boxmg-extras library
#                 ( default: $(BOXMGdist)/lib )
#  EXTRASlibi   - location to install the boxmg-extras library
#                 ( default: none )
#            
#  OS           - Operating system (e.g., Linux, SunOS )
#  CPU          - CPU (e.g., intel, sparc, alpha, mips)
#
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
  MAKE=make
endif

ifndef BOXMGdist
  BOXMGdist = ..
endif

ifndef EXTRASdist
  EXTRASdist = $(BOXMGdist)/extras
endif

ifndef SOLVdist
  SOLVdist = $(BOXMGdist)/..
endif

# ============================================================================

# -------------------------
#  Linear Solvers Macros:
# -------------------------x

include $(SOLVdist)/Makefile.inc

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(BOXMGdist)/make/global

# -------------------------
#  Path Definitions:
# -------------------------

RAND = $(EXTRASdist)/random
MACH = $(EXTRASdist)/machine/$(MACHINE)

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
#  Select Components
# -------------------------

PIECES = BUILD_MACHINE

ifndef RANDOM
  RANDOM=yes
endif

ifeq ($(RANDOM),yes)
  PIECES += BUILD_RANDOM
endif

# ========================================================================= 
# -------------------------
#   TARGETS:
# -------------------------

library: $(PIECES)

BUILD_RANDOM:
	@( $(CD) $(RAND); $(MAKE) )

BUILD_MACHINE:
	@( $(CD) $(MACH); $(MAKE) )

install: 
	@echo "Installing $(EXTRASlibp)/$(EXTRASlibn) with ... "
	@echo "$(MV) $(EXTRASlibp)/$(EXTRASlibn) $(EXTRASlibi)"
	@$(MV) $(EXTRASlibp)/$(EXTRASlibn) $(EXTRASlibi)

uninstall: 
	@echo "Uninstalling $(EXTRASlibn) from $(EXTRASlibi) with ... "
	@echo "$(MV) $(EXTRASlibi)/$(EXTRASlibn) $(EXTRASlibp)"
	@$(MV) $(EXTRASlibi)/$(EXTRASlibn) $(EXTRASlibp)

clean:
	@echo "Removing objects ... "
	@( $(CD) $(RAND); $(MAKE) $(MFLAGS) clean )
	@( $(CD) $(MACH); $(MAKE) $(MFLAGS) clean )

distclean: clean
	@echo "Removing library: " $(EXTRASlibn)
	@/bin/rm -f $(EXTRASlibp)/$(EXTRASlibn)

# ===========================================================================

