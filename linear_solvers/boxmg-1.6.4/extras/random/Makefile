# ============================================================================
#
#  Makefile:  Psuedo-Random Extras
#
#  Author:  J. David Moulton
#  Date:    Februrary 29, 2000.
#
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
  MAKE=make
endif

ifndef BOXMGdist
  BOXMGdist = ../..
endif

ifndef EXTRASdist
  EXTRASdist = $(BOXMGdist)
endif


ifndef SOLVdist
  SOLVdist = $(BOXMGdist)/..
endif

# ============================================================================

# -------------------------
#  Linear Solvers Macros:
# -------------------------

include $(SOLVdist)/Makefile.inc

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(BOXMGdist)/make/global

# ============================================================================

OBJS = rand.o

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: BEGIN $(foreach obj, $(OBJS), $(BOXMGlibp)/$(BOXMGlibn)($(obj)))

BEGIN:
	@echo " - building/checking random objects in $(BOXMGlibn) ... "

clean:
	@echo " - the source subdirectory: random"
	@$(RM) -f $(OBJS)

