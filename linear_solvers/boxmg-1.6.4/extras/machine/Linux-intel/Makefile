# ============================================================================
#
#  Makefile:  Linux-Intel Machine Extras
#
#  Author:  J. David Moulton
#  Date:    January 8, 2000.
#
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
  MAKE=make
endif

ifndef BOXMGdist
  BOXMGdist = ../../..
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

# ============================================================================

OBJS = i1mach.o d1mach.o BMG_timer.o

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: BEGIN $(foreach obj, $(OBJS), $(BOXMGlibp)/$(BOXMGlibn)($(obj)))

BEGIN:
	@echo " - building machine objects in $(BOXMGlibn) ... "

clean:
	@echo " - the source subdirectory: machine/$(MACHINE)"
	@$(RM) -f $(OBJS)


