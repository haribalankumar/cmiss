# ============================================================================
#
#  Makefile:  IRIX64 6.5.x IP27 Machine Extras
#
#  Author:  J. David Moulton
#  Date:    February 15, 2000.
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

# -------------------------
#  Linear Solvers Macros:
# -------------------------

include $(SOLVdist)/Makefile.inc

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(BOXMGdist)/make/global

# ============================================================================

OBJS = i1mach.o d1mach.o BMG_timer.o

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: BEGIN $(foreach obj, $(OBJS), $(BOXMGlibp)/$(BOXMGlibn)($(obj)))

BEGIN:
	@echo " - building/checking machine objects in $(BOXMGlibn) ... "

clean:
	@echo " - the source subdirectory: machine/$(MACHINE)"
	@$(RM) -f $(OBJS)


