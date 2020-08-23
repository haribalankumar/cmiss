# ============================================================================
#  
#    Makefile: Two-Dimensional Symmetric Black Box Multigrid
#
#    Original Code:    Joel Dendy and Victor Bandy
#    Port and updates: David Moulton and Joel Dendy
# 
# ===========================================================================
# ***************************  FILE LIST: ***********************************
# ===========================================================================
#
#
# ===========================================================================

# ============================================================================

SHELL=/bin/sh

BUILDING_BOXMG=true

ifndef BOXMGdist
  BOXMGdist = ../..
endif

ifndef SOLVdist
  SOLVdist = ../../..
endif

BMG_LIB = $(BOXMGlibp)/$(BOXMGlibn)

# ============================================================================

# ----------------------------
#  CMISS Linear Solver Macros:
# ----------------------------

include $(SOLVdist)/Makefile.inc

# -------------------------
#  Standard Global Macros:
# -------------------------

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

# ---------------------
#  Includes:
# ---------------------

ifndef FINCDIR
  FINCDIR += -I$(BOXMGdist)/include
else
  FINCDIR += -I$(BOXMGdist)/include
endif

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

# ------------------
#  OBJECT MODULES:
# ------------------

O_BMG  =  BMG2_SymStd_SOLVE_boxmg.o        \
          BMG2_SymStd_SOLVE_cg.o           \
          BMG2_SymStd_ncycle.o             \
          BMG2_SymStd_updown.o             \
          BMG2_SymStd_GET_pointers.o       \
          BMG2_SymStd_SETUP_parms.o        \
          BMG2_SymStd_SETUP_PtrWork.o      \
          BMG2_SymStd_SETUP_PtrGrid.o      \
          BMG2_SymStd_SETUP_space.o        \
          BMG2_SymStd_SETUP_recip.o        \
          BMG2_SymStd_SETUP_parts.o        \
          BMG2_SymStd_SETUP_interp_OI.o    \
          BMG2_SymStd_SETUP_interp_BI.o    \
          BMG2_SymStd_SETUP_ITLI_ex.o      \
          BMG2_SymStd_SETUP_ITLI_bl.o      \
          BMG2_SymStd_SETUP_cg_LU.o        \
          BMG2_SymStd_SETUP_relax.o        \
          BMG2_SymStd_SETUP_lines_x.o      \
          BMG2_SymStd_SETUP_lines_y.o      \
          BMG2_SymStd_restrict.o           \
          BMG2_SymStd_interp_add.o         \
          BMG2_SymStd_relax.o              \
          BMG2_SymStd_relax_GS.o           \
          BMG2_SymStd_relax_lines_x.o      \
          BMG2_SymStd_relax_lines_y.o      \
          BMG2_SymStd_residual.o           \
          BMG2_SymStd_ErrTrap.o

# No time to support periodic stuff (2000/07/24 - JDM)
#
O_BMGp  = BMG2_PerSymStd_SETUP_cofp.o      \
          BMG2_PerSymStd_SETUP_cg_LU.o     \
          BMG2_PerSymStd_interp_add.o      \
          BMG2_PerSymStd_residual.o        \
          BMG2_PerSymStd_relax.o           \
          BMG2_PerSymStd_SOLVE_cg.o        

O_PCG   = BMG2_SymStd_SOLVE_pcg.o          \
          BMG2_SymStd_PRECON_diag.o        \
          BMG2_SymStd_PRECON_boxmg.o


O_UTILS = BMG2_SymStd_UTILS_rV_zero.o      \
          BMG2_SymStd_UTILS_norm_l2.o      \
          BMG2_SymStd_UTILS_dot_l2.o       \
          BMG2_SymStd_UTILS_matvec.o       \
          BMG2_SymStd_UTILS_daxpy.o        \
          BMG2_SymStd_UTILS_dscal.o        \
          BMG2_SymStd_UTILS_dcopy.o        \
          BMG2_SymStd_UTILS_dxpby.o        \
          BMG2_SymStd_UTILS_reswt.o        \
          BMG2_SymStd_UTILS_zero_mean.o    \
          BMG2_SymStd_UTILS_compute_anorm.o

O_DUMP  = BMG2_SymStd_DUMP_parms.o         \
          BMG2_SymStd_DUMP_stencil.o       \
          BMG2_SymStd_DUMP_restrict.o      \
          BMG2_SymStd_DUMP_vector.o


OBJS = $(O_BMG) $(O_BMGp) $(O_PCG) $(O_UTILS) $(O_DUMP)

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

$(BMG_LIB): BEGIN $(OBJS)
	$(AR) $(AFLAGS) $(BMG_LIB) $(OBJS)
	$(RANLIB) $(BMG_LIB)

# Remove objects so that they aren't reused for a different build.
.INTERMEDIATE: $(OBJS)

BEGIN:
	@echo " - building boxmg-sym-std-2D in $(BOXMGlibp)/$(BOXMGlibn) ... "

clean:
	@echo " - the source subdirectory: boxmg-sym-std-2D"
	@$(RM) -rf $(OBJS) *~ rii_files







