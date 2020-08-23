# ============================================================================
#  
#    Makefile: Three-Dimensional Black Box Multigrid
#
#    Original Code:    Joel Dendy
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

# -------------------------
#  Linear_Solver Macros:
# -------------------------

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
  FINCDIR = -I$(BOXMGdist)/include
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

O_BMG   = BMG3_SymStd_SOLVE_boxmg.o          \
          BMG3_SymStd_SOLVE_cg.o             \
          BMG3_SymStd_ncycle.o               \
          BMG3_SymStd_updown.o               \
          BMG3_SymStd_GET_pointers.o         \
          BMG3_SymStd_SETUP_parms.o          \
          BMG3_SymStd_SETUP_PtrWork.o        \
          BMG3_SymStd_SETUP_PtrGrid.o        \
          BMG3_SymStd_SETUP_space.o          \
          BMG3_SymStd_SETUP_parts.o          \
          BMG3_SymStd_SETUP_interp_OI.o      \
          BMG3_SymStd_SETUP_interp_TI.o      \
          BMG3_SymStd_SETUP_ITLI07_bl.o      \
          BMG3_SymStd_SETUP_ITLI27_bl.o      \
          BMG3_SymStd_SETUP_ITLI_Izyx.o      \
          BMG3_SymStd_SETUP_relax.o          \
          BMG3_SymStd_SETUP_recip.o          \
          BMG3_SymStd_SETUP_cg_recip.o       \
          BMG3_SymStd_SETUP_planes_xy.o      \
          BMG3_SymStd_SETUP_planes_yz.o      \
          BMG3_SymStd_SETUP_planes_xz.o      \
          BMG3_SymStd_SETUP_cg_LU.o          \
          BMG3_SymStd_restrict.o             \
          BMG3_SymStd_interp_add.o           \
          BMG3_SymStd_relax.o                \
          BMG3_SymStd_relax_CG.o             \
          BMG3_SymStd_relax_GS.o             \
          BMG3_SymStd_relax_ptwise_GS.o      \
          BMG3_SymStd_relax_planes_xy.o      \
          BMG3_SymStd_relax_planes_yz.o      \
	  BMG3_SymStd_relax_planes_xz.o      \
          BMG3_SymStd_residual.o             \
          BMG3_SymStd_injection.o            \
          BMG3_SymStd_ErrTrap.o              \

O_2D_3D = BMG3_SymStd_COPY_SO_xy.o           \
          BMG3_SymStd_COPY_SO_yz.o           \
          BMG3_SymStd_COPY_SO_xz.o           \
          BMG3_SymStd_COPY_RHS_xy.o          \
          BMG3_SymStd_COPY_RHS_yz.o          \
          BMG3_SymStd_COPY_RHS_xz.o          \
          BMG3_SymStd_COPY_rV_23_xy.o        \
          BMG3_SymStd_COPY_rV_23_yz.o        \
          BMG3_SymStd_COPY_rV_23_xz.o        \
          BMG3_SymStd_COPY_rV_32_xy.o        \
          BMG3_SymStd_COPY_rV_32_yz.o        \
          BMG3_SymStd_COPY_rV_32_xz.o

O_PCG   = BMG3_SymStd_SOLVE_pcg.o            \
          BMG3_SymStd_PRECON_boxmg.o         \
          BMG3_SymStd_PRECON_diag.o

O_UTILS = BMG3_SymStd_UTILS_rV_zero.o        \
          BMG3_SymStd_UTILS_rV_rand.o        \
          BMG3_SymStd_UTILS_norm_l2.o        \
          BMG3_SymStd_UTILS_daxpy.o          \
          BMG3_SymStd_UTILS_dcopy.o          \
          BMG3_SymStd_UTILS_dot_l2.o         \
          BMG3_SymStd_UTILS_dscal.o          \
          BMG3_SymStd_UTILS_dxpby.o          \
          BMG3_SymStd_UTILS_matvec.o         \
          BMG3_SymStd_UTILS_reswt.o          \
          BMG3_SymStd_UTILS_diag_scal.o      \
          BMG3_SymStd_UTILS_zero_mean.o      \
          BMG3_SymStd_UTILS_compute_anorm.o

O_DUMP  = BMG3_SymStd_DUMP_parms.o           \
          BMG3_SymStd_DUMP_stencil.o         \
          BMG3_SymStd_DUMP_restrict.o        \
          BMG3_SymStd_DUMP_vector.o          \
          BMG3_SymStd_PRINT_matrix.o        

OBJS = $(O_BMG) $(O_PCG) $(O_2D_3D) $(O_UTILS) $(O_DUMP)

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
	@echo " - building boxmg-sym-std-3D objects in $(BOXMGlibn) ... "

clean:
	@echo " - the source subdirectory: boxmg-sym-std-3D"
	@$(RM) -rf $(OBJS) *~ rii_files





