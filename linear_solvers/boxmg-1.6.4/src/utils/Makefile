# ===========================================================================
#  
#    Makefile: Symmetric Black Box Multigrid
#
# ===========================================================================
# ***************************  FILE LIST: ***********************************
# ===========================================================================
#      
#    ~~~~~~~~~~~~~~~~~~~~~~~ 
#    BBMG Wrappers:
#    ~~~~~~~~~~~~~~~~~~~~~~~
#
#    ???
#
#    -----------------------
#    Supported Wrappers:
#    -----------------------
#
#
# ===========================================================================

# ===========================================================================

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

# --------------------------------
#  Additional BOXMG Global Macros:
# --------------------------------

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
  FINCDIR += -I$(BOXMGdist)/include -I$(CMISS_SRC)/source
else
  FINCDIR += -I$(BOXMGdist)/include -I$(CMISS_SRC)/source
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

# ==========================================================================

# ------------------
#  OBJECT MODULES:
# ------------------

O_CONVERT = BMG_CONVERT_csr_matrix.o      \
            BMG_CONVERT_rhs.o             \
            BMG_CONVERT_soln.o            \
            BMG2_CONVERT_csr_matrix.o     \
            BMG3_CONVERT_csr_matrix.o     \
            BMG2_CONVERT_rhs.o            \
            BMG3_CONVERT_rhs.o            \
            BMG2_CONVERT_soln.o           \
            BMG3_CONVERT_soln.o           

O_SETUP   = BMG_SETUP.o                   \
            BMG_PCG_SETUP.o               \
            BMG_memory.o                  \
            BMG_SET_parms.o               \
            BMG_PCG_SET_parms.o           \
            BMG_SymStd_SETUP_PtrWork.o    \
            BMG_SymStd_GET_pointers.o     \
            BMG_INIT_RES.o

O_UTILS   = BMG_UTILS_rV_zero.o           \
            BMG_UTILS_iV_zero.o           \
            BMG_CONSTRUCT_CM_IJK_map.o    \
            BMG_CONSTRUCT_BMG_IJK_map.o   \
            BMG_IJK_MAP_SIZE.o

O_SOLVE   = BMG_SOLVE.o                   \
            BMG_PCG_SOLVE.o


OBJS = $(O_CONVERT) $(O_SETUP) $(O_UTILS) $(O_SOLVE)

# ==========================================================================
# -----------------
#   TARGETS:
# -----------------

$(BMG_LIB): BEGIN $(OBJS)
	$(AR) $(AFLAGS) $(BMG_LIB) $(OBJS)
	$(RANLIB) $(BMG_LIB)

BEGIN:
	@echo " - building util objects in $(FFLAGS) $(BOXMGlibn) ... "

BMG_memory.o: BMG_memory.F
	$(FC) $(FFLAGS) $(FPPFLAGS) -c $<

# Remove objects so that they aren't reused for a different build.
.INTERMEDIATE: $(OBJS) BMG_memory.o

clean:
	@echo " - the source subdirectory: utils"
	@$(RM) -rf $(OBJS) *~ rii_files





