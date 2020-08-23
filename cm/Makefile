# -*- makefile -*-
#
# For use with GNU make.

#----------------------------------------------------------------------------
# Makefile for compiling CMISS with individual feuser_cell.f model.
#
# Original by Karl Tomlinson 2000-10-01
# Changes:
#	
# 2000-10-03 Chris Penland	Compile parallel binary and switched compile
#				flags to optimised (from debug). Added
#				targets for clean, allclean, help.
# 2000-10-05 Karl Tomlinson     Added option OPT for debug (false) or
#				optimized.
# 2000-10-13 Karl		Added tidying of old generated files,
#				preprocessing of .cmn files,
#				dynamic linking and STATIC,
#				splitting fortran files into routines.
# 2001-03-16 Karl               Modifications for compilation of a global
#				version and for linux compilation.
# 2002-05-09 Carey              Added AIX targets
#----------------------------------------------------------------------------

# no builtin implicit rules
MAKEFLAGS = --no-builtin-rules --warn-undefined-variables

#-----------------------------------------------------------------------------

ifndef SYSNAME
  SYSNAME := $(shell uname)
  ifeq ($(SYSNAME),)
    $(error error with shell command uname)
  endif
endif

ifndef NODENAME
  NODENAME := $(shell uname -n)
  ifeq ($(NODENAME),)
    $(error error with shell command uname -n)
  endif
endif

ifndef MACHNAME
  MACHNAME := $(shell uname -m)
  ifeq ($(MACHNAME),)
    $(error error with shell command uname -m)
  endif
endif

ifndef PROF
  PROF := false
endif

ifndef DEBUG
  ifneq ($(PROF),false)
    OPT := true
  endif
  ifndef OPT
    OPT := false
  endif
  ifeq ($(OPT),false)
    DEBUG := true
  else
    DEBUG := false
  endif
endif

ifndef MP
  MP = false
endif

# set architecture dependent directories and default options

# defaults
INSTRUCTION=$(MACHNAME)
BIN_ARCH_DIR = $(INSTRUCTION)-$(OPERATING_SYSTEM)
LIB_ARCH_DIR = $(INSTRUCTION)-$(ABI)-$(OPERATING_SYSTEM)
# OLD_ARCH_DIR = $(SYSNAME)-$(ABI)

ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
  # Specify what application binary interface (ABI) to use i.e. 32, n32 or 64
  ifndef ABI
    ifdef SGI_ABI
      ABI := $(patsubst -%,%,$(SGI_ABI))
    else
      ABI = n32
    endif
  endif
  # Specify which instruction set to use i.e. -mips#
  ifndef MIPS
    # Using mips3 for most basic version on esu* machines
    # as there are still some Indys around.
    # Although mp versions are unlikely to need mips3 they are made this way
    # because it makes finding library locations easier.
    MIPS = 4
    ifeq ($(filter-out esu%,$(NODENAME)),)
      ifeq ($(ABI),n32)
        ifneq ($(DEBUG),false)
          MIPS=3
        endif
      endif
    endif
  endif
#  OLD_ARCH_DIR := mips$(MIPS)-$(ABI)
  INSTRUCTION := mips
  OPERATING_SYSTEM := irix
endif
ifeq ($(SYSNAME),Linux)
  OPERATING_SYSTEM := linux
  ifeq ($(filter-out i%86,$(MACHNAME)),)
#    OLD_ARCH_DIR := linux86
    INSTRUCTION := i686
  endif
  ifneq ($(filter $(INSTRUCTION),i686 ia64 x86_64),)# i686, ia64, x86_64 
    LIB_ARCH_DIR = $(INSTRUCTION)-$(OPERATING_SYSTEM)# no ABI required
  endif
  ifndef ABI
    ABI=32
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      ABI=32
    endif
    ifneq (,$(filter $(MACHNAME),ia64 x86_64))# ia64 or x86_64
      ABI=64
    endif
  endif
endif
ifeq ($(SYSNAME),win32)
  LIB_ARCH_DIR = $(INSTRUCTION)-$(OPERATING_SYSTEM)# no ABI
#  OLD_ARCH_DIR := win32
  ifndef DYNAMIC
    DYNAMIC = false
  endif
  ABI=32
  INSTRUCTION := i386
  OPERATING_SYSTEM := win32
endif
ifeq ($(SYSNAME),SunOS)
  OPERATING_SYSTEM := solaris
endif
ifeq ($(SYSNAME),AIX)
  ifndef ABI
    ifdef OBJECT_MODE
      ifneq ($(OBJECT_MODE),32_64)
        ABI = $(OBJECT_MODE)
      endif
    endif
  endif
  ifndef ABI
    ABI = 32
  endif
#  OLD_ARCH_DIR = aix-$(ABI)
  INSTRUCTION := rs6000
  OPERATING_SYSTEM := aix
endif

ifndef DYNAMIC
  # STATIC was the old user option but was a bit misleading because the
  # executable was still not completely static.  Check for users who are still
  # using this.
  ifndef STATIC
    ifeq ($(DEBUG),false)
      STATIC := true
    else
      STATIC := false
    endif
  endif
  ifeq ($(STATIC),false)
    DYNAMIC := true
  else
    DYNAMIC := false
  endif
endif

EXE_ABI_SUFFIX =
ifeq ($(ABI),64)
  ifeq ($(filter $(MACHNAME),ia64 x86_64),)# not ia64, x86_64 
    EXE_ABI_SUFFIX = 64
  endif
endif

ifeq ($(MP),false)
  MT_SUFFIX =
  MP_SUFFIX =
else
  MT_SUFFIX = -mt
  MP_SUFFIX = _mp
#   ifeq ($(SYSNAME),SunOS)
#     MP_SUFFIX = _mt
#   endif
endif

ifeq ($(PROF),false)
  PROF_SUFFIX =
else
  PROF_SUFFIX = -prof
endif

ifneq ($(DEBUG),false)
  DEBUG_SUFFIX = -debug
  OLD_OPT_SUFFIX = -debug
else
  DEBUG_SUFFIX = 
  OLD_OPT_SUFFIX = -opt
endif

# Exporting SYSNAME saves a call to uname here and in preprocess_fortran.
export SYSNAME NODENAME MACHNAME

SOURCE_DIR = source
# OLD_WORKING_DIR := generated/$(OLD_ARCH_DIR)$(MP_SUFFIX)$(OLD_OPT_SUFFIX)
WORKING_DIR := generated/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)
SOLVER_DIR = ${CMISS_ROOT}/linear_solvers/solver
BOXMG_DIR = ${CMISS_ROOT}/linear_solvers/boxmg-1.6.4/include
ifneq ($(DYNAMIC),false)# dynamic linking
#  OLD_EXE_DIR := dso/$(OLD_ARCH_DIR)$(MP_SUFFIX)$(OLD_OPT_SUFFIX)
  DSO_DIR := dso/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)
  EXE_DIR := $(DSO_DIR)
  EXE_NAME := cm
#   EXE_RECORD := $(EXE_DIR)/.date-$(EXE_NAME)
else
#  OLD_EXE_DIR := bin/$(OLD_ARCH_DIR)
  EXE_DIR := bin/$(BIN_ARCH_DIR)
  EXE_NAME := cm$(EXE_ABI_SUFFIX)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)
#   EXE_RECORD :=
  DSO_DIR :=
endif
EXECUTABLE := $(EXE_DIR)/$(EXE_NAME)

ifndef CMISS_ROOT
  CMISS_ROOT := ..
  GLOBAL_ROOT := $(CURDIR)
#   $(error CMISS_ROOT not defined)
else
  GLOBAL_ROOT := ${CMISS_ROOT}/cm
endif

ifndef USE_CELLML
  USE_CELLML = false
  # I think CellML is now all working under IRIX, so we might
  # as well have it on by default for the IRIX builds.
  ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
    USE_CELLML = true
  endif
  # Now also works with Linux - although need to be careful that the
  # CellML libraries are compiled with compatible versions of compilers.
  ifeq ($(SYSNAME),Linux)
    USE_CELLML = true
  endif
  # Now also works with AIX
  ifeq ($(SYSNAME),AIX)
    USE_CELLML = true
  endif
endif

ifndef CELLML_ROOT
  CELLML_ROOT := $(CMISS_ROOT)/cellml
endif

ifndef CELLML_DEBUG
  CELLML_DEBUG = false
endif

ARCHIVE_DIRS := \
  lib/$(LIB_ARCH_DIR) \
  $(GLOBAL_ROOT)/lib/$(LIB_ARCH_DIR)
C_INCLUDE_DIRS := $(SOURCE_DIR) \
  $(GLOBAL_ROOT)/$(SOURCE_DIR) \
  $(CMISS_ROOT)/perl_interpreter/source
CMN_INCLUDE_DIRS := $(WORKING_DIR) $(GLOBAL_ROOT)/$(WORKING_DIR)
INC_INCLUDE_DIRS := $(SOURCE_DIR) $(SOLVER_DIR) $(BOXMG_DIR) $(GLOBAL_ROOT)/$(SOURCE_DIR) \
  ${CMISS_ROOT}/gx/source
# cmn must be before inc because files exist in both
F_INCLUDE_DIRS := $(CMN_INCLUDE_DIRS) $(INC_INCLUDE_DIRS)

#-----------------------------------------------------------------------------
# compiling commands

CC = cc
FC = f77
EXE_LINK = $(FC)
DSO_LINK = ld
# LD_RELOCATABLE = ld -r $(CLD_FLGS) $(L_FLGS)
DBGCF_FLGS = -g#OPT=false flags for C and fortran
# Option lists
# (suboption lists become more specific so that later ones overrule previous)
CFLAGS = $(strip $(CFL_FLGS) $(CFE_FLGS) $(CF_FLGS) $(C_FLGS))
FFLAGS = $(strip $(CFL_FLGS) $(CFE_FLGS) $(CF_FLGS) $(F_FLGS))
CPPFLAGS := $(addprefix -I, $(C_INCLUDE_DIRS) )
FPPFLAGS := $(addprefix -I, $(F_INCLUDE_DIRS) )
ELFLAGS = $(strip $(CFL_FLGS) $(L_FLGS) $(CFE_FLGS))
DLFLAGS = $(strip $(CFL_FLGS) $(L_FLGS) $(D_FLGS))
ifneq ($(DEBUG),false)
  CFLAGS += $(strip $(DBGCF_FLGS) $(DBGC_FLGS))
  FFLAGS += $(strip $(DBGCF_FLGS) $(DBGF_FLGS))
  CPPFLAGS += -DDEBUG
else
  CFLAGS += $(strip $(OPTCFE_FLGS) $(OPTCF_FLGS) $(OPTC_FLGS))
  FFLAGS += $(strip $(OPTCFE_FLGS) $(OPTCF_FLGS) $(OPTF_FLGS))
  ELFLAGS += $(OPTCFE_FLGS)
endif
ifneq ($(MP),false)
  CFLAGS += $(MP_FLGS)
  FFLAGS += $(MP_FLGS)
  ifneq ($(filter-out IRIX%,$(SYSNAME)),)
    #On IRIX we are going to use ld so we don't want the normal flags
    ELFLAGS += $(MP_FLGS)
  endif
endif
# suboption lists
CFL_FLGS =#	flags for C fortran and linking
L_FLGS =#	flags for linking only
CFE_FLGS =#	flags for C fortran and linking executables only
CF_FLGS = -c#	flags for C and fortran only
C_FLGS =#       flags for C only
F_FLGS =#       flags for fortran only
D_FLGS = -shared#     for linking dynamic shared objects only
# DBGCF_FLGS = -g#OPT=false flags for C and fortran
DBGC_FLGS =#	OPT=false flags for C only
DBGF_FLGS =#	OPT=false flags for fortran only
OPTCFE_FLGS =#	OPT=true flags for C and fortran and linking executables
OPTCF_FLGS = -O#OPT=true flags for C and fortran only
OPTC_FLGS =#	OPT=true flags for C only
OPTF_FLGS =#	OPT=true flags for fortran only

# The list of objects may be too long for the operating system limit on
# argument length so the list of objects is stored in a file.  This linker
# arguments for this file depend on the linker.  If the linker cannot
# interpret such a file then try to use the shell and hope the list isn't too
# long.
olist_args = `cat $1`

DSO_LIBS =# libraries required to build DSOs
# Don't overwrite dsos as they may be in use.
# Remove the executable permission as the dsos are not executable.
# (But HPUX requires dsos to be executable.)
# -w is merely a safety to protect the dso from being overwritten
define dso_build
	$(DSO_LINK) -o $1 $(DLFLAGS) $2 $(filter %.a, $^) $(addprefix -l, $(DSO_LIBS))
	chmod -wx $1
endef

ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
  #SAB.  On IRIX -hiddens-file and -mp conflict when they get to ld due to some
  #created flag -EB.  We cannot find out what this does and if we link with ld
  #explictly we have no problems, so I have changed all IRIX to use ld to link.
  EXE_LINK = ld

  # Remarks ignored because of declarations in included files:
  # remark 2272 `routine' `declared and never referenced'
  # remark 2274 `named constant' `declared and never referenced'

  # Warnings for things that we choose to have:
  # remark cc-3201: parameter never referenced.

  # WARNING 131: Multiply defined weak symbol - get this with the CellML
  # libraries because of the way IRIX does C++ templates...seems acceptable
  # to simply turn this warning off.

  # The following warning means that the execution of the program is seriously
  # different from that intended:
  # cc-1999 cc: WARNING File = zle_tricky.c, Line = 2145
  # "jumping out of a block containing VLAs" is not currently implemented
  CFLAGS += -DEBUG:error=1999 -woff 3201
  F_FLGS += -nocpp -G 0 -noextend_source \
	-woff 2272 -woff 2274 -woff 131
  CF_FLGS += -use_readonly_const -fullwarn
  ifeq ($(MP),false)
    DBGCF_FLGS += -DEBUG:trap_uninitialized:subscript_check:verbose_runtime
  else
    # subscript_check causes compiler errors with -mp:
    # ### Compiler Error in file ASSEMBLE1_DYNAM.f during MP Lowering phase:
    # ### out of context pragma (ATOMIC) in MP {standalone pragma} processing
    # f77 INTERNAL ERROR:  /usr/lib32/cmplrs/be returned non-zero status 1
    DBGCF_FLGS += -DEBUG:trap_uninitialized:verbose_runtime
  endif
  CFL_FLGS = -$(ABI) -mips$(MIPS)
  # do we need to change all these from defaults?
  # -LNO:opt=0 works around bug in f77 7.3.1.3m loop nest optimizer that
  # has trouble with loops around calls to subroutines with alternate returns.
  # There seems to be a small improvement in speed when this is turned off too.
  OPTF_FLGS += -LNO:opt=0 \
		-OPT:const_copy_limit=15100 \
		-OPT:memcpy_cannot_overlap=ON \
		-OPT:bcopy_cannot_overlap=ON \
		-OPT:memmove_cannot_overlap=ON
  OPTCF_FLGS = -O3 -OPT:Olimit=8000
  # warning 47 : This module contains branch instruction(s) that might
  # degrade performance on an older version (rev. 2.2) R4000 processor.
  # warning 158 : Expecting MIPS3 objects: ... MIPS4.
  L_FLGS += -rdata_shared -DEBUG:error=158 -woff 47 -woff 131
  MP_FLGS = -mp -MP:old_mp=OFF
  # Use the system so_locations to indicate where in memory the shared object
  # should be placed.  rld will move most of these as they all want to be in
  # the same place.
  # If we don't give any indication, ld Version 7.30 on hpc2 places them at
  # 0x5fe60000 which is in the middle of our heap space.
  # If we use -update_registry, the build fails when the required space of a
  # DSO increases.  rqsall bus errors and seg faults.  
  # The time rld spends moving things seems only of the order of 0.1s anyway.
  ifeq ($(ABI),n32)
    D_FLGS += -check_registry /usr/lib32/so_locations
  else
    D_FLGS += -check_registry /usr/lib64/so_locations
  endif
  # Specify soname as otherwise it is derived from the temporary file
  D_FLGS += -soname $(patsubst %.build,%,$(notdir $1))
  ifneq ($(DYNAMIC),false)# dynamic linking
    # dso symbols not remembered
    ELFLAGS += -rpath $(CURDIR)/$(DSO_DIR) -ignore_unresolved
  endif
#  ifeq ($(DYNAMIC),false)
#    OPT_FFLAGS += -IPA
#    OPT_LFLAGS += -IPA
#  endif
  olist_args = -objectlist $1
  # Use the BSD timers
  CPPFLAGS += -DBSD_TIMERS
endif
ifeq ($(SYSNAME),Linux)
  CC = gcc
  FC = g77
  OPTCF_FLGS =# Use separate flags for fortran and c
  olist_args = $1
#   ifneq ($(filter $(INSTRUCTION),ia64 x86_64 i686),)# ia64, x86_64, or i686
  ifneq ($(filter $(INSTRUCTION),ia64),)# ia64 only
    # Use Intel compilers if available
    # (icc -V sends output to STDERR and exits with error).
    ifneq (,$(shell icc -V 2>&1 | grep -i intel))
      CC = icc
    endif
    ifneq (,$(shell ifort -V 2>&1 | grep -i intel))
      FC = ifort
    endif
  endif
  ifeq ($(MACHNAME),x86_64)
    # pathf90 can't handle %VAL(ARRAY(index))
    # It seems to just pass as if %VAL were not there.
#     # pathf90 sends output to STDERR but exits success
#     ifneq (,$(shell pathf90 --version 2>/dev/null && echo 1))
#       FC = pathf90
#     else
    # gcc and Portland fortran as used by AMD
    ifneq (,$(shell pgf90 -V 2>/dev/null))
      FC = pgf90
      # pgcc is slightly slower for 21d and 8113 but slightly faster for b211u
      ifneq ($(MP),false)
        # but for multithreading use a compiler with openmp support.
        # (I don't think we currently use openmp in C but may do one day.)
        CC = pgcc
      endif
    endif
#     endif
  endif
  ifeq ($(MACHNAME),ppc64)
    # use xlf if available
    ifneq (,$(shell xlf95 2>/dev/null))
      # _r is required for linking with essl (but not for compiling)
#       ifeq ($(MP),false)
#         FC = xlf95
#       else
        FC = xlf95_r
#       endif
    endif
    ifneq (,$(shell xlc_r 2>/dev/null))
      CC = xlc_r
    endif
# use fdlibm pow:
# xlf95_r -qsmp=omp -q64 -L/usr/X11R6/lib generated/ppc64-64-linux-mt/object.list /home/tomlinso/cmiss/inear_solvers/lib/ppc64-64-linux/libsolver-mt.a /home/tomlinso/cmiss/linear_solvers/lib/ppc64-64-linux/libxblas-mt.a /home/tomlinso/cmiss/linear_solvers/lib/ppc64-64-linux/libsuperlu-mt.a /home/tomlinso/cmiss/linear_solvers/lib/ppc64-64-linux/libumfpack-4.0.a /home/tomlinso/cmiss/perl_interpreter/lib/ppc64-64-linux/libperlinterpreter.a /usr/lib64/liblapack.a /usr/lib64/libblas.a /usr/lib64/libcblas.a /usr/lib64/libatlas.a ~/fdlibm/w__pow.o ~/fdlibm/e_pow.o -ldl -lcrypt -lpthread -o bin/ppc64-linux/cm64-mt.link
  endif

  ifeq ($(CC),gcc)# gcc
    C_FLGS += -pipe
#    ELFLAGS += -pipe
    # Position independent code is actually only required for objects
    # in shared libraries but debug version may be built as shared libraries.
    DBGC_FLGS += -fPIC
    ifeq ($(filter $(INSTRUCTION),i686 ia64),)# not i686 nor ia64
      C_FLGS += -m$(ABI)
    endif
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        C_FLGS += -march=nocona
      endif
    endif
#     ifeq ($(filter-out i%86,$(MACHNAME)),)
#       ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
#         C_FLGS += -march=pentium4# -mfpmath=sse
#       endif
#     endif
    ifeq ($(MACHNAME),ppc64)
      C_FLGS += -mcpu=power5
    endif
    DBGC_FLGS += -O0 -fbounds-check
    # ? -ffastmath
    OPTC_FLGS = -O3# -funroll-all-loops
    ifeq ($(PROF),false)
      ifneq ($(filter $(INSTRUCTION),i686 x86_64),)# i686 or x86_64
        OPTC_FLGS += -momit-leaf-frame-pointer
      endif
    else
      C_FLGS += -g -pg# -fprofile-arcs -ftest-coverage
    endif
  endif
  ifeq ($(FC),g77)
    # -fstatck-check
    F_FLGS += -pipe -fno-second-underscore -Wunused
    # Position independent code is actually only required for objects
    # in shared libraries but debug version may be built as shared libraries.
    DBGF_FLGS += -fPIC
    ifeq ($(filter $(INSTRUCTION),i686 ia64),)# not i686 nor ia64
      F_FLGS += -m$(ABI)
      ELFLAGS += -m$(ABI)
    endif
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        F_FLGS += -march=nocona
      endif
    endif
#     ifeq ($(filter-out i%86,$(MACHNAME)),)
#       F_FLGS += -malign-double
#       ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
#         C_FLGS += -march=pentium4# -mfpmath=sse
#       endif
#     endif
    ifeq ($(MACHNAME),ppc64)
      F_FLGS += -mcpu=power5
    endif
    DBGF_FLGS += -O0 -fbounds-check
    OPTF_FLGS = -O3 -Wuninitialized# -funroll-all-loops
    ifeq ($(PROF),false)
      ifneq ($(filter $(INSTRUCTION),i686 x86_64),)# i686 or x86_64
        OPTF_FLGS += -momit-leaf-frame-pointer
      endif
    else
      F_FLGS += -g -pg# -fprofile-arcs -ftest-coverage
      ELFLAGS += -pg
    endif
  endif
  ifeq ($(FC),pathf90)
    F_FLGS  = -pipe -fno-second-underscore \
	-ff2c-abi /opt/pathscale/etc/f2c-abi/acml-2.0
    # Bounds checking does not appear to understand index ranges (-NIM:NIM)
#    DBGF_FLGS += -C# bounds checking
    OPTF_FLGS = -O3 -LNO:opt=0# -OPT:Ofast -LNO:prefetch=0
  endif
  ifeq ($(CC),icc)
    # Turn on remarks.
    # No headers for C externals because they are called from FORTRAN.
    # remark #1419: external declaration in primary source file
    # remark #891 or 981: operands are evaluated in unspecified order
    C_FLGS += -w2 -wd1419# -Wall
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        C_FLGS += -xP# for sse3 (90nm Pentium 4 series)
      else
        C_FLGS += -xW# Pentium4 compatible (?sse2)
      endif
    endif
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
        C_FLGS += -xN# for Pentium 4
      endif
    endif
    DBGC_FLGS += -O0
    OPTC_FLGS = -O3 -ansi_alias
  endif
  ifeq ($(FC),ifort)
    # turn on warnings,
    # suppress messages about non-standard fortran (including REAL*8,
    # more than 19 continuation lines),
    # suppress comment messages (including obsolescent alternate return),
    # avoid beyond 72 column warnings which are given even for comments,
    # auto arrays (as well as scalars), no aliasing,
    # temporary arrays on stack if possible, no preprocessing.
    F_FLGS += -W1 -w95 -cm -80 -auto -stack_temps# -fpp0
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        F_FLGS += -xP# for sse3 (90nm Pentium 4 series)
      else
        F_FLGS += -xW# Pentium4 compatible (?sse2)
      endif
    endif
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
        F_FLGS += -xN# for Pentium 4
      endif
    endif
    # Possibly also -nbs treat backslash as a normal character (not escape).
# This doesn't seem to do anything
#     ifeq ($(ABI),64)
#       CF_FLGS += -size_lp64
#     endif
    DBGF_FLGS += -traceback -O0 -check format# -check bounds
# This doesn't work around this bug with check bounds (with openmp atleast):
#     ifneq ($(DEBUG),false)
#       ELFLAGS += -Wl,-G0
#     endif
    # ld: bin/ia64-linux/cm-mt-debug.link: short data segment overflowed (0x59d9c0 >= 0x400000)
    # ld: can't relax section: No such file or directory
    ifneq ($(MACHNAME),ia64)
      DBGF_FLGS += -check bounds
    endif
    OPTF_FLGS = -O3
# I can't see any improvement with -fno-alias so
# it doesn't seem worth the risk.
#     OPTF_FLGS += -fno-alias
    MP_FLGS = -openmp
    ELFLAGS += -nofor_main
# ifort version 8.0 with -ipo
# IPO: performing multi-file optimizations
# fortcom(22486): unaligned access to 0x2000000009400cbe, ip=0x4000000000dc4750
# fortcom(22486): unaligned access to 0x2000000009400cd6, ip=0x4000000000dc4750
# fortcom(22486): unaligned access to 0x2000000009400cee, ip=0x4000000000dc4750
# fortcom(22486): unaligned access to 0x2000000009400d06, ip=0x4000000000dc4750
# ifort: error: Fatal error in /opt/intel_fc_80/bin/fortcom, terminated by segmentation violation
# ifort: error: problem during multi-file optimization compilation (code 1)
# ifort: error: empty multi-file optimizations object list
# ifort version 8.1 with -ipo
# ld: can't relax section: No such file or directory
#    OPTCFE_FLGS = -ipo
#    olist_args = @$1
  endif
  ifeq ($(CC),pgcc)
    OPTC_FLGS = -fastsse# -Mipa=fast
  endif
  ifeq ($(FC),pgf90)
    # All functions must be compiled with ipa in pgi52 for the linker to use
    # the info.
    OPTF_FLGS = -fastsse# -Mipa=fast
    MP_FLGS = -mp
    ELFLAGS += -Mnomain -fastsse# -Mipa=fast
#    # Making objects available to pgf90 for ipa.
#    # Strip linker script commands.
#    olist_args = `perl -wne '/^(INPUT\(|\))\n/ or print' $1`
  endif
  ifeq ($(filter-out xlc%,$(CC)),)# xlc* C compiler
    CFLAGS += -qinfo=gen:ini:por:pro:trd:tru:use
    C_FLGS += -q$(ABI) -qarch=auto -qhalt=e
    ifeq ($(ABI),64)
      # 1506-743 (I) 64-bit portability: possible change of result through conversion ...
      # These don't seem to serious.  Truncations are reported separately.
      C_FLGS += -qwarn64 -qsuppress=1506-743
    endif
    # -qinitauto for C is bytewise: 7F gives large integers.
    # -qinitauto gives linker message (and seg fault on cm start up)
    # /usr//bin/ld: BFD 2.15.90.0.1.1 20040303 (SuSE Linux) assertion fail ../../bfd/elf64-ppc.c:6272
    DBGC_FLGS += -qfullpath -C -qflttrap=inv:en -qinitauto=7F
    # OPT flags like these would be nice too,
    # but they seem to cause compilation problems
    # -qipa -qlibansi -qlibessl -Q
    # 1500-036 (I) The NOSTRICT option (default at OPT(3)) has the potential
    # to alter the semantics of a program.  Please refer to documentation on
    # the STRICT/NOSTRICT option for more information.
    # -qhot -qcache=auto sometimes makes things about 5% slower
    OPTC_FLGS += -O3 -qstrict -qfloat=fltint:rsqrt -qmaxmem=12000 -qtune=auto -qcache=auto -qnoignerrno
    # for profiling (really need to profile blas etc too though)
#   OPTCFE_FLGS += -pg -g
    # for trailing _ on fortran symbols
    CPPFLAGS += -Dunix
  endif
  ifeq ($(filter-out xlf%,$(FC)),)# xlf* fortran compiler
    F_FLGS += -q$(ABI) -qarch=auto -qhalt=e -qfixed -qextname -qrecur
    ELFLAGS += -q$(ABI)

#     CFLAGS += -qinfo=gen:ini:por:pro:trd:tru:use
#     C_FLGS += -q$(ABI) -qarch=auto -qhalt=e
    ifeq ($(ABI),64)
      # 1506-743 (I) 64-bit portability: possible change of result through conversion ...
      # These don't seem to serious.  Truncations are reported separately.
      # FD_SET in sys/time.h does this
      F_FLGS += -qwarn64 -qsuppress=1506-743
#       C_FLGS += -qwarn64 -qsuppress=1506-743
#       # arbitrary < max heap size 0x06FFFFFFFFFFFFF8
#       ELFLAGS += -bmaxdata:0x0001000000000000
    else
#       # max possible is 0x80000000.  I don't think we need space for mmap etc
#       ELFLAGS += -bmaxdata:0x80000000/dsa
    endif
    ifeq ($(DEBUG),false)
      MP_FLGS = -qsmp=omp
    else
      MP_FLGS = -qsmp=omp:noopt
    endif
    # lapack's dlamch performs an underflow so we don't check that.
    # -qextchk doesn't seem to like %VAL(0) at least
    # 7FF7FFFF is a large integer or NaNQ real*4 or NaNS real*8
    # -qinitauto for C is bytewise: 7F gives large integers.
    # -qinitauto gives linker message (and seg fault on cm start up)
    # /usr//bin/ld: BFD 2.15.90.0.1.1 20040303 (SuSE Linux) assertion fail ../../bfd/elf64-ppc.c:6272
    DBGF_FLGS += -qfullpath -C -qflttrap=inv:en# -qinitauto=7FF7FFFF
#     DBGC_FLGS += -qfullpath -C -qflttrap=inv:en -qinitauto=7F
    # OPT flags like these would be nice too,
    # but they seem to cause compilation problems
    # -qipa -qlibansi -qlibessl -Q
    # 1500-036 (I) The NOSTRICT option (default at OPT(3)) has the potential
    # to alter the semantics of a program.  Please refer to documentation on
    # the STRICT/NOSTRICT option for more information.
    # -qhot -qcache=auto sometimes makes things about 5% slower
    OPTF_FLGS += -O3 -qstrict -qfloat=fltint:rsqrt -qmaxmem=12000 -qtune=auto
    # This sometimes makes things about 5% slower, and with version 9.1 on AIX
    # causes a SIGSEGV with example 118 (and others: see GLOBALH below) and in
    # setting up XID in EVELEC in some of Travis's work with 64-bit binaries.
    # (-qcache=auto only has an affect if -qhot is used.)
#    OPTF_FLGS += -qhot -qcache=auto
    # for profiling (really need to profile blas etc too though)
#   OPTCFE_FLGS += -pg -g
#  olist_args = -f $1
  endif
  # --export-dynamic (as used by perl 5.8.6) exports symbols for perl binary
  # modules but these symbols then obscure any from a shared libperl
  # subsequently loaded
#  ELFLAGS += -Wl,--export-dynamic
  ifneq ($(DYNAMIC),false)# dynamic linking
    ELFLAGS += -Wl,-rpath,$(CURDIR)/$(DSO_DIR)
  endif
  # Avoid versioning problems with libgcc_s by linking statically.

  # libgcc2.c from gcc 3.4.4 says:
  # In addition to the permissions in the GNU General Public License, the
  # Free Software Foundation gives you unlimited permission to link the
  # compiled version of this file into combinations with other programs,
  # and to distribute those combinations without any restriction coming
  # from the use of this file.

  # (With dynamic version, should copy libgcc_s.so.N if copying libstdc++.so.N)
  ELFLAGS += -static-libgcc
  # Work-around for linking objects built against glibc headers before 2.3.5.
#  ELFLAGS += -Wl,--warn-unresolved-symbol

#   ifeq ($(MACHNAME),ppc64)
#     ifeq ($(ABI),64)
#       L_FLGS += -m elf64ppc
#     # note for ABI=32 -m elf64ppclinux is probably appropriate (but default)
#     endif
#   endif
  # Use the BSD timers
  CPPFLAGS += -DBSD_TIMERS
  # Using ld here with gcc-3.3.3 gives the following error:
  # /product/cmiss/cross-compile/i386-glibc21-linux/bin/../lib/gcc-lib/i686-pc-linux-gnu/3.3.3/../../../../i686-pc-linux-gnu/bin/ld: dso/i686-linux-debug/cm.link: hidden symbol `__dso_handle' in /product/cmiss/cross-compile/i386-glibc21-linux/bin/../lib/gcc-lib/i686-pc-linux-gnu/3.3.3/crtbegin.o is referenced by DSO
  DSO_LINK = $(CC)
endif
ifeq ($(SYSNAME),win32)
  FC = g77
  F_FLGS += -fno-second-underscore
  OPTCF_FLGS = -O2
#    FC = pgf90
#    OPT_FFLAGS = -fast
#    EXE_LINK := pgf90 -Mnomain
  ifneq ($(DYNAMIC),false)# dynamic linking
    ELFLAGS += -Wl,-rpath,$(CURDIR)/$(DSO_DIR)
  else
    # For portability?? don't rely on any shared objects being available.
    # With this flag all libraries should be in ARCHIVES so that
    # the executable is updated if these libraries change - not done yet.
    ELFLAGS += -Wl,-static
  endif
  # Use the ANSI C timers
  CPPFLAGS += -DANSI_C_TIMERS
  olist_args = $1
endif
ifeq ($(SYSNAME),SunOS)
  FC = f90
  # READMEs/fortran_95 says turn this option on so that a string is
  # padded with blanks when it is assigned a shorter string.
  F_FLGS += -Qoption f90comp -fephoenix
  ELFLAGS += -xlic_lib=sunperf
  DBGF_FLGS += -C# check subscript bounds
  # need -xarch=native after -fast
  OPTCFE_FLGS += -fast $(CFE_FLGS)
  ifeq ($(ABI),64)
    CFE_FLGS += -xarch=native64
  endif
  MP_FLGS = -openmp
  # Use the POSIX timers
  CPPFLAGS += -DPOSIX_TIMERS
endif
ifeq ($(SYSNAME),AIX)
  ifeq ($(MP),false)
    FC = xlf95
    CC = xlc
  else
    FC = xlf95_r
    CC = xlc_r
  endif
  F_FLGS += -qfixed -qextname -qrecur
  CFLAGS += -qinfo=gen:ini:por:pro:trd:tru:use
  # Do we want ELFLAGS += -bstatic? can't find half the symbols with this
  # -brtl is for loading Perl 5.8.0 dynamic libraries.
  ELFLAGS += -bbigtoc -brtl $(addprefix -bE:,$(EXPORTS))
# -brename:dtrsv,dtrsv_ -brename:dtrsm,dtrsm_ -brename:dgemv,dgemv_ -brename:dgemm,dgemm_ -brename:dnrm2,dnrm2_ -brename:dscal,dscal_ -brename:dger,dger_ -brename:dtrmv,dtrmv_ -brename:dcopy,dcopy_ -brename:dtrmm,dtrmm_ -brename:drot,drot_ -brename:dswap,dswap_ -brename:ddot,ddot_ -brename:idamax,idamax_ -brename:dsymv,dsymv_ -brename:daxpy,daxpy_ -brename:dsyr2k,dsyr2k_ -brename:dsyr2,dsyr2_ -brename:dsyrk,dsyrk_ -brename:dasum,dasum_ -brename:dsymm,dsymm_ -brename:dsyr,dsyr_
# -brename:dtrsv_,dtrsv -brename:dtrsm_,dtrsm -brename:dgemv_,dgemv -brename:dgemm_,dgemm -brename:dnrm2_,dnrm2 -brename:dscal_,dscal -brename:dger_,dger -brename:dtrmv_,dtrmv -brename:dcopy_,dcopy -brename:dtrmm_,dtrmm -brename:drot_,drot -brename:dswap_,dswap -brename:ddot_,ddot -brename:idamax_,idamax -brename:dsymv_,dsymv -brename:daxpy_,daxpy -brename:dsyr2k_,dsyr2k -brename:dsyr2_,dsyr2 -brename:dsyrk_,dsyrk -brename:dasum_,dasum -brename:dsymm_,dsymm -brename:dsyr_,dsyr

  ifneq ($(DYNAMIC),false)# dynamic linking
    ELFLAGS += -blibpath:$(CURDIR)/$(DSO_DIR):/usr/lib
  else
    ELFLAGS += -bnolibpath
  endif
  # May want -qsrcmsg
  CFE_FLGS += -q$(ABI) -qarch=auto -qhalt=e
  L_FLGS += -b$(ABI)
  D_FLGS = -G -bexpall -bnoentry
  # ? mass xlopt m (in /usr/lpp/xlf/lib)
  DSO_LIBS = xlf90 c
  ifeq ($(ABI),32)
    # Without -bmaxdata, the only one 256M virtual segment is available for
    # data.
    # In AIX 5.3, 0xAFFFFFFF is the largest value we can use here and still
    # use global shared libraries. (see aixprggd/genprogc/lrg_prg_support.htm)
    # However, 0xAFFFFFFF/dsa causes the system to crash on loading of perl
    # modules (File::Find and Time::HiRes atleast).  0x80000000 seems to work.
    # dsa allows segments to be allocated dynamically for shmat/mmap or data
    # as required.
    ELFLAGS += -bmaxdata:0x80000000/dsa
  else
    # 1506-743 (I) 64-bit portability: possible change of result through conversion ...
    # These don't seem to serious.  Truncations are reported separately.
    # FD_SET in sys/time.h does this
    CF_FLGS += -qwarn64 -qsuppress=1506-743
    # It seems that somewhere between AIX 5.1 and 5.3 the kernel loader
    # started modifying a process's soft data resource limit to match to match
    # its maxdata value (if non-zero).  As 32-bit applications need a non-zero
    # maxdata value to access more than 256M of data many applications
    # (including perl) will cause the data limit to be lowered to a 32-bit
    # addressable value.  As cmiss is likely to be a child of such 32-bit
    # processes, to access more than 32-bit addressable memory, it either
    # needs to raise its data limit or use its own maxdata value.
    # max heap size is 0x06FFFFFFFFFFFFF8
    # Arbitrary.  0x0000100000000000 should provide ~16TB.
    ELFLAGS += -bmaxdata:0x0000100000000000
  endif
  ifeq ($(DEBUG),false)
    MP_FLGS = -qsmp=omp
  else
    MP_FLGS = -qsmp=omp:noopt
  endif
  # lapack's dlamch performs an underflow so we don't check that.
  # -qextchk doesn't seem to like %VAL(0) at least
  # Should -qflttrap=nans be used as well or instead of -qflttrap=inv:en?
  DBGCF_FLGS += -qfullpath -C -qflttrap=inv:en
  # 7FF7FFFF is a large integer or NaNQ real*4 or NaNS real*8
  # -qinitauto for C is bytewise: 7F gives large integers.
  DBGF_FLGS += -qinitauto=7FF7FFFF
  DBGC_FLGS += -qinitauto=7F
  # OPT flags like these would be nice too,
  # but they seem to cause compilation problems
  # -qipa -qlibansi -qlibessl -Q
  # 1500-036 (I) The NOSTRICT option (default at OPT(3)) has the potential
  # to alter the semantics of a program.  Please refer to documentation on
  # the STRICT/NOSTRICT option for more information.
  OPTCF_FLGS = -O3 -qstrict -qfloat=fltint:rsqrt -qmaxmem=12000 -qtune=auto
  # This sometimes makes things about 5% slower, and with version 9.1
  # causes a SIGSEGV with example 118 (and others: see GLOBALH below) and in
  # setting up XID in EVELEC in some of Travis's work with 64-bit binaries.
  # (-qcache=auto only has an affect if -qhot is used.)
#  OPTF_FLGS += -qhot -qcache=auto
  OPTC_FLGS += -qnoignerrno
  # for profiling (really need to profile blas etc too though)
#   OPTCFE_FLGS += -pg -g
  olist_args = -f $1
  # Use the BSD timers
  CPPFLAGS += -DBSD_TIMERS
endif

PERL = perl -w

# This returns an empty string if not found
searchdirs = $(firstword $(wildcard $(addsuffix /$(strip $2),$1)))
# This still returns the name of the desired file if not found and so is
# useful for error checking and reporting.
searchdirsforce = $(firstword $(wildcard $(addsuffix /$(strip $2),$1)) $2)

# Check that call function works (for searchdirs, olist_args, etc.)
ifeq ($(call olist_args,test),)
  $(error call function not available.  Use GNU make version 3.78 or newer)
endif

REV_TIME_SOURCE = revision_time.c

SOURCE_PATTERNS = \
	cm.c \
	fe00/*.f \
	fe01/*.f \
	fe02/*.f \
	fe03/*.f \
	fe04/*.f \
	fe05/*.f \
	fe06/*.f \
	fe07/*.f \
	fe08/*.f \
	fe09/*.f \
	fe11/*.f \
	fe12/*.f \
	fe13/*.f \
	fe14/*.f \
	fz16.f \
	fz18.f \
	fe19/*.f \
	fe20/*.f \
	fe21/*.f \
	fe22/*.f \
	fe23/*.f \
	fe24/*.f \
	fe25/*.f \
	fe26/*.f \
	fe27/*.f \
	fe28/*.f \
	fe29/*.f \
	fe30/*.f \
	fe31/*.f \
	fe40/*.f \
	fe50/*.f \
	fe60/*.f \
	fz70/*.f \
	fe90/*.f \
	feinterpreter/*.f \
	fecutils.c \
	fefutils/*.f \
	data_objects/constant/*.f \
	data_objects/core/*.f \
	data_objects/equation/*.f \
	data_objects/field/*.f \
	data_objects/maths/*.f \
	example_path.c \
	message.c \
	debug.c \
	fetriangle.c \
	feuser/*.f \
	feuser_cell/*.f \
	$(REV_TIME_SOURCE) \
	fzminos.f \
	fzsocket.f

ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
  SOURCE_PATTERNS += fe10_sgi/*.f
else
  ifeq ($(filter-out xlf%,$(FC)),)# xlf* fortran compiler
    SOURCE_PATTERNS += fe10_ibm/*.f
  else# Linux
    SOURCE_PATTERNS += fe10_linux/*.f
  endif
endif

ARCHIVES :=# archives required by the main executable
# MODULE_ARCHIVES := #archives required only by one (or few) module
REQUIRED_BEGIN_OBJECTS=
REQUIRED_END_OBJECTS=

QUADPACK_ARCHIVE := $(call searchdirs, \
		 $(ARCHIVE_DIRS) $(CMISS_ROOT)/quadpack/lib/$(LIB_ARCH_DIR), \
		 libquadpack$(DEBUG_SUFFIX).a)
ifeq ($(QUADPACK_ARCHIVE),)
  SOURCE_PATTERNS += fzquadpack.f
else
  ARCHIVES += $(QUADPACK_ARCHIVE)
endif

SQP_ARCHIVE := $(call searchdirs, \
		 $(ARCHIVE_DIRS) $(CMISS_ROOT)/sqp/lib/$(LIB_ARCH_DIR), \
		 libsqp$(DEBUG_SUFFIX).a)
ifeq ($(SQP_ARCHIVE),)
  SOURCE_PATTERNS += fzsqp.f
else
  ARCHIVES += $(SQP_ARCHIVE)
endif

LSODA_ARCHIVE := $(call searchdirs, \
	           $(ARCHIVE_DIRS) $(CMISS_ROOT)/lsoda/lib/$(LIB_ARCH_DIR), \
		   liblsoda$(MT_SUFFIX)$(DEBUG_SUFFIX).a)
ifeq ($(LSODA_ARCHIVE),)
  SOURCE_PATTERNS += fzlsoda.f
else
  ARCHIVES += $(LSODA_ARCHIVE)
endif

SOLVER_ARCHIVES :=
SOLVER_ARCHIVE_DIRS := \
	$(ARCHIVE_DIRS) $(CMISS_ROOT)/linear_solvers/lib/$(LIB_ARCH_DIR)

SOLVER_ARCHIVE := $(call searchdirs,$(SOLVER_ARCHIVE_DIRS),libsolver$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX).a)
ifeq ($(SOLVER_ARCHIVE),)
  SOURCE_PATTERNS += fzsolver.f
else
  SOURCE_PATTERNS += fesolver/*.f
  SOLVER_ARCHIVES += $(SOLVER_ARCHIVE)
endif

XBLAS_ARCHIVE := $(call searchdirs,$(SOLVER_ARCHIVE_DIRS),libxblas$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX).a)
ifeq ($(XBLAS_ARCHIVE),)
  SOURCE_PATTERNS += fzxblas.f
else
  SOLVER_ARCHIVES += $(XBLAS_ARCHIVE)
endif

SUPERLU_ARCHIVE := $(call searchdirs,$(SOLVER_ARCHIVE_DIRS),libsuperlu$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX).a)
ifeq ($(SUPERLU_ARCHIVE),)
  SOURCE_PATTERNS += fzsuperlu.c
else
  SOLVER_ARCHIVES += $(SUPERLU_ARCHIVE)
endif

UMFPACK_ARCHIVE := $(call searchdirs,$(SOLVER_ARCHIVE_DIRS),libumfpack-4.0$(DEBUG_SUFFIX)$(PROF_SUFFIX).a)
ifeq ($(UMFPACK_ARCHIVE),)
  SOURCE_PATTERNS += fzumfpack.c
else
  SOLVER_ARCHIVES += $(UMFPACK_ARCHIVE)
endif

AMG_ARCHIVE := $(call searchdirs,$(SOLVER_ARCHIVE_DIRS),libamg$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX).a)
ifeq ($(AMG_ARCHIVE),)
  SOURCE_PATTERNS += fzamg.f
else
  SOLVER_ARCHIVES += $(AMG_ARCHIVE)
endif

BMG_ARCHIVE := $(call searchdirs,$(SOLVER_ARCHIVE_DIRS),libboxmg$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX).a)
ifeq ($(BMG_ARCHIVE),)
  SOURCE_PATTERNS += fzbmg.f
else
  SOLVER_ARCHIVES += $(BMG_ARCHIVE)
endif

ARCHIVES += $(SOLVER_ARCHIVES) 

UNEMAP_ARCHIVE := $(call searchdirs,$(ARCHIVE_DIRS),libunemap.a)
ifeq ($(UNEMAP_ARCHIVE),)
  SOURCE_PATTERNS += fzunemap.c
else
  ARCHIVES += $(UNEMAP_ARCHIVE)
  SOURCE_PATTERNS += feunemap.c
endif

# Set up the cellml archives and paths
XERCES_ARCHIVE_NAME := xerces-c$(DEBUG_SUFFIX)
CELLML_ARCHIVE_NAME := cellml$(DEBUG_SUFFIX)
CMISS_CELLML_ARCHIVE_NAME := CMISSCellML$(DEBUG_SUFFIX)
# A check for the existance of CMISS_CELLML_ARCHIVE is an additional check
# when using CellML
ifeq ($(USE_CELLML),true)
  CELLML_ARCHIVE_DIR := $(CELLML_ROOT)/lib/$(LIB_ARCH_DIR)
  CMISS_CELLML_ARCHIVE := \
	$(wildcard $(CELLML_ARCHIVE_DIR)/lib$(CMISS_CELLML_ARCHIVE_NAME).a)
  ifeq ($(CMISS_CELLML_ARCHIVE),)
    USE_CELLML=false
  endif
endif
ifeq ($(USE_CELLML),true)
  CELLML_ARCHIVES := $(CMISS_CELLML_ARCHIVE) \
	$(patsubst %,$(CELLML_ARCHIVE_DIR)/lib%.a, \
	  $(CELLML_ARCHIVE_NAME) $(XERCES_ARCHIVE_NAME))
  ifeq ($(SYSNAME),Linux)
#    DYNAMIC_LIBRARIES += pthread

# For the dynamic version, copying the shared lib (with the version suffix) to
# DSO_DIR would reduce the executable size by almost 2M.
# But g++ 3.4.3 uses libstdc++.so.6 and g++ 3.3.6 uses libstdc++.so.5 so
# linking and compiling of cellml needs to be done with the same compiler.
# Should the libstdc++.a be with libcellml or should libcellml be a shared
# library?

# libstdc++ is distributed under GPL with the so-called `runtime exception'
# (from gcc-3.4.4/libstdc++-v3/docs/html/17_intro/license.html):
# if other files instantiate
# templates or use macros or inline functions from this file, or you compile
# this file and link it with other files to produce an executable, this
# file does not by itself cause the resulting executable to be covered by
# the GNU General Public License.

    # Need a static libstdc++.
    # Use a stdc++ archive file from the compiler used to build cellml.

    # gcc 3.3.3 (cross-compile) doesn't seem able to link this into the shared
    # object on i686 and gcc 4.0.1-4ubuntu9 on x86_64 has libstdc++ built
    # without -fPIC so link into the executable.
# /product/cmiss/cross-compile/i386-glibc21-linux/bin/../lib/gcc-lib/i686-pc-linux-gnu/3.3.3/../../../../i686-pc-linux-gnu/bin/ld: /tmp/cmiss/!hpc!cmiss!cm!dso!i686-linux-debug!libfecellml.so: undefined versioned symbol name std::time_get_c@@GLIBCPP_3.2
# /product/cmiss/cross-compile/i386-glibc21-linux/bin/../lib/gcc-lib/i686-pc-linux-gnu/3.3.3/../../../../i686-pc-linux-gnu/bin/ld: failed to set dynamic section sizes: Bad value
    ARCHIVES += $(CELLML_ARCHIVE_DIR)/libstdc++.a
  endif
  SOURCE_PATTERNS += fecellml/*.f
else
  CELLML_ARCHIVES :=
  SOURCE_PATTERNS += fzcellml/*.f
endif

# CellML library depends on Perl interpreter, so need this library after cellml
PERL_INTERPRETER_ARCHIVE := $(call searchdirsforce, $(ARCHIVE_DIRS) \
                     $(CMISS_ROOT)/perl_interpreter/lib/$(LIB_ARCH_DIR), \
              libperlinterpreter$(DEBUG_SUFFIX).a)
ARCHIVES += $(PERL_INTERPRETER_ARCHIVE)

ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
  DYNAMIC_LIBRARIES = \
	complib.sgimath$(MP_SUFFIX)
  ifneq ($(DEBUG),false)
#     ARCHIVES += $(foreach name, \
# 	lapack$(MT_SUFFIX)$(DEBUG_SUFFIX) blas$(MT_SUFFIX)$(DEBUG_SUFFIX), \
# 	$(call searchdirsforce, $(SOLVER_ARCHIVE_DIRS),lib$(name).a))
    DYNAMIC_LIBRARIES += \
	exc \
	malloc_ss \
	fpe
  else
    DYNAMIC_LIBRARIES += \
	fastm
  endif
endif
ifeq ($(SYSNAME),Linux)
  # for perl; do we also want nsl and util?
  DYNAMIC_LIBRARIES = dl crypt pthread
  ifdef MKLPATH
    # Use Intel libraries
    # For MKL version 7.0, threading in blas and lapack routines is
    # controlled through the OMP_NUM_THREADS environment variable.  The
    # default is 1 if the variable is not set.  Parallelization in MKL can
    # also be disabled by setting MKL_SERIAL to yes.
    MKL_LIBRARIES = mkl_lapack
    ifeq ($(MACHNAME),x86_64)
      MKL_LIBRARIES += mkl_em64t
    else
    ifeq ($(MACHNAME),ia64)
      MKL_LIBRARIES += mkl_ipf
    else
      MKL_LIBRARIES += mkl_ia32.a
    endif
    endif
    ARCHIVES += \
	$(patsubst %,$(MKL_ARCHIVE_DIR)/lib%.a, $(MKL_LIBRARIES)) \
    DYNAMIC_LIBRARIES += guide# pthread
  else
  ifdef ACMLPATH
    ifeq ($(FC),pgf90)
      ACML_ARCHIVE_DIR = $(ACMLPATH)/pgi64$(MP_SUFFIX)/lib
    else
      ACML_ARCHIVE_DIR = $(ACMLPATH)/gnu64/lib
      ifeq ($(FC),ifort)
        DYNAMIC_LIBRARIES += g2c
      endif
    endif
    ARCHIVES += $(ACML_ARCHIVE_DIR)/libacml.a
  else
  ifeq ($(FC),pgf90)
#  ifneq (,$(wildcard $(PGI)/linux86-64/5.2/lib/libacml.a))
    ARCHIVES += $(PGI)/linux86-64/5.2/lib/libacml.a
  else
  ATLAS_ARCHIVE := $(shell gcc -m$(ABI) -print-file-name=libatlas.a)
  ifneq (libatlas.a,$(ATLAS_ARCHIVE))
    ATLAS_ARCHIVE_DIR := $(dir $(ATLAS_ARCHIVE))
    # The system has a good blas library - assume lapack is good too if present
    ARCHIVES += \
	$(call searchdirsforce, $(ATLAS_ARCHIVE_DIR) $(SOLVER_ARCHIVE_DIRS), \
	  liblapack.a)
    # Atlas usually installs blas as libf77blas but
    # Gentoo renames this to libblas.
    ARCHIVES += \
	$(firstword $(wildcard $(ATLAS_ARCHIVE_DIR)libf77blas.a) \
		    $(ATLAS_ARCHIVE_DIR)libblas.a) \
	$(ATLAS_ARCHIVE_DIR)libcblas.a \
	$(ATLAS_ARCHIVE)

    ifeq ($(FC),ifort)
      DYNAMIC_LIBRARIES += g2c
#      ARCHIVE_DIRS += $(dir $(shell gcc --print-file-name=libg2c.a))
    endif
  else
  ifeq ($(filter-out xlf%,$(FC)),)# xlf* fortran compiler
#  ifneq (libessl.so,$(shell gcc -m$(ABI) -print-file-name=libessl.so))
    ARCHIVES += $(call searchdirsforce, $(SOLVER_ARCHIVE_DIRS), \
	          lapack$(MT_SUFFIX)$(DEBUG_SUFFIX).a)
    ifeq ($(MP),false)
      DYNAMIC_LIBRARIES += essl# mass
    else
      DYNAMIC_LIBRARIES += esslsmp# mass
    endif
  else
  ifneq ($(wildcard /usr/lib/libscs$(MP_SUFFIX).so),)# SGI math libraries
    DYNAMIC_LIBRARIES += scs$(MP_SUFFIX)
  else
    # Using cmiss libraries if present, otherwise system libraries.
    LAPACK_ARCHIVE = $(call searchdirs, $(SOLVER_ARCHIVE_DIRS),liblapack$(MT_SUFFIX)$(DEBUG_SUFFIX).a)
    ifeq ($(LAPACK_ARCHIVE),)
      DYNAMIC_LIBRARIES += lapack
    else
      ARCHIVES += $(LAPACK_ARCHIVE)
    endif
    BLAS_ARCHIVE = $(call searchdirs, $(SOLVER_ARCHIVE_DIRS),libblas$(MT_SUFFIX)$(DEBUG_SUFFIX).a)
    ifeq ($(BLAS_ARCHIVE),)
      DYNAMIC_LIBRARIES += blas
    else
      ARCHIVES += $(BLAS_ARCHIVE)
    endif
# acml 2.1.0 is significantly faster than pgi 5.2 libraries for examples
# 21d 8113 and e918 but slightly slower for b211v
#       ifdef PGI
#         ARCHIVE_DIRS += $(PGI)/linux86-64/5.2/lib
#       endif
  endif
  endif
  endif
  endif
#       STATIC_LIBRARIES += \
# 	goto$(MT_SUFFIX) \
# 	lapack$(MT_SUFFIX)
  endif
  endif
# I think these are probably in xlf libraries anyway
#   ifeq ($(filter-out xlf%,$(FC)),)# xlf* fortran compiler
#     ifeq ($(ABI),32)
#       STATIC_LIBRARIES += massvp4
#       ARCHIVE_DIRS += /opt/ibmcmp/xlmass/4.1/lib
#     else
#       STATIC_LIBRARIES += massvp4_$(ABI)
#       ARCHIVE_DIRS += /opt/ibmcmp/xlmass/4.1/lib$(ABI)
#     endif
#   endif
endif
ifeq ($(SYSNAME),win32)
  ARCHIVES += $(foreach name, lapack$(MT_SUFFIX) blas$(MT_SUFFIX), \
	$(call searchdirsforce, $(SOLVER_ARCHIVE_DIRS),lib$(name).a))
  DYNAMIC_LIBRARIES =
endif
ifeq ($(SYSNAME),SunOS)
  DYNAMIC_LIBRARIES = \
        socket \
	nsl \
        dl
endif
ifeq ($(SYSNAME),AIX)
  # example b221 runs 10% slower with essl than our lapack and blas.
  # daxpy is much slower with essl, but dgemm is much faster.
  # essl 4.1 is not a full lapack implementation.
  # mass 4.1 makes b221u,b211v,h3 60-70% faster.
  ARCHIVES += $(call searchdirsforce,$(SOLVER_ARCHIVE_DIRS), \
	        liblapack$(MT_SUFFIX).a)
#   ifneq ($(DEBUG),false)
#     ARCHIVES += $(call searchdirsforce, $(SOLVER_ARCHIVE_DIRS), \
# 		libblas$(MT_SUFFIX)$(DEBUG_SUFFIX).a)
#   else
#   STATIC_LIBRARIES += blas cblas atlas
    ifeq ($(MP),false)
      DYNAMIC_LIBRARIES += essl
    else
      DYNAMIC_LIBRARIES += esslsmp
    endif
#   endif
  DYNAMIC_LIBRARIES += mass
  # ld,bsd required for perl (C may be needed for versions earlier than 5.8.0).
  # m must be before bsd (and C) to use pow from m.
  # Others seg fault when called from perl.
  DYNAMIC_LIBRARIES += ld m pthreads bsd
#  DYNAMIC_LIBRARIES += ld bsd
endif

WORMHOLE_ARCHIVE := \
	$(call searchdirs, \
	  $(ARCHIVE_DIRS) $(CMISS_ROOT)/wormhole/lib/$(BIN_ARCH_DIR), \
	  libwormhole$(EXE_ABI_SUFFIX).a)
ifeq ($(WORMHOLE_ARCHIVE),)
  SOURCE_PATTERNS += fzwormhole.c
else
  ARCHIVES += $(WORMHOLE_ARCHIVE)
endif

GX_ARCHIVES :=
GX_ARCHIVE := $(call searchdirs, \
	        $(ARCHIVE_DIRS) $(CMISS_ROOT)/gx/lib/$(BIN_ARCH_DIR), \
	        libgx$(EXE_ABI_SUFFIX).a)
ifeq ($(GX_ARCHIVE),)
  SOURCE_PATTERNS += fzgx.f
else
  SOURCE_PATTERNS += fegx/*.f
  GX_ARCHIVES += $(GX_ARCHIVE)
  ifeq ($(SYSNAME),Linux)
    # For portability use a static set of Motif libraries (Mrm Xm)
    # consistent with those used to build gx.
    GX_ARCHIVE_DIR := $(dir $(GX_ARCHIVE))
    GX_ARCHIVES += $(patsubst %,$(GX_ARCHIVE_DIR)lib%.a, Mrm Xm)
    # Xm requires Xp, but libXp.so.6 is not in Ubuntu's ia32-libs amd64
    # package, so build in Xp for these systems.
    # libXp.so is more likely installed on the build machine than libXp.a.
    XP_ARCHIVE =
    ifeq ($(INSTRUCTION),i686)
      ifeq ($(DYNAMIC),false)
        XP_ARCHIVE = $(shell gcc -print-file-name=libXp.a)
        ifeq ($(XP_ARCHIVE),libXp.a) # libXp.a not in system library path
          XP_ARCHIVE = \
            $(call searchdirs, /usr/X11R6/lib$(ABI) /usr/X11R6/lib, libXp.a)
        endif
      endif
    endif
    ifneq ($(XP_ARCHIVE),)
      GX_ARCHIVES += $(XP_ARCHIVE)
    else
      DYNAMIC_LIBRARIES += Xp
    endif
# !!! Do we need to explicitly include libraries referenced by DSOs?
    # Xmu.so references Xt SM ICE Xext X11
    # Xt.so references X11 SM ICE
    X11_DSO := $(shell gcc -m$(ABI) -print-file-name=libX11.so)
    ifeq ($(X11_DSO),libX11.so) # libX11 not in system library path
      ELFLAGS += -L$(dir $(firstword \
	$(call searchdirs, /usr/X11R6/lib$(ABI) /usr/X11R6/lib, libX11.so) \
	X11LIBDIR/libX11.so))
    endif
    DYNAMIC_LIBRARIES += \
	Xmu \
	Xt \
	Xext \
	X11 \
	SM \
	ICE
  else 
    DYNAMIC_LIBRARIES += \
	Mrm \
	Xmu \
	Xm \
	Xt \
	X11
  endif
  ARCHIVES += $(GX_ARCHIVES)
endif
ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
   #For IRIX MP we are going to use ld so that we can hide symbols
    ifeq ($(ABI),n32)
      REQUIRED_BEGIN_OBJECTS += /usr/lib32/mips$(MIPS)/crt1.o
      REQUIRED_END_OBJECTS += /usr/lib32/mips$(MIPS)/crtn.o
    endif
    ifeq ($(ABI),64)
      REQUIRED_BEGIN_OBJECTS += /usr/lib64/mips$(MIPS)/crt1.o
      REQUIRED_END_OBJECTS += /usr/lib64/mips$(MIPS)/crtn.o
    endif
   DYNAMIC_LIBRARIES += ftn m c
   ifneq ($(MP),false)
      DYNAMIC_LIBRARIES += mp
   endif
endif

ifeq ($(USE_CELLML),true)
  ifeq ($(filter-out IRIX%,$(SYSNAME)),)
    DYNAMIC_LIBRARIES += C Cio c
    ifeq ($(ABI),n32)
      REQUIRED_BEGIN_OBJECTS += \
        /usr/lib32/c++init$(MP_SUFFIX).o
    endif
    ifeq ($(ABI),64)
      REQUIRED_BEGIN_OBJECTS += \
        /usr/lib64/c++init$(MP_SUFFIX).o
    endif
  endif
  ifeq ($(SYSNAME),AIX)
    DYNAMIC_LIBRARIES += dl pthread C c
  endif
endif
# libg2c contains code under LGPL.  Link dynamically so we can distribute the
# binary under the terms of this licence without needing to supply objects and
# build scripts.  Unfortunately this means that libg2c.so must be installed to
# run the binary but at least this library is widely available and has a
# stable interface.
# ifeq ($(DYNAMIC),false)
#   ifeq ($(SYSNAME),Linux)
#     ifeq ($(FC),g77)
#       # Use a static libg2c as the shared library is only installed
#       # with gcc if fortran support is enabled.
#       # Use gcc instead of g77 to link so that we can specify the static
#       # library explicitly but then need to also explicitly specify frtbegin
#       # and m.
#       EXE_LINK = gcc# or would $(CC) be better
#       ARCHIVES += $(foreach name, frtbegin g2c, \
# 	            $(shell $(FC) --print-file-name=lib$(name).a))
#       DYNAMIC_LIBRARIES += m
#     endif
#   endif
# endif

EXPORTS = $(wildcard $(patsubst %.a,%.exp, $(ARCHIVES)))

DEFINITIONS := $(WORKING_DIR)/source.list

EXISTING_TMPS :=

.PHONY : main help

main :

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

ifndef TASK
  TASK =#
endif
ifeq ($(TASK),)
#-----------------------------------------------------------------------------

  .NOTPARALLEL:

  GET_DEFINITIONS_SCRIPT := make/get-dependencies.pl
  ifeq ($(wildcard $(GET_DEFINITIONS_SCRIPT)),)
    GET_DEFINITIONS_SCRIPT := $(GLOBAL_ROOT)/$(GET_DEFINITIONS_SCRIPT)
  endif

#   # subdirectories that may need making.
#   SUBDIRS := 
#   # only make subdirectories if they are there
#   SUBDIRS := $(filter $(wildcard *),$(SUBDIRS))

  .PHONY : definitions clean allclean \
	all allsp allmp \
	debug debugmp debug_mp opt optmp opt_mp \
	debug64 debug64mp opt64 opt64mp

  main : definitions $(WORKING_DIR) $(EXE_DIR)

#   ifneq ($(SUBDIRS),)
# 	@for subdir in $(SUBDIRS); do \
# 	  $(MAKE) --directory=$$subdir DEBUG=$(DEBUG) ABI=$(ABI) || exit 1; \
# 	done
#   endif
	$(MAKE) --no-print-directory TASK=source
	$(MAKE) --no-print-directory TASK=objects
  ifneq ($(DYNAMIC),false)# dynamic linking
	$(MAKE) --no-print-directory TASK=executable
  endif

definitions : $(WORKING_DIR)
	@$(PERL) $(GET_DEFINITIONS_SCRIPT) \
	  $(SOURCE_DIR) '$(SOURCE_PATTERNS)' \
	  $(WORKING_DIR) $(GLOBAL_ROOT) $(DEFINITIONS)

#   ifneq ($(DYNAMIC),false)# dynamic linking
#     # if someone has removed the executable remove the record
#     ifeq ($(wildcard $(EXECUTABLE)),)
#       ifneq ($(wildcard $(EXE_RECORD)),)
# 	rm $(EXE_RECORD)
#       endif
#     endif
#   endif

  # additional dependencies for dynamic
  ifneq ($(DYNAMIC),false)# dynamic linking
    ifneq ($(EXE_DIR),$(DSO_DIR))
      main : $(DSO_DIR)

      $(DSO_DIR):
	mkdir -p $@
    endif
  endif

  $(WORKING_DIR) :
# 	@if [ -d $(OLD_WORKING_DIR) ]; then \
# 	  echo 'rm -rf $(OLD_WORKING_DIR)' && rm -rf $(OLD_WORKING_DIR); \
# 	fi
	mkdir -p $@

  $(EXE_DIR) :
# 	@if [ -d $(OLD_EXE_DIR) ]; then \
# 	  echo mv $(OLD_EXE_DIR) $@ && mv $(OLD_EXE_DIR) $@; \
# 	else \
# 	  echo mkdir -p $@ && mkdir -p $@; \
# 	fi
	mkdir -p $@; \

clean:
	@echo "Cleaning house ..."
	rm -rf $(WORKING_DIR) $(EXECUTABLE) $(DSO_DIR)

allclean:
	@echo "Cleaning house ..."
	rm -rf generated/* bin/* dso/*

  # DEBUG,ABI,MP could be available from export but being verbose
  # so user can see what is happening
debug debugmp opt optmp debug64 debug64mp opt64 opt64mp:
	$(MAKE) --no-print-directory DEBUG=$(DEBUG) ABI=$(ABI) MP=$(MP)

  debug debugmp debug64 debug64mp: DEBUG=true
  opt optmp opt64 opt64mp: DEBUG=false
  debug opt debug64 opt64: MP=false
  debugmp optmp debug64mp opt64mp: MP=true
  ifneq (,$(filter $(MACHNAME),ia64 x86_64))# ia64 or x86_64
    debug debugmp opt optmp: ABI=64
  else
  ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
    debug debugmp opt optmp: ABI=n32
  else
    debug debugmp opt optmp: ABI=32
  endif
  endif
  debug64 debug64mp opt64 opt64mp: ABI=64

  # these targets are provided for backward compatibility
  debug_mp: debugmp
  opt_mp: optmp

allsp: debug opt
allmp: debugmp optmp
all: allsp
  ifneq ($(SYSNAME),Linux)
    allsp: debug64 opt64
    allmp: debug64mp opt64mp
    all: allmp
  else
    ifeq ($(MACHNAME),ia64)
      allsp: allmp
    endif
  endif

#-----------------------------------------------------------------------------
endif

ifneq ($(TASK),)

   include $(DEFINITIONS)

   LOCAL_FILES := $(notdir $(wildcard $(SOURCE_DIR)/*.*) )
   LOCAL_CMN := $(filter %.cmn, $(LOCAL_FILES))

# NB: if fe%.f is local and fz%.f is wanted then use local fz%.f
   LOCAL_UNITS := $(basename $(SOURCES_LOCAL))
   GLOBAL_UNITS := $(basename $(SOURCES_GLOBAL))

   LOCAL_DEPEND_FILES := $(patsubst %,$(WORKING_DIR)/%.d, $(LOCAL_UNITS) )

# Object lists: lists of global objects have different names to ensure
# a different list is used if a local file is checked out or removed.

   REV_TIME_BASE := $(basename $(REV_TIME_SOURCE))

endif

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

ifeq ($(TASK),source)
#-----------------------------------------------------------------------------

PREPROCESS_FORTRAN_SCRIPT := make/preprocess-fortran.pl
STUB_FORTRAN_SCRIPT = make/fe2fz.pl

ifeq ($(wildcard $(PREPROCESS_FORTRAN_SCRIPT)),)
  PREPROCESS_FORTRAN_SCRIPT := $(GLOBAL_ROOT)/$(PREPROCESS_FORTRAN_SCRIPT)
endif
PREPROCESS_FORTRAN := ABI=$(ABI) DEBUG=$(DEBUG) FC=$(FC) $(PERL) $(PREPROCESS_FORTRAN_SCRIPT)
ifneq ($(wildcard $(STUB_FORTRAN_SCRIPT)),)
  STUB_FORTRAN := $(PERL) $(STUB_FORTRAN_SCRIPT)
else
  STUB_FORTRAN := $(PERL) $(GLOBAL_ROOT)/$(STUB_FORTRAN_SCRIPT)
endif

  main : $(LOCAL_DEPEND_FILES) \
    $(addprefix $(WORKING_DIR)/,$(LOCAL_CMN))

  # Because of the recursive nature of C included files the .d files depend on
  # each included file so we must include the .d files (if they exist) for the
  # list of their dependencies.
  LOCAL_C_DEPEND_FILES := \
	$(patsubst %.c,$(WORKING_DIR)/%.d, $(filter %.c, $(SOURCES)))
  ifneq ($(LOCAL_C_DEPEND_FILES),)
    sinclude $(LOCAL_C_DEPEND_FILES)
  endif

  # Local object lists do not contain full paths, because even the global
  # version may be install in a different location.

  # Implicit rules for generating the preprocessed source files and
  # making the object dependencies.
  # To remake the preprocessed files remove either the .d file.
  $(WORKING_DIR)/%.d : $(SOURCE_DIR)/%.f
	$(PREPROCESS_FORTRAN) $< $(WORKING_DIR)/$(dir $(patsubst $(SOURCE_DIR)/%,%,$<))

  # Stubs for unused modules when there is no fz file in source dir.
  # To remake the fz stubs remove the .d file.
  $(WORKING_DIR)/fz%.d : $(SOURCE_DIR)/fe%.f
	$(STUB_FORTRAN) $< > $(WORKING_DIR)/fz$*.f
	echo '$(WORKING_DIR)/fz$*.o: $(WORKING_DIR)/fz$*.f' > $@

  $(WORKING_DIR)/fz%.d : $(GLOBAL_ROOT)/$(SOURCE_DIR)/fe%.f
	$(STUB_FORTRAN) $< > $(WORKING_DIR)/fz$*.f
	echo '$(WORKING_DIR)/fz$*.o: $(WORKING_DIR)/fz$*.f' > $@

  # wildcard in the .d file removes files that no longer exist
  # (that may no longer be needed because they were included in .h files that
  #  have changed).
  # For each global prerequisite include the local version as a possible
  # prerequisite as this may be checked out into the local directory in
  # the future.
  $(WORKING_DIR)/%.d : $(SOURCE_DIR)/%.c
	( makedepend $(CPPFLAGS) -f- $< | \
	  sed -e 's%\($(GLOBAL_ROOT)/\)\([^ ]*.h\)%\1\2 \2%g' \
	    -e 's%^$(SOURCE_DIR)\([^ ]*\).o: \(.*\)%$(WORKING_DIR)\1.o $(WORKING_DIR)\1.d: $$(wildcard \2)%' && \
	  echo '$(WORKING_DIR)/$*.o: $<' ) > $@

  # implicit rule for preprocessing common blocks.
  $(WORKING_DIR)/%.cmn : $(SOURCE_DIR)/%.cmn
	$(PREPROCESS_FORTRAN) $< $(WORKING_DIR)

#-----------------------------------------------------------------------------
endif

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

ifeq ($(TASK),objects)
#-----------------------------------------------------------------------------

  OBJECTS := \
	$(patsubst %,$(WORKING_DIR)/%.o, $(LOCAL_UNITS) ) \
	$(patsubst %,$(GLOBAL_ROOT)/$(WORKING_DIR)/%.o, $(GLOBAL_UNITS) )

  OBJECTS_LIST := $(WORKING_DIR)/object.list

  vpath %.cmn $(CMN_INCLUDE_DIRS)
  %.cmn:
	$(warning $@ not found in any of $(CMN_INCLUDE_DIRS))
	@exit 1
  # This line is for physiomers who include .f files.
#  vpath %.f $(INC_INCLUDE_DIRS)
  vpath %.inc $(INC_INCLUDE_DIRS)
  %.inc:
	$(warning $@ not found in any of $(INC_INCLUDE_DIRS))
	@exit 1

#
# New %.h rule for BOXMG include files.
#
  vpath %.h $(INC_INCLUDE_DIRS)
  %.h:
	$(warning $@ not found in any of $(INC_INCLUDE_DIRS))
	@exit 1

  vpath %.a $(ARCHIVE_DIRS)
  %.a:
	$(warning $@ not found.)
	@exit 1

#   # use global objects for non-existing local objects
#   vpath %.o $(GLOBAL_ROOT)

#-----------------------------------------------------------------------------
# include the object (and .mod) dependencies)

  ifneq ($(DYNAMIC),false)# dynamic linking
    TARGETS := $(DSOS)
    TMPS := $(foreach target,$(TARGETS),$(call get_build_filename,$(target)))
  else
    TARGETS := $(EXECUTABLE)
    TMP_EXE := $(call get_build_filename,$(EXECUTABLE))
    TMPS := $(TMP_EXE)
  endif

  # Ask make to remove temporary files if the build fails or if make is
  # interrupted.
  .INTERMEDIATE : $(TMPS)

  EXISTING_TMPS := $(wildcard $(TMPS))

  main : $(TARGETS)

ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
  #SAB  We want to hide all the perl symbols and in the static version we can
  #do this by hiding everything by exporting nothing.
  EXPORT_SYMBOLS_FLAGS = -exported_symbol ""
else
  EXPORT_SYMBOLS_FLAGS =
endif

  ifeq ($(DYNAMIC),false)# not dynamic linking
    # explicit rule for the executable
    # the static executable must be remade if objects change.

    # Unlinking the old executable and replacing (rather than overwriting)
    # means that any process based on the old executable can continue as it
    # will not be overwritten.
    $(EXECUTABLE) : $(TMP_EXE)
	mv -f $< $@

    $(TMP_EXE) : $(OBJECTS) $(CELLML_ARCHIVES) $(ARCHIVES) $(EXPORTS) $(OBJECTS_LIST)
	$(EXE_LINK) $(ELFLAGS) $(REQUIRED_BEGIN_OBJECTS) \
	  $(EXPORT_SYMBOLS_FLAGS) $(call olist_args,$(OBJECTS_LIST)) \
	  $(filter %.a, $^) \
	  $(addprefix -l, $(DYNAMIC_LIBRARIES)) $(REQUIRED_END_OBJECTS) \
	  -o $@
  endif

  ifneq ($(LOCAL_DEPEND_FILES),)
    include $(LOCAL_DEPEND_FILES)
  endif

#-----------------------------------------------------------------------------

  $(DSO_DIR)/lib%.so : $(call get_build_filename,$(DSO_DIR)/lib%.so)
	mv -f $< $@

  $(call get_build_filename,$(DSO_DIR)/lib%.so) : $(WORKING_DIR)/%.olist
	$(call dso_build,$@,$(call olist_args,$<))

  TMP_FECELLML = $(call get_build_filename,$(DSO_DIR)/libfecellml.so)
  $(TMP_FECELLML) : $(CELLML_ARCHIVES)
  ifeq ($(SYSNAME),AIX)
    $(TMP_FECELLML) : DSO_LIBS += C
# ld -o dso/rs6000-64-aix-debug/libfecellml.so.build -b64 -G -bexpall -bnoentry -f generated/rs6000-64-aix-debug/fecellml.olist /hpc/cmiss/cellml/lib/rs6000-64-aix/libCMISSCellML-debug.a /hpc/cmiss/cellml/lib/rs6000-64-aix/libcellml-debug.a /hpc/cmiss/cellml/lib/rs6000-64-aix/libxerces-c-debug.a -lxlf90 -lc -lC
# ld: 0711-781 ERROR: TOC overflow. TOC size: 67208       Maximum size: 65536
    $(TMP_FECELLML) : D_FLGS += -bbigtoc
  endif

# convert_sparse_5to1 referenced from fe07
#   $(DSO_DIR)/libfesolver.so : $(SOLVER_ARCHIVES)

# GX also referenced from fe14 and fe01
#   $(DSO_DIR)/libfegx.so : $(GX_ARCHIVES)

# Can't do this currently as the static perl interpreter is not built with
# -fPIC and if it is the symbols are exported and conflict with those of a
# dlopened perl library.
#   $(DSO_DIR)/libfeinterpreter.so : $(PERL_INTERPRETER_ARCHIVE)

  # implicit rules for making objects from C or fortran source
  $(WORKING_DIR)/%.o : $(WORKING_DIR)/%.f
	$(FC) -o $@ $(FFLAGS) $(FPPFLAGS) $<
  $(WORKING_DIR)/%.o : $(SOURCE_DIR)/%.c
	$(CC) -o $@ $(CPPFLAGS) $(CFLAGS) $<

  # implicit rules for making objects from C or fortran source
  $(WORKING_DIR)/%.o : $(GLOBAL_ROOT)/$(SOURCE_DIR)/%.f
	$(FC) -o $@ $(FFLAGS) $(FPPFLAGS) $<
  $(WORKING_DIR)/%.o : $(GLOBAL_ROOT)/$(SOURCE_DIR)/%.c
	$(CC) -o $@ $(CPPFLAGS) $(CFLAGS) $<

  # special rule for revision time object which
  # depends on all other objects and archives
  ifneq ($(filter $(REV_TIME_BASE),$(LOCAL_UNITS)),)
    $(WORKING_DIR)/$(REV_TIME_BASE).o : \
      $(SOURCE_DIR)/$(REV_TIME_SOURCE) $(ARCHIVES) \
      $(filter-out $(WORKING_DIR)/$(REV_TIME_BASE).o, $(OBJECTS))
	$(CC) -o $@ -D REVISION_TIME_STRING=\""`date`"\" $(CFLAGS) \
	  $(SOURCE_DIR)/$(REV_TIME_SOURCE)
  endif

  # objects that should not be produced with optimization
# Hopefully these are OK now that we have -LNO:opt=0
#   ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
#     # optimization causes runtime bug
#     $(WORKING_DIR)/fe13/IPMESH3.o \
#     $(WORKING_DIR)/fe26/IMSIGN.o \
# 	: OPTCF_FLGS =
#   endif
  ifeq ($(SYSNAME),SunOS)
    # these take too long to compile
    $(WORKING_DIR)/fe00/IOGEOM.o \
    $(WORKING_DIR)/fe11/IPPARA.o \
    $(WORKING_DIR)/fe13/IPMESH2.o \
	: OPTCF_FLGS =
    # optimization causes runtime seg fault (something with ERROR)
    $(WORKING_DIR)/fe20/FEM_DYNAM.o \
	: OPTCF_FLGS =
    # optimization is turned on if mp is specified
    $(WORKING_DIR)/fe00/IOGEOM.o \
    $(WORKING_DIR)/fe11/IPPARA.o \
    $(WORKING_DIR)/fe13/IPMESH2.o \
    $(WORKING_DIR)/fe20/FEM_DYNAM.o \
	: MP_FLGS =
  endif
  ifeq ($(filter-out xlf%,$(FC)),)# xlf* fortran compiler
#     # Work around run-time seg fault with xlf 9.1 (example 118 et al)
#     $(WORKING_DIR)/fe02/GLOBALH.o : OPTF_FLGS += -qnohot
    # Work around run-time seg fault with xlf 9.1
    $(WORKING_DIR)/fe02/CALC_NY_GRID_DEP.o : OPTF_FLGS += -O2
  endif
  ifeq ($(SYSNAME),AIX)
    # compiler requests this
#     $(WORKING_DIR)/fe07/MARCH8.o : F_FLGS += -qspillsize=2048
    $(WORKING_DIR)/fe08/lsoda_wrapper.o : F_FLGS += -qspillsize=1024
    $(WORKING_DIR)/fe20/FEM_DYNAM.o : F_FLGS += -qspillsize=10000
    $(WORKING_DIR)/fe21/SOLVE.o : F_FLGS += -qspillsize=2048
    $(WORKING_DIR)/fe21/FIT.o : F_FLGS += -qspillsize=1024
    $(WORKING_DIR)/fe24/DEINIT.o : F_FLGS += -qspillsize=1024
#     $(WORKING_DIR)/fe24/DEMESH.o : F_FLGS += -qspillsize=1024
#     $(WORKING_DIR)/fe29/EVOBJE.o : F_FLGS += -qspillsize=1024
    $(WORKING_DIR)/fe29/EVRESI.o : F_FLGS += -qspillsize=1024	
#     $(WORKING_DIR)/fe07/GENSOL.o : F_FLGS += -qspillsize=1024
#     # seg fault without this
#     $(WORKING_DIR)/fe00/GINOUT.o : OPTF_FLGS += -qstrict
#     $(WORKING_DIR)/fe00/IPMESH3.o : OPTF_FLGS += -qstrict

# Seems to be fixed in xlf 09.01.0000.0006
#     # With 09.01.0000.0005 and bounds checking, local variables (OPFILE) and
#     # dummy arguments (ERROR) are corrupted in SOLVE on return from GENSOL in
#     # cm64-debug.
# Seems to be fixed (or changed) with xlf 09.01.0000.0005
#     # XL Fortran 9.1 bounds checking seems to result in some pointers to
#     # arguments to other subprograms (XG to MARCH6) being set to strange
#     # values (7ff7ffff or nil).  Changing the unspecified dimension of
#     # CSEG(*)*(*) and RESJAC(NREM,*) to 1 in gensol seemed to work around the
#     # problem:
#     # And bounds checking makes ipa hang:
#     $(WORKING_DIR)/fe07/GENSOL.o : DBGCF_FLGS += -qnocheck
    ifneq ($(MP),false)
      # ATOMIC directives and -qcheck/-C don't seem to mix well with
      # xlf95_r 10.01.0000.0002.
      $(WORKING_DIR)/fe07/ZPRP_DYNAM.o : DBGCF_FLGS += -qnocheck
      $(WORKING_DIR)/fe07/ASSEMBLE5_DYNAM.o : DBGCF_FLGS += -qnocheck
      $(WORKING_DIR)/fe07/SOLVE5.o : DBGCF_FLGS += -qnocheck
# Seems to be fixed with xlf 09.01.0000.0005
#       # bounds checking makes ipa hang:
#       $(WORKING_DIR)/fe29/EVRESI.o : DBGCF_FLGS += -qnocheck
      # -qinitauto and -qsmp=omp:noopt seem incompatible for this object
      # with no auto variables
      # 1586-494 (U) INTERNAL COMPILER ERROR: Signal 11.
      $(WORKING_DIR)/revision_time.o : MP_FLGS =
    endif
#     # These modules have or call modules with inconsistent types.
#     $(WORKING_DIR)/fe00/BINREADFILE.o : DBGCF_FLGS += -qnoextchk
#     $(WORKING_DIR)/fe00/BINWRITEFILE.o : DBGCF_FLGS += -qnoextchk

  endif
  ifeq ($(FC),pathf90)
    # With -O3 this is still running after 2.5 hours in CGLOOP phase.
    $(WORKING_DIR)/fe11/IPCORN.o : OPTF_FLGS = -O2
  endif
  ifeq ($(FC),pgf90)
    # pgf90 5.2 seg faults with -mp even though there are no openmp directives.
    $(WORKING_DIR)/fe29/EVELEC.o : MP_FLGS =
  endif
  ifeq ($(FC),ifort)
    ifeq ($(MACHNAME),ia64)
      # ifort 8.1 seg faults with -O3
      $(WORKING_DIR)/fe12/OPTCAP.o : OPTF_FLGS = -O2
      # SIGILL with ifort 8.1 -O3
      $(WORKING_DIR)/fe26/EXNODE.o : OPTF_FLGS = -O2
    endif
  endif

#-----------------------------------------------------------------------------
endif

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

ifeq ($(TASK),executable)
#-----------------------------------------------------------------------------

# This is done separately for dynamic linking because dsos must exist
# but the executable does not need to be recompiled when dsos are
# updated.

  # It would be nice to test that the executable works if any DSOs have
  # changed as they may require new symbols that may not be linked into
  # the executable from the archives.
  # This could also replace the cautious dependency of the executable on
  # the Makefile changes, which is there in case the list of DSOs changes.

  TMP_EXE := $(call get_build_filename,$(EXECUTABLE))
  EXISTING_TMPS := $(wildcard $(TMP_EXE))

  main : $(EXECUTABLE)

  vpath %.a $(ARCHIVE_DIRS)
  %.a:
	$(error $@ not found.)

#   # touch the executable if dsos have changed to prevent SGI's cvd's
#   # `source newer than executable message'.

#   $(EXECUTABLE) : $(DSOS) $(EXE_RECORD)
# ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
# 	touch $@
# endif

  $(EXECUTABLE) : $(TMP_EXE)
	mv -f $< $@

  # Ask make to remove this temporary file if the build fails or if make is
  # interrupted.
  .INTERMEDIATE : $(TMP_EXE)

ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
  #Exclude all these symbols from the SGI DSO version
  EXPORT_SYMBOLS_FILE = $(WORKING_DIR)/perl.symbols
  EXPORT_SYMBOLS_FLAGS = -hiddens_file $(EXPORT_SYMBOLS_FILE)
else
  EXPORT_SYMBOLS_FLAGS =
endif

  # explicit rule for the executable
  # Initally build as $(*).link and then move to $(*), since AIX's
  #   ld will fail if the old executable is in use.
  $(TMP_EXE) : $(ARCHIVES) Makefile
ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
   #SAB We no longer want to export the Perl symbols but for the DSO cm
   #to work we need to export everything else so we need to specify an 
   #exclusion list.
	elfdump -Dt $(PERL_INTERPRETER_ARCHIVE) | grep '^\[' | cut -f3 | grep -vE '^interpret|^[a-z].*_interpreter_$$' > $(EXPORT_SYMBOLS_FILE);
endif
	$(EXE_LINK) $(ELFLAGS) -L$(DSO_DIR) \
	  $(EXPORT_SYMBOLS_FLAGS) $(REQUIRED_BEGIN_OBJECTS) $(addprefix -l, $(DSO_UNITS)) \
	  $(filter %.a, $^) \
	  $(addprefix -l, $(DYNAMIC_LIBRARIES)) \
	  $(REQUIRED_END_OBJECTS) -o $(@)

#-----------------------------------------------------------------------------
endif

# Remove any existing temporary files incase someone else created a file of
# this name and gave its parent directory to us.  This also reduces the
# chance of two processes owned by the same user working with the same file.
.PHONY : existing_dummy clear_existing_tmps

# existing_dummy is a place filler for the case when EXISTING_TMPS is empty
$(EXISTING_TMPS) existing_dummy : clear_existing_tmps

clear_existing_tmps :
	rm $(EXISTING_TMPS)


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

help:
	@echo "			Compile a version of cm"
	@echo "			======================="
	@echo
	@echo "Examples of usage:   "
	@echo
	@echo "	gmake"
	@echo "	gmake OPT= ABI=n32 MP="
	@echo "	gmake debug64mp"
	@echo " gmake CELLML_ROOT=/path/to/cool/stuff/cellml"
	@echo " gmake USE_CELLML=no_thanks"
	@echo
	@echo "Options: (The former is the default unless specified.)"
	@echo
	@echo "	(DEBUG=|OPT=)"
	@echo "	ABI=(32|n32|64)"
	@echo "	MP=(false|)"
	@echo "		These specify the version to compile."
	@echo "		ABI=n32 is only available on IRIX."
	@echo "		On IRIX, ABI=32 is not supported and"
	@echo "		ABI=n32 is the default."
	@echo "		On IA64 machines ABI=32 is not supported and"
	@echo "		ABi=64 is the default."
	@echo
	@echo "	DYNAMIC=(|false)"
	@echo "		A dynamic executable requires much "
	@echo "		less time to link, but relies on the existence of"
	@echo "		shared objects for execution.  Changing these shared"
	@echo "		objects can change the behaviour of the executable."
	@echo
	@echo "		When DYNAMIC=false is specified, the behaviour of the"
	@echo "		executable will not change unless the executable"
	@echo "		is overwritten.  A copy of the executable may be "
	@echo "		made, so that execution of the copy is unaffected by "
	@echo "		future compilations."
	@echo
	@echo "         DYNAMIC=false is the default with OPT=."
	@echo
	@echo " USE_CELLML=(true|)"
	@echo "         Set this variable to anything other than true and "
	@echo "         you will not get a CellML enabled executable."
	@echo
	@echo " CELLML_ROOT=(CMISS_ROOT/cellml|path to CellML root)"
	@echo "         Defaults to the cellml directory in CMISS_ROOT but "
	@echo "         can be overridden with a different path if required."
	@echo "         The CellML libraries are found in CELLML_ROOT/lib."
	@echo
	@echo " CELLML_DEBUG=(|true)"
	@echo "         By default cm is linked with the optimised versions "
	@echo "         of the CellML libraries, but the debug libraries "
	@echo "         can be used by setting this variable to true."
	@echo
	@echo "Available targets:                            "
	@echo
	@echo "	clean"
	@echo "		Remove generated files associated with a single"
	@echo "		version."
	@echo
	@echo "	allclean"
	@echo "		Remove all generated files."
	@echo
	@echo "	help"
	@echo "		Display this message."
	@echo
	@echo "	debug opt debug64 opt64 debugmp optmp debug64mp opt64mp"
	@echo "		Compile the specified version with automatic setting"
	@echo "		of DEBUG, ABI, and MP."
	@echo
	@echo "	all"
	@echo "		Compile all versions."
	@echo
	@echo "	allmp"
	@echo "		Compile all multiprocessor versions."
	@echo
	@echo "	allsp"
	@echo "		Compile all single processor versions."
	@echo

