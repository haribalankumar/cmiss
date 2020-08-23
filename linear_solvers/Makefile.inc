# -*- Makefile -*-
# 
# For use with GNU make.

#-----------------------------------------------------------------------------
# Makefile include file for building the CMISS solver/linear algebra libraries.
# The solver, xblas, SuperLU and Umfpack libraries Makefile all include this 
# file.
#
# Original by Stuart Norris 2002-10-07
# Changes:
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Set the Operating System name, Machine name, and flags as to whether the
# code is compiled with Optimisation, Debug, or MultiProcessor flags.
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
    DEBUG = true
  else
    DEBUG = false
  endif
endif

ifndef MP
  MP = false
endif

#-----------------------------------------------------------------------------
# Set suffixes used in naming libraries, and include directories for the
# header files from the SuperLU and Umfpack libraries.
#-----------------------------------------------------------------------------
ifeq ($(MP),false)
  MT_SUFF =
  MP_SUFF =
else
  MT_SUFF = -mt
  MP_SUFF = _mp
endif
ifeq ($(OPT),false)
  DEBUG_SUFF = -debug
else
  # DEBUG_SUFF = -opt
  DEBUG_SUFF =
endif

#-----------------------------------------------------------------------------
# Set options for each platform. Currently AIX, IRIX, OSF, Linux and Win32
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# AIX Options
#-----------------------------------------------------------------------------
ifeq ($(SYSNAME),AIX)

  ifndef ABI
    ABI = 32
  endif

  ifeq ($(ABI),32)
    POINTER_LENGTH = 4
  else
    POINTER_LENGTH = 8
  endif

  #
  # T.Austin: 26/4/07
  #
  #    Changed XLF compiler to be xlf95 in line with cm and 
  #    also because xlf_r was giving problems to boxmg that 
  #    were fixed by xlf95_r.  Hopefully long term this does
  #    not present any major problems upon recompilation of 
  #    linear_solvers.
  #

  ifeq ($(MP),false)
    XLF      = xlf
    SCILIB   = -lessl -llapack$(SUFFIX) -lm
    MP_FLAGS =
  else
    XLF      = xlf_r
    SCILIB   = -lesslsmp -llapack$(DEBUG_SUFF)

    ifeq ($(BUILDING_BLAS)$(BUILDING_LAPACK),true)
      # No OpenMP in blas and lapack, and lapack contains lines that look
      # like conditional compilation lines.
      # Turn off the omp to turn off conditional compilation
      # and try automatic parallelization.
      ifeq ($(DEBUG),false)
        MP_FLAGS = -qsmp
      else
        MP_FLAGS = -qsmp=noopt
      endif
    else
      ifeq ($(DEBUG),false)
        MP_FLAGS = -qsmp=omp:opt
      else
        MP_FLAGS = -qsmp=omp:noopt
      endif
    endif
  endif

  # 1512-047 (E) The I FORMAT edit descriptor, in literal FMT specifier, is 
  # missing a field width. A default value based on the edit descriptor and 
  # the type of the data item being processed is assumed.
  CFLAGS_ABI   = -q$(ABI) -qarch=auto
  FFLAGS_ABI   = -q$(ABI) -qextname -qfixed -qsuppress=1512-047 -qarch=auto

  # 1500-036 (I) The NOSTRICT option (default at OPT(3)) has the potential
  # to alter the semantics of a program.  Please refer to documentation on
  # the STRICT/NOSTRICT option for more information.
  FFLAGS_OPT   = $(FFLAGS_ABI) -O3 -qstrict -qsuppress=1500-036 -qmaxmem=12000 -qtune=auto
  CFLAGS_OPT   = $(CFLAGS_ABI) -O3 -qstrict -qsuppress=1500-036 -qtune=auto

  # 7FF7FFFF is a large integer or NaNQ real*4 or NaNS real*8
  FFLAGS_DEBUG = $(FFLAGS_ABI) -g -qfullpath -C -qflttrap=ov:zero:inv:en -qinitauto=7FF7FFFF
  CFLAGS_DEBUG = $(CFLAGS_ABI) -g -qfullpath -qinitauto=7F

  ifeq ($(DEBUG),false)
    FFLAGS = $(FFLAGS_OPT) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_OPT) $(MP_FLAGS) $(CINCDIR)
  else
    FFLAGS = $(FFLAGS_DEBUG) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_DEBUG) $(MP_FLAGS) $(CINCDIR)
  endif

  FC       = $(XLF)
  FPPFLAGS = 
  CC       = xlc
  CPPFLAGS = 
  LD       = $(FC)
  LDFLAGS  = $(FFLAGS)
  AR       = ar
  AFLAGS   = -rv -X $(ABI)
  RANLIB   = ranlib

  DIRNAM   = rs6000-$(ABI)-aix

  # Library dependant rules
  # BOXMG library
  ifeq ($(BUILDING_BOXMG),true)
    FPPFLAGS += -WF,-D'POINTER=INTEGER*$(POINTER_LENGTH)'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -D_OMP
  endif

  # Solver library
  ifeq ($(BUILDING_SOLVER),true)
    FPPFLAGS += -WF,-D'POINTER=INTEGER*$(POINTER_LENGTH)'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -D_OMP -DBSD_TIMERS
  endif

  # Config for Umfpack 4.0
  ifeq ($(BUILDING_UMFPACK),true)
    # 1: with the vendor-supplied BLAS.  This is the default.
    CONFIG = -DNCBLAS -bmaxdata:2000000000
  endif

  # Config for SuperLU
  ifeq ($(BUILDING_SUPERLU),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = $(SCILIB)
    MATHLIB = -lm
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
  endif

  # Config for SuperLU_omp
  ifeq ($(BUILDING_SUPERLU_MP),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = $(SCILIB)
    MATHLIB = -lm
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
    # Type of parallelisation: _SGI, _OMP, and type of platform
    CFLAGS += -D_OMP
    ifeq ($(DEBUG),false)
      BLASOPTS = -qalias=allp
    else
      BLASOPTS = -g
    endif
  endif

endif

#-----------------------------------------------------------------------------
# SGI Options
#-----------------------------------------------------------------------------
ifeq ($(SYSNAME:IRIX%=),)

  ifndef ABI
    ABI = n32
  endif

  ifeq ($(ABI),n32)
    POINTER_LENGTH = 4
  else
    POINTER_LENGTH = 8
  endif

  ifndef MIPS
    # Using mips3 for n32 versions on esu* machines so they will run on Indies.
    MIPS = 4
    ifeq ($(NODENAME:esu%=),)
      ifeq ($(ABI),n32)
        ifneq ($(DEBUG),false)
          MIPS = 3
        endif
      endif
    endif
  endif

  ifeq ($(MP),false)
    MP_FLAGS =
  else
    MP_FLAGS = -mp
  endif

  CFLAGS_ABI   = -mips$(MIPS) -$(ABI)
  CFLAGS_OPT   = $(CFLAGS_ABI) -O3 -OPT:alias=RESTRICT
  CFLAGS_DEBUG = $(CFLAGS_ABI) -O0 -g
  FFLAGS_ABI   = -mips$(MIPS) -$(ABI)
  FFLAGS_OPT   = $(FFLAGS_ABI) -O3 -OPT:alias=parm
  FFLAGS_DEBUG = $(FFLAGS_ABI) -O0 -g

  ifeq ($(DEBUG),false)
    FFLAGS = $(FFLAGS_OPT) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_OPT) $(MP_FLAGS) $(CINCDIR)
  else
    FFLAGS = $(FFLAGS_DEBUG) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_DEBUG) $(MP_FLAGS) $(CINCDIR)
  endif

  FC       = f77
  FPPFLAGS = 
  CC       = cc
  CPPFLAGS = 
  LD       = $(FC)
  LDFLAGS  = $(FFLAGS)
  AR       = ar
  AFLAGS   = -rvf
  RANLIB   = true

  SCILIB   = -lblas$(MP_SUFF) -lscs$(MP_SUFF) -lm
  DIRNAM   = mips-$(ABI)-irix

  # Library dependant rules
  # BOXMG library
  ifeq ($(BUILDING_BOXMG),true)
    FPPFLAGS += -D'POINTER=INTEGER*$(POINTER_LENGTH)'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -D_SGI 
  endif

  # Solver library
  ifeq ($(BUILDING_SOLVER),true)
    FPPFLAGS = -D'POINTER=INTEGER*$(POINTER_LENGTH)'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -D_SGI -DBSD_TIMERS
  endif

  # Config for Umfpack 4.0
  ifeq ($(BUILDING_UMFPACK),true)
    # 1: with the vendor-supplied BLAS.  This is the default.
    CONFIG = -DNCBLAS
    # 2: with no BLAS
    # CONFIG = -DNBLAS
    # 3: with no BLAS, in 64-bit mode
    # CONFIG = -DNBLAS -DLP64 -64
    # 4: with 32-bit int's, 64-bit long's, and the 64-byte SCSL BLAS
    # Compare the output with umfpack_**_demo.out64
    # CONFIG = -DNCBLAS -DLP64 -64
    # 5: with the Fortran interface to the vendor-supplied BLAS
    # CONFIG = -DNCBLAS -DNSCSL
    # 6: with the C-BLAS interface to the vendor-supplied BLAS
    # NOTE: with gcc and the 07/09/98 version of SGI's cblas.h, spurious warnings
    # about complex BLAS arguments will be reported.  This is because this older
    # version of SGI's cblas.h does not make use of the SCSL_VOID_ARGS parameter,
    # which is present in the 12/6/01 version of SGI's cblas.h.  You can safely
    # ignore these warnings.  See umf_config.h for more information.
    # CONFIG = -DNSCSL
  endif

  # Config for SuperLU
  ifeq ($(BUILDING_SUPERLU),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = -lblas$(MP_SUFF)
    MATHLIB = -lm
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
  endif

  # Config for SuperLU_omp
  ifeq ($(BUILDING_SUPERLU_MP),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = -lblas$(MP_SUFF)
    MATHLIB = -lm
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
    # Type of parallelisation: _SGI, _OMP, and type of platform
    CFLAGS += -D_SGI
  endif

endif

#-----------------------------------------------------------------------------
# Linux Options: use g77, or the native Digital and Intel compilers.
#-----------------------------------------------------------------------------
ifeq ($(SYSNAME),Linux)

  ifndef INSTRUCTION
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      INSTRUCTION := i686
    else
      INSTRUCTION = $(MACHNAME)
    endif
  endif

  CC = gcc
  FC = g77
#  ifneq ($(filter $(INSTRUCTION),ia64 x86_64 i686),)# ia64, x86_64, or i%86
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
    # pathf90 --version is successful even if the subscription is not valid.
    # Try preprocessing a small file.
    ifneq (,$(shell pathf90 -E /dev/null 2>/dev/null && echo 1))
      FC = pathf90
    else
    # gcc and Portland fortran as used by AMD
    ifneq (,$(shell pgf90 -V 2>/dev/null))
      FC = pgf90
      ifneq ($(MP),false)
        # need compiler with mp
        CC = pgcc
      endif
    endif
    endif
  endif
  ifeq ($(MACHNAME),ppc64)
    # use xlf_r if available (_r so libraries are thread safe)
    ifneq (,$(shell xlf_r 2>/dev/null))
      FC = xlf_r
    endif
    ifneq (,$(shell xlc_r 2>/dev/null))
      CC = xlc_r
    endif
  endif

  ifndef ABI
    ifeq (,$(filter $(MACHNAME),ia64 x86_64))# not ia64 or x86_64
      ABI = 32
    else
      ABI = 64
    endif
  endif
  ifeq ($(ABI),64)
    POINTER_LENGTH = 8
  else
    POINTER_LENGTH = 4
  endif

  #CFLAGS_ABI  = -ansi -pedantic -W -Wall -Wno-parentheses -Wshadow\
  #              -Wcast-align -Winline -Wstrict-prototypes -fPIC

  ifeq ($(CC),gcc)
    CFLAGS_ABI   = -pipe -ansi -pedantic -W -Wall# -fPIC# -march=pentium4# -mfpmath=sse
    ifeq ($(filter $(INSTRUCTION),i686 ia64),)# not i686 nor ia64
      CFLAGS_ABI += -m$(ABI)
    endif
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        CFLAGS_ABI += -march=nocona
      endif
    endif
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      CFLAGS_ABI += -malign-double
#       ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
#         CFLAGS_ABI += -march=pentium4# -mfpmath=sse
#       endif
    endif
    ifeq ($(MACHNAME),ppc64)
      CFLAGS_ABI += -mcpu=power5
    endif
    CFLAGS_OPT = $(CFLAGS_ABI) -O3# -funroll-all-loops
    ifeq ($(PROF),false)
      ifneq ($(filter $(INSTRUCTION),i686 x86_64),)# i686 or x86_64
	CFLAGS_OPT += -momit-leaf-frame-pointer
      endif
    else
      CFLAGS_ABI += -pg -g
    endif
  endif
  ifeq ($(CC),pathcc)
    CFLAGS_ABI   = -pipe -ansi -pedantic -W -Wall
    CFLAGS_OPT = $(CFLAGS_ABI) -O3 -OPT:alias=typed# -OPT:Ofast -LNO:prefetch=0
  endif
  ifeq ($(CC),icc)
    CFLAGS_ABI = -ansi -w2
    # -ipo0 
    CFLAGS_OPT   = $(CFLAGS_ABI) -O3 -ansi_alias
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        CFLAGS_ABI += -xP# for sse3 (90nm Pentium 4 series)
      else
        CFLAGS_ABI += -xW# Pentium4 compatible (?sse2)
      endif
    endif
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
        CFLAGS_ABI += -xN# for Pentium 4
      endif
    endif
    ifeq ($(MACHNAME),ia64)
      ifeq ($(BUILDING_UMFPACK),true)
        # Work around a run-time seg fault with icc 8.1 on ia64 in
        # init_rows_cols(umf_colamd.c) with example 21d.
        CFLAGS_OPT   = $(CFLAGS_ABI) -O2 -ansi_alias
      endif
    endif
  endif
  ifeq ($(CC),pgcc)
    CFLAGS_ABI   = -ansi -pedantic -W -Wall# -fPIC
    CFLAGS_OPT = -fastsse# -Mipa=fast
  endif
  CFLAGS_DEBUG = $(CFLAGS_ABI) -O0 -g

  # G77 compiler
  ifeq ($(FC),g77)

    # Switch off globals warning (with -Wno-globals) to avoid being told
    # about %VAL(P) being a problem. 
    FFLAGS_ABI   = -pipe -fno-second-underscore -W -Wall -Wsurprising -Wunused\
                   -Wno-globals# -fPIC
    ifeq ($(filter $(INSTRUCTION),i686 ia64),)# not i686 nor ia64
      FFLAGS_ABI += -m$(ABI)
    endif
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        FFLAGS_ABI += -march=nocona
      endif
    endif
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      FFLAGS_ABI += -malign-double
#       ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
#         FFLAGS_ABI += -march=pentium4# -mfpmath=sse
#       endif
    endif
    ifeq ($(MACHNAME),ppc64)
      FFLAGS_ABI += -mcpu=power5
    endif

    FFLAGS_OPT = $(FFLAGS_ABI) -O3# -funroll-all-loops
    ifeq ($(PROF),false)
      ifneq ($(filter $(INSTRUCTION),i686 x86_64),)# i686 or x86_64
	FFLAGS_OPT += -momit-leaf-frame-pointer
      endif
    else
      FFLAGS_ABI += -pg -g
    endif
    FFLAGS_DEBUG = $(FFLAGS_ABI) -O0 -g
    MP_FLAGS     =

  # DEC Compiler
  else
  ifeq ($(INSTRUCTION),alpha)

    FFLAGS_ABI   = -assume underscore -assume no2underscores
    FFLAGS_OPT   = $(FFLAGS_ABI) -O4 -automatic -fast
    FFLAGS_DEBUG = $(FFLAGS_ABI) -O0 -g

    ifeq ($(MP),false)
      MP_FLAGS = 
    else
      MP_FLAGS = -check omp_bindings -omp
    endif

    FC       = fort

  endif
  endif
  ifeq ($(FC),pathf90)
    FFLAGS_ABI   = -pipe -fno-second-underscore -W -Wall -Wsurprising -Wunused\
                   -Wno-globals# -fPIC# -march=pentium4# -mfpmath=sse
    FFLAGS_OPT = $(FFLAGS_ABI) -O3# -OPT:Ofast -LNO:prefetch=0
    FFLAGS_DEBUG = $(FFLAGS_ABI) -g# -O0 implied by -g
  endif
  ifeq ($(FC),ifort)
    # Intel compiler

    # turn on warnings,
    # suppress messages about non-standard fortran (including REAL*8,
    # more than 19 continuation lines),
    # suppress comment messages (including obsolescent alternate return),
    # auto arrays (as well as scalars), no aliasing,
    # 72 columns, temporary arrays on stack if possible.
    FFLAGS_ABI   = -W1 -w95 -cm  -auto -72 -stack_temps
#     FFLAGS_ABI   = -Vaxlib -w95
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        FFLAGS_ABI += -xP# for sse3 (90nm Pentium 4 series)
      else
        FFLAGS_ABI += -xW# Pentium4 compatible (?sse2)
      endif
    endif
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
        FFLAGS_ABI += -xN# for Pentium 4
      endif
    endif
    # I can't see any improvement with -fno-alias so
    # it doesn't seem worth the risk.
    FFLAGS_OPT   = $(FFLAGS_ABI) -O3
#     FFLAGS_OPT   = $(FFLAGS_ABI) -O3 -fno-alias -pad -mp1 -unroll
    FFLAGS_DEBUG = $(FFLAGS_ABI) -O0 -g -traceback -check format# -check bounds
    # Problem with -check bounds
    # ld: bin/ia64-linux/cm-mt-debug.link: short data segment overflowed (0x59d9c0 >= 0x400000)
    # ld: can't relax section: No such file or directory
    ifneq ($(MACHNAME),ia64)
      FFLAGS_DEBUG += -check bounds
    endif

    ifeq ($(MP),false)
      MP_FLAGS = 
    else
      MP_FLAGS = -openmp
    endif

  endif
  ifeq ($(FC),pgf90)
    FFLAGS_OPT = -fastsse# -Mipa=fast
    ifeq ($(MP),false)
      MP_FLAGS = 
    else
      MP_FLAGS = -mp
    endif
  endif

  ifeq ($(CC),xlc_r)

    CFLAGS_ABI   = -q$(ABI) -qarch=auto

    # 1500-036 (I) The NOSTRICT option (default at OPT(3)) has the potential
    # to alter the semantics of a program.  Please refer to documentation on
    # the STRICT/NOSTRICT option for more information.
    CFLAGS_OPT   = $(CFLAGS_ABI) -O3 -qcache=auto -qsuppress=1500-036 -qtune=auto

    CFLAGS_DEBUG = $(CFLAGS_ABI) -g -qfullpath

#   # Config for Umfpack 4.0
#   ifeq ($(BUILDING_UMFPACK),true)
#     # 1: with the vendor-supplied BLAS.  This is the default.
#     CONFIG = -DNCBLAS -bmaxdata:2000000000
#   endif
  endif
  ifeq ($(FC),xlf_r)
    ifeq ($(MP),false)
      MP_FLAGS =
    else
      ifeq ($(BUILDING_BLAS)$(BUILDING_LAPACK),true)
        # No OpenMP in blas and lapack, and lapack contains lines that look
        # like conditional compilation lines.
        # Turn off the omp to turn off conditional compilation
        # and try automatic parallelization.
        ifeq ($(DEBUG),false)
          MP_FLAGS = -qsmp
        else
          MP_FLAGS = -qsmp=noopt
        endif
      else
        ifeq ($(DEBUG),false)
          MP_FLAGS = -qsmp=omp:opt
        else
          MP_FLAGS = -qsmp=omp:noopt
        endif
      endif
    endif

    # 1512-047 (E) The I FORMAT edit descriptor, in literal FMT specifier, is 
    # missing a field width. A default value based on the edit descriptor and 
    # the type of the data item being processed is assumed.
    FFLAGS_ABI   = -q$(ABI) -qextname -qfixed -qsuppress=1512-047 -qarch=auto

    # 1500-036 (I) The NOSTRICT option (default at OPT(3)) has the potential
    # to alter the semantics of a program.  Please refer to documentation on
    # the STRICT/NOSTRICT option for more information.
    FFLAGS_OPT   = $(FFLAGS_ABI) -O3 -qcache=auto -qsuppress=1500-036 -qmaxmem=12000 -qtune=auto

    # 7FF7FFFF is a large integer or NaNQ real*4 or NaNS real*8
    FFLAGS_DEBUG = $(FFLAGS_ABI) -g -qfullpath -C -qflttrap=ov:zero:inv:en -qinitauto=7FF7FFFF

    # Library dependant rules
    # Solver library
    ifeq ($(BUILDING_SOLVER),true)
      FPPFLAGS = -WF,-D'POINTER=INTEGER*$(POINTER_LENGTH)'
    endif

  endif

  ifeq ($(DEBUG),false)
    FFLAGS = $(FFLAGS_OPT) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_OPT) $(MP_FLAGS) $(CINCDIR)
  else
    FFLAGS = $(FFLAGS_DEBUG) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_DEBUG) $(MP_FLAGS) $(CINCDIR)
  endif

  CPPFLAGS =
  LD       = $(FC)
  LDFLAGS  = $(FFLAGS)
  AR       = ar
  AFLAGS   = -rvf
  RANLIB   = ranlib

  SCILIB   = -llapack$(SUFFIX) -lblas$(SUFFIX) -lm
  ifneq ($(filter $(INSTRUCTION),i686 ia64 x86_64),)# i686, ia64, x86_64 
    DIRNAM = $(INSTRUCTION)-linux# no ABI required
  else
    DIRNAM = $(INSTRUCTION)-$(ABI)-linux
  endif

  # Library dependant rules
  # BOXMG library
  ifeq ($(BUILDING_BOXMG),true)
    ifeq ($(FC),xlf_r)
      FPPFLAGS += -WF,-D'POINTER=INTEGER*$(POINTER_LENGTH)'
    else
      FPPFLAGS += -D'POINTER=INTEGER*$(POINTER_LENGTH)'
    endif
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -D_OMP 
  endif

  # Solver library
  ifeq ($(BUILDING_SOLVER),true)
    FPPFLAGS += -D'POINTER=INTEGER*$(POINTER_LENGTH)'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -D_OMP  -DBSD_TIMERS
  endif

  # Config for Umfpack 4.0
  ifeq ($(BUILDING_UMFPACK),true)
    # 1: with no BLAS (this will be slow)
    # CONFIG = -DNBLAS
    # 2: with the C-BLAS
    CONFIG = -DNCBLAS
  endif

  # Config for SuperLU
  ifeq ($(BUILDING_SUPERLU),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = -lblas$(SUFFIX)
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
  endif

  # Config for SuperLU_omp
  ifeq ($(BUILDING_SUPERLU_MP),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = -lblas$(SUFFIX)
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
    # Type of parallelisation: _SGI, _OMP, and type of platform
    CFLAGS  += -D_OMP -D'ARCH=LINUX'
    BLASOPTS = 
  endif

endif

#-----------------------------------------------------------------------------
# win32 Options: use g77.
#-----------------------------------------------------------------------------
ifeq ($(SYSNAME),win32)
  EXPORTS += ; export CFLAGS="-mno-cygwin"; \
        export CXXFLAGS="-mno-cygwin"; \
        export CPPFLAGS=""; \
	export LDFLAGS="-mno-cygwin"
  export $EXPORTS

  MP_FLAGS     =

  POINTER_LENGTH = 4

  # In gcc 3.2 and later this is now -mms-bitfields, before that it 
  # was -fnative-struct
  CFLAGS_ABI   = -ansi -W -Wall -mms-bitfields -mno-cygwin
  #CFLAGS_ABI  = -ansi -pedantic -W -Wall -Wno-parentheses -Wshadow\
  #              -Wcast-align -Winline -Wstrict-prototypes -mms-bitfields
  CFLAGS_OPT   = $(CFLAGS_ABI) -O3
  CFLAGS_DEBUG = $(CFLAGS_ABI) -O0 -g

  # Switch off globals warning (with -Wno-globals) to avoid being told
  # about %VAL(P) being a problem. 
  FFLAGS_ABI   = -fno-second-underscore -Wall -Wsurprising -Wunused\
                 -Wno-globals -mms-bitfields -mno-cygwin
  FFLAGS_OPT   = $(FFLAGS_ABI) -O3
  FFLAGS_DEBUG = $(FFLAGS_ABI) -O0 -g

  ifeq ($(DEBUG),false)
    FFLAGS = $(FFLAGS_OPT) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_OPT) $(MP_FLAGS) $(CINCDIR)
  else
    FFLAGS = $(FFLAGS_DEBUG) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_DEBUG) $(MP_FLAGS) $(CINCDIR)
  endif

  CC           = gcc
  CPPFLAGS     = -DWIN32
  FC           = g77
  FPPFLAGS     = 
  LD           = $(FC)
  LDFLAGS      = $(FFLAGS)
  AR           = ar
  AFLAGS       = -rvf
  RANLIB       = ranlib

  SCILIB       = -llapack$(SUFFIX) -lblas$(SUFFIX) -lm
  DIRNAM       = i386-win32

  # Library dependant rules
  # BOXMG library
  ifeq ($(BUILDING_BOXMG),true)
    FPPFLAGS += -D'POINTER=INTEGER*4'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += 
  endif

  # Solver library
  ifeq ($(BUILDING_SOLVER),true)
    FPPFLAGS += -D'POINTER=INTEGER*4'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -DANSI_C_TIMERS
  endif

  # Config for Umfpack 4.0
  ifeq ($(BUILDING_UMFPACK),true)
    # 1: with no BLAS (this will be slow)
    CONFIG = -DNBLAS
    # 2: with the C-BLAS
    # CONFIG = -DNCBLAS
  endif

  # Config for SuperLU
  ifeq ($(BUILDING_SUPERLU),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = -lblas$(SUFFIX)
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
  endif

  # Config for SuperLU_omp
  ifeq ($(BUILDING_SUPERLU_MP),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = -lblas$(SUFFIX)
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
    # Type of parallelisation: _SGI, _OMP, and type of platform
    CFLAGS += -D_OMP
    BLASOPTS = 
  endif

endif

#-----------------------------------------------------------------------------
# Digital OSF Options
#-----------------------------------------------------------------------------
ifeq ($(SYSNAME),OSF1)

  ifeq ($(MP),false)
    MP_FLAGS =
  else
    MP_FLAGS = -check omp_bindings -omp
  endif

  CFLAGS_ABI   =
  CFLAGS_OPT   = $(CFLAGS_ABI) -O4 -automatic -fast
  CFLAGS_DEBUG = $(CFLAGS_ABI) -O0 -g
  FFLAGS_ABI   = -assume underscore
  FFLAGS_OPT   = $(FFLAGS_ABI) -O4 -automatic -fast
  FFLAGS_DEBUG = $(FFLAGS_ABI) -O0 -g

  ifeq ($(DEBUG),false)
    FFLAGS = $(FFLAGS_OPT) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_OPT) $(MP_FLAGS) $(CINCDIR)
  else
    FFLAGS = $(FFLAGS_DEBUG) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_DEBUG) $(MP_FLAGS) $(CINCDIR)
  endif

  FC       = f77
  FPPFLAGS = 
  CC       = cc
  CPPFLAGS =
  LD       = $(FC)
  LDFLAGS  = $(FFLAGS)
  AR       = ar
  AFLAGS   = uvr
  RANLIB   = ranlib

  SCILIB   = -ldxml -lm
  DIRNAM   = alpha-osf

  # Library dependant rules
  # BOXMG library
  ifeq ($(BUILDING_BOXMG),true)
    FPPFLAGS += -D'POINTER=INTEGER*8'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -D_OMP 
  endif

  # Solver library
  ifeq ($(BUILDING_SOLVER),true)
    FPPFLAGS += -D'POINTER=INTEGER*8'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -D_OMP -DBSD_TIMERS
  endif

  # Config for Umfpack 4.0
  ifeq ($(BUILDING_UMFPACK),true)
    # 1: with the vendor-supplied BLAS
    CONFIG = -DNCBLAS -LP64
    # 2: with no BLAS (this will be slow)
    # CONFIG = -DNBLAS -LP64
  endif

  # Config for SuperLU
  ifeq ($(BUILDING_SUPERLU),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = -ldxml
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
  endif

  # Config for SuperLU_omp
  ifeq ($(BUILDING_SUPERLU_MP),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = -ldxml
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
    # Type of parallelisation: _SGI, _OMP, and type of platform
    CFLAGS += -D_OMP
  endif

endif



#-----------------------------------------------------------------------------
# Darwin Options:
#-----------------------------------------------------------------------------
ifeq ($(SYSNAME),Darwin)

  INSTRUCTION = ppc-32-darwin

  CC = gcc

  POINTER_LENGTH = 4

  FC = g77

  CFLAGS_ABI   = -ansi -pedantic -W -Wall -fPIC
  CFLAGS_OPT   = $(CFLAGS_ABI) -O3
  CFLAGS_DEBUG = $(CFLAGS_ABI) -O0 -g

  # G77 compiler
  ifeq ($(FC),g77)

    # Switch off globals warning (with -Wno-globals) to avoid being told
    # about %VAL(P) being a problem. 
    FFLAGS_ABI   = -fno-second-underscore -W -Wall -Wsurprising -Wunused\
                   -Wno-globals -fPIC
    FFLAGS_OPT   = $(FFLAGS_ABI) -O3
    FFLAGS_DEBUG = $(FFLAGS_ABI) -O0 -g
    MP_FLAGS     =

    FPPFLAGS     = -I.

  endif

  ifeq ($(DEBUG),false)
    FFLAGS = $(FFLAGS_OPT) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_OPT) $(MP_FLAGS) $(CINCDIR)
  else
    FFLAGS = $(FFLAGS_DEBUG) $(MP_FLAGS) $(FINCDIR)
    CFLAGS = $(CFLAGS_DEBUG) $(MP_FLAGS) $(CINCDIR)
  endif

  # SuperLU uses malloc.h directly which isn't in the base directories on Darwin 10.4
  CFLAGS += -I/usr/include/malloc

  CPPFLAGS = 
  LD       = $(FC)
  LDFLAGS  = $(FFLAGS)
  AR       = ar
  AFLAGS   = -rv
  RANLIB   = ranlib

  SCILIB   = -llapack$(SUFFIX) -lblas$(SUFFIX) -lm
  DIRNAM   = $(INSTRUCTION)

  # Library dependant rules
  # BOXMG library
  ifeq ($(BUILDING_BOXMG),true)
    FPPFLAGS += -D'POINTER=INTEGER*$(POINTER_LENGTH)'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -D_OMP 
  endif

  # Solver library
  ifeq ($(BUILDING_SOLVER),true)
    FPPFLAGS += -D'POINTER=INTEGER*$(POINTER_LENGTH)'
    # Type of parallelisation in link to superLU: _SGI, _OMP
    CPPFLAGS += -D_OMP  -DBSD_TIMERS
  endif

  # Config for Umfpack 4.0
  ifeq ($(BUILDING_UMFPACK),true)
    # 1: with no BLAS (this will be slow)
    # CONFIG = -DNBLAS
    # 2: with the C-BLAS
    CONFIG = -DNCBLAS
  endif

  # Config for SuperLU
  ifeq ($(BUILDING_SUPERLU),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = -lblas$(SUFFIX)
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
  endif

  # Config for SuperLU_omp
  ifeq ($(BUILDING_SUPERLU_MP),true)
    # BLAS library to use
    BLASDEF = -DUSE_VENDOR_BLAS
    BLASLIB = -lblas$(SUFFIX)
    # Fortran namespace mapping: -DNoChange, -DAdd_, or -DUpCase
    CDEFS   = -DAdd_
    # Type of parallelisation: _SGI, _OMP, and type of platform
    CFLAGS  += -D_OMP -D'ARCH=DARWIN'
    BLASOPTS = 
  endif

endif



#-----------------------------------------------------------------------------
# Options dependant on the library being built
#-----------------------------------------------------------------------------
ifeq ($(BUILDING_BOXMG),true)
  CMISS_SRC     = $(CMISS_ROOT)/cm/source
  FINCDIR       = -I$(CMISS_SRC)
endif

ifeq ($(BUILDING_SOLVER),true)
  # Set include directories for the libsolver.a library
  ifeq ($(MP),false)
    SUPERLU_DIR = ../superLU-2.0.3/SRC
  else
    SUPERLU_DIR = ../superLU-1.0_omp/SRC
  endif
  UMFPACK4_DIR  = ../umfpack-4.0
  CMISS_SRC     = $(CMISS_ROOT)/cm/source

  FINCDIR       = -I$(CMISS_SRC) -I.
  CINCDIR       = -I$(SUPERLU_DIR) -I$(UMFPACK4_DIR)
endif

# Config for SuperLU
ifeq ($(BUILDING_SUPERLU),true)
  # Compiler options when compiling non-optimised code
  NOOPTS = $(CFLAGS_DEBUG)
endif

# Config for SuperLU_omp
ifeq ($(BUILDING_SUPERLU_MP),true)
  # Compiler options when compiling non-optimised code
  NOOPTS = $(CFLAGS_DEBUG)
endif
 
# Config for LAPACK
ifeq ($(BUILDING_LAPACK),true)
  # Compiler options when compiling non-optimised code
  NOOPTS = $(FFLAGS_DEBUG)
endif

#-----------------------------------------------------------------------------
# End of platform dependant rules
#-----------------------------------------------------------------------------
SUFFIX = $(MT_SUFF)$(DEBUG_SUFF)$(PROF_SUFF)
EXESUF = $(DIRNAM)$(SUFFIX)
SERIAL = $(DEBUG_SUFF)
LIBDIR = ../lib/$(DIRNAM)

ifneq ($(PROF),false)
  LDFLAGS   += -pg
  PROF_SUFF  = -prof
endif

#-----------------------------------------------------------------------------
# Objects to delete in a clean up
#-----------------------------------------------------------------------------
CLEAN_OBJS = *.o *~ rii_files $(LIB)

#-----------------------------------------------------------------------------
# Rule to compile .F files
#-----------------------------------------------------------------------------
COMPILE.F = $(FC) $(FFLAGS) $(FPPFLAGS) $(TARGET_ARCH) -c

