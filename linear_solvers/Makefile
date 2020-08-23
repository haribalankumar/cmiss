# -*- Makefile -*-
# 
# For use with GNU make.

#-----------------------------------------------------------------------------
# Makefile for building the CMISS solver library, which contains iterative
# linear solvers, and the LAPACK, Umfpack and SuperLU direct solvers.
#
# Original by Stuart Norris 2002-10-08
# Changes:
#-----------------------------------------------------------------------------

include Makefile.inc

#-----------------------------------------------------------------------------
# End of initialisation: set targets and rules.
#-----------------------------------------------------------------------------

# libamg is not built by default as the terms of its Licence are unclear.
DIR_BUILD = libboxmg libxblas libsolver libumfpack
ifeq ($(MP),false)
  DIR_BUILD += libsuperlu
else
  DIR_BUILD += libsuperlu_omp
endif

ifeq ($(SYSNAME),AIX)
  DIR_BUILD += liblapack
endif
ifeq ($(SYSNAME),Linux)
  ifndef ACMLPATH
    ifndef MKLPATH
      ifeq (,$(wildcard /usr/lib*/libblas.* /usr/lib/libscs.*))
        DIR_BUILD += libblas
      endif
      ifeq (,$(wildcard /usr/lib*/libcblas.* /usr/lib/libscs.*))
        DIR_BUILD += libcblas
      endif
      ifeq (,$(wildcard /usr/lib*/liblapack.* /usr/lib/libscs.*))
        DIR_BUILD += liblapack
      endif
    endif
  endif
endif
ifeq ($(SYSNAME),win32)
  DIR_BUILD += libblas libcblas liblapack
endif
ifeq ($(SYSNAME),Darwin)
  DIR_BUILD += libblas libcblas liblapack
endif

default: $(DIR_BUILD)

# Can't make different versions in parallel because they use the same object
# files.
.NOTPARALLEL:

  # DEBUG,ABI,MP could be available from export but being verbose
  # so user can see what is happening
debug debugmp opt optmp debug64 debug64mp opt64 opt64mp:
	$(MAKE) --no-print-directory OPT=$(OPT) ABI=$(ABI) MP=$(MP) PROF=$(PROF)

debug debugmp debug64 debug64mp: OPT=false
opt optmp opt64 opt64mp: OPT=true
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
ifeq (,$(filter $(SYSNAME),Linux win32))
  allsp: debug64 opt64
  allmp: debug64mp opt64mp
  all: allmp
endif
ifeq ($(MACHNAME),ia64)
  all: allmp
endif

libsolver:
	( cd solver &&          $(MAKE) MP=$(MP) OPT=$(OPT) ABI=$(ABI) )

libsuperlu_omp:
	( cd superLU-1.0_omp && $(MAKE) MP=true OPT=$(OPT) ABI=$(ABI) )

libsuperlu:
	( cd superLU-2.0.3 &&   $(MAKE) MP=false OPT=$(OPT) ABI=$(ABI) )

libumfpack:
	( cd umfpack-4.0 &&     $(MAKE) MP=false OPT=$(OPT) ABI=$(ABI) )

libamg:
	( cd amg1r6 &&          $(MAKE) MP=$(MP) OPT=$(OPT) ABI=$(ABI) )

libboxmg:
	( cd boxmg-1.6.4 &&     $(MAKE) MP=$(MP) OPT=$(OPT) ABI=$(ABI) )

libxblas:
	( cd xblas &&           $(MAKE) MP=$(MP) OPT=$(OPT) ABI=$(ABI) )

libblas:
	( cd blas &&            $(MAKE) MP=$(MP) OPT=$(OPT) ABI=$(ABI) )

libcblas:
	( cd cblas &&           $(MAKE) MP=$(MP) OPT=$(OPT) ABI=$(ABI) )

liblapack:
	( cd lapack-3.0 &&      $(MAKE) MP=$(MP) OPT=$(OPT) ABI=$(ABI) )

#-----------------------------------------------------------------------------
# Make directories and clean up
#-----------------------------------------------------------------------------

LIBDIR = lib/$(DIRNAM)

libsolver libsuperlu_omp libsuperlu libumfpack libamg libboxmg: $(LIBDIR)
libxblas libcblas libblas liblapack: $(LIBDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

clean:
	( cd solver &&          $(MAKE) clean )
	( cd superLU-1.0_omp && $(MAKE) clean )
	( cd superLU-2.0.3 &&   $(MAKE) clean )
	( cd umfpack-4.0 &&     $(MAKE) clean )
	( cd amg1r6 &&          $(MAKE) distclean )
	( cd xblas &&           $(MAKE) clean )
	( cd blas &&            $(MAKE) clean )
	( cd cblas &&           $(MAKE) clean )
	( cd lapack-3.0 &&      $(MAKE) clean )

clobber: clean
	rm -rf $(LIBDIR)/lib*.a
