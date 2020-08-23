# For use with GNU make.
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

ifndef DEBUG
  ifndef OPT
    OPT := false
  endif
  ifeq ($(OPT),false)
    DEBUG := true
  else
    DEBUG := false
  endif
endif

# set architecture dependent directories and default options

# defaults
INSTRUCTION=$(MACHNAME)
BIN_ARCH_DIR = $(INSTRUCTION)-$(OPERATING_SYSTEM)
LIB_ARCH_DIR = $(INSTRUCTION)-$(ABI)-$(OPERATING_SYSTEM)

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
  INSTRUCTION := mips
  OPERATING_SYSTEM := irix
endif
ifeq ($(SYSNAME),Linux)
  OPERATING_SYSTEM := linux
  LIB_ARCH_DIR = $(INSTRUCTION)-$(OPERATING_SYSTEM)# no ABI
  ifeq ($(filter-out i%86,$(MACHNAME)),)
    INSTRUCTION := i686
  endif
  ifndef ABI
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      ABI=32
    endif
    ifeq ($(MACHNAME),ia64)
      ABI=64
    endif
  endif
endif
ifeq ($(SYSNAME:CYGWIN%=),)# CYGWIN
  SYSNAME := win32
endif
ifeq ($(SYSNAME),win32)
  LIB_ARCH_DIR = $(INSTRUCTION)-$(OPERATING_SYSTEM)# no ABI
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
    else
      ABI = 32
    endif
  endif
  INSTRUCTION := rs6000
  OPERATING_SYSTEM := aix
endif

ifneq ($(DEBUG),false)
  DEBUG_SUFFIX = -debug
else
  DEBUG_SUFFIX =
endif

#Make architecture directory names and lib name
SOURCE_DIR = source
WORKING_DIR := generated/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)

BIN_DIR := ../bin/$(BIN_ARCH_DIR)
BIN_NAME := bin2base64h$(DEBUG_SUFFIX)
BIN := $(BIN_DIR)/$(BIN_NAME)

C_SOURCES := bin2base64h.c
C_UNITS := $(basename $(C_SOURCES) )
DEPEND_FILES := $(foreach unit, $(C_UNITS), $(WORKING_DIR)/$(unit).d )


#-----------------------------------------------------------------------------
# compiling commands

CC = cc
LD = $(CC)
# Option lists
# (suboption lists become more specific so that later ones overrule previous)
CFLAGS = $(strip $(CFL_FLGS) $(CFE_FLGS) $(CF_FLGS))
CPPFLAGS := $(addprefix -I, $(C_INCLUDE_DIRS) ) '-DABI_ENV="$(ABI_ENV)"'
ARFLAGS = -cr
ifneq ($(DEBUG),false)
  CFLAGS += $(strip $(DBGCF_FLGS) $(DBGC_FLGS))
else
  CFLAGS += $(strip $(OPTCFE_FLGS) $(OPTCF_FLGS) $(OPTC_FLGS))
endif
# suboption lists
CFL_FLGS =#	flags for C fortran and linking
L_FLGS =#	flags for linking only
CFE_FLGS =#	flags for C fortran and linking executables only
CF_FLGS = -c#	flags for C and fortran only
DBGCF_FLGS = -g#OPT=false flags for C and fortran
DBGC_FLGS =#	OPT=false flags for C only
OPTCFE_FLGS =#	OPT=true flags for C and fortran and linking executables
OPTCF_FLGS = -O#OPT=true flags for C and fortran only
OPTC_FLGS =#	OPT=true flags for C only

ifeq ($(filter-out IRIX%,$(SYSNAME)),)# SGI
  CF_FLGS += -use_readonly_const
  DBGCF_FLGS += -DEBUG:trap_uninitialized:subscript_check:verbose_runtime
  # warning 158 : Expecting MIPS3 objects: ... MIPS4.
  L_FLGS += -rdata_shared -DEBUG:error=158 -woff 47
  CFL_FLGS = -$(ABI) -mips$(MIPS)
  OPTCF_FLGS = -O3 -OPT:Olimit=8000
endif
ifeq ($(SYSNAME),Linux)
  ifeq ($(MACHNAME),ia64)
    CC = icc
    CFLAGS += -w2# -Wall
# This doesn't seem to do anything
#     ifeq ($(ABI),64)
#       CF_FLGS += -size_lp64
#     endif
  else
    CPPFLAGS += -Dbool=char -DHAS_BOOL
  endif
  OPTCF_FLGS = -O2
  SHARED_LINK_LIBRARIES += -lcrypt -ldl
endif
ifeq ($(SYSNAME),SunOS)
  # need -xarch=native after -fast
  OPTCFE_FLGS += -fast $(CFE_FLGS)
  ifeq ($(ABI),64)
    CFE_FLGS += -xarch=native64
  endif
endif
ifeq ($(SYSNAME),win32)
  CC = gcc -fnative-struct
endif
ifeq ($(SYSNAME),AIX)
  CC = xlc
  # 1506-743 (I) 64-bit portability: possible change of result through conversion ...
  # FD_SET in sys/time.h does this
  # no -qinfo=gen because perl redefines many symbols
  CFLAGS += -qinfo=ini:por:pro:trd:tru:use -qsuppress=1506-743
  ARFLAGS += -X$(ABI)
  # may want -qsrcmsg
  CF_FLGS += -qfullpath
  CFE_FLGS += -q$(ABI) -qarch=auto
  L_FLGS += -b$(ABI)
  ifeq ($(ABI),64)
    CF_FLGS += -qwarn64
  endif
  OPTCF_FLGS = -O3 -qmaxmem=12000 -qtune=auto
  OPTC_FLGS += -qnoignerrno
endif
ifeq ($(SHARED_OBJECT), true)
  CPPFLAGS += -DSHARED_OBJECT
endif
ifeq ($(USE_DYNAMIC_LOADER), true)
  CPPFLAGS += -DUSE_DYNAMIC_LOADER
endif
ifneq ($(SHARED_OBJECT), true)
  ifeq ($(STATIC_PERL_LIB),)
    CPPFLAGS += -DNO_STATIC_FALLBACK
  endif
endif
SHARED_LINK_LIBRARIES += -lc
.PHONY : main

vpath $(PERL) $(subst :, ,$(PATH))

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

  main : $(WORKING_DIR) $(BIN_DIR) $(DEPEND_FILES) $(BIN)

  $(WORKING_DIR) :
	mkdir -p $@

  $(BIN_DIR) :
	mkdir -p $@

clean:
	@echo "Cleaning house ..."
	-rm -rf $(WORKING_DIR) $(BIN_DIR)

  # include the depend file dependencies
  ifneq ($(DEPEND_FILES),)
    sinclude $(DEPEND_FILES)
  endif

  # implicit rules for making the dependency files

  # KAT I think Solaris needed nawk rather than awk, but nawk is not usually
  # avaiable on Mandrake.  I don't have a Sun to try this out so I'll get it
  # working with awk on the machines I have.
$(WORKING_DIR)/%.d : $(SOURCE_DIR)/%.c $(WORKING_DIR)
	makedepend $(CPPFLAGS) -f- -Y $< 2> $@.tmp | sed -e 's%^source\([^ ]*\).o%$$(WORKING_DIR)\1.o $$(WORKING_DIR)\1.d%' > $@

$(BIN) : $(foreach unit, $(C_UNITS), $(WORKING_DIR)/$(unit).o )
	$(LD) -o $@ $^


# implicit rules for making the objects
$(WORKING_DIR)/%.o : $(SOURCE_DIR)/%.c
	$(CC) -o $@ $(CPPFLAGS) $(CFLAGS) $<

