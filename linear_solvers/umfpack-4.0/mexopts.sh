#
# mexopts.sh	Shell script for configuring MEX-file creation script,
#               mex.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_MEX_OPT: Template Options file for building MEX-files via the system ANSI compiler
#
# Copyright 1984-2000 The MathWorks, Inc.
# $Revision: 1.1.1.1 $  $Date: 2002-03-09 06:12:51 $
#----------------------------------------------------------------------------
#
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The cmex script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        alpha)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS='-ieee -std1'
            CLIBS=''
            COPTIMFLAGS='-O2 -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS='-shared'
            FLIBS='-lUfor -lfor -lFutil'
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld' 
            LDFLAGS="-expect_unresolved '*' -shared -hidden -exported_symbol $ENTRYPOINT -exported_symbol mexVersion"
	    LDOPTIMFLAGS=''
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        hpux)
#----------------------------------------------------------------------------
            CC='cc'
#
# -Wp,-H65535 - works around a compiler limitation so we can compile the
#               MEX version of standalone/compiler/messages.c
#
            CFLAGS='+z -D_POSIX_C_SOURCE=199506L -Wp,-H65535 -Ae'
            CLIBS=""
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f90'
            FFLAGS='+z'
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-b +e $ENTRYPOINT +e mexVersion"
	    LDOPTIMFLAGS=''
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        hp700)
#----------------------------------------------------------------------------
            CC='cc'
#
# -Wp,-H65535 - works around a compiler limitation so we can compile the
#               MEX version of standalone/compiler/messages.c
# +DAportable - remove from CFLAGS if you wish to optimize for target machine
#
            CFLAGS='+z -D_HPUX_SOURCE -Wp,-H65535 -Ae +DAportable'
            CLIBS=''
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f90'
            FFLAGS='+z +DAportable'
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-b +e $ENTRYPOINT +e mexVersion"
	    LDOPTIMFLAGS=''
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        ibm_rs)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS='-qlanglvl=ansi -DIBM_RS'
            CLIBS="-L$MATLAB/bin/$Arch -lmx -lmex -lmatlbmx -lm -lmat"
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS=''
            FLIBS="-lmx -lmex $MATLAB/extern/lib/ibm_rs/fmex1.o -lm -lmat"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='cc'
            LDFLAGS="-L$MATLAB/bin/$Arch -bE:$MATLAB/extern/lib/ibm_rs/$MAPFILE -bM:SRE -e $ENTRYPOINT"
            LDOPTIMFLAGS='-s'
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        glnx86)
#----------------------------------------------------------------------------
            CC='gcc'
            CFLAGS='-fPIC'
            CLIBS=''
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='g77'
            FFLAGS='-fPIC'
            FLIBS='-lf2c'
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD=$CC
            LDFLAGS='-shared -Wl,--version-script,$MATLAB/extern/lib/$Arch/$MAPFILE'
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        sgi)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS='-n32 -mips3'
            CLIBS=''
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS='-n32 -mips3'
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-shared -exported_symbol $ENTRYPOINT -exported_symbol mexVersion"
	    LDOPTIMFLAGS=''
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        sol2)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS='-KPIC -dalign'
            CLIBS=''
	    COPTIMFLAGS='-xO5 -native -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS='-dalign -KPIC'
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='/usr/ccs/bin/ld'
            LDFLAGS="-G -M $MATLAB/extern/lib/sol2/$MAPFILE"
	    LDOPTIMFLAGS=''
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
