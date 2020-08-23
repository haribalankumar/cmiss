#! /bin/sh

# Generates, for emacs, tags on cmiss source below the current directory.

curdir=`pwd -P`

# Make sure CMISS_ROOT is an absolute path
if [ "$CMISS_ROOT" ]; then
  CMISS_ROOT=`cd "$CMISS_ROOT" && pwd -P`
fi

case "$curdir" in
  ${CMISS_ROOT:-}*) # Below CMISS_ROOT; do nothing.
    ;;
  # Not below CMISS_ROOT; include global file.
  */cm/source)
    include=$CMISS_ROOT/cm/source/TAGS
    ;;
esac

# See if we have Exuberant Ctags
if command -v exuberant-ctags >/dev/null 2>&1
  then
    exuberant_ctags=exuberant-ctags
elif ctags --version 2>/dev/null | grep 'Exuberant Ctags' >/dev/null
  then
    exuberant_ctags=ctags
fi

if [ "$exuberant_ctags" ]
  then # Exuberant Ctags
    {
      # Listing .cmn and .inc files last so that, if there are documentation
      # tags in the .f files (where the variables are set up), they are found
      # before the definition of variables in the included files.
      find . -name '*.[cf]' | sort && \
	  find . -name '*.h' -o -name '*.inc' -o -name '*.cmn' | sort
    } | "$exuberant_ctags" -e --file-scope=no --langmap=fortran:+.cmn.inc \
	--regex-fortran='/^C#### *(Command|Comment|Module|Variable): *([^(]+)/\2/' \
	${include+"--etags-include=$include"} ${1+"$@"} -L -

else # use GNU etags

    # C: may want --member
    find . -name '*.c' | sort \
	| etags --no-globals --no-defines ${1+"$@"} - && \
	find . -name '*.h' | sort \
	| etags --append ${1+"$@"} - && \
	find . -name '*.f' -o -name '*.inc' -o -name '*.cmn' | sort \
	| etags --append --language=fortran \
	--regex='/C#### *\(Command\|Comment\|Module\|Variable\): *\([^(]+\)/\2/' \
	${include+"--include=$include"} ${1+"$@"} -

fi
