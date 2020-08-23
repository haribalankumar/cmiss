#! /bin/sh
#
#  Parse Fortran Check shell script
#


# -- Resource Locations (and defaults) ----------------------------------------
#

if [ "$CMISS_ROOT" = "" ]; then
  echo \$CMISS_ROOT variable unbound.
  exit 1
else
  PARSE_PROGRAM=$CMISS_ROOT/cmiss_utils/www/FCheck/parse-fortran-check
  COMMON_LIST=$CMISS_ROOT/cmiss_utils/www/MasterLists/common-list
  HTML_DIR=$CMISS_ROOT/www/help/errors
  HTML_OUTPUT=$HTML_DIR/common-errors.html
  FCHECK_OUT=$CMISS_ROOT/cmiss_utils/www/fortran_check/checkcmiss.out
fi



# -- Access/Permission Checking -----------------------------------------------
#

if [ ! -x $PARSE_PROGRAM ]; then
  echo Error: $PARSE_PROGRAM not found or not executable
  exit 1
fi
  
if [ ! -r $COMMON_LIST ]; then
  echo Error: $COMMON_LIST not found or not readable
  exit 1
fi

if [ ! -r $FCHECK_OUT ]; then
  echo Error: $FCHECK_OUT not found or not readable
  exit 1
fi

if [ ! -e $HTML_OUTPUT ]; then
	touch $HTML_OUTPUT
fi

if [ ! -w $HTML_OUTPUT ]; then
  if [ ! -d $HTML_DIR ]; then
    echo Error: directory $HTML_DIR not found
    exit 1
  else
    if [ ! -w $HTML_DIR ]; then
      echo Error: directory $HTML_DIR not writable
      exit 1
    else
      echo Error: $HTML_OUTPUT not writable
      exit 1
    fi
  fi
fi



# -- Run Program --------------------------------------------------------------
#

$PARSE_PROGRAM -c $COMMON_LIST -f $FCHECK_OUT -o $HTML_OUTPUT



# -- Return To Caller ---------------------------------------------------------
#

exit
