#! /bin/sh
#
#  Index generation script for the CMISS routine browser CGI program
#


# -- Resource Locations (and defaults) ----------------------------------------
#

if [ "$CMISS_ROOT" = "" ]; then
  echo \$CMISS_ROOT variable unbound.
  exit 1
else
  INDEX_PROGRAM=$CMISS_ROOT/cmiss_utils/www/RoutineBrowser/make-index
  MODULE_LIST=$CMISS_ROOT/cmiss_utils/www/MasterLists/module-list-browser
  INDEX_DIR=$CMISS_ROOT/www/help/cm/routines
  INDEX_OUPUT=$INDEX_DIR/routine-index
fi



# -- Access/Permission Checking -----------------------------------------------
#

if [ ! -x $INDEX_PROGRAM ]; then
  echo Error: $INDEX_PROGRAM not found or not executable
  exit 1
fi

if [ ! -r $MODULE_LIST ]; then
  echo Error: $MODULE_LIST not found or not readable
  exit 1
fi

if [ ! -w $INDEX_OUPUT ]; then
  if [ ! -d $INDEX_DIR ]; then
    echo Error: directory $INDEX_DIR not found
    exit 1
  else
    if [ ! -w $INDEX_DIR ]; then
      echo Error: directory $INDEX_DIR not writable
      exit 1
    else
      if [ -f $INDEX_OUPUT ]; then
        echo Error: $INDEX_OUPUT not writable
        exit 1
      fi
    fi
  fi
fi
  


# -- Run Program --------------------------------------------------------------
#

$INDEX_PROGRAM -m $MODULE_LIST -o $INDEX_OUPUT



# -- Return To Caller ---------------------------------------------------------
#

exit
