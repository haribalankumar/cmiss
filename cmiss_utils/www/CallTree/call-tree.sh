#! /bin/sh
#
#  Call Tree shell script
#


# -- Resource Locations (and defaults) ----------------------------------------
#

if [ "$CMISS_ROOT" = "" ]; then
  echo \$CMISS_ROOT variable unbound.
  exit 1
else
  CALLTREE_PROGRAM=$CMISS_ROOT/cmiss_utils/www/CallTree/call-tree
  MODULE_LIST=$CMISS_ROOT/cmiss_utils/www/MasterLists/module-list
fi



# -- Access/Permission Checking -----------------------------------------------
#

if [ ! -x $CALLTREE_PROGRAM ]; then
  echo Error: $CALLTREE_PROGRAM not found or not executable
  exit 1
fi
  
if [ ! -r $MODULE_LIST ]; then
  echo Error: $MODULE_LIST not found or not readable
  exit 1
fi



# -- Run Program --------------------------------------------------------------
#

$CALLTREE_PROGRAM -l $MODULE_LIST



# -- Return To Caller ---------------------------------------------------------
#

exit
