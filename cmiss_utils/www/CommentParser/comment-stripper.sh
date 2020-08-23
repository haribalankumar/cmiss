#! /bin/sh
#
#  Comment Stripper shell script
#


# -- Resource Locations (and defaults) ----------------------------------------
#

if [ "$CMISS_ROOT" = "" ]; then
  echo \$CMISS_ROOT variable unbound.
  exit 1
else
  STRIPPER_PROGRAM=$CMISS_ROOT/cmiss_utils/www/CommentParser/comment-stripper
  MODULE_LIST=$CMISS_ROOT/cmiss_utils/www/MasterLists/module-list
  C_MODULE_LIST=$CMISS_ROOT/cmiss_utils/www/MasterLists/c-module-list
  RULE_FILE=${0%/*}/cmiss-comment.rules
fi



# -- Access/Permission Checking -----------------------------------------------
#

if [ ! -x $STRIPPER_PROGRAM ]; then
  echo Error: $STRIPPER_PROGRAM not found or not executable
  exit 1
fi
  
if [ ! -r $MODULE_LIST ]; then
  echo Error: $MODULE_LIST not found or not readable
  exit 1
fi

if [ ! -r $C_MODULE_LIST ]; then
  echo Error: $C_MODULE_LIST not found or not readable
  exit 1
fi

if [ ! -r $RULE_FILE ]; then
  echo Error: $RULE_FILE not found or not readable
  exit 1
fi

# -- Run Program --------------------------------------------------------------
#

$STRIPPER_PROGRAM -r $RULE_FILE -l $MODULE_LIST -c $C_MODULE_LIST



# -- Return To Caller ---------------------------------------------------------
#

exit
