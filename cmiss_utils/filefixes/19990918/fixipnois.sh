#! /bin/sh
#
# The is a shell script to alter ipnois files
#
# Usage:
#   ${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixipnois.sh
#
# Created:
#   John Bodley , 29 May 2000
#
# Updates:
#  
#

#
# Script name (no .sh extension)
#
SCRIPTNAME=fixipnois

#
# String to indicate the types of files to change
# Use a "?" to cover both *.ir* and .*ip* files
#
FILESTOCHANGE=*.i?nois

#
# The awk script to use to change the file
#
AWKSCRIPT=${CMISS_ROOT}/cmiss_utils/filefixes/19990918/$SCRIPTNAME.awk

# ----------------------------------------
# Below this line should not need changing
# ----------------------------------------

#
# Location of the list of files changed
#
LISTOFFILES=/tmp/${SCRIPTNAME}_${LOGNAME}.files

#
# file for temporary file operations
#
tmpfile=/tmp/${SCRIPTNAME}_$LOGNAME.tmp

# Check if necessary files (don't) exist
#
if [ ! -s  $AWKSCRIPT ]; then
  echo "ERROR : Cannot find awk script"
  exit 1
fi
if [ -a  $tmpfile ]; then
  echo "ERROR : Remove $tmpfile"
  exit 1
fi
if [ -a  $LISTOFFILES ]; then
  echo "ERROR : Remove $LISTOFFILES"
  exit 1
fi

# Find the files (excluding links) needing to be changed
#
find . -type f -name "$FILESTOCHANGE" -print > $LISTOFFILES

# Convert the files
#
for file in `cat $LISTOFFILES`
do
  if nawk -f $AWKSCRIPT $file > $tmpfile; then
    if [ -s $tmpfile ]; then
      if cmp -s $file $tmpfile; then
        rm $tmpfile
      else
        echo "Converting file "$file  
        mv $tmpfile $file
      fi
    else 
      echo "ERROR : $tmpfile is of zero size"
      exit 2   
    fi
  else
    echo "ERROR : Awkscript failed.. reverting"
    exit 2
  fi
done
rm $LISTOFFILES
