#! /bin/sh
#
# A shell script to add another option to CMISS .ipmate files
#
# Usage:
#   ${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixipmate.sh
#
# Created:
#   Carey Stevens, Wed Mar  1 18:44:12 NZDT 2000
#
# Updates:
#  
#

#
# String to indicate the types of files to change
# Use a "?" to cover both *.ir* and .*ip* files
#
FILESTOCHANGE=*.i?mate

#
# The awk script to use to change the file
#
AWKSCRIPT=${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixipmate.awk
AWKSCRIPT2=${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixipmate2.awk

#
# Location of the list of files changed
#
LISTOFFILES=/tmp/fixipmate_${USER}.files



# ----------------------------------------
# Below this line should not need changing
# ----------------------------------------

#
# file for temporary file operations
#
tmpfile=/tmp/fixipbase_$USER.tmp

# Check if necessary files (don't) exist
#
if [ ! -s  $AWKSCRIPT ]; then
  echo "ERROR : Cannot find awk script"
  exit 1
fi
if [ ! -s  $AWKSCRIPT2 ]; then
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
# KAT $ACKSCRIPT cannot be run more than once
#    nawk -f $AWKSCRIPT $file > $tmpfile

#    if [ -s $tmpfile ]
#    then
#      if cmp -s $file $tmpfile
#      then
#        rm $tmpfile
#      else
#        echo "Converting file "$file phase 1
#        mv $tmpfile $file
#      fi
#    else 
#      echo "ERROR : $tmpfile is of zero size"
#      exit 2   
#    fi

  if nawk -f $AWKSCRIPT2 $file > $tmpfile; then
    if [ -s $tmpfile ]; then
      if cmp -s $file $tmpfile; then
        rm $tmpfile
      else
        echo "Converting file "$file phase 2
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
