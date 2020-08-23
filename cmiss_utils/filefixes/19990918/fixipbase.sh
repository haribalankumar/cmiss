#! /bin/sh
#
# A shell script to add another question to CMISS .ipbase files
#
# Usage:
#   ${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixipbase.sh
#
# Created:
#   Karl Tomlinson, 14 December 1999
#
# Updates:
#  
#

#
# String to indicate the types of files to change
# Use a "?" to cover both *.ir* and .*ip* files
#
FILESTOCHANGE=*.i?base

#
# The awk script to use to change the file
#
AWKSCRIPT=${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixipbase.awk

#
# Location of the list of files changed
#
LISTOFFILES=/tmp/fixipbase_${USER}.files



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
if [ -a  $tmpfile ]; then
  echo "ERROR : Remove $tmpfile"
  exit 1
fi

# Find the files (excluding links) needing to be changed
#
find . -type f -name "$FILESTOCHANGE" -print > $LISTOFFILES

# Convert the files
#
for file in `cat $LISTOFFILES`
do
  nawk -f $AWKSCRIPT $file > $tmpfile

  if [ -s $tmpfile ]
  then
    if cmp -s $file $tmpfile
    then
      rm $tmpfile
    else
      echo "Converting file "$file  
      mv $tmpfile $file
    fi
  else 
    echo "ERROR : $tmpfile is of zero size"
    exit 2   
  fi
done
rm $LISTOFFILES