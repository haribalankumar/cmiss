#! /bin/sh
#
# The is a shell script to add another question to ipopti files
#
# Usage:
#   ${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixipopti.sh
#
# Created:
#   Leo Cheng, 17 November 1999
#
# Updates:
#  
#

#
# String to indicate the types of files to change
# Use a "?" to cover both *.ir* and .*ip* files
#
FILESTOCHANGE=*.i?opti

#
# The awk script to use to change the file
#
AWKSCRIPT=${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixipopti.awk

#
# Location of the list of files changed
#
LISTOFFILES=/tmp/fixipopti_${LOGNAME}.files



# ----------------------------------------
# Below this line should not need changing
# ----------------------------------------

echo "List of files to change : $LISTOFFILES"

# Check if necessary files exist
#
if [ ! -s  $AWKSCRIPT ]; then
  echo "ERROR : Cannot find awk script"
  exit 1
fi

# Find the files (excluding links) needing to be changed
#
find . -type f -name "$FILESTOCHANGE" -print > $LISTOFFILES

# Convert the files
#
for file in `cat $LISTOFFILES`
do
  echo "Converting file "$file  
  nawk -f $AWKSCRIPT $file > $file"_tmp"

  if [ -s $file"_tmp" ]
  then  mv $file"_tmp" $file
  else 
    echo "ERROR : $file"_tmp" is of zero size"
    exit 2   
  fi
done
