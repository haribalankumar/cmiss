#!/bin/sh
#
# The is a shell script to add another question to ipequa files
#
# Usage:
#   ${CMISS_ROOT}/cmiss_utils/filefixes/19990107/fixipequa.sh
#
# Created:
#   Chris Bradley, 3/2/99
#
# Updates:

#TVK 31/05/2000 added option to finite elasticity problem type question
#KAT 31/05/2000 included commit to CVS repository.
#KAT 13/06/2000 removed CVS stuff


#
# String to indicate the types of files to change
# Use a "?" to cover both *.ir* and .*ip* files
#
FILESTOCHANGE=*.i?equa

#
# The awk script to use to change the file
#
AWKSCRIPT=${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixipequa.awk
# Location of the list of files changed
#
LISTOFFILES=/tmp/fixipequa_${LOGNAME}.files

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

  if [ -s $file"_tmp" ]; then
    mv $file"_tmp" $file
# KAT: I think it's better to let the user decide when to commit.
# For updating examples, a better plan might be to find and checkout,
# then run this script.
#      if [ -d `dirname $file`/CVS ] && cvs update $file; then
#        cvs commit -m "fixipequa.sh: added option to finite elasticity problem type question" $file
#      fi
  else 
    echo "ERROR : $file"_tmp" is of zero size"
    exit 2   
  fi
done
