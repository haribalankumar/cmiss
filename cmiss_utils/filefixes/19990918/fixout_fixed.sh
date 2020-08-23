#! /bin/sh
#
# The is a shell script to alter ***** files
#
# Usage:
#   ${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixipout_fixed.sh
#
# Created:
#   *NAME* , *DATE*
#
# Updates:
#  
#

#
# String to indicate the types of files to change
# Use a "?" to cover both *.ir* and .*ip* files
#
FILESTOCHANGE=cmiss_test.out_fixed

#
# Location of the list of files changed
#
LISTOFFILES=/tmp/fixipout_fixed_${USER}.files

#
# file for temporary file operations
#
tmpfile=/tmp/fixipout_fixed_$USER.tmp

# Check if necessary files (don't) exist
#
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
#  if grep -q '^[^!#]* !' $file; then
  if [ -s $file ]; then
    if sed -e 's%^\([^!#]* \)!%\1#%' -e 's%^ #%> #%' -e 's% >%>%' $file > $tmpfile; then
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
      echo "ERROR : sed failed.. reverting"
      exit 2
    fi
  fi
done
rm $LISTOFFILES
