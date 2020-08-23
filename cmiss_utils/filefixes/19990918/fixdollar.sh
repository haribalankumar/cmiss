#! /bin/sh
#
# The is a shell script to remove the $ character from the first column of *.i* files
# generated under VMS
#
# Usage:
#   ${CMISS_ROOT}/cmiss_utils/filefixes/19990107/fixdollar.sh
#
# Created:
#   Carey Stevens, Wed May 19 14:53:31 NZT 1999
#
# Updates:
#   Karl Tomlinson, 10 December 1999
#

#
# String to indicate the types of files to change
# Use a "?" to cover both *.ir* and .*ip* files
#
#FILESTOCHANGE=*.i*

#
# Location of the list of files changed
#
LISTOFFILES=/tmp/fixdollar_${USER}.files
#
# Find the files (excluding links) needing to be changed
#
  find . -type f -name "*.ip*" -print > $LISTOFFILES
  find . -type f -name "*.ir*" -print >> $LISTOFFILES

echo "List of file to change : $LISTOFFILES"

#
# file for temporary file operations
#
tmpfile=/tmp/fixdollar_$USER.tmp

#
# Convert the files
#
for file in `cat $LISTOFFILES`
do
  if grep -q '^ \$ ' $file
  then
    echo "Converting file "$file  
    sed 's%^ \$ %   %' $file > $tmpfile
    if [ -s $tmpfile ]
    then  mv $tmpfile $file
    else 
      echo "ERROR : result for $file is of zero size"
      exit 2   
    fi
  fi
  if grep -q '^\$' $file
  then
    echo "Converting file "$file  
    sed 's%^\$% %' $file > $tmpfile
    if [ -s $tmpfile ]
    then  mv $tmpfile $file
    else 
      echo "ERROR : result for $file is of zero size"
      exit 2   
    fi
  fi
done
