#!/bin/sh
#
# The is a shell script to add another question to ipequa files
#
# Usage:
#   ${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixcom.sh
#
# Created:
#   Shane Blackett 7/3/99
#

FILESTOCHANGE=*.com
FILESTOCHANGE2=*_com.cmiss

#
# The awk script to use to change the file
#
AWKSCRIPT=${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixcom.awk
#
# Location of the list of files changed
#
LISTOFFILES=/tmp/fixcom_${LOGNAME}.files

#
# Location of the log file recording script usage
#
LOGFILE=${CMISS_ROOT}/cmiss_utils/filefixes/19990918/fixcom.log 
#
# The name of this script as written into the comfiles when changes
# are made.
#
SCRIPTNAME=fixcom.sh

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
if [ -a  $LISTOFFILES ]; then
  echo "ERROR : Remove $LISTOFFILES"
  exit 1
fi

# Find the files (excluding links) needing to be changed
#
find . -type f -name "$FILESTOCHANGE" -o -type f -name "$FILESTOCHANGE2" > $LISTOFFILES

# Convert the files
#
for file in `cat $LISTOFFILES`
do
  echo -n "Converting file "$file" ...."

  if [ ! -s $file ]; then
    echo " empty - no change."
  else if nawk -f $AWKSCRIPT -v DATE="`date`" -v PROGRAM_NAME=$SCRIPTNAME $file > $file.tmp; then

    if [ ! -s $file.tmp ]; then
      echo "\nERROR : Corrected $file is zero size.. reverting"
      rm $file.tmp
      exit 2
    else
#        diff $file $file.tmp > $file.diff ;
#        if [ ! -s $file.diff ]; then
      if cmp -s $file $file.tmp; then
        rm $file.tmp
        echo " no changes required."
        if [ -w $LOGFILE ]; then
	  echo `whoami` ": Parsed  " $PWD"/"$file " " `date` >> $LOGFILE
        fi
      else
        mv $file $file.old
        mv $file.tmp $file
        echo " done."
        if [ -w $LOGFILE ]; then
	  echo `whoami` ": Updated " $PWD"/"$file " " `date` >> $LOGFILE
        fi
      fi
#      rm $file.diff
    fi
  else
    echo "\nERROR : Awkscript failed.. reverting"
    exit 2
  fi fi
done
rm $LISTOFFILES
