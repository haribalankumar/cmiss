#! /bin/sh
#
# Shell script to find mismatched enters and exits calls
#
# Usage:
#   findenterexit.sh
#
# Created:
#   Chris Bradley and Glen Harris, 31 October 1996
#
# Updates:
#   Karl Tomlinson, 17 September 1999:
#     Removed <UL> of files for error counting.
#   Karl Tomlinson, 28 April 2000:
#     Allowed ROUTINENAME parameter for name.

outfile=$CMISS_ROOT/www/help/errors/findenterexit.html

{
  echo "<HTML><HEAD><TITLE>CMISS Enters/Exits errors</TITLE></HEAD>"
  echo "<H1>CMISS Enters/Exits errors</H1>"
  echo "<UL>"

  module_list=$CMISS_ROOT/cmiss_utils/www/MasterLists/module-list

  for file in `cat $module_list`; do
      egrep -i '^ *[0-9]* +(SUBROUTINE|CALL ENTERS|CALL EXITS|RETURN|END *(\!|$)|PARAMETER\(ROUTINENAME)' $file | awk -f ${CMISS_ROOT}/cmiss_utils/www/findenterexit/findenterexit.awk
  done

  echo "</UL>"
} > "$outfile"

${CMISS_ROOT}/html_utils/addfooter.sh programmer $outfile

