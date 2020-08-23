#!/bin/sh
#
#  Run mail-user.sh
#
#    Syntax: mail-user.sh <username> <email>
#


# -- Module Constants ---------------------------------------------------------
#
ROOT=${HOME}/webcmiss
HOST=`hostname`
OS=`uname`

case $OS in
  "AIX"|"Linux")
    MAIL=/usr/bin/Mail
    ;;
  "IRIX"|"IRIX64")
    MAIL=/usr/sbin/Mail
    ;;
  *)
    MAIL=/usr/bin/Mail
    ;;
esac

# -- Module Entry Point -------------------------------------------------------
#
if [ $# -ne 2 ] ; then
  echo "Syntax: mail-user.sh <username> <email>"
  exit 1
fi

# Read parameters into variables
#
USERNAME=$1
EMAIL=$2

# Generate a list of the output files.
#
cd $ROOT/people/$USERNAME/Output
OUTPUT=`ls`

# Mail the user telling them that their job is complete. Include
# a list of the output files generated.
$MAIL -s "CMISS Run Completed on $HOST" $EMAIL <<BLORT
Your CMISS job on $HOST has completed it's run.

It produced the following output files:
$OUTPUT

You can view the files via the CMISS submission system at:

  http://www.bioeng.auckland.ac.nz/cmiss/hpc_access.php
BLORT
exit 0
