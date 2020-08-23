#! /bin/sh 
#
# Shell file to run cmgui fronm the HPC Internet Submission System
#
# Usage: 
#   hpc_cmgui.sh tar_filename
# Created: 
#   Carey Stevens, 23 July 1997
# Updates: 
#

IDENT=`basename $1 .cgi`
DIR=/usr/tmp/$IDENT
TARFILE=$1
RUNFILE="$IDENT.com"

mkdir $DIR
cd $DIR

tar xf $TARFILE
COMFILE=`ls *.com`
echo "set directory /usr/tmp/$IDENT" >| $RUNFILE
cat $COMFILE >> $RUNFILE
cmgui $RUNFILE > /dev/null

cd /usr/tmp
rm -rf $DIR

