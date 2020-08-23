#!/bin/sh
#
# CMISS final daily crontab commands.  Executed in the morning
# hopefully after all tests have run so that it should reflect the
# state of all the overnight testing including linux and hpc.
#

# Set up PATH and other environment as for login shell
. $HOME/.profile

if [ -z "$CMISS_ROOT" ]; then
    echo CMISS_ROOT not set 1>&2
    exit 1
fi

# Day of the week
weekday=`date +%A`

# All STDOUT is sent to a log file.
# Only STDERR that is expected is sent to this file.
# In this way, if there is an unexpected error, a mail message is sent.

logdir=$CMISS_ROOT/cmiss_utils/cronjobs
name=final
daylog=$logdir/$name.day.out
dayerr=$logdir/$name.day.err
weeklog=$logdir/$name.week.log
lastweeklog=$logdir/$name.lastweek.log
 
# Back up week's log
[ "$weekday" = Sunday ] && mv "$weeklog" "$lastweeklog"
# Back up log file in week's log
cat "$daylog" >> "$weeklog"
# Include any stderr from the previous day in the week log.
# Use the first line of the stdout as a record of when stderr was generated.
if [ -s "$dayerr" ]; then
    {
      echo STDERR from `head -1 "$daylog"`
      cat "$dayerr"
      echo END OF STDERR
      echo
    } >> "$weeklog"
fi

# If there is stderr, we want to record it in a log file but also mail
# it.  To do this first save it in a log file then send it to stderr,
# which will be mailed by the cron daemon.

{

startdate=`date`
# This first line is used when stderr is recorded.
echo Cronjob started "$startdate":
echo

# Check cmgui and unemap builds are up to date

#
# CMGUI compilation jobs
#

echo "Ensuring cmgui is up to date everywhere ..."

make -C "$CMISS_ROOT/cmgui" cronjob 2>&1
# Generate changes list
${CMISS_ROOT}/bin/cmisschanges.pl -tree cmgui/source -out html -since week > ${CMISS_ROOT}/cmgui/utilities/cmgui_cvs_changes.log

# note that cmgui executables are only updated if all examples passed.

#
# Zinc compilation jobs
#

echo "Ensuring zinc is up to date everywhere ..."

make -C "$CMISS_ROOT/zinc" cronjob 2>&1

#   
# UNEMAP compilation jobs
# 
echo "Ensuring unemap is up to date..."
make -C "$CMISS_ROOT/unemap" cronjob 2>&1

# update unemap executables
echo "Copying unemap..." 2>&1

(
  cd $CMISS_ROOT || exit 1
  update_if_newer () {
      # Copy $1 to $2 if it is newer
      if [ "$1" -nt "$2" ]; then
          mv "$2" "$2.save"
          cp "$1" "$2"
      fi
  } 

  update_if_newer unemap/bin/Unemap bin/Unemap
   
  update_if_newer unemap/bin/mips-irix/unemap bin/mips-irix/unemap
  update_if_newer unemap/bin/mips-irix/unemap-debug bin/mips-irix/unemap-debug
  update_if_newer unemap/bin/mips-irix/unemap-3d bin/mips-irix/unemap-3d
  update_if_newer unemap/bin/mips-irix/unemap64 bin/mips-irix/unemap64

  update_if_newer unemap/bin/i686-linux/unemap bin/i686-linux/unemap
  update_if_newer unemap/bin/i686-linux/unemap-3d bin/i686-linux/unemap-3d

  update_if_newer unemap/bin/x86_64-linux/unemap bin/x86_64-linux/unemap
  update_if_newer unemap/bin/x86_64-linux/unemap-3d bin/x86_64-linux/unemap-3d

  update_if_newer unemap/bin/rs6000-aix/unemap bin/rs6000-aix/unemap
  update_if_newer unemap/bin/rs6000-aix/unemap-3d bin/rs6000-aix/unemap-3d

  date
)


echo
echo "--------------------------------------------------------"
echo "Process examples"
echo "--------------------------------------------------------"

# Generate example changes list
# PJB: not sure this works anymore
${CMISS_ROOT}/bin/cmisschanges.pl -tree examples -out html > ${CMISS_EXAMPLES}/example_cvs_changes.log



# This needs to run on a machine that can has an smtp client
echo "Sending mail ..."
mailout=/tmp/mail-$$.out
mailerr=/tmp/mail-$$.err
{
  ${CMISS_EXAMPLES}/common/cmiss_type1/mail_tree.pl -p ".*" -f ${CMISS_EXAMPLES}/notify_all_mail.html
  date
} > $mailout 2> $mailerr &

echo "Updating web ..."
${CMISS_EXAMPLES}/common/cmiss_type1/html_tree.pl -f ${CMISS_EXAMPLES}/list_all.html
${CMISS_EXAMPLES}/common/cmiss_type1/generate_index_thumbs.pl
date
#
# Update the executables if the a examples were successful.
# I would test on all the cmgui examples but that would require changes to mail_tree.pl
#
${CMISS_EXAMPLES}/common/cmiss_type1/mail_tree.pl -e a -f ${CMISS_EXAMPLES}/notify_a_mail.html 2>&1 ;
if grep Success ${CMISS_EXAMPLES}/notify_a_mail.html > /dev/null; then
    echo "Updating cmgui executables..."
	 (cd ${CMISS_ROOT}/cmgui/bin/i686-linux && cp --update --backup --suffix=.previous --remove-destination cmgui cmgui-debug cmgui-motif cmgui-motif-debug cmgui-wx cmgui-wx-debug cmgui-gtk ${CMISS_ROOT}/bin/i686-linux/)
	 (cd ${CMISS_ROOT}/cmgui/bin/rs6000-aix && cp --update  --backup --suffix=.previous --remove-destination cmgui cmgui-debug cmgui-motif cmgui-motif-debug cmgui64 ${CMISS_ROOT}/bin/rs6000-aix/)
	 (cd ${CMISS_ROOT}/cmgui/bin/mips-irix && cp --update  --backup --suffix=.previous --remove-destination cmgui cmgui-debug cmgui-motif cmgui-motif-debug cmgui64 ${CMISS_ROOT}/bin/mips-irix/)
	 (cd ${CMISS_ROOT}/cmgui/bin/x86_64-linux && cp --update  --backup --suffix=.previous --remove-destination cmgui cmgui-debug cmgui-motif cmgui-motif-debug cmgui-wx cmgui-wx-debug ${CMISS_ROOT}/bin/x86_64-linux/)
	 #Update the command help to match this version
	 perl ${CMISS_ROOT}/cmgui/source/utilities/commands2html/commands2html.pl cmgui > ${CMISS_ROOT}/www/help/cmgui/gfx/index.html
else
    echo "Skipping the executable update as examples failed."
fi
date

wait # for mail tree
echo STDOUT from mail_tree
cat $mailout
echo END OF STDOUT from mail_tree
rm $mailout
if [ -s $mailerr ]; then
    {
      echo STDERR from mail_tree
      cat $mailerr
      echo END OF STDERR from mail_tree
    } >&2
fi
rm $mailerr

echo "-----------------------------------------------------"

# End of jobs
#
echo cronjob started "$startdate" ended `date`.
echo

# record stdout and stderr in log files.
} >"$daylog" 2>"$dayerr"

if [ -s $dayerr ]; then
    # stderr was recorded in the log file;
    # send it to the old stderr to be mailed.
    cat $dayerr 1>&2
else
    # successful; clean out empty file
    rm $dayerr
fi
