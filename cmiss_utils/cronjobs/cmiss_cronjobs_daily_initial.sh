#!/bin/sh
#

# CMISS initial daily crontab commands.  Executed just prior to midnight
# before any comilations or tests have run.  Responsible for ensuring
# source and example tree is up to date. Also cleans source weekly

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
name=initial
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

# Do a make clean-all once a week, so that we are building
# everything from scratch at least once a week.
# We don't clean out the uidh directory as these can be tricky
# to generate and often require manual intervention to get working
if [ "$weekday" = Saturday ]; then
    {
      i386-glibc23-linux-cross-make -C "$CMISS_ROOT/cmgui" clean-all-except-uidh 2>&1
    }
fi

# update cmgui source
cd $CMISS_ROOT/cmgui
svn update

# update cmgui examples
cd $CMISS_EXAMPLES/a
svn update

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
