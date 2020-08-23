#!/bin/sh
#
# CMISS cm documentation daily crontab commands.  
# 
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
name=cm-docs
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

#
# CMISS documentation
#
echo "Generating ctags TAG files for cm source"
( cd "$CMISS_ROOT/cm/source" && cmissetags )

# Generate data for cmiss-lookup and routine-browser
( cd "$CMISS_ROOT/cmiss_utils/www" && ./parse-source.sh )



#
# cm error checking
#

# Note that cm should be made before running this as it uses the parsed cm
# source files.
echo "Running ftnchek on CMISS source files ... `date`"
(cd "$CMISS_ROOT/cmiss_utils/www" && ./ftnchek-and-parse.sh)


#
# more cm error checking
#
echo "Checking dimensions errors in CMISS"
(cd "$CMISS_ROOT/cm/source" && $CMISS_ROOT/cmiss_utils/www/check_dimension/check_dimension.sh)

echo "Checking enters and exit calls"
"$CMISS_ROOT/cmiss_utils/www/findenterexit/findenterexit.sh"

echo "Counting errors in error pages"
"$CMISS_ROOT/cmiss_utils/www/count_errors/count_errors.sh"


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
