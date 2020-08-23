#!/bin/sh
#
# CMISS daily master crontab commands
# Also spawns jobs on bioeng85, bioeng22
#  (AJC-20080418: bioeng69 has been replaced by bioeng85 and esu20 nolonger exists)
#

# Set up PATH and other environment as for login shell
. $HOME/.profile

if [ -z "$CMISS_ROOT" ]; then
    echo CMISS_ROOT not set 1>&2
    exit 1
fi

# Day of the week
weekday=`date +%A`

# Short hostname
host=`hostname -s`

# All STDOUT is sent to a log file.
# Only STDERR that is expected is sent to this file.
# In this way, if there is an unexpected error, a mail message is sent.

logdir=$CMISS_ROOT/cmiss_utils/cronjobs

name=master
dayout="$logdir/$name.day.out"
dayerr="$logdir/$name.day.err"
weeklog="$logdir/$name.week.log"
lastweeklog="$logdir/$name.lastweek.log"

# Back up week's log
[ "$weekday" = Sunday ] && mv "$weeklog" "$lastweeklog"
# Back up output file in week's log
cat "$dayout" >> "$weeklog"
# Include any stderr from the previous day in the week log.
# Use the first line of the stdout as a record of when stderr was generated.
if [ -s "$dayerr" ]; then
    {
      echo STDERR from `head -1 "$dayout"`
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

# Number of cpus to flog
threads=1


# Both cm and cmgui require a perl interpreter and linear solvers.

#echo "Ensuring linear solvers are up to date..." 

# Linear solvers do not make on more than one machine using
# esu1:/product/cmiss concurrently as they share the same temporary files, so
# make these on esu8 before spawning other jobs.

#gmake -j "$threads" -k -C "$CMISS_ROOT/linear_solvers" all 2>&1 || {
#    echo Error making linear_solvers
#    gmake -k -C "$CMISS_ROOT/linear_solvers" all
#} 1>&2

# Start bioeng85 jobs.
# There should be no stdout as the script stores output in a log file.
bioeng85err=/tmp/bioeng85-$$.err

ssh bioeng85 exec '$CMISS_ROOT/cmiss_utils/cronjobs/cmiss_cronjobs_daily_bioeng85.sh' \
    >$bioeng85err 2>&1 &

# Start bioeng1031 jobs.
# There should be no stdout as the script stores output in a log file.
bioeng1031err=/tmp/bioeng1031-$$.err

ssh bioeng1031 exec '$CMISS_ROOT/cmiss_utils/cronjobs/cmiss_cronjobs_daily_bioeng1031.sh' \
    >$bioeng1031err 2>&1 &

# Start bioeng22 jobs.
# There should be no stdout as the script stores output in a log file.
#bioeng22err=/tmp/bioeng22-$$.err

#ssh bioeng22 exec '$CMISS_ROOT/cmiss_utils/cronjobs/cmiss_cronjobs_daily_bioeng22.sh' \
#    >$bioeng22err 2>&1 &

#echo "Ensuring perl interpreter is up to date..." 

#gmake -j "$threads" -k -C "$CMISS_ROOT/perl_interpreter" all 2>&1 || {
#    echo Error making perl_interpreter
#    gmake -k -C "$CMISS_ROOT/perl_interpreter" all
#} 1>&2



#
# End of jobs
#
echo cronjob started "$startdate" ended `date`.
echo

# record stdout and stderr in log files.
} >"$dayout" 2>"$dayerr"

if [ -s $dayerr ]; then
    # stderr was recorded in the log file;
    # send it to the old stderr to be mailed.
    cat $dayerr 1>&2
else
    # successful; clean out empty file
    rm $dayerr
fi

wait # for jobs on other machines
for machine in bioeng85 bioeng1031; do
    eval "errfile=\$${machine}err"
    if [ -s $errfile ]; then
	{
	  echo STDERR from $machine jobs
	  cat $errfile
	  echo END OF STDERR from $machine
	} >&2
    fi
    rm $errfile
done
