#!/bin/sh
#
# CMISS daily crontab commands
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
name=`hostname -s`
daylog=$logdir/$name.day.out
dayerr=$logdir/$name.day.err
weeklog=$logdir/$name.week.log
lastweeklog=$logdir/$name.lastweek.log
cmmakelog=$logdir/i686-linux_cm_make_fail.log
 
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

# Number of cpus to flog
threads=2

# Cross compilation against older glibc
MAKE=i386-glibc23-linux-cross-make

# Record the system load.
# Use two displays to try to record info over a longer period between displays.
# The first display is not recorded.

# top from procps 3.1.6 wants TERM set even in batch mode.
# This is fixed by procps 3.2.4.

TERM=dumb top -b -n2 -d10 | perl -n -we 'BEGIN{ $/ = "" }; $. > 2 and print'

echo "-----------------------------------------------------"
echo

#
# CMISS compilation jobs
#

echo "Ensuring perl interpreter is up to date..."

$MAKE -j "$threads" -k -C "$CMISS_ROOT/perl_interpreter" all 2>&1 || {
    echo Error making perl_interpreter
    $MAKE -k -C "$CMISS_ROOT/perl_interpreter" all
} 1>&2

echo "Ensuring linear solvers are up to date..." 

$MAKE -j "$threads" -k -C "$CMISS_ROOT/linear_solvers" all 2>&1 || {
    echo Error making linear_solvers
    $MAKE -k -C "$CMISS_ROOT/linear_solvers" all
} 1>&2

echo
echo "-----------------------------------------------------"
echo "Making the CMISS backend..."
echo "-----------------------------------------------------"

: > "$cmmakelog"

$MAKE -j "$threads" -k -C "$CMISS_ROOT/cm" all 2>&1 || {
    date
    $MAKE -k -C "$CMISS_ROOT/cm" all
} >> "$cmmakelog" 2>&1

echo
echo "-----------------------------------------------------"
echo "Compiling cmgui executables"
echo "-----------------------------------------------------"

# Ensure cmgui executables are up to date before running tests.
# The make update job does this too.
$MAKE -j "$threads" -k -C "$CMISS_ROOT/cmgui" all 2>&1 

echo
echo "-----------------------------------------------------"
echo "Compiling zinc and zinc-npruntime"
echo "-----------------------------------------------------"
$MAKE -k -C "$CMISS_ROOT/zinc" all 2>&1 
$MAKE -k -C "$CMISS_ROOT/zinc-npruntime" all 2>&1 

echo
date
echo
echo "--------------------------------------------------------"
echo "Testing examples"
echo "--------------------------------------------------------"

testdir=${TMPDIR:-/tmp}/cmiss_testing
mkdir "$testdir"
(
    if cd "$testdir"; then
        # Test optimized cm on Sunday; debug cm other days, all cmgui every day
        #SAB Added a new flag that allows us to test all versions that were broken
        # in addition to the ones that are selected in the -v flag
	if [ "$weekday" = Sunday ]; then
	    "$CMISS_EXAMPLES/common/cmiss_type1/test_tree.pl" \
		-v '.*/(cmgui.*|cm)' -z "$threads" --verbose --check-all-broken-versions -M time=18000
	else
	    "$CMISS_EXAMPLES/common/cmiss_type1/test_tree.pl" \
		-v '.*/(cmgui.*|cm-debug)' -z "$threads" --verbose --check-all-broken-versions -M time=18000
	fi
    fi
)
rm -rf "$testdir"

echo "-----------------------------------------------------"

#
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
