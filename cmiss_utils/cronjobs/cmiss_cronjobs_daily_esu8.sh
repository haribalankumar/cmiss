#!/bin/sh
#
# CMISS daily crontab commands for esu8
# Also spawns jobs on bioeng85, bioeng22
#  (AJC-20080418: bioeng69 has been replaced by bioeng85 and esu20 nolonger exists)
#

# Set up PATH and other environment as for login shell
. $HOME/.profile

# Day of the week
weekday=`date +%A`
# Short hostname
host=`hostname -s`

# All STDOUT is sent to a log file.
# Only STDERR that is expected is sent to this file.
# In this way, if there is an unexpected error, a mail message is sent.

logdir=$CMISS_ROOT/cmiss_utils/cronjobs

dayout="$logdir/$host.day.out"
dayerr="$logdir/$host.day.err"
weeklog="$logdir/$host.week.log"
lastweeklog="$logdir/$host.lastweek.log"
# cmmakelog=$logdir/mips-irix_cm_make_fail.log

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
threads=2

# Record the system load and processes using the most cpu
# and virtual and physical memory.

# CPU: need two displays to record info between displays.
# The first display is not recorded.
top -b -d2 -s10 | perl -n -e '$empties == 2 ? print : m/^\n/ && $empties++'
# virtual and physical memory: the header is not rerecorded.
top -b -I -osize 10 | tail +6
top -b -I -ores 10 | tail +6

# Both cm and cmgui require a perl interpreter and linear solvers.

echo "Ensuring linear solvers are up to date..." 

# Linear solvers do not make on more than one machine using
# esu1:/product/cmiss concurrently as they share the same temporary files, so
# make these on esu8 before spawning other jobs.

gmake -j "$threads" -k -C "$CMISS_ROOT/linear_solvers" all 2>&1 || {
    echo Error making linear_solvers
    gmake -k -C "$CMISS_ROOT/linear_solvers" all
} 1>&2

echo "Ensuring perl interpreter is up to date..." 

gmake -j "$threads" -k -C "$CMISS_ROOT/perl_interpreter" all 2>&1 || {
    echo Error making perl_interpreter
    gmake -k -C "$CMISS_ROOT/perl_interpreter" all
} 1>&2



# CMGUI compilation jobs
#
# Note cmgui executables need to be built from up to date source before running
# the tests. We cannot do an svn update on esu8 so it needs to be done by 
# another cronjob before compilation occurs.

# !!! There could be some issues if the same compilation is happening due to
# other cronjobs. Currently the only other cron job which should be building
# cmgui is cmiss_cronjobs_daily_final.sh, which runs at the end of the 
# overnight build/testing and ensures cmgui is up to date everywhere.

echo
echo "-----------------------------------------------------"
echo "Compiling cmgui executables"
echo "-----------------------------------------------------"

gmake -j "$threads" -k -C "$CMISS_ROOT/cmgui" all 2>&1 || {
    echo Error making cmgui
    gmake -k -C "$CMISS_ROOT/cmgui" all
} 1>&2

# Now wait for successful example testing before updating executables 
# executables are updated as part of the cmiss_cronjobs_daily_final.sh cronjob


# CM compilation jobs 
#
# PJB: esu20 does not exist anymore so I have added the cm IRIX build
# here

echo
echo "--------------------------------------------------------"
echo "Compiling cm executables"
echo "--------------------------------------------------------"

 gmake -j "$threads" -k -C "$CMISS_ROOT/cm" all 2>&1 || {
   echo Error making cm
   gmake -k -C "$CMISS_ROOT/cm" all
 } 1>&2


# Generate changes list
# PJB: I think this should be done as part of the final job, will try it there
#${CMISS_ROOT}/bin/cmisschanges.pl -tree cmgui/source -out html -since week > ${CMISS_ROOT}/cmgui/utilities/cmgui_cvs_changes.log


echo
echo "--------------------------------------------------------"
echo "Testing examples"
echo "--------------------------------------------------------"


testdir=${TMPDIR:-/tmp}/cmiss_testing$$
mkdir "$testdir"
(
    if cd "$testdir"; then

   #SAB Added a new flag that allows us to test all versions that were broken
   # in addition to the ones that are selected in the -v flag
	$CMISS_EXAMPLES/common/cmiss_type1/test_tree.pl -p cmgui -v '.*/cmgui.*' -z "$threads" --verbose --check-all-broken-versions

    fi
)
rm -rf "$testdir"


#
# CMISS data (update lists)
#
echo "Updating all CMISS data lists ... `date`"
$CMISS_ROOT/cmiss_utils/www/data/viewing/cmiss_data.sh

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

