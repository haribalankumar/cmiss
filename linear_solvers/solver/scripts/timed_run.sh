#!/bin/sh

OMP_NUM_THREADS=1
export OMP_NUM_THREADS

XLSMPOPTS="delays=100 : yields=0 : spins=0"
export XLSMPOPTS

function timed_run {

  # Make tmpdirectory
  tmpfile=tmp.$$
  tmpdir=tmp
  if [ ! -d $tmpdir ]
  then
    mkdir $tmpdir
  fi

  if [ $nrep = 0 ]
  then
    nrep=15
  fi

  # Run once so we can check the solution is right
  echo
  echo "***************************************************************"
  echo
  echo "Program = $exe, Solver = $solver, NCPU = $ncpu"
  echo "Problem = $problem"
  echo "Nrep = $nrep, Tmpfile = $tmpfile, Tmpdir = $tmpdir"
  echo "Maxtime = $maxtime, Maxit = $maxit"
  echo "Other options = $other_options"
  echo
  $exe -q $problem -solver $solver -ncpu $ncpu -niter $maxit -maxtime $maxtime $other_options
  exit_code=$?
  echo

  if [ $exit_code -ne 0 ]
  then
    return 1
  fi

  # Run a few reps, so we can check the fastest run
  i=0
  while [ $i -lt $nrep ]
  do
    $exe -q $problem -solver $solver -ncpu $ncpu -niter $maxit -maxtime $maxtime $other_options 2> /dev/null
    i=` echo "$i + 1" | bc `
  done > $tmpfile

  # Get the fastest runs
  for flag in "Factorisation" "Solution time" "Total time"
  do
    min=`grep "$flag" $tmpfile | w2w | cut -f2 -d: | cut -f3 -d' ' | sort -ni | head -1`
    echo " $flag	:	$min"
  done
  echo

  mv $tmpfile $tmpdir

  return 0
}

