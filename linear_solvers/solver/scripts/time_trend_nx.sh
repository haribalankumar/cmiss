#!/bin/sh

function timed_trend {

  maxtime=30
  maxit=20000

	if [ "$exe" = "" ] ; then
  	exe=../bin/soltest_32
	fi
	if [ "$solver" = "" ] ; then
	  solver=jacobi
	fi
	if [ "$nlow" = "" ] ; then
	  nlow=10
	fi
	if [ "$nhigh" = "" ] ; then
	  nhigh=100
	fi
	if [ "$nstep" = "" ] ; then
	  nstep=10
	fi
	if [ "$sys" = "" ] ; then
	  sys="-fv3d"
	fi
	if [ "$nrep" = "" ] ; then
	  nrep=15
	fi
	if [ "$ncpu" = "" ] ; then
	  ncpu=1
	fi
	if [ "$slu_co" = "" ] ; then
	  slu_co=0
	fi

	n=$nlow
	not_too_long=0
	while [ $not_too_long -eq 0 ] && [ $n -le $nhigh ]
	do
	  problem="$sys $n,$n,$n $trans -slu_co $slu_co"
	  timed_run
	  not_too_long=$?

	  let n+=$nstep
	done
}

