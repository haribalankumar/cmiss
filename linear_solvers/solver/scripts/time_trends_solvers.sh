#!/bin/sh

. ./timed_run.sh
. ./time_trend_nx.sh

exe=../bin/soltest_64
nrep=9
nlow=10
nhigh=100
nstep=10
maxtime=35

for trans in -trans -notrans
do

  for solver in jacobi sor ilu0 ilu1
  do
    sys=-fv3d
    nlow=10
    nhigh=100
    nstep=10
    echo "timed_trend > ${solver}${sys}${trans}.out"
    timed_trend > ${solver}${sys}${trans}.out

    sys=-fv2d
    nlow=100
    nhigh=1000
    nstep=100
    echo "timed_trend > ${solver}${sys}${trans}.out"
    timed_trend > ${solver}${sys}${trans}.out
  done

  for cg in bicgstab gmres
  do
    for precon in none jacobi sor ilu0 ilu1
    do
      solver="$cg -precon $precon"
      outfile=${cg}_${precon}

      sys=-fv3d
      nlow=10
      nhigh=100
      nstep=10
      echo "timed_trend > ${solver}${sys}${trans}.out"
      timed_trend > ${outfile}${sys}${trans}.out

      sys=-fv2d
      nlow=100
      nhigh=3000
      nstep=100
      echo "timed_trend > ${solver}${sys}${trans}.out"
      timed_trend > ${outfile}${sys}${trans}.out
    done
  done
done

for trend in -trans
do

  solver=superlu

  for slu_co in 0 1 2 3
  do
    sys=-fv3d
    nlow=10
    nhigh=100
    nstep=10
    echo "timed_trend > ${solver}${slu_co}${sys}${trans}.out"
    timed_trend > ${solver}${slu_co}${sys}${trans}.out

    sys=-fv2d
    nlow=100
    nhigh=3000
    nstep=100
    echo "timed_trend > ${solver}${slu_co}${sys}${trans}.out"
    timed_trend > ${solver}${slu_co}${sys}${trans}.out
  done

  for solver in umfpack harwell
  do
    sys=-fv3d
    nlow=10
    nhigh=100
    nstep=10
    other_options="-umf_narr 100 -har_narr 50,50"
    echo "timed_trend > ${solver}${sys}${trans}.out"
    timed_trend > ${solver}${sys}${trans}.out

    sys=-fv2d
    nlow=100
    nhigh=3000
    nstep=100
    other_options="-umf_narr 100 -har_narr 50,50"
    echo "timed_trend > ${solver}${sys}${trans}.out"
    timed_trend > ${solver}${sys}${trans}.out
  done
done
