#! /bin/sh
#
#  shell script for invoking the subdirectory shell scripts that run ftnchek
#  and parse its output.
#

export CM_PARSE_FTNCHK=${CM_PARSE_FTNCHK:-`cd ${0%/*}; pwd`/fortran_check}

$CM_PARSE_FTNCHK/fortran_check_cm.pl > $CM_PARSE_FTNCHK/checkcmiss.out \
    || exit 1

echo "Parsing ftnchek errors in CMISS"
$CM_PARSE_FTNCHK/parse_ftnchek_errors.sh

echo "Parsing common block errors"
${0%/*}/FCheck/parse-fortran-check.sh
