#! /bin/sh
#
# Shell file to strip out calls with varying numbers of arguments.
#
# Usage: 
#   parse_ftnchek_errors.sh
# Created: 
#   Glen Harris 7 May 1996
# Updates: 
#   Martyn Nash 28 May 96: Added error counts for html pages
#   Glen Harris 29 May 96: Removed error counts for html pages - done 
#     at the referring page level (ie in ${CMISS_ROOT}/www/help/errors)
#   Glen Harris 11 June 96: Now works on all awk files in this directory
#   Glen Harris 24 October 96: Passing in a 'omit' list for each awk 
#     file.  The first file to each awk script is a list of lines
#     not to output.
#   Glen Harris 25 October 96: Added an option to just run one script.
#

export CM_PARSE_FTNCHK=${CM_PARSE_FTNCHK:-`cd ${0%/*}; pwd`}

if [ ! -d "$CM_PARSE_FTNCHK" ]; then
    echo $0: no parse ftnchk directory found 1>&2
    exit 1
fi

if [ -z "$1" ]; then
    # run perl parser
    ${CM_PARSE_FTNCHK}/parse_ftnchek_errors.pl $CM_PARSE_FTNCHK/checkcmiss.out > ${CM_PARSE_FTNCHK}/unprocessed.out
    if [ -s ${CM_PARSE_FTNCHK}/unprocessed.out ]; then
	{
	  echo "$0: unprocessed ftnchek output":
	  cat ${CM_PARSE_FTNCHK}/unprocessed.out
	} >&2
    fi
    data_files="${CM_PARSE_FTNCHK}/*.dat"
else
	data_files=$1
fi
echo $data_files

tempfile=parse_ftncheck.tmp$$
for full_name in ${data_files}; do
    tail=${full_name%.*}
    name=${tail##*/}
    # Check for the omit file - if it doesn't exist, make it
    if [ ! -e ${CM_PARSE_FTNCHK}/${name}.omit ]; then
	touch ${CM_PARSE_FTNCHK}/${name}.omit
    fi
    $CM_PARSE_FTNCHK/update_date.pl $name
    echo -n Processing omits for $name...
    awk -f	${CM_PARSE_FTNCHK}/omit.awk_util ${CM_PARSE_FTNCHK}/${name}.omit \
		${CM_PARSE_FTNCHK}/$name.date > $tempfile
    echo done
    mv $tempfile ${CM_PARSE_FTNCHK}/$name.trim
    ${CM_PARSE_FTNCHK}/make_html.sh ${name} ${CM_PARSE_FTNCHK}/$name.trim
done
