#!/bin/sh
#
# Shell file for finding errors in CMISS fortran dimensioning of arrays
#
# Usage: check_dimension.sh fe*.f
# Created: Glen Harris
# Updates: Chris Bradley 29-Feb-1996
#          Martyn Nash 28 May 96: added error count
#

script_dir=${0%/*}
check_dimension_dir=${CM_CHECK_DIMENSION:-`cd "$script_dir" && pwd`}

if [ ! -d "$check_dimension_dir" ]; then
    echo $0: no check_dimension directory found 1>&2
    exit 1
fi

{
    if [ ${1+set} ]; then
	for filename; do
	    echo $filename
	done
    else
	module_list=$CMISS_ROOT/cmiss_utils/www/MasterLists/module-list
	cat $module_list | sed s%^.*/cm/source/%%
# 	find . -name '*.f' | sed s%^.*/%%
    fi
} | {
  cd $CMISS_ROOT/cm/source && xargs awk -f $check_dimension_dir/step1.awk
} > $check_dimension_dir/all_variables

awk -f $check_dimension_dir/step2.awk $check_dimension_dir/variables.check $check_dimension_dir/all_variables > $check_dimension_dir/error_list

rm $check_dimension_dir/all_variables

# Use cat instead of cp so as not to preserve read-only permissions.
cat $check_dimension_dir/html_dimension_template > $check_dimension_dir/dimension_error.html

awk -f $check_dimension_dir/step3.awk $check_dimension_dir/error_list >> $check_dimension_dir/dimension_error.html

$CMISS_ROOT/html_utils/addfooter.sh programmer $check_dimension_dir/dimension_error.html

mv $check_dimension_dir/dimension_error.html $CMISS_ROOT/www/help/errors/dimension_error.html
