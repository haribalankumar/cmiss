#! /bin/sh
#
# The is a shell script to fix change #/name to #s/name in ip files
#
# Usage:
#   ipcellipmate.sh
#
# Created:
#   Martin Buist, 13 January 2000
#
# Updates:
#

find . -name '*.ipmate' -print > /tmp/ipcellipmate.files
sed s/.ipmate/.ip/ /tmp/ipcellipmate.files > /tmp/ipcellipmate.files_no_ext

for file in `cat /tmp/ipcellipmate.files_no_ext`
do
	echo "Converting file "$file"mate"
	if [ ! -s $file"mate_cell_old" ]; then
        	cp $file"mate" $file"mate_cell_old"

		nawk -f ${CMISS_ROOT}/cmiss_utils/filefixes/19990918/splitipmate.awk $file"mate_cell_old" > $file"mate_tmp"
		nawk -f ${CMISS_ROOT}/cmiss_utils/filefixes/19990918/splitipcell.awk $file"mate_cell_old" > $file"cell_tmp"

		mv $file"mate_tmp" $file"mate"
		if [ ! -s $file"cell" ]; then 
			mv $file"cell_tmp" $file"cell"
		fi
		if [ -s $file"cell_tmp" ]; then 
			rm $file"cell_tmp"
		fi
	fi
done

rm /tmp/ipcellipmate.files
rm /tmp/ipcellipmate.files_no_ext
