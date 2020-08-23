#!/bin/tcsh
set data_name1=${1:t}
set data_name=${data_name1:r}
sed -e "s%cmiss_www_url_goes_here%${CMISS_WWW_URL}%g" -e "s%rel_dir_name_goes_here%$1%" -e "/name_file_goes_here/r ${CMISS_DATA}/$1/name.txt" -e "/name_file_goes_here/r ${CMISS_DATA}/$1/node_name.txt" -e "/name_file_goes_here/d" -e "s%daname%$data_name%" ${CM_VIEW_DATA}/html_name_template >> ${CM_VIEW_DATA}/html_names
unset data_name
