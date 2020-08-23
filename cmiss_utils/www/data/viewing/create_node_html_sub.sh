#!/bin/tcsh
ls -1 $1/*/generated/*.html | sed -e "s%\.dir%%g" > file_names 
sed -e "s%${CMISS_DATA}/%%g" -e "s%/generated/data_page.html%%g" ${CM_VIEW_DATA}/file_names > rel_file_names
./create_html_names.sh `cat rel_file_names`
set data_name=data_${1:t}
sed s%daname%$data_name% create_node_html.sed > create_node_html_custom_1.sed
sed s%dirname%$1% create_node_html_custom_1.sed > create_node_html_custom_2.sed
unset data_name
rm -f $1/generated/data_page.html
if( { test -f $1/description.html } ) then
#Change the path inside the description files
sed -e 's%\${CMISS_DATA_URL}%'${CMISS_WWW_URL}help/data/data_files%g -e 's%\${CMISS_WWW_URL}%'${CMISS_WWW_URL}%g $1/description.html > $1/description_temp.html
else
echo "No description for this data. <BOLD>Invent a description and add it!</BOLD>" > $1/description_temp.html
endif
sed -f create_node_html_custom_2.sed data_node_template.html > $1/generated/data_page.html
rm $1/description_temp.html
${CMISS_ROOT}/html_utils/addfooter.sh user $1/generated/data_page.html
rm create_node_html_custom_1.sed
rm create_node_html_custom_2.sed
