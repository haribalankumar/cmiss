#!/bin/tcsh
if ( ! -d $1/generated ) then
        mkdir $1/generated
endif
ls -1 $1 > file_names
sed s%dirname%${CMISS_WWW_URL}help/data/data_files$1% create_data_html.awk > create_data_html_custom_1.awk
nawk -f create_data_html_custom_1.awk file_names > data_file_names
set data_name=data_${1:t}
sed s%daname%$data_name% create_data_html.sed > create_data_html_custom_1.sed
sed s%dirname%$1% create_data_html_custom_1.sed > create_data_html_custom_2.sed
unset data_name
rm -f $1/generated/data_page.html
if( { test -f $1/description.html } ) then
#Change the path inside the description files
sed -e 's%\${CMISS_DATA_URL}%'${CMISS_WWW_URL}help/data/data_files%g -e 's%\${CMISS_WWW_URL}%'${CMISS_WWW_URL}%g $1/description.html > $1/description_temp.html
else
echo "No description for this data. <BOLD>Invent a description and add it!</BOLD>" > $1/description_temp.html
endif
sed -f create_data_html_custom_2.sed data_html_template.html > $1/generated/data_page.html
rm $1/description_temp.html
${CMISS_ROOT}/html_utils/addfooter.sh user $1/generated/data_page.html
rm create_data_html_custom_1.awk
rm create_data_html_custom_1.sed
rm create_data_html_custom_2.sed
