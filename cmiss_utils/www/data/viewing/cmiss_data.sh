#!/bin/sh 
#
# Shell file for creating html views of the data files
#
# Usage: 
#   cmiss_data.sh <data_subset_dir>
# Created: 
#   Carey Stevens 4-Feb-1997
# Updates: 
#

export CM_VIEW_DATA=${CM_VIEW_DATA:-${0%/*}}

if [ ! -d "$CM_VIEW_DATA" ]; then
    echo $0: no view data directory found 1>&2
    exit 1
fi

export CMISS_DATA=${CMISS_DATA:-${CMISS_ROOT:-/product/cmiss}/data}

if [ ! -d "$CMISS_DATA" ]; then
    echo $0: no cmiss data directory found 1>&2
    exit 1
fi

echo "   Updating CMISS data html files..."

#Start off in the correct directory
cd ${CM_VIEW_DATA}

#Find all data directories that contain the file name.txt
find ${CMISS_DATA}/$1 -follow -name name.txt -print > ${CM_VIEW_DATA}/dir_name_list.txt 

#Remove file part from the list - just get directory name
sed -f ${CM_VIEW_DATA}/remove_name.sed ${CM_VIEW_DATA}/dir_name_list.txt > ${CM_VIEW_DATA}/dir_list.txt

#For each directory, create the html page for the data
${CM_VIEW_DATA}/create_data_html.sh `cat ${CM_VIEW_DATA}/dir_list.txt`

#Find all data directories that contain the file node_name.txt
find ${CMISS_DATA}/$1 -follow -name node_name.txt -print > ${CM_VIEW_DATA}/dir_name_list.txt 

#Remove file part from the list - just get directory name
sed -f ${CM_VIEW_DATA}/remove_node_name.sed ${CM_VIEW_DATA}/dir_name_list.txt > ${CM_VIEW_DATA}/dir_list.txt

#For each directory, setup the data_page.html
${CM_VIEW_DATA}/setup_node.sh `cat ${CM_VIEW_DATA}/dir_list.txt`

#For each directory, create the html page for the data
${CM_VIEW_DATA}/create_node_html.sh `cat ${CM_VIEW_DATA}/dir_list.txt`
