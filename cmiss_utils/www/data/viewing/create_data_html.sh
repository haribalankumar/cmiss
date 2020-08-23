#!/bin/sh
#
# Creates the HTML page for each directory
#
# Usage: 
#   create_data_html.sh <list of directories>
# Created: 
#   Carey Stevens 4-Feb-1997
# Updates: 
#
for dirname
do
  ${CM_VIEW_DATA}/create_data_html_sub.sh $dirname
done
