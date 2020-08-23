#!/bin/sh
#
# Shell file for creating the HTML list of data which
# goes in the main data list document.
#
# Usage: 
# Created: 
#   Carey Stevens 4-Feb-1997
# Updates: 
#
rm -f ${CM_VIEW_DATA}/html_names
for reldirname
do
  ${CM_VIEW_DATA}/create_html_names_sub.sh $reldirname
done
