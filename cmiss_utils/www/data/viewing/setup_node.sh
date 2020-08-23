#!/bin/sh
#
# Sets up the structure for each directory.
#
# Usage: 
#   setup_node.sh <list of directories>
# Created: 
#   Carey Stevens
# Updates: 
#
for dirname
do
  ${CM_VIEW_DATA}/setup_node_sub.sh $dirname
done
