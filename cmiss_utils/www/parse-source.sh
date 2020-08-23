#! /bin/sh
#
#  shell script for invoking the subdirectory shell scripts that parse and
#  index cm source.
#
#  NB. All the scripts are run *relative* to the location of this script. 
#


# -- File and Directory Locations ---------------------------------------------
#

MASTERLISTS_SCRIPT=MasterLists/make-lists.sh
COMMENTSTRIPPER_SCRIPT=CommentParser/comment-stripper.sh
ROUTINEINDEXER_SCRIPT=RoutineBrowser/make-index.sh
#  KAT 2Dec99:  This doesn't seem appropriate anymore
#  HTMLIZER_SCRIPT=HTMLizeUnusedVars/htmlize-unused-vars.sh
CALLTREE_SCRIPT=CallTree/call-tree.sh


# -- Batch Invocation ---------------------------------------
#

echo "Making Master Lists ..."
$MASTERLISTS_SCRIPT

echo "Stripping Comments ..."
# Sends stuff to stderr that has not been handled.  Redirect for cronjobs.
$COMMENTSTRIPPER_SCRIPT 2>&1 || \
    echo $COMMENTSTRIPPER_SCRIPT failed to strip comments 1>&2

echo "Making Routine Index ..."
$ROUTINEINDEXER_SCRIPT

echo "Parsing Call Tree ..."
# Sends stuff to stderr that has not been handled.  Redirect for cronjobs.
$CALLTREE_SCRIPT 2>&1 || \
    echo $CALL_TREE_SCRIPT failed to generate tree 1>&2
