#! /bin/sh
#
# Make Lists.sh
#
#   Make the two master lists of modules (e.g. fe10.f) and common files
#   (e.g. b01.cmn) for use with the other programs
#


# -- Script Constants ---------------------------------------------------------
#

MODULE_PATTERN=*.f
# MODULE_EXCLUDE_PATTERN='fe(10\.f|(09|10)_[^sl])'
# MODULE_EXCLUDE_PATTERN_BR='fe(09|10)_linux.f'
if [ "$CMISS_LOCALE" == OXFORD ]; then
  ARCHIVE_PATTERN=cmiss_archive*.f
else
  ARCHIVE_PATTERN=../archive/cmiss_archive*.f
fi
C_MODULE_PATTERN=[fc]*.c
COMMON_PATTERN=*.cmn
LIST=/bin/ls
LIST_OPTIONS=-1



# -- Resource Locations (and defaults) ----------------------------------------
#

if [ "$CMISS_ROOT" = "" ]; then
  echo \$CMISS_ROOT variable unbound.
  exit 1
else
  INPUT_DIR=$CMISS_ROOT/cm/source
  OUTPUT_DIR=${0%/*}
  COMMON_LIST=$OUTPUT_DIR/common-list
  f_module_list=$OUTPUT_DIR/module-list
  BROWSER_MODULE_LIST=$OUTPUT_DIR/module-list-browser #for routine browser index
  C_MODULE_LIST=$OUTPUT_DIR/c-module-list
  ARCHIVE_LIST=$OUTPUT_DIR/archive-list
fi

# -- Read/Write/Execute Tests -------------------------------------------------
#

if [ ! -d $INPUT_DIR ]; then
  echo $0: input directory $INPUT_DIR does not exist \(aborting\).
  exit 1
fi

if [ ! -r $INPUT_DIR ]; then
  echo $0: input directory $INPUT_DIR is not readable \(aborting\).
  exit 1
fi

if [ ! -d $OUTPUT_DIR ]; then
  echo $0: output directory $OUTPUT_DIR does not exist \(aborting\).
  exit 1
fi

if [ ! -w $f_module_list ]; then
  if [ ! -w $OUTPUT_DIR ]; then
    echo ERROR [$0]: Directory $OUTPUT_DIR is not writable
    exit 1
  else
    if [ -f $f_module_list ]; then
      echo ERROR [$0]: File $f_module_list is not writable
      exit 1
    fi
  fi
fi

if [ ! -w $COMMON_LIST ]; then
  if [ ! -w $OUTPUT_DIR ]; then
    echo ERROR [$0]: Directory $OUTPUT_DIR is not writiable
    exit 1
  else
    if [ -f $COMMON_LIST ]; then
      echo ERROR [$0]: File $COMMON_LIST is not writable
      exit 1
    fi
  fi
fi

if [ ! -x $LIST ]; then
  echo $0: cannot run $LIST to list the files \(aborting\).
  exit 1
fi



# -- Generate Lists -----------------------------------------------------------
#

find $INPUT_DIR -name "$MODULE_PATTERN" > $f_module_list
cp $f_module_list $BROWSER_MODULE_LIST
# $LIST $LIST_OPTIONS $INPUT_DIR/$MODULE_PATTERN1   > ${MODULE_LIST}_tmp
# egrep -v "$MODULE_EXCLUDE_PATTERN" ${MODULE_LIST}_tmp > $MODULE_LIST
# rm ${MODULE_LIST}_tmp
# egrep -v "$MODULE_EXCLUDE_PATTERN_BR" $MODULE_LIST > $BROWSER_MODULE_LIST
$LIST $LIST_OPTIONS $INPUT_DIR/$ARCHIVE_PATTERN   > $ARCHIVE_LIST
$LIST $LIST_OPTIONS $INPUT_DIR/$C_MODULE_PATTERN > $C_MODULE_LIST
$LIST $LIST_OPTIONS $INPUT_DIR/$COMMON_PATTERN   > $COMMON_LIST



# -- Return -------------------------------------------------------------------
#

exit 0
