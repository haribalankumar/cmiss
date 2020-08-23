#! /bin/csh
#
# Shell file to output html levels as separate files.
#
# Usage: 
#   make_html.sh name data_file
# Created: 
#   Glen Harris 25 October 1996
# Updates: 
#   Richard Boyes 4th August 1999
#   Getting rid of html files as there may be a build up of them
#
set name=$1
set data_file=$2

if x$name == x then
    exit 1
endif

#Richard Boyes - getting rid of old html files.
#  echo "Removing the "${name}" html files..."
#  if { cd ${CM_PARSE_FTNCHK} } then
#    /usr/bin/ls -1 | grep "^${name}.*\.html$" | xargs rm -f
#  endif
if -d ${CM_PARSE_FTNCHK}/${name} then
    #Richard Boyes - getting rid of old html files
    echo "Removing the "${name}" html files..."
    # rm -r and find -exec has trouble over nfs
    find ${CM_PARSE_FTNCHK}/${name} -type f | xargs -l rm -f
    rmdir ${CM_PARSE_FTNCHK}/${name}
endif
mkdir ${CM_PARSE_FTNCHK}/${name}

echo "Generating the ${name} html files..."
(cd ${CM_PARSE_FTNCHK}/${name}; ../format_level.pl ${name} ${data_file})

# This foreach statement needs to be replaced with something a little
# trickier as if there are more than 256 files the argument list will
# be too long - Richard Boyes
#  foreach file (${CM_PARSE_FTNCHK}/${name}*.html)
#  	echo "Adding footer to "$file
#  	${CMISS_ROOT}/html_utils/addfooter.sh programmer ${file}
#  end
cd ${CM_PARSE_FTNCHK}/${name} 
echo "Adding footer to ${name}/*.html"
/bin/ls -1 | grep ^"${name}"'.*\.html$' | xargs -l ${CMISS_ROOT}/html_utils/addfooter.sh programmer

#  cp ${CM_PARSE_FTNCHK}/${name}*.html ${CMISS_ROOT}/www/help/errors/${name}
if ! -d ${CMISS_ROOT}/www/help/errors/${name} then
  ln -s ${CM_PARSE_FTNCHK}/${name}/ ${CMISS_ROOT}/www/help/errors/
endif
#find ${CM_PARSE_FTNCHK} -name ${name}\*.html -print | sed -e s%${CM_PARSE_FTNCHK}%"."%g | cpio -puvm ${CMISS_ROOT}/www/help/errors/${name}

#  echo "Finished copying files....."

#Allow the error counting program to get at the base html file
ln -sf ${CMISS_ROOT}/www/help/errors/${name}/${name}.html ${CMISS_ROOT}/www/help/errors/${name}.html

echo "Finished making symbolic link....."

