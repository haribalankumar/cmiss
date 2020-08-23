#! /bin/csh -f
#
# This file sets up all cmiss programmer commands and aliases. 
#
# Note: This script relies on (and hence checks that) the CMISS_ROOT 
# environment variable is set.
#
# Usage: 
#   source cmiss_programmer.sh
# Created: 
#   Chris Bradley, 3 June 1996
# Updates: 
#
#

if(${CMISS_ROOT} == "") then
	echo "Must have CMISS_ROOT enviroment variable set"
else

#
#       Enviroment variables
#


	setenv CM_ROOT ${CMISS_ROOT}/cm
	setenv CMGUI_ROOT ${CMISS_ROOT}/cmgui
	setenv CMISS_DATA ${CMISS_ROOT}/data

# DDOT in the sgi 7.2.1 blas library seems to perform INVALID
# operations that do not affect the result.  To avoid these being trapped in the
# debug version, libfpe is used to quietly set the results of fpe's involving
# NaN to NaN.
# Removed as esu35 now uses 7.2.0 blas
#	setenv TRAP_FPE "INVALID=APPROPRIATE"

#
#       Other script files
#

# KAT out of date
# 	source ${CM_ROOT}/scripts/runaliases.sh
# 	source ${CM_ROOT}/scripts/sourcealiases.sh
# 	source ${CM_ROOT}/scripts/otheraliases.sh
# 	source ${CMGUI_ROOT}/scripts/runaliases.sh
#        source ${CMGUI_ROOT}/scripts/sourcealiases.sh
# 	source ${CMGUI_ROOT}/scripts/otheraliases.sh

#	alias cmiss ${CMISS_ROOT}/cmiss/cmiss
	alias cmiss '${CMISS_ROOT}/bin/cmgui -cm'
	alias mycmiss '${CMISS_ROOT}/bin/cmgui -mycm'
endif
