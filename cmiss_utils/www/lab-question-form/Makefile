#
# Makefile for the lab-question-form program
#

CGI_PROGS	=	lab-question-form.cgi 

WEB_PATH	=	/hpc/cmissweb/

CC_OPTS		=	 \
			-DCMISS_WWW_URL=\"http://www.cmiss.org/cm/wiki/CMISSTutorials\" \
			-DCMISS_CGI_URL=\"http://cmiss.bioeng.auckland.ac.nz/\" \
			-DCMISS_WWW_SERVERPATH=\"$(WEB_PATH)\" \
			-DQUESTION_PATH=\"$(WEB_PATH)/lab_questions\" \
			-I/hpc/cmiss/cmiss_utils/www/CGIDecode

#CC_OPTS		=	-n32 \
#			-DCMISS_WWW_URL=\"${CMISS_WWW_URL}\" \
#			-DCMISS_CGI_URL=\"${CMISS_WWW_URL}\" \
#			-DCMISS_WWW_SERVERPATH=\"${CMISS_WWW_SERVERPATH}\" \
#			-DQUESTION_PATH=\"${CMISS_WWW_SERVERPATH}/lab_questions\" \
#			-I/usr/local/include
OBJECTS		=	
LIBRARY_PATH	=	/hpc/cmiss/cmiss_utils/www/CGIDecode/
LIBRARIES	=	-lcgi

SHELL		=	/bin/sh

.SUFFIXES:	.cgi .c

.c.cgi:
	cc $(CC_OPTS) -o $@ $*.c $(OBJECTS) -L$(LIBRARY_PATH) $(LIBRARIES)
	cp $@ $(WEB_PATH)/cgi-bin/
	chmod +s $(WEB_PATH)/cgi-bin/$@
#	cp $@ ${CMISS_ROOT}/www/scripts/.
#	chmod +s ${CMISS_ROOT}/www/scripts/$@

all :	$(CGI_PROGS)

lab-question-form.cgi	:	lab-question-form.c

clean : 
	rm -f $(CGI_PROGS)

clobber : clean
