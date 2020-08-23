#
#  Makefile for CMISS Master File Lists
#
#    This makefile does not compile anything - rather it is here to
#    support the same Utility functions that the other CMISS
#    information management makefiles do.
#


# -- Define Shell -------------------------------------------------------------
#

SHELL=/bin/sh



# -- Files --------------------------------------------------------------------
#

SCRIPTS      = make-lists.sh
DOCS         = README
OUTPUT_FILES = module-list common-list



# -- Default Target -----------------------------------------------------------
#
 
all: 



# -- Utility Targets ----------------------------------------------------------
#

tar:
	-rm -f source.tar.gz 
	tar cfv source.tar Makefile $(DOCS) $(SCRIPTS)
	gzip source.tar

untar: source.tar.gz
	zcat source.tar.gz | tar xfv -

clean: 
	-rm -f Makefile.bak source.tar.gz $(OUTPUT_FILES)
