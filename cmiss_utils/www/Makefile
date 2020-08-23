#
#  Makefile for CMISS information managment programs
#
#


# -- Define Shell -------------------------------------------------------------
#

SHELL=/bin/sh



# -- Files --------------------------------------------------------------------
#

DIRECTORIES = CallTree CommentParser FCheck HTMLizeUnusedVars MasterLists CGIDecode ExampleViewer RoutineBrowser RegEx Utilities
TARFILES = $(DIRECTORIES:=/source.tar.gz)
MAKEFILES = $(DIRECTORIES:=/Makefile)
SCRIPTS = run-batch-jobs.sh



# -- Default Target -----------------------------------------------------------
#

all: 
	for dir in $(DIRECTORIES); do \
	  ( cd $$dir ; make ) \
	done



# -- Utility Targets ----------------------------------------------------------
#

tar:
	for dir in $(DIRECTORIES); do \
          ( cd $$dir ; make tar ) \
        done
	-rm -f source.tar.gz 
	tar cfv source.tar $(TARFILES) $(MAKEFILES) Makefile README $(SCRIPTS)
	gzip source.tar

untar: source.tar.gz
	zcat source.tar.gz | tar xfv -
	for dir in $(DIRECTORIES); do \
          ( cd $$dir ; make untar ) \
        done

clean: 
	-rm -f Makefile.bak source.tar.gz
	for dir in $(DIRECTORIES); do \
          ( cd $$dir ; make clean ) \
        done
