#
#  Makefile for the CGI routine browser
#  and it's related cron indexer
#


# -- Define Shell -------------------------------------------------------------
#

SHELL=/bin/sh



# -- Directory Paths ----------------------------------------------------------
#

COMMENTPARSER_DIR = ../CommentParser
UTILITY_DIR = ../Utilities/

CGIINCLUDE_DIR = ../CGIDecode
CGILIBRARY_DIR = ../CGIDecode

INSTALL_DIR = /usr/httpd/www.esc.auckland.ac.nz/Groups/Bioengineering/CMISS/scripts



# -- Additional Files ---------------------------------------------------------
#

SCRIPTS = make-index.sh
DOCS = README
OUTPUT_FILES = routine-index



# -- Processor Directives -----------------------------------------------------
#

CC        = cc
CFLAGS    = -ansi -pedantic -O0 -g
CDEFINES  = -DCMISS_ROOT=\"${CMISS_ROOT}\"
CINCLUDES = -I$(COMMENTPARSER_DIR) -I$(CGIINCLUDE_DIR) -I$(UTILITY_DIR)



# -- Rules --------------------------------------------------------------------
#

.c.o:
	$(CC) -c $(CFLAGS) $(CDEFINES) $(CINCLUDES) -o $@ $<



# -- Target Files -------------------------------------------------------------
#

UTILITY_SOURCES = $(UTILITY_DIR)htstring.c $(UTILITY_DIR)simple-parser.c
UTILITY_HEADERS = $(UTILITY_SOURCES:.c=.h)
UTILITY_OBJECTS = $(UTILITY_SOURCES:.c=.o)

BROWSER_SOURCES = routine-browser.c
BROWSER_HEADERS = 
BROWSER_OBJECTS = $(BROWSER_SOURCES:.c=.o)

INDEXER_SOURCES = make-index.c
INDEXER_HEADERS = 
INDEXER_OBJECTS = $(INDEXER_SOURCES:.c=.o)

BROWSER = routine-browser.cgi
INDEXER = make-index

LOCAL_SOURCES = $(BROWSER_SOURCES) $(INDEXER_SOURCES)
LOCAL_HEADERS = $(BROWSER_HEADERS) $(INDEXER_HEADERS)
LOCAL_OBJECTS = $(BROWSER_OBJECTS) $(INDEXER_OBJECTS)

SOURCES = $(UTILITY_SOURCES) $(LOCAL_SOURCES) 
HEADERS = $(UTILITY_HEADERS) $(LOCAL_HEADERS) 
OBJECTS = $(UTILITY_OBJECTS) $(LOCAL_OBJECTS) 

CGI_LIBRARY = $(CGILIBRARY_DIR)/libcgi.a



# -- Default Target -----------------------------------------------------------
#

all: $(BROWSER) $(INDEXER)



# -- Specific Targets ---------------------------------------------------------
#

$(BROWSER): $(BROWSER_OBJECTS) $(UTILITY_OBJECTS) $(CGI_LIBRARY)
	cc -o $@ $(BROWSER_OBJECTS) $(UTILITY_OBJECTS) $(CGI_LIBRARY)

$(INDEXER): $(INDEXER_OBJECTS) $(UTILITY_OBJECTS)
	cc -o $@ $(INDEXER_OBJECTS) $(UTILITY_OBJECTS)

$(CGI_LIBRARY):
	( cd $(CGILIBRARY_DIR); make )


# -- Utility Targets ----------------------------------------------------------
#

clean: 
	rm -f $(BROWSER) $(INDEXER) $(OBJECTS) tags Makefile.bak

tar:
	-rm -f source.tar.gz
	tar cfv source.tar $(LOCAL_SOURCES) $(LOCAL_HEADERS) Makefile $(DOCS) $(SCRIPTS)
	gzip source.tar

untar: source.tar.gz
	zcat source.tar.gz | tar xfv -



# -- Installation Targets -----------------------------------------------------
#

install: $(INSTALL_DIR)/$(BROWSER)

$(INSTALL_DIR)/$(BROWSER): $(BROWSER)
	cp $? $@



# -- Automatic Dependancy List ------------------------------------------------
# 

