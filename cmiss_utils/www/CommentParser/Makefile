#
# Makefile for Comment Stripping Routines...
#


# -- Shell Directives ---------------------------------------------------------
#

SHELL = /bin/sh



# -- Directories --------------------------------------------------------------
#

UTILITY_DIR = ../Utilities
CGILIBRARY_DIR = ../CGIDecode
REGEXP_DIR = ../RegEx



# -- Processor Directives -----------------------------------------------------
#

CC        = cc
CFLAGS    = -ansi -pedantic -O0 -g
CDEFINES  = -DCMISS_ROOT=\"${CMISS_ROOT}\"
CINCLUDES = -I$(UTILITY_DIR) -I$(CGILIBRARY_DIR) -I$(REGEXP_DIR)



# -- Rules --------------------------------------------------------------------
#

.c.o:
	$(CC) -c $(CFLAGS) $(CDEFINES) $(CINCLUDES) -o $@ $<



# -- Resources ----------------------------------------------------------------
#

CGI_LIBRARY = $(CGILIBRARY_DIR)/libcgi.a
REGEXP_LIBRARY = $(REGEXP_DIR)/libregexp.a

UTILITY_SOURCES = $(UTILITY_DIR)/cmiss-lineio.c  $(UTILITY_DIR)/groups.c \
                  $(UTILITY_DIR)/htstring.c      $(UTILITY_DIR)/lineio.c \
                  $(UTILITY_DIR)/rule-parse.c    $(UTILITY_DIR)/rules.c  \
                  $(UTILITY_DIR)/simple-parser.c $(UTILITY_DIR)/tags.c   \
                  $(UTILITY_DIR)/utilities.c
UTILITY_HEADERS = $(UTILITY_SOURCES:.c=.h)
UTILITY_OBJECTS = $(UTILITY_SOURCES:.c=.o)

PARSE_SOURCES = output.c comment-stripper.c
PARSE_HEADERS = output.h
PARSE_OBJECTS = $(PARSE_SOURCES:.c=.o)

LOOKUP_SOURCES = cmiss-lookup.c 
LOOKUP_HEADERS = 
LOOKUP_OBJECTS = $(LOOKUP_SOURCES:.c=.o)

DATABASE_SOURCES = reducers.c
DATABASE_HEADERS = $(DATABASE_SOURCES:.c=.h)
DATABASE_OBJECTS = $(DATABASE_SOURCES:.c=.o)

OUTPUTCMDS_SOURCES = outputcommands.c
OUTPUTCMDS_HEADERS = $(OUTPUTCMDS_SOURCES:.c=.h)
OUTPUTCMDS_OBJECTS = $(OUTPUTCMDS_SOURCES:.c=.o)

TARGETS = comment-stripper cmiss-lookup.cgi

LOCAL_SOURCES = $(PARSE_SOURCES) $(LOOKUP_SOURCES) $(DATABASE_SOURCES) $(OUTPUTCMDS_SOURCES)
LOCAL_HEADERS = $(PARSE_HEADERS) $(LOOKUP_HEADERS) $(DATABASE_HEADERS) $(OUTPUTCMDS_HEADERS)
LOCAL_OBJECTS = $(PARSE_OBJECTS) $(LOOKUP_OBJECTS) $(DATABASE_OBJECTS) $(OUTPUTCMDS_OBJECTS)

HEADERS = $(UTILITY_HEADERS) $(LOCAL_SOURCES)
SOURCES = $(UTILITY_SOURCES) $(LOCAL_HEADERS)
OBJECTS = $(UTILITY_OBJECTS) $(LOCAL_OBJECTS)

DOCS = README
SCRIPTS = comment-stripper.sh



# -- Default Target -----------------------------------------------------------
#

all: $(TARGETS)



# -- Specific Targets ---------------------------------------------------------
#

comment-stripper: $(UTILITY_OBJECTS) $(PARSE_OBJECTS) $(OUTPUTCMDS_OBJECTS) $(DATABASE_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(UTILITY_OBJECTS) $(PARSE_OBJECTS) $(OUTPUTCMDS_OBJECTS) $(DATABASE_OBJECTS) $(CGI_LIBRARY) $(REGEXP_LIBRARY)

cmiss-lookup.cgi: $(LOOKUP_OBJECTS) $(DATABASE_OBJECTS) $(UTILITY_OBJECTS) $(CGI_LIBRARY) $(REGEXP_LIBRARY)
	$(CC) $(CFLAGS) -o $@ $(UTILITY_OBJECTS) $(LOOKUP_OBJECTS) $(OUTPUTCMDS_OBJECTS) $(DATABASE_OBJECTS) $(CGI_LIBRARY) $(REGEXP_LIBRARY)

$(CGI_LIBRARY):
	( cd $(CGILIBRARY_DIR); make )

$(REGEXP_LIBRARY):
	( cd $(REGEXP_DIR); make )



# -- Utility Targets ----------------------------------------------------------
#

clean:
	rm -f $(TARGETS) $(LOCAL_OBJECTS) tags Makefile.bak

tar: 
	-rm -f source.tar.gz
	tar cf source.tar $(LOCAL_SOURCES) $(LOCAL_HEADERS) $(DOCS) $(SCRIPTS) Makefile
	gzip source.tar

untar: source.tar.gz
	zcat source.tar.gz | tar xfv -

tags: $(SOURCES) $(HEADERS) 
	ctags -T $(SOURCES) $(HEADERS)

depend: $(SOURCES) $(HEADERS)
	makedepend $(CDEFINES) $(CINCLUDES) $(SOURCES)



# -- Automatic Dependancy List ------------------------------------------------
#

# DO NOT DELETE
