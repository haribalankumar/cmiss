#
# Makefile for Fortran-check re-parser
#


# -- Shell Directives ---------------------------------------------------------
#

SHELL = /bin/sh



# -- Directory Paths ----------------------------------------------------------
#

UTILITY_DIR = ../Utilities/



# -- Processor Directives -----------------------------------------------------
#

CC        = cc
CFLAGS    = -ansi -pedantic -O0 -g
CDEFINES  = -DCMISS_ROOT=\"${CMISS_ROOT}\" \
            -DCMISS_WWW_URL=\"${CMISS_WWW_URL}\" \
            -DCMISS_WWW_SERVERPATH=\"${CMISS_WWW_SERVERPATH}\"
CINCLUDES = -I$(UTILITY_DIR)



# -- Rules --------------------------------------------------------------------
#

.c.o:
	$(CC) -c $(CFLAGS) $(CDEFINES) $(CINCLUDES) -o $@ $<



# -- Targets ------------------------------------------------------------------
#

UTILITY_SOURCES = $(UTILITY_DIR)lineio.c $(UTILITY_DIR)simple-parser.c
UTILITY_HEADERS = $(UTILITY_SOURCES:.c=.h)
UTILITY_OBJECTS = $(UTILITY_SOURCES:.c=.o)

LOCAL_SOURCES = parser.c parse-fortran-check.c
LOCAL_HEADERS = parser.h
LOCAL_OBJECTS = $(LOCAL_SOURCES:.c=.o)

TARGET = parse-fortran-check

HEADERS = $(UTILITY_HEADERS) $(LOCAL_HEADERS)
SOURCES = $(UTILITY_HEADERS) $(LOCAL_SOURCES)
OBJECTS = $(UTILITY_OBJECTS) $(LOCAL_OBJECTS)



# -- Additional Files ---------------------------------------------------------
#

SCRIPTS = parse-fortran-check.sh
DOCS    = README



# -- Default Target -----------------------------------------------------------
#

all: $(TARGET)



# -- Specific Targets ---------------------------------------------------------
#

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS)



# -- Utility Targets ----------------------------------------------------------
#

clean:
	rm -f $(TARGET) $(OBJECTS) tags Makefile.bak source.tar.gz

tar: 
	-rm -f source.tar.gz
	tar cfv source.tar $(LOCAL_SOURCES) $(LOCAL_HEADERS) $(SCRIPTS) $(DOCS) Makefile
	gzip source.tar

untar: source.tar.gz
	zcat source.tar.gz | tar xfv -

tags: $(SOURCES) $(HEADERS) 
	ctags -T $(SOURCES) $(HEADERS)

depend: $(SOURCES) $(HEADERS)
	makedepend $(CDEFINES) $(CINCLUDES) $(SOURCES)



# -- Automatic Dependancy List ------------------------------------------------
#

