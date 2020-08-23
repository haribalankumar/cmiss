#
# Makefile for the CMISS Unused Var converter
#


# -- Shell Directives ---------------------------------------------------------
#

SHELL = /bin/sh



# -- Processor Directives -----------------------------------------------------
#

CC        = cc
CFLAGS    = -ansi -pedantic -O0
CDEFINES  = 
CINCLUDES = 



# -- Rules --------------------------------------------------------------------
#

.c.o:
	$(CC) -c $(CFLAGS) $(CDEFINES) $(CINCLUDES) -o $@ $<



# -- Targets ------------------------------------------------------------------
#

TARGET = libregexp.a

SOURCES = regexp.c regsub.c regerror.c
HEADERS = regexp.h regmagic.h
OBJECTS = $(SOURCES:.c=.o)



# -- Additional Files ---------------------------------------------------------
#

DOCS    = README



# -- Default Target -----------------------------------------------------------
#

all: $(TARGET)



# -- Specific Targets ---------------------------------------------------------
#

$(TARGET): $(OBJECTS)
	-rm -f $@
	ar cq $@ $(OBJECTS)



# -- Utility Targets ----------------------------------------------------------
#

clean:
	rm -f $(TARGET) $(OBJECTS) tags Makefile.bak source.tar.gz

tar: 
	-rm -f source.tar.gz
	tar cfv source.tar $(SOURCES) $(HEADERS) $(SCRIPTS) $(DOCS) Makefile
	gzip source.tar

untar: source.tar.gz
	zcat source.tar.gz | tar xfv -

tags: $(SOURCES) $(HEADERS) 
	ctags -T $(SOURCES) $(HEADERS)

depend: $(SOURCES) $(HEADERS)
	makedepend $(CDEFINES) $(CINCLUDES) $(SOURCES)



# -- Automatic Dependancy List ------------------------------------------------
#

