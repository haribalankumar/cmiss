#
# Makefile for Call Tree program
#


# -- Shell Directives ---------------------------------------------------------
#

SHELL = /bin/sh



# -- Directories --------------------------------------------------------------
#

UTILITY_DIR = ../Utilities
CGI_DIR = ../CGIDecode
REDUCERS_DIR = ../CommentParser



# -- Processor Directives -----------------------------------------------------
#

CC        = cc
CFLAGS    = -ansi -pedantic -O0
CDEFINES  =
CINCLUDES = -I$(UTILITY_DIR) -I$(CGI_DIR) -I$(REDUCERS_DIR)



# -- Rules --------------------------------------------------------------------
#

.c.o:
	$(CC) -c $(CFLAGS) $(CDEFINES) $(CINCLUDES) -o $@ $<



# -- Targets ------------------------------------------------------------------
#

UTILITY_SOURCES = $(UTILITY_DIR)/cmiss-lineio.c  $(UTILITY_DIR)/groups.c \
                  $(UTILITY_DIR)/htstring.c      $(UTILITY_DIR)/lineio.c \
                  $(UTILITY_DIR)/rule-parse.c    $(UTILITY_DIR)/rules.c  \
                  $(UTILITY_DIR)/simple-parser.c $(UTILITY_DIR)/tags.c   \
                  $(UTILITY_DIR)/utilities.c
UTILITY_HEADERS = $(UTILITY_SOURCES:.c=.h)
UTILITY_OBJECTS = $(UTILITY_SOURCES:.c=.o)

LOCAL_SOURCES = call-tree.c collection.c
LOCAL_HEADERS = collection.h
LOCAL_OBJECTS = $(LOCAL_SOURCES:.c=.o)

REDUCER_SOURCES = 
REDUCER_HEADERS = $(REDUCERS_DIR)/reducers.h
REDUCER_OBJECTS = $(REDUCER_HEADERS:.h=.o)

TARGET = call-tree

SOURCES = $(UTILITY_SOURCES) $(LOCAL_SOURCES) $(REDUCER_SOURCES)
HEADERS = $(UTILITY_HEADERS) $(LOCAL_HEADERS) $(REDUCER_HEADERS)
OBJECTS = $(UTILITY_OBJECTS) $(LOCAL_OBJECTS) $(REDUCER_OBJECTS)

DOCS = README

SCRIPTS = call-tree.sh



# -- Default Target -----------------------------------------------------------
#

all: $(TARGET)



# -- Specific Targets ---------------------------------------------------------
#

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS)

$(REDUCER_OBJECTS):
	(cd $(REDUCERS_DIR); make $(REDUCER_OBJECTS))



# -- Utility Targets ----------------------------------------------------------
#

clean:
	rm -f $(TARGET) $(OBJECTS) tags Makefile.bak source.tar.gz

tar: 
	-rm -f source.tar.gz
	tar cf source.tar $(LOCAL_SOURCES) $(LOCAL_HEADERS) Makefile $(DOCS) $(SCRIPTS)
	gzip source.tar

untar: source.tar.gz
	zcat source.tar.gz | tar xf -

tags: $(SOURCES) $(HEADERS) 
	ctags -T $(SOURCES) $(HEADERS)

depend: $(SOURCES) $(HEADERS)
	makedepend $(CDEFINES) $(CINCLUDES) $(SOURCES)



# -- Automatic Dependancy List ------------------------------------------------
#

