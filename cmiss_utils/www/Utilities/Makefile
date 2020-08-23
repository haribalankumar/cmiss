#
# Makefile for various Utility Functions/Class/Modules
#

CC        = cc
CFLAGS    = -ansi -pedantic -Wall  -O0
AR        = ar
AFLAGS    = rvf


.c.o:
	$(CC) -c $(CFLAGS) $(CDEFINES) $(CINCLUDES) -o $@ $<

LIB = libUte.a
OBJ = lineio.o simple-parser.o tags.o cmiss-lineio.o utilities.o htstring.o\
	rules.o rule-parse.o groups.o 
TEST = test-lineio
TOBJ = test-lineio.o


$(LIB): $(OBJ)
	$(AR) $(AFLAGS) $@ $(OBJ)

$(TEST): $(TOBJ) $(LIB)
	cc -o $@ $(CFLAGS) $(TOBJ) -L. -lUte

clean:
	rm -f $(OBJ) *~ *.a $(TEST)

# -- Automatic Dependancy List ------------------------------------------------
#

# DO NOT DELETE

