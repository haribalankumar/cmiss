#
# Makefile for various Utility Functions/Class/Modules
#

CC        = cc
CFLAGS    = -ansi -pedantic -O0
AR        = ar
AFLAGS    = rvf


.c.o:
	$(CC) -c $(CFLAGS) $(CINCLUDES) -o $@ $<

LIB = libcgi.a
OBJ = cgi-decode.o

TEST = test-cgi
TOBJ = test-cgi.o

$(LIB): $(OBJ)
	$(AR) $(AFLAGS) $(LIB) $(OBJ)

$(TEST): $(TOBJ) $(LIB)
	$(CC) -o $@ $(TOBJ) $(CFLAGS) -L. -lcgi

clean:
	rm -f *.o *.a $(TEST) *~

