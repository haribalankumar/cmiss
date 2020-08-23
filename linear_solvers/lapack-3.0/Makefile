# -*- Makefile -*-

#
#  Top Level Makefile for LAPACK
#  Version 3.0
#  June 30, 1999
#

include make.inc

lib: $(LIBDIR) lapacklib

$(LIBDIR):
	mkdir -p $(LIBDIR)

lapacklib:
	( cd SRC; $(MAKE) )

clean:
	( cd SRC; $(MAKE) clean )

