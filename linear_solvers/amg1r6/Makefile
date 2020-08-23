# Makefile for AMG code library

BUILDING_AMG = true

include ../Makefile.inc

# Object Modules

TIME = amgtime.o
NDRV = ndrv.o
AMG  = amg1r6.o
WRAP = amg_wrapper.o
PRMS = amg_set_parms.o
CPY  = amg_copy_matrix.o

OBJS = $(TIME) $(NDRV) $(AMG) $(WRAP) $(PRMS) $(CPY) 

AMGLIB = libamg$(SUFFIX).a

LIB = $(LIBDIR)/$(AMGLIB)

# ============================================================================
# -----------------
#   TARGETS:
# -----------------

$(LIB): $(OBJS)
	@echo "Executing $(AR) $(AFLAGS) $(LIB) on objects";\
	$(AR) $(AFLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

# Remove objects so that they aren't reused for a different build.
.INTERMEDIATE: $(OBJS)

clean:
	@echo "Removing objects for $(EXE) ... " ;\
	$(RM) -f $(OBJS) *~

distclean:
	@echo "Removing objects and library ";\
	$(RM) -f $(OBJS) $(LIB)

# ============================================================================
