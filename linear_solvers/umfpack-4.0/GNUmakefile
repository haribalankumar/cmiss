# -*- Makefile -*-
#-------------------------------------------------------------------------------
# UMFPACK Makefile for compiling on Unix systems (for GNU Make)
#-------------------------------------------------------------------------------

# Include generic Make rules
BUILDING_UMFPACK = true

include ../Makefile.inc

LIBUMF = $(LIBDIR)/libumfpack-4.0$(SUFFIX).a

# C compiler rule (CFLAGS and CONFIG are set in the Make.* file above):
C = $(CC) $(CFLAGS) $(CONFIG) -DNDEBUG

#-------------------------------------------------------------------------------
# no modifications should be necessary below this line.
#-------------------------------------------------------------------------------

# non-user-callable umf_*.[ch] files:
UMF = umf_analyze umf_apply_order umf_assemble umf_blas3_update \
	umf_build_tuples umf_build_tuples_usage umf_colamd umf_create_element \
	umf_dump umf_extend_front umf_free umf_garbage_collection \
	umf_get_memory umf_init_front umf_is_permutation umf_kernel \
	umf_kernel_init umf_kernel_init_usage umf_kernel_wrapup \
	umf_local_search umf_lsolve umf_ltsolve umf_malloc \
	umf_mem_alloc_element umf_mem_alloc_head_block \
	umf_mem_alloc_tail_block umf_mem_free_tail_block \
	umf_mem_init_memoryspace umf_order_front_tree umf_report_perm \
	umf_realloc umf_report_vector umf_row_search umf_scale_column \
	umf_set_stats umf_solve umf_symbolic_usage umf_transpose \
	umf_tuple_lengths umf_usolve umf_utsolve umf_valid_numeric \
	umf_valid_symbolic

# non-user-callable, created from umf_ltsolve, umf_utsolve, and umf_triplet:
UMFH = umf_lhsolve umf_uhsolve umf_triplet_map_x umf_triplet_map_nox \
	umf_triplet_nomap_x umf_triplet_nomap_nox

INTERNAL = $(UMF) $(UMFH)

# user-callable umfpack_*.[ch] files:
UMFPACK = umfpack_col_to_triplet umfpack_defaults umfpack_free_numeric \
	umfpack_free_symbolic umfpack_get_numeric umfpack_get_lunz \
	umfpack_get_symbolic umfpack_numeric umfpack_qsymbolic \
	umfpack_report_control umfpack_report_info umfpack_report_matrix \
	umfpack_report_numeric umfpack_report_perm umfpack_report_status \
	umfpack_report_symbolic umfpack_report_triplet \
	umfpack_report_vector umfpack_solve umfpack_symbolic \
	umfpack_transpose umfpack_triplet_to_col

# user-callable, created from umfpack_solve.c (umfpack_wsolve.h exists, though):
UMFPACKW = umfpack_wsolve

USER = $(UMFPACK) $(UMFPACKW)

# user-callable, only one version for int/long, real/complex, *.[ch] files:
GENERIC = umfpack_timer


#-------------------------------------------------------------------------------
# include files:
#-------------------------------------------------------------------------------

INC = umfpack.h umf_config.h umf_version.h umf_internal.h umf_triplet.h \
	$(addsuffix .h,$(UMF)) $(addsuffix .h,$(USER)) \
	$(addsuffix .h,$(GENERIC))


#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------

DI = $(subst umf_,umf_di_,$(INTERNAL)) $(subst umfpack_,umfpack_di_,$(USER))
DL = $(subst umf_,umf_dl_,$(INTERNAL)) $(subst umfpack_,umfpack_dl_,$(USER))
ZI = $(subst umf_,umf_zi_,$(INTERNAL)) $(subst umfpack_,umfpack_zi_,$(USER))
ZL = $(subst umf_,umf_zl_,$(INTERNAL)) $(subst umfpack_,umfpack_zl_,$(USER))

# defaults, for Unix:
ifndef OBJEXT
OBJEXT = .o
endif
ifndef OUTPUT
OUTPUT = -o
endif

# add the appropriate .o or .obj extension
DI := $(addsuffix $(OBJEXT),$(DI))
DL := $(addsuffix $(OBJEXT),$(DL))
ZI := $(addsuffix $(OBJEXT),$(ZI))
ZL := $(addsuffix $(OBJEXT),$(ZL))
GN := $(addsuffix $(OBJEXT),$(GENERIC))


#-------------------------------------------------------------------------------
# compile the C-callable library (umfpack.a), the demos, and the mexFunction
#-------------------------------------------------------------------------------

lib: $(LIBDIR) $(LIBUMF)

all: $(LIBUMF) \
	umfpack_di_demo umfpack_zi_demo umfpack_dl_demo umfpack_zl_demo \
	umfpack_simple $(UMFPACKMEX)

$(LIBDIR):
	mkdir -p $(LIBDIR)

#-------------------------------------------------------------------------------
# compile each routine in the DI version
#-------------------------------------------------------------------------------

umf_di_%$(OBJEXT): umf_%.c $(INC)
	$(C) -DDINT -c $< $(OUTPUT) $@

umf_di_%hsolve$(OBJEXT): umf_%tsolve.c $(INC)
	$(C) -DDINT -DCONJUGATE_SOLVE -c $< $(OUTPUT) $@

umf_di_triplet_map_x$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DDINT -DDO_MAP -DDO_VALUES -c $< $(OUTPUT) $@

umf_di_triplet_map_nox$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DDINT -DDO_MAP -c $< $(OUTPUT) $@

umf_di_triplet_nomap_x$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DDINT -DDO_VALUES -c $< $(OUTPUT) $@

umf_di_triplet_nomap_nox$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DDINT -c $< $(OUTPUT) $@

umfpack_di_wsolve$(OBJEXT): umfpack_solve.c $(INC)
	$(C) -DDINT -DWSOLVE -c $< $(OUTPUT) $@

umfpack_di_%$(OBJEXT): umfpack_%.c $(INC)
	$(C) -DDINT -c $< $(OUTPUT) $@


#-------------------------------------------------------------------------------
# compile each routine in the DL version
#-------------------------------------------------------------------------------

umf_dl_%$(OBJEXT): umf_%.c $(INC)
	$(C) -DDLONG -c $< $(OUTPUT) $@

umf_dl_%hsolve$(OBJEXT): umf_%tsolve.c $(INC)
	$(C) -DDLONG -DCONJUGATE_SOLVE -c $< $(OUTPUT) $@

umf_dl_triplet_map_x$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DDLONG -DDO_MAP -DDO_VALUES -c $< $(OUTPUT) $@

umf_dl_triplet_map_nox$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DDLONG -DDO_MAP -c $< $(OUTPUT) $@

umf_dl_triplet_nomap_x$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DDLONG -DDO_VALUES -c $< $(OUTPUT) $@

umf_dl_triplet_nomap_nox$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DDLONG -c $< $(OUTPUT) $@

umfpack_dl_wsolve$(OBJEXT): umfpack_solve.c $(INC)
	$(C) -DDLONG -DWSOLVE -c $< $(OUTPUT) $@

umfpack_dl_%$(OBJEXT): umfpack_%.c $(INC)
	$(C) -DDLONG -c $< $(OUTPUT) $@


#-------------------------------------------------------------------------------
# compile each routine in the ZI version
#-------------------------------------------------------------------------------

umf_zi_%$(OBJEXT): umf_%.c $(INC)
	$(C) -DZINT -c $< $(OUTPUT) $@

umf_zi_%hsolve$(OBJEXT): umf_%tsolve.c $(INC)
	$(C) -DZINT -DCONJUGATE_SOLVE -c $< $(OUTPUT) $@

umf_zi_triplet_map_x$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DZINT -DDO_MAP -DDO_VALUES -c $< $(OUTPUT) $@

umf_zi_triplet_map_nox$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DZINT -DDO_MAP -c $< $(OUTPUT) $@

umf_zi_triplet_nomap_x$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DZINT -DDO_VALUES -c $< $(OUTPUT) $@

umf_zi_triplet_nomap_nox$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DZINT -c $< $(OUTPUT) $@

umfpack_zi_wsolve$(OBJEXT): umfpack_solve.c $(INC)
	$(C) -DZINT -DWSOLVE -c $< $(OUTPUT) $@

umfpack_zi_%$(OBJEXT): umfpack_%.c $(INC)
	$(C) -DZINT -c $< $(OUTPUT) $@


#-------------------------------------------------------------------------------
# compile each routine in the ZL version
#-------------------------------------------------------------------------------

umf_zl_%$(OBJEXT): umf_%.c $(INC)
	$(C) -DZLONG -c $< $(OUTPUT) $@

umf_zl_%hsolve$(OBJEXT): umf_%tsolve.c $(INC)
	$(C) -DZLONG -DCONJUGATE_SOLVE -c $< $(OUTPUT) $@

umf_zl_triplet_map_x$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DZLONG -DDO_MAP -DDO_VALUES -c $< $(OUTPUT) $@

umf_zl_triplet_map_nox$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DZLONG -DDO_MAP -c $< $(OUTPUT) $@

umf_zl_triplet_nomap_x$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DZLONG -DDO_VALUES -c $< $(OUTPUT) $@

umf_zl_triplet_nomap_nox$(OBJEXT): umf_triplet.c $(INC)
	$(C) -DZLONG -c $< $(OUTPUT) $@

umfpack_zl_wsolve$(OBJEXT): umfpack_solve.c $(INC)
	$(C) -DZLONG -DWSOLVE -c $< $(OUTPUT) $@

umfpack_zl_%$(OBJEXT): umfpack_%.c $(INC)
	$(C) -DZLONG -c $< $(OUTPUT) $@


#-------------------------------------------------------------------------------
# Create the generic routines (GN) using a generic rule
#-------------------------------------------------------------------------------

%$(OBJEXT): %.c $(INC)
	$(C) -c $< $(OUTPUT) $@


#-------------------------------------------------------------------------------
# Clean-up
#-------------------------------------------------------------------------------

clean:
	rm -rf *$(OBJEXT) $(CLEAN_OBJS) $(LIBUMF)


#-------------------------------------------------------------------------------
# Create the umfpack.a library
#-------------------------------------------------------------------------------

OBJ := $(ZI) $(ZL) $(DI) $(DL) $(GN)

$(LIBUMF): $(OBJ)
	$(AR) $(AFLAGS) $(LIBUMF) $^
	@ echo "If ranlib isn't found then it's not needed."
	- $(RANLIB) $(LIBUMF)
	- cat License

# Remove objects so that they aren't reused for a different build.
.INTERMEDIATE: $(OBJ)

#-------------------------------------------------------------------------------
# Create the five demo programs, run them, and compare the output
#-------------------------------------------------------------------------------

DEMO = umfpack_di_demo.c umfpack_dl_demo.c umfpack_zi_demo.c umfpack_zl_demo.c

%_demo.c:  umfpack_xx_demo.c %_demo.sed
	- sed -f $*_demo.sed < umfpack_xx_demo.c > $*_demo.c

%_demo: %_demo.c $(INC) $(LIBUMF)
	$(C) $(OUTPUT) $*_demo $*_demo.c $(LIBUMF) $(SCILIB)
	./$*_demo > my_$*_demo.out
	- diff $*_demo.out my_$*_demo.out

umfpack_simple: umfpack_simple.c $(INC) $(LIBUMF)
	$(C) $(OUTPUT) umfpack_simple umfpack_simple.c $(LIBUMF) $(SCILIB)
	- cat License
	./umfpack_simple

demo: $(DEMO)


#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

# Keep the umfpack_*_demo.c files, in case the user doesn't have sed.

purge: clean
	- rm $(LIBUMF) umfpack_simple 
	- rm umfpack_di_demo umfpack_dl_demo
	- rm umfpack_zi_demo umfpack_zl_demo
	- rm umfpack*.mex* umfpack.dll umfpackmex.a
	- rm UserGuide.aux UserGuide.bbl UserGuide.blg
	- rm UserGuide.log UserGuide.toc
	- rm my_umfpack_di_demo.out my_umfpack_dl_demo.out
	- rm my_umfpack_zi_demo.out my_umfpack_zl_demo.out

# This removes some files in the distribution, but which can be regenerated
allpurge: purge
	- rm umfpack_di_demo.c umfpack_dl_demo.c
	- rm umfpack_zi_demo.c umfpack_zl_demo.c
	- rm UserGuide.ps UserGuide.pdf UserGuide.dvi UserGuide.tex


#-------------------------------------------------------------------------------
# Create the umfpack mexFunction for MATLAB.
#-------------------------------------------------------------------------------

MEXDI = $(subst umf_,umf_md_,$(INTERNAL)) $(subst umfpack_,umfpack_md_,$(USER))
MEXZI = $(subst umf_,umf_mz_,$(INTERNAL)) $(subst umfpack_,umfpack_mz_,$(USER))
MEXDI := $(addsuffix $(OBJEXT),$(MEXDI))
MEXZI := $(addsuffix $(OBJEXT),$(MEXZI))

# Note that mex has no "-o" option, thus the need for "mv" commands ...

umf_md_%$(OBJEXT): umf_%.c $(INC)
	$(MEX) -DDINT -c $<
	- mv umf_$*$(OBJEXT) $@

umf_md_%hsolve$(OBJEXT): umf_%tsolve.c $(INC)
	$(MEX) -DDINT -DCONJUGATE_SOLVE -c $<
	- mv umf_$*tsolve$(OBJEXT) $@

umf_md_triplet_map_x$(OBJEXT): umf_triplet.c $(INC)
	$(MEX) -DDINT -DDO_MAP -DDO_VALUES -c $<
	- mv umf_triplet$(OBJEXT) $@

umf_md_triplet_map_nox$(OBJEXT): umf_triplet.c $(INC)
	$(MEX) -DDINT -DDO_MAP -c $<
	- mv umf_triplet$(OBJEXT) $@

umf_md_triplet_nomap_x$(OBJEXT): umf_triplet.c $(INC)
	$(MEX) -DDINT -DDO_VALUES -c $<
	- mv umf_triplet$(OBJEXT) $@

umf_md_triplet_nomap_nox$(OBJEXT): umf_triplet.c $(INC)
	$(MEX) -DDINT -c $<
	- mv umf_triplet$(OBJEXT) $@

umfpack_md_wsolve$(OBJEXT): umfpack_solve.c $(INC)
	$(MEX) -DDINT -DWSOLVE -c $<
	- mv umfpack_solve$(OBJEXT) $@

umfpack_md_%$(OBJEXT): umfpack_%.c $(INC)
	$(MEX) -DDINT -c $<
	- mv umfpack_$*$(OBJEXT) $@

umf_mz_%$(OBJEXT): umf_%.c $(INC)
	$(MEX) -DZINT -c $<
	- mv umf_$*$(OBJEXT) $@

umf_mz_%hsolve$(OBJEXT): umf_%tsolve.c $(INC)
	$(MEX) -DZINT -DCONJUGATE_SOLVE -c $<
	- mv umf_$*tsolve$(OBJEXT) $@

umf_mz_triplet_map_x$(OBJEXT): umf_triplet.c $(INC)
	$(MEX) -DZINT -DDO_MAP -DDO_VALUES -c $<
	- mv umf_triplet$(OBJEXT) $@

umf_mz_triplet_map_nox$(OBJEXT): umf_triplet.c $(INC)
	$(MEX) -DZINT -DDO_MAP -c $<
	- mv umf_triplet$(OBJEXT) $@

umf_mz_triplet_nomap_x$(OBJEXT): umf_triplet.c $(INC)
	$(MEX) -DZINT -DDO_VALUES -c $<
	- mv umf_triplet$(OBJEXT) $@

umf_mz_triplet_nomap_nox$(OBJEXT): umf_triplet.c $(INC)
	$(MEX) -DZINT -c $<
	- mv umf_triplet$(OBJEXT) $@

umfpack_mz_wsolve$(OBJEXT): umfpack_solve.c $(INC)
	$(MEX) -DZINT -DWSOLVE -c $<
	- mv umfpack_solve$(OBJEXT) $@

umfpack_mz_%$(OBJEXT): umfpack_%.c $(INC)
	$(MEX) -DZINT -c $<
	- mv umfpack_$*$(OBJEXT) $@

umfpack_mtimer$(OBJEXT):  umfpack_timer.c $(INC)
	- cp -f umfpack_timer.c umfpack_mtimer.c
	$(MEX) -v -c umfpack_mtimer.c
	- rm umfpack_mtimer.c

umfpackmex.a: $(MEXZI) $(MEXDI) $(INC) umfpack_mtimer$(OBJEXT)
	- ar cr umfpackmex.a $(MEXZI) $(MEXDI) umfpack_mtimer$(OBJEXT)

umfpack: umfpackmex.c umfpackmex.a
	$(MEX) -output umfpack umfpackmex.c umfpackmex.a

# Try the no-transpose version of x=b/A:
umfpack2: umfpackmex.c umfpackmex.a
	$(MEX) -DNO_TRANSPOSE_FORWARD_SLASH -v -output umfpack2 umfpackmex.c umfpackmex.a


#-------------------------------------------------------------------------------
# Create the User Guide.  This assumes you have sed, latex, and ghostscript.
# The User Guide comes pre-made (UserGuide.ps and UserGuide.pdf) so you
# shouldn't need to create it here.
#-------------------------------------------------------------------------------

UserGuide.tex:  UserGuide.stex UserGuide.sed1 UserGuide.sed2 $(INC) umfpack_simple.c
	sed -f UserGuide.sed1 < UserGuide.stex | sed -f UserGuide.sed2 | expand -8 > UserGuide.tex

UserGuide.ps: UserGuide.tex UserGuide.bib
	latex UserGuide
	bibtex UserGuide
	latex UserGuide
	latex UserGuide
	dvips UserGuide

UserGuide.pdf:  UserGuide.ps
	ps2pdf UserGuide.ps UserGuide.pdf

dist:  UserGuide.pdf $(DEMO)
	- rm UserGuide.aux UserGuide.bbl UserGuide.blg UserGuide.log UserGuide.toc
