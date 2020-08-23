
.SUFFIXES: .o .f .F .f90

FC=f90
#FFLAGS=-O0 -g -check_bounds
LFLAGS=
OBJ=constants.o fileio.o thickness.o slicefile.o spline.o\
    ipelemfile.o ipnodefile.o createnodes.o create_thickness.o 
EXE=create_thickness

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

.f90.o:
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *~ *.o *.mod
