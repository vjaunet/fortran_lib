#construction des modules
MODULES=$(shell grep -l 'module' SOURCE/*.f90)
MOD=$(MODULES:.f90=.mod)

#source : rech de .f90
SRC=$(wildcard SOURCE/*.f90)

#construction des objets
OBJ=$(SRC:.f90=.o)

#-----------------------------------------------------------
#compilation :
CC=ifort
LIBS= -lfftw3 -llapack

#option de vectorisation et parallelisation:
OPT_para= -vec-report0 -parallel -par-report0
OPT_OMP= #-openmp
#option de debugage :
OPT_Debug= #-traceback -CB -warn alignment -ftrapuv -mp1

ALL:$(MOD) $(OBJ)

%.mod:%.f90
	$(CC) -O3 $(OPT_OMP) $(OPT_para) $(OPT_Debug) -c $^	   $(LIBS)

%.o: %.f90
	$(CC) -O3 $(OPT_OMP) $(OPT_para) $(OPT_Debug) -o $@ -c $^  $(LIBS)

clean:
	rm -rf *.o *~ *.mod

install:
	ar rc libtecplot.a tecplot_io.o
	mv *.mod MOD/.
	mv *.a LIB/.