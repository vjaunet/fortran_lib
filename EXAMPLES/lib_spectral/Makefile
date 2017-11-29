#nom de l'executable
EXEC = spectral_prog

#construction des modules
MODULES=$(shell grep -l 'module' *.f90)
MOD=$(MODULES:.f90=.mod)

#source : rech de .f90
SRC=$(wildcard *.f90)

#construction des objets
OBJ=$(SRC:.f90=.o)

#-----------------------------------------------------------
#compilation :

# those folders need to be changed according to your
# fftw3 install path
IPATH=-I../../MOD -I/usr/local/include
LPATH=-L../../LIB -L/usr/local/lib

CC=gfortran $(IPATH)
LIBS= -lspectral -lfftw3 $(LPATH) #

#option de vectorisation et parallelisation:
# OPT_OMP= -fopenmp

#option de debugage :
#OPT_Debug= -g -Wall

ALL:$(MOD) $(EXEC)

$(EXEC):$(OBJ)
	$(CC) -O3 $(OPT_OMP) $(OPT_para) $(OPT_Debug) -o $@ $^     $(LIBS)

%.mod:%.f90
	$(CC) -O3 $(OPT_OMP) $(OPT_para) $(OPT_Debug) -c $^	   $(LIBS)

%.o: %.f90
	$(CC) -O3 $(OPT_OMP) $(OPT_para) $(OPT_Debug) -o $@ -c $^  $(LIBS)

clean:
	rm -rf *.o *~ *.mod
	rm -rf $(EXEC)
