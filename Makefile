############################################################

#source : rech de .f90
SRC=$(wildcard SOURCE/*.f90)

#construction des objets
OBJ=$(SRC:.f90=.mod)

#-----------------------------------------------------------

#option de vectorisation et parallelisation:
OPT_para= -vec-report0 -parallel -par-report0
OPT_OMP= -openmp

#option de debugage :
Debug= -traceback -CB -warn alignment -ftrapuv -mp1

IPATH=-I./
LPATH=-L./
LIBS= -lfftw3 -mkl -lstat -lspectral -lpod $(LPATH)

CC = ifort -O2 $(IPATH)
CC += $(OPT_para)
CC += $(OPT_OMP)

ALL:spectral pod spectralpod stat interpol tecplot qsort piv_data press_data

debug: CC += $(Debug)
debug: ALL

spectral:lib_spectral.mod
lib_spectral.mod:SOURCE/lib_spectral.f90
	$(CC) -c $^ $(LIBS);

stat:lib_stat.mod
lib_stat.mod:SOURCE/lib_stat.f90
	$(CC) -c $^ $(LIBS);

interpol:lib_interpol.mod
lib_interpol.mod:SOURCE/lib_interpol.f90
	$(CC) -c $^ $(LIBS);

tecplot:lib_tecplot_io.mod
lib_tecplot_io.mod:SOURCE/lib_tecplot_io.f90
	$(CC) -c $^ $(LIBS);

pod:lib_pod.mod
lib_pod.mod:SOURCE/lib_pod.f90
	$(CC) -c $^ $(LIBS);

spectralpod:spectral pod lib_spectralpod.mod
lib_spectralpod.mod:SOURCE/lib_spectralPOD.f90
	$(CC) -c $^ $(LIBS);

qsort:qsort_c.mod
qsort_c.mod:SOURCE/qsort_c.f90
	$(CC) -c $^ $(LIBS);

piv_data:stat pod qsort lib_piv_data.mod
lib_piv_data.mod:SOURCE/lib_piv_data.f90
	$(CC) -c $^ $(LIBS);

press_data:lib_press_data.mod
lib_press_data.mod:SOURCE/lib_press_data.f90
	$(CC) -c $^ $(LIBS);


install:
	ar rc libtecplot.a lib_tecplot_io.o
	ar rc libinterpol.a lib_interpol.o
	ar rc libstat.a lib_stat.o
	ar rc libpod.a lib_pod.o
	ar rc libspectralpod.a lib_spectralPOD.o lib_spectral.o
	ar rc libpivdata.a \
	lib_piv_data.o lib_pod.o lib_stat.o qsort_c.o
	ar rc libspectral.a lib_spectral.o
	ar rc libpressdata.a lib_press_data.o
	cp *.mod MOD/.
	mv *.a LIB/.

clean:
	rm -rf *.o *.mod *.a SOURCE/*~
