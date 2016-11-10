############################################################

#source : rech de .f90
SRC=$(wildcard SOURCE/*.f90)

#construction des objets
OBJ=$(SRC:.f90=.mod)

#-----------------------------------------------------------

#option de vectorisation et parallelisation:
# OPT_para= -vec-report0 -parallel -par-report0
OPT_OMP= -fopenmp

#option de debugage :
Debug= -Wall -fcheck=bounds

IPATH=-I./ -I/usr/local/include
LPATH=-L./
LIBS= -lfftw3 -llapack -lstat -lspectral -lpod $(LPATH)

CC = gfortran -O2 $(IPATH)
CC += $(OPT_para)
CC += $(OPT_OMP)

ALL:spectral pod spectralpod stat interpol tecplot qsort piv_data press_data netcdf

debug: CC += $(Debug)
debug: ALL

netcdf:lib_netcdf.o lib_netcdf.mod
lib_netcdf.mod:SOURCE/lib_netcdf.f90 lib_netcdf.o
	@true
lib_netcdf.o:SOURCE/lib_netcdf.f90
	$(CC) -c $^ -lnetcdff -lnetcdff

spectral:lib_spectral.mod
lib_spectral.mod:SOURCE/lib_spectral.f90
	$(CC) -c $^ -lfftw3;

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
	$(CC) -c $^;

piv_data:stat pod qsort netcdf lib_piv_data.mod
lib_piv_data.mod:SOURCE/lib_piv_data.f90
	$(CC) -c $^ $(LIBS) -llibnetcdf -lnetcdff -lnetcdff;

press_data:lib_press_data.mod
lib_press_data.mod:SOURCE/lib_press_data.f90
	$(CC) -c $^ $(LIBS);


install:
	ar rc liblibnetcdf.a lib_netcdf.o
	ar rc libtecplot.a lib_tecplot_io.o
	ar rc libinterpol.a lib_interpol.o
	ar rc libstat.a lib_stat.o
	ar rc libpod.a lib_pod.o
	ar rc libspectralpod.a lib_spectralPOD.o lib_spectral.o
	ar rc libpivdata.a \
	lib_piv_data.o lib_pod.o lib_stat.o qsort_c.o
	ar rc libspectral.a lib_spectral.o
	ar rc libpressdata.a lib_press_data.o
	mkdir -p MOD
	mkdir -p LIB
	cp *.mod MOD/.
	mv *.a LIB/.
#for file in $(ls MOD/*); do echo ln -s $(pwd)/$file /usr/local/include/$file; done

clean:
	rm -rf *.o *.mod *.a SOURCE/*~
