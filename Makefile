############################################################

#source : rech de .f90
SRC=$(wildcard SOURCE/*.f90)

#construction des objets
OBJ=$(SRC:.f90=.mod)

#-----------------------------------------------------------

#option de vectorisation et parallelisation:
OPT_OMP= -fopenmp

#option de debugage :
Debug= -Wall -fcheck=all -g

IPATH=-I./ -I/usr/include -I/usr/local/Cellar/fftw/3.3.8/include/
LPATH=-L./ -L/usr/lib/x86_64-linux-gnu

CC = gfortran -O2 $(IPATH)
CC += $(OPT_OMP)

ALL:static install

static:spectral pod spectralpod stat interpol tecplot qsort piv_data press_data #netcdf

debug: CC += $(Debug)
debug: ALL

netcdf:lib_netcdf.o lib_netcdf.mod
lib_netcdf.mod:SOURCE/lib_netcdf.f90 lib_netcdf.o
	@true
lib_netcdf.o:SOURCE/lib_netcdf.f90
	$(CC) -c $^

spectral:lib_spectral.mod
lib_spectral.mod:SOURCE/lib_spectral.f90 lib_spectral.o
	@true
lib_spectral.o:SOURCE/lib_spectral.f90
	$(CC) -c $^

stat:lib_stat.mod
lib_stat.mod:SOURCE/lib_stat.f90 lib_stat.o
	@true
lib_stat.o:SOURCE/lib_stat.f90
	$(CC) -c $^

interpol:lib_interpol.mod
lib_interpol.mod:SOURCE/lib_interpol.f90 lib_interpol.o
	@true
lib_interpol.o:SOURCE/lib_interpol.f90
	$(CC) -c $^

tecplot:lib_tecplot_io.mod
lib_tecplot_io.mod:SOURCE/lib_tecplot_io.f90 lib_tecplot_io.o
	@true
lib_tecplot_io.o:SOURCE/lib_tecplot_io.f90
	$(CC) -c $^

pod:lib_pod.mod
lib_pod.mod:SOURCE/lib_pod.f90 lib_pod.o
	@true
lib_pod.o:SOURCE/lib_pod.f90
	$(CC) -c $^

spectralpod:spectral pod lib_spectralpod.mod
lib_spectralpod.mod:SOURCE/lib_spectralPOD.f90 lib_spectralPOD.o
	@true
lib_spectralPOD.o:SOURCE/lib_spectralPOD.f90
	$(CC) -c $^

qsort:qsort_c.mod
qsort_c.mod:SOURCE/qsort_c.f90 qsort_c.o
	@true
qsort_c.o:SOURCE/qsort_c.f90
	$(CC) -c $^

piv_data:stat pod qsort lib_piv_data.mod #netcdf
lib_piv_data.mod:SOURCE/lib_piv_data.f90 lib_piv_data.o
	@true
lib_piv_data.o:SOURCE/lib_piv_data.f90
	$(CC) -c $^

press_data:lib_press_data.mod
lib_press_data.mod:SOURCE/lib_press_data.f90 lib_press_data.o
	@true
lib_press_data.o:SOURCE/lib_press_data.f90
	$(CC) -c $^


install:
	ar rc libtecplot.a lib_tecplot_io.o
	ar rc libinterpol.a lib_interpol.o
	ar rc libstat.a lib_stat.o
	ar rc libpod.a lib_pod.o
	ar rc libspectralpod.a lib_spectralPOD.o lib_spectral.o
	ar rc libpivdata.a lib_piv_data.o qsort_c.o lib_pod.o lib_stat.o #lib_netcdf.o
	ar rc libspectral.a lib_spectral.o
	ar rc libpressdata.a lib_press_data.o
	mkdir -p MOD
	mkdir -p LIB
	cp *.mod MOD/.
	mv *.a LIB/.

        #ar rc liblibnetcdf.a lib_netcdf.o

clean:
	rm -rf *.o *.mod *.a SOURCE/*~
