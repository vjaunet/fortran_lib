############################################################

#source : rech de .f90
SRC=$(wildcard SOURCE/*.f90)

#construction des objets
OBJ=$(SRC:.f90=.o)

#-----------------------------------------------------------

#option de vectorisation et parallelisation:
OPT_para= -vec-report0 -parallel -par-report0
OPT_OMP= -openmp

#option de debugage :
Debug= -traceback -CB -warn alignment -ftrapuv -mp1

IPATH=-I/calcul/DATA_BV/VINCENT/fortran-lib/SOURCE
LPATH=-L/calcul/DATA_BV/VINCENT/fortran-lib/LIB
LIBS= -lfftw3 -mkl -lstat -lspectral -lpod $(LPATH)

CC = ifort -O2 $(IPATH)
CC += $(OPT_para)
CC += $(OPT_OMP)

ALL:$(OBJ)

debug:
	@cd SOURCE/; \
	$(CC) $(Debug) -c lib_spectral.f90  $(LIBS); \
	$(CC) $(Debug) -c lib_pod.f90  $(LIBS); \
	$(CC) $(Debug) -c lib_spectralPOD.f90  $(LIBS); \
	$(CC) $(Debug) -c lib_stat.f90  $(LIBS); \
	$(CC) $(Debug) -c interpol.f90  $(LIBS); \
	$(CC) $(Debug) -c lib_tecplot_io.f90  $(LIBS); \
	$(CC) $(Debug) -c qsort.f90  $(LIBS); \
	$(CC) $(Debug) -c lib_piv_data.f90  $(LIBS);


%.o: %.f90
	@cd SOURCE/; \
	$(CC) -c lib_spectral.f90  $(LIBS); \
	$(CC) -c lib_pod.f90  $(LIBS); \
	$(CC) -c lib_spectralPOD.f90  $(LIBS); \
	$(CC) -c lib_stat.f90  $(LIBS); \
	$(CC) -c interpol.f90  $(LIBS); \
	$(CC) -c lib_tecplot_io.f90  $(LIBS); \
	$(CC) -c qsort.f90  $(LIBS); \
	$(CC) -c lib_piv_data.f90  $(LIBS);

clean:
	rm -rf SOURCE/*.o SOURCE/*~

spectral:
	@cd SOURCE/; \
	$(CC) -c lib_spectral.f90  $(LIBS);

stat:
	@cd SOURCE/; \
	$(CC) -c lib_stat.f90  $(LIBS);

interpol:
	@cd SOURCE/; \
	$(CC) -c interpol.f90  $(LIBS);

tecplot:
	@cd SOURCE/; \
	$(CC) -c lib_tecplot_io.f90  $(LIBS);

pod:
	@cd SOURCE/; \
	$(CC) -c lib_pod.f90  $(LIBS);

spectralpod: spectral
	@cd SOURCE/; \
	$(CC) -c lib_spectralPOD.f90  $(LIBS);

qsort:
	@cd SOURCE/; \
	$(CC) -c qsort.f90  $(LIBS);

piv_data: stat pod qsort
	@cd SOURCE/; \
	$(CC) -c lib_piv_data.f90  $(LIBS);


install:
	ar rc SOURCE/libtecplot.a SOURCE/lib_tecplot_io.o
	ar rc SOURCE/libinterpol.a SOURCE/interpol.o
	ar rc SOURCE/libstat.a SOURCE/lib_stat.o
	ar rc SOURCE/libpod.a SOURCE/lib_pod.o
	ar rc SOURCE/libspectralpod.a SOURCE/lib_spectralPOD.o SOURCE/lib_spectral.o
	ar rc SOURCE/libpivdata.a \
	SOURCE/lib_piv_data.o SOURCE/lib_pod.o SOURCE/lib_stat.o SOURCE/qsort.o
	ar rc SOURCE/libspectral.a SOURCE/lib_spectral.o
	mv SOURCE/*.mod MOD/.
	mv SOURCE/*.a LIB/.
