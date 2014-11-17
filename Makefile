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

IPATH=-I/home/vjaunet/VJT/FORTRAN-LIB/SOURCE
LPATH=-L/home/vjaunet/VJT/FORTRAN-LIB/LIB
LIBS= -lfftw3 -llapack -lstat -lpod $(LPATH)

CC = ifort -O2 $(IPATH)
CC += $(OPT_para)
CC += $(OPT_OMP)
#CC += $(Debug)

ALL:$(OBJ)

%.o: %.f90
	@cd SOURCE/; \
	$(CC) -c lib_spectral.f90  $(LIBS); \
	$(CC) -c lib_pod.f90  $(LIBS); \
	$(CC) -c lib_stat.f90  $(LIBS); \
	$(CC) -c interpol.f90  $(LIBS); \
	$(CC) -c tecplot_io.f90  $(LIBS); \
	$(CC) $(IPATH) -c lib_piv_data.f90  $(LIBS);

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
	$(CC) -c tecplot_io.f90  $(LIBS);

pod:
	@cd SOURCE/; \
	$(CC) -c lib_pod.f90  $(LIBS);


piv_data: stat pod
	@cd SOURCE/; \
	$(CC) -c lib_piv_data.f90  $(LIBS);


install:
	ar rc SOURCE/libtecplot.a SOURCE/tecplot_io.o
	ar rc SOURCE/libinterpol.a SOURCE/interpol.o
	ar rc SOURCE/libstat.a SOURCE/lib_stat.o
	ar rc SOURCE/libpod.a SOURCE/lib_pod.o
	ar rc SOURCE/libpivdata.a \
	SOURCE/lib_piv_data.o SOURCE/lib_pod.o SOURCE/lib_stat.o
	ar rc SOURCE/libspectral.a SOURCE/lib_spectral.o
	mv SOURCE/*.mod MOD/.
	mv SOURCE/*.a LIB/.
