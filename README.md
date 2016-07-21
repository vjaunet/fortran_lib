Fortran Librairies
------------------

Installation :
make;

make install;

In the program make file include the installation path "IPATH" as the following:

ifort -IIPATH/MOD

To call a librairy, compile with the following options :

ifort -IIPATH/MOD -c myprog.f90 -o exec -LIPATH/LIB -lmylib

and "use module_name" has to be place before the implicit none of the main
program file.

mylib has to be raplaced with one the following librairy names :
- tecplot
  a tecplot_io module for read and writing acsii tecplot files

- interp
  a interpolation module containing 1d and 2 routines

- sepctral
  psd, coherence, cross-psd, auto/cross-correlation of real*8 and complex*16 varialbes. And more..

- even more to be described

TODO : put some examples, do a real wiki
