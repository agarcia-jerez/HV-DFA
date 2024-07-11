
HV-DFA
(c) 2016 Antonio Garcia Jerez, José Piña Flores

This program is mainly written in FORTRAN 90. Nevertheless some FORTRAN 2003 procedures are used to read the command line.
Alternative subroutines for command line reading are provided for compatibility with some FORTRAN 90/f95 compilers, 
which have to be manually enabled in HV.f90

The code has been checked under windows with:
gfortran compiler (recommended)
ifort
ftn95 v. 7.1
 
Digital Visual Fortran Optimizing Compiler v. 5.0 (Microsoft Developer Studio 5.0). Requires choosing suitable subroutines in HV.f90 (MISCELLANEA section)

The program takes advantage of basic parallel computation features (via OpenMP).


COMPILATION (examples)

gfortran HV.f90 -o HVf.exe -static -Ofast -fopenmp     (with gcc 5.2.0 from www.equation.com). Note that -static may be incompatible with -fopenmp in official gcc releases

ifort HV.f90 -o iHVf.exe -static -fast -openmp -Qparallel (Intel Fortran)

ftn95 HV.f90 /OPTIMISE /LINK                   (with ftn95 7.1 at http://silverfrost.com)


USAGE(examples)

HV.exe -f model.txt -nf 100 -fmin 0.1 -fmax 10 -logsam -nmr 20 -nml 20 -ph -hv  > HV.dat

Type HV.exe -h for help.

 