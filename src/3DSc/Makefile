####################################################################
# We don't want to hear anything
#.SILENT:
######################################################################
#
# Some global definitions
# - Fortran and C compilers, and C preprocessor 
#FCOMP = g77
FCH=gfortran
FC=gfortran
CC=gcc
%FFLAGS= -C -g
FFLAGS= -O 


# - Include files
INCLUDE = bounds.h param.h grid.h path.h num.h


# - Name of the executable (gets suffixed by '.exe')
# - Default rules

####################################################################
#
# Source programs


OBJS = mainprogram.o\
         Greentens.o\
         green_halfspace.o\
         mksource.o \
         mkmatD.o \
         mkmatDh.o\
         mkmatE.o \
         readpar.o \
         initial.o\
         zgemm.o\
         getpar.o\
         files.o\
         zgemv.o xerbla.o lsame.o matproloog.o calcfield.o

       
####################################################################   
#
# Creation of the executables, linking by means of Fortran-compiler
Plane: $(OBJS)
	$(FCH) $(FFLAGS) -o green.exe $(OBJS) 

.f.o:
	$(FC) $(FFLAGS) -c $<

clean:
	-rm -f *.o green.exe

mainprogram.o : mainprogram.f  bounds.h param.h record.h num.h
Greentens.o : Greentens.f  bounds.h param.h num.h
green_halfspace.o : green_halfspace.f bounds.h param.h num.h embank.h
mksource.o : mksource.f  bounds.h param.h num.h
mkmatD.o : mkmatD.f  bounds.h param.h num.h
mkmatE.o : mkmatE.f  bounds.h param.h num.h
initial.o : initial.f  bounds.h param.h num.h record.h
readpar.o : readpar.f  bounds.h param.h record.h path.h num.h 
matproloog.o : matproloog.f  bounds.h param.h  num.h 
mkmatRuRdTuTd.o :mkmatRuRdTuTd.f  bounds.h param.h num.h 
calcfield.o : calcfield.f  bounds.h param.h num.h 
mkmatDh.o : mkmatDh.f  bounds.h param.h num.h 
initmsg.o : initmsg.c 
files.o : files.c 
getpar.o : getpar.c
