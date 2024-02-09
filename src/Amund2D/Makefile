LIBDIR=Lib

FCH=gfortran
FC=gfortran
FFLAGS= -O
TRG = demul.e

#-Libraries for linking
ALIBS = $(LIBDIR)/lapackLinux.a $(LIBDIR)/tmglibLinux.a $(LIBDIR)/blasLinux.a

OBJS3D7 = main_demult.o \
cgeco.o\
cgesl.o\
zgeco.o\
zgesl.o\
zgemm.o\
zgemv.o\
zgedi.o\
zswap.o\
matproloog.o\
xerbla.o\
lsame.o

$(TRG): $(OBJS3D7)
	$(FCH) $(FFLAGS) -o $(TRG) $(OBJS3D7)
	
	
.f.o:	
	$(FC) $(FFLAGS) -c $<
	
.c.o:	
	$(CC) $(FFLAGS) -c $<
	
clean:
	-rm *.o 
