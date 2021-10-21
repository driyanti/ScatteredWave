LN=mpif77
FC=mpif77
CC=mpicc
CFLAGS= -O3 -Wall
FFLAGS= -x f77-cpp-input -O
LFLAGS= -O
TRG= hlm3d.e

include OBJ.in
### ./phdims.h

$(TRG): $(OBJS3D7) cntrl.h prms3d.h grids.h
	$(LN) $(LFLAGS) -o $(TRG) $(OBJS3D7)

.f.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	-rm -f *.tempt *.o *.e *.loger
