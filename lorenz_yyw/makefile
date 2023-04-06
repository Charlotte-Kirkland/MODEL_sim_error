F90   = gfortran
OBJS  =	Lorenz63.o test.o

.SUFFIXES: .f90 .o
FC=gfortran

all:$(OBJS)
	$(FC) -o test.exe  $(OBJS) 
.f90.o:
	$(F90) -c  $(NCINC) $<
.f.o:
	$(F90) -c $(NCINC) $<
clean:
	rm *.o *.mod	


