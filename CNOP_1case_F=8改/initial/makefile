F90   = ifort
OBJS  =	initial.o

.SUFFIXES: .f90 .o
FC=ifort

all:$(OBJS)
	$(FC) -o initial.exe  $(OBJS) 
.f90.o:
	$(F90) -c  $(NCINC) $<
.f.o:
	$(F90) -c $(NCINC) $<
clean:
	rm *.o *.mod	


