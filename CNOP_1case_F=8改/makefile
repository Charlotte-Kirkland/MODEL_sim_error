F90   = gfortran
#OBJS  =	Grid.o Adirsol.o Af1model.o Ajcb.o analy.o di2.o dirsol.o costfun.o EvalF.o EvalFG.o EvalG.o F1model.o FFT.o jcb.o l2norm.o Lf1model.o LineSearch2.o Ljcb.o project.o spg2.o weiwo.o weiwob.o cnop-F.o
OBJS  =	module.o ADJ_Lorenz96.o costandgrad.o spg.o TLM_Lorenz96.o Lorenz96.o test.o

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


