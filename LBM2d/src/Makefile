FC= /opt/local/bin/gfortran-mp-6
#FC= gfortran -fbounds-check
FO=-c -g
default:
	$(FC) $(FO) sub.f90 -o sub.o 
	$(FC) $(FO) cdata.f90 -o cdata.o
	$(FC) $(FO) messages.f90 -o messages.o
	$(FC) $(FO) Avg.f90 -o Avg.o
	$(FC) $(FO) BoundCond.f90  -o BoundCond.o
	$(FC) $(FO) io.f90 -o io.o
	$(FC) $(FO) InitCond.f90 -o InitCond.o 
	$(FC) $(FO) diagnostic.f90 -o diagnostic.o 
	$(FC) $(FO) Force.f90  -o Force.o 
	$(FC) $(FO) Evolve.f90  -o Evolve.o 
	$(FC) $(FO) ShearStress.f90  -o ShearStress.o
	$(FC) lbm.f90 io.o InitCond.o BoundCond.o diagnostic.o Evolve.o cdata.o sub.o Avg.o messages.o Force.o ShearStress.o -o lbm.x 
clean :
	-rm *.o *.mod *.out 

