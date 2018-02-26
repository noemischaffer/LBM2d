#FC= /opt/local/bin/gfortran-mp-6
FC= gfortran -fbounds-check
FO=-c -g
default:
	$(FC) $(FO) sub.f90 -o sub.o 
	$(FC) $(FO) cdata.f90 sub.o -o cdata.o 
	$(FC) $(FO) messages.f90 -o messages.o 
	$(FC) $(FO) Avg.f90 cdata.o sub.o messages.o -o Avg.o
	$(FC) $(FO) InitCond.f90 cdata.o messages.o Avg.o -o InitCond.o 
	$(FC) $(FO) BoundCond.f90 cdata.o cdata.o Avg.o -o BoundCond.o
	$(FC) $(FO) Force.f90 cdata.o sub.o messages.o -o Force.o 
	$(FC) $(FO) Evolve.f90 cdata.o sub.o Avg.o messages.o Force.o -o Evolve.o 
	$(FC) $(FO) ShearStress.f90 cdata.o Evolve.o -o ShearStress.o
	$(FC) lbm.f90 InitCond.o BoundCond.o Evolve.o cdata.o sub.o Avg.o messages.o Force.o ShearStress.o -o lbm.x 
clean :
	-rm *.o *.mod *.out 

