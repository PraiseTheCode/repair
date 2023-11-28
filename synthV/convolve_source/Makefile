#
#Linux with Intel Fortran compiler
#
FF = ifort
FLAG = -cm 
BDIR = .

all: $(BDIR)/convolve Makefile

$(BDIR)/convolve:  main.o disc_integ.o convolution.o  map1.o Makefile
	$(FF) $(FLAG) -o $(BDIR)/convolve *.o -static 

main.o: main.f90 
	$(FF) -c $(FLAG) main.f90

disc_integ.o: disc_integ.f90 
	$(FF) -c $(FLAG) disc_integ.f90

convolution.o: convolution.f90 
	$(FF) -c $(FLAG) convolution.f90

map1.o: map1.f90 
	$(FF) -c $(FLAG) map1.f90
