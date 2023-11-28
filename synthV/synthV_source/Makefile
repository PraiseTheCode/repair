#
#Linux with Intel Fortran compiler
#
FF = ifort
FLAG = -c
BDIR = .

all: $(BDIR)/SynthV_26 Makefile

$(BDIR)/SynthV_26:  data_modules.o modules.o arrays_allocation.o calculate_spectrum.o continuum_opacity.o effective_depths.o Lemkesdat.o\
       lines_opacity.o first_record.o HeI_prof.o Hydrogen.o Josh.o kapp.o load_model.o Math_subroutines.o \
       number_density.o synthV_input.o SynthV_main.o Makefile
	$(FF)  -o $(BDIR)/SynthV_26 *.o

data_modules.o: data_modules.f90 
	$(FF)  $(FLAG) data_modules.f90

modules.o: modules.f90 
	$(FF)  $(FLAG) modules.f90

arrays_allocation.o: arrays_allocation.f90 
	$(FF)  $(FLAG) arrays_allocation.f90

calculate_spectrum.o: calculate_spectrum.f90 
	$(FF)  $(FLAG) calculate_spectrum.f90

continuum_opacity.o: continuum_opacity.f90 
	$(FF)  $(FLAG) continuum_opacity.f90

lines_opacity.o: lines_opacity.f90 
	$(FF)  $(FLAG) lines_opacity.f90

effective_depths.o: effective_depths.f90 
	$(FF)  $(FLAG) effective_depths.f90

first_record.o: first_record.f90 
	$(FF)  $(FLAG) first_record.f90

Lemkesdat.o: Lemkesdat.f90 
	$(FF)  $(FLAG) Lemkesdat.f90

Hydrogen.o: Hydrogen.f90  
	$(FF)  $(FLAG) Hydrogen.f90

HeI_prof.o: HeI_prof.f90 
	$(FF)  $(FLAG) HeI_prof.f90

load_model.o: load_model.f90 
	$(FF)  $(FLAG) load_model.f90

Josh.o: Josh.f90 
	$(FF)  $(FLAG) Josh.f90

kapp.o: kapp.f90
	$(FF)  $(FLAG) kapp.f90

Math_subroutines.o: Math_subroutines.f90 
	$(FF)  $(FLAG) Math_subroutines.f90

number_density.o: number_density.f90 
	$(FF)  $(FLAG) number_density.f90

synthV_input.o: synthV_input.f90 
	$(FF)  $(FLAG) synthV_input.f90

SynthV_main.o: SynthV_main.f90 
	$(FF)  $(FLAG) SynthV_main.f90


clean:
	@rm -f *.o
	@rm -f *.mod
	@rm -f SynthV_26
