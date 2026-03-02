###### Makefile for ellipFor driver programs ######
## Instructions ##
# GNU Fortran: make gfortran
# Intel Fortran: make ifx
### Note for ifx compiler accessed via Intel oneAPI ###
# The appropriate script must be run in the terminal prior to using GNU Make
# e.g., source /opt/intel/oneapi/setvars.sh

# define directory containing ellipFor modules
mod_dir := ellipFor_modules

# GNU Fortran
gfortran: FC := gfortran
gfortran: FOPTS := -flto -O3 -std=f2018
gfortran: ellipFor_test_driver test_material_driver

# GNU Fortran Debug Mode
gfortran_debug: FC := gfortran
gfortran_debug: FOPTS := -Wall -std=f2018 -ffree-line-length-132 -g -fbacktrace -ffpe-trap='invalid,overflow,zero,underflow' -fbounds-check
gfortran_debug: ellipFor_test_driver test_material_driver

# Intel Fortran
ifx: FC := ifx
ifx: FOPTS := -ipo -O3 -fp-model strict
ifx: ellipFor_test_driver test_material_driver

# Note: Fortran module file names must be lowercase in GNU Make rules

test_material_driver: kind_parameters.o xelbdj2_all_routines.o xgscd_routines.o elliptic.o test_material_driver.o
	$(FC) -o test_material_driver $(FOPTS) kind_parameters.o xelbdj2_all_routines.o xgscd_routines.o elliptic.o test_material_driver.o

test_material_driver.o: kind_parameters.mod elliptic.mod
	$(FC) -c $(FOPTS) test_material_driver.f90

ellipFor_test_driver: kind_parameters.o xelbdj2_all_routines.o xgscd_routines.o elliptic.o ellipFor_test_driver.o
	$(FC) -o ellipFor_test_driver $(FOPTS) kind_parameters.o xelbdj2_all_routines.o xgscd_routines.o elliptic.o ellipFor_test_driver.o

ellipFor_test_driver.o: kind_parameters.mod elliptic.mod
	$(FC) -c $(FOPTS) ellipFor_test_driver.f90

elliptic.o elliptic.mod: kind_parameters.mod xelbdj2_all_routines.mod xgscd_routines.mod
	$(FC) -c $(FOPTS) $(mod_dir)/elliptic.f90

xgscd_routines.o xgscd_routines.mod: kind_parameters.mod
	$(FC) -c $(FOPTS) $(mod_dir)/xgscd_routines.f90

xelbdj2_all_routines.o xelbdj2_all_routines.mod: kind_parameters.mod
	$(FC) -c $(FOPTS) $(mod_dir)/xelbdj2_all_routines.f90

kind_parameters.o kind_parameters.mod:
	$(FC) -c $(FOPTS) $(mod_dir)/kind_parameters.f90

clean:
	rm -f *.mod *.o
