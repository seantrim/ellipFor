###### Makefile for ellipFor driver programs ######
## Instructions ##
# GNU Fortran: make gfortran
# Intel Fortran: make ifx
### Note for ifx compiler accessed via Intel oneAPI ###
# The appropriate script must be run in the terminal prior to using GNU Make
# e.g., source /opt/intel/oneapi/setvars.sh

# GNU Fortran
gfortran: FC = gfortran
gfortran: FOPTS = -flto -O3 -std=f2018
gfortran: ellipFor_test_driver test_material_driver
gfortran:
	echo $(FC) $(FOPTS)

# GNU Fortran Debug Mode
gfortran_debug: FC = gfortran
gfortran_debug: FOPTS = -Wall -std=f2018 -ffree-line-length-132 -g -fbacktrace -ffpe-trap='invalid,overflow,zero,underflow' -fbounds-check
gfortran_debug: ellipFor_test_driver test_material_driver
gfortran_debug:
	echo $(FC) $(FOPTS)

# Intel Fortran
ifx: FC = ifx
ifx: FOPTS = -ipo -O3
ifx: ellipFor_test_driver test_material_driver
ifx:
	echo $(FC) $(FOPTS)

test_material_driver: kind_parameters.o xelbdj2_all_routines.o xgscd_routines.o elliptic.o elliptic.mod test_material_driver.o kind_parameters.mod elliptic.mod xelbdj2_all_routines.mod xgscd_routines.mod
	$(FC) -o test_material_driver $(FOPTS) test_material_driver.f90 kind_parameters.f90 xelbdj2_all_routines.f90 xgscd_routines.f90 elliptic.f90

test_material_driver.o: test_material_driver.f90 kind_parameters.mod elliptic.mod
	$(FC) -c $(FOPTS) test_material_driver.f90

ellipFor_test_driver: kind_parameters.o xelbdj2_all_routines.o xgscd_routines.o elliptic.o elliptic.mod ellipFor_test_driver.o kind_parameters.mod elliptic.mod xelbdj2_all_routines.mod xgscd_routines.mod
	$(FC) -o ellipFor_test_driver $(FOPTS) ellipFor_test_driver.f90 kind_parameters.f90 xelbdj2_all_routines.f90 xgscd_routines.f90 elliptic.f90

ellipFor_test_driver.o: ellipFor_test_driver.f90 kind_parameters.mod elliptic.mod
	$(FC) -c $(FOPTS) ellipFor_test_driver.f90

elliptic.o elliptic.mod: elliptic.f90 kind_parameters.mod xelbdj2_all_routines.mod xgscd_routines.mod
	$(FC) -c $(FOPTS) elliptic.f90

xgscd_routines.o xgscd_routines.mod: xgscd_routines.f90 kind_parameters.mod
	$(FC) -c $(FOPTS) xgscd_routines.f90

xelbdj2_all_routines.o xelbdj2_all_routines.mod: xelbdj2_all_routines.f90 kind_parameters.mod
	$(FC) -c $(FOPTS) xelbdj2_all_routines.f90

kind_parameters.o kind_parameters.mod: kind_parameters.f90
	$(FC) -c $(FOPTS) kind_parameters.f90

clean:
	rm -f *.mod *.o
