#!/bin/bash
#### Script for compiling and linking test driver programs for the ellipFor library ####

# ---------------- Input for Driver and Compiler Selection --------------------

# Select driver program
DRIVER_FILE="ellipFor_test_driver.f90" # options are "ellipFor_test_driver.f90" and "test_material_driver.f90"

# Select Fortran compiler
COMPILER="gfortran" # options are "gfortran" and "ifx"

# If using ifx compiler, specify the Intel oneAPI directory
oneAPI_dir=/opt/intel/oneapi # Intel oneAPI

# ----------------- No User Input Required Past This Point --------------------

# source files containing algorithm implementation and the driver
SOURCE_FILES="kind_parameters.f90 xelbdj2_all_routines.f90 xgscd_routines.f90 elliptic.f90 ${DRIVER_FILE}"

# determine object file names
OBJECT_FILES=( $(basename -s .f90 $(basename -a $SOURCE_FILES)) )
OBJECT_FILES="${OBJECT_FILES[@]/%/.o}" 

if [ $COMPILER = "gfortran" ] # GNU Fortran
then
  # Optimized for speed and verify that the Fortran 2018 standard is met
  gfortran -flto -O3 -std=f2018 -c $SOURCE_FILES #complile object files
  gfortran -flto -O3 -std=f2018 $OBJECT_FILES -o $(basename -s .f90 ${DRIVER_FILE}) #link object files
elif [ $COMPILER = "ifx" ]    # Intel Fortran
then
  # Optimized for speed	
  source $oneAPI_dir/setvars.sh # initialize environment variables for Intel oneAPI
  ifx -ipo -O3 -c $SOURCE_FILES #complile object files
  ifx -ipo -O3 $OBJECT_FILES -o $(basename -s .f90 ${DRIVER_FILE}) #link object files
else
  echo "Error: the requested Fortran compiler is not available."
fi

##Clean up object files
rm -f *.o *.mod
