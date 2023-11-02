# Select driver program
DRIVER_FILE="ellipFor_test_driver.f90" #options are "ellipFor_test_driver.f90" "test_material_driver.f90"

# ----------------- No User Input Required Past This Point --------------------

# source files containing algorithm implementation and the driver
SOURCE_FILES="elliptic.f90 xelbdj2_all_routines.f90 xgscd_routines.f90 ${DRIVER_FILE}"

# determine object file names
OBJECT_FILES=( $(basename -s .f90 $(basename -a $SOURCE_FILES)) )
OBJECT_FILES="${OBJECT_FILES[@]/%/.o}" 

#Optimized for speed
gfortran -flto -O3 -std=f2018 -c $SOURCE_FILES #complile object files
gfortran -flto -O3 -std=f2018 $OBJECT_FILES -o $(basename -s .f90 ${DRIVER_FILE}) #link object files

##Clean up object files
rm -f *.o *.mod
