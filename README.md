# ellipFor [![DOI](https://zenodo.org/badge/710577422.svg)](https://zenodo.org/doi/10.5281/zenodo.10071175)

Fortran software for the evaluation of Legendre elliptic integrals and Jacobi elliptic functions for generalized input parameters.

## Contents
1. [Background](#background)
2. [Variable Definitions](#variable-definitions)
3. [File Description](#file-description)
    1. [Fortran](#fortran)
    2. [SageMath](#sagemath)
4. [Main ellipFor Subroutines](#main-ellipfor-subroutines)
5. [How to Use the Fortran Routines](#how-to-use-the-fortran-routines)
    1. [Standalone](#standalone)
    2. [With Another Code](#with-another-code)
6. [Legal](#legal)

## Background

This repository contains files and data for the ellipFor library supporting the article "A Fortran software library for Legendre elliptic integrals and Jacobi elliptic
functions with generalized input arguments" by S.J. Trim and R.J. Spiteri. <? Computer algebra scripts for the exact solution are provided in SageMath and Maple. Symbolic computation of the internal heating rate is performed using Maple, which has been translated into Fortran. The Fortran routines can be used to calculate quantities from the exact solution both independently and within an existing convection code. ?>

## Variable Definitions
* $m$ = parameter
* $\phi$ = Jacobi amplitude
* $u$ = first argument of Jacobi elliptic functions
* $K(m)$ = complete Legendre elliptic integral of the first kind
* $E(m)$ = complete Legendre elliptic integral of the second kind
* $F(\phi|m)$ = incomplete Legendre elliptic integral of the first kind
* $E(\phi|m)$ = incomplete Legendre elliptic integral of the second kind
* $\text{sn}(u|m)$ = elliptic sine Jacobi elliptic function
* $\text{cn}(u|m)$ = elliptic cosine Jacobi elliptic function
* $\text{dn}(u|m)$ = delta amplitude Jacobi elliptic function

## File Description

### [Fortran](/Fortran)

[elliptic.f90](/Fortran/elliptic.f90)
* Contains routines that evaluate Legendre elliptic integrals and Jacobi elliptic functions for generalized input parameter ranges
* Input ranges were generalized by combining transformations with calls to routines from [xelbdj2_all_routines.f90](/Fortran/xelbdj2_all_routines.f90) and [xgscd_routines.f90](/Fortran/xgscd_routines.f90)

[xelbdj2_all_routines.f90](/Fortran/xelbdj2_all_routines.f90)
* Contains routines for evaluation of associate incomplete elliptic integrals of first, second, and third kinds
* Assumes standard input parameter ranges
* Adapted from routines by Toshio Fukushima (see [Legal](#legal))

[xgscd_routines.f90](/Fortran/xgscd_routines.f90)
* Contains routines for the evaluation of the Jacobi elliptic functions sn, cn, and dn
* Assumes standard input parameter ranges
* Adapted from routines by Toshio Fukushima (see [Legal](#legal))

[compile_script.sh](/Fortran/compile_script.sh)
* Terminal script for compiling the [standalone](#standalone) version of the Fortran routines using gfortran

[ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90)
* Test driver program for [standalone](#standalone) version of the Fortran routines
* Evaluates $K(m)$, $E(m)$, $F(\phi|m)$, $E(\phi|m)$, $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$ for a given $m$, $\phi$, and $u$

[ellipFor_test_driver](/Fortran/ellipFor_test_driver)
* Sample exexcutable for the [standalone](#standalone) driver program based on [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90)
* Results from running [compile_script](/Fortran/compile_script) in the terminal using the .f90 files in the [Fortran](/Fortran) folder
    * for this sample, the script shell variables `DRIVER_FILE="ellipFor_test_driver.f90"` and `COMPILER="gfortran"`were used
* gfortran 11.4.0 was used

[ellipFor_test_driver.dat](/Fortran/ellipFor_test_driver.dat)
* Output file produced by executing [ellipFor_test_driver](/Fortran/ellipFor_test_driver)
* Contains data for $K(m)$, $E(m)$, $F(\phi|m)$, $E(\phi|m)$, $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$

[test_material_driver.f90](/Fortran/test_material_driver.f90)
* Test material driver program for the Fortran routines 
* Evaluates $K(m)$, $E(m)$, $F(\phi|m)$, $E(\phi|m)$, $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$ for $0 \leq m \leq 100$, $\phi = 3 \pi /4$, and $u=1+i$

[test_material_driver](/Fortran/test_material_driver)
* Sample exexcutable for the test material driver program based on [test_material_driver.f90](/Fortran/test_material_driver.f90)
* Results from running [compile_script](/Fortran/compile_script) in the terminal using the .f90 files in the [Fortran](/Fortran) folder
    * for this sample, the script shell variables `DRIVER_FILE="test_material_driver.f90"` and `COMPILER="gfortran"`were used
* gfortran 11.4.0 was used

[complete_elliptic_integrals.dat](/Fortran/complete_elliptic_integrals.dat)
* Output file produced by executing [test_material_driver](/Fortran/test_material_driver)
* Contains data for $K(m)$ and $E(m)$
* Used in section 6.1 in the article (see [Background](#background))

[incomplete_elliptic_integrals.dat](/Fortran/incomplete_elliptic_integrals.dat)
* Output file produced by executing [test_material_driver](/Fortran/test_material_driver)
* Contains data for $F(\phi|m)$, $E(\phi|m)$
* Used in section 6.2 in the article (see [Background](#background))

[Jacobi_elliptic_functions.dat](/Fortran/Jacobi_elliptic_functions.dat)
* Output file produced by executing [test_material_driver](/Fortran/test_material_driver)
* Contains data for $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$
* Used in section 6.3 in the article (see [Background](#background))

### [SageMath](/SageMath)
[test_problem_solutions.sage](/SageMath/test_problem_solutions.sage)
* SageMath script that can be used to calculate test problem data corresponding to the data produced by [test_material_driver](/Fortran/test_material_driver)
* Output is used in section 6 of the article (see [Background](#background))

[CAS_complete.dat](/SageMath/CAS_complete.dat)
* Output file produced by [test_problem_solutions.sage](/SageMath/test_problem_solutions.sage)
* Contains test problem data for complete Legendre elliptic integrals
* Used in section 6.1 of the article (see [Background](#background))

[CAS_incomplete.dat](/SageMath/CAS_incomplete.dat)
* Output file produced by [test_problem_solutions.sage](/SageMath/test_problem_solutions.sage)
* Contains test problem data for incomplete Legendre elliptic integrals
* Used in section 6.2 of the article (see [Background](#background))

[CAS_functions.dat](/SageMath/CAS_functions.dat)
* Output file produced by [test_problem_solutions.sage](/SageMath/test_problem_solutions.sage)
* Contains test problem data for Jacobi elliptic functions
* Used in section 6.3 of the article (see [Background](#background))

## Main ellipFor Subroutines
`complete_elliptic_integrals(m,Fc,Ec)`
* Evaluate $K(m)$ and $E(m)$
* Input: `m` = $m$ for $m \geq 0$
* Output: `Fc` = $K(m)$ and `Ec` = $E(m)$

`incomplete_elliptic_integrals(phi,m,F,E)`
* Evaluate $F(\phi|m)$ and $E(\phi|m)$
* Input: `phi` = $\phi$ for $\phi \in \mathbb{R}$, and `m` = $m$ for $m \geq 0$
* Output: `F` = $F(\phi|m)$ and `E` = $E(\phi|m)$ 

`Jacobi_elliptic_functions(u,m,sn,cn,dn)`
* Evaluate $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$
* Input: `u` = $u$ for $u \in \mathbb{C}$, and `m` = $m$ for $m \geq 0$
* Output: `sn` = $\text{sn}(u|m)$, `cn` = $\text{cn}(u|m)$, and `dn` = $\text{dn}(u|m)$

Note that the source code for the above routines is in [elliptic.f90](/Fortran/elliptic.f90) and examples for calling the subroutines are in [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90). 

## How to Use the Fortran Routines

### Standalone

Can be used to generate values for $K(m)$, $E(m)$, $F(\phi|m)$, $E(\phi|m)$, $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$ using the driver program [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90).

1. Specify the desired $m$, $\phi$, and $u$ in [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90)
    * Examples are shown in [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90)
2. Compile the code by running [compile_script.sh](/Fortran/compile_script.sh) from the command line.
    * set the driver file shell variable `DRIVER_FILE` in [compile_script.sh](/Fortran/compile_script.sh) as `DRIVER_FILE="ellipFor_test_driver.f90"`
    * set the Fortran compiler shell variable `COMPILER` in [compile_script.sh](/Fortran/compile_script.sh) as `COMPILER="gfortran"` for GNU Fortran or `COMPILER="ifx"` for Intel Fortran
    * Linux/Mac: `$ source compile_script.sh`
        * This produces the executable [ellipFor_test_driver](/Fortran/ellipFor_test_driver)
    * gfortran 11.4.0 or later is recommended
    * for Intel Fortran it is assumed that the ifx compiler is accessed through the Intel oneAPI software package (tested with version 2023.2.1-16 for Linux)
    * Other compilers may be possible but results should be tested
3. Run [ellipFor_test_driver](/Fortran/ellipFor_test_driver) from the command line
    * Linux/Mac: `$ ./ellipFor_test_driver`
    * This will produce data for $K(m)$, $E(m)$, $F(\phi|m)$, $E(\phi|m)$, $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$ in the output file [ellipFor_test_driver.dat](/Fortran/ellipFor_test_driver.dat)
    * Example subroutine calls are shown in [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90)

### With Another Code

Can be used to calculate $K(m)$, $E(m)$, $F(\phi|m)$, $E(\phi|m)$, $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$ from within another code. These instructions presume that the other code is written in Fortran. The following steps are guidelines only. The precise procedure may depend on the particular code used.

#### Fortran

1. Insert calls to the subroutines for Legendre elliptic integrals and Jacobi elliptic functions within the source of the other code where necessary
    * Examples of how to call the subroutines are shown in [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90)  
2. Link the f90 files from the [Fortran](/Fortran) folder named [elliptic.f90](/Fortran/elliptic.f90), [xelbdj2_all_routines.f90](/Fortran/xelbdj2_all_routines.f90), and [xgscd_routines.f90](/Fortran/xgscd_routines.f90) to the source for the other code
    * Example: `gfortran -flto -O3 other_code.f90 elliptic.f90 xelbdj2_all_routines.f90 xgscd_routines.f90 -o other_code`
    * In the above example, the source for the other code is `other_code.f90` and the resulting executable is `other_code`
        * Modify these names as needed
    * gfortran 11.4.0 or later is recommended
    * Other compilers (and compiler options) may be possible but results should be tested
    * Warning: Duplicate variable/routine names may occur
        * Resolve any related compiler errors
        * Verify that the arguments of subroutine calls correspond to the correct values and data types 
3. Run the code execuatable as usual

## Legal

This repository is subject to the GPLv3 [license](/LICENSE).

The routines contained within [xelbdj2_all_routines.f90](/Fortran/xelbdj2_all_routines.f90) and [xgscd_routines.f90](/Fortran/xgscd_routines.f90) are adapted from routines by Toshio Fukushima available under the CC BY-SA 4.0 license. Original versions of these routines can be found at http://dx.doi.org/10.13140/RG.2.2.27011.66085 and https://www.researchgate.net/publication/233903220_xgscdtxt_Fortran_program_package_to_compute_the_Jacobian_elliptic_functions_snum_cnum_dnum.
