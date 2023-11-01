# ellipFor <? [![DOI](https://zenodo.org/badge/534176632.svg)](https://zenodo.org/badge/latestdoi/534176632) ?>

Fortran software for the evaluation of Legendre elliptic integrals and Jacobi elliptic functions for generalized input parameters.

## Contents
1. [Background](#background)
2. [Variable Definitions](#variable-definitions)
3. [File Description](#file-description)
    1. [Fortran](#fortran)
    2. [SageMath](#sagemath)
    3. [Data](#data)
4. [How to Use the Fortran Routines](#how-to-use-the-fortran-routines)
    1. [Standalone](#standalone)
    2. [With Another Code](#with-another-code)
5. [Legal](#legal)

## Background

This repository contains files and data supporting the article "A Fortran software library for Legendre elliptic integrals and Jacobi elliptic
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

<?
### [SageMath](/SageMath)
[exact_solution.sage](/SageMath/exact_solution.sage)
* SageMath script that can be used to:
    * symbolically compute the formulas for $T$ and $v_{RMS}$
    * generate data and plots for $C$ and $T$
    * note: does not calculate $H$ (see [Maple](/Maple))

### [Maple](/Maple)
[maple_analytic_solution_include_exterior.mw](/Maple/maple_analytic_solution_include_exterior.mw)
* Maple worksheet for symbolic computation of $H$ (everywhere except at the domain boundaries)
* Results stored in [foo_exterior.m](/Maple/foo_exterior.m) 

[foo_exterior.m](/Maple/foo_exterior.m)
* Results from running [maple_analytic_solution_include_exterior.mw](/Maple/maple_analytic_solution_include_exterior.mw)
* Once loaded into a Maple worksheet, functions such as $C(x,z,t)$, $T(x,z,t)$, and $H(x,z,t)$ become available
    * Saves time compared to rerunning [maple_analytic_solution_include_exterior.mw](/Maple/maple_analytic_solution_include_exterior.mw) 

[Fortran_code_generation.mw](/Maple/Fortran_code_generation.mw)
* Maple worksheet for translating Maple's $H$ formula into Fortran 77 code
* The result was adapted to free form Fortran for use in [H_func.f90](/Fortran/H_func.f90)

### [Fortran](/Fortran)
[exact_solution_main.f90](/Fortran/exact_solution_main.f90)
* Main program for [standalone](#standalone) version of the Fortran routines

[exact_solution_routines.f90](/Fortran/exact_solution_routines.f90)
* Contains routines for calculating physical quantities including:
    * $C$, $T$, $H$, $v_{RMS}$, and $E$

[H_func.f90](/Fortran/H_func.f90)
* Contains routine for evaluating $H(x,z,t)$

[H_helper_routines.f90](/Fortran/H_helper_routines.f90)
* Contains Fortran equivalents of Maple functions referenced in [H_func.f90](/Fortran/H_func.f90)

[input_functions.f90](/Fortran/input_functions.f90)
* Contains user defined function $f(t)$, with its integral and derivatives

[elliptic.f90](/Fortran/elliptic.f90)
* Contains routines that evaluate elliptic integrals and Jacobi elliptic functions for the input parameter ranges needed
* Input ranges were generalized by combining identities with calls to routines from [xelbdj2_all_routines.f90](/Fortran/xelbdj2_all_routines.f90) and [xgscd_routines.f90](/Fortran/xgscd_routines.f90)

[xelbdj2_all_routines.f90](/Fortran/xelbdj2_all_routines.f90)
* Contains routines for evaluation of associate incomplete elliptic integrals of first, second, and third kinds
* Assumes standard input parameter ranges
* Adapted from routines by Toshio Fukushima (see [Legal](#legal))

[xgscd_routines.f90](/Fortran/xgscd_routines.f90)
* Contains routines for the evaluation of the Jacobi elliptic functions sn, cn, and dn
* Assumes standard input parameter ranges
* Adapted from routines by Toshio Fukushima (see [Legal](#legal))

[compile_script](/Fortran/compile_script)
* Terminal script for compiling the [standalone](#standalone) version of the Fortran routines using gfortran

[exact_solution_code](/Fortran/exact_solution_code)
* Sample exexcutable for the [standalone](#standalone) version of the Fortran routines
* Results from running [compile_script](/Fortran/compile_script) in the terminal using the .f90 files in the [Fortran](/Fortran) folder
* gfortran 11.3.0 was used  

### [Data](/Data)
[entrainment_sample_1_401x401.dat](/Data/entrainment_sample_1_401x401.dat)
* $E$ time series data for temporally periodic case in "Sample Results" section 
* Computed using [Fortran](/Fortran) routines

[entrainment_sample_2_751x501.dat](/Data/entrainment_sample_2_751x501.dat)
* $E$ time series data for approaching steady state case in "Sample Results" section
* Computed using [Fortran](/Fortran) routines
?>

## How to Use the Fortran Routines

### Standalone

Can be used to generate values for $K(m)$, $E(m)$, $F(\phi|m)$, $E(\phi|m)$, $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$ using the driver program [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90).

1. Specify the desired $m$, $\phi$, and $u$ in [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90)
    * Examples are shown in [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90)
2. Compile the code by running [compile_script.sh](/Fortran/compile_script.sh) from the command line.
    * set the driver file shell variable in the first line of [compile_script.sh](/Fortran/compile_script.sh) as `DRIVER_FILE="ellipFor_test_driver.f90"`
    * Linux: `$ source compile_script.sh`
        * This produces the executable [ellipFor_test_driver](/Fortran/ellipFor_test_driver)
    * gfortran 11.3.0 or later is recommended
    * Other compilers may be possible but results should be tested
3. Run [ellipFor_test_driver](/Fortran/ellipFor_test_driver) from the command line
    * Linux: `$ ./ellipFor_test_driver`
    * This will produce data for $K(m)$, $E(m)$, $F(\phi|m)$, $E(\phi|m)$, $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$ in the output file [ellipFor_test_driver.dat](/Fortran/ellipFor_test_driver.dat)
    * Example subroutine calls are shown in [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90)

### With Another Code

Can be used to calculate $K(m)$, $E(m)$, $F(\phi|m)$, $E(\phi|m)$, $\text{sn}(u|m)$, $\text{cn}(u|m)$, and $\text{dn}(u|m)$ from within another code. These instructions presume that the other code is written in Fortran. The following steps are guidelines only. The precise procedure may depend on the particular code used.

#### Fortran

1. Insert calls to the subroutines for Legendre elliptic integrals and Jacobi elliptic functions within the source of the other code where necessary
    * Examples of how to call the subroutines are shown in [ellipFor_test_driver.f90](/Fortran/ellipFor_test_driver.f90)
        * The H value is returned in the rightmost argument  
2. Link the f90 files from the [Fortran](/Fortran) folder named [elliptic.f90](/Fortran/elliptic.f90), [xelbdj2_all_routines.f90](/Fortran/xelbdj2_all_routines.f90), and [xgscd_routines.f90](/Fortran/xgscd_routines.f90) to the source for the other code
    * Example: `gfortran -flto -O3 other_code.f90 elliptic.f90 xelbdj2_all_routines.f90 xgscd_routines.f90 -o other_code`
    * In the above example, the source for the other code is `other_code.f90` and the resulting executable is `other_code`
        * Modify these names as needed
    * gfortran 11.3.0 or later is recommended
    * Other compilers (and compiler options) may be possible but results should be tested
    * Warning: Duplicate variable/routine names may occur
        * Resolve any related compiler errors
        * Verify that the arguments of subroutine calls correspond to the correct values and data types 
3. Run the code execuatable as usual

## Legal

This repository is subject to the GPLv3 [license](/LICENSE).

The routines contained within [xelbdj2_all_routines.f90](/Fortran/xelbdj2_all_routines.f90) and [xgscd_routines.f90](/Fortran/xgscd_routines.f90) are adapted from routines by Toshio Fukushima available under the CC BY-SA 4.0 license. Original versions of these routines can be found at http://dx.doi.org/10.13140/RG.2.2.27011.66085 and https://www.researchgate.net/publication/233903220_xgscdtxt_Fortran_program_package_to_compute_the_Jacobian_elliptic_functions_snum_cnum_dnum.
