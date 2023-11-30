program ellipFor_test_driver
! *** Driver program demonstrating the basic use of ellipFor featuring the evaluation of 1) incomplete Legendre elliptic integrals, 2) complete Legendre elliptic integrals, and 3) the primary Jacobi elliptic functions ***
 use kind_parameters
 use elliptic, only: complete_elliptic_integrals,incomplete_elliptic_integrals,Jacobi_elliptic_functions
implicit none

!!Useful parameters
real(dp), parameter :: pii=3.1415926535897932d0

!!variables for elliptic functions and integrals
real(dp) :: m           !!elliptic parameter
real(dp) :: phi         !!Jacobi amplitude
complex(dp) :: u        !!complex argument
complex(dp) :: Fc,Ec    !!complete elliptic integrals of first and second kinds
complex(dp) :: Fi,Ei    !!incomplete elliptic integrals of first and second kinds 
complex(dp) :: sn,cn,dn !!Jacobi elliptic function values

open(unit=101,file="ellipFor_test_driver.dat") !! file for driver output

!!! 1) complete Legendre elliptic integrals
m=100.d0
call complete_elliptic_integrals(m,Fc,Ec)
write(101,*) "Complete Legendre Elliptic Integrals:"
write(101,*) m,Fc,Ec

!!! 2) incomplete Legendre elliptic integrals
phi=0.75*pii 
m=100.d0     
call incomplete_elliptic_integrals(phi,m,Fi,Ei)
write(101,*) "Incomplete Legendre Elliptic Integrals:"
write(101,*) phi,m,Fi,Ei

!!! 3) Jacobi elliptic functions
u=(1.d0,1.d0) ! complex-valued u
m=100.d0      
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
write(101,*) "Jacobi Elliptic Functions:"
write(101,*) u,m,sn,cn,dn

close(101) !! close output file
end program ellipFor_test_driver

