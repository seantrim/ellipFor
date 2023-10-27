program exact_solution
 use OMP_LIB
 use kind_parameters
implicit none

!!Useful parameters
real(dp), parameter :: pii=3.1415926535897932d0

!!variables for elliptic functions and integrals
real(dp) :: m            !!elliptic parameter
real(dp) :: phi          !!Jacobi amplitude
complex(dp) :: u        !!complex argument
complex(dp) :: Fc,Ec    !!complete elliptic integrals of first and second kinds
complex(dp) :: Fi,Ei    !!incomplete elliptic integrals of first and second kinds 
complex(dp) :: sn,cn,dn !!Jacobi elliptic function values

integer(isp) :: i !!loop counter
integer(isp) :: N_m
real(dp) :: m_min,m_max,delta_m

!!define range of parameter m for test problems
m_min=0.01d0   !!minimum
m_max=100.01d0  !!maximum
N_m=100       !! # of m values
delta_m=(m_max-m_min)/N_m !!interval between m values
!!end define range of parameter m for test problems


!!!elliptic integrals: loop over m
write(*,*) "Elliptic Integrals: Loop over parameter m."
open(unit=101,file="complete_elliptic_integrals.dat")
open(unit=102,file="incomplete_elliptic_integrals.dat")
phi=0.75*pii!1.6d0
delta_m=(m_max-m_min)/N_m
do i=1,N_m+1
 m=m_min+delta_m*(i-1)
 if (m.ne.1.d0) then 
  call complete_elliptic_integrals(m,Fc,Ec)
  write(101,*) m,Fc%RE,Fc%IM,Ec%RE,Ec%IM
  call incomplete_elliptic_integrals(phi,m,Fi,Ei)
  write(102,*) m,Fi%RE,Fi%IM,Ei%RE,Ei%IM
 end if
end do
close(101)
close(102)
write(*,*) ""
!!!end elliptic integrals: loop over m

!!!!!!!!!!!!!!!!!!!!!!!!!!!Jacobi elliptic functions
write(*,*) "Jacobi Elliptic Functions: Loop over parameter m."
open(unit=103,file="Jacobi_elliptic_functions.dat")
u=(1.d0,1.d0)              ! complex-valued u
delta_m=(m_max-m_min)/N_m
do i=1,N_m+1
 m=m_min+delta_m*(i-1)
 call Jacobi_elliptic_functions(u,m,sn,cn,dn)
 write(103,*) m,sn%RE,sn%IM,cn%RE,cn%IM,dn%RE,dn%IM
end do
close(103)
write(*,*) ""
!!!!!!!!!!!!!!!!!!!!!!!!!!!End Jacobi elliptic functions

end program exact_solution

