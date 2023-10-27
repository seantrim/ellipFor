module kind_parameters
 use iso_fortran_env, only : int32,int64 
implicit none
integer, parameter :: sp = kind(0.0e0)
integer, parameter :: dp = kind(0.0d0)
integer, parameter :: qp = kind(0.0q0)
integer, parameter :: isp = int32
integer, parameter :: idp = int64
end module kind_parameters

subroutine Jacobi_elliptic_functions(u,m,sn,cn,dn)
!!computes the Jacobi elliptic functions sn, cn, and dn
!!argument u can be any complex value
!!real elliptic parameter must satisfy 0<=m (m>1 is allowed)
 use kind_parameters
implicit none

!!input
complex(dp),intent(in) :: u !!argument 
real(dp),intent(in) :: m !!elliptic parameter

!!output
complex(dp),intent(out) :: sn,cn,dn !!Jacobi elliptic functions

!!internal variables
real(dp) :: mr !!reciprocal parameter
real(dp) :: k  !!modulus
complex(dp) :: u_temp !!temporary argument
complex(dp) :: sn_temp,cn_temp,dn_temp !!temporary Jacobi elliptic function values

if (m.gt.1.d0) then
 k=sqrt(m); mr=1.d0/m
 u_temp=k*u
 !write(*,*) "JEF u,m,k,u_temp,mr=",u,m,k,u_temp,mr
 call Jacobi_elliptic_functions_complex_argument_standard_parameter(u_temp,real(mr,qp),sn_temp,cn_temp,dn_temp)
 sn=sn_temp/k; dn=cn_temp; cn=dn_temp
else
 call Jacobi_elliptic_functions_complex_argument_standard_parameter(u,real(m,qp),sn,cn,dn)
end if
end subroutine Jacobi_elliptic_functions

subroutine Jacobi_elliptic_functions_complex_argument_standard_parameter(u,m,sn,cn,dn)
!!u can be any complex value
!!assumes 0<=m<=1
 use kind_parameters
implicit none

!!input
complex(dp),intent(in) :: u !!argument 
real(qp),intent(in) :: m !!elliptic parameter

!!output
complex(dp),intent(out) :: sn,cn,dn !!Jacobi elliptic functions

!!internal variables
real(dp) :: u_temp !!temporary argument
real(dp) :: u_temp_c !!temporary argument associated with complimentary parameter
real(dp) :: sn_temp,cn_temp,dn_temp !!temporary Jacobi elliptic
real(dp) :: sn_temp_c,cn_temp_c,dn_temp_c !!temporary Jacobi elliptic using complimentary parameters
real(qp) :: mc !!compliment of parameter
real(dp) :: delta !!divisor
real(dp) :: sn_temp_dn_temp_c,cn_temp_cn_temp_c,dn_temp_cn_temp_c !!products of temporary variables
real(dp) :: k
real(dp) :: m_dp !!double precision value of m

if (u%IM.ne.0.d0) then !!complex argument
 m_dp=real(m,dp)
 k=sqrt(m_dp)
 mc=1.q0-m
 !mc=(1.d0+k)*(1.d0-k)
 u_temp_c=u%IM
 !write(*,*) "JEFCargSParam: u,m=",u,m
 !write(*,*) "JEFCargSParam: u_temp_c,mc=",u_temp_c,mc
 call Jacobi_elliptic_functions_standard_input_range(u_temp_c,mc,sn_temp_c,cn_temp_c,dn_temp_c)
 !write(*,*) "JEFCargSParam: sn_temp_c,cn_temp_c,dn_temp_c=",sn_temp_c,cn_temp_c,dn_temp_c
 !write(*,*) ""
 u_temp=u%RE
 !write(*,*) "JEFCargSParam: u_temp,m=",u_temp,m
 call Jacobi_elliptic_functions_standard_input_range(u_temp,m,sn_temp,cn_temp,dn_temp)
 !write(*,*) "JEFCargSParam: sn_temp,cn_temp,dn_temp=",sn_temp,cn_temp,dn_temp
 delta=cn_temp_c**2+m_dp*(sn_temp*sn_temp_c)**2 !!reduces truncation error (this form avoids subtraction)

 sn_temp_dn_temp_c=sn_temp*dn_temp_c; cn_temp_cn_temp_c=cn_temp*cn_temp_c; dn_temp_cn_temp_c=dn_temp*cn_temp_c

 sn=complex(sn_temp_dn_temp_c,sn_temp_c*cn_temp_cn_temp_c*dn_temp)/delta
 cn=complex(cn_temp_cn_temp_c,-sn_temp_dn_temp_c*dn_temp*sn_temp_c)/delta
 dn=complex(dn_temp_cn_temp_c*dn_temp_c,-m_dp*sn_temp*cn_temp*sn_temp_c)/delta

else !!real argument
 u_temp=u%RE
 call Jacobi_elliptic_functions_standard_input_range(u_temp,m,sn_temp,cn_temp,dn_temp)
 sn=complex(sn_temp,0.d0); cn=complex(cn_temp,0.d0); dn=complex(dn_temp,0.d0)
end if
end subroutine Jacobi_elliptic_functions_complex_argument_standard_parameter

subroutine Jacobi_elliptic_functions_standard_input_range(u,m,sn,cn,dn)
!!computes the Jacobi elliptic functions sn, cn, and dn
!!argument u can be any real value
!!elliptic parameter must be within 0<=m<=1  **** m is quadruple precision ****
  use kind_parameters
implicit none

!!input
real(dp),intent(in) :: u !!argument and elliptic parameter
real(qp),intent(in) :: m

!!output
real(dp),intent(out) :: sn,cn,dn !!Jacobi elliptic functions

!!internal variables
real(qp) :: mc !!compliment of elliptic parameter
!real(dp) :: k !!elliptic modulus

!k=sqrt(m)
!mc=(1.d0+k)*(1.d0-k)
mc=1.q0-m
call gscd(u,mc,sn,cn,dn)
end subroutine Jacobi_elliptic_functions_standard_input_range

subroutine complete_elliptic_integrals(m,Fc,Ec)
!!!!SJT: computes the complete elliptic integrals of first and second kind given parameter m
!!!!SJT: assumes 0<=m
!!!!SJT: allows m>1 by applying the reciprocal-modulus transformation for complete elliptic integrals
 use kind_parameters
implicit none

!!input
real(dp),intent(in) :: m !!elliptic characteristic and parameter

!!output
complex(dp),intent(out) :: Fc,Ec !!complete elliptic integrals of the first and  second kinds

!!internal variables
real(dp) :: mc !!compliment of elliptic parameter
real(dp) :: mr,mrc !!reciprocal parameter and its compliment 
real(dp) :: k !!elliptic modulus and reciprocal
real(dp) :: n !!characteristic
real(dp) :: Fc_temp,Ec_temp,Pc_temp !!temporary complete elliptic integral values from standard input range
real(dp) :: Fc_r,Ec_r,Pc_r,Fc_rc,Ec_rc,Pc_rc !!complete elliptic integral values based on reciprocal parameter and its compliment
real(qp) :: m_qp,mc_qp,mr_qp,mrc_qp !!quad precision variables

n=0.d0 !!arbitrary characteristic value -- integrals of the third kind are not used

if (m.gt.1.d0) then
 k=sqrt(m)
 !mc=1.d0-m; mr=1.d0/m; mrc=1.d0-mr
 m_qp=real(m,qp); mr_qp=real(1.d0/m,qp); mrc_qp=1.q0-mr_qp; mc_qp=1.q0-m_qp; mc=real(mc_qp,dp)
 call complete_elliptic_integrals_standard_input_range(n,mr_qp,Fc_r,Ec_r,Pc_r)
 call complete_elliptic_integrals_standard_input_range(n,mrc_qp,Fc_rc,Ec_rc,Pc_rc)
 Fc=complex(Fc_r,-Fc_rc)/k; Ec=complex(m*Ec_r+mc*Fc_r,-Fc_rc+m*Ec_rc)/k !!OG -- cancellation error possible for real part of Ec
else
 call complete_elliptic_integrals_standard_input_range(n,real(m,qp),Fc_temp,Ec_temp,Pc_temp)
 Fc=complex(Fc_temp,0.d0); Ec=complex(Ec_temp,0.d0)
end if
end subroutine complete_elliptic_integrals

subroutine complete_elliptic_integrals_standard_input_range(n,m,Fc,Ec,Pc)
!!!!SJT: computes the complete elliptic integrals of kinds first to third given parameter m and characteristic n
!!!!SJT: assumes 0<=m<=1 and 0<=n<=1
 use kind_parameters
implicit none

!!input
real(dp),intent(in) :: n !!elliptic characteristic and parameter
real(qp),intent(in) :: m

!!output
real(dp),intent(out) :: Fc,Ec,Pc !!complete elliptic integrals of the first kind, second kind, and third kind

!!internal variables
real(dp), parameter :: pii=3.1415926535897932d0
real(dp) :: zero,piio2 !!constants
real(dp) :: bc,dc,jc !!associate complete elliptic integrals
real(qp) :: mc !!complimentary parameter
!real(dp) :: k !!elliptic modulus

zero=0.d0
piio2=pii/2.d0

!mc=1.d0-m !!complimentary parameter
mc=1.q0-m
!k=sqrt(m); mc=(1.d0+k)*(1.d0-k)
call elbdj2(piio2,zero,n,mc,bc,dc,jc) !!returns the associate complete elliptic integrals
Fc=bc+dc; Ec=bc+real(mc,dp)*dc; Pc=Fc+n*jc !!build standard elliptic integrals from associate elliptic integrals
end subroutine complete_elliptic_integrals_standard_input_range

subroutine incomplete_elliptic_integrals_standard_input_range(phi,n,m,F,E,P)
!!!!SJT: computes the incomplete elliptic integrals of kinds first to third given parameter m and characteristic n
!!assumes 0<=phi<=pi/2, 0<=m<=1, 0<=n<=1
 use kind_parameters
implicit none

!!input
real(dp),intent(in) :: phi,n !!elliptic amplitude, characteristic, and parameter
real(qp),intent(in) :: m

!!output
real(dp),intent(out) :: F,E,P !!incomplete elliptic integrals of the first kind, second kind, and third kind

!!internal variables
real(dp), parameter :: pii=3.1415926535897932d0
real(dp) :: piio2      !!constants
real(dp) :: phic    !!complimentary variables
real(qp) :: mc
!real(dp) :: k          !!elliptic modulus
real(dp) :: b,d,j      !!associate incomplete elliptic integrals

piio2=pii/2.d0

phic=piio2-phi;
!k=sqrt(m); mc=(1.d0+k)*(1.d0-k)
!mc=1.d0-m
mc=1.q0-m
call elbdj2(phi,phic,n,mc,b,d,j) !!associate incomplete elliptic integrals
F=b+d; E=b+real(mc,dp)*d; P=F+n*j !!build standard elliptic integrals from associate elliptic integrals
end subroutine incomplete_elliptic_integrals_standard_input_range

subroutine incomplete_elliptic_integrals_standard_amp_large_parameter(phi,m,F,E)
!!!!SJT: reduce parameter using reciprocal-modulus transformation (if needed)
!!!!SJT: computes the incomplete elliptic integrals of first and second kind given amplitude phi and parameter m
!!assumes 0<=phi<=pi/2 and 0<=m (m>1 is allowed)
 use kind_parameters
implicit none

!!input
real(dp),intent(in) :: phi,m !!elliptic amplitude and parameter

!!output
complex(dp),intent(out) :: F,E !!incomplete elliptic integrals of the first and second kind

!!internal variables
real(dp), parameter :: pii=3.1415926535897932d0
real(qp), parameter :: pii_qp=3.14159265358979323846264338327950288q0
real(dp) :: piio2 !!constants
real(dp) :: phic,mc    !!complimentary variables
real(dp) :: phi_temp,phic_temp,b_temp,d_temp,j_temp !!temporary variables
real(dp) :: u !!argument of arcsine
real(dp) :: k !!elliptic modulus
real(dp) :: kr,mr !!reciprocal modulus and parameter
real(dp) :: mrc !!compliment of the reciprocal parameter
complex(dp) :: Fc,Ec !!complete elliptic integral values computed for m>1
real(dp) :: Fc_temp,Ec_temp,Pc_temp !!temporary complete elliptic integral values
real(dp) :: F_temp,E_temp,P_temp !!temporary incomplete elliptic integral values
real(dp) :: amp,ampc !!temporary elliptic amplitude variable and compliment
real(dp) :: sin_amp !!sine amplitude
real(dp) :: n !!characteristic
real(qp) :: mr_qp,mrc_qp,m_qp,mc_qp
real(qp) :: u_qp,phi_temp_qp,amp_qp,sin_amp_qp

piio2=pii/2.d0
n=0.d0 !!arbitrary characteristic -- integrals of the third kind are not used

if (phi.eq.piio2) then
 call complete_elliptic_integrals(m,F,E)
else !!0<=phi<pi/2
 if (m.gt.1.d0) then !!large parameter m>1
  k=sqrt(m)
  !mc=1.d0-m; mr=1.d0/m !!OG
  mc=real(1.q0-real(m,qp),dp); mr_qp=1.q0/real(m,qp)
  !u=k*sin(phi) !!argument of arcsine
  u_qp=real(k,qp)*sin(real(phi,qp))
  if (u_qp.le.1.q0) then
   !amp=asin(u)
   amp=real(asin(u_qp),dp)
   call incomplete_elliptic_integrals_standard_input_range(amp,n,mr_qp,F_temp,E_temp,P_temp)
   F=complex(F_temp/k,0.d0); E=complex(k*E_temp+mc*F%RE,0.d0)
  else
   !mrc=1.d0-mr; !!OG
   mrc_qp=1.q0-mr_qp; mrc=real(mrc_qp,dp)
   !sin_amp=(sqrt(u**2.d0-1.d0))/(u*sqrt(mrc))
   sin_amp_qp=(sqrt(u_qp**2-1.q0))/(u_qp*sqrt(mrc_qp))
   !amp=asin(min(sin_amp,1.d0)) !!it is possible for sin_amp to be slightly greater than unity due to round off error
   amp_qp=asin(min(sin_amp_qp,1.q0)); amp=real(amp_qp,dp) !!it may be possible for sin_amp to be slightly greater than unity due to round off error
   call complete_elliptic_integrals_standard_input_range(n,mr_qp,Fc_temp,Ec_temp,Pc_temp)
   call incomplete_elliptic_integrals_standard_input_range(amp,n,mrc_qp,F_temp,E_temp,P_temp)
   F=complex(Fc_temp,-F_temp)/k
   E=complex(m*Ec_temp+mc*Fc_temp,&
     &-F_temp+m*E_temp+real(((1.q0-m)*sin_amp_qp*cos(amp_qp))/sqrt(1.q0-mrc_qp*sin_amp_qp**2),dp))/k
  end if
 else !!standard input ranges
  m_qp=real(m,qp)
  call incomplete_elliptic_integrals_standard_input_range(phi,n,m_qp,F_temp,E_temp,P_temp)
  F=complex(F_temp,0.d0); E=complex(E_temp,0.d0)
 end if
end if

end subroutine incomplete_elliptic_integrals_standard_amp_large_parameter

subroutine incomplete_elliptic_integrals(phi,m,F,E)
!!!!SJT: computes the incomplete elliptic integrals of first and second kind given amplitude phi and parameter m
!!assumes phi is real
!!assumes 0<=m (m>1 is allowed but m must be real)
 use kind_parameters
implicit none

!!input
real(dp),intent(in) :: phi,m !!elliptic amplitude and parameter

!!output
complex(dp),intent(out) :: F,E !!incomplete elliptic integrals of the first and second kind

!!internal variables
real(dp), parameter :: pii=3.1415926535897932d0
real(dp) :: piio2 !!constants
complex(dp) :: Fc,Ec !!complete elliptic integral values computed for m>1
real(dp) :: phi_std
real(dp) :: N0
integer(isp) :: N,i_sign

piio2=pii/2.d0

if ((0.d0.le.phi).and.(phi.le.piio2)) then !!standard amplitude range
 call incomplete_elliptic_integrals_standard_amp_large_parameter(phi,m,F,E)
else !!amplitude outside of standard range
 call complete_elliptic_integrals(m,Fc,Ec)
 N0=phi/pii;
 N=nint(N0,isp) !!integer multiplier
 
 phi_std=phi-N*pii !!find appropriate phi in standard range
 i_sign=1
 if (phi_std.lt.0.d0) then
  i_sign=-1
  phi_std=abs(phi_std)
 end if

 phi_std=abs(phi_std)
 call incomplete_elliptic_integrals_standard_amp_large_parameter(phi_std,m,F,E)
 F=2*N*Fc+i_sign*F
 E=2*N*Ec+i_sign*E
end if

end subroutine incomplete_elliptic_integrals

subroutine incomplete_elliptic_integral_trapezoidal_rule(phi,k,Ntheta)
!!compute the incomplete elliptic integral of the first kind via the trapezoidal rule in quad precision
!!used to compute reference values by brute force that may be used for testing
!!not intended to be used in production runs
 use kind_parameters
implicit none

!!input
real(qp),intent(in) :: phi !!Jacobi amplitude 
real(qp),intent(in) :: k !!elliptic modulus
integer(isp),intent(in) :: Ntheta !!# of mesh points

!!output (via print statement)
complex(qp) :: F !!elliptic integral of the first kind

!!internal variable
real(qp) :: theta !!variable of integration
real(qp) :: dtheta !!mesh spacing
integer(isp) :: Nindex !!index counter

dtheta=phi/real(Ntheta-1,qp)

F=(0.q0,0.q0)
theta=0.q0

do Nindex=1,Ntheta-1
 theta=real(Nindex-1,qp)*dtheta
 !F=F+dtheta*(integrand(theta,k)+integrand(theta+dtheta,k))/2.q0
 !F=F+dtheta*(integrand(theta,k)+4.q0*integrand(theta+dtheta/2.q0,k)+integrand(theta+dtheta,k))/6.q0
 F=F+dtheta*(integrand(theta,k)+3.q0*integrand(theta+dtheta/3.q0,k)+3.q0*integrand(theta+dtheta*2.q0/3.q0,k)&
            &+integrand(theta+dtheta,k))/8.d0
 !write(*,*) theta,theta+dtheta
 !write(*,*) integrand(theta,m),integrand(theta+dtheta,m)
end do
write(*,*) "IEITR: F",F

contains

 complex(qp) function integrand(theta,k)
 implicit none
 real(qp),intent(in) :: theta,k
  integrand=(1.q0,0.q0)/sqrt(complex(1.q0-(k*sin(theta))**2,0.q0))
 end function integrand

end subroutine incomplete_elliptic_integral_trapezoidal_rule
