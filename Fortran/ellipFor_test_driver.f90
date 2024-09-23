program ellipFor_test_driver
! *** Driver program demonstrating the basic use of ellipFor featuring the evaluation of:
! 1) complete Legendre elliptic integrals, 2) incomplete Legendre elliptic integrals, 
! and 3) the principal Jacobi elliptic functions ***
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
m=100._dp
call complete_elliptic_integrals(m,Fc,Ec)
write(101,'(a37)') "Complete Legendre Elliptic Integrals:"
write(101,'(a5,g26.17)') "m=   ",m
write(101,'(a5,g26.17,a2,g26.17,a1)') "Fc= (",Fc % re," ,",Fc % im,")" 
write(101,'(a5,g26.17,a2,g26.17,a1)') "Ec= (",Ec % re," ,",Ec % im,")"
write(101,'(a1)') " "

!!! 2) incomplete Legendre elliptic integrals
phi=0.75_dp*pii 
m=100._dp     
call incomplete_elliptic_integrals(phi,m,Fi,Ei)
write(101,'(a39)') "Incomplete Legendre Elliptic Integrals:"
write(101,'(a5,g26.17)') "phi= ",phi
write(101,'(a5,g26.17)') "m=   ",m
write(101,'(a5,g26.17,a2,g26.17,a1)') "Fi= (",Fi % re," ,",Fi % im,")" 
write(101,'(a5,g26.17,a2,g26.17,a1)') "Ei= (",Ei % re," ,",Ei % im,")" 
write(101,'(a1)') " "

!!! 3) Jacobi elliptic functions
u=(1._dp,1._dp) ! complex-valued u
m=100._dp      
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
write(101,'(a26)') "Jacobi Elliptic Functions:"
write(101,'(a5,g26.17,a2,g26.17,a1)') "u=  (",u % re," ,",u % im,")" 
write(101,'(a5,g26.17)') "m=   ",m
write(101,'(a5,g26.17,a2,g26.17,a1)') "sn= (",sn % re," ,",sn % im,")" 
write(101,'(a5,g26.17,a2,g26.17,a1)') "cn= (",cn % re," ,",cn % im,")" 
write(101,'(a5,g26.17,a2,g26.17,a1)') "dn= (",dn % re," ,",dn % im,")" 

close(101) !! close output file

!validate output by comparing with reference data
call validate_output

contains

 subroutine validate_output
  !! **** Compare driver output to reference data ****
  !!tolerance for validation
  real(dp),parameter :: tol=1.e-15_dp
  !!flag for test status
  logical :: test_flag
  logical :: test_flag_complete
  logical :: test_flag_incomplete
  logical :: test_flag_functions
  !!variables for reference values of elliptic functions and integrals
  real(dp) :: m_ref                   !!elliptic parameter
  real(dp) :: phi_ref                 !!Jacobi amplitude
  complex(dp) :: u_ref                !!complex argument
  complex(dp) :: Fc_ref,Ec_ref        !!complete elliptic integrals of first and second kinds
  complex(dp) :: Fi_ref,Ei_ref        !!incomplete elliptic integrals of first and second kinds 
  complex(dp) :: sn_ref,cn_ref,dn_ref !!Jacobi elliptic function values
  !!variables for errors
  real(dp) :: m_err                   !!elliptic parameter
  real(dp) :: phi_err                 !!Jacobi amplitude
  complex(dp) :: u_err                !!complex argument
  complex(dp) :: Fc_err,Ec_err        !!complete elliptic integrals of first and second kinds
  complex(dp) :: Fi_err,Ei_err        !!incomplete elliptic integrals of first and second kinds 
  complex(dp) :: sn_err,cn_err,dn_err !!Jacobi elliptic function values
  !!strings for reference data file
  character(1)  :: space,bracket
  character(2)  :: comma
  character(37) :: complete_integral_label
  character(39) :: incomplete_integral_label
  character(26) :: elliptic_functions_label
  character(5)  :: m_label,phi_label,Fc_label,Ec_label,Fi_label,Ei_label
  character(5)  :: u_label,sn_label,cn_label,dn_label 
  
  open(unit=102,file="expected_data/ellipFor_test_driver_OG.dat") !! file for reference data

  !!!!!! 1) complete Legendre elliptic integrals
  test_flag_complete=.true. !!initialize test flag

  !!read reference data
  read(102,'(a37)') complete_integral_label
  read(102,'(a5,g26.17)') m_label,m_ref
  read(102,'(a5,g26.17,a2,g26.17,a1)') Fc_label,Fc_ref % re,comma,Fc_ref % im,bracket 
  read(102,'(a5,g26.17,a2,g26.17,a1)') Ec_label,Ec_ref % re,comma,Ec_ref % im,bracket
  read(102,'(a1)') space

  !!compute relative errors (compared to reference data) and check against tolerance
  m_err=abs((m-m_ref)/m_ref)                        ; if (m_err.gt.tol)       test_flag_complete=.false.
  Fc_err % re=abs((Fc % re-Fc_ref % re)/Fc_ref % re); if (Fc_err % re.gt.tol) test_flag_complete=.false. 
  Fc_err % im=abs((Fc % im-Fc_ref % im)/Fc_ref % im); if (Fc_err % im.gt.tol) test_flag_complete=.false. 
  Ec_err % re=abs((Ec % re-Ec_ref % re)/Ec_ref % re); if (Ec_err % re.gt.tol) test_flag_complete=.false. 
  Ec_err % im=abs((Ec % im-Ec_ref % im)/Ec_ref % im); if (Ec_err % im.gt.tol) test_flag_complete=.false. 

  !!print relative errors if necessary
  if (test_flag_complete.eqv..false.) then
   print '(a37)', complete_integral_label
   print '(a11,a14)',   "Variable    ","Relative Error"
  end if
  if (m_err.gt.tol) print '(a8,g26.17)',       "m       ",m_err
  if (Fc_err % re.gt.tol) print '(a8,g26.17)', "Re[Fc]  ",Fc_err % re
  if (Fc_err % im.gt.tol) print '(a8,g26.17)', "Im[Fc]  ",Fc_err % im
  if (Ec_err % re.gt.tol) print '(a8,g26.17)', "Re[Ec]  ",Ec_err % re
  if (Ec_err % im.gt.tol) print '(a8,g26.17)', "Im[Ec]  ",Ec_err % im

  !!!!!! 2) incomplete Legendre elliptic integrals
  test_flag_incomplete=.true. !!initialize test flag

  !!read reference data
  read(102,'(a39)') incomplete_integral_label
  read(102,'(a5,g26.17)') phi_label,phi_ref
  read(102,'(a5,g26.17)') m_label,m_ref
  read(102,'(a5,g26.17,a2,g26.17,a1)') Fi_label,Fi_ref % re,comma,Fi_ref % im,bracket 
  read(102,'(a5,g26.17,a2,g26.17,a1)') Ei_label,Ei_ref % re,comma,Ei_ref % im,bracket 
  read(102,'(a1)') space 

  !!compute relative errors (compared to reference data) and check against tolerance
  phi_err=abs((phi-phi_ref)/phi_ref)                ; if (phi_err.gt.tol)     test_flag_incomplete=.false.
  m_err=abs((m-m_ref)/m_ref)                        ; if (m_err.gt.tol)       test_flag_incomplete=.false.
  Fi_err % re=abs((Fi % re-Fi_ref % re)/Fi_ref % re); if (Fi_err % re.gt.tol) test_flag_incomplete=.false. 
  Fi_err % im=abs((Fi % im-Fi_ref % im)/Fi_ref % im); if (Fi_err % im.gt.tol) test_flag_incomplete=.false. 
  Ei_err % re=abs((Ei % re-Ei_ref % re)/Ei_ref % re); if (Ei_err % re.gt.tol) test_flag_incomplete=.false. 
  Ei_err % im=abs((Ei % im-Ei_ref % im)/Ei_ref % im); if (Ei_err % im.gt.tol) test_flag_incomplete=.false. 

  !!print relative errors if necessary
  if (test_flag_incomplete.eqv..false.) then
   print '(a39)', incomplete_integral_label
   print '(a11,a14)',   "Variable    ","Relative Error"
  end if
  if (phi_err.gt.tol) print '(a8,g26.17)',     "phi     ",phi_err
  if (m_err.gt.tol) print '(a8,g26.17)',       "m       ",m_err
  if (Fi_err % re.gt.tol) print '(a8,g26.17)', "Re[Fi]  ",Fi_err % re
  if (Fi_err % im.gt.tol) print '(a8,g26.17)', "Im[Fi]  ",Fi_err % im
  if (Ei_err % re.gt.tol) print '(a8,g26.17)', "Re[Ei]  ",Ei_err % re
  if (Ei_err % im.gt.tol) print '(a8,g26.17)', "Im[Ei]  ",Ei_err % im

  !!!!!! 3) Jacobi elliptic functions
  test_flag_functions=.true. !!initialize test flag

  !!read reference data
  read(102,'(a26)') elliptic_functions_label
  read(102,'(a5,g26.17,a2,g26.17,a1)') u_label,u_ref % re,comma,u_ref % im,bracket 
  read(102,'(a5,g26.17)') m_label,m_ref
  read(102,'(a5,g26.17,a2,g26.17,a1)') sn_label,sn_ref % re,comma,sn_ref % im,bracket 
  read(102,'(a5,g26.17,a2,g26.17,a1)') cn_label,cn_ref % re,comma,cn_ref % im,bracket 
  read(102,'(a5,g26.17,a2,g26.17,a1)') dn_label,dn_ref % re,comma,dn_ref % im,bracket 

  !!compute relative errors (compared to reference data) and check against tolerance
  u_err % re=abs((u % re-u_ref % re)/u_ref % re)    ; if (u_err % re.gt.tol)  test_flag_functions=.false. 
  u_err % im=abs((u % im-u_ref % im)/u_ref % im)    ; if (u_err % im.gt.tol)  test_flag_functions=.false. 
  m_err     =abs((m-m_ref)/m_ref)                   ; if (m_err.gt.tol)       test_flag_functions=.false.
  sn_err % re=abs((sn % re-sn_ref % re)/sn_ref % re); if (sn_err % re.gt.tol) test_flag_functions=.false. 
  sn_err % im=abs((sn % im-sn_ref % im)/sn_ref % im); if (sn_err % im.gt.tol) test_flag_functions=.false. 
  cn_err % re=abs((cn % re-cn_ref % re)/cn_ref % re); if (cn_err % re.gt.tol) test_flag_functions=.false. 
  cn_err % im=abs((cn % im-cn_ref % im)/cn_ref % im); if (cn_err % im.gt.tol) test_flag_functions=.false. 
  dn_err % re=abs((dn % re-dn_ref % re)/dn_ref % re); if (dn_err % re.gt.tol) test_flag_functions=.false. 
  dn_err % im=abs((dn % im-dn_ref % im)/dn_ref % im); if (dn_err % im.gt.tol) test_flag_functions=.false. 

  !!print relative errors if necessary
  if (test_flag_functions.eqv..false.) then
   print '(a26)', elliptic_functions_label
   print '(a11,a14)',   "Variable    ","Relative Error"
  end if
  if (u_err % re.gt.tol) print '(a8,g26.17)',  "Re[u]   ",u_err % re
  if (u_err % im.gt.tol) print '(a8,g26.17)',  "Im[u]   ",u_err % im
  if (m_err.gt.tol) print '(a8,g26.17)',       "m       ",m_err
  if (sn_err % re.gt.tol) print '(a8,g26.17)', "Re[sn]  ",sn_err % re
  if (sn_err % im.gt.tol) print '(a8,g26.17)', "Im[sn]  ",sn_err % im
  if (cn_err % re.gt.tol) print '(a8,g26.17)', "Re[cn]  ",cn_err % re
  if (cn_err % im.gt.tol) print '(a8,g26.17)', "Im[cn]  ",cn_err % im
  if (dn_err % re.gt.tol) print '(a8,g26.17)', "Re[dn]  ",dn_err % re
  if (dn_err % im.gt.tol) print '(a8,g26.17)', "Im[dn]  ",dn_err % im

  close(102) !! close file for reference data

  test_flag=all([test_flag_complete,test_flag_incomplete,test_flag_functions])

  if (test_flag) then
   print '(a64)', "Driver output was successfully validated against reference data."
  else
   print '(a80)', "Driver output significantly differs from reference data. See error report above."
  end if
 end subroutine validate_output

end program ellipFor_test_driver

