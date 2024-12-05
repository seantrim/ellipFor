program test_material_driver
! *** Driver program for generating ellipFor data for test problems featuring:
! 1) complete Legendre elliptic integrals, 2) incomplete Legendre elliptic integrals, 
! and 3) the principal Jacobi elliptic functions ***
! * test data is automatically compared with SageMath values *
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_is_negative,ieee_value,ieee_positive_inf
 use kind_parameters
 use elliptic, only: complete_elliptic_integrals,incomplete_elliptic_integrals,Jacobi_elliptic_functions 
 implicit none

 
 real(dp),parameter :: tol_complete  =5.e-15_dp ! tolerance for relative errors of complete Legendre elliptic integrals
 real(dp),parameter :: tol_incomplete=5.e-15_dp ! tolerance for relative errors of incomplete Legendre elliptic integrals
 real(dp),parameter :: tol_functions =9.e-14_dp ! tolerance for relative errors of Jacobi elliptic functions
 real(dp),parameter :: tol_special   =5.e-15_dp ! tolerance for relative errors of special values of all ellipFor functions

 logical :: special_value_tests_passed      ! flag to indicate if all special value tests were passed
 logical :: complete_integral_test_passed   ! flag to indicate if general complete elliptic integral test passed
 logical :: incomplete_integral_test_passed ! flag to indicate if general incomplete elliptic integral test passed
 logical :: functions_test_passed           ! flag to indicate if general Jacobi elliptic function test passed
 logical :: general_testing_passed          ! flag to indicate if all general value tests were passed
 
 call test_special_values(tol_special,special_value_tests_passed)
 if (special_value_tests_passed) then
  print '(a39)', "- All special value tests were passed -"
 else
  print '(a63)', "------ Warning: Failures detected in special value tests ------"
 end if
 print '(a1)', " "


 print '(a46)', "Testing general values for ellipFor functions:"
 call test_complete_elliptic_integrals(tol_complete,complete_integral_test_passed)
 call test_incomplete_elliptic_integrals(tol_incomplete,incomplete_integral_test_passed)
 call test_Jacobi_elliptic_functions(tol_functions,functions_test_passed)

 general_testing_passed=all([complete_integral_test_passed,incomplete_integral_test_passed,functions_test_passed])
 if (general_testing_passed) then
  print '(a39)', "- All general value tests were passed -"
 else
  print '(a63)', "------ Warning: Failures detected in general value tests ------"
 end if
 print '(a1)', " "

 if (special_value_tests_passed.and.general_testing_passed) then
  print '(a47)', "* All ellipFor tests were successfully passed *"
 else
  print '(a76)', "* Warning: ellipFor test failures detected -- see above output for details *"
 end if
contains
 
 subroutine test_special_values(tol,all_tests_passed)
  ! *** Test ellipFor for special values of Legendre elliptic integrals and Jacobi elliptic functions ***
  ! input
  real(dp),intent(in) :: tol              ! tolerance for automatic test verification

  ! output
  logical,intent(out) :: all_tests_passed ! flag to indicate if all tests passed successfully

  ! local variables
  ! parameters
  real(dp), parameter :: pii=3.1415926535897932_dp ! pi
  
  ! variables for elliptic functions and integrals
  real(dp) :: m           !!elliptic parameter
  real(dp) :: mc          !!complimentary parameter
  real(dp) :: k           !!elliptic modulus
  real(dp) :: phi         !!Jacobi amplitude
  complex(dp) :: kc       !!complimentary modulus
  complex(dp) :: u        !!complex argument
  complex(dp) :: Fc,Ec    !!complete elliptic integrals of first and second kinds
  complex(dp) :: Fcc,Ecc  !!complimentary complete elliptic integrals of first and second kinds
  complex(dp) :: Fi,Ei    !!incomplete elliptic integrals of first and second kinds 
  complex(dp) :: sn,cn,dn !!Jacobi elliptic function values

  ! variables for testing identities
  integer(isp) :: i           ! loop index
  integer(isp) :: nsize       ! # of values to be tested for identities requiring a range of arguments
  real(dp) :: m_min,m_max     ! min/max m values
  real(dp) :: phi_min,phi_max ! min/max phi values
  real(dp),allocatable :: m_array(:),phi_array(:) ! arrays for storing randomly selected argument values
  logical :: verified                      ! flag for verify_complex_values subroutine
  logical :: Fi_verified,Fi_verified_total ! verification flags
  logical :: Ei_verified,Ei_verified_total
  logical :: sn_verified,sn_verified_total
  logical :: cn_verified,cn_verified_total
  logical :: dn_verified,dn_verified_total

  ! variable for positive infinity (assigned below)
  real(dp) :: positive_inf 

  ! assign positive infinity variable
  positive_inf=ieee_value(0._dp,ieee_positive_inf)

  ! initialize test flag
  all_tests_passed=.true. ! for one or more failures, this flag shall be set to .false.

  print '(a46)', "Testing special values for ellipFor functions:"
  ! test K(0)=pi/2, E(0)=pi/2
  m=0._dp
  call complete_elliptic_integrals(m,Fc,Ec)
  print '(a23)', " *** Test K(0)=pi/2 ***"
  print '(a15,a1,f18.16,a1,f18.16,a1)', " EllipFor K(0)=","(",Fc%re    ,",",Fc%im,")"
  print '(a15,a1,f18.16,a1,f18.16,a1)', " Exact    K(0)=","(",pii/2._dp,",",0._dp,")"
  call verify_complex_values(Fc,cmplx(pii/2._dp,0._dp,dp),tol,.true.,verified)
  if (.not.verified) all_tests_passed=.false. ! if test fails, set test flag accordingly
  print '(a1)', " "
  print '(a23)', " *** Test E(0)=pi/2 ***"
  print '(a15,a1,f18.16,a1,f18.16,a1)', " EllipFor E(0)=","(",Ec%re    ,",",Ec%im,")"
  print '(a15,a1,f18.16,a1,f18.16,a1)', " Exact    E(0)=","(",pii/2._dp,",",0._dp,")"
  call verify_complex_values(Ec,cmplx(pii/2._dp,0._dp,dp),tol,.true.,verified)
  if (.not.verified) all_tests_passed=.false. ! if test fails, set test flag accordingly
  print '(a1)', " "
 
  ! test K(1)=Infinity and E(1)=1
  m=1._dp
  call complete_elliptic_integrals(m,Fc,Ec)
  print '(a27)', " *** Test K(1)=Infinity ***"
  print '(a43)', " Note: the definition from SageMath is used"
  print '(a15,a1,f18.16,a1,f18.16,a1)', " EllipFor K(1)=","(",Fc%re               ,",",Fc%im,")"
  print '(a15,a1,a18,a1,f18.16,a1)'   , " Exact    K(1)=","(","          Infinity",",",0._dp,")"
  call verify_complex_values(Fc,cmplx(positive_inf,0._dp,dp),tol,.true.,verified) 
  if (.not.verified) all_tests_passed=.false. ! if test fails, set test flag accordingly
  print '(a1)', " "
  print '(a20)', " *** Test E(1)=1 ***"
  print '(a15,a1,f18.16,a1,f18.16,a1)', " EllipFor E(1)=","(",Ec%re               ,",",Ec%im,")"
  print '(a15,a1,f18.16,a1,f18.16,a1)', " Exact    E(1)=","(",1._dp               ,",",0._dp,")"
  call verify_complex_values(Ec,cmplx(1._dp,0._dp,dp),tol,.true.,verified)
  if (.not.verified) all_tests_passed=.false. ! if test fails, set test flag accordingly
  print '(a1)', " "
 
  ! test F(pi/2|1)=Infinity, E(pi/2|1)=1
  phi=pii/2._dp
  m=1._dp
  call incomplete_elliptic_integrals(phi,m,Fi,Ei)
  print '(a32)', " *** Test F(pi/2|1)=Infinity ***"
  print '(a43)', " Note: the definition from SageMath is used"
  print '(a20,a1,f18.16,a1,f18.16,a1)', " EllipFor F(pi/2|1)=","(",Fi%re               ,",",Fi%im,")"
  print '(a20,a1,a18,a1,f18.16,a1)'   , " Exact    F(pi/2|1)=","(","          Infinity",",",0._dp,")"
  call verify_complex_values(Fi,cmplx(positive_inf,0._dp,dp),tol,.true.,verified) 
  if (.not.verified) all_tests_passed=.false. ! if test fails, set test flag accordingly
  print '(a1)', " "
  print '(a25)', " *** Test E(pi/2|1)=1 ***"
  print '(a20,a1,f18.16,a1,f18.16,a1)', " EllipFor E(pi/2|1)=","(",Ei%re               ,",",Ei%im,")"
  print '(a20,a1,f18.16,a1,f18.16,a1)', " Exact    E(pi/2|1)=","(",1._dp               ,",",0._dp,")"
  call verify_complex_values(Ei,cmplx(1._dp,0._dp,dp),tol,.true.,verified)
  if (.not.verified) all_tests_passed=.false. ! if test fails, set test flag accordingly
  print '(a1)', " "
 
  ! F(phi|0)=phi, E(phi|0)=phi
  print '(a72)', " *** Test F(phi|0)=phi and E(phi|0)=phi for -10**12 <= phi <= 10**12 ***"
  m=0._dp
  phi_min=-1.e12_dp; phi_max=1.e12_dp ! range of m
  nsize=40000_isp ! # of m values tested
  allocate(phi_array(1:nsize))
  phi_array(1)=phi_min; phi_array(nsize)=phi_max ! ensure min/max m values are tested
  call random_number(phi_array(2:nsize-1))                    ! random values from [0,1)
  phi_array(2:nsize-1)=(phi_max-phi_min)*phi_array(2:nsize-1)+phi_min ! use random values to populate remaining argument values
  Fi_verified_total=.true. ! initialize verification flags
  Ei_verified_total=.true. 
  do i=1,nsize
   call incomplete_elliptic_integrals(phi_array(i),m,Fi,Ei)
   call verify_complex_values(Fi,cmplx(phi_array(i),0._dp,dp),tol,.false.,Fi_verified)
   call verify_complex_values(Ei,cmplx(phi_array(i),0._dp,dp),tol,.false.,Ei_verified)
   Fi_verified_total=Fi_verified_total.and.Fi_verified ! update verification flags
   Ei_verified_total=Ei_verified_total.and.Ei_verified
  end do
  print '(a7,i6,a41)', " Note: ",nsize," randomly selected phi values were tested"
  if (Fi_verified) then
   print '(a25)', " F(phi|0)=phi test passed"
  else
   print '(a48)', " ------ Warning: F(phi|0)=phi test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  if (Ei_verified) then
   print '(a25)', " E(phi|0)=phi test passed"
  else
   print '(a48)', " ------ Warning: E(phi|0)=phi test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  print '(a1)', " "

  ! F(0|m)=0, E(0|m)=0
  print '(a56)', " *** Test F(0|m)=0 and E(0|m)=0 for 0 <= m <= 10**12 ***"
  phi=0._dp
  m_min=0._dp; m_max=1.e12_dp ! range of m
  nsize=20000_isp ! # of m values tested
  allocate(m_array(1:nsize))
  m_array(1)=m_min; m_array(nsize)=m_max ! ensure min/max m values are tested
  call random_number(m_array(2:nsize-1))                    ! random values from [0,1)
  m_array(2:nsize-1)=(m_max-m_min)*m_array(2:nsize-1)+m_min ! use random values to populate remaining argument values
  Fi_verified_total=.true. ! initialize verification flags
  Ei_verified_total=.true. 
  do i=1,nsize
   call incomplete_elliptic_integrals(phi,m_array(i),Fi,Ei)
   call verify_complex_values(Fi,cmplx(0._dp,0._dp,dp),tol,.false.,Fi_verified)
   call verify_complex_values(Ei,cmplx(0._dp,0._dp,dp),tol,.false.,Ei_verified)
   Fi_verified_total=Fi_verified_total.and.Fi_verified ! update verification flags
   Ei_verified_total=Ei_verified_total.and.Ei_verified
  end do
  deallocate(m_array)
  print '(a7,i6,a39)', " Note: ",nsize," randomly selected m values were tested"
  if (Fi_verified) then
   print '(a21)', " F(0|m)=0 test passed"
  else
   print '(a44)', " ------ Warning: F(0|m)=0 test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  if (Ei_verified) then
   print '(a21)', " E(0|m)=0 test passed"
  else
   print '(a44)', " ------ Warning: E(0|m)=0 test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  print '(a1)', " "


  ! sn(0|m)=0, cn(0|m)=1, dn(0|m)=1
  print '(a70)', " *** Test sn(0|m)=0, cn(0|m)=1, and dn(0|m)=1 for 0 <= m <= 10**12 ***"
  u=(0._dp,0._dp)
  m_min=0._dp; m_max=1.e12_dp ! range of m
  nsize=20000_isp ! # of m values tested
  allocate(m_array(1:nsize))
  m_array(1)=m_min; m_array(nsize)=m_max ! ensure min/max m values are tested
  call random_number(m_array(2:nsize-1))                    ! random values from [0,1)
  m_array(2:nsize-1)=(m_max-m_min)*m_array(2:nsize-1)+m_min ! use random values to populate remaining argument values
  sn_verified_total=.true. ! initialize verification flags
  cn_verified_total=.true. 
  dn_verified_total=.true. 
  do i=1,nsize
   call Jacobi_elliptic_functions(u,m_array(i),sn,cn,dn)
   call verify_complex_values(sn,cmplx(0._dp,0._dp,dp),tol,.false.,sn_verified)
   call verify_complex_values(cn,cmplx(1._dp,0._dp,dp),tol,.false.,cn_verified)
   call verify_complex_values(dn,cmplx(1._dp,0._dp,dp),tol,.false.,dn_verified)
   sn_verified_total=sn_verified_total.and.sn_verified ! update verification flags
   cn_verified_total=cn_verified_total.and.cn_verified
   dn_verified_total=dn_verified_total.and.dn_verified
  end do
  deallocate(m_array)
  print '(a7,i6,a39)', " Note: ",nsize," randomly selected m values were tested"
  if (sn_verified) then
   print '(a22)', " sn(0|m)=0 test passed"
  else
   print '(a45)', " ------ Warning: sn(0|m)=0 test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  if (cn_verified) then
   print '(a22)', " cn(0|m)=0 test passed"
  else
   print '(a45)', " ------ Warning: cn(0|m)=0 test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  if (dn_verified) then
   print '(a22)', " dn(0|m)=0 test passed"
  else
   print '(a45)', " ------ Warning: dn(0|m)=0 test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  print '(a1)', " "

   
  ! sn(K(m)|m)=1, cn(K(m)|m)=0, dn(K(m)|m)=sqrt(1-m)
  print '(a87)', " *** Test sn(K(m)|m)=1, cn(K(m)|m)=0, and dn(K(m)|m)=(1-m)**0.5 for 0 <= m <= 10**7 ***"
  m_min=0._dp; m_max=1.e7_dp ! range of m
  nsize=16000_isp ! # of m values tested
  allocate(m_array(1:nsize))
  m_array(1)=m_min; m_array(nsize)=m_max ! ensure min/max m values are tested
  call random_number(m_array(2:nsize-1))                    ! random values from [0,1)
  m_array(2:nsize-1)=(m_max-m_min)*m_array(2:nsize-1)+m_min ! use random values to populate remaining argument values
  sn_verified_total=.true. ! initialize verification flags
  cn_verified_total=.true. 
  dn_verified_total=.true. 
  do i=1,nsize
   call complete_elliptic_integrals(m_array(i),Fc,Ec)
   u=Fc
   mc=1._dp-m_array(i); kc=sqrt(cmplx(mc,0._dp,dp))
   call Jacobi_elliptic_functions(u,m_array(i),sn,cn,dn)
   call verify_complex_values(sn,cmplx(1._dp,0._dp,dp),tol,.false.,sn_verified)
   call verify_complex_values(cn,cmplx(0._dp,0._dp,dp),tol,.false.,cn_verified)
   call verify_complex_values(dn,kc,tol,.false.,dn_verified)
   sn_verified_total=sn_verified_total.and.sn_verified ! update verification flags
   cn_verified_total=cn_verified_total.and.cn_verified
   dn_verified_total=dn_verified_total.and.dn_verified
  end do
  deallocate(m_array)
  print '(a7,i6,a39)', " Note: ",nsize," randomly selected m values were tested"
  if (sn_verified) then
   print '(a25)', " sn(K(m)|m)=1 test passed"
  else
   print '(a48)', " ------ Warning: sn(K(m)|m)=1 test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  if (cn_verified) then
   print '(a25)', " cn(K(m)|m)=0 test passed"
  else
   print '(a48)', " ------ Warning: cn(K(m)|m)=0 test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  if (dn_verified) then
   print '(a34)', " dn(K(m)|m)=(1-m)**0.5 test passed"
  else
   print '(a57)', " ------ Warning: dn(K(m)|m)=(1-m)**0.5 test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  print '(a1)', " "


  ! sn(K(m)+i*K(m')|m)=1/k, cn(K(m)+i*K(m')|m)=-i*k'/k, dn(K(m)+i*K(m')|m)=0
  print '(a104)', " *** Test sn(K(m)+i*K(m')|m)=1/k, cn(K(m)+i*K(m')|m)=-i*k'/k, and dn(K(m)+i*K(m')|m)=0 for 0 < m < 1 ***"
  m_min=0._dp+tol; m_max=1.e0_dp-tol ! range of m
  nsize=10000_isp ! # of m values tested
  allocate(m_array(1:nsize))
  m_array(1)=m_min; m_array(nsize)=m_max ! ensure min/max m values are tested
  call random_number(m_array(2:nsize-1))                    ! random values from [0,1)
  m_array(2:nsize-1)=(m_max-m_min)*m_array(2:nsize-1)+m_min ! use random values to populate remaining argument values
  sn_verified_total=.true. ! initialize verification flags
  cn_verified_total=.true. 
  dn_verified_total=.true. 
  do i=1,nsize
   call complete_elliptic_integrals(m_array(i),Fc,Ec)
   mc=1._dp-m_array(i); kc=sqrt(cmplx(mc,0._dp,dp)); k=sqrt(m_array(i))
   call complete_elliptic_integrals(mc,Fcc,Ecc)
   u=Fc+(0._dp,1._dp)*Fcc
   call Jacobi_elliptic_functions(u,m_array(i),sn,cn,dn)
   call verify_complex_values(sn,cmplx(1._dp/k,0._dp,dp),tol,.false.,sn_verified)
   call verify_complex_values(cn,-(0._dp,1._dp)*kc/k,tol,.false.,cn_verified)
   call verify_complex_values(dn,(0._dp,0._dp),tol,.false.,dn_verified)
   sn_verified_total=sn_verified_total.and.sn_verified ! update verification flags
   cn_verified_total=cn_verified_total.and.cn_verified
   dn_verified_total=dn_verified_total.and.dn_verified
  end do
  deallocate(m_array)
  print '(a7,i6,a39)', " Note: ",nsize," randomly selected m values were tested"
  if (sn_verified) then
   print '(a35)', " sn(K(m)+i*K(m')|m)=1/k test passed"
  else
   print '(a58)', " ------ Warning: sn(K(m)+i*K(m')|m)=1/k test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  if (cn_verified) then
   print '(a39)', " cn(K(m)+i*K(m')|m)=-i*k'/k test passed"
  else
   print '(a62)', " ------ Warning: cn(K(m)+i*K(m')|m)=-i*k'/k test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  if (dn_verified) then
   print '(a33)', " dn(K(m)+i*K(m')|m)=0 test passed"
  else
   print '(a56)', " ------ Warning: dn(K(m)+i*K(m')|m)=0 test failed ------"
   all_tests_passed=.false. ! if test fails, set test flag accordingly
  end if
  print '(a1)', " "

 end subroutine test_special_values

 function negative_infinity_test(x) result(neg_inf_test)
  ! *** Test for negative infinity ***
  ! input
  real(dp),intent(in) :: x
  ! output
  logical :: neg_inf_test
  neg_inf_test=(.not.ieee_is_finite(x)).and.(ieee_is_negative(x))
 end function

 function positive_infinity_test(x) result(pos_inf_test)
  ! *** Test for negative infinity ***
  ! input
  real(dp),intent(in) :: x
  ! output
  logical :: pos_inf_test
  pos_inf_test=(.not.ieee_is_finite(x)).and.(.not.ieee_is_negative(x))
 end function

 subroutine verify_complex_values(z,z_ref,tol,verbose,verified)
  ! *** Verify that the complex variable z agrees with z_ref within the specified tolerance ***
  ! note: z_ref can include infinite real and/or imaginary parts
  ! input
  complex(dp),intent(in) :: z     ! trial value
  complex(dp),intent(in) :: z_ref ! reference value
  real(dp)   ,intent(in) :: tol   ! tolerance (max around 1.e-15_dp for precise agreement)
  logical    ,intent(in) :: verbose ! flag to indicate verbose output (silent if .false.) 

  ! output
  logical    ,intent(out) :: verified ! flag to indicate if trial and reference values agree within tolerance

  ! local variables
  logical :: verify_real_part ! flag for testing the real part
  logical :: verify_imag_part ! flag for testing the imaginary part
  logical :: pos_infinity_ref ! flag for positive infinity in reference value
  logical :: neg_infinity_ref ! flag for positive infinity in reference value

  ! initialize flags
  verify_real_part=.false.
  verify_imag_part=.false.

  ! test real part
  pos_infinity_ref=positive_infinity_test(z_ref%re) ! check for positive infinity in reference value
  neg_infinity_ref=negative_infinity_test(z_ref%re) ! check for negative infinity in reference value
  if (pos_infinity_ref) then       ! if reference value is +Infinity 
   verify_real_part=positive_infinity_test(z%re)
  else if (neg_infinity_ref) then  ! if reference value is -Infinity
   verify_real_part=negative_infinity_test(z%re)
  else if (z_ref%re.eq.0._dp) then ! absolute difference if reference value is zero
   if (abs(z%re-z_ref%re).lt.abs(tol)) verify_real_part=.true.
  else                             ! relative difference otherwise
   if (abs((z%re-z_ref%re)/z_ref%re).lt.abs(tol)) verify_real_part=.true.
  end if

  ! test imaginary part
  pos_infinity_ref=positive_infinity_test(z_ref%im) ! check for positive infinity in reference value
  neg_infinity_ref=negative_infinity_test(z_ref%im) ! check for negative infinity in reference value
  if (pos_infinity_ref) then       ! if reference value is +Infinity 
   verify_real_part=positive_infinity_test(z%im)
  else if (neg_infinity_ref) then  ! if reference value is -Infinity
   verify_real_part=negative_infinity_test(z%im)
  else if (z_ref%im.eq.0._dp) then ! absolute difference if reference value is zero
   if (abs(z%im-z_ref%im).lt.abs(tol)) verify_imag_part=.true.
  else                             ! relative difference otherwise
   if (abs((z%im-z_ref%im)/z_ref%im).lt.abs(tol)) verify_imag_part=.true.
  end if

  ! verification result
  verified=verify_real_part.and.verify_imag_part

  ! output of results
  if (verbose) then
   if (verified) then 
    print '(a12)', " Test Passed"
   else if (verify_real_part) then
    print '(a31)', " Test Failed for Imaginary Part"
   else if (verify_imag_part) then
    print '(a26)', " Test Failed for Real Part"
   else
    print '(a12)', " Test Failed"
   end if
  end if
 end subroutine verify_complex_values

 subroutine test_complete_elliptic_integrals(tol,test_passed)
  ! *** Compute errors of complete elliptic integral values relative to reference CAS (SageMath) values ***
  ! input
  real(dp),intent(in) :: tol ! tolerance for tests

  ! output
  logical,intent(out) :: test_passed             ! flag to indicate if test was passed within tolerance

  ! local variables
  logical      :: file_exists   ! flag for the existence of files
  logical      :: Fc_test_passed,Ec_test_passed ! flags for tests of Fc and Ec values
  real(dp)     :: m             ! the parameter
  complex(dp)  :: Fc_CAS,Ec_CAS ! CAS reference values for complete Legendre elliptic intergrals of first and second kinds
  complex(dp)  :: Fc,Ec         ! ellipFor values for complete Legendre elliptic integrals of first and second kinds
  complex(dp)  :: Fc_err,Ec_err ! error values for complete Legendre elliptic integrals relative to CAS values
  complex(dp)  :: Fc_err_max,Ec_err_max   ! max error values
  integer(isp) :: nargs                   ! # of m values tested
  integer(isp) :: io                      ! for read iostat argument
  ! introduce test in terminal: range of m should match randomized list
  print '(a75)', " *** Test K(m) and E(m) relative to SageMath values for 0 <= m <= 10**7 ***" 
  open(unit=101,file="expected_data/CAS_complete.dat",status="old",action="read") ! CAS reference file
  inquire(file="error_complete.dat",exist=file_exists)
  if (file_exists) then
   open(unit=202,file="error_complete.dat",status="replace",action="write")       ! file for error data
  else
   open(unit=202,file="error_complete.dat",status="new",action="write")           ! file for error data
  end if
  ! initialize stat variables
  Fc_err_max=(0._dp,0._dp)
  Ec_err_max=(0._dp,0._dp)
  nargs=0 ! intialize counter for randomized argument list
  do
   ! read in CAS reference values 
   read(101,*,iostat=io) m,Fc_CAS%re,Fc_CAS%im,Ec_CAS%re,Ec_CAS%im
   if (io.lt.0_isp) exit ! exit loop if end of file reached

   ! compute ellipFor values
   call complete_elliptic_integrals(m,Fc,Ec)

   ! compute errors
   Fc_err%re=error(Fc%re,Fc_CAS%re)
   Fc_err%im=error(Fc%im,Fc_CAS%im)
   Ec_err%re=error(Ec%re,Ec_CAS%re)
   Ec_err%im=error(Ec%im,Ec_CAS%im)

   ! update max error values
   if (Fc_err%re.gt.Fc_err_max%re) Fc_err_max%re=Fc_err%re
   if (Fc_err%im.gt.Fc_err_max%im) Fc_err_max%im=Fc_err%im
   if (Ec_err%re.gt.Ec_err_max%re) Ec_err_max%re=Ec_err%re
   if (Ec_err%im.gt.Ec_err_max%im) Ec_err_max%im=Ec_err%im

   ! print error data to file
   write(202,'(5(g27.17))') m,Fc_err%re,Fc_err%im,Ec_err%re,Ec_err%im
   nargs=nargs+1
  end do
  close(101)
  close(202)
  print '(a7,i6,a39)', " Note: ",nargs," randomly selected m values were tested"
  if (max(Fc_err_max%re,Fc_err_max%im) .lt. tol) then
   print '(a17)', " K(m) test passed"
   Fc_test_passed=.true.
  else
   print '(a63,g25.17,a7)', " ------ Warning: K(m) test failed with a max relative error of ",&
   &max(Fc_err_max%re,Fc_err_max%im)," ------"
   Fc_test_passed=.false.
  end if 
  if (max(Ec_err_max%re,Ec_err_max%im) .lt. tol) then
   print '(a17)', " E(m) test passed"
   Ec_test_passed=.true.
  else
   print '(a63,g25.17,a7)', " ------ Warning: E(m) test failed with a max relative error of ",&
   &max(Ec_err_max%re,Ec_err_max%im)," ------"
   Ec_test_passed=.false.
  end if
  if (Fc_test_passed.and.Ec_test_passed) then
   test_passed=.true.
  else
   test_passed=.false.
  end if
  print '(a1)', " " 
 end subroutine test_complete_elliptic_integrals

 subroutine test_incomplete_elliptic_integrals(tol,test_passed)
  ! *** Compute errors of incomplete elliptic integral values relative to reference CAS (SageMath) values ***

  ! input
  real(dp),intent(in) :: tol ! tolerance for tests

  ! output
  logical,intent(out) :: test_passed             ! flag to indicate if test was passed within tolerance

  ! local variables
  logical      :: file_exists   ! flag for the existence of files
  logical      :: Fi_test_passed,Ei_test_passed ! flags for tests of Fc and Ec values
  real(dp)     :: phi           ! Jacobi amplitude
  real(dp)     :: m             ! the parameter
  complex(dp)  :: Fi_CAS,Ei_CAS ! CAS reference values for complete Legendre elliptic intergrals of first and second kinds
  complex(dp)  :: Fi,Ei         ! ellipFor values for complete Legendre elliptic integrals of first and second kinds
  complex(dp)  :: Fi_err,Ei_err ! error values for complete Legendre elliptic integrals relative to CAS values
  complex(dp)  :: Fi_err_max,Ei_err_max   ! max error values
  complex(dp)  :: Fi_err_mean,Ei_err_mean ! mean error values
  integer(isp) :: nargs                 ! # of argument combinations tested
  integer(isp) :: io                    ! for read iostat argument 

  ! introduce test in terminal: range of m should match randomized list
  print '(a110)', " *** Test F(phi|m) and E(phi|m) relative to SageMath values for -10**9 <= phi <= 10**9 and 0 <= m <= 10**7 ***" 
  ! open files
  open(unit=101,file="expected_data/CAS_incomplete.dat",status="old",action="read") ! CAS reference file
  inquire(file="error_incomplete.dat",exist=file_exists)
  if (file_exists) then
   open(unit=202,file="error_incomplete.dat",status="replace",action="write")       ! file for error data
  else
   open(unit=202,file="error_incomplete.dat",status="new",action="write")           ! file for error data
  end if
  ! initialize stat variables
  Fi_err_max =(0._dp,0._dp)
  Ei_err_max =(0._dp,0._dp)
  Fi_err_mean=(0._dp,0._dp)
  Ei_err_mean=(0._dp,0._dp)
  nargs=0 ! intialize counter for randomized argument list
  do
   ! read in CAS reference values 
   read(101,*,iostat=io) phi,m,Fi_CAS%re,Fi_CAS%im,Ei_CAS%re,Ei_CAS%im
   if (io.lt.0_isp) exit ! exit loop if end of file reached

   ! compute ellipFor values
   call incomplete_elliptic_integrals(phi,m,Fi,Ei)

   ! compute errors
   Fi_err%re=error(Fi%re,Fi_CAS%re)
   Fi_err%im=error(Fi%im,Fi_CAS%im)
   Ei_err%re=error(Ei%re,Ei_CAS%re)
   Ei_err%im=error(Ei%im,Ei_CAS%im)

   ! update max error values
   if (Fi_err%re.gt.Fi_err_max%re) Fi_err_max%re=Fi_err%re
   if (Fi_err%im.gt.Fi_err_max%im) Fi_err_max%im=Fi_err%im
   if (Ei_err%re.gt.Ei_err_max%re) Ei_err_max%re=Ei_err%re
   if (Ei_err%im.gt.Ei_err_max%im) Ei_err_max%im=Ei_err%im

   ! compute sum used for mean
   Fi_err_mean=Fi_err_mean+Fi_err
   Ei_err_mean=Ei_err_mean+Ei_err

   ! print error data to file
   write(202,'(6(g27.17))') phi,m,Fi_err%re,Fi_err%im,Ei_err%re,Ei_err%im
   nargs=nargs+1
  end do
  close(101)
  close(202)
  ! finish mean calculation
  Fi_err_mean=Fi_err_mean/nargs
  Ei_err_mean=Ei_err_mean/nargs
  print '(a7,i6,a56)', " Note: ",nargs," randomly selected combinations of phi and m were tested"
  if (max(Fi_err_max%re,Fi_err_max%im) .lt. tol) then
   print '(a21)', " F(phi|m) test passed"
   Fi_test_passed=.true.
  else
   print '(a44)', " ------ Warning: F(phi|m) test failed ------" 
   print '(a23,g25.17,a1,g25.17,a1)', "  Relative error max =(",Fi_err_max %re,",",Fi_err_max %im,")"
   print '(a23,g25.17,a1,g25.17,a1)', "  Relative error mean=(",Fi_err_mean%re,",",Fi_err_mean%im,")"
   Fi_test_passed=.false.
  end if 
  if (max(Ei_err_max%re,Ei_err_max%im) .lt. tol) then
   print '(a21)', " E(phi|m) test passed"
   Ei_test_passed=.true.
  else
   print '(a44)', " ------ Warning: E(phi|m) test failed ------" 
   print '(a23,g25.17,a1,g25.17,a1)', "  Relative error max =(",Ei_err_max %re,",",Ei_err_max %im,")"
   print '(a23,g25.17,a1,g25.17,a1)', "  Relative error mean=(",Ei_err_mean%re,",",Ei_err_mean%im,")"
   Ei_test_passed=.false.
  end if
  if (Fi_test_passed.and.Ei_test_passed) then
   test_passed=.true.
  else
   test_passed=.false.
  end if
  print '(a1)', " " 

 end subroutine test_incomplete_elliptic_integrals

 subroutine test_Jacobi_elliptic_functions(tol,test_passed)
  ! *** Compute errors of Jacobi elliptic function values relative to reference CAS (SageMath) values ***

  ! input
  real(dp),intent(in) :: tol              ! tolerance for tests

  ! output
  logical,intent(out) :: test_passed      ! flag to indicate if test was passed within tolerance

  ! local variables
  logical      :: file_exists             ! flag for the existence of files
  logical      :: sn_test_passed,cn_test_passed,dn_test_passed ! test flags
  complex(dp)  :: u                       ! first input parameter for the Jacobi elliptic functions
  real(dp)     :: m                       ! second input parameter for the Jacobi elliptic functions 
  real(dp),parameter :: CAS_tol=5.e-16_dp ! tolerance for CAS agreement test
  real(dp)     :: CAS_test_array(1:6)     ! array for testing agreement between CAS packages
  logical      :: CAS_test                ! flag for CAS agreement test
  complex(dp)  :: sn_SM,cn_SM,dn_SM       ! SageMath values for complete Legendre elliptic intergrals of first and second kinds
  complex(dp)  :: sn_Math,cn_Math,dn_Math ! Mathematica values for complete Legendre elliptic intergrals of first and second kinds
  complex(dp)  :: sn,cn,dn                ! ellipFor values for complete Legendre elliptic integrals of first and second kinds
  complex(dp)  :: sn_err,cn_err,dn_err    ! error values for complete Legendre elliptic integrals relative to CAS values
  complex(dp)  :: sn_err_max,cn_err_max,dn_err_max    ! max errors for complete Legendre elliptic integrals relative to CAS values
  complex(dp)  :: sn_err_mean,cn_err_mean,dn_err_mean ! max errors for complete Legendre elliptic integrals relative to CAS values
  integer(isp) :: nargs                   ! # of arguments
  integer(isp) :: io                      ! for read iostat argument

  ! introduce test in terminal: range of m should match randomized list
  print '(a71)', " *** Test sn(u|m), cn(u|m), and dn(u|m) relative to SageMath values ***"
  print '(a74)', "  Argument ranges: -1 <= Re[u] <= 1, -1 <= Im[u] <= 1, and 0 <= m <= 10**2" 
  open(unit=101,file="expected_data/CAS_functions.dat",status="old",action="read") ! CAS reference file
  inquire(file="error_functions.dat",exist=file_exists)
  if (file_exists) then
   open(unit=202,file="error_functions.dat",status="replace",action="write")       ! file for error data
  else
   open(unit=202,file="error_functions.dat",status="new",action="write")           ! file for error data
  end if
  ! initialize stat variables
  sn_err_max =(0._dp,0._dp)
  cn_err_max =(0._dp,0._dp)
  dn_err_max =(0._dp,0._dp)
  sn_err_mean=(0._dp,0._dp)
  cn_err_mean=(0._dp,0._dp)
  dn_err_mean=(0._dp,0._dp)
  nargs=0 ! intialize counter for randomized argument list
  do
   ! read in CAS reference values 
   read(101,*,iostat=io) u%re,u%im,m,&
                        &sn_SM  %re,sn_SM  %im,cn_SM  %re,cn_SM  %im,dn_SM  %re,dn_SM  %im,&
                        &sn_Math%re,sn_Math%im,cn_Math%re,cn_Math%im,dn_Math%re,dn_Math%im
   if (io.lt.0_isp) exit ! exit loop if end of file reached

   ! determine if SageMath and Mathematica values agree within tolerance
   CAS_test_array=[&
   &error(sn_SM%re,sn_Math%re),error(sn_SM%im,sn_Math%im),&
   &error(cn_SM%re,cn_Math%re),error(cn_SM%im,cn_Math%im),&
   &error(dn_SM%re,dn_Math%re),error(dn_SM%im,dn_Math%im)]

   CAS_test=all(CAS_test_array.lt.CAS_tol)

   if (CAS_test) then
    ! compute ellipFor values
    call Jacobi_elliptic_functions(u,m,sn,cn,dn)

    ! compute errors
    sn_err%re=error(sn%re,sn_SM%re)
    sn_err%im=error(sn%im,sn_SM%im)
    cn_err%re=error(cn%re,cn_SM%re)
    cn_err%im=error(cn%im,cn_SM%im)
    dn_err%re=error(dn%re,dn_SM%re)
    dn_err%im=error(dn%im,dn_SM%im)

    ! update max error values
    if (sn_err%re.gt.sn_err_max%re) sn_err_max%re=sn_err%re
    if (sn_err%im.gt.sn_err_max%im) sn_err_max%im=sn_err%im
    if (cn_err%re.gt.cn_err_max%re) cn_err_max%re=cn_err%re
    if (cn_err%im.gt.cn_err_max%im) cn_err_max%im=cn_err%im
    if (dn_err%re.gt.dn_err_max%re) dn_err_max%re=dn_err%re
    if (dn_err%im.gt.dn_err_max%im) dn_err_max%im=dn_err%im

    ! compute sum used for mean
    sn_err_mean=sn_err_mean+sn_err
    cn_err_mean=cn_err_mean+cn_err
    dn_err_mean=dn_err_mean+dn_err

    ! print error data to file
    write(202,'(9(g27.17))') u%re,u%im,m,sn_err%re,sn_err%im,cn_err%re,cn_err%im,dn_err%re,dn_err%im
    nargs=nargs+1
   end if ! if CAS agreement test passed
  end do
  close(101)
  close(202)

  ! finish mean calculation
  sn_err_mean=sn_err_mean/nargs
  cn_err_mean=cn_err_mean/nargs
  dn_err_mean=dn_err_mean/nargs

  print '(a7,i6,a54)', " Note: ",nargs," randomly selected combinations of u and m were tested"

  if (max(sn_err_max%re,sn_err_max%im) .lt. tol) then
   print '(a20)', " sn(u|m) test passed"
   sn_test_passed=.true.
  else
   print '(a43)', " ------ Warning: sn(u|m) test failed ------" 
   print '(a23,g25.17,a1,g25.17,a1)', "  Relative error max =(",sn_err_max %re,",",sn_err_max %im,")"
   print '(a23,g25.17,a1,g25.17,a1)', "  Relative error mean=(",sn_err_mean%re,",",sn_err_mean%im,")"
   sn_test_passed=.false.
  end if 
  if (max(cn_err_max%re,cn_err_max%im) .lt. tol) then
   print '(a20)', " cn(u|m) test passed"
   cn_test_passed=.true.
  else
   print '(a43)', " ------ Warning: cn(u|m) test failed ------" 
   print '(a23,g25.17,a1,g25.17,a1)', "  Relative error max =(",cn_err_max %re,",",cn_err_max %im,")"
   print '(a23,g25.17,a1,g25.17,a1)', "  Relative error mean=(",cn_err_mean%re,",",cn_err_mean%im,")"
   cn_test_passed=.false.
  end if 
  if (max(dn_err_max%re,dn_err_max%im) .lt. tol) then
   print '(a20)', " dn(u|m) test passed"
   dn_test_passed=.true.
  else
   print '(a43)', " ------ Warning: dn(u|m) test failed ------" 
   print '(a23,g25.17,a1,g25.17,a1)', "  Relative error max =(",dn_err_max %re,",",dn_err_max %im,")"
   print '(a23,g25.17,a1,g25.17,a1)', "  Relative error mean=(",dn_err_mean%re,",",dn_err_mean%im,")"
   dn_test_passed=.false.
  end if 
  if (all([sn_test_passed,cn_test_passed,dn_test_passed])) then
   test_passed=.true.
  else
   test_passed=.false.
  end if
  print '(a1)', " "

 end subroutine test_Jacobi_elliptic_functions

 function error(val,val_ref) result(e)
  ! *** Compute relative error based on reference value ***
  ! note: absolute error is computed if either input value is zero
  ! input
  real(dp),intent(in) :: val     ! value
  real(dp),intent(in) :: val_ref ! reference value

  ! output
  real(dp)            :: e       ! error

  if ((val_ref.ne.0._dp).and.(val.ne.0._dp)) then
   e=abs((val-val_ref)/val_ref) ! relative error
  else
   e=abs(val-val_ref) ! absolute error
  end if
 end function error
end program test_material_driver

