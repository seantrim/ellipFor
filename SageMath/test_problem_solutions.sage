#SageMath script for generating reference CAS values for the test_material_driver.f90 program
m_min=0.01                # minimum
m_max=100.01              # maximum
N_m=100                   # number of m parameter values
delta_m=(m_max-m_min)/N_m # interval between m values
phi=3*pi/4                # argument for incomplete Legendre elliptic integrals
u=1+I                     # argument for Jacobi elliptic functions

print("Computing Data for Complete Elliptic Integrals")
with open('CAS_complete.dat', 'w') as f: 
 for i in range(1,N_m+2):
  m=m_min+delta_m*(i-1)
  Fc=elliptic_kc(m).n()
  Ec=elliptic_ec(m).n()
  print(m,Fc.real(),Fc.imag(),Ec.real(),Ec.imag(), file=f)

print("Computing Data for Incomplete Elliptic Integrals")
with open('CAS_incomplete.dat', 'w') as f:
 for i in range(1,N_m+2):
  m=m_min+delta_m*(i-1)
  F=elliptic_f(phi,m).n()
  E=elliptic_e(phi,m).n()
  print(m,F.real(),F.imag(),E.real(),E.imag(), file=f)

print("Computing Data for Jacobi Elliptic Functions")
with open('CAS_functions.dat', 'w') as f:
 for i in range(1,N_m+2):
  m=m_min+delta_m*(i-1)
  sn=jacobi_sn(u,m).n()
  cn=jacobi_cn(u,m).n()
  dn=jacobi_dn(u,m).n()
  print(m,sn.real(),sn.imag(),cn.real(),cn.imag(),dn.real(),dn.imag(), file=f)
