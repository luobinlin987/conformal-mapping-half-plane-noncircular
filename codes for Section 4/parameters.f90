Module parameters
  ! This module contains the parameters of the solution
  ! n  : the items of ck, Eq. (1)
  ! m  : the quantity of sample points, Eq. (6)
  ! h  : buried depth the elliptical cavity, Eq. (12)
  ! A  : half length of the major axis of the elliptical cavity, Eq. (12)
  ! B  : half length of the minor axis of the elliptical cavity, Eq. (12)
  ! pi : ratio of the circumference of a circle to the diameter
  ! OMEGA :  the \varOmega, Eq. (8)
  integer,parameter::n=15,m=180
  real(8),parameter::h=5,A=3,B=2
  real(8),parameter::pi=3.1415926535897932384626433832795
  real(8),parameter::OMEGA=1d10
End Module parameters

