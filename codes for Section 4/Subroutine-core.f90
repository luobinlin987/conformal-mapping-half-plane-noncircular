Subroutine core(F,alpha,C,theta)
  use parameters
  implicit none
  ! ******** variable introduction ********
  ! F    : F({\bm C},{\bm \theta}), Eq. (8)
  ! alpha: \alpha, Eq. (6)
  ! C    : c_{k} (k=0,1,2,3,\cdots,n), Eq. (6)
  ! theta: \theta_{i} (i=1,2,3,\cdots,m), Eq. (6)
  ! vt   : \phi_{i} (i=1,2,3,\cdots,m), Eq. (13)
  ! x0,y0: x_{i}^{0} and y_{i}^{0}, Eq. (13)
  ! x,y,xc,yc
  !      : X_{i} and Y_{i} (i=1,2,3,\cdots,m), Eq. (6)
  real(8),intent(out)::F
  real(8),intent(in)::alpha
  real(8),dimension(n+1),intent(in)::C
  real(8),dimension(m),intent(in)::theta
  real(8),dimension(m)::vt,x0,y0
  real(8),dimension(m)::x,y,xc,yc
  integer::i,k

  ! always setting intial values of the variable to zero for stablility
  x=0
  y=0
  xc=0
  yc=0
  vt=0
  x0=0
  y0=0
  F=0

  ! computing X_{i} and Y_{i}, Eq. (6)
  do i=1,m
     do k=1,n
        xc(i)=xc(i)+C(1+k)*(alpha**(k)+alpha**(-k))*sin(real(k)*theta(i))
        yc(i)=yc(i)+C(1+k)*(alpha**(k)-alpha**(-k))*cos(real(k)*theta(i))
     enddo
     x(i)=-C(1)*2*alpha*sin(theta(i))/(1d0+alpha**2-2d0*alpha*cos(theta(i)))-xc(i)
     y(i)=C(1)*(1d0-alpha**2)/(1d0+alpha**2-2d0*alpha*cos(theta(i)))+yc(i)
  enddo

  ! computing x_{i}^{0} and y_{i}^{0}, Eq. (13)
  do i=1,m
     vt(i)=pi*(i-1)/(m-1)-pi/2
     x0(i)=A*cos(vt(i))
     y0(i)=B*sin(vt(i))-H
  enddo

  ! computing the penalty function, Eq. (8)
  do i=2,m-1
     F=F+(x(i)-x0(i))**2
  enddo

  do i=1,m
     F=F+(y(i)-y0(i))**2
  enddo
  
  F=F+(max(0d0,alpha-1))**2+(max(0d0,-alpha))**2

  do i=1,m-1
     F=F+(max(0d0,theta(i)-theta(i+1)))**2
  enddo

  F=F*OMEGA

End Subroutine core
