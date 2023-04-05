include 'parameters.f90'
Program non_convex_verification
  use parameters
  implicit none
  ! ******** variable introduction ********
  ! theta: \phi_{i} (i=1,2,3,\cdots,m), Eq. (13)
  ! x0,y0: x_{i}^{0} and y_{i}^{0}, Eq. (13)
  ! C    : {\bm C}, Algorithm 1
  ! th   : {\bm \theta}, Algorithm 1
  ! x,y,xc,yc
  !      : X_{i} and Y_{i} (i=1,2,3,\cdots,m), Eq. (6)
  ! G    : G, Eq. (10)
  ! F    : F({\bm C},[:,1],{\bm \theta}[:,1]), 
  !        F({\bm C},[:,2],{\bm \theta}[:,2]),
  !        F({\bm C},[:,3],{\bm \theta}[:,3])
  ! t    : random numeric t within range [0,1), Lines 10,11,13,17, Algorithm 1
  ! s    : random numeric \lambda, Eqs. (10) and (11)
  ! num  : iteration rep, N, Algorithm 1
  real(8),dimension(m)::theta,x0,y0
  real(8),dimension(m,3)::th,x,y,xc,yc
  real(8),dimension(n+2,3)::C
  real(8)::G,F(3),t,s,temp
  integer::i,j,k,num
  integer::start,finish

  call system_clock(start)
  call random_seed()
  ! computing rectangular coordinates in Eq. (13), Line 1, Algorithm 1
  do i=1,m
     theta(i)=pi*(i-1)/(m-1)-pi/2
     x0(i)=A*cos(theta(i))
     y0(i)=B*sin(theta(i))-H
  enddo

  ! Lines 2 and 3, Algorithm 1
  G=1d0
  num=1
  ! starting iteration
  do while (G > 0)
  ! Line 6, Algorithm 1
     C=0
     th=0
     ! Lines 10-14, Algorithm 1
     do j=1,2
        call random_number(t)
        C(1,j)=t
        call random_number(t)
        C(2,j)=-(t)*2d1
        do i=3,n+2
           call random_number(t)
           C(i,j)=(t*2-1)*(1d1)**(2-i)
        enddo
     enddo
     ! Lines 15-19, Algorithm 1
     do j=1,2
        th(1,j)=0
        do i=2,m-1
           call random_number(t)
           th(i,j)=pi*t
        enddo
        th(m,j)=pi
     enddo
     ! Line 20, Algorithm 1
     do j=1,2
        do i=2,m
           do k=i,m
              if (abs(th(k,j)) < abs(th(i,j))) then
                 temp=th(i,j)
                 th(i,j)=th(k,j)
                 th(k,j)=temp
              endif
           enddo
        enddo
     enddo
     ! Lines 22-23, Algorithm 1
     call random_number(s)
     do i=1,n+2
        C(i,3)=s*C(i,1)+(1-s)*C(i,2)
     enddo
     do i=1,m
        th(i,3)=s*th(i,1)+(1-s)*th(i,2)
     enddo
     ! Line 24, Algorithm 1
     xc=0
     yc=0
     x=0
     y=0
     do j=1,3
        do i=1,m
           do k=3,n+2
              xc(i,j)=xc(i,j)+C(k,j)*(C(1,j)**(k-2)+C(1,j)**(-k+2))*sin((k-2)*th(i,j))
              yc(i,j)=yc(i,j)+C(k,j)*(C(1,j)**(k-2)-C(1,j)**(-k+2))*cos((k-2)*th(i,j))
           enddo
           x(i,j)=-2d0*C(2,j)*C(1,j)*sin(th(i,j))/(1+C(1,j)**2-2d0*C(1,j)*cos(th(i,j)))-xc(i,j)
           y(i,j)=C(2,j)*(1-C(1,j)**2)/(1+C(1,j)**2-2d0*C(1,j)*cos(th(i,j)))+yc(i,j)
        enddo
        x(1,j)=0
        x(m,j)=0
     enddo
     ! Line 25, Algorithm 1
     F=0
     do j=1,3
        do i=2,m-1
           F(j)=F(j)+(x(i,j)-x0(i))**2
        enddo
        do i=1,m
           F(j)=F(j)+(y(i,j)-y0(i))**2
        enddo
        F(j)=F(j)+(max(0d0,C(1,j)-1d0))**2+(max(0d0,-C(1,j)))**2
        do i=1,m-1
           F(j)=F(j)+(max(0d0,th(i,j)-th(i+1,j)))**2
        enddo
        F(j)=F(j)*OMEGA
     enddo
     G=0
     G=s*F(1)+(1-s)*F(2)-F(3)
     ! Line 26, Algorithm 1
     num=num+1
     write(*,*) num,',',G
  enddo

  call system_clock(finish)

  write(*,*) '*************************************'

  ! Lines 29-30, Algorithm 1
  write(*,*) 'Iteration reps =',num,',','G_minus =',G
  write(*,*) 'Elapsing time is',(real(finish)-real(start))*1d-3,'seconds'

End Program non_convex_verification
