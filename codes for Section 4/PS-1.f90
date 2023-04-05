Subroutine PS1(lgF,alpha0,c0,nummax,size)
  use parameters
  implicit none
  ! the five variables below are can be found in file 'main.f90'
  real(8),intent(in)::alpha0,c0
  integer,intent(in)::size,nummax
  real(8),intent(out)::lgF
  ! ******** variable introduction *********
  ! ******** the variables below are related to the algorithm
  ! 
  ! alhpa  : {\bm C}^{*}[1], Line 26, Algorithm 2
  ! C  : {\bm C}^{*}[2,n+1], Line 26, Algorithm 2
  ! theta  : {\bm \varTheta}^{*}, Line 27, Algorithm 2
  ! VTR  : value to reach, Input, Algorithm 2
  ! num  : current iteration rep, N, Lines 2 and 24, Algorithm 2
  ! winit,wend  : w_{\rm init} and w_{\rm end}, Input, Algorithm 2
  ! w  : w, Line 1, Algorithm 4
  ! d1init,d1end  : d_{1\rm init} and d_{1\rm end}, Input, Algorithm 2
  ! d2init,d2end  : d_{2\rm init} and d_{2\rm end}, Input, Algorithm 2
  ! d1,d2  : d_{1} and d_{2}, Lines 2 and 3, Algorithm 4
  ! r12  : r_{1} and r_{2}, Lines 2 and 18, Algorithm 5
  ! r    : r, Line 16, Algorithm 3
  ! temp : temporary data
  ! beta : \beta, Line 14, Algorithm 3
  ! valphamininit, valphaminend, valphamaxinit, valphamaxend 
  !      : velocity ranges of the upper and lower limits from initialization to the end
  !        during the iteration for \alpha,
  !        Input, Algorithm 2
  ! vCmininit, vCminend, vCmaxinit, vCmaxend
  !      : velocity ranges of the upper and lower limits from initialization to the end
  !        during the iteration for c_{k} (k=0,1,2,3,\cdots,n),
  !        Input, Algorithm 2
  ! vthetamininit, vthetaminend, vthetamaxinit, vthetamaxend
  !      : velocity ranges of the upper and lower limits from initialization to the end
  !        during the iteration for \theta_{i} (i=1,2,3,\cdots,m),
  !        Input, Algorithm 2
  ! valphamin,valphamax
  !      : ranges of the upper and lower limits for the velocity of \alpha,
  !        Lines 5 and 6, Algorithm 4
  ! alphamin,alphamax
  !      : ranges of the upper and lower limits for \alpha
  !        Eqs. (14) and (15)
  ! vCmin, vCmax
  !      : ranges of the upper and lower limits for the velocity of c_{k} (k=1,2,3,\cdots,n),
  !        Lines 5 and 6, Algorithm 4
  ! Cmin, Cmax
  !      : ranges of the upper and lower limits for c_{k} (k=1,2,3,\cdots,n),
  !        Eqs. (14) and (15)
  ! vthetamin, vthetamax
  !      : ranges of the upper and lower limits for the velocity of \theta_{i} (i=1,2,3,\cdots,m),
  !        Lines 9 and 10, Algorithm 4
  ! thetamin, thetamax
  !      : ranges of the upper and lower limits for \theta_{i} (i=1,2,3,\cdots,m),
  !        between 0 and pi
  ! valpha, xalpha
  !      : velocity and position for \alpha, respectively,
  !        Lines 2 and 10, respectively
  ! vC, xC
  !      : velocity and position for c_{k} (k=1,2,3,\cdots,n), respectively
  !        Lines 2 and 10, respectively
  ! vtheta, xtheta
  !      : velocity and position for theta_{i} (i=1,2,3,\cdots,m), respectively,
  !        Lines 18 and 26, respectively
  ! palpha, pC, ptheta
  !      : individual best position for \alpha, c_{k} (k=1,2,3,\cdots,n),
  !        and \theta_{i} (i=1,2,3,\cdots,m), respectively,
  !        Lines 12-15, Algorithm 2
  ! pgalpha, pgC, pgtheta
  !      : global best postion for \alpha, c_{k} (k=1,2,3,\cdots,n),
  !        and \theta_{i} (i=1,2,3,\cdots,m), respectively,
  !        Lines 16-20, Algorithm 2
  ! F1, F2
  !      : F({\bm C}[:,j],{\bm \theta}[:,j]),
  !        and F({\bm C}_{\rm ibest}[:,j],{\bm \theta}_{\rm ibest}[:,j]), respectively,
  !        Lines 8-10, Algorithm 2
  ! Fbest: F_{N}({\bm C}_{\rm gbest},{\bm \varTheta}_{\rm gbest}),
  !        and F_{N+1}({\bm C}_{\rm gbest},{\bm \varTheta}_{\rm gbest}),
  !        Lines 3, 4, 5 and 23, Algorithm 2
  !******** variable introduction finished ********
  real(8)::alpha
  real(8),dimension(n+1)::C
  real(8),dimension(m)::theta
  real(8)::VTR
  integer::num
  integer::i,j,k
  real(8)::winit,wend,w
  real(8)::d1,d2,d1init,d1end,d2init,d2end
  real(8)::r12(2),r,temp,beta
  real(8)::valphamininit,valphaminend,valphamaxinit,valphamaxend
  real(8),allocatable,dimension(:)::vCmininit,vCminend,vCmaxinit,vCmaxend
  real(8),allocatable,dimension(:)::vthetamininit,vthetaminend,vthetamaxinit,vthetamaxend
  real(8)::alphamin,alphamax,valphamin,valphamax
  real(8),allocatable,dimension(:)::vthetamin,vthetamax,vCmin,vCmax
  real(8),allocatable,dimension(:)::Cmin,Cmax,thetamin,thetamax
  real(8),allocatable,dimension(:)::xalpha,valpha
  real(8),allocatable,dimension(:,:)::xC,vC,xtheta,vtheta
  real(8),allocatable,dimension(:)::palpha
  real(8),allocatable,dimension(:,:)::pC,ptheta
  real(8)::pgalpha
  real(8),allocatable,dimension(:)::pgC,pgtheta
  real(8),allocatable,dimension(:)::F1,F2
  real(8),dimension(2)::Fbest
  ! ******** variables below are not related to algorithm
  ! x, y, xcc, ycc
  !       : rectangular coordinates of cavity boundary in each iteration,
  !         which are recorded for data visualization to monitor the iteration trends
  ! ********
  real(8),allocatable,dimension(:)::x,y,xcc,ycc
  ! ******** path for data saving
  character(*),parameter::datapath=&
       '/home/luobinlin/researchpaper/paper3/code/data/particle-swarm/PS1/'
  character(5)::var

  open(20000,file=datapath//'20000.csv')

  ! storage allocation
  allocate(vCmininit(n+1),vCminend(n+1),vCmaxinit(n+1),vCmaxend(n+1))
  allocate(vthetamininit(m),vthetaminend(m),vthetamaxinit(m),vthetamaxend(m))
  allocate(vCmin(n+1),vCmax(n+1),vthetamin(m),vthetamax(m))
  allocate(Cmin(n+1),Cmax(n+1),thetamin(m),thetamax(m))
  allocate(xalpha(size),valpha(size))
  allocate(xC(n+1,size),vC(n+1,size))
  allocate(xtheta(m,size),vtheta(m,size))
  allocate(palpha(size),pC(n+1,size),ptheta(m,size),F1(size),F2(size))
  allocate(pgC(n+1),pgtheta(m))
  allocate(x(m),y(m),xcc(m),ycc(m))

  ! value assignment of the ranges of variable speeds, Eqs. (17) and (18)
  valphamininit=-0.02
  valphaminend=-0.01
  valphamaxinit=-valphamininit
  valphamaxend=-valphaminend

  vCmininit=-0.02
  vCminend=-0.01
  vCmaxinit=-vCmininit
  vCmaxend=-vCminend

  vthetamininit=-0.01
  vthetaminend=-0.005
  vthetamaxinit=-vthetamininit
  vthetamaxend=-vthetaminend

  ! value assignment of variable ranges, Eqs. (14) and (15)
  alphamin=alpha0
  alphamax=1d0
  Cmin(1)=2d0*c0
  Cmax(1)=0d0
  do i=2,n+1
     Cmin(i)=-1d0*exp(real(2-i))
     Cmax(i)=-Cmin(i)
  enddo

  ! value assigment of the range of \theta
  thetamin=0
  thetamax=pi

  ! parameters of modified Particle Swarm Method, Eq. (19)
  winit=0.7
  wend=0.4

  d1init=2.5
  d1end=0.5
  d2init=0.9
  d2end=2.25

  ! value of \beta, Line 14, Algorithm 3 and Eq. (20)
  beta=1d-2

  ! static state of modified Particle Swarm Method, Lines 2-7, Algorithm 3
  valpha=0
  vC=0
  vtheta=0

  ! initial values of variable position
  call random_seed()
  do j=1,size
     ! Lines 10-12, Algorithm 3
     call random_number(r)
     xalpha(j)=alphamin+r*(alphamax-alphamin)
     !
     do i=1,n+1
        call random_number(r)
        xC(i,j)=0+(Cmin(i)+r*(Cmax(i)-Cmin(i)))
     enddo
     ! Line 16, Algorithm 3, activated for Sets b and d
     ! do i=1,m
     !    call random_number(r)
     !    xtheta(i,j)=r*pi
     ! enddo
     ! Line 14, Algorithm 3, activated for Sets a and c
     do i=1,m
        call random_number(r)
        xtheta(i,j)=pi*(i-1)/(m-1)+pi/(m-1)*beta*(2d0*r-1d0)
     enddo
     ! Line 18, Algorithm 3
     xtheta(1,j)=0
     xtheta(m,j)=pi
     ! Line 19, Algorithm 3, activated for Sets a and b
     do i=1,m
        do k=i,m
           if (xtheta(i,j) > xtheta(k,j)) then
              temp=xtheta(i,j)
              xtheta(i,j)=xtheta(k,j)
              xtheta(k,j)=temp
           endif
        enddo
     enddo
  enddo

  ! iteration begins
  palpha=xalpha
  pC=xC
  ptheta=xtheta
  pgalpha=OMEGA
  pgC=OMEGA
  pgtheta=OMEGA

  num=0
  VTR=1d-6
  Fbest=(/OMEGA**3,OMEGA**2/)
  do while (Fbest(1)-Fbest(2) >= VTR .or. num <= nummax)
     write(var,'(I5)') 10000+num
     open(10000+num, file=datapath//var//'.csv')
     ! Line 7, Algorithm 2
     Fbest(1)=Fbest(2)
     ! Line 1-11, Algorithm 4
     w=winit+(wend-winit)*num/nummax
     d1=d1init+(d1end-d1init)*num/nummax
     d2=d2init+(d2end-d2init)*num/nummax
     valphamin=valphamininit+(valphaminend-valphamininit)*num/nummax
     valphamax=valphamaxinit+(valphamaxend-valphamaxinit)*num/nummax
     vCmin=vCmininit+(vCminend-vCmininit)*num/nummax
     vCmax=vCmaxinit+(vCmaxend-vCmaxinit)*num/nummax
     vthetamin=vthetamininit+(vthetaminend-vthetamininit)*num/nummax
     vthetamax=vthetamaxinit+(vthetamaxend-vthetamaxinit)*num/nummax

     ! Lines 8-10, Algorithm 2, parallel computation accelarated via OPENACC and CUDA
     !$acc data copyin(xalpha,palpha,xC,xtheta,pC,ptheta), copyout(F1,F2)
     !$acc parallel loop 
     do j=1,size
        call core(F1(j),xalpha(j),xC(:,j),xtheta(:,j))
        call core(F2(j),palpha(j),pC(:,j),ptheta(:,j))
     enddo
     !$acc end parallel loop
     !$acc end data

     ! Line 11-20, Algorithm 2
     do j=1,size
        if (F1(j) < F2(j)) then
           palpha(j)=xalpha(j)
           pC(:,j)=xC(:,j)
           ptheta(:,j)=xtheta(:,j)
        endif
        if (F2(j) < Fbest(1)) then
           Fbest(2)=F2(j)
           pgalpha=palpha(j)
           pgC=pC(:,j)
           pgtheta=ptheta(:,j)
        endif
        ! Line 1-32, Algorithm 4
        call random_number(r12)
        valpha(j)=w*valpha(j)&
             +d1*r12(1)*(palpha(j)-xalpha(j))&
             +d2*r12(2)*(pgalpha-xalpha(j))
        xalpha(j)=xalpha(j)+valpha(j)
        if (valpha(j) > valphamax) then
           valpha(j)=valphamax
        else if (valpha(j) < valphamin) then
           valpha(j)=valphamin
        endif
        if (xalpha(j) > alphamax) then
           xalpha(j)=alphamax
        else if (xalpha(j) < alphamin) then
           xalpha(j)=alphamin
        endif
        do i=1,n+1
           call random_number(r12)
           vC(i,j)=w*vC(i,j)&
                +d1*r12(1)*(pC(i,j)-xC(i,j))&
                +d2*r12(2)*(pgC(i)-xC(i,j))
           xC(i,j)=xC(i,j)+vC(i,j)
        enddo
        do i=1,n+1
           if (vC(i,j) > vCmax(i)) then
              vC(i,j)=vCmax(i)
           else if (vC(i,j) < vCmin(i)) then
              vC(i,j)=vCmin(i)
           endif
        enddo
        do i=1,n+1
           if (xC(i,j) < Cmin(i)) then
              xC(i,j)=Cmin(i)
           else if (xC(i,j) > Cmax(i)) then
              xC(i,j)=Cmax(i)
           endif
        enddo
        do i=1,m
           call random_number(r12)
           vtheta(i,j)=w*vtheta(i,j)&
                +d1*r12(1)*(ptheta(i,j)-xtheta(i,j))&
                +d2*r12(2)*(pgtheta(i)-xtheta(i,j))
           xtheta(i,j)=xtheta(i,j)+vtheta(i,j)
        enddo
        do i=1,m
           if (vtheta(i,j) < vthetamin(i)) then
              vtheta(i,j)=vthetamin(i)
           else if (vtheta(i,j) > vthetamax(i)) then
              vtheta(i,j)=vthetamax(i)
           endif
        enddo
        do i=1,m
           if (xtheta(i,j) < thetamin(i)) then
              xtheta(i,j) = thetamin(i)
           else if (xtheta(i,j) > thetamax(i)) then
              xtheta(i,j) = thetamax(i)
           endif
        enddo
        ! Line 34, Algorithm 5, activated in Sets a and b
        do i=1,m
           do k=i,m
              if (xtheta(k,j) < xtheta(i,j)) then
                 temp=xtheta(i,j)
                 xtheta(i,j)=xtheta(k,j)
                 xtheta(k,j)=temp
              endif
           enddo
        enddo
        ! Line 33, Algorithm 5
        xtheta(1,j)=0
        xtheta(m,j)=pi
     enddo
     ! rectangular coordinate computation for data visulization and convergence monitor
     x=0
     y=0
     xcc=0
     ycc=0
     do i=1,m
        do k=1,n
           xcc(i)=xcc(i)+pgC(1+k)*(pgalpha**(k)+pgalpha**(-k))*sin(real(k)*pgtheta(i))
           ycc(i)=ycc(i)+pgC(1+k)*(pgalpha**(k)-pgalpha**(-k))*cos(real(k)*pgtheta(i))
        enddo
        x(i)=-pgC(1)*2d0*pgalpha*sin(pgtheta(i))/(1d0+pgalpha**2-2d0*pgalpha*cos(pgtheta(i)))-xcc(i)
        y(i)=pgC(1)*(1d0-pgalpha**2)/(1d0+pgalpha**2-2d0*pgalpha*cos(pgtheta(i)))+ycc(i)
     enddo
     do i=1,m
        write(10000+num,1001) x(i),',',y(i)
     enddo
1001 format(F20.15,A,F20.15)
     close(10000+num)
     write(20000,2001) num,',',log10(Fbest(2))
2001 format(I4,A,F30.10)
     ! Line 23, Algorithm 2
     write(*,*) 'PS-1',',','N=',num,',','Fbest=',Fbest(2)
     ! Line 24, Algorithm 2
     num=num+1
  enddo
  ! Lines 26 and 27, Algorithm 2
  alpha=pgalpha
  C=pgC
  theta=pgtheta

  lgF=log10(Fbest(2))

  ! final rectangular coordinates for data visulization in Fig. 2
  open(1001,file=datapath//'1001.csv')
  xcc=0
  ycc=0
  x=0
  y=0
  do i=1,m
     do k=1,n
        xcc(i)=xcc(i)+C(1+k)*(alpha**(k)+alpha**(-k))*sin(real(k)*theta(i))
        ycc(i)=ycc(i)+C(1+k)*(alpha**(k)-alpha**(-k))*cos(real(k)*theta(i))
     enddo
     x(i)=-C(1)*2*alpha*sin(theta(i))/(1d0+alpha**2-2d0*alpha*cos(theta(i)))-xcc(i)
     y(i)=C(1)*(1d0-alpha**2)/(1d0+alpha**2-2d0*alpha*cos(theta(i)))+ycc(i)
  enddo

  do i=1,m
     write(1001,3001) x(i),',',y(i)
  enddo
3001 format(F20.15,A,F20.15)
  close(1001)

  write(*,*) '********************************'
  write(*,*) 'alpha0=',alpha0,'c0=',c0
  write(*,*) 'RESULT'
  write(*,*) 'Iteration reps =',num,',','Fbest =',Fbest(2)
  ! Line 28, Algorithm 2
  do i=1,n+1
     write(*,*) C(i)
  enddo

End Subroutine PS1
