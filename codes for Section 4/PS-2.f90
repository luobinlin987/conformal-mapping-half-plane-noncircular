Subroutine PS2(lgF,alpha0,c0,nummax,size)
  use parameters
  implicit none
  real(8),intent(in)::alpha0,c0
  integer,intent(in)::size,nummax
  real(8),intent(out)::lgF
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
  real(8),allocatable,dimension(:)::x,y,xcc,ycc
  real(8),dimension(2)::Fbest
  character(*),parameter::datapath=&
       '/home/luobinlin/researchpaper/paper3/code/data/particle-swarm/PS2/'
  character(5)::var

  open(20000,file=datapath//'20000.csv')

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

  alphamin=alpha0
  alphamax=1d0
  Cmin(1)=2d0*c0
  Cmax(1)=0d0
  do i=2,n+1
     Cmin(i)=-1d0*exp(real(2-i))
     Cmax(i)=-Cmin(i)
  enddo

  thetamin=0
  thetamax=pi

  winit=0.7
  wend=0.4

  d1init=2.5
  d1end=0.5
  d2init=0.9
  d2end=2.25

  beta=1d-2

  valpha=0
  vC=0
  vtheta=0

  call random_seed()
  do j=1,size
     call random_number(r)
     xalpha(j)=alphamin+r*(alphamax-alphamin)
     do i=1,n+2
        call random_number(r)
        xC(i,j)=0+(Cmin(i)+r*(Cmax(i)-Cmin(i)))
     enddo
     do i=1,m
        call random_number(r)
        xtheta(i,j)=r*pi
     enddo
     ! do i=1,m
     !    call random_number(r)
     !    xtheta(i,j)=pi*(i-1)/(m-1)+pi/(m-1)*beta*(2d0*r-1d0)
     ! enddo
     xtheta(1,j)=0
     xtheta(m,j)=pi
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
     Fbest(1)=Fbest(2)
     w=winit+(wend-winit)*num/nummax
     d1=d1init+(d1end-d1init)*num/nummax
     d2=d2init+(d2end-d2init)*num/nummax
     valphamin=valphamininit+(valphaminend-valphamininit)*num/nummax
     valphamax=valphamaxinit+(valphamaxend-valphamaxinit)*num/nummax
     vCmin=vCmininit+(vCminend-vCmininit)*num/nummax
     vCmax=vCmaxinit+(vCmaxend-vCmaxinit)*num/nummax
     vthetamin=vthetamininit+(vthetaminend-vthetamininit)*num/nummax
     vthetamax=vthetamaxinit+(vthetamaxend-vthetamaxinit)*num/nummax

     !$acc data copyin(xalpha,palpha,xC,xtheta,pC,ptheta), copyout(F1,F2)
     !$acc parallel loop 
     do j=1,size
        call core(F1(j),xalpha(j),xC(:,j),xtheta(:,j))
        call core(F2(j),palpha(j),pC(:,j),ptheta(:,j))
     enddo
     !$acc end parallel loop
     !$acc end data

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
        do i=1,m
           do k=i,m
              if (xtheta(k,j) < xtheta(i,j)) then
                 temp=xtheta(i,j)
                 xtheta(i,j)=xtheta(k,j)
                 xtheta(k,j)=temp
              endif
           enddo
        enddo
        xtheta(1,j)=0
        xtheta(m,j)=pi
     enddo
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
     write(*,*) 'PS-2',',','N=',num,',','Fbest=',Fbest(2)
     num=num+1
  enddo
  alpha=pgalpha
  C=pgC
  theta=pgtheta

  lgF=log10(Fbest(2))

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
  do i=1,n+1
     write(*,*) C(i)
  enddo

End Subroutine PS2
