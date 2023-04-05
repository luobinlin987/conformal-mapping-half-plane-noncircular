! distributed encoding ideology is applied, thus, other .f90 files are included
  include 'parameters.f90'
  include 'Subroutine-core.f90'
  include 'PS-1.f90'
  include 'PS-2.f90'
  include 'PS-3.f90'
  include 'PS-4.f90'
  include 'PS-final.f90'
  ! 'parameters.f90' contains the parameters for computation
  ! 'Subroutine-core.f90' is the core subroutine to compute the penalty function, Eq. (8)
  ! 'PS-i.f90' (i=1,2,3,4) are subroutines to the four sets of strategy combiniations, Table 1
  ! 'PS-final.f90' is the subroutine for the discussion of solution and efficiency

  Program main
    use parameters
    implicit none
    ! variable description
    ! vt  : values of phi_{i}, Eq. (13)
    ! x0,y0  : values of x_{i}^{0} and y_{i}^{0}, Eq. (13)
    ! r0  : the distance between arbitrary point and cavity center, Eq. (16)
    ! r0min  : the radius of incribed circle of the cavity, Eq. (16)
    ! alpha0 : {\bm C}_{\min}[1], Eq. (14)
    ! c0  : {\bm C}_{\min}[2], Eq. (14)
    ! lgF : final value of penalty function for discussion of solution accuracy and efficiency
    !       in Section 4.4
    ! lgFi (i=1,2,3,4)  : final value of penalty functions for the four sets of strategy
    !                     combinations
    ! nummax : N_{\max} of the Input parameters in Algorithm 2
    ! size : the population size of the Particle Swarm Method,
    !        S of the Input parameters in Algorithm 2
    ! 
    real(8),dimension(m)::vt,x0,y0
    real(8),dimension(m)::r0
    real(8)::r0min
    integer::i,k
    real(8)::alpha0,c0,lgF,lgF1,lgF2,lgF3,lgF4
    integer::nummax,size

    ! the variables below are related to computing ellapsing time of the solutions
    character(10)::start,finish
    character(10)::start1,start2,start3,start4,start5,startfinal
    character(10)::finish1,finish2,finish3,finish4,finish5,finishfinal
    integer::start1hrs,finish1hrs,start1mins,finish1mins
    real(8)::start1secs,finish1secs
    integer::start2hrs,finish2hrs,start2mins,finish2mins
    real(8)::start2secs,finish2secs
    integer::start3hrs,finish3hrs,start3mins,finish3mins
    real(8)::start3secs,finish3secs
    integer::start4hrs,finish4hrs,start4mins,finish4mins
    real(8)::start4secs,finish4secs
    integer::starthrs,finishhrs,startmins,finishmins
    real(8)::startsecs,finishsecs
    real(8)::elapsingtime

    ! the variable below is related to data storage
    character(*),parameter::datapath=&
         '/home/luobinlin/researchpaper/paper3/code/data/'

    
    call date_and_time(TIME=start)

    open(1000,file=datapath//'1000.csv')

    ! computing the coordinates of the right half boundary of the elliptical cavity, Eq. (13)
    do i=1,m
       vt(i)=pi*(i-1)/(m-1)-pi/2
       x0(i)=A*cos(vt(i))
       y0(i)=B*sin(vt(i))-H
    enddo

    do i=1,m
       write(1000,1001) x0(i),',',y0(i)
    enddo
1001 format(F20.15,A,F20.15)

    ! computing the distance between the sample points and cavity center
    do i=1,m
       r0(i)=sqrt((x0(i)-0)**2+(y0(i)+H)**2)
    enddo

    ! selecting the minimum distance,
    ! in other words, the radius of the incribed circle of the cavity, Eq. (16)
    r0min=minval(r0)

    ! Eq. (14)
    alpha0=r0min/(H+sqrt(H**2-r0min**2))
    c0=-H*(1d0-alpha0**2)/(1d0+alpha0**2)

    ! Particle Swarm Method begins
    write(*,*) 'Particle Swarm Method'

    ! maximum iteration rep and population size
    nummax=3000
    size=10*(n+2+m)

    ! computing the four sets of strategy combinations
    ! alpha0, c0, nummax, size are sent into the subroutine PSi (i=1,2,3,4)
    ! while lgFi (i=1,2,3,4) are recieved as results from the subroutines
    call date_and_time(TIME=start1)
    call PS1(lgF1,alpha0,c0,nummax,size)
    call date_and_time(TIME=finish1)

    call date_and_time(TIME=start2)
    call PS2(lgF2,alpha0,c0,nummax,size)
    call date_and_time(TIME=finish2)

    call date_and_time(TIME=start3)
    call PS3(lgF3,alpha0,c0,nummax,size)
    call date_and_time(TIME=finish3)

    call date_and_time(TIME=start4)
    call PS4(lgF4,alpha0,c0,nummax,size)
    call date_and_time(TIME=finish4)

    ! computing ellpasing time for the four sets
    write(*,*) 'start time 1 : ',start1,',','finish time 1 : ',finish1
    read(start1(1:2),*) start1hrs
    read(start1(3:4),*) start1mins
    read(start1(5:10),*) start1secs
    read(finish1(1:2),*) finish1hrs
    read(finish1(3:4),*) finish1mins
    read(finish1(5:10),*) finish1secs
    write(*,*) 'Lg(F1) = ', lgF1
    write(*,*) 'Elapsing time of set 1 = ', 3600.0*(finish1hrs-start1hrs)&
         +60.0*(finish1mins-start1mins)&
         +(finish1secs-start1secs), 'seconds'

    write(*,*) 'start time 2 : ',start2,',','finish time 2 : ',finish2
    read(start2(1:2),*) start2hrs
    read(start2(3:4),*) start2mins
    read(start2(5:10),*) start2secs
    read(finish2(1:2),*) finish2hrs
    read(finish2(3:4),*) finish2mins
    read(finish2(5:10),*) finish2secs
    write(*,*) 'Lg(F2) = ', lgF2
    write(*,*) 'Elapsing time of set 2 = ', 3600.0*(finish2hrs-start2hrs)&
         +60.0*(finish2mins-start2mins)&
         +(finish2secs-start2secs), 'seconds'

    write(*,*) 'start time 3 : ',start3,',','finish time 3 : ',finish3
    read(start3(1:2),*) start3hrs
    read(start3(3:4),*) start3mins
    read(start3(5:10),*) start3secs
    read(finish3(1:2),*) finish3hrs
    read(finish3(3:4),*) finish3mins
    read(finish3(5:10),*) finish3secs
    write(*,*) 'Lg(F3) = ', lgF3
    write(*,*) 'Elapsing time of set 3 = ', 3600.0*(finish3hrs-start3hrs)&
         +60.0*(finish3mins-start3mins)&
         +(finish3secs-start3secs), 'seconds'

    write(*,*) 'start time 4 : ',start4,',','finish time 4 : ',finish4
    read(start4(1:2),*) start4hrs
    read(start4(3:4),*) start4mins
    read(start4(5:10),*) start4secs
    read(finish4(1:2),*) finish4hrs
    read(finish4(3:4),*) finish4mins
    read(finish4(5:10),*) finish4secs
    write(*,*) 'Lg(F4) = ', lgF4
    write(*,*) 'Elapsing time of set 4 = ', 3600.0*(finish4hrs-start4hrs)&
         +60.0*(finish4mins-start4mins)&
         +(finish4secs-start4secs), 'seconds'

    ! ******************solution procedure for the four sets are finished



    ! ******************solution procedure for the discussion of solution accuracy and efficiency
    !                   are below

    ! call date_and_time(TIME=startfinal)
    ! call PS(lgF,alpha0,c0,nummax,size)
    ! call date_and_time(TIME=finishfinal)

    ! write(*,*) 'start time final : ',startfinal,',','finish time final : ',finishfinal

    
    ! write(*,*) 'Lg(F) = ',lgF

    
    ! write(*,*) 'Particle Swarm Method finished'

    ! call date_and_time(TIME=finish)

    ! write(*,*) 'start time 0 : ',start, ',','finish time 0 : ',finish

    ! read(start(1:2),*) starthrs
    ! read(start(3:4),*) startmins
    ! read(start(5:10),*) startsecs
    ! read(finish(1:2),*) finishhrs
    ! read(finish(3:4),*) finishmins
    ! read(finish(5:10),*) finishsecs

    ! write(*,*) 'Elapsing time of all sets = ', 3600.0*(finishhrs-starthrs)&
    !      +60.0*(finishmins-startmins)&
    !      +(finishsecs-startsecs), 'seconds'

  End Program main
  
