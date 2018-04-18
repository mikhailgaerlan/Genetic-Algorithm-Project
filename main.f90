program main
  use rands
  use convert
  use gen_tools
  implicit none
  
  integer :: popsize,indsize,objectives,ngen
  logical :: prints,seeded
  real(kind=8) :: cxpb,distance,unifpb
  real(kind=8) :: mutpb,flipprob,mu,sigma
  real(kind=8),dimension(:,:),allocatable :: pop,fitness,offspring,offfitness,total,totalfitness
  real(kind=8) :: rand
  integer :: i,j,k,n,prntrows
  character(len=32) :: numstr,prntformat
  character(len=32) :: funcname,crossover,mutation
  namelist /parameters/prints,seeded,popsize,ngen,funcname,&
       crossover,cxpb,distance,unifpb,mutation,mutpb,flipprob,mu,sigma
  
  !=============================
  !  Parameter Initialization
  !=============================
  open(10,file='parameters.txt',status='old')
  read(10,nml=parameters)

  if (seeded) call set_seed
  
  if (funcname.eq.'ackley') then
     indsize = 2
     objectives = 1
  elseif (funcname.eq.'binhkorn') then
     indsize = 2
     objectives = 2
  elseif (funcname.eq.'kursawe') then
     indsize = 3
     objectives = 2
  endif

  allocate(pop(popsize,indsize),fitness(popsize,objectives+3))
  allocate(offspring(popsize,indsize),offfitness(popsize,objectives+3))
  allocate(total(2*popsize,indsize),totalfitness(2*popsize,objectives+3))
  prntrows = indsize+objectives+3
  write(numstr,'(i4)') prntrows
  prntformat = '('//trim(adjustl(numstr))//'f10.4)'
  
  !=============================
  !  Population Initialization
  !=============================
  do j=1,indsize
     do i=1,popsize
        call random_number(rand)
        pop(i,j) = 10.0d0*(rand-0.5d0)
     enddo
  enddo
  do i=1,popsize
     fitness(i,1) = f1(pop(i,:))
     fitness(i,2) = f2(pop(i,:))
     fitness(i,objectives+3) = 1.0d0
  enddo
  
  call rank(fitness)
  call distancing(fitness)
  call sort(pop,fitness)
  do i=1,popsize
     write(*,prntformat) pop(i,:),fitness(i,:)
  enddo
  !write(*,*) ''

  !=============================
  !      Begin Algorithm
  !=============================
  do n=1,ngen
     offspring = pop
     offfitness = fitness
     
     !=============================
     !    Crossover Operation
     !=============================
     if (prints) write(*,*) "Crossing..."
     call distancing(offfitness)
     do i=1,popsize
        do j=1,popsize
           call random_number(rand)
           if ((rand<cxpb).and.(abs(offfitness(i,objectives+1)-offfitness(j,objectives+1)).ge.distance)) then
              do k=1,indsize
                 if (crossover.eq.'oneptcross') call oneptcross(offspring(i,k),offspring(j,k),offspring(i,k),offspring(j,k))
                 if (crossover.eq.'twoptcross') call twoptcross(offspring(i,k),offspring(j,k),offspring(i,k),offspring(j,k))
                 if (crossover.eq.'unifcross') call unifcross(offspring(i,k),offspring(j,k),offspring(i,k),offspring(j,k),unifpb)
                 if (offspring(i,k)<-5.0d0) offspring(i,k)=-5.0d0
                 if (offspring(i,k)>5.0d0) offspring(i,k)=5.0d0
                 if (offspring(j,k)<-5.0d0) offspring(j,k)=-5.0d0
                 if (offspring(j,k)>5.0d0) offspring(j,k)=5.0d0
              enddo
              offfitness(i,4) = 0.0d0
              offfitness(j,4) = 0.0d0
           endif
        enddo
     enddo
     if (prints) write(*,*) "Crossing done."
     
     !=============================
     !           Mutation
     !=============================
     do i=1,popsize
        call random_number(rand)
        if (rand<mutpb) then
           do k=1,indsize
              if (mutation.eq.'mutflipbit') call mutgauss(offspring(i,k),mu,sigma)
              if (mutation.eq.'mutgauss') call mutflipbit(offspring(i,k),flipprob)
              if (offspring(i,k)<-5.0d0) offspring(i,k)=-5.0d0
              if (offspring(i,k)>5.0d0) offspring(i,k)=5.0d0
           enddo
           offfitness(i,4) = 0.0d0
        endif
     enddo
     
     !=============================
     !     Calculate Fitness
     !=============================
     if (prints) write(*,*) "Evaluating fitnesses..."
     do i=1,popsize
        if (offfitness(i,4).eq.0.0d0) then
           offfitness(i,1) = f1(offspring(i,:))
           offfitness(i,2) = f2(offspring(i,:))
           offfitness(i,4) = 1.0d0
        endif
     enddo
     if (prints) write(*,*) "Evaluating done."
     
     !=============================
     !          Sorting
     !=============================
     total(:popsize,:) = pop
     totalfitness(:popsize,:) = fitness
     total(popsize+1:,:) = offspring
     totalfitness(popsize+1:,:) = offfitness
     if (prints) write(*,*) "Ranking..."
     call rank(totalfitness)
     if (prints) write(*,*) "Ranking done"
     if (prints) write(*,*) "Sorting..."
     call dominsort(total,totalfitness)
     if (prints) write(*,*) "Sortind done."
     pop = total(:popsize,:)
     fitness = totalfitness(:popsize,:)
     if (prints) write(*,*) ""
     do i=1,popsize
        write(*,prntformat) pop(i,:),fitness(i,:)
     enddo
  enddo
  
contains

  function f1(x)
    real(kind=8) :: f1
    real(kind=8),intent(in) :: x(:)
    real(kind=8),parameter :: pi=4.0d0*atan(1.0d0)
    integer :: i
    
    f1=0.0d0
    if (funcname.eq.'ackley') then
       f1 = -20.0d0*exp(-0.2d0*sqrt(0.5d0*(x(1)**2.0d0+x(2)**2.0d0)))-&
            exp(0.5d0*(cos(2*pi*x(1))+cos(2*pi*x(2))))+exp(1.0d0)+20
    elseif (funcname.eq.'binhkorn') then
       f1 = 4.0d0*(x(1)**2.0d0)+4.0d0*(x(2)**2.0d0) !Binh-Korn 1
    elseif (funcname.eq.'kursawe') then
       do i=1,2
          f1 = f1-10.0d0*exp(-0.2d0*sqrt(x(i)**2.0d0+x(i+1)**2.0d0)) !Kursawe 1
       enddo
    endif
  end function f1

  function f2(x)
    real(kind=8) :: f2
    real(kind=8),intent(in) :: x(:)
    integer :: i
    
    f2=0.0d0
    if (funcname.eq.'binhkorn') then
       f2 = (x(1)-5.0d0)**2.0d0+(x(2)-5.0d0)**2.0d0 !Binh-Korn 2
    elseif (funcname.eq.'kursawe') then
       do i=1,3
          f2 = f2+abs(x(i))**(0.8d0)+5.0d0*sin(x(i)**3.0d0) !Kursawe 2
       enddo
    endif
  end function f2

end program main
