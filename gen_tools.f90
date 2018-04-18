module gen_tools
  use rands
  use convert
  implicit none

  private
  public :: oneptcross,twoptcross,unifcross,mutflipbit,mutgauss,rank,sort,dominsort,distancing

contains

  !==========================
  !     One Point Cross
  !==========================
  subroutine oneptcross(parent1,parent2,offspring1,offspring2)
    real(kind=8),intent(in) :: parent1,parent2
    real(kind=8),intent(out) :: offspring1,offspring2
    character(len=64) :: parent1bin,parent2bin,offspring1bin,offspring2bin
    integer :: randint
    
    call dec2bin(parent1,parent1bin)
    call dec2bin(parent2,parent2bin)
    call random_integer(13,62,randint)
    
    offspring1bin(:randint)=parent1bin(:randint)
    offspring1bin(randint+1:)=parent2bin(randint+1:)
    offspring2bin(:randint)=parent2bin(:randint)
    offspring2bin(randint+1:)=parent1bin(randint+1:)
    
    call bin2dec(offspring1bin,offspring1)
    call bin2dec(offspring2bin,offspring2)
  end subroutine oneptcross

  !==========================
  !     Two Point Cross
  !==========================
  subroutine twoptcross(parent1,parent2,offspring1,offspring2)
    real(kind=8),intent(in) :: parent1,parent2
    real(kind=8),intent(out) :: offspring1,offspring2
    character(len=64) :: parent1bin,parent2bin,offspring1bin,offspring2bin
    integer :: randint1,randint2
        
    call dec2bin(parent1,parent1bin)
    call dec2bin(parent2,parent2bin)
    call random_integer(13,59,randint1)
    call random_integer(randint1+2,62,randint2)
    
    offspring1bin(:randint1)=parent1bin(:randint1)
    offspring1bin(randint1+1:randint2)=parent2bin(randint1+1:randint2)
    offspring1bin(randint2+1:)=parent1bin(randint2+1:)
    offspring2bin(:randint1)=parent2bin(:randint1)
    offspring2bin(randint1+1:randint2)=parent1bin(randint1+1:randint2)
    offspring2bin(randint2+1:)=parent2bin(randint2+1:)
    
    call bin2dec(offspring1bin,offspring1)
    call bin2dec(offspring2bin,offspring2)
  end subroutine twoptcross
  
  !==========================
  !   Uniform Crossover
  !==========================
  subroutine unifcross(parent1,parent2,offspring1,offspring2,prob)
    real(kind=8),intent(in) :: parent1,parent2,prob
    real(kind=8),intent(out) :: offspring1,offspring2
    character(len=64) :: parent1bin,parent2bin,offspring1bin,offspring2bin
    real(kind=8) :: rand
    integer :: i
    
    call dec2bin(parent1,parent1bin)
    call dec2bin(parent2,parent2bin)
    
    do i=1,12
       offspring1bin(i:i) = parent1bin(i:i)
       offspring2bin(i:i) = parent2bin(i:i)
    enddo
    do i=13,64
       call random_number(rand)
       if (rand.lt.prob) then
          offspring1bin(i:i) = parent2bin(i:i)
          offspring2bin(i:i) = parent1bin(i:i)
       else
          offspring1bin(i:i) = parent1bin(i:i)
          offspring2bin(i:i) = parent2bin(i:i)
       endif
    enddo
    
    call bin2dec(offspring1bin,offspring1)
    call bin2dec(offspring2bin,offspring2)
  end subroutine unifcross

  !==========================
  !    Mutation Flip Bit
  !==========================
  subroutine mutflipbit(individual,prob)
    real(kind=8),intent(inout) :: individual
    real(kind=8),intent(in) :: prob
    real(kind=8) :: rand
    character(len=64) :: indbin
    integer :: i
    
    call dec2bin(individual,indbin)
    do i=1,64
       call random_number(rand)
       if (rand < prob) then
          if (indbin(i:i)=="1") then
             indbin(i:i) = "0"
          else
             indbin(i:i) = "1"
          endif
       endif
    enddo
    call bin2dec(indbin,individual)    
  end subroutine mutflipbit

  !==========================
  !    Mutation Gaussian
  !==========================
  subroutine mutgauss(individual,mu,sigma)
    real(kind=8),intent(inout) :: individual
    real(kind=8),intent(in) :: mu,sigma
    real(kind=8) :: rand
    
    call random_gauss(rand,mu,sigma)
    individual = individual+rand
  end subroutine mutgauss
  
  !==========================
  !     Rank Assignment
  !==========================
  subroutine rank(fitness)
    real(kind=8),intent(inout) :: fitness(:,:)
    real(kind=8),allocatable :: fitnesses(:)
    real(kind=8) :: tier
    logical :: dominating
    integer,dimension(2) :: fitshape
    integer :: i,j,k,ranked
    fitshape = shape(fitness)
    allocate(fitnesses(fitshape(1)))
    fitnesses = 0.0d0
    do i=1,fitshape(1)
       fitness(i,fitshape(2)-1) = 0.0d0
    enddo
    fitnesses = 0.0d0
    tier = 1.0d0
    ranked = fitshape(1)
    do while(ranked.gt.0)
       ranked = 0
       do i=1,fitshape(1)
          if (fitness(i,fitshape(2)-1).eq.0.0d0) then
             do j=1,fitshape(1)
                if ((i.ne.j).and.(fitness(j,fitshape(2)-1).eq.0.0d0)) then
                   dominating = .true.
                   do k=1,fitshape(2)-3
                      if (fitness(i,k).le.fitness(j,k)) then
                         dominating = .false.
                         exit
                      endif
                   enddo
                   if (dominating) then
                      fitnesses(i) = 0.0d0
                      exit
                   else
                      fitnesses(i) = tier
                      ranked = ranked+1
                   endif
                endif
                if (ranked.eq.0) then
                   fitnesses(i) = tier
                   ranked = ranked+1
                endif
             enddo
          endif
       enddo
       do i=1,fitshape(1)
          fitness(i,fitshape(2)-1) = fitnesses(i)
       enddo
       tier = tier+1.0d0
       ranked=0
       do i=1,fitshape(1)
          if (fitness(i,fitshape(2)-1).eq.0) ranked = ranked+1
       enddo
    enddo
  end subroutine rank
  
  !==========================
  !    Tier Bubble Sort
  !==========================
  subroutine sort(pop,fitness)
    real(kind=8),intent(inout) :: pop(:,:)
    real(kind=8),intent(inout) :: fitness(:,:)
    real(kind=8),allocatable :: swapfit(:)
    real(kind=8),allocatable :: swappop(:)
    integer,dimension(2) :: popshape,fitshape
    integer :: i,j
    
    popshape = shape(pop)
    fitshape = shape(fitness)
    allocate(swapfit(fitshape(2)),swappop(popshape(2)))
    do j=1,fitshape(2)-1
       if (j.ne.fitshape(2)-2) then
          i = 1
          do while(i<popshape(1))
             if (fitness(i+1,j)<fitness(i,j)) then
                swapfit(:) = fitness(i,:)
                swappop(:) = pop(i,:)
                fitness(i,:) = fitness(i+1,:)
                pop(i,:) = pop(i+1,:)
                fitness(i+1,:) = swapfit(:)
                pop(i+1,:) = swappop(:)
                i = 1
             else
                i = i+1
             endif
          enddo
       endif
    enddo
  end subroutine sort

  !==========================
  !     Dominating Sort
  !==========================
  subroutine dominsort(pop,fitness)
    real(kind=8),intent(inout) :: pop(:,:)
    real(kind=8),intent(inout) :: fitness(:,:)
    real(kind=8),allocatable :: swappop(:,:)
    real(kind=8),allocatable :: swapfit(:,:)
    integer,dimension(2) :: popshape,fitshape
    real(kind=8) :: tier
    integer :: i,j
    
    popshape = shape(pop)
    fitshape = shape(fitness)
    allocate(swapfit(fitshape(1),fitshape(2)),swappop(popshape(1),popshape(2)))
    swappop=0.0d0
    swapfit=0.0d0
    
    tier=0.0d0
    i = 1
    do while(i.le.popshape(1))
       tier=tier+1.0d0
       do j=1,popshape(1)
          if (fitness(j,fitshape(2)-1).eq.tier) then
             swappop(i,:) = pop(j,:)
             swapfit(i,:) = fitness(j,:)
             pop(j,:) = 0.0d0
             fitness(j,:) = 0.0d0
             i = i+1
          endif
       enddo
    enddo
    pop = swappop
    fitness = swapfit
  end subroutine dominsort
  
  !==========================
  !   Distance Function
  !==========================
  subroutine distancing(fitness)
    real(kind=8),intent(inout) :: fitness(:,:)
    real(kind=8),allocatable :: avgfit(:)
    real(kind=8) :: maxfit
    integer,dimension(2) :: fitshape
    integer :: i
    
    fitshape = shape(fitness)
    allocate(avgfit(fitshape(2)-2))
    fitness(:,fitshape(2)-2) = 0.0d0
    do i=1,fitshape(2)-3
       avgfit(i) = sum(fitness(:,i))/fitshape(1)
    enddo
    do i=1,fitshape(1)
       fitness(i,fitshape(2)-2) = sum((avgfit-fitness(i,:fitshape(2)-3))**2.0d0)
    enddo
    maxfit = maxval(fitness(:,fitshape(2)-2))
    fitness(:,fitshape(2)-2) = fitness(:,fitshape(2)-2)/maxfit
  end subroutine distancing
  
end module gen_tools
