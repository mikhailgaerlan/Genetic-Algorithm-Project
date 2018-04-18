program tester
  use rands
  use convert
  use gen_tools
  implicit none
  
  real(kind=8), dimension(2) :: pop,off
  real(kind=8) :: rand
  integer :: i
  call set_seed

  do i=1,2
     call random_number(rand)
     pop(i) = (10.d0)*(rand-0.5d0)
  enddo
  write(*,*) pop
  call unifcross(pop(1),pop(2),pop(1),pop(2),0.5d0)
  write(*,*) pop
  
end program tester
