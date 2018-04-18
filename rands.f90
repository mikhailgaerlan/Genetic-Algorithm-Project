module rands
  implicit none

  private
  public :: set_seed,random_integer,random_gauss

contains

  !==========================
  !    Random Integer
  !==========================
  subroutine random_integer(min,max,randint)
    integer,intent(in) :: min, max
    integer,intent(out) :: randint
    real(kind=8) :: rand
    call random_number(rand)
    
    randint = min+floor((max-min+1)*rand)
  end subroutine random_integer
  
  !==========================
  !    Random Integer
  !==========================
  subroutine random_gauss(rand,mu,sigma)
    real(kind=8),intent(in) :: mu,sigma
    real(kind=8),intent(out) :: rand
    real(kind=8),parameter :: twopi = 8.0d0*atan(1.0d0)
    real(kind=8) :: t,phi,x1,x2,x3
    
    call random_number(t)
    call random_number(phi)
    x1 = sqrt(-2.0d0*log(t))*cos(twopi*phi)
    x1 = sqrt(-2.0d0*log(t))*sin(twopi*phi)
    
    call random_number(x3)
    if (x3.lt.0.5d0) then
       rand = x1*sigma+mu
    else
       rand = x2*sigma+mu
    endif
  end subroutine random_gauss

  !==========================
  !       Set Seed
  !==========================
  subroutine set_seed
    integer,allocatable :: seed(:)
    integer :: n_seed,clock,i
    real(kind=8) :: x

    call random_seed(size=n_seed)
    allocate(seed(n_seed))
    call system_clock(count=clock)
    do i=1,n_seed
       seed(i)=clock/i
    enddo
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine set_seed

end module rands
