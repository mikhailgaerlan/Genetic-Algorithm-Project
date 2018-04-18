module convert
  implicit none

  private
  public :: dec2bin, bin2dec

contains

  !==========================
  !    Decimal to Binary
  !==========================
  subroutine dec2bin(num,binstr)
    real(kind=8),intent(in) :: num
    character(len=64),intent(out) :: binstr
    character(len=64) :: sign, exponent, mantissa

    write(binstr, '(b64.64)') num
  end subroutine dec2bin

  !==========================
  !    Binary to Decimal
  !==========================
  subroutine bin2dec(binstr,num)
    character(len=64),intent(in) :: binstr
    real(kind=8),intent(out) :: num
    character(len=64) :: signstr,exponentstr,mantissastr
    real(kind=8) :: mantissa
    integer :: temp,exponent,sign
    integer :: i
    
    signstr = binstr(1:1)
    exponentstr = binstr(2:12)
    mantissastr = binstr(13:)
    
    read(signstr, *) sign
    sign = (-1)**sign

    exponent = 0
    do i=1,11
       read(exponentstr(i:i), *) temp
       exponent = exponent+temp*(2**(11-i))
    enddo
    exponent = exponent-1023
    
    mantissa = 0.0
    do i=1,52
       read(mantissastr(i:i), *) temp
       mantissa = mantissa+temp*(2.0d0**(-i))
    enddo
    mantissa = 1.0d0+mantissa

    num = sign*(2.0d0**exponent)*mantissa
  end subroutine bin2dec

end module convert
