module sinc_mod
  implicit none
  private
  public :: sincinterp1d!,sincinterp2d
  real,parameter   :: pi=3.14159265359

  contains

  ! Whitaker-Shannon Interpolator ('sinc interpolation')
  ! From regular to irregular
  subroutine sincinterp1d(n1,f,g,s,n)
    real :: f(:),g(:),s(:)
    integer :: ii,kk,n1,n
    g=0.
    do ii=1,n1
      do kk=-n,n
        if ((ii+kk).lt. 1) cycle
        if ((ii+kk).gt.n1) cycle
        g(ii)=g(ii)+sinc(s(ii)+real(kk))*f(ii+kk)
      end do
    end do
  end subroutine
  !subroutine sincinterp2d(n1,n2,f,g,s,n)
  !  real :: f(:),g(:),s(:)
  !  integer :: ii,jj,kk,n1,n2,n
  !  g=0.
  !  do jj=1,n2
  !  do ii=1,n1
  !    do kk=-n,n
  !      if ((ii+kk).lt. 1) cycle
  !      if ((ii+kk).gt.n1) cycle
  !      g((jj-1)*n1+ii)=g((jj-1)*n1+ii)+sinc(s((jj-1)*n1+ii)+kk)*f((jj-1)*n1+ii+kk)
  !    end do
  !  end do
  !  end do
  !end subroutine

  ! Simple sinc function
  real function sinc(x)
    real :: x
    sinc=1.
    if (x.ne.0.) sinc=sin(pi*x)/(pi*x)
  end function

end module
