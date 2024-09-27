module gradients2d_mod ! Lou edit
  use lsqrDataModule, only: dp
  use processing_mod

  implicit none

  private
  public :: gradients2d_mod_init,gradients2d_mod_close
  public :: gradient1_2d_scaled_for,gradient1_2d_scaled_adj
  public :: gradient2_2d_scaled_for,gradient2_2d_scaled_adj
  public :: laplacian_2d_scaled_for,laplacian_2d_scaled_adj
  public :: gradient1_2d_for,gradient2_2d_for
  public :: gradient1_2d_for_dd,gradient1_2d_adj_dd
  public :: gradient1_2d_scaled_for_x, gradient1_2d_scaled_adj_x

  integer :: n1,n2
  real    :: d1

  contains

 !--------------------------------------------------------------
  subroutine gradients2d_mod_init(n1_in,d1_in,n2_in)
    integer :: n1_in,n2_in
    real    :: d1_in
    call gradients2d_mod_close()
    n1=n1_in
    n2=n2_in
    d1=d1_in
  end subroutine 

  subroutine gradients2d_mod_close()
  end subroutine

 !--------------------------------------------------------------
  subroutine gradient1_2d_for_dd(m,d)
    real(dp):: m(n1,n2),d(n1,n2)
    integer :: ii,jj
    do jj=1,n2
    do ii=1,n1-1
      d(ii,jj)=d(ii,jj)+m(ii+1,jj)/d1
      d(ii,jj)=d(ii,jj)-m(ii  ,jj)/d1
    end do
    d(n1,jj)=d(n1,jj)+m(n1  ,jj)/d1
    d(n1,jj)=d(n1,jj)-m(n1-1,jj)/d1
    end do
  end subroutine

 !--------------------------------------------------------------
  subroutine gradient1_2d_adj_dd(m,d)
    real(dp):: m(n1,n2),d(n1,n2)
    integer :: ii,jj
    do jj=1,n2
    do ii=1,n1-1
      m(ii+1,jj)=m(ii+1,jj)+d(ii,jj)/d1
      m(ii  ,jj)=m(ii  ,jj)-d(ii,jj)/d1
    end do
    m(n1  ,jj)=m(n1  ,jj)+d(n1,jj)/d1
    m(n1-1,jj)=m(n1-1,jj)-d(n1,jj)/d1
    end do
  end subroutine

 !--------------------------------------------------------------
  subroutine gradient1_2d_for(m,d)
    real    :: m(n1,n2),d(n1,n2)
    integer :: ii,jj
    do jj=1,n2
    do ii=1,n1-1
      d(ii,jj)=d(ii,jj)+m(ii+1,jj)/d1
      d(ii,jj)=d(ii,jj)-m(ii  ,jj)/d1
    end do
    d(n1,jj)=d(n1,jj)+m(n1  ,jj)/d1
    d(n1,jj)=d(n1,jj)-m(n1-1,jj)/d1
    end do
  end subroutine

 !--------------------------------------------------------------
  subroutine gradient2_2d_for(m,d)
    real    :: m(n1,n2),d(n1,n2)
    integer :: ii,jj
    do jj=1,n2
    d(1,jj)=d(1,jj)+  m(1,jj)/d1/d1
    d(1,jj)=d(1,jj)-2*m(2,jj)/d1/d1
    d(1,jj)=d(1,jj)+  m(3,jj)/d1/d1
    do ii=2,n1-1
      d(ii,jj)=d(ii,jj)+  m(ii-1,jj)/d1/d1
      d(ii,jj)=d(ii,jj)-2*m(ii  ,jj)/d1/d1
      d(ii,jj)=d(ii,jj)+  m(ii+1,jj)/d1/d1
    end do
    d(n1,jj)=d(n1,jj)+  m(n1-2,jj)/d1/d1
    d(n1,jj)=d(n1,jj)-2*m(n1-1,jj)/d1/d1
    d(n1,jj)=d(n1,jj)+  m(n1  ,jj)/d1/d1
    end do
  end subroutine

 !--------------------------------------------------------------
  subroutine gradient2_2d_scaled_for(m,d,eps)
    real(dp)    :: m(n1,n2),d(n1,n2)
    real        :: eps
    integer :: ii,jj
    do jj=2,n2
    d(1,jj)=d(1,jj)+  m(1,jj)*eps
    d(1,jj)=d(1,jj)-2*m(2,jj)*eps
    d(1,jj)=d(1,jj)+  m(3,jj)*eps
    do ii=2,n1-1
      d(ii,jj)=d(ii,jj)+  m(ii-1,jj)*eps
      d(ii,jj)=d(ii,jj)-2*m(ii  ,jj)*eps
      d(ii,jj)=d(ii,jj)+  m(ii+1,jj)*eps
    end do
    d(n1,jj)=d(n1,jj)+  m(n1-2,jj)*eps
    d(n1,jj)=d(n1,jj)-2*m(n1-1,jj)*eps
    d(n1,jj)=d(n1,jj)+  m(n1  ,jj)*eps
    end do
  end subroutine

 !--------------------------------------------------------------
  subroutine gradient2_2d_scaled_adj(m,d,eps)
    real(dp)    :: m(n1,n2),d(n1,n2)
    real        :: eps
    integer :: ii,jj

    do jj=2,n2
    m(n1  ,jj)=m(n1  ,jj)+  d(n1,jj)*eps
    m(n1-1,jj)=m(n1-1,jj)-2*d(n1,jj)*eps
    m(n1-2,jj)=m(n1-2,jj)+  d(n1,jj)*eps

    do ii=2,n1-1
      m(ii+1,jj)=m(ii+1,jj)+  d(ii,jj)*eps
      m(ii  ,jj)=m(ii  ,jj)-2*d(ii,jj)*eps
      m(ii-1,jj)=m(ii-1,jj)+  d(ii,jj)*eps
    end do

    m(3,jj)=m(3,jj)+  d(1,jj)*eps
    m(2,jj)=m(2,jj)-2*d(1,jj)*eps
    m(1,jj)=m(1,jj)+  d(1,jj)*eps
    end do
  end subroutine

 !--------------------------------------------------------------
  subroutine gradient1_2d_scaled_for(m,d,eps)
    real(dp)    :: m(n1,n2),d(n1,n2)
    real        :: eps
    integer     :: ii,jj
    do jj=1,n2
    do ii=1,n1-1
      d(ii,jj)=d(ii,jj)+m(ii+1,jj)*eps
      d(ii,jj)=d(ii,jj)-m(ii  ,jj)*eps
    end do
    d(n1,jj)=d(n1,jj)+m(n1  ,jj)*eps
    d(n1,jj)=d(n1,jj)-m(n1-1,jj)*eps
    end do
  end subroutine

!---------------------------------------------------------------
  subroutine gradient1_2d_scaled_adj(m,d,eps)
    real(dp)  :: m(n1,n2),d(n1,n2)
    real      :: eps
    integer   :: ii,jj
    
    do jj=1,n2
    do ii=1,n1-1
      m(ii  ,jj)=m(ii  ,jj)-d(ii,jj)*eps
      m(ii+1,jj)=m(ii+1,jj)+d(ii,jj)*eps
    end do
    m(n1-1,jj) = m(n1-1,jj) - d(n1,jj)*eps
    m(n1  ,jj) = m(n1  ,jj) + d(n1,jj)*eps
    end do
  end subroutine
 !--------------------------------------------------------------

  subroutine laplacian_2d_scaled_for(m,d,eps)
    real(dp)    :: m(n1,n2),d(n1,n2)
    real        :: eps
    integer :: ii,jj
    d(1,1)=d(1,1)+m(1,1)*eps
    do ii=2,n1-1
      d(ii,1 )=d(ii,1 )+  m(ii-1,1  )*eps
      d(ii,1 )=d(ii,1 )+  m(ii+1,1  )*eps
      d(ii,1 )=d(ii,1 )+  m(ii  ,1+1)*eps
      d(ii,1 )=d(ii,1 )-3*m(ii  ,1  )*eps
    end do
    d(n1,1)=d(n1,1)+m(n1,1)*eps
    do jj=2,n2-1
      d(1 ,jj)=d(1 ,jj)+  m(1  ,jj-1)*eps
      d(1 ,jj)=d(1 ,jj)+  m(1+1,jj  )*eps
      d(1 ,jj)=d(1 ,jj)+  m(1  ,jj+1)*eps
      d(1 ,jj)=d(1 ,jj)-3*m(1  ,jj  )*eps
      do ii=2,n1-1
        d(ii,jj)=d(ii,jj)+  m(ii,jj-1)*eps
        d(ii,jj)=d(ii,jj)+  m(ii-1,jj)*eps
        d(ii,jj)=d(ii,jj)-4*m(ii  ,jj)*eps
        d(ii,jj)=d(ii,jj)+  m(ii+1,jj)*eps
        d(ii,jj)=d(ii,jj)+  m(ii,jj+1)*eps
      end do
      d(n1,jj)=d(n1,jj)-3*m(n1  ,jj  )*eps
      d(n1,jj)=d(n1,jj)+  m(n1  ,jj-1)*eps
      d(n1,jj)=d(n1,jj)+  m(n1-1,jj  )*eps
      d(n1,jj)=d(n1,jj)+  m(n1  ,jj+1)*eps
    end do
    d(1,n2)=d(1,n2)+m(1,n2)*eps
    do ii=2,n1-1
      d(ii,n2)=d(ii,n2)-3*m(ii  ,n2)*eps
      d(ii,n2)=d(ii,n2)+  m(ii,n2-1)*eps
      d(ii,n2)=d(ii,n2)+  m(ii-1,n2)*eps
      d(ii,n2)=d(ii,n2)+  m(ii+1,n2)*eps
    end do
    d(n1,n2)=d(n1,n2)+m(n1,n2)*eps
  end subroutine
!---------------------------------------------------------------

  subroutine laplacian_2d_scaled_adj(m,d,eps)
    real(dp)    :: m(n1,n2),d(n1,n2)
    real        :: eps
    integer :: ii,jj

    m(1,1)=m(1,1)+d(1,1)*eps
    do ii=2,n1-1
      m(ii-1,1  )=m(ii-1,1  ) +  d(ii,1 )*eps
      m(ii+1,1  )=m(ii+1,1  ) +  d(ii,1 )*eps
      m(ii  ,1+1)=m(ii  ,1+1) +  d(ii,1 )*eps
      m(ii  ,1  )=m(ii  ,1  ) -3*d(ii,1 )*eps
    end do
    m(n1,1)=m(n1,1)+d(n1,1)*eps
    do jj=2,n2-1
      m(1  ,jj-1) = m(1  ,jj-1) +  d(1 ,jj)*eps
      m(1+1,jj  ) = m(1+1,jj  ) +  d(1 ,jj)*eps
      m(1  ,jj+1) = m(1  ,jj+1) +  d(1 ,jj)*eps
      m(1  ,jj  ) = m(1  ,jj  ) -3*d(1 ,jj)*eps
      do ii=2,n1-1
        m(ii,jj-1)=m(ii,jj-1) +  d(ii,jj)*eps
        m(ii-1,jj)=m(ii-1,jj) +  d(ii,jj)*eps
        m(ii  ,jj)=m(ii  ,jj) -4*d(ii,jj)*eps
        m(ii+1,jj)=m(ii+1,jj) +  d(ii,jj)*eps
        m(ii,jj+1)=m(ii,jj+1) +  d(ii,jj)*eps
      end do
      m(n1  ,jj  )=m(n1  ,jj  ) -3*d(n1,jj)*eps
      m(n1  ,jj-1)=m(n1  ,jj-1) +  d(n1,jj)*eps
      m(n1-1,jj  )=m(n1-1,jj  ) +  d(n1,jj)*eps
      m(n1  ,jj+1)=m(n1  ,jj+1) +  d(n1,jj)*eps
    end do
    m(1,n2)=m(1,n2)+d(1,n2)*eps
    do ii=2,n1-1
      m(ii  ,n2)=m(ii  ,n2) - 3*d(ii,n2)*eps
      m(ii,n2-1)=m(ii,n2-1) +   d(ii,n2)*eps
      m(ii-1,n2)=m(ii-1,n2) +   d(ii,n2)*eps
      m(ii+1,n2)=m(ii+1,n2) +   d(ii,n2)*eps
    end do
    m(n1,n2)=m(n1,n2)+d(n1,n2)*eps

  end subroutine



!----------------------------------------------07/05/2024
  subroutine gradient1_2d_scaled_for_x(m,d,eps)
    real(dp)    :: m(n1,n2),d(n1,n2)
    real        :: eps
    integer     :: ii,jj

    do ii=1,n1
    do jj=1,n2-1
      d(ii,jj)=d(ii,jj)+m(ii,jj+1)*eps
      d(ii,jj)=d(ii,jj)-m(ii  ,jj)*eps
    end do
    d(ii,n2)=d(ii,n2)+m(ii  ,n2)*eps
    d(ii,n2)=d(ii,n2)-m(ii,n2-1)*eps
    end do
  end subroutine

!---------------------------------------------------------------
  subroutine gradient1_2d_scaled_adj_x(m,d,eps)
    real(dp)  :: m(n1,n2),d(n1,n2)
    real      :: eps
    integer   :: ii,jj
    

    do ii=1,n1
    do jj=1,n2-1
      m(ii  ,jj)=m(ii  ,jj)-d(ii,jj)*eps
      m(ii,jj+1)=m(ii,jj+1)+d(ii,jj)*eps
    end do
    m(ii, n2-1) = m(ii, n2-1) - d(ii, n2)*eps
    m(ii, n2) = m(ii, n2) + d(ii, n2)*eps
    end do
  end subroutine
!---------------------------------------------------------------

end module
