module warp2d_mod ! Lou edit for 1D & inversion - passes dot test
  use sinc_mod
  use gradients2d_mod
  use lsqrModule, only: lsqr
  use lsqrDataModule, only: dp
  implicit none

  private
  public   :: warp_mod_init,warp_dottest,warp_mod_close,warp_update_model,warp_res
  public   :: warp_inv_for,warp_inv_adj,warp_inv

  integer  :: n1,n2
  real     :: eps_tik
  real(dp) :: eps_mod

  real(dp) :: atol,btol,conlim
  integer  :: niter,itnlim
  logical  :: wantse,strainreg

  real,   allocatable,dimension(:) :: grads
  real,   allocatable,dimension(:) :: modelA,modelB,modelW,modelBASE
  real,   allocatable,dimension(:) :: shifts,res
  real(dp),allocatable             :: pshifts(:),b(:),t(:)


  contains

  !!!!!!!! warp module initialisation !!!!!!!!
  subroutine warp_mod_init(n1_in,n2_in,modelA_in,modelB_in,eps_tik_in,eps_mod_in,niter_in,itnlim_in,strainreg_in)
    integer :: n1_in,n2_in,niter_in,itnlim_in,n
    real :: eps_tik_in,eps_mod_in
    real :: modelA_in(:),modelB_in(:)
    logical :: strainreg_in
    call warp_mod_close()
    !!! external parameters
    n1=n1_in
    n2=n2_in
    n=n1*n2
    strainreg=strainreg_in
    eps_tik=eps_tik_in; eps_mod=eps_mod_in
    niter=niter_in; itnlim=itnlim_in
    call gradients2d_mod_init(n1,1.,n2)
    !!! warp parameters
    allocate(modelA(n)); modelA=modelA_in
    allocate(modelBASE(n),modelB(n)); modelBASE=modelB_in; modelB=0.
    allocate(modelW(n)); modelW=0.
    allocate(shifts(n)); shifts=0.
    allocate( grads(n)); grads=0.
    !!! inversion parameters
    allocate(pshifts(n)); pshifts=0.
    allocate(res(n)); res=0.
    allocate(t(n)); t=0.
    allocate(b(2*n)); b=0.
    wantse=.false.
    atol=1.0e-18; btol=atol
    conlim=1.0e+8
  end subroutine

  !!!!!!!! warp module closing !!!!!!!!
  subroutine warp_mod_close
    call gradients2d_mod_close()
    if (allocated(modelA))    deallocate(modelA)
    if (allocated(modelB))    deallocate(modelB)
    if (allocated(modelBASE)) deallocate(modelBASE)
    if (allocated(modelW))    deallocate(modelW)
    if (allocated(shifts))    deallocate(shifts)
    if (allocated(grads))     deallocate(grads)
    if (allocated(pshifts))   deallocate(pshifts)
    if (allocated(b))         deallocate(b)
    if (allocated(t))         deallocate(t)
    if (allocated(res))       deallocate(res)
  end subroutine

  !!!!!!!! warp dottest (for adjointedness) !!!!!!!!
  subroutine warp_dottest()
    real(dp),allocatable           :: x1(:),x2(:),y1(:),y2(:)
    real                           :: sum1,sum2,rand,drand
    integer                        :: ii,n
    n=n1*n2
    !!! Fill the shifts with random numbers !!!
    do ii=1,n; call random_number(rand); shifts(ii)=rand; end do
    call warp_update_model(shifts)
    !!! Fill x1 and y2 with random numbers !!!
    allocate(x1(n),x2(n),y1(n*2),y2(n*2)); x1=0.0_dp; x2=0.0_dp; y1=0.0_dp;y2=0.0_dp
    do ii=1,n;   call random_number(drand); x1(ii)=drand; end do
    do ii=1,n*2; call random_number(drand); y2(ii)=drand; end do
    !!! Call forward and adjoint !!!
    call warp_inv_for(n*2,n,x1,y1)
    call warp_inv_adj(n*2,n,x2,y2)
    !!! Evaluate inner products !!!
    sum1=0.0_dp
    sum2=0.0_dp
    do ii=1,n;   sum1=sum1+x1(ii)*x2(ii); end do
    do ii=1,n*2; sum2=sum2+y1(ii)*y2(ii); end do
    write(0,*) "dot test results: ",real(sum1),real(sum2)
    !!! Reset shifts !!!
    shifts=0.
    deallocate(x1,x2,y1,y2)
  end subroutine

  !!!!!!!! warp inversion !!!!!!!!
  subroutine warp_inv(modelA_in,shifts_out,modelW_out)
    real,allocatable,dimension(:)     :: modelA_in,shifts_out,modelW_out
    real(dp),allocatable              :: se(:)
    real(dp)                          :: Anorm,Acond,rnorm,Arnorm,xnorm ! outputs
    integer                           :: istop,itn ! outputs
    integer                           :: iter,n
    modelA=modelA_in
    modelB=modelBASE
    !write(0,*) 'maxval(shifts_out),mean(shifts_out)',maxval(shifts_out),sum(shifts_out)/size(shifts_out)
    !write(0,*) 'maxval(modelA),mean(modelA)', maxval(modelA),sum(modelA)/size(modelA)
    !write(0,*) 'maxval(modelB),mean(modelB)', maxval(modelB),sum(modelB)/size(modelB)
    shifts=shifts_out; grads=0.; res=0.
    n=n1*n2
    !call warp_dottest()
    allocate(se(n))
    do iter=1,niter
      !!! calculate residual, modelW and grads
      call warp_res(shifts,res)
      !write(0,*) iter,sum(abs(shifts)),sum(abs(res))
      b(1:n)=dble(res);b(n+1:n*2)=0.
      pshifts=0.0_dp
      !!! compute shift perturbation
      !write(0,*) iter,sum(abs(modelW))
      call lsqr(2*n,n,warp_inv_for,warp_inv_adj,b,eps_mod,wantse,pshifts,se,atol,btol,conlim, itnlim,-1,istop,itn,Anorm,Acond,rnorm,Arnorm,xnorm)
      !!! update shifts
      shifts=shifts+real(pshifts)
    end do
    deallocate(se)
    shifts_out=shifts
    modelW_out=modelW
  end subroutine

  !!!!!!!! warp residual calculation (inc. model update) !!!!!!!!
  subroutine warp_res(d,r)
    real :: d(:),r(:)
    integer :: n
    call warp_update_model(d)
    n=n1*n2
    r(1:n)=modelA-modelW
  end subroutine

  !!!!!!!! warp model update (inc. interpolation and gradient) !!!!!!!!
  subroutine warp_update_model(shifts)
    real :: shifts(:)
    integer :: n
    n=n1*n2
    !write(0,*) sum(abs(modelW)),sum(abs(grads))
    call sincinterp1d(n,modelB,modelW,shifts,20)
    !write(0,*) sum(abs(modelW)),sum(abs(grads))
    call gradient1_2d_for(modelW,grads)
    !write(0,*) sum(abs(modelW)),sum(abs(grads))
  end subroutine

  !!!!!!!! warp inversion forward operator !!!!!!!!
  subroutine warp_inv_for(nd,nm,m,d) ! leave format - probably needed for lsqr???
    integer,intent(in)     :: nd,nm
    real(dp),intent(in)    :: m(nm)
    real(dp),intent(inout) :: d(nd)
    integer :: n
    ! Evaluate the block column forward operator
    n=n1*n2
    call warp_grad_wop_for(m,d(1:n))
    !call gradient2_2d_scaled_for(m,d(n+1:n*2),eps_tik)
    t=0.
    !call gradient1_2d_for_dd(m,t)
    call laplacian_2d_scaled_for(m,d(n+1:n*2),eps_tik)
  end subroutine

  !!!!!!!! warp inversion adjoint operator !!!!!!!!
  subroutine warp_inv_adj(nd,nm,m,d) ! leave format - probably needed for lsqr???
    integer,intent(in)     :: nd,nm
    real(dp),intent(in)    :: d(nd)
    real(dp),intent(inout) :: m(nm)
    integer :: n
    n=n1*n2
    ! Evaluate the block column adjoint operator (a block row)
    t=0.
    call laplacian_2d_scaled_adj(m,d(n+1:n*2),eps_tik)
    !call gradient1_2d_adj_dd(m,t)
    !call gradient2_2d_scaled_adj(m,d(n+1:n*2),eps_tik)
    call warp_grad_wop_adj(m,d(1:n))
  end subroutine

  !!!!!!!! warp forward gradient calculation !!!!!!!!
  subroutine warp_grad_wop_for(m,d)
    real(dp) :: m(:),d(:)
    d=d-grads*m
  end subroutine

  !!!!!!!! warp adjoint gradient calculation !!!!!!!!
  subroutine warp_grad_wop_adj(m,d)
    real(dp) :: m(:),d(:)
    m=m-grads*d
  end subroutine

end module
