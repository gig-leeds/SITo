module processing_mod
  use rsf

  implicit none
  private

  public :: processing_mod_init, processing_mod_close
  public :: mask_topmute_lop_for !, mask_seabed_lop_for
  public :: write_rsf!, softclip_2, softclip, softclip_1
  public :: taper_rtm!, clip_rtm

  !!!!!!!!  !!!!!!!!
  integer                                :: nx,nz,nt
  integer                                :: nr,ns,ne,nrt
  real                                   :: dz,dx,dt,oz,ox,ot
  real                                   :: vp
  integer                                :: d,n,l
  integer,allocatable,dimension(:,:,:)   :: rcv_zx,src_zx
  real,allocatable,dimension(:)          :: mask
  logical                                :: seam
  integer :: ibatch, nbatch

contains

  !-------------------------------------------------------------
  !!!!!!!! initialise and allocate !!!!!!!!
  subroutine processing_mod_init(nz_in, nx_in, nt_in, dz_in, dx_in, dt_in, &
                                 oz_in, ox_in, ot_in, src_zx_in, rcv_zx_in, &
                                 vp_in, ibatch_in, nbatch_in)
    integer :: ibatch_in, nbatch_in
    real,allocatable,dimension(:,:,:) :: src_zx_in,rcv_zx_in
    real,allocatable,dimension(:)     :: vp_in
    real                              :: dz_in,dx_in,dt_in,oz_in,ox_in,ot_in
    integer                           :: nz_in,nx_in,nt_in
    integer                           :: ie,is,ir,ix,iz
    integer                           :: d0,d1
    call processing_mod_close()
    seam=.false.
    nbatch=nbatch_in
    ibatch=ibatch_in
    
    !!! read in dimensions
    nz=nz_in; nx=nx_in; nt=nt_in
    dz=dz_in; dx=dx_in; dt=dt_in
    oz=oz_in; ox=ox_in; ot=ot_in

    !!! calculate a velocity to use for the direct arrivals
    vp=0.
    do ix=1,nx
      vp=vp+1./sqrt(vp_in((ix-1)*nz+1))
    end do
    vp=vp/nx
    if (seam) vp=1400

    !!! calculate depth to first reflector
    d0=nz;
    do ix=1,nx
      iz=1;d1=nz
      do while ((d1.eq.nz).and.(iz.lt.nz))
        if (1./sqrt(vp_in((ix-1)*nz+iz)).gt.vp) d1=iz ! if vp_in is greater than & 
                                                      !the value of the top layer &
                                                      !and the current value for d is greater than iz
        iz=iz+1
      end do
      if (d0.gt.d1) d0=d1
    end do
    d=d0

    !!! set up geometries 
    nr=size(rcv_zx_in,1)
    ns=size(src_zx_in,1)
    ne=size(rcv_zx_in,3)

    allocate(src_zx(ns,2,ne))
    allocate(rcv_zx(nr,2,ne))
    nrt=nr*nt
    do ie=1,ne
      do is=1,ns
        src_zx(is,1,ie)=1.5+(src_zx_in(is,1,ie)-oz)/dz
        src_zx(is,2,ie)=1.5+(src_zx_in(is,2,ie)-ox)/dx
      end do
      do ir=1,nr
        rcv_zx(ir,1,ie)=1.5+(rcv_zx_in(ir,1,ie)-oz)/dz
        rcv_zx(ir,2,ie)=1.5+(rcv_zx_in(ir,2,ie)-ox)/dx
      end do
    end do

    !!!create mask
    allocate(mask(nr*nt*ne));
    do ie=ibatch,ne,nbatch ! for each experiment
      call mask_topmute_init(src_zx(:,:,ie),rcv_zx(:,:,ie),ie)
    end do
  end subroutine

!---------------------------------------------------------------
  !!!!!!!! initialise top mute mask !!!!!!!!
  subroutine mask_topmute_init(szx,rzx,ie)
    integer,dimension(:,:)            :: szx,rzx
    integer                           :: is,ir,it,ie,tramp
    integer,allocatable,dimension(:)  :: t
    real                              :: l,ds,xsr,t_src
    integer                           :: tind
    allocate(t(nr));t=1
    t_src=0.25
    do is=1,ns
      do ir=1,nr
        tramp=nint(0.2/dt) ! ramp length (should probably be external)

        !! calculate source receiver separation - assumes s
        xsr=abs(szx(is,2)-rzx(ir,2))*dx

        !! calculate reflection distance
        ds=(d-szx(is,1))*dz
        l=sqrt((2*ds)**2+xsr**2)

        !! calculate mute time
        xsr=xsr+t_src*vp
        if (l.lt.xsr) l=xsr
        t(ir)=nint(l/vp/dt)
        if (t(ir).lt.1) t(ir)=1
        
        
        !! make mask
        do it=1,t(ir)
          tind=nrt*(ie-1)+nr*(it-1)+ir
          if (tind.lt.1         ) cycle
          if (tind.gt.size(mask)) cycle
          mask(tind)=0
        end do
        do it=t(ir)+1,t(ir)+tramp
          tind=nrt*(ie-1)+nr*(it-1)+ir
          if (tind.lt.1         ) cycle
          if (tind.gt.size(mask)) cycle
          mask(tind)=real(it-t(ir)-1)/tramp
        end do
        do it=t(ir)+tramp+1,nt
          tind=nrt*(ie-1)+nr*(it-1)+ir
          if (tind.lt.1         ) cycle
          if (tind.gt.size(mask)) cycle
          mask(tind)=1
        end do
      end do
    end do
    deallocate(t)

  end subroutine

  !-------------------------------------------------------------
  !!!!!!!!      apply top mute mask !!!!!!!!
  subroutine mask_topmute_lop_for(m,d)
    real,allocatable,dimension(:) :: m,d
    integer :: nrt,nst,ie,ieb

    !d=m*mask
    d=0.
    nrt=nr*nt
    nst=ns*nt
    do ie=ibatch,ne,nbatch
      ieb=ie/nbatch+1
      if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
      !write(0,*) 'ie',ie,(ieb-1)*nrt+1,ieb*nrt,size(m),size(d)
      d((ieb-1)*nrt+1:ieb*nrt) = m((ieb-1)*nrt+1:ieb*nrt) * &
                               mask((ie-1)*nrt+1:ie*nrt)
    end do
  end subroutine
!------------------------------------------------

!  !!!!!!!! apply seabed mask !!!!!!!!
!  subroutine mask_seabed_lop_for(m,d,nmask,nramp)
!    real,allocatable,dimension(:) :: m,d
!    integer                       :: ix,iz,nmask,nramp
!    d=0.
!    do ix=1,nx
!      do iz=nmask-nramp,nmask-1
!        d((ix-1)*nz+iz)=m((ix-1)*nz+iz)*real(nramp+iz-nmask)/nramp
!      end do
!      do iz=nmask,nz
!        d((ix-1)*nz+iz)=m((ix-1)*nz+iz)
!      end do
!    end do
!  end subroutine

!!---------------------------------------------------------------
!  !!!!!!!!              apply well mask !!!!!!!!
!  subroutine mask_well_lop_for(d, nmask)
!    real,allocatable,dimension(:) :: d
!    integer                       :: iw, nmask, nw
!
!    nw=size(d);
!    if (nmask.ne.0) then
!      do iw=1,nmask
!        d(iw)=0
!        d(nw+1-iw)=0
!      end do
!    end if
!  end subroutine

!!---------------------------------------------------------------
!  subroutine agc_nlop_for(f,g,n,l) ! applies
!    real,allocatable,dimension(:) :: f,g
!    integer                       :: i,id,n,l,li
!    real                          :: mean
!
!    g=0.
!    if (mod(l,2).eq.0) l=l+1
!
!    !!! calculate desired mean !!!
!    mean=sum(abs(f))/n
!    !!! first l/2 data !!!
!    do i=1,ceiling(l/2.)
!      li=i*2-1
!      g(i)=mean*li/sum(abs(f(1:li)))
!      if (sum(abs(f(1:li))).eq.0) g(i)=1
!    end do
!    !!! normal window !!!
!    do i=2,n-l
!      id=i+floor(l/2.)
!      g(id)=mean*l/sum(abs(f(i:i+l-1)))
!      if (sum(abs(f(1:li))).eq.0) g(i)=1
!    end do
!    !!! last l/2 data !!!
!    do i=ceiling(l/2.),1,-1
!      li=i*2-1
!      id=n-i+1
!      g(id)=mean*li/sum(abs(f(n-li+1:n)))
!      if (sum(abs(f(1:li))).eq.0) g(i)=1
!    end do
!    call agc_lop(f,g)
!  end subroutine

!!---------------------------------------------------------------
!  subroutine agc_lop(f,g) ! adjoint and forward based on linearisation
!    real,allocatable,dimension(:) :: f
!    real,dimension(:)             :: g
!    f=f*g
!  end subroutine

!---------------------------------------------------------------
!  subroutine rollmean_for(m,d,n,l)
!    real,allocatable,dimension(:) :: m,d
!    integer                       :: i,id,n,l,li
!    if (mod(l,2).eq.0) l=l+1
!    !!! first l/2 data !!!
!    do i=1,ceiling(l/2.)
!      li=i*2-1
!      !write(0,*) li
!      d(i)=sum(m(1:li))/li
!    end do
!    !!! normal window !!!
!    do i=2,n-l
!      id=i+floor(l/2.)
!      d(id)=sum(m(i:i+l-1))/l
!    end do
!    !!! last l/2 data !!!
!    do i=ceiling(l/2.),1,-1
!      li=i*2-1
!      id=n-i+1
!      d(id)=sum(m(n-li+1:n))/li
!    end do
!  end subroutine
!
!!---------------------------------------------------------------
!  subroutine xcorr(m1,m2,d,n,L) ! assumes equal length
!    real,allocatable,dimension(:)  :: m1,m2,d
!    integer                        :: n,L,i,j,id
!    write(0,*) size(m1),size(m2),size(d)
!    do j=1,n
!      do i=0,L ! loop through output indexes
!        id=L+1+i
!        if ((j+i).le.n) d(id)=d(id)+m1(j)*m2(j+i)
!      end do
!    end do
!    do i=1,L
!      d(L+1-i)=d(L+1+i)
!    end do
!  end subroutine

!!---------------------------------------------------------------
!  subroutine conv(m1,m2,d,n1,n2)
!    real,allocatable,dimension(:)  :: m1,m2,d
!    integer                        :: n1,n2,i,j
!    do i=1,n2 ! loop through output indexes
!      do j=1,n1
!        if (j.lt.i) d(i+j-1)=d(i+j-1)+m1(j)*m2(i)
!      end do
!    end do
!  end subroutine

!---------------------------------------------------------------
  !!!!!!!!           close and deallocate !!!!!!!!
  subroutine processing_mod_close()
    if (allocated(mask)) deallocate(mask)
    if (allocated(src_zx)) deallocate(src_zx)
    if (allocated(rcv_zx)) deallocate(rcv_zx)
  end subroutine

!------------------------------------------------softclip
!  subroutine softclip_1(d, clip_scalar, clipped_d)
!    real,allocatable,dimension(:)  :: d, clipped_d
!    real :: clip_scalar ! g in the softclip eq of Jon Claerbout
!    integer :: n
!    real :: percent_val, gval  
!
!    n = size(d)
!
!    call find_percentile(d, clip_scalar, percent_val)
!    write(0,*) percent_val
!    !percent_val = 4.0 * 1.0e09
!    gval = 1./ percent_val
!
!    !do i = 1, n
!    !  clipped_d(i) = gval * d(i) /sqrt(1. + gval*gval * d(i)*d(i))
!    !end do
!    clipped_d = gval * d /sqrt(1. + gval*gval * d*d)
!
!  end subroutine
!
!!-----------------------------------------------------
!subroutine find_percentile(d, percent, percent_val)
!  real,allocatable,dimension(:)  :: d
!  integer :: i, n, krow
!  real    ::  percent
!  real :: percent_val, buff  
!  !real :: tmp
!
!  d = abs(d)
!  n =  size(d)
!  
!
!  ! sorting d
!  do i = 1, n
!    krow = minloc(d(i:n), dim=1) + i - 1
! 
!    buff = d(i)
!    d(i) = d(krow)
!    d(krow) = buff
!  end do
!  !do i =1, (n * percent)
!  !  tmp = tmp + d(i) 
!  !end do
!  !write(0,*) n * percent, tmp/(n * percent)
!  ! finding percentile
!  percent_val = d(n * percent)
!  write(0,*) 'percent_val = ', percent_val
!  
!end subroutine
!-----------------------------------------------------
!subroutine softclip_tanh(x, percent, y)
!  
!  ! Calculates y(x) = tanh(kx) 
!  ! x: input
!  ! k: scalar value
!  real,allocatable,dimension(:)  :: x, y
!  real :: k, percent
!
!  !call find_percentile(x, percent, k)
!  k = percent
!  !write(0,*) percent
!  y = tanh(k * x)
!
!end subroutine   
!
! !----------------------------------------------------
!subroutine hardclip(x, percent, y)
!
!  ! Calculates y(x) = clip(x, -T, T)
!  real,allocatable,dimension(:)  :: x, y
!  real :: percent, threshold
!  integer :: i
!  call find_percentile(x, percent, threshold)
!  write(0,*) 'percent', percent
!
!  do i = 1, size(x)
!    !y(i) = max(-threshold, min(x(i), threshold))
!    y(i) = max(0., min(x(i), threshold))
!  end do
!
!end subroutine
!
! !----------------------------------------------------
subroutine write_rsf(out, n3, d3, o3, n1, d1, o1, &
                       n2, d2, o2)
  type(file) :: out
  integer    :: n1, n2, n3
  real       :: d1, d2, d3, o1, o2, o3
 
  call to_par(out,'n3',n3)
  call to_par(out,'d3',d3)  
  call to_par(out,'o3',o3) 
 
  call to_par(out,'n1',n1)
  call to_par(out,'d1',d1)  
  call to_par(out,'o1',o1) 


  call to_par(out,'n2',n2)
  call to_par(out,'d2',d2)  
  call to_par(out,'o2',o2) 

end subroutine

!---------------------------------------------------------------
!subroutine laplacian(img, nx, nz, nxyz, img_lap)
!
!  integer :: nx, nz, nxyz, ii, jj
!  real,allocatable,dimension(:) :: img, img_lap
!
!  do jj=2,nx-1
!    do ii=2,nz-1
!      img_lap(nxyz+(jj-1)*nz+ii) = img_lap(nxyz+(jj-1)*nz+ii) + &
!                                   1. * img(nxyz+(jj-0)*nz+ii+0)
!      img_lap(nxyz+(jj-1)*nz+ii) = img_lap(nxyz+(jj-1)*nz+ii) - &
!                                   1. *img(nxyz+(jj-1)*nz+ii  )
!      img_lap(nxyz+(jj-1)*nz+ii) = img_lap(nxyz+(jj-1)*nz+ii) + &
!                                   1. * img(nxyz+(jj-2)*nz+ii+0)
!      img_lap(nxyz+(jj-1)*nz+ii) = img_lap(nxyz+(jj-1)*nz+ii) - &
!                                   1. * img(nxyz+(jj-1)*nz+ii  )
!      img_lap(nxyz+(jj-1)*nz+ii) = img_lap(nxyz+(jj-1)*nz+ii) + &
!                                   1. * img(nxyz+(jj-1)*nz+ii+1)
!      img_lap(nxyz+(jj-1)*nz+ii) = img_lap(nxyz+(jj-1)*nz+ii) - &
!                                   1. * img(nxyz+(jj-1)*nz+ii  )
!      img_lap(nxyz+(jj-1)*nz+ii) = img_lap(nxyz+(jj-1)*nz+ii) + &
!                                   1. * img(nxyz+(jj-1)*nz+ii-1)
!      img_lap(nxyz+(jj-1)*nz+ii) = img_lap(nxyz+(jj-1)*nz+ii) - &
!                                   1. * img(nxyz+(jj-1)*nz+ii  )
!    end do 
!  end do
!end subroutine

!!---------------------------------------------------------------
!  subroutine softclip_2(d, clip_scalar, clipped_d)
!    real,allocatable,dimension(:)  :: d, clipped_d
!    real :: clip_scalar, percentile
! 
!    !call find_percentile(d, clip_scalar, percentile)  not a good choice
!    percentile = 1.0
!    clipped_d = percentile * d / sqrt( 1 + percentile * d &
!                * percentile * d)
!
!  end subroutine 
!!---------------------------------------------------------------
!  subroutine softclip(d, clip, clipped_d)
!    real,allocatable,dimension(:)  :: d, clipped_d
!    real :: clip, k, c
!    integer :: ii
!
!
!    k = 10.
!    do ii = 1, size(d)
!    !c = clip * sign(d(ii)); ! matlab
!    c = sign(clip, d(ii)) !fortran
!    !clipped_d(ii) = -sign(d(ii)) * log(1 + exp(k * sign(d(ii)) * (c - d(ii)))) / k + c  ! matlab
!     clipped_d(ii) = -sign(log(1 + exp(sign(k, d(ii)) * (c - d(ii))) / k), d(ii)) + c  
!                      
!  end do
!
!  end subroutine
!------ubroutine taper_rtm(rtm, n1, n2, rtm_taper)-----------------------------------------------

  subroutine taper_rtm(rtm, nx, nz, rtm_taper)
    
    integer :: nzeros, ntaper, nones, nz, nx, iz, ix
    real, allocatable, dimension(:) :: rtm_taper, taper
    real, allocatable, dimension(:) :: rtm

    real, parameter :: pi = 3.14159265359

    nzeros = 20
    ntaper = 180
    nones = 101

    !allocate(rtm_taper(nx, nz), taper(nx, nz))
    allocate(taper(nx*nz))
    rtm_taper = 0.
    taper = 0.
  !do ix = 1, nx
  !  do iz =1, nz
  !     rtm_taper(ix, iz) = rtm((ix-1) * nz + iz)
  !  end do
  !end do

    do ix = 1, nx
      do iz = 1, nzeros
        !taper(ix, iz) = 0.
        if (iz.lt.nz) taper((ix-1)*nz+iz)= 0.
      end do
      !do iz = nzeros + 1, nzeros + ntaper
             !taper(ix, iz) = 0.5 * (1 - cos((iz - 1)/ ntaper * pi))
      !end do
      !do iz = nzeros + ntaper + 1, nz 
      !   taper(ix, iz) = 1.  
      !end do
      do iz=1, ntaper
           taper((ix-1)* nz + iz + nzeros) = &
                 0.5 * (1.0 - cos((real(iz) -1.0)/ real(ntaper) * pi))
           write(0,*) 'taper', taper((ix-1)* nz + iz + nzeros)
      end do
      do iz=1, nones
        taper((ix-1) * nz + iz + nzeros + ntaper) = 1.
      end do
    
    end do

    !do ix = 1, nx
      !do iz = 1, nz
      !  rtm_taper(ix, iz) = rtm_taper(ix, iz) * taper(ix, iz)
      !end do
    !end do 
    rtm_taper = rtm * taper
  
    !write(0,*) 'taper', rtm(200), taper(200), rtm_taper(200)
    
    deallocate(taper)

  end subroutine
!!------------------------------------------------------
!  subroutine clip_rtm(rtm, clip_val)
!    real, allocatable, intent(inout), dimension(:) :: rtm
!    real, intent(in)  :: clip_val
!    integer :: i
! 
!    do i = 1, size(rtm)
!      !clipped_rtm = rtm
!      if (rtm(i).gt.clip_val) rtm(i) =  clip_val
!      if (rtm(i).lt.-clip_val) rtm(i) =  -clip_val
!    end do
!  
!  end subroutine
!
!------------------------------------------------------
!  subroutine filter_rtm(rtm_in, n1, n2, &
!                        f1, f2, f3, f4, rtm_out)
!    integer :: n1, n2, ic, ix
!    real    :: f1, f2, f3, f4
!    real :: rtm_in(:)
!    complex, allocatable :: arr_fft(:,:)
!    real, allocatable, dimension(:) :: filter, rtm_out
!
!      
!      !real_trace = real_trace/maxval(real_trace) 
!      call fft(n1, n2, rtm_in, arr_fft)
!      call bandpass_filter(n1, n2, &
!!                             0.01, 0.05, 0.1, 0.15, filter)
!                             f1, f2, f3, f4, filter)
!      write(0,*) 'sizes: ', size(arr_fft, 2), size(filter)
!      !arr_fft = arr_fft * cmplx(1.0)
!      do ix =1, n2
!        do ic = 1, size(arr_fft,2)
!          arr_fft(ix, ic) = arr_fft(ix,ic) * filter(ic)
!        end do
!      end do
!      call ifft(n1, n2, arr_fft, rtm_out)
!      write(0,*) 'ifft done!'
!      !call sf_init()
!      !out = rsf_output('fft')
!      !call  write_rsf(out, 0., 0., 0., n1, 10., 0., &
!      !                 n2, 10., 0.)
!
!  end subroutine filter_rtm
!!------------------------------------------------------

!  subroutine fft(n1, n2, arr_in, arr_fft_out)
!    include "fftw3.f"
!    real    :: arr_in(:)
!    real, allocatable :: arr_real(:)
!    complex, allocatable :: arr_fft(:), arr_fft_out(:,:)
!    
!    integer*8           :: planf
!    integer             :: n, nc, ic, ix, iz, n1, n2
! 
!    n = n1
!    nc = n/2 + 1
!
!    if (.not. allocated(arr_fft) ) allocate(arr_fft(nc))
!    allocate(arr_fft_out(nx, nc))
!    allocate(arr_real(n1))
!    
!    arr_real =0. 
!
!
!    do ix = 1, n2
!      do iz = 1, n1
!          arr_real(iz) = arr_in((ix-1) * n1 + iz)
!      end do
!      call sfftw_plan_dft_r2c_1d(planf, n, arr_real, arr_fft, FFTW_ESTIMATE)
!      call sfftw_execute(planf, arr_real, arr_fft)
!      !if (ix <= 100) write(0, *) arr_fft
!      do ic = 1, nc
!        arr_fft_out(ix, ic)  = arr_fft(ic)
!      end do
!      call sfftw_destroy_plan(planf, arr_real, arr_fft) 
!    end do
!
!  end subroutine fft
!
!------------------------------------------------------
!  subroutine ifft(n1, n2, arr_in, arr_ifft_out)
!
!    include "fftw3.f"
!    real,allocatable    :: arr_ifft(:), arr_ifft_out(:)
!    complex, allocatable :: arr_in(:, :), arr_cmplx(:)
!    integer*8           :: plani
!    integer             :: i, n, n1, n2, nc, ix, iz, ic
!
!   
!    n = n1
!    nc = n/2 + 1
!
!    if (.not. allocated(arr_ifft)) allocate(arr_ifft(n1))
!    if (.not. allocated(arr_cmplx)) allocate(arr_cmplx(nc))
!    if (.not. allocated(arr_ifft_out)) allocate(arr_ifft_out(n1*n2))
!
!    do ix = 1, n2
!      arr_cmplx = arr_in(ix,:)
!      call sfftw_plan_dft_c2r_1d(plani, n, arr_cmplx, arr_ifft, FFTW_ESTIMATE)
!      call sfftw_execute(plani, arr_cmplx, arr_ifft)
!
!      do iz = 1, n1
!        arr_ifft_out((ix-1) * n1 + iz) = real(arr_ifft(iz))
!      end do
!      call sfftw_destroy_plan(plani, arr_cmplx, arr_ifft) 
!    end do
!
!
!  end subroutine ifft
!!!---------------------------------------------------------------
!  subroutine bandpass_filter(n1, n2, &
!                             f1, f2, f3, f4, filter)
!
!    integer  :: n1, n2
!    integer  :: f1_ind, f2_ind, f3_ind, f4_ind
!    real :: f1, f2, f3, f4, dkz
!    real, allocatable, dimension(:) :: filter   
!    real, parameter :: pi = 3.14159265359
!
!
!
!    if (.not. allocated(filter)) allocate(filter(n1/2+1))
!    dkz = 1./(n1) ! Depth wavenumber unit
!
!    f1_ind = f1/dkz + 1.5 ; f1_ind = int(f1_ind)
!    f2_ind = f2/dkz + 1.5 ; f2_ind = int(f2_ind)
!    f3_ind = f3/dkz + 1.5 ; f3_ind = int(f3_ind)
!    f4_ind = f4/dkz + 1.5 ; f4_ind = int(f4_ind)
!    filter = 0.
!    write(0,*) 'filter index ', n1, dkz, f1_ind, f2_ind, f3_ind, f4_ind
!    do i = 1, n1
!        if ( i.ge.f1_ind .and. i.le.f2_ind) filter(i) = &
!                                            0.5 * (1.0 -cos((real(i) - 1.0)/&
!                                            (real(f2_ind-f1_ind)-1.0) * pi))
!        if (i.gt.f2_ind .and. i.le.f3_ind) filter(i) = 1.
!        if (i.gt.f3_ind .and. i.le.f4_ind) filter(i) = &
!                                           0.5 * (1.0 + cos((real(i) - 1.0)/ &
!                                           (real(f4_ind - f3_ind) - 1.0) * pi))
!      !write(0,*) 'filter: ', i, filter(i)
!    end do
!    write(0,*) 'filter calculated: ' , 'sum(filter)', sum(filter)
!
!  end subroutine
!!------------------------------------------------------

end module

