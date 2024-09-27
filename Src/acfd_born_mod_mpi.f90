  ! Module For Acoustic Finite Difference Modelling

module acfd_born_mod_mpi 
  use rsf
  !use omp_lib
  use sdr_tmp_io_mpi
  !use mpi_mod
  use processing_mod
  implicit none
  private

  public  :: acfd_born_init,acfd_born_close
  public  :: acfd_born_lop,acfd_born_fmd,acfd_fmd_adjsrc,acfd_born_res,acfd_born_grd,acfd_born_res_grd
  public  :: acfd_born_slo_update,acfd_born_slo2_update,acfd_born_vel_update
  public  :: itoa

  logical :: freesrf

  integer :: nz,nx,nt,nzz,nxx
  real    :: dz,dx,dt
  real    :: oz,ox,ot
  integer :: nr,ns,ne,np,i,neb

  integer,allocatable,dimension(:,:,:) :: rcv_zx,src_zx
  real,allocatable,dimension(:)        :: slo2,vel2,damp
  real,allocatable,dimension(:)        :: p1,p2,p3,lp,s,srcdat,r1,r2,r3,mm
  real                                 :: ddamp
  integer ibatch, nbatch
  ! For OMP
  integer :: thread_num

  ! For SDR scratch io
  character(50) :: path,filename
  integer       :: jslices,nslices

  ! Spatial FD coefficients
  real,parameter   :: c0=-1.0
  real,parameter   :: c1=+1.0

  integer,parameter :: dp=8
  integer(dp) :: io_buff_dptr
  integer     :: io_buff_sptr


contains

  ! Initialize and Allocate
  subroutine acfd_born_init(np_in,nz_in,nx_in,nt_in,dz_in,dx_in,dt_in,oz_in,ox_in,ot_in,srcdat_in,src_zx_in,rcv_zx_in,frsrf,jslices_in,path_in,filename_in,neb_in,ibatch_in,nbatch_in)
    real,allocatable,dimension(:,:,:) ::  src_zx_in,rcv_zx_in
    real,allocatable,dimension(:)     ::  srcdat_in
    real                              ::  dz_in,dx_in,dt_in,oz_in,ox_in,ot_in
    integer                           ::  np_in,nz_in,nx_in,nt_in
    character(50) :: path_in,filename_in
    integer       :: jslices_in
    integer       :: is,ir,ie
    logical       :: frsrf
    integer,optional :: neb_in,ibatch_in,nbatch_in
    i=0 ! allows different filenames for different 'mpi_reduce' locks
    freesrf=frsrf ! Free Surface (or not)
    path = path_in
    ibatch=1
    nbatch=1
    if (present(ibatch_in)) ibatch=ibatch_in
    if (present(nbatch_in)) nbatch=nbatch_in
!    write(0,*) ibatch,'inside 1'

    call acfd_born_close()
    !nw=1; 
    np=np_in
    nz=nz_in
    nx=nx_in
    nt=nt_in
    dz=dz_in
    dx=dx_in
    dt=dt_in
    oz=oz_in
    ox=ox_in
    ot=ot_in
    nzz=nz+2*np
    nxx=nx+2*np
    nr=size(rcv_zx_in,1)
    ns=size(src_zx_in,1)
    ne=size(rcv_zx_in,3)
    allocate(src_zx(ns,2,ne))
    allocate(rcv_zx(nr,2,ne))
    do ie=1,ne
    do is=1,ns
      src_zx(is,1,ie)=1.5+(src_zx_in(is,1,ie)-oz)/dz + np
      src_zx(is,2,ie)=1.5+(src_zx_in(is,2,ie)-ox)/dx + np
    end do
    do ir=1,nr
      rcv_zx(ir,1,ie)=1.5+(rcv_zx_in(ir,1,ie)-oz)/dz + np
      rcv_zx(ir,2,ie)=1.5+(rcv_zx_in(ir,2,ie)-ox)/dx + np
    end do
    end do
    allocate(p1(nzz*nxx),p2(nzz*nxx),p3(nzz*nxx),lp(nzz*nxx),s(nzz*nxx))
    allocate(r1(nr),r2(nr),r3(nr),mm(nzz*nxx))
    ddamp=minval((/real(np)/200,1.0/))

    ! Initialize the wavefield scratchspace
    filename=filename_in
    jslices=jslices_in
    nslices=floor(real(nt-1)/jslices)+1
    !write(0,*) nslices,nt,nt+1,jslices
    !write(0,*) 'path: ', path,ibatch, nbatch
    call tmp_io_init(path, ibatch, nbatch)

    ! Store sourcing wavefield
    allocate(srcdat(ns*nt*ne))
    srcdat=srcdat_in

    ! OMP
    !thread_num = omp_get_max_threads ( )
    !write(0,*), ''
    !write(0,*), 'The number of processors available = ', omp_get_num_procs()
    !write(0,*), 'The number of threads available    = ', thread_num
    !write(0,*), ''
    !write(0,*) 'nzz,nxx',nzz,nxx
    neb=1
    if (present(neb_in)) neb=neb_in

  end subroutine

  subroutine acfd_born_close()
  !  character(50)  :: mlockfile
    !write(0,*)'here 5.1'
    if (allocated(rcv_zx)) deallocate(rcv_zx)
    if (allocated(src_zx)) deallocate(src_zx)
    if (allocated(slo2))   deallocate(slo2)
    if (allocated(vel2))   deallocate(vel2)
    if (allocated(damp))   deallocate(damp)
    if (allocated(p1))     deallocate(p1)
    if (allocated(p2))     deallocate(p2)
    if (allocated(p3))     deallocate(p3)
    if (allocated(lp))     deallocate(lp)
    if (allocated(s))      deallocate(s)
    if (allocated(mm))     deallocate(mm)
    if (allocated(srcdat)) deallocate(srcdat)
    !write(0,*)'here 5.5'
  !  mlockfile='born_end_mlock.dat'
    !call mpi_mlock(mlockfile)
    !call mpi_mrelease(mlockfile)
    call tmp_io_close()
    !write(0,*)'here 5.999'
  end subroutine


  ! Set or Update the Velocity or Slowness
  subroutine  acfd_born_slo_update(slo,updatewfld)
    real,allocatable,dimension(:) :: slo
    integer :: ix,iz
    logical,optional :: updatewfld
    logical :: update
    !call mpi_cast('slo',slo)
!    write(0,*) ibatch,'inside 1a'
    update=.false.
    if (present(updatewfld)) update=updatewfld
    if (.not.allocated(slo2)) allocate(slo2(nzz*nxx)); slo2=0.
    if (.not.allocated(vel2)) allocate(vel2(nzz*nxx)); vel2=0.
    if (.not.allocated(damp)) allocate(damp(nzz*nxx)); damp=0.
    call abl_prep(np,dz,dx,dt,nzz,nxx,minval(1./slo))
!    write(0,*) ibatch,'inside 1b'
    ! Square slowness
    do ix=1,nx
    do iz=1,nz
      slo2((ix-1+np)*nzz+iz+np)=slo((ix-1)*nz+iz)*slo((ix-1)*nz+iz)
    end do
    end do
!    write(0,*) ibatch,'inside 1c'
    call padd_model(slo2)
    vel2=dt*dt/dz/dx/slo2
!    write(0,*) ibatch,'inside 1d'
    if (update) then
      call acfd_born_update_background_wfld()
    end if
!    write(0,*) ibatch,'inside 1e'
  end subroutine

!------------------------------------------------
  subroutine  acfd_born_slo2_update(slo,updatewfld)
    real,allocatable,dimension(:) :: slo
    integer :: ix,iz
    logical,optional :: updatewfld
    logical :: update
    !call mpi_cast('slo2',slo)
    update=.false.
    if (present(updatewfld)) update=updatewfld
    if (.not.allocated(slo2)) allocate(slo2(nxx*nzz)); slo2=0.
    if (.not.allocated(vel2)) allocate(vel2(nxx*nzz)); vel2=0.
    if (.not.allocated(damp)) allocate(damp(nxx*nzz)); damp=0.
    !!!
    !!!
    call abl_prep(np,dz,dx,dt,nzz,nxx,1500.)
    !call abl_prep(np,dz,dx,dt,nzz,nxx,minval(1./sqrt(slo)))
    ! Square slowness
    do ix=1,nx
    do iz=1,nz
      slo2((ix-1+np)*nzz+iz+np)=slo((ix-1)*nz+iz)
    end do
    end do
    call padd_model(slo2)
    vel2=dt*dt/dz/dx/slo2
    if (update) then
      call acfd_born_update_background_wfld()
    end if
    i=0
  end subroutine
!------------------------------------------------

  subroutine  acfd_born_vel_update(vel,updatewfld)
    real,allocatable,dimension(:) :: vel
    integer :: ix,iz
    logical,optional :: updatewfld
    logical :: update
    !call mpi_cast('vel',vel)
    update=.false.
    if (present(updatewfld)) update=updatewfld
    if (.not.allocated(slo2)) allocate(slo2((2*np+nz)*(2*np+nx))); slo2=0.
    if (.not.allocated(vel2)) allocate(vel2(nxx*nzz)); vel2=0.
    if (.not.allocated(damp)) allocate(damp(nxx*nzz)); damp=0.
    call abl_prep(np,dz,dx,dt,nzz,nxx,minval(vel))
    ! Square slowness
    do ix=1,nx
    do iz=1,nz
      slo2((ix-1+np)*nzz+iz+np)=1./(vel((ix-1)*nz+iz)*vel((ix-1)*nz+iz))
    end do
    end do
    call padd_model(slo2)
    vel2=dt*dt/dz/dx/slo2
    if (update) then
      call acfd_born_update_background_wfld()
    end if
  end subroutine

  subroutine padd_model(m)
    real,allocatable,dimension(:) :: m
    integer :: ip,ix,iz
    ! Padd the edges (top & bottom)
    do ip=1,np
    do ix=1,nxx
      m((ix-1)*nzz+      ip)=m((ix-1)*nzz+np+   1)
      m((ix-1)*nzz+np+nz+ip)=m((ix-1)*nzz+np+nz-1)
    end do
    end do
    ! Padd the edges (left & right)
    do ip=1,np
    do iz=1,nzz
      m((      ip-1)*nzz+iz)=m((np+1   )*nzz+  iz)
      m((np+nx+ip-1)*nzz+iz)=m((np+nx-1)*nzz+  iz)
    end do
    end do
  end subroutine

  ! Born Forward Modelling
  subroutine  acfd_born_fmd         (m,d)
    real,allocatable,dimension(:) :: m,d
    ! Model Background of Each Experiment (and store wavefields)
    call acfd_fmd(d,.true.)
    ! Add perturbation (using stored background wavefields)
    call acfd_born_lop(.true.,.false.,m,d)
  end subroutine

  ! Born residual (from background - optionally store wavefields)
  subroutine  acfd_born_res         (d,r,storewfld_in)
    real,allocatable,dimension(:) :: d,r
    logical,optional :: storewfld_in
    logical          :: storewfld
    storewfld=.false.
    if (present(storewfld_in)) storewfld=storewfld_in
    r=-d
    call acfd_fmd(r,storewfld)
  end subroutine

  ! Born gradient
  subroutine  acfd_born_grd         (r,g)
    real,allocatable,dimension(:) :: r,g
    call acfd_born_lop(.true.,.true.,g,r)
  end subroutine

  ! Born residual and gradient
  subroutine  acfd_born_res_grd     (d,r,g)
    real,allocatable,dimension(:) :: d,r,g
    call acfd_born_res(d,r,.true.) ! Store background wavefields for gradient computation
    call acfd_born_grd(r,g)
  end subroutine

!------------------------------------------------
  ! Linear Operator Handle
  subroutine acfd_born_lop(add,adj,m,d,storewfld_in,rtm_in,prestack_in,fd_in)
    logical                       :: add,adj
    real,allocatable,dimension(:) :: m,d
    integer                       :: ie,ieb,nrt,nxz,nxt
    logical,optional              :: storewfld_in,rtm_in,prestack_in
    logical                       :: storewfld,rtm,prestack
    real,optional,dimension(:)    :: fd_in
    real,allocatable,dimension(:) :: fd
    !type(file)                    :: mout
    !character(50)                 :: filename,ie_string

    i=i+1

    storewfld=.false.
    rtm=.false.
    prestack=.false.
    if (present(rtm_in)) rtm=rtm_in 
    if (present(prestack_in)) prestack=prestack_in 
    if (present(storewfld_in)) storewfld=storewfld_in 
    if (present(fd_in)) then
      allocate(fd(nx*nz))
      fd=fd_in
    end if
    call adjadd(add,adj,m,d)
    nrt=nr*nt

    if (.not.adj) then
    ! forward
!      call mpi_mlock('born_lop_for')
      do ie=ibatch,ne,nbatch
        ieb=ie/nbatch+1
        if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
        nxz=(ieb-1)*nx*nz
        nxt=(ieb-1)*nrt
        if (.not.prestack) nxz=0.
        if (.not.prestack) nxt=0.
!        write(0,*) 'acfd_born_lop: FORWARD',ie,nxz+1,nxz+nx*nz,nxt+1,nxt+nrt,size(m),size(d)
        call acfd_born_lop_for(m(nxz+1:nxz+nx*nz),d(nxt+1:nxt+nrt),rcv_zx(:,:,ie),ie)
      end do
!      call mpi_mrelease('born_lop_for')
    else

      ! adjoint
      if (ibatch.ne.nbatch) m=0.
      do ie=ibatch,ne,nbatch
!        write(0,*) 'acfd_born_lop: adjoint',ie
        ieb=ie/nbatch+1
        if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
        nxz=(ieb-1)*nx*nz
        nxt=(ieb-1)*nrt
        if (.not.prestack) nxz=0.
        !if (.not.prestack) nxt=0.
!        write(0,*) 'adjoint: ADJOINT',ie,nxz+1,nxz+nx*nz,nxt+1,nxt+nrt,size(m),size(d),ibatch,ne,nbatch
!        write(0,*) ' 1 ===',sum(abs(m(nxz+1:nxz+nx*nz))),sum(abs(d(nxt+1:nxt+nrt)))
        if (present(fd_in)) then
          call acfd_born_lop_adj(m(nxz+1:nxz+nx*nz),d(nxt+1:nxt+nrt),rcv_zx(:,:,ie),ie,storewfld,rtm,fd)
        else 
          call acfd_born_lop_adj(m(nxz+1:nxz+nx*nz),d(nxt+1:nxt+nrt),rcv_zx(:,:,ie),ie,storewfld,rtm)
        end if
!        write(0,*) ' 2 ===',sum(abs(m(nxz+1:nxz+nx*nz))),sum(abs(d(nxt+1:nxt+nrt)))

      end do
      
      !if (.not. prestack) then
      !  write(0,*) 'Reducing - m'
      !  call mpi_reduce('m' // itoa(i),m)
      !  write(0,*) 'Finished Reducing - m'
      !end if

      if (present(fd_in)) then
!        call mpi_reduce('fd',fd)
        fd_in=fd
      end if
    end if
    if (allocated(fd)) deallocate(fd)
  end subroutine

!------------------------------------------------

  subroutine  acfd_born_lop_for(m,d,rzx,ie)
    real,dimension(:)  :: m,d
    integer,dimension(:,:) :: rzx
    integer :: iz,ix,it,ir,rind,ind,ie
    integer :: in_tmp
    character(50) :: filenamenr

    ! Open wavefield input file
    filenamenr=trim(filename) // itoa(ie) // '.dat'
    in_tmp=sdr_input(filenamenr)

    ! Initalize Padded model space and Padd
    mm=0;
    do ix=1,nx
    do iz=1,nz
      mm((ix-1+np)*nzz+iz+np)=m((ix-1)*nz+iz)
    end do
    end do
!    write(0,*) 'Called Born Forward 2'

    ! Include some useful factors
    mm = mm*vel2*dz*dx*jslices/dt/dt

    p1=0.;p2=0.;p3=0.
    r1=0.;r2=0.;r3=0.
    do it=2,nt+1,1

      ! Explode model perturbations
      if (mod(it-2,jslices).eq.0) then
        io_buff_dptr=4_dp*dble(floor(real(it-2)/jslices))*dble(nzz*nxx)+1_dp
        call sdr_tmp_read_dptr(in_tmp,s,io_buff_dptr)
        do ix=1,nxx
        do iz=1,nzz
          ind=(ix-1)*nzz+iz
          p3(ind)=p3(ind)-s(ind)*mm(ind)
        end do
        end do
      end if

      ! Compute Laplacian
      lp=0.
      do ix=2,nxx-1
      do iz=2,nzz-1
        ind=(ix-1)*nzz+iz
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1)*nzz+iz-1)
        lp(ind)=lp(ind) + vel2(ind)*c0*p2((ix-1)*nzz+iz+0)*2
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1)*nzz+iz+1)
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1+1)*nzz+iz)
        lp(ind)=lp(ind) + vel2(ind)*c0*p2((ix-1+0)*nzz+iz)*2
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1-1)*nzz+iz)
      end do
      end do

      ! Time-Step
      do ix=1,nxx
      do iz=1,nzz
        ind=(ix-1)*nzz+iz
        ! 2nd derivative wrt time
        p3(ind)=p3(ind)+p2(ind)
        p3(ind)=p3(ind)+p2(ind)
        p3(ind)=p3(ind)-p1(ind)
        ! Laplacian
        p3(ind)=p3(ind)+lp(ind)
      end do
      end do

      ! Absorbing BC's
      call abl_for(p1,p2,p3)

      ! Set Free Surface BC (pressure=0)
      if (freesrf) then
      do ix=1,nxx
      do iz=1,np
        p3((ix-1)*nzz+iz)=0.
      end do
      end do
      end if

      ! Cycle receivers
      do ir=1,nr
        r1(ir)=r2(ir)
        r2(ir)=r3(ir)
        r3(ir)=0.
      end do

      ! Extract Receivers
      do ir=1,nr
        rind=(rzx(ir,2)-1)*nzz+rzx(ir,1)
        r3(ir)=r3(ir)+p3(rind)
      end do

      ! Second order derivative on the fly
      do ir=1,nr
        rind=(it-2)*nr+ir
        d(rind)=d(rind)+   r1(ir)
        d(rind)=d(rind)-2.*r2(ir)
        d(rind)=d(rind)+   r3(ir)
      end do

      ! Cycle
      do ix=1,nxx
      do iz=1,nzz
        ind=(ix-1)*nzz+iz
        p1(ind)=p2(ind)
        p2(ind)=p3(ind)
        p3(ind)=0
      end do
      end do
    end do

    ! Close wavefield input file
    call sdr_input_close(in_tmp)
  end subroutine


  subroutine  acfd_born_lop_adj(m,d,rzx,ie,storewfld,rtm,f_in)
    real,dimension(:)      :: m,d
    integer,dimension(:,:) :: rzx
    logical                :: storewfld,stepsize,rtm
    integer                :: iz,ix,it,ir,rind,ind,ie,iw,iww
    integer                :: in_tmp,in_tmp_wlls,in_tmp_wllr,out_tmp,out_tmp_wll
    real,allocatable,dimension(:) :: tmp,sw,p0,s0,f
    real,optional,dimension(:)    :: f_in
    character(50)                 :: filenamenr
    write(0,*) 'acfd_born_lop_adj'
    
    ! Open wavefield input files
    filenamenr=trim(filename) // itoa(ie) // '.dat'
    in_tmp=sdr_input(filenamenr)
    stepsize=.false.
    if (present(f_in)) then
      stepsize=.true.
      allocate(f(nx*nz)); f=0
      allocate(s0(nx*nz)); s0=0
      allocate(p0(nx*nz)); p0=0
      filenamenr=trim(filename) // itoa(ie) // '_wllfor.dat'
      in_tmp_wlls=sdr_input(filenamenr)
      filenamenr=trim(filename) // itoa(ie) // '_wlladj.dat'
      in_tmp_wllr=sdr_input(filenamenr)
    end if

    ! Open wavefield output file
    if (storewfld) then
      filenamenr=trim(filename) // itoa(ie) // '_adj.dat'
      out_tmp=sdr_output(filenamenr)
      filenamenr=trim(filename) // itoa(ie) // '_wlladj.dat'
      out_tmp_wll=sdr_output(filenamenr)
      allocate(tmp(nx*nz));tmp=0
    end if

    allocate(sw(nx*nz)); sw=0

    ! Padded model space
    mm=0;

    p1=0.;p2=0.;p3=0.
    r1=0.;r2=0.;r3=0.

    do it=nt+1,2,-1

      ! Cycle
      do ix=1,nxx
      do iz=1,nzz
        ind=(ix-1)*nzz+iz
        p3(ind)=p2(ind)
        p2(ind)=p1(ind)
        p1(ind)=0
      end do
      end do

      ! Inject Receivers
      do ir=1,nr
        rind=(it-2)*nr+ir
        if (rtm) then
          r2(ir) = r2(ir) +  d(rind) !?!?!?!?!?!?!?!?!
        else
          r3(ir) = r3(ir) +  d(rind)
          r2(ir) = r2(ir) -2*d(rind)
          r1(ir) = r1(ir) +  d(rind)
        end if
      end do
      do ir=1,nr
        rind=(rzx(ir,2)-1)*nzz+rzx(ir,1)
        p3(rind) = p3(rind)+r3(ir)
      end do

      ! Cycle receivers
      do ir=1,nr
        r3(ir)=r2(ir)
        r2(ir)=r1(ir)
        r1(ir)=0.
      end do

      ! Set Free Surface BC (pressure=0)
      if (freesrf) then
      do ix=1,nxx
      do iz=1,np
        p3((ix-1)*nzz+iz)=0.
      end do
      end do
      end if

      ! Absorbing BC's
      call abl_adj(p1,p2,p3)

      ! Compute Laplacian
      lp=0
      do ix=2,nxx-1
      do iz=2,nzz-1
        ind=(ix-1)*nzz+iz
        lp((ix-1)*nzz+iz-1)=lp((ix-1)*nzz+iz-1) + c1*vel2(ind)*p3(ind)
        lp((ix-1)*nzz+iz+0)=lp((ix-1)*nzz+iz+0) + c0*vel2(ind)*p3(ind)*2
        lp((ix-1)*nzz+iz+1)=lp((ix-1)*nzz+iz+1) + c1*vel2(ind)*p3(ind)
        lp((ix-1+1)*nzz+iz)=lp((ix-1+1)*nzz+iz) + c1*vel2(ind)*p3(ind)
        lp((ix-1+0)*nzz+iz)=lp((ix-1+0)*nzz+iz) + c0*vel2(ind)*p3(ind)*2
        lp((ix-1-1)*nzz+iz)=lp((ix-1-1)*nzz+iz) + c1*vel2(ind)*p3(ind)
      end do
      end do

!      ! Time-Step
      do ix=1,nxx
      do iz=1,nzz
        ind=(ix-1)*nzz+iz
        p2(ind)=p2(ind)+lp(ind)
        p1(ind)=p1(ind)-p3(ind)
        p2(ind)=p2(ind)+p3(ind)
        p2(ind)=p2(ind)+p3(ind)
      end do
      end do
    
      if (storewfld) then
        if (mod(it-2,jslices).eq.0) call sdr_tmp_write(out_tmp,p3)
        do ix=1,nx
        do iz=1,nz
          iw=(ix-1)*nz+iz
          iww=(ix-1+np)*nz+iz+np
          tmp(iw)=p3(iww)
        end do
        end do
        call sdr_tmp_write(out_tmp_wll,tmp)
      end if

      ! Implode reflectors
      if (mod(it-2,jslices).eq.0) then
        io_buff_dptr=4_dp*dble(floor(real(it-2)/jslices))*dble(nzz*nxx)+1_dp
        call sdr_tmp_read_dptr(in_tmp,s,io_buff_dptr)
        do ix=1,nxx
        do iz=1,nzz
          ind=(ix-1)*nzz+iz
          mm(ind) = mm(ind) - s(ind)*p3(ind)
        end do
        end do
      end if

      ! Calculate du_r/dm.du_s/dm
      if (stepsize) then
        io_buff_dptr=4_dp*dble(nx*nz)*(dble(it)-2_dp)+1_dp
        !io_buff_sptr=4*nx*nz*(it-2)+1
        call sdr_tmp_read_dptr(in_tmp_wllr,p0,io_buff_dptr)
        call sdr_tmp_read_dptr(in_tmp_wlls,s0,io_buff_dptr)
        do ix=1,nx
        do iz=1,nz
          iw=(ix-1)*nz+iz
          iww=(ix+np-1)*nzz+iz+np
          !wind=(iw+nw-1+np)*nzz+iw+np
          f(iw)=f(iw) + (s(iww)-s0(iw))*(p3(iww)-p0(iw))
        end do
        end do
      end if
    end do

    ! Include some useful factors
    mm = mm*vel2*dz*dx*jslices/dt/dt

    ! Truncate the Model space
    do ix=1,nx
    do iz=1,nz
      m((ix-1)*nz+iz)=m((ix-1)*nz+iz) + mm((ix-1+np)*nzz+iz+np)
    end do
    end do

    if (stepsize) f_in=f

    ! Close wavefield files
    call sdr_input_close(in_tmp)
    call sdr_input_close(in_tmp_wlls)
    call sdr_input_close(in_tmp_wllr)
    call sdr_output_close(out_tmp)
    call sdr_output_close(out_tmp_wll)
    if (allocated(tmp)) deallocate(tmp)
    if (allocated(sw)) deallocate(sw)
    if (allocated(s0)) deallocate(s0)
    if (allocated(p0)) deallocate(p0)
    if (allocated(f)) deallocate(f)
    write(0,*) 'acfd_born_lop_adj_end'
  end subroutine


  ! Regular FD Kernel to create background wavefields and record receivers
  subroutine acfd_fmd(d,storewfld)
    real,allocatable,dimension(:) :: d
    integer                       :: ie,ieb,nrt,nst
    logical                       :: storewfld
    nrt=nr*nt
    nst=ns*nt

    do ie=ibatch,ne,nbatch
      ieb=ie/nbatch+1
      if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
      call acfd_lop_for(srcdat((ie-1)*nst+1:ie*nst),d((ieb-1)*nrt+1:ieb*nrt),src_zx(:,:,ie),rcv_zx(:,:,ie),ie,storewfld)
    end do
  end subroutine

  subroutine  acfd_lop_for(m,d,szx,rzx,ie,storewfld)
    real,dimension(:)      :: m,d
    integer,dimension(:,:) :: szx,rzx
    integer                :: iz,ix,it,is,ir,sind,rind,ind,ie,iw,iww
    integer                :: out_tmp,out_tmp_wll
    real,allocatable,dimension(:) :: tmp
    character(50)          :: filenamenr
    logical                :: storewfld
    write(0,*) 'acfd_lop_for'

    ! Create wavefield output files
    if (storewfld) then
      filenamenr=trim(filename) // itoa(ie) // '.dat'
      out_tmp=sdr_output(filenamenr)
      filenamenr=trim(filename) // itoa(ie) // '_wllfor.dat'
      out_tmp_wll=sdr_output(filenamenr)
      allocate(tmp(nx*nz));tmp=0
    end if

    p1=0.;p2=0.;p3=0.
    do it=2,nt+1
     
      ! Inject Sources
      do is=1,ns
        ind = (it-2)*ns+is
        sind=(szx(is,2)-1)*nzz+szx(is,1)
        p3(sind)=p3(sind)+vel2(sind)*m(ind)*dz*dx
      end do

      ! Compute Laplacian
      lp=0.
      do ix=2,nxx-1
      do iz=2,nzz-1
        ind=(ix-1)*nzz+iz
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1)*nzz+iz-1)
        lp(ind)=lp(ind) + vel2(ind)*c0*p2((ix-1)*nzz+iz+0)*2
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1)*nzz+iz+1)
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1-1)*nzz+iz)
        lp(ind)=lp(ind) + vel2(ind)*c0*p2((ix-1+0)*nzz+iz)*2
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1+1)*nzz+iz)
      end do
      end do

      ! Time-Step
      do ix=1,nxx
      do iz=1,nzz
        ind=(ix-1)*nzz+iz
        p3(ind)=p3(ind)+p2(ind)
        p3(ind)=p3(ind)+p2(ind)
        p3(ind)=p3(ind)-p1(ind)
        p3(ind)=p3(ind)+lp(ind)
      end do
      end do

      ! Absorbing BC's
      call abl_for(p1,p2,p3)

      ! Set Free Surface BC (pressure=0)
      if (freesrf) then
      do ix=1,nxx
      do iz=1,np
        p3((ix-1)*nzz+iz)=0.
      end do
      end do
      end if

      ! Extract Receivers
      do ir=1,nr
        if (it.gt.nt) cycle
        ind =(it-1)*nr+ir
        rind=(rzx(ir,2)-1)*nzz+rzx(ir,1)
        d(ind)=d(ind)+p3(rind)
      end do
  
      if (storewfld) then
        if (mod(it-2,jslices).eq.0) call sdr_tmp_write(out_tmp,p3)
        do ix=1,nx
        do iz=1,nz
          iww=(ix-1+np)*nzz+iz+np
          iw=(ix-1)*nz+iz
          tmp(iw)=p3(iww)
        end do
        end do
        call sdr_tmp_write(out_tmp_wll,tmp)
      end if

      ! Cycle
      do ix=1,nxx
      do iz=1,nzz
        ind=(ix-1)*nzz+iz
        p1(ind)=p2(ind)
        p2(ind)=p3(ind)
        p3(ind)=0
      end do
      end do

    end do
    ! Close wavefield output file
    call sdr_output_close(out_tmp)
    call sdr_output_close(out_tmp_wll)
    if(allocated(tmp)) deallocate(tmp)
    write(0,*) 'acfd_lop_for_end'
  end subroutine

  subroutine acfd_fmd_adjsrc(phi,m,srcterm,prestack)
    real,allocatable,dimension(:) :: m,phi
    integer                       :: ie,nrt,nst,ieb,nxz
    logical                       :: srcterm,prestack
    nrt=nr*nt
    nst=ns*nt
    !call mpi_cast('phi',phi) ! Already done in main program
    do ie=ibatch,ne,nbatch
      ieb=ie/nbatch+1
      if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
      nxz=(ieb-1)*nx*nz
      if (.not. prestack) nxz=0
      call acfd_lop_for_adjsrc(phi(nxz+1:nxz+nx*nz),m(nxz+1:nxz+nx*nz),ie,srcterm)
    end do
!    call mpi_slock('m_image')
!    if (ibatch.eq.nbatch) call mpi_srelease('m_image')
    !if (.not. prestack) then
      !if (srcterm) call mpi_reduce('grdwrp_src',m)
      !if (.not.srcterm) call mpi_reduce('grdwrp_rcv',m)
    !end if
  end subroutine

  subroutine  acfd_lop_for_adjsrc(phi,m,ie,srcterm)
    real,dimension(:)             :: phi
    real,dimension(:) :: m
    integer                       :: iz,ix,it,ie,iw,ind,iww
    integer                       :: in_tmp_src,in_tmp_wvf,jout
    character(50)                 :: filenamenr
    logical                       :: srcterm
    real,allocatable,dimension(:) :: sw
    real                          :: val
    jout=50
    ! Open wavefield input files
    if (srcterm) then
      filenamenr=trim(filename) // itoa(ie) // '_wlladj.dat'
      in_tmp_src=sdr_input(filenamenr)
      filenamenr=trim(filename) // itoa(ie) // '.dat'
      in_tmp_wvf=sdr_input(filenamenr)
    else
      filenamenr=trim(filename) // itoa(ie) // '_wllfor.dat'
      in_tmp_src=sdr_input(filenamenr)
      filenamenr=trim(filename) // itoa(ie) // '_adj.dat'
      in_tmp_wvf=sdr_input(filenamenr)
    end if
    allocate(sw(nx*nz));sw=0
    p1=0.;p2=0.;p3=0.

    do it=2,nt+1
      ! Inject Sources (make this Laplacian adjoint injection)
      io_buff_dptr=4_dp*dble(nx*nz)*(dble(it)-2_dp)+1_dp ! Read file starting from start
      !io_buff_sptr=4*nx*nz*(it-2)+1 ! Read file starting from start
      !write(0,*) 'it,start_int',it,io_buff_sptr
      call sdr_tmp_read_dptr(in_tmp_src,sw,io_buff_dptr)
      do ix=1,nx
      do iz=1,nz
        iw=(ix-1)*nz+iz
    !    wind=(wzx(iw+nw)-1+np)*nzz+wzx(iw)+np
    !    p3(wind)=p3(wind)-sw(iw)*phi(iw)*vel2(wind)*dz*dx

    !    ! Order reverted & left and right - now pull apart the above coordinates
    !    rtm1((jj-1)*nz+ii  ) = rtm1((jj-1)*nz+ii  ) - 1*rtm((jj-1)*nz+ii)
    !    rtm1((jj-1)*nz+ii-1) = rtm1((jj-1)*nz+ii-1) + 1*rtm((jj-1)*nz+ii) ! 1
    !    rtm1((jj-1)*nz+ii  ) = rtm1((jj-1)*nz+ii  ) - 1*rtm((jj-1)*nz+ii)
    !    rtm1((jj-1)*nz+ii+1) = rtm1((jj-1)*nz+ii+1) + 1*rtm((jj-1)*nz+ii) ! 2
    !    rtm1((jj-1)*nz+ii  ) = rtm1((jj-1)*nz+ii  ) - 1*rtm((jj-1)*nz+ii)
    !    rtm1((jj-2)*nz+ii+0) = rtm1((jj-2)*nz+ii+0) + 1*rtm((jj-1)*nz+ii) ! 3
    !    rtm1((jj-1)*nz+ii  ) = rtm1((jj-1)*nz+ii  ) - 1*rtm((jj-1)*nz+ii)
    !    rtm1((jj-0)*nz+ii+0) = rtm1((jj-0)*nz+ii+0) + 1*rtm((jj-1)*nz+ii) ! 4
        
        !wind=(iw+nw-1+np)*nzz+iw+np

        iw=(ix-1)*nz+iz
        iww=(ix-1+np)*nzz+iz+np
        val=-sw(iw)*phi(iw)*vel2(iww)*dz*dx

        iww=(ix-1+np)*nzz+iz+np
        p3(iww) = p3(iww) - 1*val
        iww=(ix-1+np)*nzz+iz+np-1
        p3(iww) = p3(iww) + 1*val

        iww=(ix-1+np)*nzz+iz+np
        p3(iww) = p3(iww) - 1*val
        iww=(ix-1+np)*nzz+iz+np+1
        p3(iww) = p3(iww) + 1*val

        iww=(ix-1+np)*nzz+iz+np
        p3(iww) = p3(iww) - 1*val
        iww=(ix-2+np)*nzz+iz+np
        p3(iww) = p3(iww) + 1*val

        iww=(ix-1+np)*nzz+iz+np
        p3(iww) = p3(iww) - 1*val
        iww=(ix-0+np)*nzz+iz+np
        p3(iww) = p3(iww) + 1*val

      end do
      end do

      ! Compute Laplacian
      lp=0.
      do ix=2,nxx-1
      do iz=2,nzz-1
        ind=(ix-1)*nzz+iz
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1)*nzz+iz-1)
        lp(ind)=lp(ind) + vel2(ind)*c0*p2((ix-1)*nzz+iz+0)*2
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1)*nzz+iz+1)
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1-1)*nzz+iz)
        lp(ind)=lp(ind) + vel2(ind)*c0*p2((ix-1+0)*nzz+iz)*2
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1+1)*nzz+iz)
      end do
      end do

      ! Time-Step
      do ix=1,nxx
      do iz=1,nzz
        ind=(ix-1)*nzz+iz
        p3(ind)=p3(ind)+p2(ind)
        p3(ind)=p3(ind)+p2(ind)
        p3(ind)=p3(ind)-p1(ind)
        p3(ind)=p3(ind)+lp(ind)
      end do
      end do

      ! Absorbing BC's
      call abl_for(p1,p2,p3)

      ! Set Free Surface BC (pressure=0)
      if (freesrf) then
      do ix=1,nxx
        do iz=1,np
          p3((ix-1)*nzz+iz)=0.
        end do
      end do
      end if

      ! Add Contribution to Gradient
      if (mod(it-2,jslices).eq.0) then
        lp=0.
        lp=p1
        lp=lp-2.*p2
        lp=lp+p3
        io_buff_dptr=4_dp*dble(nzz*nxx)*dble(nslices-((it-2)/jslices+1))+1_dp
        call sdr_tmp_read_dptr(in_tmp_wvf,s,io_buff_dptr)
        lp=lp*s
        do ix=1,nx
        do iz=1,nz  
          m((ix-1)*nz+iz)=m((ix-1)*nz+iz)+lp((ix-1+np)*nzz+iz+np)
        end do
        end do
      end if

      ! Cycle
      do ix=1,nxx
      do iz=1,nzz
        ind=(ix-1)*nzz+iz
        p1(ind)=p2(ind)
        p2(ind)=p3(ind)
        p3(ind)=0
      end do
      end do

    end do
    ! Close wavefield output file
    call sdr_input_close(in_tmp_src)
    call sdr_input_close(in_tmp_wvf)
    if (allocated(sw)) deallocate(sw)
    write(0,*) 'acfd_lop_for_adjsrc_end'
  end subroutine

!------------------------------------------------
  ! Regular FD Kernel to create background wavefields
  subroutine acfd_born_update_background_wfld()
    integer                       :: ie,nst
    nst=ns*nt
    do ie=ibatch,ne,nbatch
      call update_background_wfld(srcdat((ie-1)*nst+1:ie*nst),src_zx(:,:,ie),ie)
    end do
  end subroutine

!------------------------------------------------
  subroutine  update_background_wfld(m,szx,ie)
    real,dimension(:)  :: m
    integer,dimension(:,:) :: szx
    integer :: iz,ix,it,is,sind,ind,ie
    integer :: out_tmp
    character(50) :: filenamenr

    ! Create wavefield output file
    filenamenr=trim(filename) // itoa(ie) // '.dat'
    out_tmp=sdr_output(filenamenr)

    p1=0.;p2=0.;p3=0.
    do it=2,nt+1
     
      ! Inject Sources
      do is=1,ns
        ind = (it-2)*ns+is
        sind=(szx(is,2)-1)*nzz+szx(is,1)
        p3(sind)=p3(sind)+vel2(sind)*m(ind)*dz*dx
      end do

      ! Compute Laplacian
      lp=0.
      do ix=2,nxx-1
      do iz=2,nzz-1
        ind=(ix-1)*nzz+iz
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1)*nzz+iz-1)
        lp(ind)=lp(ind) + vel2(ind)*c0*p2((ix-1)*nzz+iz+0)*2
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1)*nzz+iz+1)
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1-1)*nzz+iz)
        lp(ind)=lp(ind) + vel2(ind)*c0*p2((ix-1+0)*nzz+iz)*2
        lp(ind)=lp(ind) + vel2(ind)*c1*p2((ix-1+1)*nzz+iz)
      end do
      end do

      ! Time-Step
      do ix=1,nxx
      do iz=1,nzz
        ind=(ix-1)*nzz+iz
        p3(ind)=p3(ind)+p2(ind)
        p3(ind)=p3(ind)+p2(ind)
        p3(ind)=p3(ind)-p1(ind)
        p3(ind)=p3(ind)+lp(ind)
      end do
      end do

      ! Absorbing BC's
      call abl_for(p1,p2,p3)

      ! Set Free Surface BC (pressure=0)
      if (freesrf) then
      do ix=1,nxx
      do iz=1,np
        p3((ix-1)*nzz+iz)=0.
      end do
      end do
      end if

      if (mod(it-2,jslices).eq.0) then
        call sdr_tmp_write(out_tmp,p3)
      end if

      ! Cycle
      do ix=1,nxx
      do iz=1,nzz
        ind=(ix-1)*nzz+iz
        p1(ind)=p2(ind)
        p2(ind)=p3(ind)
        p3(ind)=0
      end do
      end do

    end do
    ! Close wavefield output file
    call sdr_output_close(out_tmp)
  end subroutine
!------------------------------------------------

  subroutine  abl_for(p1,p2,p3)
    real,dimension(:) :: p1,p2,p3
    integer :: iz,ix,ind
    do ix=1,nxx
    do iz=1,nzz
      ind=(ix-1)*nzz+iz
      p3(ind) = (p3(ind)-damp(ind)*(damp(ind)*p2(ind)-p1(ind)))/(1+damp(ind))
    end do
    end do
  end subroutine

  subroutine  abl_adj(p1,p2,p3)
    real,dimension(:) :: p1,p2,p3
    integer :: iz,ix,ind
    do ix=1,nxx
    do iz=1,nzz
      ind=(ix-1)*nzz+iz
      p1(ind) = p1(ind) +           damp(ind)*p3(ind)/(1+damp(ind))
      p2(ind) = p2(ind) - damp(ind)*damp(ind)*p3(ind)/(1+damp(ind))
      p3(ind) = p3(ind)/(1+damp(ind))
    end do
    end do
  end subroutine

  subroutine abl_prep(np,dz,dx,dt,nzz,nxx,vmin)
    real           ::    dz,dx,dt,        vmin,kappa,val
    integer        :: np,         nzz,nxx, ind,ix,iz
    damp=0.
    kappa=dt*25*vmin/((np-1)*sqrt(dz*dx))*ddamp/np**2.
    ! Sides
    do iz=np+1,nzz-np
      do ix=1,np
        val=kappa*(real(ix-1))**2
        damp((np-ix+1 -1)*nzz+iz)=val
        damp((np+nx+ix-1)*nzz+iz)=val
      end do
    end do
    do ix=np+1,nxx-np
      do iz=1,np
        val=kappa*(real(iz-1))**2.
        ind=(ix-1)*nzz+np-iz+1
        damp((ix-1)*nzz+np-iz+1 )=val
        damp((ix-1)*nzz+np+nx+iz)=val
      end do
    end do
    ! Corners
    do iz=1,np
    do ix=1,np
        val=kappa*( (sqrt(real(iz-1)**2.+real(ix-1)**2)) )**2.
        damp((np   -ix  )*nzz+np-iz+1 )=val
        damp((nx+np+ix-1)*nzz+nz+np+iz)=val
        damp((nx+np+ix-1)*nzz+np-iz+1 )=val
        damp((np   -ix  )*nzz+nz+np+iz)=val
    end do
    end do
  end subroutine

  ! Add or not to add
  subroutine  adjadd(add,adj,m,d)
    logical       :: add,adj
    real,dimension(:) ::     m,d
    if (.not.add) then
      if (.not.adj) d(:)=0.
      if (     adj) m(:)=0.
    end if
  end subroutine

  ! Integer to String
  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmptmptp
    write(tmptmptp,'(i0)') i
    res = trim(tmptmptp)
  end function

end module

