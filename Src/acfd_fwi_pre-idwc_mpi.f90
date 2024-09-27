!-------------------------------------------------------------------
!   Filename:  acfd_fwi_pre-idwc_mpi
!   Purpose:   Pre-Stack Image Domain Wavefield Tomography in Parallel
!   Developed by: Sjoerd de Ridder
!   Modified by:  Afsaneh Mohammadzaheri
!   Email:     s.deridder@leeds.ac.uk
!   License:   MIT
!----------------------------------------------------

program acdf_fwi_pre_idwc_mpi
  use rsf
  use MPI
  use acfd_born_mod_mpi
  use processing_mod
  use warp2d_mod
  use gradients2d_mod
  implicit none

  !!!!!!!! Declare experiment parameters !!!!!!!!
  type(file)                        :: in,in_rcv_zx,in_src_zx,in_rcv_dat,in_src_dat
  type(file)                        :: out,out_rtm_dat,in_rtm_dat
  type(file)                        :: out_grdwrp,out_grdfwi,out_img,out_grdtot
  integer                           :: nx,nz,nt,ns,nr,ne,ne1,ne2,im,id,iw,ieb,ie,neb,ntaper,ii,jj,iz,ix
  real                              :: dx,ox,dz,oz,dt,ot,dr,or,mute_angle

  real,allocatable,dimension(:)     :: model,pmodel,taper,muter
  real,allocatable,dimension(:)     :: rtm,fd_wflds,rtm1
  real,allocatable,dimension(:)     :: rcvgrd,wrpgrd,total_wrpgrd,totgrd0,totgrd1,search,total_rcvgrd
  real,allocatable,dimension(:)     :: srcdat,rcvdat,rtmdat,rtmdat0,rtmdat1
  real,allocatable,dimension(:,:,:) :: src_zx,rcv_zx
  real,allocatable,dimension(:)     :: rcvres0,rcvres1
  real,allocatable,dimension(:)     :: wrp0,wrp1,wrp00,wrp11, total_wrp1
  real,allocatable,dimension(:)     :: rtmimg0,rtmimg1
  real,allocatable,dimension(:)     :: rtmwrp0,rtmwrp1
  real,allocatable,dimension(:)     :: rtmwrpgrd1,rtmwrpgrd2,fimg1,fimg0,fwrp0,fwrp1
  real,allocatable,dimension(:)     :: phi0,phi1,phi00
  real                              :: rcvres, total_rcvres, Mtotal_rcvres,rtmres,total_rtmres
  real                              :: fit,newres,oldres
  real                              :: dguess
  integer                           :: nguess,iguess,ntrial,nxyz,ntop,nbot

  !!!!!!!! Declare FWI parameters !!!!!!!!
  integer                           :: niter,iter,npad,logfile
  double precision                  :: beta1,beta2,beta,epsscale
  double precision                  :: alpha,alpha0,alpha1,alpha2,alpha1a,alpha1b,alpha2a,alpha2b
  double precision                  :: total_alpha1a,total_alpha1b,total_alpha2a,total_alpha2b
  real                              :: eps_rtm,eps_rcv,perturbp,shiftguess
  character(50)                     :: wflpath,mpipath,filename,fid,ie_string,logfilename,rid,fidm

  !!!!!!!! Declare warp parameters !!!!!!!!
  real                              :: eps_tik,eps_mod,water,waterval
  integer                           :: niter_wrp,itnlim_wrp,nw_mute
  
  ! For muter
  real  :: minz,delz
  integer :: delzn,minzn

  !!!!!!!! Processing parameters !!!!!!!!
  logical                           :: agc,strainreg
  integer                           :: na,nsmooth,voffset

  !real,parameter :: pi=3.14159265359
  real :: pi
  integer :: ibatch, nbatch
  !----------------------------------------------MPI init
  integer :: rank, size_of_cluster, ierror

  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size_of_cluster, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
  !------------------------------------------------------
  ibatch = rank + 1
  nbatch = size_of_cluster
  pi = 3.14159265359


  !------------------------------------ Setup RSF !!!!!!!!
  call sf_init()

  call from_par('fid',fid,' ') ! Data-file identifier
  call from_par('fidm',fidm,' ') ! Data-file identifier for monitor data
  call from_par('rid',rid,' ') ! Run identifier

  !!!!!!!! I/O locations for MPI & Wavefields !!!!!!!!
  mpipath = trim('./scratch/idcfwi-mpi-'//trim(rid)//'/')
  wflpath = trim('./scratch/idcfwi-wfl-'//trim(rid)//'/')

  !!!!!!!! Setup I/O !!!!!!!!
  in  = rsf_input('init_model')
  if (ibatch.eq.nbatch) then
    out = rsf_output('output')
    out_grdwrp = rsf_output('grd-wrp')
    out_grdfwi = rsf_output('grd-fwi')
    out_grdtot = rsf_output('grd-tot')
    out_img = rsf_output('current_image')
  end if
  in_rcv_zx = rsf_input('rcv_zx')
  in_src_zx = rsf_input('src_zx')
  in_src_dat = rsf_input('src_dat')


  !!!!!!!! Read in model dimensions !!!!!!!!
  call from_par(in,'n1',nz,1 ) ! nz : number of samples in z dimension
  call from_par(in,'d1',dz,1.) ! dz : sample spacing in z dimension (m)
  call from_par(in,'o1',oz,0.) ! oz : origin location in z dimension (m)
  call from_par(in,'n2',nx,1 ) ! nx : number of samples in x dimension
  call from_par(in,'d2',dx,1.) ! dx : sample spacing in x dimension (m)
  call from_par(in,'o2',ox,0.) ! ox : origin location in x dimension (m)
  
  !!!!!!!! Read in experiment dimensions !!!!!!!!
  call from_par(in_rcv_zx,'n1',nr,1 ) ! nr : number of receivers
  call from_par(in_rcv_zx,'d1',dr,1.) ! dr : receiver spacing (m)
  call from_par(in_rcv_zx,'o1',or,0.) ! or : location of first receiver in x dimension (m)
  call from_par(in_src_zx,'n1',ns,1 ) ! ns : number of sources
  call from_par(in_rcv_zx,'n3',ne1,1 )! ne1 : number of experiments in source array
  call from_par(in_src_zx,'n3',ne2,1 )! ne2 : number of experiments in receiver array
  if (ne1.ne.ne2) then; write(0,*) "Number of experiments in source and receiver arrays differ"; else; ne=ne1; end if

  !!!!!!!! Read in source signature dimensions !!!!!!!!
  call from_par(in_src_dat,'n2',nt,1 ) ! nt : number of time samples
  call from_par(in_src_dat,'d2',dt,1.) ! dt : sample time period
  call from_par(in_src_dat,'o2',ot,0.) ! ot : start time

  !!!!!!!! Modelling parameters !!!!!!!!
  call from_par('npad',npad,0)  ! npad : thickness of absorbing boundary condition
  call from_par('nsmooth',nsmooth,1)  ! nsmooth : smoothing window for image space gradient
  write(0,*) 'nsmooth',nsmooth

  call from_par('mute_angle',mute_angle,0.)  ! angle with horizon

  !!! For Trial updates
  call from_par('ntrial',ntrial,0)  ! ntrial : number of trial updates

  call from_par('ntop',ntop,0)
  call from_par('nbot',nbot,0)

  call from_par('minz',minz,60.) ! dx : sample spacing in x dimension (m)
  call from_par('delz',delz,60.) ! ox : origin location in x dimension (m)

  delzn=delz/dz
  minzn=minz/dz

  !!!!!!!! Read in FWI parameters !!!!!!!!
  call from_par('niter',  niter,     1   )  ! niter : number of iterations
  call from_par('perturbp',perturbp,1.0)   ! perturbp : trial step-length in search direction (% of model/gradient)
  call from_par('eps_rtm',eps_rtm,   1.  ) ! eps_rtm : weighting of IDWC
  call from_par('eps_rcv',eps_rcv,   1.  ) ! eps_rcv : weighting of receiver residual
  call from_par('agc',agc,.true.)          ! agc : images AGC'ed pre warp?
  call from_par('agclength',na,51)         ! na : agc window length in number of indexes
  call from_par('waterlevel',water,1.)    ! water : waterlevel for phi
  call from_par('dguess',dguess,1.)
  call from_par('nguess',nguess,0)
  call from_par('shiftguess',shiftguess,0.)    ! water : waterlevel for phi
  call from_par('nw_mute',nw_mute,10)      ! nw_mute : number of data muted at each end of well post-warp

  call from_par('strainreg',strainreg,.true.)      ! 

  !!!!!!!! Read in warp/LSQR parameters !!!!!!!!
  call from_par('eps_tik',eps_tik,0.)     ! eps_tik : 2nd order TK regularisation parameter for warp inversion
  call from_par('eps_model',eps_mod,0.)     ! eps_mod : 0th order TK regularisation parameter for warp inversion (damp)
  call from_par('itnlim_warp',itnlim_wrp,25) ! itnlin_wrp : limit on iterations in LSQR in absence of convergence
  call from_par('niter_warp',niter_wrp,5)    ! niter_wrp : number of iterations

  call from_par("ntaper",ntaper,0) ! Length of taper for well adjoint sources (phi)

  if (ibatch.le.mod(ne,nbatch)) then ! defines how many experiments are dealt
    neb=ne/nbatch+1
  else
    neb=ne/nbatch
  end if

  !!!!!!!! Write out FWI parameters !!!!!!!!
  if (ibatch.eq.nbatch) then
    call to_par(out, 'n3',niter+1)
    call to_par(out, 'd3',1.)
    call to_par(out, 'o3',0.)
    call to_par(out,'n1',nz) ! nz : number of samples in z dimension
    call to_par(out,'d1',dz) ! dz : sample spacing in z dimension (m)
    call to_par(out,'o1',oz) ! oz : origin location in z dimension (m)
    call to_par(out,'n2',nx) ! nx : number of samples in x dimension
    call to_par(out,'d2',dx) ! dx : sample spacing in x dimension (m)
    call to_par(out,'o2',ox) ! ox : origin location in x dimension (m)

    call to_par(out_grdfwi, 'n3',niter)
    call to_par(out_grdfwi, 'd3',1.)
    call to_par(out_grdfwi, 'o3',0.)
    call to_par(out_grdfwi,'n1',nz) ! nz : number of samples in z dimension
    call to_par(out_grdfwi,'d1',dz) ! dz : sample spacing in z dimension (m)
    call to_par(out_grdfwi,'o1',oz) ! oz : origin location in z dimension (m)
    call to_par(out_grdfwi,'n2',nx) ! nx : number of samples in x dimension
    call to_par(out_grdfwi,'d2',dx) ! dx : sample spacing in x dimension (m)
    call to_par(out_grdfwi,'o2',ox) ! ox : origin location in x dimension (m)

    call to_par(out_grdwrp, 'n3',niter)
    call to_par(out_grdwrp, 'd3',1.)
    call to_par(out_grdwrp, 'o3',0.)
    call to_par(out_grdwrp,'n1',nz) ! nz : number of samples in z dimension
    call to_par(out_grdwrp,'d1',dz) ! dz : sample spacing in z dimension (m)
    call to_par(out_grdwrp,'o1',oz) ! oz : origin location in z dimension (m)
    call to_par(out_grdwrp,'n2',nx) ! nx : number of samples in x dimension
    call to_par(out_grdwrp,'d2',dx) ! dx : sample spacing in x dimension (m)
    call to_par(out_grdwrp,'o2',ox) ! ox : origin location in x dimension (m)

    call to_par(out_grdtot, 'n3',niter)
    call to_par(out_grdtot, 'd3',1.)
    call to_par(out_grdtot, 'o3',0.)
    call to_par(out_grdtot,'n1',nz) ! nz : number of samples in z dimension
    call to_par(out_grdtot,'d1',dz) ! dz : sample spacing in z dimension (m)
    call to_par(out_grdtot,'o1',oz) ! oz : origin location in z dimension (m)
    call to_par(out_grdtot,'n2',nx) ! nx : number of samples in x dimension
    call to_par(out_grdtot,'d2',dx) ! dx : sample spacing in x dimension (m)
    call to_par(out_grdtot,'o2',ox) ! ox : origin location in x dimension (m)

  end if

  !!!!!!!! Allocations using read in dimensions !!!!!!!!
  allocate(rtmdat(nz*nx*neb))
  rtmdat=0. ! rtmdat : I0 image from well data, read in
  allocate(srcdat(ns*nt*ne))
  srcdat=0. ! srcdat : source signal
  allocate(rcv_zx(nr,2,ne))
  rcv_zx=0. ! rcv_zx : receiver geometry
  allocate(src_zx(ns,2,ne))
  src_zx=0. ! src_zx : source geometry
  allocate(rtm(nz*nx*neb),rtm1(nz*nx*neb))
  rtm=0.;rtm1=0.         ! rtm : RTM image
  allocate(fd_wflds(nz*nx))
  fd_wflds=0.    ! fd_wflds : finite difference of d(us)/ds.d(ur)/ds
  allocate(model(nz*nx),pmodel(nz*nx))
  model=0.; pmodel=0.  ! model/pmodel : model in slowness-squared, initial/perturbed
  allocate(wrp0(nz*nx*neb),wrp00(nz*nx),wrp1(nz*nx*neb),wrp11(nz*nx),fwrp0(nz*nx),fwrp1(nz*nx))
  wrp0=0.; wrp00=0.; wrp1=0.              ! wrp0/1 : warp, initial/perturbed
  allocate(total_wrp1(nz*nx*neb)); total_wrp1=0.
  allocate(phi0(nz*nx*neb),phi1(nz*nx),phi00(nz*nx))
  phi0=0.; phi00=0.; phi1=0.              ! phi0/1 : phi, initial/perturbed
  allocate(rtmimg0(nz*nx),rtmimg1(nz*nx),fimg1(nz*nx),fimg0(nz*nx))
  rtmimg0=0.; rtmimg1=0.; fimg1=0.;fimg0=1;  ! rtmimg0/1 : I1 RTM-derived well image, initial/perturbed
  allocate(rtmwrp0(nz*nx),rtmwrp1(nx*nz))
  rtmwrp0=0.; rtmwrp1=0.  ! rtmwrpgrd1/2 : warped well-derived image
  allocate(rtmwrpgrd1(nx*nz),rtmwrpgrd2(nx*nz))
  rtmwrpgrd1=0.; rtmwrpgrd2=0. ! rtmwrpgrd1/2 : 1st/2nd derivative of warped well-derived image
  allocate(search(nz*nx))
  search=0.                                    ! search : search directon
  allocate(rcvgrd (nz*nx),wrpgrd(nz*nx*neb),total_wrpgrd(nz*nx*neb),total_rcvgrd(nz*nx*neb));
  rcvgrd= 0.; wrpgrd= 0.; total_wrpgrd=0.;total_rcvgrd=0.       ! rcvgrd/wrpgrd : receiver and idwc gradients
  allocate(totgrd0(nz*nx),totgrd1(nz*nx))
  totgrd0=0.; totgrd1=0.       ! totgrd0/1 : previous/current iteration's total gradients
  allocate(taper(nz*nx*neb))
  taper=1.
  allocate(rtmdat0(nz*nx),rtmdat1(nz*nx))
  rtmdat0=0.;rtmdat1=0.
  allocate(muter(nz*nx))

  !!!!!!!! Set up model for FWI !!!!!!!!
  call rsf_read(in,model) ! standard in is velocity model
  model=1./model**2       ! converts model to slowness squared
  if (ibatch.eq.nbatch) then
     call rsf_write(out,1./sqrt(model)) ! model writes out as velocity
     write(0,*) 'Output written'
  end if
  epsscale=0.

  !!!!!!!! Read in geometries and data !!!!!!!!
  call rsf_read(in_rcv_zx,rcv_zx);
  call rsf_read(in_src_zx,src_zx);
  call rsf_read(in_src_dat,srcdat)

  do ie=ibatch,ne,nbatch ! cycles only experiments in current batch
    !write(0,*) src_zx(1,2,ie),mute_angle,dz,dx
    ieb=ie/nbatch+1
    if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
    nxyz=(ieb-1)*nz*nx
  do ix=1,nx
    voffset=(abs(ox+(ix-1)*dx-src_zx(1,2,ie))*tand(mute_angle))/dz+0.5
    !write(0,*) ix,voffset
  do iz=1,voffset
    if (iz.lt.nz) taper(nxyz+(ix-1)*nz+iz)=0.
  end do
  do iz=1,nbot
    if (iz+voffset.lt.nz) taper(nxyz+(ix-1)*nz+iz+voffset)=0.
  end do
  do iz=1,ntaper
    if (iz+voffset+nbot.lt.nz) taper(nxyz+(ix-1)*nz+iz+nbot+voffset)=sind(real(iz-1)*90./real(ntaper))**2.
  end do
  do iz=1,ntaper
    taper(nxyz+(ix-1)*nz+nz-iz+1-ntop)=sind(real(iz-1)*90./real(ntaper))**2.*taper(nxyz+(ix-1)*nz+nz-iz+1-ntop)
  end do
  do iz=1,ntop
    taper(nxyz+(ix-1)*nz+nz-iz+1)=0.*taper(nxyz+(ix-1)*nz+nz-iz+1)
  end do
  end do
  end do

  muter=1.
  write(0,*) 'minzn,delzn'
  write(0,*)  minzn,delzn
  do ix=1,nx
  do iz=1,minzn
    muter((ix-1)*nz+iz)=0.
  end do
  do iz=1,delzn
    muter((ix-1)*nz+iz+minzn)=sind(real(iz-1)*90./real(delzn))**2.
  end do
  end do


  do ie=ibatch,ne,nbatch ! cycles only experiments in current batch
    ieb=ie/nbatch+1
    if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
    write(ie_string,*) ie
    write(0,*) 'reading in rtmdat '
    filename=trim(trim(fid)//'image_e'//trim(adjustl(ie_string))//'.rsf')
    write(0,*) 'reading in rtmdat ',filename
    in_rtm_dat=rsf_input(filename)
    call rsf_read(in_rtm_dat,rtmdat((ieb-1)*nx*nz+1:ieb*nx*nz))
  end do

  !!!!!!!! Setup 'MPI' and 'MPI' parameters !!!!!!!!
  allocate(rcvdat (nr*nt*neb)) ! rcvdat : observed receiver data
  allocate(rcvres0(nr*nt*neb)) ! rcvres0 : inital receiver residual
  allocate(rcvres1(nr*nt*neb)) ! rcvres1 : perturbed receiver residual (used temporarily as masked rcvdat)
  rcvdat=0.

  !!!!!!!! Setup log/progress output !!!!!!!!
  write(0,*) 'ibatch',ibatch
  write(ie_string,*) ibatch
  logfilename=trim('./logs/'//trim(rid)//'_log'//trim(adjustl(ie_string))//'.txt')
  write(0,*) logfilename
  logfile=77+ibatch ! different log file for each batch
  open(logfile,file=logfilename,status='unknown')


  !!!!!!!! Read in receiver data !!!!!!!!
  do ie=ibatch,ne,nbatch ! cycles only experiments in current batch
    ieb=ie/nbatch+1
    if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
    write(ie_string,*) ie
    filename=trim(trim(fidm)//'e'//trim(adjustl(ie_string))//'.rsf')
    write(0,*) 'reading in rcvdat ',filename
    in_rcv_dat=rsf_input(filename)
    call rsf_read(in_rcv_dat,rcvdat((ieb-1)*nr*nt+1:ieb*nr*nt)) ! reads in the receiver data
  end do

  !!!!!!!! Initialize modules !!!!!!!!
  write(0,*) 'wflpath',wflpath
  call MPI_BARRIER(MPI_COMM_WORLD, ierror) ! Do we need this to ensure safe cleaning of the tmp space?
  call acfd_born_init(npad,nz,nx,nt,dz,dx,dt,oz,ox,ot,srcdat,src_zx,rcv_zx,.false.,10,wflpath,filename,neb,ibatch,nbatch)
  call MPI_BARRIER(MPI_COMM_WORLD, ierror) ! Do we need this to ensure safe cleaning of the tmp space?
  call processing_mod_init(nz,nx,nt,dz,dx,dt,oz,ox,ot,src_zx,rcv_zx,model,ibatch,nbatch)
  deallocate(srcdat,src_zx,rcv_zx) ! deallocates parameters only used in modules

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FWI loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rtmres=0.
  oldres=0.
  newres=0.

  do iter=1,niter
    if (ibatch.eq.nbatch) write(logfile,*) '=== Iteration ===',iter

    !!! Update velocity model estimate !!!
    call MPI_BCAST(model, size(model), MPI_REAL, nbatch-1, MPI_COMM_WORLD, ierror)
    call acfd_born_slo2_update(model) ! casts the model, then calculates and stores u_s
    
    !!! Cycle and Clear variables !!!
    totgrd0=totgrd1
    totgrd1=0.
    rcvgrd=0.
    wrpgrd=0.

    !!! Compute and back project receiver residual !!!
    rcvres0=0.
    total_rcvgrd=0.
    call acfd_born_res_grd(rcvdat,rcvres0,rcvgrd) ! stores u_r(residual), outputs receiver residual and gradient
    call MPI_REDUCE(rcvgrd,total_rcvgrd,size(rcvgrd), MPI_REAL, MPI_SUM,nbatch-1,MPI_COMM_WORLD, ierror)
    rcvgrd=total_rcvgrd

    ! event if eps_rcv=0; this has been computed - I may aswell output it
    if (ibatch.eq.nbatch) call rsf_write(out_grdfwi,rcvgrd)

    !!! Computing the receiver residual for monitoring !!!
    rcvres=0.5*eps_rcv*sum(rcvres0**2)
  
    call MPI_REDUCE(rcvres, Mtotal_rcvres, 1, MPI_REAL, MPI_SUM,nbatch-1,MPI_COMM_WORLD, ierror)
    rcvres = Mtotal_rcvres

    if (ibatch.eq.nbatch) write(logfile,*) iter,' ---- Receiver residual m0   :', rcvres !, sum(abs(rcvgrd))
    if (rcvres.ne.rcvres) then
      write(0,*) iter,' --- propagator no longer stable'
    end if

    !!! Produce RTM image !!! (using rcvres1 temporarily as masked rcvdat) 
    rcvres1=0.;rtm1=0.
    call mask_topmute_lop_for(rcvdat,rcvres1)            ! processes the receiver data (and temporarily stores it in rcvres1)
    call acfd_born_lop(.true.,.true.,rtm1,rcvres1,.true.,.true.,.true.) ! ADJOINT, reads u_s, stores u_r(full), outputs RTM image

    ! ADAPT
    do ie=ibatch,ne,nbatch ! cycles only experiments in current batch

      write(0,*) "batch",ie
      ieb=ie/nbatch+1
      if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
      nxyz=(ieb-1)*nz*nx
      write(ie_string,*) ie

      !!!!!!!! Laplacian Filter
      rtm(nxyz+1:nxyz+nz*nx)=0.
      ! rtm=rtm1
      do jj=2,nx-1
      do ii=2,nz-1
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) + 1.*rtm1(nxyz+(jj-0)*nz+ii+0)
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) - 1.*rtm1(nxyz+(jj-1)*nz+ii  )
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) + 1.*rtm1(nxyz+(jj-2)*nz+ii+0)
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) - 1.*rtm1(nxyz+(jj-1)*nz+ii  )
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) + 1.*rtm1(nxyz+(jj-1)*nz+ii+1)
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) - 1.*rtm1(nxyz+(jj-1)*nz+ii  )
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) + 1.*rtm1(nxyz+(jj-1)*nz+ii-1)
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) - 1.*rtm1(nxyz+(jj-1)*nz+ii  )
      end do
      end do
      !call rsf_write(out_rtm,rtm) ! Write the RTM image

      write(0,*) "Write out RTM"
      filename=trim(trim(fid)//'rtmdat2_e'//trim(adjustl(ie_string))//'.rsf')
      out_rtm_dat=rsf_output(filename)
      call to_par(out_rtm_dat,'n3',3 )
      call to_par(out_rtm_dat,'d3',1.)
      call to_par(out_rtm_dat,'o3',0.)
      call to_par(out_rtm_dat,'n1',nz) ! nz : number of samples in z dimension
      call to_par(out_rtm_dat,'d1',dz) ! dz : sample spacing in z dimension (m)
      call to_par(out_rtm_dat,'o1',oz) ! oz : origin location in z dimension (m)
      call to_par(out_rtm_dat,'n2',nx) ! nx : number of samples in x dimension
      call to_par(out_rtm_dat,'d2',dx) ! dx : sample spacing in x dimension (m)
      call to_par(out_rtm_dat,'o2',ox) ! ox : origin location in x dimension (m)
      call rsf_write(out_rtm_dat,rtmdat(nxyz+1:nxyz+nz*nx)*taper(nxyz+1:nxyz+nz*nx))
      call rsf_write(out_rtm_dat,rtm(nxyz+1:nxyz+nz*nx)*taper(nxyz+1:nxyz+nz*nx))
      write(0,*) 'Done!'


      !!! Sample and process well image !!!
      rtmimg0=rtm(nxyz+1:nxyz+nz*nx)*taper(nxyz+1:nxyz+nz*nx); ! Update the monitor image
      rtmdat0=rtmdat(nxyz+1:nxyz+nz*nx)*taper(nxyz+1:nxyz+nz*nx)
      call warp_mod_init(nz,nx,rtmdat0*0,rtmimg0,eps_tik,eps_mod,niter_wrp,itnlim_wrp,strainreg) ! initialises the warp module

      !!! Invert for warp !!!
      rtmwrp0=0.
      fimg0=0.; fit=huge(1.)
      do iguess=1,nguess+1
        wrp00=shiftguess+(iguess-1)*dguess
        call warp_inv(rtmdat0,wrp00,rtmwrp0)
        if (sum(abs(rtmimg0-rtmwrp0)**2.).lt.fit) then
          fit=sum(abs(rtmimg0-rtmwrp0)**2.)
          fimg0=rtmwrp0
          fwrp0=wrp00
        end if
      end do
      rtmwrp0=fimg0
      wrp0(nxyz+1:nxyz+nz*nx)=fwrp0
      call rsf_write(out_rtm_dat,rtmwrp0)

      !write(0,*) "before"
      filename=trim(trim(fid)//'rtmwrp_e'//trim(adjustl(ie_string))//'.rsf')
      out_rtm_dat=rsf_output(filename)
      call to_par(out_rtm_dat,'n3',1 )
      call to_par(out_rtm_dat,'d3',1.)
      call to_par(out_rtm_dat,'o3',0.)
      call to_par(out_rtm_dat,'n1',nz) ! nz : number of samples in z dimension
      call to_par(out_rtm_dat,'d1',dz) ! dz : sample spacing in z dimension (m)
      call to_par(out_rtm_dat,'o1',oz) ! oz : origin location in z dimension (m)
      call to_par(out_rtm_dat,'n2',nx) ! nx : number of samples in x dimension
      call to_par(out_rtm_dat,'d2',dx) ! dx : sample spacing in x dimension (m)
      call to_par(out_rtm_dat,'o2',ox) ! ox : origin location in x dimension (m)
      call rsf_write(out_rtm_dat,wrp0(nxyz+1:nxyz+nz*nx))

      if (sum(wrp0(nxyz+1:nxyz+nz*nx)**2).eq.0) write(logfile,*) 'warp shift is zero, panic!'
      
      !!! Compute gradients of images !!!
      rtmwrpgrd1=0.; rtmwrpgrd2=0.
      call gradient1_2d_for(rtmwrp0,rtmwrpgrd1) ! calculates the 1st derivative of wllwrp w.r.t the well coordinate
      call gradient2_2d_for(rtmwrp0,rtmwrpgrd2) ! calculates the 2nd derivative of wllwrp w.r.t the well coordinate
      
      phi00=rtmwrpgrd1**2-rtmwrpgrd2*(rtmimg0-rtmwrp0)        ! calculates phi's
      waterval=water*sum(abs(phi00))/(nx*nz)
      if (waterval.eq.0.) waterval=1.
      do im=1,nx*nz
         if (phi00(im).gt.0.) then
            phi00(im)=(wrp0(nxyz+im)*rtmwrpgrd1(im))/(phi00(im)+waterval) ! calculates
         else
            phi00(im)=(wrp0(nxyz+im)*rtmwrpgrd1(im))/(phi00(im)-waterval) ! calculates
         end if
      end do
      phi0(nxyz+1:nxyz+nz*nx)=phi00*taper(nxyz+1:nxyz+nz*nx)

      filename=trim(trim(fid)//'rtmphi_e'//trim(adjustl(ie_string))//'.rsf')
      out_rtm_dat=rsf_output(filename)
      call to_par(out_rtm_dat,'n3',1 )
      call to_par(out_rtm_dat,'d3',1.)
      call to_par(out_rtm_dat,'o3',0.)
      call to_par(out_rtm_dat,'n1',nz) ! nz : number of samples in z dimension
      call to_par(out_rtm_dat,'d1',dz) ! dz : sample spacing in z dimension (m)
      call to_par(out_rtm_dat,'o1',oz) ! oz : origin location in z dimension (m)
      call to_par(out_rtm_dat,'n2',nx) ! nx : number of samples in x dimension
      call to_par(out_rtm_dat,'d2',dx) ! dx : sample spacing in x dimension (m)
      call to_par(out_rtm_dat,'o2',ox) ! ox : origin location in x dimension (m)
      call rsf_write(out_rtm_dat,phi0(nxyz+1:nxyz+nz*nx))

    end do

    if (iter.eq.1) then
      if (rcvres.ne.0..and.eps_rtm.ne.0.) then
        eps_rtm=eps_rtm*rcvres/(0.5*sum(wrp0(1:nz*nx)**2))
        call MPI_BCAST(eps_rtm, 1, MPI_REAL, nbatch-1,MPI_COMM_WORLD,ierror)
      end if
    end if
    rtmres=eps_rtm*0.5*sum(wrp0**2)

    write(0,*) 'rank (rtmres)',eps_rtm,rtmres
    call MPI_REDUCE(rtmres, total_rtmres, 1, MPI_REAL, MPI_SUM, nbatch-1,MPI_COMM_WORLD, ierror)
    rtmres = total_rtmres
    write(0,*) 'rank (total rtmres)',rtmres

    oldres=rtmres+rcvres
    if (ibatch.eq.nbatch) then
      write(logfile,*) iter,' ---- Image    residual m0   :', rtmres
      write(logfile,*) iter,' ---- Total    residual m0   :', oldres
    end if


    fwrp0=0.
    wrp1=wrp0
    !phi00=0.
    write(0,*) 'Reduce wrp'
    call MPI_REDUCE(wrp1, total_wrp1, size(wrp1), MPI_REAL, MPI_SUM,nbatch-1,MPI_COMM_WORLD, ierror)
    wrp1 = total_wrp1

    if (ibatch.eq.nbatch) then
      do ieb=1,neb ! cycles only experiments in current batch
        write(0,*) 'ieb - here',ieb,size(fwrp0),size(wrp1)
        nxyz=(ieb-1)*nz*nx
        fwrp0=fwrp0+wrp1(nxyz+1:nxyz+nz*nx)
      end do
      !call rsf_write(out_wrp,fwrp0)
    end if
    wrp1=0.

    !if (iter.eq.1) call mpi_cast_datum('eps',eps_rtm)
    
    !!! Compute wavefields for the gradient !!!
    wrpgrd=0.
    call acfd_fmd_adjsrc(phi0,wrpgrd,.true. ,.true.)                ! calculates the source term of the warp gradient
    call acfd_fmd_adjsrc(phi0,wrpgrd,.false.,.true.)                ! calculates the receiver term of the warp gradient

    do ie=ibatch,ne,nbatch ! cycles only experiments in current batch
      write(0,*) "batch",ie
      ieb=ie/nbatch+1
      if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
      nxyz=(ieb-1)*nz*nx
      filename=trim(trim(fid)//'wrpgrd_e'//trim(adjustl(ie_string))//'.rsf')
      out_rtm_dat=rsf_output(filename)
      call to_par(out_rtm_dat,'n3',1 )
      call to_par(out_rtm_dat,'d3',1.)
      call to_par(out_rtm_dat,'o3',0.)
      call to_par(out_rtm_dat,'n1',nz) ! nz : number of samples in z dimension
      call to_par(out_rtm_dat,'d1',dz) ! dz : sample spacing in z dimension (m)
      call to_par(out_rtm_dat,'o1',oz) ! oz : origin location in z dimension (m)
      call to_par(out_rtm_dat,'n2',nx) ! nx : number of samples in x dimension
      call to_par(out_rtm_dat,'d2',dx) ! dx : sample spacing in x dimension (m)
      call to_par(out_rtm_dat,'o2',ox) ! ox : origin location in x dimension (m)
      call rsf_write(out_rtm_dat,wrpgrd(nxyz+1:nxyz+nz*nx))
    end do

    total_wrpgrd=0.
    call MPI_REDUCE(wrpgrd, total_wrpgrd, size(wrpgrd), MPI_REAL,MPI_SUM,nbatch-1,MPI_COMM_WORLD, ierror)
    wrpgrd = total_wrpgrd


    if (ibatch.eq.nbatch) then
      wrpgrd=-wrpgrd
      fwrp0=0.
      do ieb=1,neb ! cycles only experiments in current batch
        nxyz=(ieb-1)*nz*nx
        fwrp0=fwrp0+wrpgrd(nxyz+1:nxyz+nz*nx)
      end do
      !if (nsmooth.gt.1) then
      !  call smooth2d(fwrp0,totgrd1,nz,nx,nsmooth)
      !  fwrp0=totgrd1
      !  totgrd1=0.
      !end if
      call rsf_write(out_grdwrp,fwrp0)
    
      !!! Compute search direction using Polak-Ribiere !!!
      beta1=0.0d0; beta2=0.0d0; beta=0.0d0
      totgrd1=eps_rcv*rcvgrd+eps_rtm*fwrp0 ! adds together receiver and idwc gradients (fwrp0 is a temp array)
      totgrd1=totgrd1*muter
      call rsf_write(out_grdtot,totgrd1)
      do im=1,nz*nx
        beta1=beta1+totgrd1(im)*(totgrd1(im)-totgrd0(im))
        beta2=beta2+totgrd0(im)* totgrd0(im)
      end do
      if (beta2.ne.0.0d0) beta = beta1/beta2
      search = beta*search-totgrd1

      alpha0=dble(perturbp/100.0)*dble(sum(abs(model)))/dble(sum(abs(search)))

      pmodel=model+alpha0*search
    end if

    !!! Update velocity model !!!
    call MPI_BCAST(pmodel, size(pmodel), MPI_REAL, nbatch-1, MPI_COMM_WORLD, ierror)
    call acfd_born_slo2_update(pmodel) ! casts pmodel only

    !!! Compute new receiver residual !!!
    rcvres1=0.
    call acfd_born_res(rcvdat,rcvres1,.true.) ! saves u_s
    
    !!! Computing the receiver residual for monitoring !!!
    rcvres=0.5*eps_rcv*sum(rcvres1**2)
    call MPI_REDUCE(rcvres, total_rcvres, 1,MPI_REAL,MPI_SUM,nbatch-1,MPI_COMM_WORLD, ierror)
    rcvres = total_rcvres


    if (ibatch.eq.nbatch) write(logfile,*) iter,' ---- Receiver residual m0+dm:', rcvres
    if (rcvres.ne.rcvres) then
      write(0,*) iter,' --- propagator no longer stable - PANIC---'
    end if

    !!! Compute finite difference of receiver residual
    rcvres1=(rcvres1-rcvres0)

    !!! Calculate receiver contribution to step-size !!!
    alpha1a=0.0d0; alpha2a=0.0d0
    alpha1b=0.0d0; alpha2b=0.0d0; alpha=0.0d0
    do id=1,nr*nt*neb
      alpha1a=alpha1a+eps_rcv*rcvres1(id)*rcvres0(id)
      alpha2a=alpha2a+eps_rcv*rcvres1(id)*rcvres1(id)
    end do

    !!! Produce new RTM image !!! (re-use rcvres1)
    rtm1=0.; rcvres1=0.
    call mask_topmute_lop_for(rcvdat,rcvres1)
    call acfd_born_lop(.true.,.true.,rtm1,rcvres1,.false.,.true.,.true.,fd_wflds) ! uses u_s as saved, calcs stepsize component

    do ie=ibatch,ne,nbatch ! cycles only experiments in current batch
      write(ie_string,*) ie
      ieb=ie/nbatch+1
      if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
      write(0,*) "batch",ie
      !!!!!!!! Laplacian Filter
      nxyz=(ieb-1)*nz*nx
      rtm(nxyz+1:nxyz+nz*nx)=0.
      do jj=2,nx-1
      do ii=2,nz-1
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) + 1.*rtm1(nxyz+(jj-0)*nz+ii+0)
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) - 1.*rtm1(nxyz+(jj-1)*nz+ii  )
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) + 1.*rtm1(nxyz+(jj-2)*nz+ii+0)
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) - 1.*rtm1(nxyz+(jj-1)*nz+ii  )
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) + 1.*rtm1(nxyz+(jj-1)*nz+ii+1)
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) - 1.*rtm1(nxyz+(jj-1)*nz+ii  )
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) + 1.*rtm1(nxyz+(jj-1)*nz+ii-1)
        rtm(nxyz+(jj-1)*nz+ii)=rtm(nxyz+(jj-1)*nz+ii) - 1.*rtm1(nxyz+(jj-1)*nz+ii  )
      end do
      end do

      !!! Sample and process well image !!!
      rtmimg1=rtm(nxyz+1:nxyz+nz*nx)*taper; ! Update the monitor image
      rtmdat1=rtmdat(nxyz+1:nxyz+nz*nx)*taper
      call warp_mod_init(nz,nx,rtmdat1*0,rtmimg1,eps_tik,eps_mod,niter_wrp,itnlim_wrp,strainreg) ! initialises the warp module

      !!! Invert for warp !!!
      rtmwrp1=0.
      fimg1=0.; fit=huge(1.)
      do iguess=1,nguess+1
        wrp11=shiftguess+(iguess-1)*dguess
        call warp_inv(rtmdat1,wrp11,rtmwrp1)
        if (sum(abs(rtmimg1-rtmwrp1)**2.).lt.fit) then
          fit=sum(abs(rtmimg1-rtmwrp1)**2.)
          fimg1=rtmwrp1
          fwrp1=wrp11
        end if
      end do
      rtmwrp1=fimg1
      wrp1(nxyz+1:nxyz+nz*nx)=fwrp1
      if (sum(wrp1(nxyz+1:nxyz+nz*nx)**2).eq.0) write(logfile,*) 'warp shift is zero, panic!'

    end do

    do ie=ibatch,ne,nbatch ! cycles only experiments in current batch
      write(ie_string,*) ie
      write(0,*) "rtm_trial batch",ie
      ieb=ie/nbatch+1
      if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
      nxyz=(ieb-1)*nz*nx
      filename=trim(trim(fid)//'rtm_trial_e'//trim(adjustl(ie_string))//'.rsf')
      out_rtm_dat=rsf_output(filename)
      call to_par(out_rtm_dat,'n3',1 )
      call to_par(out_rtm_dat,'d3',1.)
      call to_par(out_rtm_dat,'o3',0.)
      call to_par(out_rtm_dat,'n1',nz) ! nz : number of samples in z dimension
      call to_par(out_rtm_dat,'d1',dz) ! dz : sample spacing in z dimension (m)
      call to_par(out_rtm_dat,'o1',oz) ! oz : origin location in z dimension (m)
      call to_par(out_rtm_dat,'n2',nx) ! nx : number of samples in x dimension
      call to_par(out_rtm_dat,'d2',dx) ! dx : sample spacing in x dimension (m)
      call to_par(out_rtm_dat,'o2',ox) ! ox : origin location in x dimension (m)
      call rsf_write(out_rtm_dat,rtm(nxyz+1:nxyz+nz*nx))
    end do

    do ie=ibatch,ne,nbatch ! cycles only experiments in current batch
      write(ie_string,*) ie
      ieb=ie/nbatch+1
      if (mod(ie,nbatch).eq.0) ieb=ie/nbatch
      nxyz=(ieb-1)*nz*nx
      filename=trim(trim(fid)//'wrp_trial_e'//trim(adjustl(ie_string))//'.rsf')
      out_rtm_dat=rsf_output(filename)
      call to_par(out_rtm_dat,'n3',1 )
      call to_par(out_rtm_dat,'d3',1.)
      call to_par(out_rtm_dat,'o3',0.)
      call to_par(out_rtm_dat,'n1',nz) ! nz : number of samples in z dimension
      call to_par(out_rtm_dat,'d1',dz) ! dz : sample spacing in z dimension (m)
      call to_par(out_rtm_dat,'o1',oz) ! oz : origin location in z dimension (m)
      call to_par(out_rtm_dat,'n2',nx) ! nx : number of samples in x dimension
      call to_par(out_rtm_dat,'d2',dx) ! dx : sample spacing in x dimension (m)
      call to_par(out_rtm_dat,'o2',ox) ! ox : origin location in x dimension (m)
      call rsf_write(out_rtm_dat,wrp1(nxyz+1:nxyz+nz*nx))
    end do


    rtmres=0.5*eps_rtm*sum(wrp1**2)
    call MPI_REDUCE(rtmres, total_rtmres,1,MPI_REAL,MPI_SUM,nbatch-1,MPI_COMM_WORLD, ierror)
    rtmres = total_rtmres


    if (ibatch.eq.nbatch) then
    write(logfile,*) iter,' ---- Image    residual m0+dm:', rtmres
    write(logfile,*) iter,' ---- Total    residual m0+dm:', rtmres+rcvres
    end if
    newres=rtmres+rcvres
    ! (Perhaps we don't need this)
    call MPI_BCAST(newres, 1, MPI_REAL, nbatch-1,MPI_COMM_WORLD,ierror)

    wrp1=(wrp1-wrp0)

    !!! Calculate warp contribution to step-size (image terms) !!!
    
    do iw=1,nx*nz*neb
      alpha1b=alpha1b+eps_rtm*wrp1(iw)*wrp0(iw)
      alpha2b=alpha2b+eps_rtm*wrp1(iw)*wrp1(iw)
    end do
    
    !!! Reduce step length components !!!
    call MPI_REDUCE(alpha1a,total_alpha1a,1,MPI_DOUBLE_PRECISION,MPI_SUM,nbatch-1,MPI_COMM_WORLD, ierror)
    alpha1a = total_alpha1a
    call MPI_REDUCE(alpha2a,total_alpha2a,1,MPI_DOUBLE_PRECISION,MPI_SUM,nbatch-1,MPI_COMM_WORLD, ierror)
    alpha2a = total_alpha2a
    call MPI_REDUCE(alpha1b,total_alpha1b,1,MPI_DOUBLE_PRECISION,MPI_SUM,nbatch-1,MPI_COMM_WORLD, ierror)
    alpha1b = total_alpha1b
    call MPI_REDUCE(alpha2b,total_alpha2b,1,MPI_DOUBLE_PRECISION,MPI_SUM,nbatch-1,MPI_COMM_WORLD, ierror)
    alpha2b = total_alpha2b


    if (ibatch.eq.nbatch) then
      !!! Calculate step size !!!
      alpha1=-(alpha1a+alpha1b)
      alpha2=+(alpha2a+alpha2b)
      if (alpha2.ne.0.0d0) alpha=alpha1/alpha2*alpha0
      write(logfile,*) iter, ' ---- Alpha 1:                  ',real(alpha1a),real(alpha1b)
      write(logfile,*) iter, ' ---- Alpha 2:                  ',real(alpha2a),real(alpha2b)
      write(logfile,*) iter, ' ---- Step (versus trial step): ',real(alpha)  ,real(alpha0)

      model=model+alpha*search
      call rsf_write(out,1./sqrt(model))
      write(0,*) 'Written out model'
    end if
  end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FWI loop end !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!! Clean up and close !!!!!!!!
  deallocate(rcvdat)
  deallocate(model,pmodel)
  deallocate(rtm,rtm1)
  deallocate(rcvres0,rcvres1)
  deallocate(wrp0,wrp1,wrp00,wrp11)
  deallocate(rtmimg0,rtmimg1,fimg0,fimg1)
  deallocate(rtmwrp0,rtmwrp1)
  deallocate(rtmwrpgrd1,rtmwrpgrd2)
  deallocate(phi00)
  deallocate(phi0)
  deallocate(phi1)
  deallocate(rcvgrd,wrpgrd)
  deallocate(totgrd0,totgrd1)
  deallocate(search)
  deallocate(fd_wflds)
  deallocate(taper)
  deallocate(fwrp0)
  deallocate(fwrp1)
  deallocate(rtmdat0,rtmdat1)
  deallocate(muter)
  if (ibatch.eq.nbatch) deallocate(rtmdat)

  call acfd_born_close()
  call processing_mod_close()
  call warp_mod_close()
  call sf_close()
  close(logfile)

  call MPI_FINALIZE(ierror)

  call exit(0)
  
end program
