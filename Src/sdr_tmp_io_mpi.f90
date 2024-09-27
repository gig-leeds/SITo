module sdr_tmp_io_mpi
  use rsf
  implicit none

  private

  public :: tmp_io_init,tmp_io_close
  public :: sdr_input,sdr_input_close
  public :: sdr_output,sdr_output_close
  !public :: sdr_tmp_read,sdr_tmp_write
  public :: sdr_tmp_write
  public :: sdr_tmp_read_sptr,sdr_tmp_read_dptr

  character(50)  :: path
  integer        :: unit_cntr,nbatch,ibatch
  integer,parameter :: dp=8

  !interface sdr_tmp_read
  !  module procedure sdr_tmp_read_sptr
  !  module procedure sdr_tmp_read_dptr
  !end interface

  contains

!------------------------------------------------
  subroutine tmp_io_init(path_in, ibatch_in, nbatch_in)

    integer ibatch_in, nbatch_in
    character(50) :: path_in
    character(100) :: makedirectory*45,rmdirectory*45
    path=path_in
    ibatch=ibatch_in
    nbatch=nbatch_in
    !call from_par('nbatch',nbatch,1)
    !call from_par('ibatch',ibatch,1)
    call tmp_io_close()

    ! Set unit counter
    unit_cntr=900

    ! Ensure folder is made
    if (ibatch.eq.nbatch) then
      rmdirectory = 'rm -rf ' // trim(path)
      call system(rmdirectory)
      makedirectory = 'mkdir -p ' // trim(path)
      call system(makedirectory)
    end if
  end subroutine
  subroutine tmp_io_close()
    character(100) :: rmdirectory*45
    if (ibatch.eq.nbatch) then
      rmdirectory = 'rm -rf ' // trim(path)
      call system(rmdirectory)
    end if
  end subroutine

  integer function sdr_input(filename)
    character(50) :: filename
    character(100) :: fullpath
    integer :: unitnr
    ! Determine unit counter
    unit_cntr=unit_cntr+1
    unitnr=unit_cntr
    fullpath = trim(path) // trim(filename)
    !write(0,*) fullpath
    !write(0,*) "Trying to Open File",unitnr,fullpath
    !open(unit=unitnr,file=fullpath,form='unformatted',status='old',action='read',access='direct')
    open(unit=unitnr,file=fullpath,form='unformatted',status='old',action='read',access='stream')
    !write(0,*) "Opened File"
    sdr_input=unitnr
  end function
  subroutine sdr_input_close(unitnr)
    integer :: unitnr
    close(unit=unitnr)
  end subroutine

  integer function sdr_output(filename)
    character(50) :: filename
    character(100) :: fullpath
    integer :: unitnr
    ! Determine unit counter
    unit_cntr=unit_cntr+1
    unitnr=unit_cntr
    fullpath = trim(path) // trim(filename)
    !write(0,*) fullpath
    open(unit=unitnr,file=fullpath,form='unformatted',status='replace',action='write',access='stream')
    sdr_output=unitnr
  end function
  subroutine sdr_output_close(unitnr)
    integer :: unitnr
    close(unit=unitnr)
  end subroutine

  subroutine sdr_tmp_write(unitnr,buff,bytenr)
    integer           ::   unitnr
    integer,optional  ::               bytenr
    real,dimension(:) ::          buff
    if (present(bytenr)) then
      write(unitnr,pos=bytenr) buff  ! Untested
    else
    !  write(0,*) "Written",size(buff),'bytes to file'
      write(unitnr) buff
    end if
  end subroutine

  subroutine sdr_tmp_read_sptr(unitnr,buff,bytenr)
    integer           ::  unitnr
    integer,optional  ::              bytenr
    real,dimension(:) ::         buff
    if (present(bytenr)) then
      read(unitnr,pos=bytenr) buff
    else
      read(unitnr) buff
    end if
  end subroutine

  subroutine sdr_tmp_read_dptr(unitnr,buff,bytenr)
    integer           ::  unitnr
    integer(dp),optional  ::              bytenr
    real,dimension(:) ::         buff
    if (present(bytenr)) then
      read(unitnr,pos=bytenr) buff
    else
      read(unitnr) buff
    end if
  end subroutine
  
end module
