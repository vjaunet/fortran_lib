MODULE LIB_PRESS_DATA
  implicit none

  !=========================================================
  !
  ! Kulite pressure routines and container
  !*  author : Vincent Jaunet
  !*  License: GPL v3.0
  !
  !=========================================================

  integer, private :: i,j,ic,is
  integer, private, parameter ::ndim=2
  character(len=4),private ::cur_version="v0.1"

  private ::press_io_constructor,press_io_read,&
       press_io_write,press_info,press_destroy

  type PRESSdata
     character(len=4) ::version
     integer          ::nsamples
     integer          ::nsensors
     real             ::Fs=1.d0,Pa=1.d0
     real, dimension(:,:),   allocatable ::p
     real, dimension(:,:),   allocatable ::x

     character(len=500)  ::comments = ""

   contains

     procedure  :: create => press_io_constructor
     procedure  :: read_bin => press_io_read
     procedure  :: get_fluctuations => press_get_fluc
     procedure  :: write_bin => press_io_write
     procedure  :: read_data => press_io_read
     procedure  :: write_data => press_io_write
     procedure  :: print_info => press_info
     procedure  :: destroy => press_destroy

  end type PRESSdata


contains

  subroutine press_info(p)
    class(PRESSdata)              ::p

    write(06,*)"Numbers of sensors : ",P%nsensors
    write(06,*)"Numbers of samples : ",P%nsamples
    write(06,*)"Sampling frequency : ",P%Fs

    return
  end subroutine press_info

  subroutine press_io_constructor(p,nsamples,nsensors,Fs,Pa,com)
    class(PRESSdata)              ::p
    integer                       ::nsamples,nsensors
    real                          ::Fs,Pa
    character(len=*),optional     ::com
    !---------------------------------------
    p%version = cur_version
    p%nsamples = nsamples
    p%nsensors = nsensors
    allocate(p%p(nsamples,nsensors))
    allocate(p%x(nsensors,ndim))

    p%Fs=Fs
    p%Pa=Pa

    if (present(com)) then
       p%comments=trim(com)
    end if

  end subroutine press_io_constructor

  subroutine press_destroy(this)
    class(PRESSdata)              ::this
    !---------------------------------------


    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%x)) deallocate(this%x)

  end subroutine press_destroy

  subroutine press_get_fluc(this)
    class(PRESSdata)              ::this
    integer                       ::ic
    !----------------------------------------------

    do ic=1,this%nsensors
       this%p(:,ic) =  this%p(:,ic) -&
            sum(this%p(:,ic))/this%nsamples
    end do

    return
  end subroutine press_get_fluc

  subroutine press_io_read(p,ifile)
    class(PRESSdata)              ::p
    character(len=*)              ::ifile
    logical                       ::file_exists
    integer, parameter            ::fid=1515
    !----------------------------------------------

    inquire(file=trim(ifile), EXIST=file_exists)
    if (.not.file_exists) then
       write(06,'(a,a)')'PRESSdata error : cannot find ',ifile
       return
    end if

    !file exists let's read it
    open(unit=fid,file=trim(ifile),action='read',access='stream',status='old')
    read(fid)p%version
    if (p%version==cur_version) then
       read(fid)p%nsensors,p%nsamples,p%Fs,p%Pa
       read(fid)p%comments

       if (.not.allocated(p%p)) allocate(p%p(p%nsamples,p%nsensors))
       if (.not.allocated(p%x)) allocate(p%x(p%nsensors,ndim))

       read(fid)p%x
       read(fid)p%p
    else
       write(06,'(a,a)')'PRESSdata error : version unknown'
       STOP
    end if
    close(fid)

    return
  end subroutine press_io_read

  subroutine press_io_write(p,ofile)
    class(PRESSdata)              ::p
    character(len=*)              ::ofile
    integer, parameter            ::fid=1515
    !------------------------------------------------

    open(unit=fid, access='stream',file=trim(ofile),status='unknown')
    write(fid)cur_version
    write(fid)p%nsensors,p%nsamples,p%Fs,p%Pa
    write(fid)p%comments
    write(fid)p%x
    write(fid)p%p
    close(fid)

    return
  end subroutine press_io_write

end MODULE LIB_PRESS_DATA
