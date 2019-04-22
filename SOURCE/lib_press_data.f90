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

<<<<<<< HEAD
=======

  private ::press_io_constructor,press_io_read,&
       press_io_write,press_info,press_destroy


>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
  type PRESSdata
     character(len=4) ::version
     integer          ::nsamples
     integer          ::nsensors
     real             ::Fs=1.d0,Pa=1.d0
     real, dimension(:,:),   allocatable ::p
     real, dimension(:,:),   allocatable ::x

     character(len=500)  ::comments = ""

   contains

<<<<<<< HEAD
     procedure  :: create => ccor
     procedure  :: read_bin => press_io_read
     procedure  :: write_bin => press_io_write
=======
     procedure  :: create => press_io_constructor
     procedure  :: read_bin => press_io_read
     procedure  :: get_fluctuations => press_get_fluc
     procedure  :: write_bin => press_io_write
     procedure  :: read_data => press_io_read
     procedure  :: write_data => press_io_write
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
     procedure  :: print_info => press_info
     procedure  :: destroy => press_destroy

  end type PRESSdata


contains

  subroutine press_info(p)
    class(PRESSdata)              ::p

<<<<<<< HEAD
    write(06,*)"Numbers of sensors : ",P.nsensors
    write(06,*)"Numbers of samples : ",P.nsamples
    write(06,*)"Sampling frequency : ",P.Fs
=======
    write(06,*)"Numbers of sensors : ",P%nsensors
    write(06,*)"Numbers of samples : ",P%nsamples
    write(06,*)"Sampling frequency : ",P%Fs
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0

    return
  end subroutine press_info

<<<<<<< HEAD
  subroutine ccor(p,nsamples,nsensors,Fs,Pa,com)
=======
  subroutine press_io_constructor(p,nsamples,nsensors,Fs,Pa,com)
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
    class(PRESSdata)              ::p
    integer                       ::nsamples,nsensors
    real                          ::Fs,Pa
    character(len=*),optional     ::com
    !---------------------------------------
<<<<<<< HEAD
    p.version = cur_version
    p.nsamples = nsamples
    p.nsensors = nsensors
    allocate(p.p(nsamples,nsensors))
    allocate(p.x(nsensors,ndim))

    p.Fs=Fs
    p.Pa=Pa

    if (present(com)) then
       p.comments=trim(com)
    end if

  end subroutine ccor
=======
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
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0

  subroutine press_destroy(this)
    class(PRESSdata)              ::this
    !---------------------------------------


<<<<<<< HEAD
    if (allocated(this.p)) deallocate(this.p)
    if (allocated(this.x)) deallocate(this.x)

  end subroutine press_destroy

=======
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

>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
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
<<<<<<< HEAD
    read(fid)p.version
    if (p.version==cur_version) then
<<<<<<< HEAD
       read(fid)p.nsensors,p.nsamples,p.Fs
=======
       read(fid)p.nsensors,p.nsamples,p.Fs,p.Pa
>>>>>>> b8e6bcaa7a21208163788d5d7a30a1fa9b2cf1e8
       read(fid)p.comments

       if (.not.allocated(p.p)) allocate(p.p(p.nsamples,p.nsensors))
       if (.not.allocated(p.x)) allocate(p.x(p.nsensors,ndim))

       read(fid)p.x
       read(fid)p.p
    else
       write(06,'(a,a)')'PRESSdata error : version unknown'
       return
=======
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
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
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
<<<<<<< HEAD
<<<<<<< HEAD
    write(fid)p.nsensors,p.nsamples,p.Fs
=======
    write(fid)p.nsensors,p.nsamples,p.Fs,p.Pa
>>>>>>> b8e6bcaa7a21208163788d5d7a30a1fa9b2cf1e8
    write(fid)p.comments
    write(fid)p.x
    write(fid)p.p
=======
    write(fid)p%nsensors,p%nsamples,p%Fs,p%Pa
    write(fid)p%comments
    write(fid)p%x
    write(fid)p%p
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
    close(fid)

    return
  end subroutine press_io_write

end MODULE LIB_PRESS_DATA
