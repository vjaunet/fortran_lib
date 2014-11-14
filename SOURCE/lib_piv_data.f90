module lib_piv_data
  use lib_stat
  use lib_pod
  implicit none

  !=================Specification=============================
  !
  !
  ! PIV data container
  !
  !
  !
  !
  !===========================================================

  integer, private :: i,j,ic,is

  type PIVdata
     character ::typeofgrid="C"
     integer   ::nx,ny,nsamples,ncomponent,pixel_step
     real      ::dx,x0,dy,y0 !for scaling
     integer   ::ntheta,nr
     real      ::dr,dtheta
     real      ::z_pos=0,fs=1
     real, dimension(:,:,:,:), allocatable ::u
     real, dimension(:,:,:),   allocatable ::x

     !stat containers
     real, dimension(:,:,:),   allocatable ::ustat
     real, dimension(:,:,:),   allocatable ::w

   contains
     procedure :: read_bin  => piv_io_read
     procedure :: write_bin => piv_io_write
     procedure :: cal_stats => piv_stats
     procedure :: set_x     => piv_set_x
     procedure :: destroy   => piv_destroy
     procedure :: replace_outlier => piv_replace_outlier

  end type PIVdata

contains

  subroutine piv_replace_outlier(datapiv,method,wsize)
    class(PIVdata)                     ::datapiv
    character(len=*)                   ::method
    real                               ::um
    integer                            ::ii,jj,i0
    real                               ::nval
    integer    ,optional               ::wsize
    !-----------------------------------------------

    select case(method)

    case ("POD")

       call gappypod(datapiv.u,datapiv.w,int(10),real(1e-8))

    case("AVG")

       if (.not.PRESENT(wsize)) wsize = 3;

       i0=wsize/2

       do is=1,datapiv.nsamples
          do j=1+i0,datapiv.ny-i0
             do i=1+i0,datapiv.nx-i0
                if (datapiv.w(i,j,is)==0.d0) then

                   do ic=1,datapiv.ncomponent
                      um = 0.d0
                      nval = 0
                      do ii=-i0,i0,1
                         do jj=-i0,i0,1
                            um = um + datapiv.u(i+ii,j+jj,ic,is)&
                                 *datapiv.w(i+ii,j+jj,is)
                            nval = nval + datapiv.w(i+ii,j+jj,is)
                         end do
                      end do

                      if (nval > 0) um = um / nval

                      datapiv.u(i,j,ic,is) = um
                   end do

                end if
             end do
          end do
       end do

    end select

  end subroutine piv_replace_outlier

  subroutine piv_set_x(datapiv)
    class(PIVdata)                     ::datapiv
    !-----------------------------------------------

    if (datapiv.typeofgrid == "C") then

       !fill in x for cartesian coodinate systeme
       allocate(datapiv.x(datapiv.nx,&
            datapiv.ny,&
            datapiv.ncomponent))
       do i=1,datapiv.nx
          do j=1,datapiv.ny
             datapiv.x(i,j,1) = real(i*datapiv.pixel_step) * datapiv.dx + datapiv.x0
             datapiv.x(i,j,2) = real(j*datapiv.pixel_step) * datapiv.dy + datapiv.y0
             datapiv.x(i,j,3) = datapiv.z_pos
          end do
       end do

    else if (datapiv.typeofgrid == "P") then

       !fill in x for polar coodinate systeme
       allocate(datapiv.x(datapiv.nr,&
            datapiv.ntheta,&
            datapiv.ncomponent))
       do i=1,datapiv.nr
          do j=1,datapiv.ntheta
             datapiv.x(i,j,1) = (i-1) * datapiv.dr * cos((j-1)*datapiv.dtheta)
             datapiv.x(i,j,2) = (i-1) * datapiv.dr * sin((j-1)*datapiv.dtheta)
             datapiv.x(i,j,3) = datapiv.z_pos
          end do
       end do

    else
       write(06,*)"piv_io_read : impossible to define the type of grid"
       STOP
    end if

  end subroutine piv_set_x

  subroutine piv_stats(datapiv)
    class(PIVdata)                     ::datapiv
    !-----------------------------------------------
    allocate(datapiv.ustat(datapiv.nx,datapiv.ny,2*datapiv.ncomponent))

    if (.not. allocated(datapiv.w)) then
       allocate(datapiv.w(datapiv.nx,datapiv.ny,datapiv.nsamples))
       datapiv.w = 1.0
    end if

    do j=1,datapiv.ny
       do i=1,datapiv.nx
          do ic=1,datapiv.ncomponent
             call average(datapiv.u(i,j,ic,:),&
                  datapiv.ustat(i,j,ic),&
                  datapiv.w(i,j,:))
             call rms    (datapiv.u(i,j,ic,:),&
                  datapiv.ustat(i,j,ic+datapiv.ncomponent),&
                  datapiv.w(i,j,:))
          end do
       end do
    end do


  end subroutine piv_stats

  subroutine piv_io_read(datapiv,filename)
    class(PIVdata)                     ::datapiv
    character(len=*)                   ::filename
    logical                            ::file_exists
    !----------------------------------------------

    !check existence of binary data file
    inquire(file=trim(filename), EXIST=file_exists)
    if (file_exists) then

       open(unit=110,file=trim(filename),form='unformatted',&
            action='read', access='stream', status='old')

       if (datapiv.typeofgrid == "P") then
          read(110)datapiv.nr, datapiv.ntheta,&
               datapiv.ncomponent,datapiv.nsamples,&
               datapiv.dr, datapiv.dtheta,&
               datapiv.fs

          allocate(datapiv.u(datapiv.nr,&
               datapiv.ntheta,&
               datapiv.ncomponent,datapiv.nsamples))
          read(110)datapiv.u

       else if (datapiv.typeofgrid == "C") then
          !read datapiv info header
          read(110)datapiv.nx, datapiv.ny,&
               datapiv.ncomponent,datapiv.nsamples,&
               datapiv.dx, datapiv.dy,&
               datapiv.x0, datapiv.y0,&
               datapiv.pixel_step,datapiv.fs

          !read velocity samples
          allocate(datapiv.u(datapiv.nx,&
               datapiv.ny,&
               datapiv.ncomponent,datapiv.nsamples))
          read(110)datapiv.u

       else
          write(06,*)"piv_io_read : impossible to define the type of grid"
          STOP
       end if

       close(110)
    else
       write(06,*) trim(filename)," doesn't exist..."
       STOP
    end if

  end subroutine piv_io_read

  subroutine piv_io_write(datapiv,filename)
    class(PIVdata)                     ::datapiv
    character(len=*)                   ::filename
    logical                            ::file_exists
    !----------------------------------------------

    open(unit=110,file=trim(filename),form='unformatted',&
         action='write', access='stream', status='unknown')

    !write datapiv info header
    if (datapiv.typeofgrid == "P") then
       write(110)datapiv.nr, datapiv.ntheta,&
            datapiv.ncomponent,datapiv.nsamples,&
            datapiv.dr, datapiv.dtheta,&
            datapiv.fs

    else if (datapiv.typeofgrid == "C") then
       write(110)datapiv.nx, datapiv.ny,&
            datapiv.ncomponent,datapiv.nsamples,&
            datapiv.dx, datapiv.dy,&
            datapiv.x0, datapiv.y0,&
            datapiv.pixel_step,datapiv.fs
    else
       write(06,*)"piv_io_write : impossible to define the type of grid"
       STOP
    end if

    write(110)datapiv.u

    close(110)

  end subroutine piv_io_write

  subroutine piv_destroy(datapiv)
    class(PIVdata)                     ::datapiv
    !-------------------------------------------

    if (allocated(datapiv.u)) deallocate(datapiv.u)
    if (allocated(datapiv.ustat)) deallocate(datapiv.ustat)
    if (allocated(datapiv.x)) deallocate(datapiv.x)
    if (allocated(datapiv.w)) deallocate(datapiv.w)

  end subroutine piv_destroy


end module lib_piv_data
