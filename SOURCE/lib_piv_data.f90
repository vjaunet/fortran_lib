module lib_piv_data
  implicit none

  !=================Specification=============================
  !*
  !*
  !*
  !*  PIV data container and useful routines
  !*
  !*  author : Vincent Jaunet
  !*  License: GPL v3.0
  !*
  !*  -11/2016 :
  !*        v. jaunet : added support for netcdf data format
  !*
  !===========================================================

  integer, private :: i,j,ic,is

  type statistics
     !statistics containers
     real, dimension(:,:,:),   allocatable ::u_mean
     real, dimension(:,:,:),   allocatable ::u_rms
     real, dimension(:,:,:),   allocatable ::u_skew
     real, dimension(:,:,:),   allocatable ::u_flat
     real, dimension(:,:,:),   allocatable ::xmoments
  end type statistics

  type PIVdata
     !data indices and scaling
     character ::typeofgrid="C"
     integer   ::nx,ny,nsamples,ncomponent
     real      ::dx=1.d0,x0=0.d0,dy=1.d0,y0=0.d0 !for scaling
     integer   ::pixel_step=1
     real      ::z_pos=0.d0
     real      ::fs=1

     !comments
     character(len=500)  ::comments = ""

     !stagnation conditions
     integer   ::ncgen
     real ,dimension(:) ,allocatable :: cgen

     !raw data containers
     real, dimension(:,:,:,:), allocatable ::u !velocity container
     real, dimension(:,:,:),   allocatable ::x !positions
     real, dimension(:,:,:),   allocatable ::w !Flag and weigth

     !statistics
     type(statistics) ::stat

   contains
     procedure :: read_bin  => piv_io_read
     procedure :: write_bin => piv_io_write
     procedure :: read_data  => piv_io_read
     procedure :: write_data => piv_io_write
     procedure :: set_x     => piv_set_x
     procedure :: replace_outlier => piv_replace_outlier
     procedure :: detect_outlier => piv_detect_outlier
     procedure :: cal_stats => piv_stats
     procedure :: cal_avg => piv_average
     procedure :: cal_rms => piv_rms
     procedure :: get_fluctuations => piv_fluctuations
     procedure :: destroy   => piv_destroy
     procedure :: create    => piv_create
     procedure :: print_info=> piv_info

  end type PIVdata

contains

  subroutine piv_set_x(datapiv)
    class(PIVdata)                     ::datapiv
    !-----------------------------------------------

    if (datapiv%typeofgrid == "C") then

       !fill in x for cartesian coodinate systeme
       allocate(datapiv%x(datapiv%nx,&
            datapiv%ny,3))
       do i=0,datapiv%nx-1
          do j=0,datapiv%ny-1
             datapiv%x(i+1,j+1,1) = real(i*datapiv%pixel_step) * datapiv%dx + datapiv%x0
             datapiv%x(i+1,j+1,2) = real(j*datapiv%pixel_step) * datapiv%dy + datapiv%y0
             datapiv%x(i+1,j+1,3) = datapiv%z_pos
          end do
       end do

    else if (datapiv%typeofgrid == "P") then

       !fill in x for polar coodinate systeme
       allocate(datapiv%x(datapiv%nx,&
            datapiv%ny,3))
       do i=1,datapiv%nx
          do j=1,datapiv%ny
             datapiv%x(i,j,1) = (i-1) * datapiv%dx * cos((j-1)*datapiv%dy)
             datapiv%x(i,j,2) = (i-1) * datapiv%dx * sin((j-1)*datapiv%dy)
             datapiv%x(i,j,3) = datapiv%z_pos
          end do
       end do

    else
       write(06,*)"piv_set_x : impossible to define the type of grid"
       STOP
    end if

  end subroutine piv_set_x

  subroutine piv_stats(datapiv,Nstats)
    use lib_stat
    class(PIVdata)                     ::datapiv
    integer          ,optional         ::Nstats
    integer                            ::n1,n2
    integer                            ::Nmax
    !-----------------------------------------------
    n1 = datapiv%nx
    n2 = datapiv%ny

    if (PRESENT(Nstats)) then
       Nmax = Nstats
    else
       Nmax = datapiv%nsamples
    end if

    allocate(datapiv%stat%u_mean(n1,n2,datapiv%ncomponent))
    allocate(datapiv%stat%u_rms(n1,n2,datapiv%ncomponent))
    allocate(datapiv%stat%u_skew(n1,n2,datapiv%ncomponent))
    allocate(datapiv%stat%u_flat(n1,n2,datapiv%ncomponent))

    if (.not. allocated(datapiv%w)) then
       allocate(datapiv%w(n1,n2,datapiv%nsamples))
       datapiv%w = 1.d0
    end if

    do ic=1,datapiv%ncomponent
       do j=1,n2
          do i=1,n1

             call average(datapiv%u(i,j,ic,1:Nmax),&
                  datapiv%stat%u_mean(i,j,ic),&
                  datapiv%w(i,j,1:Nmax))

             call rms    (datapiv%u(i,j,ic,1:Nmax),&
                  datapiv%stat%u_rms(i,j,ic),&
                  datapiv%w(i,j,1:Nmax))

             call skewness (datapiv%u(i,j,ic,1:Nmax),&
                  datapiv%stat%u_skew(i,j,ic),&
                  datapiv%w(i,j,1:Nmax))

             call flatness (datapiv%u(i,j,ic,1:Nmax),&
                  datapiv%stat%u_flat(i,j,ic),&
                  datapiv%w(i,j,1:Nmax))

          end do
       end do
    end do

    if (datapiv%ncomponent == 2) then
       allocate(datapiv%stat%xmoments(n1,n2,1))
       do j=1,n2
          do i=1,n1
             call xmoment(datapiv%u(i,j,1,1:Nmax),&
                  datapiv%u(i,j,2,1:Nmax),&
                  datapiv%stat%xmoments(i,j,1),&
                  datapiv%w(i,j,1:Nmax))
          end do
       end do
    else if (datapiv%ncomponent == 3) then
       allocate(datapiv%stat%xmoments(n1,n2,3))
       do j=1,n2
          do i=1,n1
             call xmoment(datapiv%u(i,j,1,1:Nmax),&
                  datapiv%u(i,j,2,:),&
                  datapiv%stat%xmoments(i,j,1),&
                  datapiv%w(i,j,1:Nmax))

             call xmoment(datapiv%u(i,j,1,1:Nmax),&
                  datapiv%u(i,j,3,:),&
                  datapiv%stat%xmoments(i,j,2),&
                  datapiv%w(i,j,1:Nmax))

             call xmoment(datapiv%u(i,j,2,1:Nmax),&
                  datapiv%u(i,j,3,:),&
                  datapiv%stat%xmoments(i,j,3),&
                  datapiv%w(i,j,1:Nmax))
          end do
       end do
    end if

  end subroutine piv_stats

  subroutine piv_average(datapiv,Nmax)
    use lib_stat
    class(PIVdata)                     ::datapiv
    integer                            ::n1,n2,nn
    integer         ,optional          ::Nmax
    !-----------------------------------------------


    if (present(Nmax)) then
       nn = Nmax
    else
       nn = datapiv%nsamples
    end if

    n1 = datapiv%nx
    n2 = datapiv%ny

    if (.not. allocated(datapiv%w)) then
       allocate(datapiv%w(n1,n2,datapiv%nsamples))
       datapiv%w = 1.d0
    end if

    if (.not.allocated(datapiv%stat%u_mean)) then
       allocate(datapiv%stat%u_mean(n1,n2,datapiv%ncomponent))
    end if

    do ic=1,datapiv%ncomponent
       do j=1,n2
          do i=1,n1

             call average(datapiv%u(i,j,ic,1:nn),&
                  datapiv%stat%u_mean(i,j,ic),&
                  datapiv%w(i,j,1:nn))

          end do
       end do
    end do

  end subroutine piv_average

  subroutine piv_rms(datapiv,Nmax)
    use lib_stat
    class(PIVdata)                     ::datapiv
    integer                            ::n1,n2,nn
    integer         ,optional          ::Nmax
    !-----------------------------------------------

    if (present(Nmax)) then
       nn = Nmax
    else
       nn = datapiv%nsamples
    end if

    if (.not. allocated(datapiv%w)) then
       allocate(datapiv%w(n1,n2,datapiv%nsamples))
       datapiv%w = 1.d0
    end if

    if (.not.allocated(datapiv%stat%u_rms)) then
       n1 = datapiv%nx
       n2 = datapiv%ny

       allocate(datapiv%stat%u_rms(n1,n2,datapiv%ncomponent))
       do ic=1,datapiv%ncomponent
          do j=1,n2
             do i=1,n1

                call rms    (datapiv%u(i,j,ic,1:nn),&
                     datapiv%stat%u_rms(i,j,ic),&
                     datapiv%w(i,j,1:nn))

             end do
          end do
       end do
    end if

  end subroutine piv_rms

  subroutine piv_fluctuations(datapiv)
    use lib_stat
    class(PIVdata)                     ::datapiv
    integer                            ::n1,n2
    !-----------------------------------------------

    if (.not.allocated(datapiv%stat%u_mean)) then
       call datapiv%cal_avg()
    end if

    do is=1,datapiv%nsamples
       datapiv%u(:,:,:,is) = datapiv%u(:,:,:,is) - datapiv%stat%u_mean(:,:,:)
    end do

  end subroutine piv_fluctuations



  !*****************************************************************
  !
  !      INPUT - OUTPUTS
  !
  !
  !*****************************************************************


  subroutine piv_io_read(datapiv,filename,fmt)
    class(PIVdata)                     ::datapiv
    character(len=*)                   ::filename
    character(len=*)   ,optional       ::fmt
    !----------------------------------------------

    if (present(fmt)) then
       select case (fmt)
       case ('bin')
          call piv_io_read_bin(datapiv,filename)
       ! case ('netCDF')
       !    call piv_io_read_netcdf(datapiv,filename)
       case default
          write(06,*) 'piv_io_read : unknown data format "',fmt,'"'
          STOP
       end select
    else
       call piv_io_read_bin(datapiv,filename)
    end if

  end subroutine piv_io_read

  ! subroutine piv_io_read_netcdf(datapiv,filename)
  !   use lib_netcdf
  !   class(PIVdata)                     ::datapiv
  !   type(netcdf_data)                  ::ncdfpiv
  !   character(len=*)                   ::filename
  !   logical                            ::file_exists

  !   integer                            ::n1,n2,ic
  !   !----------------------------------------------

  !   !check existence of binary data file
  !   inquire(file=trim(filename), EXIST=file_exists)
  !   if (file_exists) then
  !      call ncdfpiv%read_data(filename)
  !   else
  !      write(06,*) trim(filename)," doesn't exist..."
  !      STOP
  !   end if

  !   !switch from netCDF to pivdata format
  !   !------------------------------------
  !   if (ncdfpiv%ndim == 1) &
  !        STOP 'piv_io_read_netcdf : wrong number of dimensions in the piv file'

  !   !beware netcdf data are column major mode (c++ style)
  !   !that is why x and y are inverted
  !   datapiv%nx = ncdfpiv%dimensions(2)%len
  !   datapiv%ny = ncdfpiv%dimensions(1)%len
  !   if (ncdfpiv%ndim == 3) then
  !      datapiv%nsamples   = ncdfpiv%dimensions(3)%len
  !      datapiv%fs         = 1.d0/(ncdfpiv%var(3)%data(2)-ncdfpiv%var(3)%data(1))
  !   else
  !      datapiv%nsamples   = 1
  !   end if
  !   datapiv%ncomponent = ncdfpiv%nvar-ncdfpiv%ndim
  !   datapiv%ncgen = 0

  !   !allocate memory
  !   call piv_create(datapiv)

  !   !set coordinates
  !   datapiv%x0 = ncdfpiv%var(2)%data(1)
  !   datapiv%y0 = ncdfpiv%var(1)%data(1)
  !   datapiv%dx = ncdfpiv%var(2)%data(2)-ncdfpiv%var(1)%data(1)
  !   datapiv%dy = ncdfpiv%var(1)%data(2)-ncdfpiv%var(2)%data(1)

  !   call piv_set_x(datapiv)

  !   !fill variable memory
  !   do ic=1,datapiv%ncomponent
  !      datapiv%u(:,:,ic,:)=reshape(ncdfpiv%var(ic+ncdfpiv%ndim)%data,&
  !           (/datapiv%nx,datapiv%ny,datapiv%nsamples/))
  !   end do

  !   !free netcdf memory use
  !   call ncdfpiv%destroy()

  ! end subroutine piv_io_read_netcdf

  subroutine piv_io_read_bin(datapiv,filename)
    class(PIVdata)                     ::datapiv
    character(len=*)                   ::filename
    logical                            ::file_exists

    integer                            ::n1,n2
    !----------------------------------------------

    !check existence of binary data f<ile
    inquire(file=trim(filename), EXIST=file_exists)
    if (file_exists) then

       open(unit=110,file=trim(filename),form='unformatted',&
            action='read', access='stream', status='old')

       read(110)datapiv%typeofgrid

       if (datapiv%typeofgrid == "C" .or. datapiv%typeofgrid == "P") then

          !read datapiv info header
          read(110)datapiv%nx, datapiv%ny,&
               datapiv%ncomponent,datapiv%nsamples,&
               datapiv%dx, datapiv%dy,&
               datapiv%x0, datapiv%y0,&
               datapiv%pixel_step,&
               datapiv%fs,&
               datapiv%z_pos

          n1 = datapiv%nx
          n2 = datapiv%ny

       else
          write(06,*)"piv_io_read : impossible to define the type of grid"
          write(06,*)"              be sure to define the correct data format"
          STOP
       end if

       !read statgnation conditions if some
       read(110)datapiv%ncgen
       if (datapiv%ncgen > 0) then
          allocate(datapiv%cgen(datapiv%ncgen))
          read(110)datapiv%cgen
       end if

       !read 500 comment characters
       read(110)datapiv%comments

       !read velocity samples<
       allocate(datapiv%u(n1,n2,datapiv%ncomponent,&
            datapiv%nsamples))
       read(110)datapiv%u

       close(110)
    else
       write(06,*) trim(filename)," doesn't exist..."
       STOP
    end if

  end subroutine piv_io_read_bin

  subroutine piv_io_write(datapiv,filename,fmt)
    class(PIVdata)                     ::datapiv
    character(len=*)                   ::filename
    character(len=*)   ,optional       ::fmt
    !----------------------------------------------

    if (present(fmt)) then
       select case (fmt)
       case ('bin')
          call piv_io_write_bin(datapiv,filename)
       ! case ('netCDF')
       !       call piv_io_write_netcdf(datapiv,filename)
       case default
          STOP 'piv_io_write : unknom data format (bin or netCDF)'
       end select
    else
       call piv_io_write_bin(datapiv,filename)
    end if

  end subroutine piv_io_write

  subroutine piv_io_write_bin(datapiv,filename)
    class(PIVdata)                     ::datapiv
    character(len=*)                   ::filename
    logical                            ::file_exists
    !----------------------------------------------

    open(unit=110,file=trim(filename),form='unformatted',&
         action='write', access='stream', status='unknown')

    write(110)datapiv%typeofgrid

    !write datapiv info header
    write(110)datapiv%nx, datapiv%ny,&
         datapiv%ncomponent,datapiv%nsamples,&
         datapiv%dx, datapiv%dy,&
         datapiv%x0, datapiv%y0,&
         datapiv%pixel_step,datapiv%fs,&
         datapiv%z_pos

    !write statgnation conditions if some
    write(110)datapiv%ncgen
    if (datapiv%ncgen > 0) then
       write(110)datapiv%cgen
    end if

    !write comments
    write(110)datapiv%comments

    !write data
    write(110)datapiv%u

    close(110)

  end subroutine piv_io_write_bin

  ! subroutine piv_io_write_netcdf(datapiv,filename)
  !   use lib_netcdf
  !   class(PIVdata)                     ::datapiv
  !   type(netcdf_data)                  ::ncdfpiv
  !   character(len=*)                   ::filename

  !   character(len=4)   ,dimension(:), allocatable  ::varname
  !   !----------------------------------------------

  !   !switch from pivdata to netCDF format
  !   !------------------------------------
  !   !beware netcdf data are column major mode (c++ style)
  !   !that is why x and y are inverted
  !   allocate(varname(datapiv%nsamples))
  !   do ic=1,datapiv%ncomponent
  !      write(varname(ic),'(a1,i1)')'u',ic
  !   end do

  !   if (datapiv%nsamples /= 1) then

  !      call ncdfpiv%create(&
  !           (/datapiv%ny,datapiv%nx, datapiv%nsamples/),&
  !           (/"y","x","t"/),&
  !           datapiv%ncomponent+3,&
  !           (/"yvar","xvar","tvar",(varname(ic),ic=1,datapiv%ncomponent)/)&
  !           )
  !   else
  !      call ncdfpiv%create(&
  !           (/datapiv%ny,datapiv%nx/),&
  !           (/"y","x"/),&
  !           datapiv%ncomponent+2,&
  !           (/"yvar","xvar",(varname(ic),ic=1,datapiv%ncomponent)/)&
  !           )
  !   end if

  !   do i=1,datapiv%nx
  !      ncdfpiv%var(2)%data(i) = datapiv%x0 + (i-1)*datapiv%dx
  !   end do
  !   do i=1,datapiv%ny
  !      ncdfpiv%var(1)%data(i) = datapiv%y0 + (i-1)*datapiv%dy
  !   end do
  !   if (datapiv%nsamples /= 1) then
  !      do i=1,datapiv%nsamples
  !         ncdfpiv%var(3)%data(i) = (i-1)/datapiv%fs
  !      end do
  !   end if
  !   do ic=1,datapiv%ncomponent
  !      ncdfpiv%var(ncdfpiv%ndim+ic)%data = reshape(datapiv%u(:,:,ic,:),&
  !           (/size(ncdfpiv%var(ncdfpiv%ndim+ic)%data)/))
  !   end do

  !   !write the data
  !   call ncdfpiv%write_data(filename)

  !   !clear some memory
  !   call ncdfpiv%destroy()

  ! end subroutine piv_io_write_netcdf

  !***********************************************************
  !
  !
  !       Constructor / Destuctor Routines
  !
  !
  !***********************************************************


  subroutine piv_destroy(datapiv)
    class(PIVdata)                     ::datapiv
    !-------------------------------------------

    if (allocated(datapiv%u)) deallocate(datapiv%u)
    if (allocated(datapiv%stat%u_mean)) deallocate(datapiv%stat%u_mean)
    if (allocated(datapiv%stat%u_rms))  deallocate(datapiv%stat%u_rms)
    if (allocated(datapiv%stat%u_skew)) deallocate(datapiv%stat%u_skew)
    if (allocated(datapiv%stat%u_flat)) deallocate(datapiv%stat%u_flat)
    if (allocated(datapiv%x)) deallocate(datapiv%x)
    if (allocated(datapiv%w)) deallocate(datapiv%w)
    if (allocated(datapiv%cgen)) deallocate(datapiv%cgen)

  end subroutine piv_destroy

  subroutine piv_create(datapiv)
    class(PIVdata)                     ::datapiv
    !-------------------------------------------

    if (.not. allocated(datapiv%u)) then
       if ((datapiv%nx /= 0) .and. (datapiv%ny /= 0)&
            &.and. (datapiv%ncomponent/=0) .and. (datapiv%nsamples /=0)) then
          allocate(datapiv%u(datapiv%nx,datapiv%ny,datapiv%ncomponent,datapiv%nsamples))
       else
          STOP "datapiv : can't allocate memory, a table size equals 0"
       end if
    end if

    if (.not. allocated(datapiv%cgen)) then
       if ((datapiv%ncgen /= 0))  allocate(datapiv%cgen(datapiv%ncgen))
    end if

  end subroutine piv_create

  subroutine piv_info(datapiv)
    class(PIVdata)                     ::datapiv
    !-------------------------------------------

    write(06,*)"File infos :"
    if (datapiv%typeofgrid=="C") then
       write(06,*)"Cartesian grid"
       write(06,'(a,i3,a,i3,a,i3,a,i5)')"  - nx = ",datapiv%nx,", ny = ",datapiv%ny,&
            ", ncompoments = ",datapiv%ncomponent,", nsamples = ", datapiv%nsamples
       write(06,'(a,f6.2,a,f6.2,a,f6.2,a,f6.2)')"  - x0 = ",datapiv%x0,", y0 = ",datapiv%y0,&
            ", dx = ",datapiv%dx,", dy = ", datapiv%dy

    end if

    if (datapiv%typeofgrid=="P") then
       write(06,*)"Polar grid"
       write(06,'(a,i3,a,i3,a,i3,a,i5)')"  - nr = ",datapiv%nx,", ntheta = ",datapiv%ny,&
            ", ncompoments = ",datapiv%ncomponent,", nsamples = ", datapiv%nsamples
       write(06,'(a,f6.2,a,f6.2,a,f6.2,a,f6.2)')"  - r0 = ",datapiv%x0,", thetha0 = ",datapiv%y0,&
            ", dr = ",datapiv%dx,", dtheta = ", datapiv%dy

    end if

    write(06,*)" - Sampling frequency :",datapiv%fs
    write(06,*)" - z position :",datapiv%z_pos
    write(06,*)" - Comments :",trim(datapiv%comments)

    if (datapiv%ncgen > 0) then
       write(06,'(a,10(f10.3,2x))')"  - Stagnation Conditions :",(datapiv%cgen(i),i=1,datapiv%ncgen)
    else
       write(06,'(a)') " No Stagnation condition stored"
    end if

  end subroutine piv_info



  !***********************************************************************
  !
  !          PIV filter routines
  !
  !***********************************************************************
  subroutine piv_detect_outlier(datapiv,method,wsize,nsigma)
    use lib_pod
    use lib_stat
    use omp_lib
    class(PIVdata)                            ::datapiv
    character(len=*)                          ::method
    integer    ,optional                      ::wsize
    real       ,optional                      ::Nsigma

    real                                      ::um,utest
    integer                                   ::ii,jj,i0
    integer                                   ::nx,ny,ns,nc
    real                                      ::nval
    integer                                   ::w_def
    real                                      ::Nsigma_def


    integer(kind=8)                           ::id,if,jd,jf,ivec,nvec
    real(kind=8) ,dimension(:) ,allocatable   ::neighbor,neighflag
    !-------------------------------------------------------------

    !check method is correct
    if ( method /= "UOD" .and. &
         method /= "AVG" .and. &
         method /= "MED" .and. &
         method /= "SIG") then

       STOP 'detect_outlier : unknown method must be &
            &"UOD", "AVG", "MED" or "SIG"'

    end if

    !put default interrogation area size
    if (.not. PRESENT(wsize)) then
       w_def = 3;
    else
       w_def = wsize
    end if

    !put default Nsigma parameter
    if (.not. PRESENT(Nsigma)) then
       Nsigma_def = 5.d0;
    else
       Nsigma_def = Nsigma
    end if

    ns = datapiv%nsamples
    nc = datapiv%ncomponent
    ny = datapiv%ny
    nx = datapiv%nx

    if (.not. allocated(datapiv%w)) then
       allocate(datapiv%w(nx,ny,ns))
       datapiv%w = 1.d0
    end if

    !--------------------------------------
    ! Nsigma detection

    if (method == "SIG") then

       !computing stats for nsigma outlier detection
       call piv_average(datapiv)
       call piv_rms(datapiv)

       !loop through the samples
       do is=1,ns
          do ic=1,nc
             do j=1,ny
                do i=1,nx
                   if (datapiv%w(i,j,is) /= 0.d0) then
                      if (abs(datapiv%u(i,j,ic,is)-datapiv%stat%u_mean(i,j,ic))&
                           .gt. Nsigma_def*datapiv%stat%u_rms(i,j,ic)) then
                         datapiv%w(i,j,is) = 0.d0
                      end if
                   end if
                end do
             end do
          end do
       end do

       deallocate(datapiv%stat%u_mean)
       deallocate(datapiv%stat%u_rms)

       return
    end if

    !------------------------------
    ! spatial detection methods

    !loop through all the samples
    call omp_set_num_threads(25)
!!!! !$OMP PARALLEL DO PRIVATE(ic,j,i,ii,jj,id,if,jd,jf,utest,nvec,neighbor,neighflag) SHARED(datapiv,method) SCHEDULE(Dynamic)
    do is=1,ns
       do ic=1,nc
          do j=1,ny
             do i=1,nx

                !test only if vector is not already an outlier
                if (datapiv%w(i,j,is) /= 0.d0) then

                   !computing the Interrogation Area size
                   if (i <= w_def/2) then
                      id = i-w_def/2 ; if = w_def/2
                   else if (i > nx-w_def/2) then
                      id = -w_def/2 ; if = nx-i
                   else
                      id = -w_def/2 ; if = w_def/2
                   end if

                   if (j <= w_def/2) then
                      jd = j-w_def/2 ; jf = w_def/2
                   else if (j > ny-w_def/2) then
                      jd = -w_def/2 ; jf = ny - j
                   else
                      jd = -w_def/2 ; jf = w_def/2
                   end if

                   !storing the Interrogation area
                   !excluding the center sample
                   utest = datapiv%u(i,j,ic,is)
                   nvec = (if-id+1)*(jf-jd+1)-1
                   allocate (neighbor(nvec))
                   allocate (neighflag(nvec))
                   ivec = 1
                   do ii=id,if
                      do jj =jd,jf
                         if (ii/=0 .or. jj/=0 ) then
                            neighbor(ivec)  = dble(datapiv%u(i+ii,j+jj,ic,is))
                            neighflag(ivec) = dble(datapiv%w(i+ii,j+jj,is))
                            ivec = ivec+1
                         end if
                      end do
                   end do

                   select case(method)

                   case ("UOD")

                      !computing the UOD on the subsample
                      call UOD_filter(utest,neighbor,neighflag,nvec)

                   case("AVG")

                      !computing the AVERAGE filter on the subsample
                      call average_filter(utest,neighbor,neighflag,nvec)

                   case("MED")

                      !computing the MEDIAN filter on the subsample
                      call median_filter(utest,neighbor,neighflag,nvec)

                   end select

                   !storing the flag :
                   !if flag = 0 an outlier has been detected
                   datapiv%w(i,j,ic) = neighflag(1)

                   deallocate(neighbor)
                   deallocate(neighflag)

                end if !w/=0

             end do
          end do
       end do
    end do
!!!!   !$OMP END PARALLEL DO

    return

  end subroutine piv_detect_outlier

  subroutine piv_replace_outlier(datapiv,method,wsize)
    use lib_pod
    class(PIVdata)                            ::datapiv
    character(len=*)                          ::method
    integer    ,optional                      ::wsize

    real                                      ::um
    integer                                   ::ii,jj,i0
    integer                                   ::nx,ny,ns,nc
    real                                      ::nval
    integer                                   ::w_def


    integer(kind=8)                           ::id,if,jd,jf,ivec,nvec
    real(kind=8) ,dimension(:) ,allocatable   ::neighbor,neighflag
    !-------------------------------------------------------------

    !check method is correct
    if (method /= "POD" .and. &
         method /= "UOD" .and. &
         method /= "AVG" .and. &
         method /= "MED" ) then

       STOP 'replace_outlier : unknown method must be &
            &"POD", "UOD", "AVG", "MED"'

    end if

    !put default interrogation area size
    if (.not. PRESENT(wsize)) then
       w_def = 3;
    else
       w_def = wsize
    end if

    if (method == "POD") then
       call gappypod(datapiv%u,datapiv%w,int(10),real(1e-8))
       return
    end if

    ns = datapiv%nsamples
    nc = datapiv%ncomponent
    ny = datapiv%ny
    nx = datapiv%nx

    !loop through all the samples
    do is=1,ns
       do ic=1,nc
          do j=1,ny
             do i=1,nx

                !Replacing only if necessary
                if (datapiv%w(i,j,is) == 0.d0) then

                   !computing the Interrogation Area size
                   if (i <= w_def/2) then
                      id = i-w_def/2 ; if = w_def/2
                   else if (i > nx-w_def/2) then
                      id = -w_def/2 ; if = nx-i
                   else
                      id = -w_def/2 ; if = w_def/2
                   end if

                   if (j <= w_def/2) then
                      jd = j-w_def/2 ; jf = w_def/2
                   else if (j > ny-w_def/2) then
                      jd = -w_def/2 ; jf = ny - j
                   else
                      jd = -w_def/2 ; jf = w_def/2
                   end if

                   !storing the Interrogation area
                   !excluding the center sample
                   nvec = (if-id+1)*(jf-jd+1)-1
                   allocate (neighbor(nvec))
                   allocate (neighflag(nvec))
                   ivec = 1
                   do ii=id,if
                      do jj =jd,jf
                         if (ii/=0 .or. jj/=0 ) then
                            neighbor(ivec)  = dble(datapiv%u(i+ii,j+jj,ic,is))
                            neighflag(ivec) = dble(datapiv%w(i+ii,j+jj,is))
                            ivec = ivec+1
                         end if
                      end do
                   end do

                   select case(method)

                   case ("UOD")

                      !computing the UOD on the subsample
                      call UOD_filter(datapiv%u(i,j,ic,is),neighbor,neighflag,nvec)

                   case("AVG")

                      !computing the AVERAGE filter on the subsample
                      call average_filter(datapiv%u(i,j,ic,is),neighbor,neighflag,nvec)

                   case("MED")

                      !computing the MEDIAN filter on the subsample
                      call median_filter(datapiv%u(i,j,ic,is),neighbor,neighflag,nvec)

                   end select

                   !storing the flag :
                   !if flag = 0 a replacement has been done
                   datapiv%w(i,j,ic) = neighflag(1)

                   deallocate(neighbor)
                   deallocate(neighflag)

                end if !end w=0

             end do
          end do
       end do
    end do

    return

  end subroutine piv_replace_outlier

  SUBROUTINE UOD_filter(input,neighbor,flag,Nneigh)
    use qsort_c
    implicit none

    !input variables
    real(kind=4)                                       ::input
    integer(kind=8)                                    ::Nneigh
    real(kind=8)        ,dimension(NNEIGH)             ::neighbor,flag

    !Intern variables
    real(kind=8)                                       ::Um,rm,r0
    real(kind=8)        ,dimension(:) ,allocatable     ::tosort
    real(kind=8)        ,dimension(:) ,allocatable     ::resid
    real(kind=8)                                       ::epsilon = 0.1
    real(kind=8)                                       ::threshold = 2.d0

    !loops
    integer(kind=8)                                    ::i,pos,nvalid
    integer(kind=8)                                    ::iv

    !------------------------------------------------------

    nvalid = 0
    do i=1,Nneigh
       if (flag(i) == 1.d0) then
          nvalid = nvalid + 1
       end if
    end do

    if (nvalid < 3) return;

    allocate(tosort(Nvalid))
    iv = 1
    do i=1,Nneigh
       if (flag(i) == 1.d0) then
          tosort(iv) = neighbor(i)
          iv = iv +1
       end if
    end do

    !calculating median of Neighbors Um
    call QsortC(tosort)
    Um = tosort(nvalid/2)

    !calculating median of residuals rm
    tosort =  abs(tosort - Um)
    call QsortC(tosort)
    rm = tosort(Nvalid/2)

    !replace sample if necessary
    r0 = abs(input - Um)/(rm + epsilon)
    if (r0 > threshold) then
       input = Um
       flag(1) = 0.d0
    end if


  END SUBROUTINE UOD_filter

  SUBROUTINE median_filter(input,neighbor,flag,Nneigh)
    use qsort_c
    implicit none

    !input variables
    real(kind=4)                                       ::input
    integer(kind=8)                                    ::Nneigh
    real(kind=8)        ,dimension(NNEIGH)             ::neighbor,flag

    !Intern variables
    real(kind=8)                                       ::median
    real(kind=8)        ,dimension(:) ,allocatable     ::tosort

    !loops
    integer(kind=8)                                    ::i,pos,nvalid

    !------------------------------------------------------

    nvalid = 0
    do i=1,Nneigh
       if (flag(i) == 1.d0) then
          nvalid = nvalid + 1
       end if
    end do

    if (nvalid < 3) return

    allocate(tosort(Nvalid))
    nvalid = 1
    do i=1,Nneigh
       if (flag(i) == 1.d0) then
          tosort(nvalid) = neighbor(i)
          nvalid = nvalid + 1
       end if
    end do

    call QsortC(tosort)
    median = neighbor(nvalid/2)

    if (abs(input) > 2.d0*median) then
       input = median
       flag(1) = 0.d0
    end if


  END SUBROUTINE median_filter

  SUBROUTINE average_filter(input,neighbor,flag,Nneigh)
    implicit none

    !input variables
    real(kind=4)                                       ::input
    integer(kind=8)                                    ::Nneigh
    real(kind=8)        ,dimension(NNEIGH)             ::neighbor,flag

    !Intern variables
    real(kind=8)                                       ::mean
    real(kind=8)                                       ::rms
    real(kind=8)                                       ::sumf

    !loops
    integer(kind=8)                                    ::i,pos

    !------------------------------------------------------

    sumf = sum(flag(:))
    if (sumf < 3. ) then
       return
    end if

    mean = 0.d0
    do i=1,Nneigh
       mean = mean + neighbor(i)*flag(i)
    end do
    mean = mean/sumf

    rms = 0.d0
    do i=1,Nneigh
       rms = rms + (((neighbor(i)-mean)**2)*flag(i))
    end do
    rms = sqrt(rms/sumf)

    !if (abs(input-mean) > 2.d0*rms) then
    input = mean
    flag(1) = 0.d0
    !end if

  END SUBROUTINE average_filter

end module lib_piv_data
