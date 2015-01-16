module lib_piv_data
  implicit none

  !=================Specification=============================
  !
  !
  !
  !
  !            PIV data container and useful routines
  !
  !
  !
  !
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
     integer   ::ntheta,nr
     real      ::dx=1,x0=0,dy=1,y0=0 !for scaling
     real      ::dr=0,dtheta=0
     integer   ::pixel_step=1
     real      ::z_pos=0
     integer   ::fs=1

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
     procedure :: set_x     => piv_set_x
     procedure :: replace_outlier => piv_replace_outlier
     procedure :: cal_stats => piv_stats
     procedure :: get_fluctuations => piv_fluctuations
     procedure :: destroy   => piv_destroy
     procedure :: print_info=> piv_info

  end type PIVdata

contains

  subroutine piv_set_x(datapiv)
    class(PIVdata)                     ::datapiv
    !-----------------------------------------------

    if (datapiv.typeofgrid == "C") then

       !fill in x for cartesian coodinate systeme
       allocate(datapiv.x(datapiv.nx,&
            datapiv.ny,3))
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
            datapiv.ntheta,3))
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
    use lib_stat
    class(PIVdata)                     ::datapiv
    integer                            ::n1,n2
    !-----------------------------------------------
    if(datapiv.typeofgrid == 'C') then
       n1 = datapiv.nx
       n2 = datapiv.ny
    else
       n1 = datapiv.nr
       n2 = datapiv.ntheta
    end if

    allocate(datapiv.stat.u_mean(n1,n2,datapiv.ncomponent))
    allocate(datapiv.stat.u_rms(n1,n2,datapiv.ncomponent))
    allocate(datapiv.stat.u_skew(n1,n2,datapiv.ncomponent))
    allocate(datapiv.stat.u_flat(n1,n2,datapiv.ncomponent))

    if (.not. allocated(datapiv.w)) then
       allocate(datapiv.w(n1,n2,datapiv.nsamples))
       datapiv.w = 1.d0
    end if

    do ic=1,datapiv.ncomponent
       do j=1,n2
          do i=1,n1

             call average(datapiv.u(i,j,ic,:),&
                  datapiv.stat.u_mean(i,j,ic),&
                  datapiv.w(i,j,:))

             call rms    (datapiv.u(i,j,ic,:),&
                  datapiv.stat.u_rms(i,j,ic),&
                  datapiv.w(i,j,:))

             call skewness (datapiv.u(i,j,ic,:),&
                  datapiv.stat.u_skew(i,j,ic),&
                  datapiv.w(i,j,:))

             call flatness (datapiv.u(i,j,ic,:),&
                  datapiv.stat.u_flat(i,j,ic),&
                  datapiv.w(i,j,:))

          end do
       end do
    end do

    if (datapiv.ncomponent == 2) then
       allocate(datapiv.stat.xmoments(n1,n2,1))
       do j=1,n2
          do i=1,n1
             call xmoment(datapiv.u(i,j,1,:),&
                  datapiv.u(i,j,2,:),&
                  datapiv.stat.xmoments(i,j,1),&
                  datapiv.w(i,j,:))
          end do
       end do
    else if (datapiv.ncomponent == 3) then
       allocate(datapiv.stat.xmoments(n1,n2,3))
       do j=1,n2
          do i=1,n1
             call xmoment(datapiv.u(i,j,1,:),&
                  datapiv.u(i,j,2,:),&
                  datapiv.stat.xmoments(i,j,1),&
                  datapiv.w(i,j,:))

             call xmoment(datapiv.u(i,j,1,:),&
                  datapiv.u(i,j,3,:),&
                  datapiv.stat.xmoments(i,j,2),&
                  datapiv.w(i,j,:))

             call xmoment(datapiv.u(i,j,2,:),&
                  datapiv.u(i,j,3,:),&
                  datapiv.stat.xmoments(i,j,3),&
                  datapiv.w(i,j,:))
          end do
       end do
    end if

  end subroutine piv_stats

  subroutine piv_fluctuations(datapiv)
    use lib_stat
    class(PIVdata)                     ::datapiv
    integer                            ::n1,n2
    !-----------------------------------------------

    if (.not.allocated(datapiv.stat.u_mean)) then
       if(datapiv.typeofgrid == 'C') then
          n1 = datapiv.nx
          n2 = datapiv.ny
       else
          n1 = datapiv.nr
          n2 = datapiv.ntheta
       end if

       allocate(datapiv.stat.u_mean(n1,n2,datapiv.ncomponent))
       do ic=1,datapiv.ncomponent
          do j=1,n2
             do i=1,n1

                call average(datapiv.u(i,j,ic,:),&
                     datapiv.stat.u_mean(i,j,ic))

             end do
          end do
       end do
    end if

    do is=1,datapiv.nsamples
       datapiv.u(:,:,:,is) = datapiv.u(:,:,:,is) - datapiv.stat.u_mean(:,:,:)
    end do

  end subroutine piv_fluctuations

  subroutine piv_io_read(datapiv,filename)
    class(PIVdata)                     ::datapiv
    character(len=*)                   ::filename
    logical                            ::file_exists

    integer                            ::n1,n2
    !----------------------------------------------

    !check existence of binary data file
    inquire(file=trim(filename), EXIST=file_exists)
    if (file_exists) then

       open(unit=110,file=trim(filename),form='unformatted',&
            action='read', access='stream', status='old')

       read(110)datapiv.typeofgrid

       if (datapiv.typeofgrid == "P") then
          read(110)datapiv.nr, datapiv.ntheta,&
               datapiv.ncomponent,datapiv.nsamples,&
               datapiv.dr, datapiv.dtheta,&
               datapiv.fs

          n1 = datapiv.nr
          n2 = datapiv.ntheta

          !read data containers
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


          n1 = datapiv.nx
          n2 = datapiv.ny

       else
          write(06,*)"piv_io_read : impossible to define the type of grid"
          STOP
       end if

       !read statgnation conditions if some
       read(110)datapiv.ncgen
       if (datapiv.ncgen > 0) then
          allocate(datapiv.cgen(datapiv.ncgen))
          read(110)datapiv.cgen
       end if

       !read 500 comment characters
       read(110)datapiv.comments

       !read velocity samples
       allocate(datapiv.u(n1,n2,datapiv.ncomponent,&
            datapiv.nsamples))
       read(110)datapiv.u

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

    write(110)datapiv.typeofgrid

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

    !write statgnation conditions if some
    write(110)datapiv.ncgen
    if (datapiv.ncgen > 0) then
       write(110)datapiv.cgen
    end if

    !write comments
    write(110)datapiv.comments


    write(110)datapiv.u

    close(110)

  end subroutine piv_io_write

  subroutine piv_destroy(datapiv)
    class(PIVdata)                     ::datapiv
    !-------------------------------------------

    if (allocated(datapiv.u)) deallocate(datapiv.u)
    if (allocated(datapiv.stat.u_mean)) deallocate(datapiv.stat.u_mean)
    if (allocated(datapiv.stat.u_rms))  deallocate(datapiv.stat.u_rms)
    if (allocated(datapiv.stat.u_skew)) deallocate(datapiv.stat.u_skew)
    if (allocated(datapiv.stat.u_flat)) deallocate(datapiv.stat.u_flat)
    if (allocated(datapiv.x)) deallocate(datapiv.x)
    if (allocated(datapiv.w)) deallocate(datapiv.w)
    if (allocated(datapiv.cgen)) deallocate(datapiv.cgen)

  end subroutine piv_destroy

  subroutine piv_info(datapiv)
    class(PIVdata)                     ::datapiv
    !-------------------------------------------

    write(06,*)"File infos :"
    if (datapiv.typeofgrid=="C") then
       write(06,*)"Cartesian grid"
       write(06,'(a,i3,a,i3,a,i3,a,i5)')"  - nx = ",datapiv.nx,", ny = ",datapiv.ny,&
            ", ncompoments = ",datapiv.ncomponent,", nsamples = ", datapiv.nsamples
    end if

    if (datapiv.typeofgrid=="P") then
       write(06,*)"Polar grid"
       write(06,'(a,i3,a,i3,a,i3,a,i5)')"  - nr = ",datapiv.nr,", ntheta = ",datapiv.ntheta,&
            ", ncompoments = ",datapiv.ncomponent,", nsamples = ", datapiv.nsamples
    end if

    write(06,*)" - Sampling frequency :",datapiv.fs
    write(06,*)" - Comments :",trim(datapiv.comments)

    if (datapiv.ncgen > 0) then
       write(06,'(a,10(f10.3,2x))')"  - Stagnation Conditions :",(datapiv.cgen(i),i=1,datapiv.ncgen)
    end if

  end subroutine piv_info



  !***********************************************************************
  !
  !          PIV filter routines
  !
  !***********************************************************************
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
       call gappypod(datapiv.u,datapiv.w,int(10),real(1e-8))
       return
    end if

    ns = datapiv.nsamples
    nc = datapiv.ncomponent
    if (datapiv.typeofgrid == 'C') then
       ny = datapiv.ny
       nx = datapiv.nx
    else if (datapiv.typeofgrid == 'P') then
       ny = datapiv.nr
       nx = datapiv.ntheta
    end if

    !loop through all the samples
    do is=1,ns
       do ic=1,nc
          do j=1,ny
             do i=1,nx

                !Replacing only if necessary
                if (datapiv.w(i,j,is) == 0.d0) then

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
                            neighbor(ivec)  = dble(datapiv.u(i+ii,j+jj,ic,is))
                            neighflag(ivec) = dble(datapiv.w(i+ii,j+jj,is))
                            ivec = ivec+1
                         end if
                      end do
                   end do

                   select case(method)

                   case ("UOD")

                      !computing the UOD on the subsample
                      call UOD_filter(datapiv.u(i,j,ic,is),neighbor,neighflag,nvec)

                   case("AVG")

                      !computing the AVERAGE filter on the subsample
                      call average_filter(datapiv.u(i,j,ic,is),neighbor,neighflag,nvec)

                   case("MED")

                      !computing the MEDIAN filter on the subsample
                      call median_filter(datapiv.u(i,j,ic,is),neighbor,neighflag,nvec)

                   end select

                   !storing the flag :
                   !if flag = 0 a replacement has been done
                   datapiv.w(i,j,ic) = neighflag(1)

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

    if (nvalid < 5) return;

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
    mean = mean/sum(flag(:))

    rms = 0.d0
    do i=1,Nneigh
       rms = rms + (((neighbor(i)-mean)**2)*flag(i))
    end do
    rms = sqrt(rms/sum(flag(:)))

    if (abs(input-mean) > 2.d0*rms) then
       input = mean
       flag(1) = 0.d0
    end if

  END SUBROUTINE average_filter

end module lib_piv_data
