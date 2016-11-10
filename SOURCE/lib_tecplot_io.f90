module lib_tecplot_IO
  implicit none
  !=================Specification=============================
  !*
  !*       Module for Reading and Writing Tecplot data
  !*
  !*  author : Vincent Jaunet
  !*  date   : 24-03-2014
  !*  License: GPL v3.0
  !*
  !* The module contains several routines and an
  !* extended type "filetype".
  !*
  !* - filetype contains all the necessary information
  !* about the file to be read or written : file name, access,
  !* exitence, ID, variable names and titles
  !*
  !* Usage :
  !===========================================================

  integer(kind=8), private ::i,j,k,ic,in,ix

  type tec_filetype
     character(len=250)             ::filename
     character(len=1)               ::RWaccess='R'
     integer                        ::exist=-1
     integer                        ::fid

     integer                        ::nx=0,ny=0,nz=0
     integer                        ::nvar=0,nc=0

     character(len=250)             ::varnames = ""
     character(len=250)             ::title = ""

   contains
     procedure :: openfile  => tec_openfile
     procedure :: closefile => tec_closefile

     generic, public    :: read_header  => tec_read_header_io,&
          tec_read_header_BSA
     procedure, private :: tec_read_header_io
     procedure, private :: tec_read_header_BSA

     generic, public    :: write_header => tec_write_header
     procedure, private :: tec_write_header

     procedure, public    :: read_zone => tec_get_zone

     generic, public    :: read_data => read_ascii_1d, &
          read_ascii_2d, read_ascii_3d
     procedure, private :: read_ascii_1d
     procedure, private :: read_ascii_2d
     procedure, private :: read_ascii_3d

     generic, public    :: write_data => write_ascii_1d,&
          write_ascii_2d, write_ascii_3d,&
          write_ascii_1d_x,&
          write_ascii_2d_x, write_ascii_3d_x
     procedure, private :: write_ascii_1d
     procedure, private :: write_ascii_1d_x
     procedure, private :: write_ascii_2d
     procedure, private :: write_ascii_2d_x
     procedure, private :: write_ascii_3d
     procedure, private :: write_ascii_3d_x



  end type tec_filetype


contains

  !=============Openning and closing files===================
  !*
  !*
  !==========================================================

  subroutine tec_openfile(filespec)
    implicit none
    class(tec_filetype)                      ::filespec
    !--------------------------------------------

    if (filespec%exist == 1) then
       if (filespec%RWaccess == 'R') then
          open (unit=filespec%fid,&
               file = trim(filespec%filename),&
               status = 'old',&
               action='read')
       elseif (filespec%RWaccess == 'W') then
          open (unit=filespec%fid,&
               file = trim(filespec%filename),&
               status = 'old',&
               action='write')
       end if
    elseif (filespec%exist == 0) then
       if (filespec%RWaccess == 'R') then
          open (unit=filespec%fid,&
               file = trim(filespec%filename),&
               status = 'new',&
               action='read')
       elseif (filespec%RWaccess == 'W') then
          open (unit=filespec%fid,&
               file = trim(filespec%filename),&
               status = 'new',&
               action='write')
       end if
    else
       open (unit=filespec%fid,&
            file = trim(filespec%filename))
    end if

  end subroutine tec_openfile

  subroutine tec_closefile(filespec)
    implicit none
    class(tec_filetype)                      ::filespec
    !--------------------------------------------

    close(filespec%fid)

  end subroutine tec_closefile

  !==== Header handling =================================
  !*
  !*
  !==========================================================

  subroutine tec_write_header(filespec)
    implicit none
    class(tec_filetype)                  ::filespec
    !-----------------------------------------------

    write(filespec%fid,'(a,a,a)')'Title = "',trim(filespec%title),'"'
    write(filespec%fid,'(a,a)')'Variables = ',trim(filespec%varnames)

  end subroutine tec_write_header

  subroutine tec_read_header_io(filespec)
    implicit none
    class(tec_filetype)                 ::filespec
    character(len=9)                    ::trash1
    character(len=12)                   ::trash2
    !--------------------------------------------

    !read title and variable names
    read(filespec%fid,'(a9,a)')trash1,filespec%title
    read(filespec%fid,'(a12,a)')trash2,filespec%varnames

    !Remove remaining " from title
    filespec%title = filespec%title(1:len_trim(filespec%title)-1)

  end subroutine tec_read_header_io

  subroutine tec_read_header_BSA(filespec,t0)
    implicit none
    class(tec_filetype)                 ::filespec
    integer                             ::nl
    character(len=9)                    ::trash1
    character(len=12)                   ::trash2
    character(len=18)                   ::trash18
    character(len=1)                    ::trash4
    integer(kind=8)                     ::h,m,s
    real(kind=8)                        ::t0
    !--------------------------------------------

    read(filespec%fid,'(a9,a)')trash1,filespec%title

    !bsa header
    read(filespec%fid,*)
    read(filespec%fid,'(a18,i2,a1,i2,a1,i2)')trash18,h,trash4,m,trash4,s
    t0 = (h*60+m)*60+s
    read(filespec%fid,*)
    read(filespec%fid,*)

    read(filespec%fid,'(a12,a)')trash2,filespec%varnames
    read(filespec%fid,*)

    !Remove remaining " from title
    filespec%title = filespec%title(1:len_trim(filespec%title)-1)

  end subroutine tec_read_header_BSA

  !====== Zone handling =====================================
  !*
  !*
  !==========================================================


  subroutine tec_get_zone(filespec,zone_title)
    class(tec_filetype)                      ::filespec

    character(len=250)                  ::zone
    integer                             ::i,ii
    integer                             ::ios
    character(len=10)                   ::trash1

    integer(kind=4)                     ::nx,ny,nc
    character(len=*), optional        ::zone_title
    !-----------------------------------------------

    !get the number of variables
    ios = 0
    i = 0
    do ii=1,len_trim(filespec%varnames)
       if (filespec%varnames(ii:ii) == '"') i=i+1
    end do
    filespec%nvar = i/2

    !read zone line and extract the data
    read(filespec%fid,'(a)')zone

    i=1
    ii=0
    do while (i<len_trim(zone)-1)

       if (zone(i:i) == '"') then
          !get zone title
          ii = 0
          do while (zone(i+ii+1:i+ii+1) /= '"' &
               .and. i+ii+1<len_trim(zone)+1)
             ii = ii+1
          end do
          if (present(zone_title)) then
             zone_title = adjustl(zone(i+1:i+ii))
          end if
          ii=ii+1

       else if (zone(i:i+1) == 'I=') then
          !get nx

          ii = 1
          do while (zone(i+ii+1:i+ii+1) /= ',' &
               .and. i+ii+1<len_trim(zone)+1)
             ii = ii+1
          end do
          trash1 = adjustl(zone(i+2:i+ii))
          read(trash1,*)filespec%nx

       else if (zone(i:i+1) == "J=") then
          !get ny

          ii = 1
          do while (zone(i+ii+1:i+ii+1) /= ',' &
               .and. i+ii+1<len_trim(zone)+1)
             ii = ii+1
          end do
          trash1 = adjustl(zone(i+2:i+ii))
          read(trash1,*)filespec%ny

       else if (zone(i:i+1) == "K=") then
          !get nz

          ii = 1
          do while (zone(i+ii+1:i+ii+1) /= ',' &
               .and. i+ii+1<len_trim(zone)+1)
             ii = ii+1
          end do
          trash1 = adjustl(zone(i+2:i+ii))
          read(trash1,*)filespec%nz

       end if

       i = i+ii+1
       ii=0
    end do

  end subroutine tec_get_zone



  !==== Read 1D data ========================================
  !*
  !*
  !==========================================================

  subroutine read_ascii_1d(filespec,tab)
    implicit none
    class(tec_filetype)                          ::filespec
    integer(kind=8)                              ::nn,nc
    class(*) ,dimension(:,:)                     ::tab
    integer(kind=4) ,dimension(:,:) ,allocatable ::tabi4
    integer(kind=8) ,dimension(:,:) ,allocatable ::tabi8
    real(kind=4) ,dimension(:,:) ,allocatable    ::tabr4
    real(kind=8) ,dimension(:,:) ,allocatable    ::tabr8
    !----------------------------------------------------

    nn = size(tab,1)
    nc = size(tab,2)

    if ((nn /= filespec%nx) .or. nc /= filespec%nvar) then
       STOP 'read_ascii_1d : wrong table size'
    end if

    select type (tab)
    type is (integer)
       allocate(tabi4(nn,nc))
       do in=1,nn
          read(filespec%fid,*)(tabi4(in,ic),ic=1,nc)
       end do
       tab = tabi4
       deallocate(tabi4)

    type is (integer(kind=8))
       allocate(tabi8(nn,nc))
       do in=1,nn
          read(filespec%fid,*)(tabi8(in,ic),ic=1,nc)
       end do
       tab = tabi8
       deallocate(tabi8)

    type is (real)
       allocate(tabr4(nn,nc))
       do in=1,nn
          read(filespec%fid,*)(tabr4(in,ic),ic=1,nc)
       end do
       tab = tabr4
       deallocate(tabr4)

    type is (real(kind=8))
       allocate(tabr8(nn,nc))
       do in=1,nn
          read(filespec%fid,*)(tabr8(in,ic),ic=1,nc)
       end do
       tab = tabr8
       deallocate(tabr8)

    end select


  end subroutine read_ascii_1d


  !==== Read 2D data ========================================
  !*
  !*
  !==========================================================

  subroutine read_ascii_2d(filespec,tab)
    implicit none
    class(tec_filetype)                             ::filespec
    integer(kind=8)                                 ::ni,nj,nc
    class(*), intent(out), dimension(:,:,:)         ::tab
    integer(kind=4) ,dimension(:,:,:) ,allocatable ::tabi4
    integer(kind=8) ,dimension(:,:,:) ,allocatable ::tabi8
    real(kind=4) ,dimension(:,:,:) ,allocatable    ::tabr4
    real(kind=8) ,dimension(:,:,:) ,allocatable    ::tabr8

    !----------------------------------------------------

    ni = size(tab,1)
    nj = size(tab,2)
    nc = size(tab,3)

    if (ni /= filespec%nx &
         .or. nj /= filespec%ny&
         .or. nc /= filespec%nvar) then
       STOP 'read_ascii_1d : wrong table size'
    end if

    select type (tab)
    type is (integer)
       allocate(tabi4(ni,nj,nc))
       do j=1,nj
          do i=1,ni
             read(filespec%fid,*)(tabi4(i,j,ic),ic=1,nc)
          end do
       end do
       tab = tabi4
       deallocate(tabi4)

    type is (integer(kind=8))
       allocate(tabi8(ni,nj,nc))
       do j=1,nj
          do i=1,ni
             read(filespec%fid,*)(tabi8(i,j,ic),ic=1,nc)
          end do
       end do
       tab = tabi8
       deallocate(tabi8)

    type is (real)
       allocate(tabr4(ni,nj,nc))
       do j=1,nj
          do i=1,ni
             read(filespec%fid,*)(tabr4(i,j,ic),ic=1,nc)
          end do
       end do
       tab = tabr4
       deallocate(tabr4)

    type is (real(kind=8))
       allocate(tabr8(ni,nj,nc))
       do j=1,nj
          do i=1,ni
             read(filespec%fid,*)(tabr8(i,j,ic),ic=1,nc)
          end do
       end do
       tab = tabr8
       deallocate(tabr8)

    end select


  end subroutine read_ascii_2d

  !==== read 3D data ====================================
  !*
  !*
  !==========================================================

  subroutine read_ascii_3d(filespec,tab)
    implicit none
    class(tec_filetype)                              ::filespec
    integer(kind=8)                                  ::ni,nj,nk,nc
    class(*), intent(out), dimension(:,:,:,:)        ::tab
    integer(kind=4) ,dimension(:,:,:,:) ,allocatable ::tabi4
    integer(kind=8) ,dimension(:,:,:,:) ,allocatable ::tabi8
    real(kind=4) ,dimension(:,:,:,:) ,allocatable    ::tabr4
    real(kind=8) ,dimension(:,:,:,:) ,allocatable    ::tabr8

    !----------------------------------------------------

    ni = size(tab,1)
    nj = size(tab,2)
    nk = size(tab,3)
    nc = size(tab,4)

    if (ni /= filespec%nx &
         .or. nj /= filespec%ny&
         .or. nj /= filespec%nz&
         .or. nc /= filespec%nvar) then
       STOP 'read_ascii_3d : wrong table size'
    end if

    select type (tab)
    type is (integer)
       allocate(tabi4(ni,nj,nk,nc))
       do k=1,nk
          do j=1,nj
             do i=1,ni
                read(filespec%fid,*)(tabi4(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
       tab = tabi4
       deallocate(tabi4)

    type is (integer(kind=8))
       allocate(tabi8(ni,nj,nk,nc))
       do k=1,nk
          do j=1,nj
             do i=1,ni
                read(filespec%fid,*)(tabi8(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
       tab = tabi8
       deallocate(tabi8)

    type is (real)
       allocate(tabr4(ni,nj,nk,nc))
       do k=1,nk
          do j=1,nj
             do i=1,ni
                read(filespec%fid,*)(tabr4(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
       tab = tabr4
       deallocate(tabr4)

    type is (real(kind=8))
       allocate(tabr8(ni,nj,nk,nc))
       do k=1,nk
          do j=1,nj
             do i=1,ni
                read(filespec%fid,*)(tabr8(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
       tab = tabr8
       deallocate(tabr8)

    end select


  end subroutine read_ascii_3d



  !==== Write 1D data ====================================
  !*
  !*
  !==========================================================
  subroutine write_ascii_1d_x(filespec,x,tab,zonetitle)
    implicit none
    class(tec_filetype)                                      ::filespec
    integer(kind=8)                                          ::ni,nc,nx
    class(*), intent(in), dimension(:,:)                     ::tab
    class(*), intent(in), dimension(:,:)                     ::x
    character(len=*)  ,optional                              ::zonetitle

    integer,   dimension(:,:), allocatable                   ::bufi4
    integer(kind=8), dimension(:,:), allocatable                   ::bufi8
    real,      dimension(:,:), allocatable                   ::bufr4
    real(kind=8),    dimension(:,:), allocatable                   ::bufr8
    !----------------------------------------------------------------------

    ni=size(tab,1)
    nc=size(tab,2)
    nx=size(x,2)

    select type (x)
    type is (integer)
       allocate(bufi4(ni,nc+nx))
       bufi4(:,1:nx) = x(:,:)
    type is (integer(kind=8))
       allocate(bufi8(ni,nc+nx))
       bufi8(:,1:nx) = x(:,:)
    type is (real)
       allocate(bufr4(ni,nc+nx))
       bufr4(:,1:nx) = x(:,:)
    type is (real(kind=8))
       allocate(bufr8(ni,nc+nx))
       bufr8(:,1:nx) = x(:,:)
    end select

    select type (tab)
    type is (integer)
       bufi4(:,nx+1:nx+nc) = tab(:,:)
       if (present(zonetitle)) then
          call write_ascii_1d(filespec,bufi4,zonetitle)
       else
          call write_ascii_1d(filespec,bufi4)
       end if
       deallocate(bufi4)

    type is (integer(kind=8))
       bufi8(:,nx+1:nx+nc) = tab(:,:)
       if (present(zonetitle)) then
          call write_ascii_1d(filespec,bufi8,zonetitle)
       else
          call write_ascii_1d(filespec,bufi8)
       end if
       deallocate(bufi8)

    type is (real(kind=4))
       bufr4(:,nx+1:nx+nc) = tab(:,:)
       if (present(zonetitle)) then
          call write_ascii_1d(filespec,bufr4,zonetitle)
       else
          call write_ascii_1d(filespec,bufr4)
       end if
       deallocate(bufr8)

    type is (real(kind=8))
       bufr8(:,nx+1:nx+nc) = tab(:,:)
       if (present(zonetitle)) then
          call write_ascii_1d(filespec,bufr8,zonetitle)
       else
          call write_ascii_1d(filespec,bufr8)
       end if
       deallocate(bufr8)

    end select

  end subroutine write_ascii_1d_x

  subroutine write_ascii_1d(filespec,tab,zonetitle)
    implicit none
    class(tec_filetype)                             ::filespec
    integer(kind=8)                                 ::nn,nc
    class(*), intent(in), dimension(:,:)            ::tab
    character(len=*) ,optional                      ::zonetitle

    integer, dimension(:,:), allocatable            ::tabi4
    integer(kind=8), dimension(:,:), allocatable    ::tabi8
    real   , dimension(:,:), allocatable            ::tabr4
    real(kind=8) , dimension(:,:), allocatable      ::tabr8

    character(len=150)                              ::fmt
    !-----------------------------------------------------------

    nn = size(tab,1)
    nc = size(tab,2)

    write(fmt,'(i0)')nc+1
    fmt=fmt//"(e15.8,2x)"
1666 format (a,i0)
1665 format (a,a,a,i0)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',nn
    else
       write(filespec%fid,1666)'Zone T="", I=',nn
    end if

    select type (tab)
    type is (integer)
       allocate(tabi4(nn,nc),source=tab)
       do in=1,nn
          write(filespec%fid,trim(fmt))(tabi4(in,ic),ic=1,nc)
       end do
       deallocate(tabi4)

    type is (integer(kind=8))
       allocate(tabi8(nn,nc),source=tab)
       do in=1,nn
          write(filespec%fid,trim(fmt))(tabi8(in,ic),ic=1,nc)
       end do
       deallocate(tabi8)

    type is (real)
       allocate(tabr4(nn,nc),source=tab)
       do in=1,nn
          write(filespec%fid,trim(fmt))(tabr4(in,ic),ic=1,nc)
       end do
       deallocate(tabr4)

    type is (real(kind=8))
       allocate(tabr8(nn,nc),source=tab)
       do in=1,nn
          write(filespec%fid,trim(fmt))(tabr8(in,ic),ic=1,nc)
       end do
       deallocate(tabr8)
    end select

  end subroutine write_ascii_1d


  !==== Write 2D data =============================================
  !*
  !*
  !================================================================
  subroutine write_ascii_2d_x(filespec,x,tab,zonetitle)
    implicit none
    class(tec_filetype)                                      ::filespec
    integer(kind=8)                                          ::ni,nj,nk,nc,nx
    class(*), intent(in), dimension(:,:,:)                   ::tab
    class(*), intent(in), dimension(:,:,:)                   ::x
    character(len=*)  ,optional                              ::zonetitle

    integer,   dimension(:,:,:), allocatable                 ::bufi4
    integer(kind=8), dimension(:,:,:), allocatable                 ::bufi8
    real,      dimension(:,:,:), allocatable                 ::bufr4
    real(kind=8),    dimension(:,:,:), allocatable                 ::bufr8
    !----------------------------------------------------------------------

    ni=size(tab,1)
    nj=size(tab,2)
    nc=size(tab,3)
    nx=size(x,3)

    select type (x)
    type is (integer)
       allocate(bufi4(ni,nj,nc+nx))
       bufi4(:,:,1:nx) = x(:,:,:)
    type is (integer(kind=8))
       allocate(bufi8(ni,nj,nc+nx))
       bufi8(:,:,1:nx) = x(:,:,:)
    type is (real)
       allocate(bufr4(ni,nj,nc+nx))
       bufr4(:,:,1:nx) = x(:,:,:)
    type is (real(kind=8))
       allocate(bufr8(ni,nj,nc+nx))
       bufr8(:,:,1:nx) = x(:,:,:)
    end select

    select type (tab)
    type is (integer)
       bufi4(:,:,nx+1:nx+nc) = tab(:,:,:)
       if (present(zonetitle)) then
          call write_ascii_2d(filespec,bufi4,zonetitle)
       else
          call write_ascii_2d(filespec,bufi4)
       end if
       deallocate(bufi4)

    type is (integer(kind=8))
       bufi8(:,:,nx+1:nx+nc) = tab(:,:,:)
       if (present(zonetitle)) then
          call write_ascii_2d(filespec,bufi8,zonetitle)
       else
          call write_ascii_2d(filespec,bufi8)
       end if
       deallocate(bufi8)

    type is (real(kind=4))
       bufr4(:,:,nx+1:nx+nc) = tab(:,:,:)
       if (present(zonetitle)) then
          call write_ascii_2d(filespec,bufr4,zonetitle)
       else
          call write_ascii_2d(filespec,bufr4)
       end if
       deallocate(bufr4)

    type is (real(kind=8))
       bufr8(:,:,nx+1:nx+nc) = tab(:,:,:)
       if (present(zonetitle)) then
          call write_ascii_2d(filespec,bufr8,zonetitle)
       else
          call write_ascii_2d(filespec,bufr8)
       end if
       deallocate(bufr8)

    end select

  end subroutine write_ascii_2d_x


  subroutine write_ascii_2d(filespec,tab,zonetitle)
    implicit none
    class(tec_filetype)                               ::filespec
    class(*), intent(in), dimension(:,:,:)            ::tab
    character(len=*)                  ,  optional     ::zonetitle
    integer(kind=8)                                   ::ni,nj,nc

    integer, dimension(:,:,:), allocatable            ::tabi4
    integer(kind=8), dimension(:,:,:), allocatable          ::tabi8
    real   , dimension(:,:,:), allocatable            ::tabr4
    real(kind=8) , dimension(:,:,:), allocatable            ::tabr8
    character(len=150)                                ::fmt
    !-------------------------------------------------------------

    ni=size(tab,1)
    nj=size(tab,2)
    nc=size(tab,3)

!1664 format (<2*nc>(e15.8,2x))
    write(fmt,'(a,i0,a)')"(",2*nc,"(e15.8,2x))"

1665 format (a,a,a,i0,a,i0)
1666 format (a,a,i0,a,i0)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',ni,', J=',nj
    else
       write(filespec%fid,1666)'Zone T="',&
            '", I=',ni,', J=',nj
    end if

    select type (tab)
    type is (integer)
       allocate(tabi4(ni,nj,nc),source=tab)
       do j=1,nj
          do i=1,ni
             write(filespec%fid,trim(fmt))(tabi4(i,j,ic),ic=1,nc)
          end do
       end do
       deallocate(tabi4)

    type is (integer(kind=8))
       allocate(tabi8(ni,nj,nc),source=tab)
       do j=1,nj
          do i=1,ni
             write(filespec%fid,trim(fmt))(tabi8(i,j,ic),ic=1,nc)
          end do
       end do
       deallocate(tabi8)

    type is (real)
       allocate(tabr4(ni,nj,nc),source=tab)
       do j=1,nj
          do i=1,ni
             write(filespec%fid,trim(fmt))(tabr4(i,j,ic),ic=1,nc)
          end do
       end do
       deallocate(tabr4)

    type is (real(kind=8))
       allocate(tabr8(ni,nj,nc),source=tab)
       do j=1,nj
          do i=1,ni
             write(filespec%fid,trim(fmt))(tabr8(i,j,ic),ic=1,nc)
          end do
       end do
       deallocate(tabr8)
    end select


  end subroutine write_ascii_2d


  !==== Write 3D data =======================================
  !*
  !*
  !==========================================================
  subroutine write_ascii_3d_x(filespec,x,tab,zonetitle)
    implicit none
    class(tec_filetype)                                      ::filespec
    integer(kind=8)                                          ::ni,nj,nk,nc,nx
    class(*), intent(in), dimension(:,:,:,:)                 ::tab
    class(*), intent(in), dimension(:,:,:,:)                 ::x
    character(len=*)  ,optional                              ::zonetitle

    integer, dimension(:,:,:,:), allocatable                 ::bufi4
    integer(kind=8), dimension(:,:,:,:), allocatable               ::bufi8
    real, dimension(:,:,:,:), allocatable                    ::bufr4
    real(kind=8), dimension(:,:,:,:), allocatable                  ::bufr8
    !----------------------------------------------------------------------

    ni=size(tab,1)
    nj=size(tab,2)
    nk=size(tab,3)
    nc=size(tab,4)
    nx=size(x,4)

    select type (x)
    type is (integer)
       allocate(bufi4(ni,nj,nk,nc+nx))
       bufi4(:,:,:,1:nx) = x(:,:,:,:)
    type is (integer(kind=8))
       allocate(bufi8(ni,nj,nk,nc+nx))
       bufi8(:,:,:,1:nx) = x(:,:,:,:)
    type is (real)
       allocate(bufr4(ni,nj,nk,nc+nx))
       bufr4(:,:,:,1:nx) = x(:,:,:,:)
    type is (real(kind=8))
       allocate(bufr8(ni,nj,nk,nc+nx))
       bufr8(:,:,:,1:nx) = x(:,:,:,:)
    end select

    select type (tab)
    type is (integer)
       bufi4(:,:,:,nx+1:nx+nc) = tab(:,:,:,:)
       if (present(zonetitle)) then
          call write_ascii_3d(filespec,bufi4,zonetitle)
       else
          call write_ascii_3d(filespec,bufi4)
       end if
       deallocate(bufi4)

    type is (integer(kind=8))
       bufi8(:,:,:,nx+1:nx+nc) = tab(:,:,:,:)
       if (present(zonetitle)) then
          call write_ascii_3d(filespec,bufi8,zonetitle)
       else
          call write_ascii_3d(filespec,bufi8)
       end if
       deallocate(bufi8)

    type is (real(kind=4))
       bufr4(:,:,:,nx+1:nx+nc) = tab(:,:,:,:)
       if (present(zonetitle)) then
          call write_ascii_3d(filespec,bufr4,zonetitle)
       else
          call write_ascii_3d(filespec,bufr4)
       end if
       deallocate(bufr4)

    type is (real(kind=8))
       bufr8(:,:,:,nx+1:nx+nc) = tab(:,:,:,:)
       if (present(zonetitle)) then
          call write_ascii_3d(filespec,bufr8,zonetitle)
       else
          call write_ascii_3d(filespec,bufr8)
       end if
       deallocate(bufr8)

    end select

  end subroutine write_ascii_3d_x

  subroutine write_ascii_3d(filespec,tab,zonetitle)
    implicit none
    class(tec_filetype)                                      ::filespec
    integer(kind=8)                                          ::ni,nj,nk,nc
    class(*), intent(in), dimension(:,:,:,:)                 ::tab
    character(len=*)  ,optional                              ::zonetitle


    integer, dimension(:,:,:,:), allocatable                 ::tabi4
    integer(kind=8), dimension(:,:,:,:), allocatable         ::tabi8
    real   , dimension(:,:,:,:), allocatable                 ::tabr4
    real(kind=8) , dimension(:,:,:,:), allocatable           ::tabr8

    character(len=150)                                       ::fmt
    !----------------------------------------------------------------------

    ni=size(tab,1)
    nj=size(tab,2)
    nk=size(tab,3)
    nc=size(tab,4)

    !1664 format (<2*nc>(e15.8,2x))
    write(fmt,'(a,i0,a)')"(",2*nc,"(e15.8,2x))"

1665 format (a,a,a,i0,a,i0,a,i0)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),&
            '", I=',ni,', J=',nj,', K=',nk
    else
       write(filespec%fid,1665)'Zone T="',&
            '", I=',ni,', J=',nj,', K=',nk
    end if

    select type (tab)
    type is (integer)
       allocate(tabi4(ni,nj,nk,nc), source=tab)
       do k=1,nk
          do j=1,nj
             do i=1,ni
                write(filespec%fid,trim(fmt))(tabi4(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
       deallocate(tabi4)

    type is (integer(kind=8))
       allocate(tabi8(ni,nj,nk,nc), source=tab)
       do k=1,nk
          do j=1,nj
             do i=1,ni
                write(filespec%fid,trim(fmt))(tabi8(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
       deallocate(tabi8)

    type is (real)
       allocate(tabr4(ni,nj,nk,nc), source=tab)
       do k=1,nk
          do j=1,nj
             do i=1,ni
                write(filespec%fid,trim(fmt))(tabr4(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
       deallocate(tabr4)

    type is (real(kind=8))
       allocate(tabr8(ni,nj,nk,nc), source=tab)
       do k=1,nk
          do j=1,nj
             do i=1,ni
                write(filespec%fid,trim(fmt))(tabr8(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
       deallocate(tabr8)
    end select

  end subroutine write_ascii_3d


  !--------------------------------------------
  !! end of module
end module lib_tecplot_IO
