module tecplot_IO

  !=================Specification=============================
  !*
  !*       Module for Reading and Writing Tecplot data
  !*
  !*  author : Vincent Jaunet
  !*  date   : 24-03-2014
  !*  License: MIT
  !*
  !* The module contains several routines and an
  !* extended type "filetype".
  !*
  !* - filetype contains all the necessary information
  !* about the file to be read or written : file name, access,
  !* exitence, ID, variable names and titles
  !*
  !* Usage :
  !*
  !* PROGRAM test
  !* use tecplot_IO
  !* implicit none
  !*
  !* type(filetype)   ::ofile
  !*
  !* !Defines the file parameters :
  !* ofile = filetype(filename = 'file.dat',&
  !*      RWaccess = 'W',& ! "W" for writting, "R" for reading, "" for both
  !*      exist = .true.,& ! if the file must exist, .false if not
  !*      fid = 10,&       ! fid number
  !*      varnames = '"x" "y" "U" "V"',& !Varaible names as a string
  !*      title = 'Title')               !File title as string
  !*
  !* !Open the file
  !* call tec_openfile(ofile)
  !*
  !* !Write the Tecplot header:
  !* !Title = title
  !* !Variables = varnames
  !* call tec_write_header(ofile)
  !*
  !* !write a Zone data + Zone header
  !* !Note that the same call can be used for different data type
  !* !Real, Double or integer values are implemented
  !* !Ex : call tec_write_ascii(filetype,& !file descrition
  !*                            x,&        !Optional spatial coordiantes
  !*                            nx,&       !Optional for 2D number spatial coordiantes
  !*                            data,&
  !*                            N_samples,&
  !*                            N_variables,&
  !*                            Zone_title)
  !*
  !* !For 1D data:
  !* call tec_write_ascii(ofile,data,'Zone_Title')
  !* call tec_write_ascii(ofile,x,data,'Zone_Title')
  !*
  !* !For 2D data:
  !* call tec_write_ascii(ofile,data,'Zone_Title')
  !* call tec_write_ascii(ofile,x,data,'Zone_Title')
  !*
  !* !For 3D data:
  !* call tec_write_ascii(ofile,data,'Zone_Title')
  !* call tec_write_ascii(ofile,x,data,'Zone_Title')
  !*
  !* !close file
  !* call tec_closefile(ofile)
  !*
  !*
  !* Reading data :
  !* iifile = filetype(filename = 'res-new-eemd-io.dat',&
  !*      RWaccess = 'R',&
  !*      exist = .true.,&
  !*      fid = 10)
  !*
  !* call tec_openfile(iifile)
  !* call tec_read_header(iifile)
  !* call tec_get_zone(iifile,nt,nc)
  !* call tec_read_ascii(iifile,data,nt,nc)
  !* call tec_close(iifile)
  !*
  !*
  !*
  !* end
  !============================================================

  integer(kind=8), private ::i,j,k,ic,in,ix

  type filetype
     character(len=250)             ::filename
     character(len=1)               ::RWaccess
     integer                        ::exist=-1
     integer                        ::fid

     character(len=250)             ::varnames = ""
     character(len=250)             ::title = ""


  end type filetype

  public  :: tec_openfile,tec_closefile,tec_write_ascii,tec_read_ascii,&
       tec_write_header,tec_read_header,tec_get_zone

  private :: fwrite_ascii_1d, dwrite_ascii_1d, iwrite_ascii_1d,&
       fwrite_ascii_1d_1c, dwrite_ascii_1d_1c, iwrite_ascii_1d_1c,&
       fwrite_ascii_2d,dwrite_ascii_2d,iwrite_ascii_2d,&
       dwrite_ascii_2d_1c,&
       fwrite_ascii_3d,dwrite_ascii_3d,iwrite_ascii_3d,&
       fread_ascii_1d, dread_ascii_1d, iread_ascii_1d,&
       fread_ascii_2d,dread_ascii_2d,iread_ascii_2d,&
       fread_ascii_3d,dread_ascii_3d,iread_ascii_3d

  interface tec_read_ascii
     module procedure fread_ascii_1d, dread_ascii_1d, iread_ascii_1d,&
          fread_ascii_2d,dread_ascii_2d,iread_ascii_2d,&
          fread_ascii_3d,dread_ascii_3d,iread_ascii_3d
  end interface tec_read_ascii

  interface tec_write_ascii
     module procedure fwrite_ascii_1d, dwrite_ascii_1d, iwrite_ascii_1d,&
          fwrite_ascii_1d_1c, dwrite_ascii_1d_1c, iwrite_ascii_1d_1c,&
          fwrite_ascii_2d,dwrite_ascii_2d,iwrite_ascii_2d,&
          dwrite_ascii_2d_1c,&
          fwrite_ascii_3d,dwrite_ascii_3d,iwrite_ascii_3d
  end interface tec_write_ascii


contains

  !=============Openning and closing files===================
  !*
  !*
  !==========================================================

  subroutine tec_openfile(filespec)
    implicit none
    type(filetype)                      ::filespec
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
    type(filetype)                      ::filespec
    !--------------------------------------------

    close(filespec%fid)

  end subroutine tec_closefile

  !==== Header handling =================================
  !*
  !*
  !==========================================================

  subroutine tec_write_header(filespec)
    implicit none
    type(filetype)                      ::filespec
    !--------------------------------------------

    write(filespec%fid,'(a,a,a)')'Title = "',trim(filespec%title),'"'
    write(filespec%fid,'(a,a)')'Variables = ',trim(filespec%varnames)

  end subroutine tec_write_header

  subroutine tec_read_header(filespec)
    implicit none
    type(filetype)                      ::filespec
    character(len=9)                    ::trash1
    character(len=12)                   ::trash2
    !--------------------------------------------

    !read title and variable names
    read(filespec%fid,'(a9,a)')trash1,filespec%title
    read(filespec%fid,'(a12,a)')trash2,filespec%varnames

    !Remove remaining " from title
    filespec%title = filespec%title(1:len_trim(filespec%title)-1)

  end subroutine tec_read_header

  subroutine tec_get_zone(filespec,ni,nc,nj,nk,zone_title)
    type(filetype)                      ::filespec

    character(len=250)                  ::zone
    integer                             ::i,ii
    integer                             ::ios
    character(len=10)                   ::trash1

    integer(kind=8)                     ::ni,nc
    integer(kind=8),    optional        ::nj,nk
    character(len=250), optional        ::zone_title
    !-----------------------------------------------

    !get the number of variables
    ios = 0
    i = 0
    do ii=1,len_trim(filespec%varnames)
       if (filespec%varnames(ii:ii) == '"') i=i+1
    end do
    nc = i/2

    !read zone line and extract the data
    read(filespec%fid,'(a)')zone

    i=1
    do while (i<len_trim(zone)-1)

       if (zone(i:i+1) == '"') then
          !get zone title

          ii = 0
          do while (zone(i+ii+1:i+ii+1) == '"' &
               .or. i+ii+1==len_trim(zone))
             ii = ii+1
          end do
          zone_title = adjustl(zone(i+2:i+ii))

       else if (zone(i:i+1) == "I=") then
          !get nx

          ii = 1
          do while (zone(i+ii+1:i+ii+1) == ',' &
               .or. i+ii+1==len_trim(zone))
             ii = ii+1
          end do
          trash1 = adjustl(zone(i+2:i+ii))
          read(trash1,*)nx

       else if (zone(i:i+1) == "J=") then
          !get ny

          ii = 1
          do while (zone(i+ii+1:i+ii+1) == ',' &
               .or. i+ii+1==len_trim(zone))
             ii = ii+1
          end do
          trash1 = adjustl(zone(i+2:i+ii))
          read(trash1,*)ny

       else if (zone(i:i+1) == "K=") then
          !get nz

          ii = 1
          do while (zone(i+ii+1:i+ii+1) == ',' &
               .or. i+ii+1==len_trim(zone))
             ii = ii+1
          end do
          trash1 = adjustl(zone(i+2:i+ii))
          read(trash1,*)nz
       end if

       i = i+ii+1
    end do

  end subroutine tec_get_zone


  !==== Read 1D data ========================================
  !*
  !*
  !==========================================================

  subroutine fread_ascii_1d(filespec,tab)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::nn,nc
    real, intent(out), dimension(:,:)          ::tab
    !----------------------------------------------------

    nn = size(tab,1)
    nc = size(tab,2)

    do in=1,nn
       read(filespec%fid,*)(tab(in,ic),ic=1,nc)
    end do

  end subroutine fread_ascii_1d

  subroutine dread_ascii_1d(filespec,tab)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::nn,nc
    real*8,     intent(out), dimension(:,:)   ::tab
    !----------------------------------------------------

    nn = size(tab,1)
    nc = size(tab,2)

    do in=1,nn
       read(filespec%fid,*)(tab(in,ic),ic=1,nc)
    end do

  end subroutine dread_ascii_1d

  subroutine iread_ascii_1d(filespec,tab)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::nn,nc
    integer,   intent(out), dimension(:,:)     ::tab
    !----------------------------------------------------

    nn = size(tab,1)
    nc = size(tab,2)

    do in=1,nn
       read(filespec%fid,*)(tab(in,ic),ic=1,nc)
    end do

  end subroutine iread_ascii_1d

  !==== Read 2D data ========================================
  !*
  !*
  !==========================================================

  subroutine fread_ascii_2d(filespec,tab)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::ni,nj,nc
    real, intent(out), dimension(:,:,:)        ::tab
    !----------------------------------------------------

    ni = size(tab,1)
    nj = size(tab,2)
    nc = size(tab,3)

    do j=1,nj
       do i=1,ni
          read(filespec%fid,*)(tab(i,j,ic),ic=1,nc)
       end do
    end do

  end subroutine fread_ascii_2d

  subroutine dread_ascii_2d(filespec,tab)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::ni,nj,nc
    real*8, intent(out), dimension(:,:,:)      ::tab
    !----------------------------------------------------

    ni = size(tab,1)
    nj = size(tab,2)
    nc = size(tab,3)

    do j=1,nj
       do i=1,ni
          read(filespec%fid,*)(tab(i,j,ic),ic=1,nc)
       end do
    end do

  end subroutine dread_ascii_2d

  subroutine iread_ascii_2d(filespec,tab)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::ni,nj,nc
    integer, intent(out), dimension(:,:,:)     ::tab
    !----------------------------------------------------

    ni = size(tab,1)
    nj = size(tab,2)
    nc = size(tab,3)

    do j=1,nj
       do i=1,ni
          read(filespec%fid,*)(tab(i,j,ic),ic=1,nc)
       end do
    end do

  end subroutine iread_ascii_2d

  !==== Read 3D data ========================================
  !*
  !*
  !==========================================================

  subroutine fread_ascii_3d(filespec,tab)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::ni,nj,nk,nc
    real, intent(out), dimension(:,:,:,:)      ::tab
    !----------------------------------------------------

    ni = size(tab,1)
    nj = size(tab,2)
    nk = size(tab,3)
    nc = size(tab,4)

    do k=1,nk
       do j=1,nj
          do i=1,ni
             read(filespec%fid,*)(tab(i,j,k,ic),ic=1,nc)
          end do
       end do
    end do

  end subroutine fread_ascii_3d

  subroutine dread_ascii_3d(filespec,tab)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::ni,nj,nk,nc
    real*8, intent(out), dimension(:,:,:,:)    ::tab
    !----------------------------------------------------

    ni = size(tab,1)
    nj = size(tab,2)
    nk = size(tab,3)
    nc = size(tab,4)

    do k=1,nk
       do j=1,nj
          do i=1,ni
             read(filespec%fid,*)(tab(i,j,k,ic),ic=1,nc)
          end do
       end do
    end do

  end subroutine dread_ascii_3d

  subroutine iread_ascii_3d(filespec,tab)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::ni,nj,nk,nc
    integer, intent(out), dimension(:,:,:,:)   ::tab
    !-----------------------------------------------------------

    ni = size(tab,1)
    nj = size(tab,2)
    nk = size(tab,3)
    nc = size(tab,4)

    do k=1,nk
       do j=1,nj
          do i=1,ni
             read(filespec%fid,*)(tab(i,j,k,ic),ic=1,nc)
          end do
       end do
    end do

  end subroutine iread_ascii_3d



  !==== Write 1D data ====================================
  !*
  !*
  !==========================================================

  subroutine fwrite_ascii_1d(filespec,t,tab,zonetitle)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::nn,nc
    real, intent(in), dimension(:,:)           ::tab
    real, intent(in), dimension(:), optional   ::t
    character(len=*) ,optional                 ::zonetitle
    !-------------------------------------------------

    nn = size(tab,1)
    nc = size(tab,2)

1664 format (<nc+1>(e15.8,2x))
1665 format (a,a,a,i0)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',nn
    else
       write(filespec%fid,1665)'Zone T="", I=',nn
    end if

    if (present(t)) then
       do in=1,nn
          write(filespec%fid,1664)t(in),(tab(in,ic),ic=1,nc)
       end do
    else
       do in=1,nn
          write(filespec%fid,1664)(tab(in,ic),ic=1,nc)
       end do
    end if

  end subroutine fwrite_ascii_1d

  subroutine dwrite_ascii_1d(filespec,t,tab,zonetitle)
    implicit none
    type(filetype)                                       ::filespec
    integer(kind=8)                                      ::nn,nc
    real(kind=8), intent(in), dimension(:,:)             ::tab
    real(kind=8), intent(in), dimension(:),  optional    ::t
    character(len=*)  ,optional                          ::zonetitle
    !----------------------------------------------------------------

    nn = size(tab,1)
    nc = size(tab,2)

1664 format (<nc+1>(e15.8,2x))
1665 format (a,a,a,i0)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',nn
    else
       write(filespec%fid,1665)'Zone T="", I=',nn
    end if

    if (present(t)) then
       do in=1,nn
          write(filespec%fid,1664)t(in),(tab(in,ic),ic=1,nc)
       end do
    else
       do in=1,nn
          write(filespec%fid,1664)(tab(in,ic),ic=1,nc)
       end do
    end if

  end subroutine dwrite_ascii_1d

  subroutine iwrite_ascii_1d(filespec,t,tab,zonetitle)
    implicit none
    type(filetype)                              ::filespec
    integer(kind=8)                             ::nn,nc
    integer, intent(in), dimension(:,:)         ::tab
    integer, intent(in), dimension(:), optional ::t
    character(len=*)   ,optional                ::zonetitle
    !-------------------------------------------------

    nn = size(tab,1)
    nc = size(tab,2)

1664 format (<nc+1>(e15.8,2x))
1665 format (a,a,a,i0)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',nn
    else
       write(filespec%fid,1665)'Zone T="", I=',nn
    end if

    if (present(t)) then
       do in=1,nn
          write(filespec%fid,1664)t(in),(tab(in,ic),ic=1,nc)
       end do
    else
       do in=1,nn
          write(filespec%fid,1664)(tab(in,ic),ic=1,nc)
       end do
    end if

  end subroutine iwrite_ascii_1d

  !=======================================================
  ! for 1d_1c

  subroutine fwrite_ascii_1d_1c(filespec,t,zonetitle,tab)
    implicit none
    type(filetype)                             ::filespec
    integer(kind=8)                            ::nn
    real, intent(in), dimension(:)             ::t
    real, intent(in), dimension(:) ,optional   ::tab
    character(len=*) ,optional                 ::zonetitle
    !-------------------------------------------------

    nn = size(t)

1664 format (2(e15.8,2x))
1665 format (a,a,a,i0)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',nn
    else
       write(filespec%fid,1665)'Zone T="", I=',nn
    end if

    if (present(tab)) then
       do in=1,nn
          write(filespec%fid,1664)t(in),tab(in)
       end do

    else
       do in=1,nn
          write(filespec%fid,1664)t(in)
       end do
    end if

  end subroutine fwrite_ascii_1d_1c

  subroutine dwrite_ascii_1d_1c(filespec,t,tab,zonetitle)
    implicit none
    type(filetype)                                    ::filespec
    real(kind=8), intent(in), dimension(:)            ::t
    real(kind=8), intent(in), dimension(:)            ::tab
    character(len=*)  ,optional                       ::zonetitle
    integer(kind=8)                                   ::nn
    !------------------------------------------------------------

    nn = size(t)

1664 format (2(e15.8,2x))
1665 format (a,a,a,i0)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',nn
    else
       write(filespec%fid,'(a,i0)')'Zone T="", I=',nn
    end if

    !    if (present(tab)) then
    do in=1,nn
       write(filespec%fid,1664)t(in),tab(in)
    end do
    ! else
    !    do in=1,nn
    !       write(filespec%fid,1664)t(in)
    !    end do
    ! end if

  end subroutine dwrite_ascii_1d_1c

  subroutine iwrite_ascii_1d_1c(filespec,t,zonetitle,tab)
    implicit none
    type(filetype)                              ::filespec
    integer(kind=8)                             ::nn
    integer, intent(in), dimension(:)           ::t
    integer, intent(in), dimension(:),optional  ::tab
    character(len=*)   ,optional                ::zonetitle
    !-------------------------------------------------

    nn = size(t)

1664 format (2(e15.8,2x))
1665 format (a,a,a,i0)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',nn
    else
       write(filespec%fid,1665)'Zone T="", I=',nn
    end if

    if (present(tab)) then
       do in=1,nn
          write(filespec%fid,1664)t(in),tab(in)
       end do

    else
       do in=1,nn
          write(filespec%fid,1664)t(in)
       end do
    end if

  end subroutine iwrite_ascii_1d_1c

  !==== Write 2D data =============================================
  !*
  !*
  !================================================================
  subroutine fwrite_ascii_2d(filespec,x,tab,zonetitle)
    implicit none
    type(filetype)                                  ::filespec
    real, intent(in), dimension(:,:,:),  optional   ::x
    real, intent(in), dimension(:,:,:)              ::tab
    character(len=*)                  ,  optional   ::zonetitle
    integer(kind=8)                                 ::ni,nj,nc,nx
    !-------------------------------------------------------------

1664 format (<2*nc>(e15.8,2x))
1665 format (a,a,a,i0,a,i0)

    ni=size(tab,1)
    nj=size(tab,2)
    nc=size(tab,3)
    nx=size(x,3)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',ni,', J=',nj
    else
       write(filespec%fid,1665)'Zone T="',&
            '", I=',ni,', J=',nj
    end if

    if (present(x)) then
       do j=1,nj
          do i=1,ni
             write(filespec%fid,1664)(x(i,j,ix),ix=1,nx),(tab(i,j,ic),ic=1,nc)
          end do
       end do
    else
       do j=1,nj
          do i=1,ni
             write(filespec%fid,1664)(tab(i,j,ic),ic=1,nc)
          end do
       end do
    end if

  end subroutine fwrite_ascii_2d

  subroutine dwrite_ascii_2d(filespec,x,tab,zonetitle)
    implicit none
    type(filetype)                                        ::filespec
    real(kind=8), intent(in), dimension(:,:,:),  optional ::x
    real(kind=8), intent(in), dimension(:,:,:)            ::tab
    character(len=*)                          ,  optional ::zonetitle

    integer(kind=8)                                       ::ni,nj,nc,nx
    !------------------------------------------------------------------

1664 format (<2*nc>(e15.8,2x))
1665 format (a,a,a,i0,a,i0)

    ni=size(tab,1)
    nj=size(tab,2)
    nc=size(tab,3)
    nx=size(x,3)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',ni,', J=',nj
    else
       write(filespec%fid,1665)'Zone T="", I=',ni,', J=',nj
    end if

    if (present(x)) then
       do j=1,nj
          do i=1,ni
             write(filespec%fid,1664)(x(i,j,ix),ix=1,nx),(tab(i,j,ic),ic=1,nc)
          end do
       end do
    else
       do j=1,nj
          do i=1,ni
             write(filespec%fid,1664)(tab(i,j,ic),ic=1,nc)
          end do
       end do
    end if
  end subroutine dwrite_ascii_2d

  subroutine iwrite_ascii_2d(filespec,x,tab,zonetitle)
    implicit none
    type(filetype)                                   ::filespec
    integer, intent(in), dimension(:,:,:)            ::tab
    integer, intent(in), dimension(:,:,:),  optional ::x
    character(len=*)   ,optional                     ::zonetitle
    integer(kind=8)                                  ::ni,nj,nn,nc,nx
    !----------------------------------------------------------------

1664 format (<2*nc>(e15.8,2x))
1665 format (a,a,a,i0,a,i0)

    ni=size(tab,1)
    nj=size(tab,2)
    nc=size(tab,3)
    nx=size(x,3)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',ni,', J=',nj
    else
       write(filespec%fid,1665)'Zone T="',&
            '", I=',ni,', J=',nj
    end if

    if (present(x)) then
       do j=1,nj
          do i=1,ni
             write(filespec%fid,1664)(x(i,j,ix),ix=1,nx),(tab(i,j,ic),ic=1,nc)
          end do
       end do
    else
       do j=1,nj
          do i=1,ni
             write(filespec%fid,1664)(tab(i,j,ic),ic=1,nc)
          end do
       end do
    end if

  end subroutine iwrite_ascii_2d

  subroutine dwrite_ascii_2d_1c(filespec,x,tab,zonetitle)
    implicit none
    type(filetype)                                        ::filespec
    real(kind=8), intent(in), dimension(:,:,:)            ::x
    real(kind=8), intent(in), dimension(:,:)              ::tab
    character(len=*)                          ,  optional ::zonetitle

    integer(kind=8)                                       ::ni,nj,nx
    !------------------------------------------------------------------

1664 format (<4>(e15.8,2x))
1665 format (a,a,a,i0,a,i0)

    ni=size(tab,1)
    nj=size(tab,2)
    nx=size(x,3)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),'", I=',ni,', J=',nj
    else
       write(filespec%fid,1665)'Zone T="", I=',ni,', J=',nj
    end if

    ! if (present(x)) then
    do j=1,nj
       do i=1,ni
          write(filespec%fid,1664)(x(i,j,ix),ix=1,nx),tab(i,j)
       end do
    end do
    ! else
    !    do j=1,nj
    !       do i=1,ni
    !          write(filespec%fid,1664)tab(i,j)
    !       end do
    !    end do
    ! end if

  end subroutine dwrite_ascii_2d_1c



  !==== Write 3D data =======================================
  !*
  !*
  !==========================================================
  subroutine fwrite_ascii_3d(filespec,x,tab,zonetitle)
    implicit none
    type(filetype)                                      ::filespec
    integer(kind=8)                                     ::ni,nj,nk,nc
    real, intent(in), dimension(:,:,:,:)                ::tab
    real, intent(in), dimension(:,:,:,:),    optional   ::x
    character(len=*)  ,optional                         ::zonetitle
    !--------------------------------------------------------------

1664 format (<2*nc>(e15.8,2x))
1665 format (a,a,a,i0,a,i0,a,i0)

    ni=size(tab,1)
    nj=size(tab,2)
    nk=size(tab,3)
    nc=size(tab,4)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),&
            '", I=',ni,', J=',nj,', K=',nk
    else
       write(filespec%fid,1665)'Zone T="',&
            '", I=',ni,', J=',nj,', K=',nk
    end if

    if (present(x)) then
       do k=1,nk
          do j=1,nj
             do i=1,ni
                write(filespec%fid,1664)(x(i,j,k,ic),ic=1,nc),(tab(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
    else
       do k=1,nk
          do j=1,nj
             do i=1,ni
                write(filespec%fid,1664)(tab(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
    end if
  end subroutine fwrite_ascii_3d

  subroutine dwrite_ascii_3d(filespec,x,tab,zonetitle)
    implicit none
    type(filetype)                                            ::filespec
    integer(kind=8)                                           ::ni,nj,nk,nc
    real(kind=8), intent(in), dimension(:,:,:,:)              ::tab
    real(kind=8), intent(in), dimension(:,:,:,:),    optional ::x
    character(len=*)   ,optional                              ::zonetitle
    !--------------------------------------------------------------------

1664 format (<nc>(e15.8,2x))
1665 format (a,a,a,i0,a,i0,a,i0)

    ni=size(tab,1)
    nj=size(tab,2)
    nk=size(tab,3)
    nc=size(tab,4)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),&
            '", I=',ni,', J=',nj,', K=',nk
    else
       write(filespec%fid,1665)'Zone T="',&
            '", I=',ni,', J=',nj,', K=',nk
    end if
    if (present(x)) then
       do k=1,nk
          do j=1,nj
             do i=1,ni
                write(filespec%fid,1664)(x(i,j,k,ic),ic=1,nc),(tab(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
    else
       do k=1,nk
          do j=1,nj
             do i=1,ni
                write(filespec%fid,1664)(tab(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
    end if

  end subroutine dwrite_ascii_3d

  subroutine iwrite_ascii_3d(filespec,x,tab,zonetitle)
    implicit none
    type(filetype)                                         ::filespec
    integer(kind=8)                                        ::ni,nj,nk,nc
    integer, intent(in), dimension(:,:,:,:)                ::tab
    integer, intent(in), dimension(:,:,:,:),    optional   ::x
    character(len=*)  ,optional                            ::zonetitle
    !-----------------------------------------------------------------

1664 format (<2*nc>(i5,2x))
1665 format (a,a,a,i0,a,i0,a,i0)

    ni=size(tab,1)
    nj=size(tab,2)
    nk=size(tab,3)
    nc=size(tab,4)

    if (present(zonetitle)) then
       write(filespec%fid,1665)'Zone T="',trim(zonetitle),&
            '", I=',ni,', J=',nj,', K=',nk
    else
       write(filespec%fid,1665)'Zone T="',&
            '", I=',ni,', J=',nj,', K=',nk
    end if

    if (present(x)) then
       do k=1,nk
          do j=1,nj
             do i=1,ni
                write(filespec%fid,1664)(x(i,j,k,ic),ic=1,nc),(tab(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
    else
       do k=1,nk
          do j=1,nj
             do i=1,ni
                write(filespec%fid,1664)(tab(i,j,k,ic),ic=1,nc)
             end do
          end do
       end do
    end if

  end subroutine iwrite_ascii_3d

end module tecplot_IO
