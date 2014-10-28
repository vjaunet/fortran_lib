module lib_stat
  implicit none

  !=================Specification=============================
  !*
  !* Module for Statisitic Routines
  !*
  !*
  !*  author : Vincent Jaunet
  !*  date   : 24-03-2014
  !*  License: MIT
  !*  contact: v.jaunet@gmail.com
  !*
  !!===========================================================
  integer(kind=8), private ::i,j,k,ic,in,nn,nc,ii,nx,ny

  public :: average,rms

  private :: d_average_1d_1c,d_rms_1d_1c,d_rms_1d_1c_moy,&
       f_average_1d_1c,f_rms_1d_1c,f_rms_1d_1c_moy

  interface average
     module procedure d_average_1d_1c,f_average_1d_1c
  end interface average

  interface rms
     module procedure d_rms_1d_1c,d_rms_1d_1c_moy
     module procedure f_rms_1d_1c,f_rms_1d_1c_moy
  end interface rms

contains

  subroutine d_average_1d_1c(var,moy,weight)
    implicit none
    real(kind=8)    ,dimension(:)                   ::var
    real(kind=8)                                    ::moy,sumw
    real(kind=8)    ,dimension(:)        ,optional  ::weight
    real(kind=8)    ,dimension(:)     ,allocatable  ::w
    !----------------------------------------------------------

    nn = size(var,1)

    allocate(w(nn))
    if (present(weight)) then
       w(:) = weight(:)
    else
       w(:) = 1.d0
    end if

    moy = 0.d0
    sumw = sum(w(:))
    if (sumw /= 0.d0) then
       moy = sum(var(:)*w(:))/sumw
    end if


  end subroutine d_average_1d_1c

  subroutine f_average_1d_1c(var,moy,weight)
    implicit none
    real(kind=4)    ,dimension(:)                   ::var
    real(kind=4)                                    ::moy,sumw
    real(kind=4)    ,dimension(:)        ,optional  ::weight
    real(kind=4)    ,dimension(:)     ,allocatable  ::w
    !----------------------------------------------------------

    nn = size(var,1)

    allocate(w(nn))
    if (present(weight)) then
       w(:) = weight(:)
    else
       w(:) = 1.d0
    end if

    moy = 0.d0
    sumw = sum(w(:))
    if (sumw /= 0.d0) then
       moy = sum(var(:)*w(:))/sumw
    end if

  end subroutine f_average_1d_1c


  subroutine d_rms_1d_1c(var,rms,weight)
    implicit none
    real(kind=8)    ,dimension(:)                   ::var
    real(kind=8)                                    ::rms
    real(kind=8)                                    ::moy, sumw
    real(kind=8)    ,dimension(:)        ,optional  ::weight
    real(kind=8)    ,dimension(:)     ,allocatable  ::w
    !----------------------------------------------------------

    nn = size(var,1)

    allocate(w(nn))
    if (present(weight)) then
       w(:) = weight(:)
    else
       w(:) = 1.d0
    end if

    call d_average_1d_1c(var,moy,w)

    rms = 0.d0
    sumw = sum(w(:))
    if (sumw /= 0.d0) then
       rms = sum((var(:)-moy)**2*w(:))/sumw
       rms = sqrt(rms)
    end if



  end subroutine d_rms_1d_1c

  subroutine f_rms_1d_1c(var,rms,weight)
    implicit none
    real(kind=4)    ,dimension(:)                   ::var
    real(kind=4)                                    ::rms
    real(kind=4)                                    ::moy, sumw
    real(kind=4)    ,dimension(:)        ,optional  ::weight
    real(kind=4)    ,dimension(:)     ,allocatable  ::w
    !----------------------------------------------------------

    nn = size(var,1)

    allocate(w(nn))
    if (present(weight)) then
       w(:) = weight(:)
    else
       w(:) = 1.d0
    end if

    call f_average_1d_1c(var,moy,w)

    rms = 0.d0
    sumw = sum(w(:))
    if (sumw /= 0.d0) then
       rms = sum((var(:)-moy)**2*w(:))/sumw
       rms = sqrt(rms)
    end if



  end subroutine f_rms_1d_1c

  subroutine d_rms_1d_1c_moy(var,moy,rms,weight)
    implicit none
    real(kind=8)    ,dimension(:)                   ::var
    real(kind=8)                                    ::rms,sumw
    real(kind=8)                                    ::moy
    real(kind=8)    ,dimension(:)        ,optional  ::weight
    real(kind=8)    ,dimension(:)     ,allocatable  ::w
    !----------------------------------------------------------

    nn = size(var,1)

    allocate(w(nn))
    if (present(weight)) then
       w(:) = weight(:)
    else
       w(:) = 1.d0
    end if

    rms = 0.d0
    sumw = sum(w(:))
    if (sumw /= 0.d0) then
       rms = sum((var(:)-moy)**2*w(:))/sumw
       rms = sqrt(rms)
    end if

  end subroutine d_rms_1d_1c_moy

  subroutine f_rms_1d_1c_moy(var,moy,rms,weight)
    implicit none
    real(kind=4)    ,dimension(:)                   ::var
    real(kind=4)                                    ::rms,sumw
    real(kind=4)                                    ::moy
    real(kind=4)    ,dimension(:)        ,optional  ::weight
    real(kind=4)    ,dimension(:)     ,allocatable  ::w
    !----------------------------------------------------------

    nn = size(var,1)

    allocate(w(nn))
    if (present(weight)) then
       w(:) = weight(:)
    else
       w(:) = 1.d0
    end if

    rms = 0.d0
    sumw = sum(w(:))
    if (sumw /= 0.d0) then
       rms = sum((var(:)-moy)**2*w(:))/sumw
       rms = sqrt(rms)
    end if

  end subroutine f_rms_1d_1c_moy



end module lib_stat
