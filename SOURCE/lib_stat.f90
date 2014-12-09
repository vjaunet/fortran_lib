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

  interface skewness
     module procedure f_skewness_1d_1c
  end interface skewness

  interface flatness
     module procedure f_flatness_1d_1c
  end interface flatness

  interface xmoment
     module procedure f_xmom_1d_1c
  end interface xmoment

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

    deallocate(w)

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

  subroutine f_skewness_1d_1c(var,skew,weight)
    implicit none
    real(kind=4)    ,dimension(:)                   ::var
    real(kind=4)                                    ::skew
    real(kind=4)                                    ::moy,rms,sumw
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
    call f_rms_1d_1c_moy(var,moy,rms,w)

    skew = 0.d0
    sumw = sum(w(:))
    if (sumw /= 0.d0) then
       skew = sum(w*((var(:)-moy)/rms)**3)
       skew = skew/sumw
    end if

  end subroutine f_skewness_1d_1c

  subroutine f_flatness_1d_1c(var,flat,weight)
    implicit none
    real(kind=4)    ,dimension(:)                   ::var
    real(kind=4)                                    ::flat
    real(kind=4)                                    ::moy,rms,sumw
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
    call f_rms_1d_1c_moy(var,moy,rms,w)

    flat = 0.d0
    sumw = sum(w(:))
    if (sumw /= 0.d0) then
       flat = sum(w*((var(:)-moy)/rms)**4) - 3.d0
       flat = flat/sumw
    end if

  end subroutine f_flatness_1d_1c

  subroutine f_xmom_1d_1c(var1,var2,xmom,weight)
    implicit none
    real(kind=4)    ,dimension(:)                   ::var1,var2
    real(kind=4)                                    ::moy1,rms1,xmom
    real(kind=4)                                    ::moy2,rms2,sumw
    real(kind=4)    ,dimension(:)        ,optional  ::weight
    real(kind=4)    ,dimension(:)     ,allocatable  ::w
    !----------------------------------------------------------

    nn = size(var1,1)

    allocate(w(nn))
    if (present(weight)) then
       w(:) = weight(:)
    else
       w(:) = 1.d0
    end if

    call f_average_1d_1c(var1,moy1,w)
    call f_average_1d_1c(var2,moy2,w)

    xmom = 0.d0
    sumw = sum(w(:))
    if (sumw /= 0.d0) then
       xmom = sum((var1-moy1)*(var2-moy2)*w)
       xmom = xmom/sumw
    end if

  end subroutine f_xmom_1d_1c


end module lib_stat
