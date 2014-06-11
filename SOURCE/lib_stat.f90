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
  integer(kind=8), private ::i,j,k,ic,in,nn,nc,ii

  public :: cal_statistics

  private :: cal_stat_1d_1c

  interface cal_statistics
     module procedure cal_stat_1d_1c, cal_stat_1d_nc
  end interface cal_statistics

contains

  subroutine cal_stat_1d_1c(var,stat,weight)
    implicit none
    real(kind=8)    ,dimension(:)               ::var
    real(kind=8)    ,dimension(4)               ::stat
    real(kind=8)    ,dimension(:)    ,optional  ::weight
    real(kind=8)    ,dimension(:) ,allocatable  ::w
    !-------------------------------------------------

    nn = size(var,1)

    allocate(w(nn))
    if (present(weight)) then
       w(:) = weight(:)
    else
       w(:) = 1.d0
    end if

    stat(1) = sum(var(:)*w(:))/sum(w(:))
    stat(2) = sum((var(:)-stat(1))**2*w(:))/sum(w(:))
    stat(3) = sum(((var(:)-stat(1))/stat(2))**3*w(:))/sum(w(:))
    stat(4) = sum(((var(:)-stat(1))/stat(2))**4*w(:))/sum(w(:))

    return

  end subroutine cal_stat_1d_1c

  subroutine cal_stat_1d_nc(var,stat,weight)
    implicit none
    real(kind=8)    ,dimension(:,:)             ::var
    real(kind=8)    ,dimension(:)               ::stat
    real(kind=8)    ,dimension(:)    ,optional  ::weight
    real(kind=8)    ,dimension(:) ,allocatable  ::w
    integer(kind=8)                             ::nstat
    !----------------------------------------------------

    nn = size(var,1)
    nc = size(var,2)
    nstat = 0
    do i=1,nc
       do j=i,nc
          nstat = nstat+1
       end do
    end do
    nstat = nstat + 2*nc

    if (size(stat) .lt. nstat) then
       print*,'size(stat) must be ',nstat
       stop
    end if

    stat = 0.d0

    allocate(w(nn))
    if (present(weight)) then
       w(:) = weight(:)
    else
       w(:) = 1.d0
    end if

    !average
    do ic=1,nc
       stat(ic) = sum(var(:,ic)*w(:))/sum(w(:))
    end do

    !center the input variable
    do ic=1,nc
       var(:,ic) = var(:,ic) - stat(ic)
    end do

    !rms
    do i=nc+1,2*nc
       stat(i) = sum(var(:,i-nc)*var(:,i-nc)*w(:))/sum(w(:))
       stat(i) = sqrt(stat(i))
    end do

    !cross-moments
    ii=2*nc+1
    do i=1,nc-1
       do j=i+1,nc
          stat(ii) = sum(var(:,i)*var(:,j)*w(:))/sum(w(:))
          ii=ii+1
       end do
    end do

    !normalize the variable
    do i=1,nc
       var(:,i) = var(:,i)/stat(i+nc)
    end do

    !skewness
    do i=1,nc
       stat(ii) = sum(var(:,i)**3*w(:))/sum(w(:))
       ii = ii+1
    end do

    !flatness
    do i=1,nc
       stat(ii) = sum(var(:,i)**4*w(:))/sum(w(:))
       ii = ii+1
    end do

    !uncenter and unormallize var
    do i=1,nc
       var(:,i) = var(:,i)*stat(i+nc) + stat(i)
    end do

    return

  end subroutine cal_stat_1d_nc


end module lib_stat
