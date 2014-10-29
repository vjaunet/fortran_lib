module interpol

  !=================Specification=============================
  !*
  !*       Module for Interpolating data
  !*
  !*  author : Vincent Jaunet
  !*  date   : 24-03-2014
  !*  License: MIT
  !*  contact: v.jaunet@gmail.com
  !*
  !*
  !=================Specification=============================

  !=================    Usage    =============================
  !*
  !* use interpol
  !*
  !*
  !* For 1d Data :
  !* call spline_inter(x(:),y(:),xint(:),yint(:))
  !*
  !* call inv_dist_inter(x(:),y(:),xint(:),yint(:),power,type)
  !*  - power is the power of the weigthing distance
  !*  - type has to be passed as a character(len=1) for a local
  !*    interpolation
  !*
  !* For 2d Data :
  !*
  !* call bicubic_inter(x(:,:,2),y(:,:),xint(:,:,2),yint(:,:))
  !*
  !*
  !*
  !=================    Usage    =============================

  !=================    TO DO    =============================
  !*
  !*          2D inverse distance
  !*
  !=================    TO DO    =============================

  integer(kind=8), private ::i,j,k,ic,in,iint,jint

  public :: bicubic_inter, &
       inv_dist_inter,&
       spline_inter

  private :: d_bicubic_interp_2d, diff, fillA, fillB, get_4neighbors,&
       d_inverse_dist_1d,&
       d_spline_inter_1d, d_spline_coef

  interface spline_inter
     module procedure d_spline_inter_1d
  end interface spline_inter

  interface inv_dist_inter
     module procedure d_inverse_dist_1d
  end interface inv_dist_inter

  interface bicubic_inter
     module procedure d_bicubic_interp_2d, f_bicubic_interp_2d
  end interface bicubic_inter

contains

  !=============Bicubic Interpolation of 2D data ===================
  !*
  !*         Computes bicubic interp of 2D data y = f(x1,x2)
  !*
  !================================================================
  subroutine d_bicubic_interp_2d(x,y,xint,yint)
    implicit none
    real(kind=8)    ,dimension(:,:,:)               ::x
    real(kind=8)    ,dimension(:,:)                 ::y

    real(kind=8)    ,dimension(:,:,:) ,allocatable  ::dy
    real(kind=8)    ,dimension(:,:)   ,allocatable  ::A
    real(kind=8)    ,dimension(:)     ,allocatable  ::B
    real(kind=8)    ,dimension(:)     ,allocatable  ::alpha

    real(kind=8)    ,dimension(:,:)   ,allocatable  ::dist
    real(kind=8)    ,dimension(2,2,2)               ::xtmp
    real(kind=8)    ,dimension(2,2,3)               ::dytmp
    real(kind=8)    ,dimension(2,2)                 ::ytmp
    integer(kind=8) ,dimension(2)                   ::itmp
    real(kind=8)                                    ::gdx,gdy

    real(kind=8)    ,dimension(:,:,:)               ::xint
    real(kind=8)    ,dimension(:,:)                 ::yint

    real(kind=8)                                    ::trash
    integer(kind=8)                                 ::nx,ny
    integer(kind=8)                                 ::nxint,nyint
    integer(kind=8)                                 ::ii
    !-------------------------------------------------------------

    nx = size(x,1)
    ny = size(x,2)
    nxint = size(xint,1)
    nyint = size(xint,2)

    allocate(dy(nx,ny,3))

    !compute gradients along x1 : dy/dx1
    call diff(y,nx,ny,dy(:,:,1),'X')
    !compute gradients along x2 ; dy/dx2
    call diff(y,nx,ny,dy(:,:,2),'Y')
    !computes cros gradients : d(dy/dx)/dy = dy/dx1dx2
    call diff(dy(:,:,1),nx,ny,dy(:,:,3),'Y')

    allocate(A(16,16))
    call fillA(A)
    allocate(B(16),alpha(16))

    allocate(dist(nx,ny))
    do iint=1,nxint
       do jint=1,nyint

          !computes distance between POI and all data
          dist = 0.d0
          do i=1,nx
             do j=1,ny

                do ic=1,2
                   dist(i,j) = dist(i,j) + &
                        (x(i,j,ic)-xint(iint,jint,ic))**2
                end do
                dist(i,j) = dsqrt(dist(i,j))

             end do
          end do

          !localization of the surrounding 4 interpolant points
          itmp = minloc(dist)
          if (minval(dist) ==0.d0) then
             !the interpolated point is loacated on an interpolant
             !get the exact value
             yint(iint,jint) = y(itmp(1),itmp(2))

          else

             !get the values of the surrounding interpolants
             call get_4neighbors(itmp,xint(iint,jint,:),&
                  x,y,dy,xtmp,ytmp,dytmp)

             !get the bicubic interpolation coefficient alpha
             call fillB(B,ytmp,dytmp(:,:,1),dytmp(:,:,2),dytmp(:,:,3))
             alpha = 0.d0
             do i=1,16
                do j=1,16
                   alpha(i) = alpha(i) + A(i,j)*B(j)
                end do
             end do

             !get distance to compute the interpolation
             gdx = xint(iint,jint,1)- xtmp(1,1,1)
             gdx = gdx/(xtmp(2,2,1) - xtmp(1,1,1))
             gdy = xint(iint,jint,2) - xtmp(1,1,2)
             gdy = gdy/(xtmp(2,2,2) - xtmp(1,1,2))

             !Interpolate
             yint(iint,jint) = 0.d0
             ii = 1
             do j=0,3
                do i=0,3
                   yint(iint,jint) = yint(iint,jint) + &
                        alpha(ii) * gdx**i * gdy**j
                   ii=ii+1
                end do
             end do

             !for debugging
             ! if (isnan(yint(iint,jint))) then
             !    print*,itmp(1),itmp(2)
             !    print*,xint(iint,jint,1),xint(iint,jint,1)
             !    ! pause
             !    yint(iint,jint) = -10.d0
             ! end if

          end if

       end do !end loop on interpolated data
    end do !end loop on interpolated data


    deallocate(dy)
    deallocate(A)
    deallocate(B)
    deallocate(alpha)
    deallocate(dist)

  end subroutine d_bicubic_interp_2d


  subroutine f_bicubic_interp_2d(x,y,xint,yint)
    implicit none
    real(kind=4)    ,dimension(:,:,:)               ::x
    real(kind=4)    ,dimension(:,:)                 ::y

    real(kind=8)    ,dimension(:,:,:) ,allocatable  ::dy
    real(kind=8)    ,dimension(:,:)   ,allocatable  ::A
    real(kind=8)    ,dimension(:)     ,allocatable  ::B
    real(kind=8)    ,dimension(:)     ,allocatable  ::alpha

    real(kind=8)    ,dimension(:,:)   ,allocatable  ::dist
    real(kind=8)    ,dimension(2,2,2)               ::xtmp
    real(kind=8)    ,dimension(2,2,3)               ::dytmp
    real(kind=8)    ,dimension(2,2)                 ::ytmp
    integer(kind=8) ,dimension(2)                   ::itmp
    real(kind=8)                                    ::gdx,gdy

    real(kind=4)    ,dimension(:,:,:)               ::xint
    real(kind=4)    ,dimension(:,:)                 ::yint

    real(kind=8)                                    ::trash
    integer(kind=8)                                 ::nx,ny
    integer(kind=8)                                 ::nxint,nyint
    integer(kind=8)                                 ::ii
    !-------------------------------------------------------------

    nx = size(x,1)
    ny = size(x,2)
    nxint = size(xint,1)
    nyint = size(xint,2)

    allocate(dy(nx,ny,3))

    !compute gradients along x1 : dy/dx1
    call diff(dble(y),nx,ny,dy(:,:,1),'X')
    !compute gradients along x2 ; dy/dx2
    call diff(dble(y),nx,ny,dy(:,:,2),'Y')
    !computes cros gradients : d(dy/dx)/dy = dy/dx1dx2
    call diff(dy(:,:,1),nx,ny,dy(:,:,3),'Y')

    allocate(A(16,16))
    call fillA(A)
    allocate(B(16),alpha(16))

    allocate(dist(nx,ny))
    do iint=1,nxint
       do jint=1,nyint

          !computes distance between POI and all data
          dist = 0.d0
          do i=1,nx
             do j=1,ny

                do ic=1,2
                   dist(i,j) = dist(i,j) + &
                        (x(i,j,ic)-xint(iint,jint,ic))**2
                end do
                dist(i,j) = dsqrt(dist(i,j))

             end do
          end do

          !localization of the surrounding 4 interpolant points
          itmp = minloc(dist)
          if (minval(dist) ==0.d0) then
             !the interpolated point is loacated on an interpolant
             !get the exact value
             yint(iint,jint) = y(itmp(1),itmp(2))

          else

             !get the values of the surrounding interpolants
             call get_4neighbors(itmp,dble(xint(iint,jint,:)),&
                  dble(x),dble(y),dy,xtmp,ytmp,dytmp)

             !get the bicubic interpolation coefficient alpha
             call fillB(B,ytmp,dytmp(:,:,1),dytmp(:,:,2),dytmp(:,:,3))
             alpha = 0.d0
             do i=1,16
                do j=1,16
                   alpha(i) = alpha(i) + A(i,j)*B(j)
                end do
             end do

             !get distance to compute the interpolation
             gdx = xint(iint,jint,1)- xtmp(1,1,1)
             gdx = gdx/(xtmp(2,2,1) - xtmp(1,1,1))
             gdy = xint(iint,jint,2) - xtmp(1,1,2)
             gdy = gdy/(xtmp(2,2,2) - xtmp(1,1,2))

             !Interpolate
             yint(iint,jint) = 0.d0
             ii = 1
             do j=0,3
                do i=0,3
                   yint(iint,jint) = yint(iint,jint) + &
                        alpha(ii) * gdx**i * gdy**j
                   ii=ii+1
                end do
             end do

             !for debugging
             ! if (isnan(yint(iint,jint))) then
             !    print*,itmp(1),itmp(2)
             !    print*,xint(iint,jint,1),xint(iint,jint,1)
             !    ! pause
             !    yint(iint,jint) = -10.d0
             ! end if

          end if

       end do !end loop on interpolated data
    end do !end loop on interpolated data


    deallocate(dy)
    deallocate(A)
    deallocate(B)
    deallocate(alpha)
    deallocate(dist)

  end subroutine f_bicubic_interp_2d


  !=============Spline Interpolation of 1D data ===================
  !*
  !*
  !================================================================

  subroutine d_spline_inter_1d(x,y,xint,yint)
    !
    !
    !cccc Routine d'interpolation
    !     En entree
    !     signal y=f(x) sur nbp points
    !     nbp coefficient
    !     nbpint abscisse xint d'interpolation
    !     En sortie :
    !     nbpint valeurs interpollees

    integer*8                            ::nbp,nbpint
    real*8                               ::x(:),y(:)
    real*8                               ::xint(:),yint(:)
    real*8 ,allocatable                  ::ypp(:)
    integer*8                            ::i,j,sens
    real*8                               ::a,b,dx

    nbp = size(y,1)
    nbpint = size(yint,1)

    allocate(ypp(nbp))

    ! Get spline coeficients
    call d_spline_coef(x,y,ypp,nbp)

    !     On suppose les tableaux x et xint sont ordonnes
    if(nbpint.eq.1) then
       sens=1
    else if(sign(1.0,x(2)-x(1)).eq.sign(1.0,xint(2)-xint(1))) then
       sens=1
    else
       sens=-1
    endif
    !
    if(sens.eq.1) then
       !       Meme orientation
       j=2
       do i=1,nbpint
          do while((x(j).le.xint(i)).and.(j.lt.nbp))
             j=j+1
          enddo
          dx=(x(j)-x(j-1))
          a=(x(j)-xint(i))/dx
          b=1.0-a
          yint(i)=a*y(j-1)+b*y(j)&
               +((a*a*a-a)*ypp(j-1)+(b*b*b-b)*ypp(j))*dx*dx/6.0
       enddo
    else if(sens.eq.-1) then
       !       Orientation opposee
       j=nbp-1
       do i=1,nbpint
          do while((x(j).ge.xint(i)).and.(j.gt.1))
             j=j-1
          enddo
          dx=(x(j+1)-x(j))
          a=(x(j+1)-xint(i))/dx
          b=1.0-a
          yint(i)=a*y(j)+b*y(j+1)&
               +((a*a*a-a)*ypp(j)+(b*b*b-b)*ypp(j+1))*dx*dx/6.0
       enddo
    endif

    deallocate(ypp)

    return

  end subroutine d_spline_inter_1d

  subroutine d_spline_coef(x,y,ypp,nbp)
    !cccc Routine de calcul des coefficients spline ypp
    !     En entree :
    !     signal y=f(x) sur nbp points
    !     En sortie :
    !     nbp coeffients ypp

    integer*8 nbp
    real*8 x(nbp),y(nbp),ypp(nbp)
    integer*8 i
    real*8,allocatable :: dx(:),dy(:),tmp(:)
    real*8 sigma,p

    allocate(dx(nbp-1),dy(nbp-1),tmp(nbp-1))
    !     Constitution du systeme a resoudre (spline naturelle : y''=0 aux bornes)
    do i=1,nbp-1
       dx(i)=x(i+1)-x(i)
       dy(i)=y(i+1)-y(i)
    enddo

    ypp(1)=0.0
    tmp(1)=0.0
    do i=2,nbp-1
       sigma=dx(i-1)/(dx(i)+dx(i-1))
       p=sigma*ypp(i-1)+2.0
       ypp(i)=(sigma-1.0)/p
       tmp(i)=(6.0*(dy(i)/dx(i)-dy(i-1)/dx(i-1))/(dx(i)+dx(i-1))&
            -sigma*tmp(i-1))/p
    enddo

    ypp(nbp)=0.0
    do i=nbp-1,1,-1
       ypp(i)=ypp(i)*ypp(i+1)+tmp(i)
    enddo
    deallocate(dx,dy,tmp);
    return
  end subroutine d_spline_coef


  !=============Inverse Distance Interpolation            ===================
  !*
  !*
  !==========================================================================

  subroutine d_inverse_dist_1d(x,y,xint,yint,pow,type)
    implicit none
    real(kind=8)     ,dimension(:)                 ::x
    real(kind=8)     ,dimension(:)                 ::y

    real(kind=8)                                   ::pow
    character(len=*) ,optional                     ::type

    real(kind=8)                                   ::wi,sumwi
    real(kind=8)     ,dimension(:)   ,allocatable  ::d
    logical                                        ::samept

    real(kind=8)     ,dimension(:)                 ::xint
    real(kind=8)     ,dimension(:)                 ::yint
    integer(kind=8)                                ::nx,nxint,j0
    !--------------------------------------------------------------

    nx = size(x,1)
    nxint = size(xint,1)

    allocate(d(nxint))

    if (.not. present(type)) then
       !Shepard(1968)
       yint = 0.d0
       do j=1,nxint
          sumwi   = 0.d0
          samept = .false.
          i=1
          do while(.not.samept .and. i<nx)
             wi = dsqrt((xint(j) - x(i))**2)
             if (wi /= 0.d0) then
                wi = wi**(-pow)
                sumwi = sumwi + wi
                yint(j) = yint(j) + wi*y(i)
                i = i+1
             else
                samept = .true.
             end if
          end do
          if (samept) then
             yint(j) = y(i)
          else
             yint(j) = yint(j)/sumwi
          end if
       end do

    else
       !Local Shepard method
       yint = 0.d0
       do j=1,nxint
          sumwi   = 0.d0
          samept = .false.
          i=1
          wi = dsqrt((xint(j) - x(i))**2)

          do while (wi > pow .and. wi /= 0.d0)
             i=i+1
             wi = dsqrt((xint(j) - x(i))**2)
          end do
          do while(wi <= pow .and. wi /= 0.d0 .and. i<nx)
             wi = ((pow - wi)/(pow*wi))**2
             sumwi = sumwi + wi
             yint(j) = yint(j) + wi*y(i)
             i = i+1
             wi = dsqrt((xint(j) - x(i))**2)
          end do

          if (wi == 0.d0) then
             yint(j) = y(i)
          else
             yint(j) = yint(j)/sumwi
          end if
       end do

    end if

  end subroutine d_inverse_dist_1d

  !-------------------------------------------------------------------
  !*
  !*                TOOLS
  !*
  !-------------------------------------------------------------------

  subroutine diff(tab,nx,ny,dtab,dir)
    implicit none
    integer(kind=8)                                    ::nx,ny
    real(kind=8)     ,dimension(nx,ny)                 ::tab,dtab
    character(len=1)                                   ::dir
    !--------------------------------------------------------------------

    !computes the second order Finite difference over a 2D domain

    if (dir == 'X') then
       do j=1,ny
          dtab(1,j)  = (-tab(3,j) + 4.d0*tab(2,j) + &
               3.d0*tab(1,j))/3.d0
          dtab(nx,j) = (tab(nx-2,j) - 4.d0*tab(nx-1,j) - &
               3.d0*tab(nx,j))/3.d0
       end do

       do i=2,nx-1
          do j=2,ny-1
             dtab(i,j)  = (tab(i+1,j) - tab(i-1,j))/2.d0
          end do
       end do

    else if (dir == "Y") then
       do i=1,nx
          dtab(i,1) = (-tab(i,3) + 4.d0*tab(i,2) + &
               3.d0*tab(i,1))/3.d0
          dtab(i,ny) = (tab(i,ny-2) - 4.d0*tab(i,ny-1) - &
               3.d0*tab(i,ny))/3.d0
       end do

       do i=2,nx-1
          do j=2,ny-1
             dtab(i,j)  = (tab(i,j+1) - tab(i,j-1))/2.d0
          end do
       end do

    else

       STOP "Differentiation direction unknown"

    end if
  end subroutine diff

  subroutine fillB(x,tab,dxtab,dytab,dxytab)
    real(kind=8) ,dimension(16)        ::x
    real(kind=8) ,dimension(2,2)       ::tab,dxtab,dytab,dxytab
    !----------------------------------------------------------

    x(1)  = tab(1,1)
    x(2)  = tab(2,1)
    x(3)  = tab(1,2)
    x(4)  = tab(2,2)
    x(5)  = dxtab(1,1)
    x(6)  = dxtab(2,1)
    x(7)  = dxtab(1,2)
    x(8)  = dxtab(2,2)
    x(9)  = dytab(1,1)
    x(10) = dytab(2,1)
    x(11) = dytab(1,2)
    x(12) = dytab(2,2)
    x(13) = dxytab(1,1)
    x(14) = dxytab(2,1)
    x(15) = dxytab(1,2)
    x(16) = dxytab(2,2)

  end subroutine fillB

  subroutine fillA(A)
    real(kind=8) ,dimension(:,:) ::A
    !-----------------------------------

    A = 0.d0
    A(1,1)  = 1.d0
    A(2,5)  = 1.d0
    A(3,1)  = -3.d0
    A(3,2)  =  3.d0
    A(3,5)  = -2.d0
    A(3,6)  = -1.d0
    A(4,1)  = 2.d0
    A(4,2)  = -2.d0
    A(4,5)  = 1.d0
    A(4,6)  = 1.d0
    A(5,9)  = 1.d0
    A(6,13) = 1.d0
    A(7,9)  = -3.d0
    A(7,10) = 3.d0
    A(7,13) = -2.d0
    A(7,14) = -1.d0
    A(8,9)  = 2.d0
    A(8,10) = -2.d0
    A(8,13) = 1.d0
    A(8,14) = 1.d0
    A(9,1)  = -3.d0
    A(9,3)  = 3.d0
    A(9,9)  = -2.d0
    A(9,11) = -1.d0
    A(10,5) = -3.d0
    A(10,7) = 3.d0
    A(10,13)= -2.d0
    A(10,15)= -1.d0
    A(11,1) = 9.d0
    A(11,2) = -9.d0
    A(11,3) = -9.d0
    A(11,4) = 9.d0
    A(11,5) = 6.d0
    A(11,6) = 3.d0
    A(11,7) = -6.d0
    A(11,8) = -3.d0
    A(11,9) = 6.d0
    A(11,10)= -6.d0
    A(11,11)= 3.d0
    A(11,12)= -3.d0
    A(11,13)= 4.d0
    A(11,14)= 2.d0
    A(11,15)= 2.d0
    A(11,16)= 1.d0
    A(12,1) = -6.d0
    A(12,2) = 6.d0
    A(12,3) = 6.d0
    A(12,4) = -6.d0
    A(12,5) = -3.d0
    A(12,6) = -3.d0
    A(12,7) = 3.d0
    A(12,8) = 3.d0
    A(12,9) = -4.d0
    A(12,10)= 4.d0
    A(12,11)= -2.d0
    A(12,12)= 2.d0
    A(12,13)= -2.d0
    A(12,14)= -2.d0
    A(12,15)= -1.d0
    A(12,16)= -1.d0
    A(13,1) = 2.d0
    A(13,3) = -2.d0
    A(13,9) = 1.d0
    A(13,11)= 1.d0
    A(14,5) = 2.d0
    A(14,7) = -2.d0
    A(14,13)= 1.d0
    A(14,15)= 1.d0

    A(15,1) = -6.d0
    A(15,2) = 6.d0
    A(15,3) = 6.d0
    A(15,4) = -6.d0
    A(15,5) = -4.d0
    A(15,6) = -2.d0
    A(15,7) = 4.d0
    A(15,8) = 2.d0
    A(15,9) = -3.d0
    A(15,10)= 3.d0
    A(15,11)= -3.d0
    A(15,12)= 3.d0
    A(15,13)= -2.d0
    A(15,14)= -1.d0
    A(15,15)= -2.d0
    A(15,16)= -1.d0

    A(16,1) = 4.d0
    A(16,2) = -4.d0
    A(16,3) = -4.d0
    A(16,4) = 4.d0
    A(16,5) = 2.d0
    A(16,6) = 2.d0
    A(16,7) = -2.d0
    A(16,8) = -2.d0
    A(16,9) = 2.d0
    A(16,10)= -2.d0
    A(16,11)= 2.d0
    A(16,12)= -2.d0
    A(16,13)= 1.d0
    A(16,14)= 1.d0
    A(16,15)= 1.d0
    A(16,16)= 1.d0

  end subroutine fillA


  subroutine get_4neighbors(itmp,xint,x,y,dy,xtmp,ytmp,dytmp)
    implicit none
    integer(kind=8)    ,dimension(:)           ::itmp
    real(kind=8)       ,dimension(:)           ::xint
    real(kind=8)       ,dimension(:,:,:)       ::x
    real(kind=8)       ,dimension(:,:)         ::y
    real(kind=8)       ,dimension(:,:,:)       ::dy
    real(kind=8)       ,dimension(:,:,:)       ::xtmp
    real(kind=8)       ,dimension(:,:)         ::ytmp
    real(kind=8)       ,dimension(:,:,:)       ::dytmp

    integer(kind=8)                            ::ideb,ifin
    integer(kind=8)                            ::jdeb,jfin
    !---------------------------------------------------------

    ideb=1
    jdeb=1
    ifin = size(y,1)
    jfin = size(y,2)

    !**************************************
    !*         handling of corners
    !**************************************

    if (itmp(1) == ideb .and. itmp(2) == jdeb)  then
       !the closest interpolant is in the lower left corner
       xtmp(1,1,1) = x(itmp(1),itmp(2),1)
       xtmp(1,1,2) = x(itmp(1),itmp(2),2)
       ytmp(1,1)   = y(itmp(1),itmp(2))
       dytmp(1,1,1)= dy(itmp(1),itmp(2),1)
       dytmp(1,1,2)= dy(itmp(1),itmp(2),2)
       dytmp(1,1,3)= dy(itmp(1),itmp(2),3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1),itmp(2)+1,1)
       xtmp (1,2,2) = x(itmp(1),itmp(2)+1,2)
       ytmp (1,2)   = y(itmp(1),itmp(2)+1)
       dytmp(1,2,1)= dy(itmp(1),itmp(2)+1,1)
       dytmp(1,2,2)= dy(itmp(1),itmp(2)+1,2)
       dytmp(1,2,3)= dy(itmp(1),itmp(2)+1,3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1)+1,itmp(2)+1,1)
       xtmp (2,2,2) = x(itmp(1)+1,itmp(2)+1,2)
       ytmp (2,2)   = y(itmp(1)+1,itmp(2)+1)
       dytmp(2,2,1)= dy(itmp(1)+1,itmp(2)+1,1)
       dytmp(2,2,2)= dy(itmp(1)+1,itmp(2)+1,2)
       dytmp(2,2,3)= dy(itmp(1)+1,itmp(2)+1,3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1)+1,itmp(2),1)
       xtmp (2,1,2) = x(itmp(1)+1,itmp(2),2)
       ytmp (2,1)   = y(itmp(1)+1,itmp(2))
       dytmp(2,1,1)= dy(itmp(1)+1,itmp(2),1)
       dytmp(2,1,2)= dy(itmp(1)+1,itmp(2),2)
       dytmp(2,1,3)= dy(itmp(1)+1,itmp(2),3)

       return

    else if (itmp(1)==ideb .and. itmp(2) == jfin) then

       !the closest interpolant is in the lower left corner
       xtmp(1,1,1) = x(itmp(1) ,itmp(2)-1,1)
       xtmp(1,1,2) = x(itmp(1) ,itmp(2)-1,2)
       ytmp(1,1)   = y(itmp(1) ,itmp(2)-1)
       dytmp(1,1,1)= dy(itmp(1),itmp(2)-1,1)
       dytmp(1,1,2)= dy(itmp(1),itmp(2)-1,2)
       dytmp(1,1,3)= dy(itmp(1),itmp(2)-1,3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1),itmp(2),1)
       xtmp (1,2,2) = x(itmp(1),itmp(2),2)
       ytmp (1,2)   = y(itmp(1),itmp(2))
       dytmp(1,2,1)= dy(itmp(1),itmp(2),1)
       dytmp(1,2,2)= dy(itmp(1),itmp(2),2)
       dytmp(1,2,3)= dy(itmp(1),itmp(2),3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1)+1,itmp(2),1)
       xtmp (2,2,2) = x(itmp(1)+1,itmp(2),2)
       ytmp (2,2)   = y(itmp(1)+1,itmp(2))
       dytmp(2,2,1)= dy(itmp(1)+1,itmp(2),1)
       dytmp(2,2,2)= dy(itmp(1)+1,itmp(2),2)
       dytmp(2,2,3)= dy(itmp(1)+1,itmp(2),3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1)+1,itmp(2)-1,1)
       xtmp (2,1,2) = x(itmp(1)+1,itmp(2)-1,2)
       ytmp (2,1)   = y(itmp(1)+1,itmp(2)-1)
       dytmp(2,1,1)= dy(itmp(1)+1,itmp(2)-1,1)
       dytmp(2,1,2)= dy(itmp(1)+1,itmp(2)-1,2)
       dytmp(2,1,3)= dy(itmp(1)+1,itmp(2)-1,3)

       return


    else if (itmp(1)==ifin .and. itmp(2) == jdeb) then

       !lower left corner
       xtmp(1,1,1) = x(itmp(1) -1,itmp(2),1)
       xtmp(1,1,2) = x(itmp(1) -1,itmp(2),2)
       ytmp(1,1)   = y(itmp(1) -1,itmp(2))
       dytmp(1,1,1)= dy(itmp(1)-1,itmp(2),1)
       dytmp(1,1,2)= dy(itmp(1)-1,itmp(2),2)
       dytmp(1,1,3)= dy(itmp(1)-1,itmp(2),3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1)-1,itmp(2)+1,1)
       xtmp (1,2,2) = x(itmp(1)-1,itmp(2)+1,2)
       ytmp (1,2)   = y(itmp(1)-1,itmp(2)+1)
       dytmp(1,2,1)= dy(itmp(1)-1,itmp(2)+1,1)
       dytmp(1,2,2)= dy(itmp(1)-1,itmp(2)+1,2)
       dytmp(1,2,3)= dy(itmp(1)-1,itmp(2)+1,3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1),itmp(2)+1,1)
       xtmp (2,2,2) = x(itmp(1),itmp(2)+1,2)
       ytmp (2,2)   = y(itmp(1),itmp(2)+1)
       dytmp(2,2,1)= dy(itmp(1),itmp(2)+1,1)
       dytmp(2,2,2)= dy(itmp(1),itmp(2)+1,2)
       dytmp(2,2,3)= dy(itmp(1),itmp(2)+1,3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1),itmp(2),1)
       xtmp (2,1,2) = x(itmp(1),itmp(2),2)
       ytmp (2,1)   = y(itmp(1),itmp(2))
       dytmp(2,1,1)= dy(itmp(1),itmp(2),1)
       dytmp(2,1,2)= dy(itmp(1),itmp(2),2)
       dytmp(2,1,3)= dy(itmp(1),itmp(2),3)

       return

    else if (itmp(1)==ifin .and. itmp(2) == jfin) then

       !lower left corner
       xtmp(1,1,1) = x(itmp(1) -1,itmp(2)-1,1)
       xtmp(1,1,2) = x(itmp(1) -1,itmp(2)-1,2)
       ytmp(1,1)   = y(itmp(1) -1,itmp(2)-1)
       dytmp(1,1,1)= dy(itmp(1)-1,itmp(2)-1,1)
       dytmp(1,1,2)= dy(itmp(1)-1,itmp(2)-1,2)
       dytmp(1,1,3)= dy(itmp(1)-1,itmp(2)-1,3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1)-1,itmp(2),1)
       xtmp (1,2,2) = x(itmp(1)-1,itmp(2),2)
       ytmp (1,2)   = y(itmp(1)-1,itmp(2))
       dytmp(1,2,1)= dy(itmp(1)-1,itmp(2),1)
       dytmp(1,2,2)= dy(itmp(1)-1,itmp(2),2)
       dytmp(1,2,3)= dy(itmp(1)-1,itmp(2),3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1),itmp(2),1)
       xtmp (2,2,2) = x(itmp(1),itmp(2),2)
       ytmp (2,2)   = y(itmp(1),itmp(2))
       dytmp(2,2,1)= dy(itmp(1),itmp(2),1)
       dytmp(2,2,2)= dy(itmp(1),itmp(2),2)
       dytmp(2,2,3)= dy(itmp(1),itmp(2),3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1),itmp(2)-1,1)
       xtmp (2,1,2) = x(itmp(1),itmp(2)-1,2)
       ytmp (2,1)   = y(itmp(1),itmp(2)-1)
       dytmp(2,1,1)= dy(itmp(1),itmp(2)-1,1)
       dytmp(2,1,2)= dy(itmp(1),itmp(2)-1,2)
       dytmp(2,1,3)= dy(itmp(1),itmp(2)-1,3)

       return

    end if


    !**************************************
    !*       handling of boundaries
    !**************************************
    if (itmp(1) == ideb) then

       if (x(itmp(1),itmp(2),2) < xint(2)) then

          !lower left corner
          xtmp(1,1,1)  = x(itmp(1),itmp(2),1)
          xtmp(1,1,2)  = x(itmp(1),itmp(2),2)
          ytmp(1,1)    = y(itmp(1),itmp(2))
          dytmp(1,1,1)= dy(itmp(1),itmp(2),1)
          dytmp(1,1,2)= dy(itmp(1),itmp(2),2)
          dytmp(1,1,3)= dy(itmp(1),itmp(2),3)

          !upper left corner
          xtmp (1,2,1) = x(itmp(1),itmp(2)+1,1)
          xtmp (1,2,2) = x(itmp(1),itmp(2)+1,2)
          ytmp (1,2)   = y(itmp(1),itmp(2)+1)
          dytmp(1,2,1)= dy(itmp(1),itmp(2)+1,1)
          dytmp(1,2,2)= dy(itmp(1),itmp(2)+1,2)
          dytmp(1,2,3)= dy(itmp(1),itmp(2)+1,3)

          !upper right corner
          xtmp (2,2,1) = x(itmp(1)+1,itmp(2)+1,1)
          xtmp (2,2,2) = x(itmp(1)+1,itmp(2)+1,2)
          ytmp (2,2)   = y(itmp(1)+1,itmp(2)+1)
          dytmp(2,2,1)= dy(itmp(1)+1,itmp(2)+1,1)
          dytmp(2,2,2)= dy(itmp(1)+1,itmp(2)+1,2)
          dytmp(2,2,3)= dy(itmp(1)+1,itmp(2)+1,3)

          !lower right corner
          xtmp (2,1,1) = x(itmp(1)+1,itmp(2),1)
          xtmp (2,1,2) = x(itmp(1)+1,itmp(2),2)
          ytmp (2,1)   = y(itmp(1)+1,itmp(2))
          dytmp(2,1,1)= dy(itmp(1)+1,itmp(2),1)
          dytmp(2,1,2)= dy(itmp(1)+1,itmp(2),2)
          dytmp(2,1,3)= dy(itmp(1)+1,itmp(2),3)

       else

          !lower left corner
          xtmp(1,1,1)  = x(itmp(1),itmp(2)-1,1)
          xtmp(1,1,2)  = x(itmp(1),itmp(2)-1,2)
          ytmp(1,1)    = y(itmp(1),itmp(2)-1)
          dytmp(1,1,1)= dy(itmp(1),itmp(2)-1,1)
          dytmp(1,1,2)= dy(itmp(1),itmp(2)-1,2)
          dytmp(1,1,3)= dy(itmp(1),itmp(2)-1,3)

          !upper left corner
          xtmp (1,2,1) = x(itmp(1),itmp(2),1)
          xtmp (1,2,2) = x(itmp(1),itmp(2),2)
          ytmp (1,2)   = y(itmp(1),itmp(2))
          dytmp(1,2,1)= dy(itmp(1),itmp(2),1)
          dytmp(1,2,2)= dy(itmp(1),itmp(2),2)
          dytmp(1,2,3)= dy(itmp(1),itmp(2),3)

          !upper right corner
          xtmp (2,2,1) = x(itmp(1)+1,itmp(2),1)
          xtmp (2,2,2) = x(itmp(1)+1,itmp(2),2)
          ytmp (2,2)   = y(itmp(1)+1,itmp(2))
          dytmp(2,2,1)= dy(itmp(1)+1,itmp(2),1)
          dytmp(2,2,2)= dy(itmp(1)+1,itmp(2),2)
          dytmp(2,2,3)= dy(itmp(1)+1,itmp(2),3)

          !lower right corner
          xtmp (2,1,1) = x(itmp(1)+1,itmp(2),1)
          xtmp (2,1,2) = x(itmp(1)+1,itmp(2),2)
          ytmp (2,1)   = y(itmp(1)+1,itmp(2))
          dytmp(2,1,1)= dy(itmp(1)+1,itmp(2),1)
          dytmp(2,1,2)= dy(itmp(1)+1,itmp(2),2)
          dytmp(2,1,3)= dy(itmp(1)+1,itmp(2),3)

       end if

       return

    end if

    if (itmp(1) == ifin) then

       if (x(itmp(1),itmp(2),2) < xint(2)) then

          !lower left corner
          xtmp(1,1,1)  = x(itmp(1)-1,itmp(2),1)
          xtmp(1,1,2)  = x(itmp(1)-1,itmp(2),2)
          ytmp(1,1)    = y(itmp(1)-1,itmp(2))
          dytmp(1,1,1)= dy(itmp(1)-1,itmp(2),1)
          dytmp(1,1,2)= dy(itmp(1)-1,itmp(2),2)
          dytmp(1,1,3)= dy(itmp(1)-1,itmp(2),3)

          !upper left corner
          xtmp (1,2,1) = x(itmp(1)-1,itmp(2)+1,1)
          xtmp (1,2,2) = x(itmp(1)-1,itmp(2)+1,2)
          ytmp (1,2)   = y(itmp(1)-1,itmp(2)+1)
          dytmp(1,2,1)= dy(itmp(1)-1,itmp(2)+1,1)
          dytmp(1,2,2)= dy(itmp(1)-1,itmp(2)+1,2)
          dytmp(1,2,3)= dy(itmp(1)-1,itmp(2)+1,3)

          !upper right corner
          xtmp (2,2,1) = x(itmp(1),itmp(2)+1,1)
          xtmp (2,2,2) = x(itmp(1),itmp(2)+1,2)
          ytmp (2,2)   = y(itmp(1),itmp(2)+1)
          dytmp(2,2,1)= dy(itmp(1),itmp(2)+1,1)
          dytmp(2,2,2)= dy(itmp(1),itmp(2)+1,2)
          dytmp(2,2,3)= dy(itmp(1),itmp(2)+1,3)

          !lower right corner
          xtmp (2,1,1) = x(itmp(1),itmp(2),1)
          xtmp (2,1,2) = x(itmp(1),itmp(2),2)
          ytmp (2,1)   = y(itmp(1),itmp(2))
          dytmp(2,1,1)= dy(itmp(1),itmp(2),1)
          dytmp(2,1,2)= dy(itmp(1),itmp(2),2)
          dytmp(2,1,3)= dy(itmp(1),itmp(2),3)

       else

          !lower left corner
          xtmp(1,1,1)  = x(itmp(1)-1,itmp(2)-1,1)
          xtmp(1,1,2)  = x(itmp(1)-1,itmp(2)-1,2)
          ytmp(1,1)    = y(itmp(1)-1,itmp(2)-1)
          dytmp(1,1,1)= dy(itmp(1)-1,itmp(2)-1,1)
          dytmp(1,1,2)= dy(itmp(1)-1,itmp(2)-1,2)
          dytmp(1,1,3)= dy(itmp(1)-1,itmp(2)-1,3)

          !upper left corner
          xtmp (1,2,1) = x(itmp(1)-1,itmp(2),1)
          xtmp (1,2,2) = x(itmp(1)-1,itmp(2),2)
          ytmp (1,2)   = y(itmp(1)-1,itmp(2))
          dytmp(1,2,1)= dy(itmp(1)-1,itmp(2),1)
          dytmp(1,2,2)= dy(itmp(1)-1,itmp(2),2)
          dytmp(1,2,3)= dy(itmp(1)-1,itmp(2),3)

          !upper right corner
          xtmp (2,2,1) = x(itmp(1),itmp(2),1)
          xtmp (2,2,2) = x(itmp(1),itmp(2),2)
          ytmp (2,2)   = y(itmp(1),itmp(2))
          dytmp(2,2,1)= dy(itmp(1),itmp(2),1)
          dytmp(2,2,2)= dy(itmp(1),itmp(2),2)
          dytmp(2,2,3)= dy(itmp(1),itmp(2),3)

          !lower right corner
          xtmp (2,1,1) = x(itmp(1),itmp(2)-1,1)
          xtmp (2,1,2) = x(itmp(1),itmp(2)-1,2)
          ytmp (2,1)   = y(itmp(1),itmp(2)-1)
          dytmp(2,1,1)= dy(itmp(1),itmp(2)-1,1)
          dytmp(2,1,2)= dy(itmp(1),itmp(2)-1,2)
          dytmp(2,1,3)= dy(itmp(1),itmp(2)-1,3)

       end if

       return

    end if

    if (itmp(2) == jdeb) then

       if (x(itmp(1),itmp(2),1) < xint(1)) then

          !lower left corner
          xtmp(1,1,1)  = x(itmp(1),itmp(2),1)
          xtmp(1,1,2)  = x(itmp(1),itmp(2),2)
          ytmp(1,1)    = y(itmp(1),itmp(2))
          dytmp(1,1,1)= dy(itmp(1),itmp(2),1)
          dytmp(1,1,2)= dy(itmp(1),itmp(2),2)
          dytmp(1,1,3)= dy(itmp(1),itmp(2),3)

          !upper left corner
          xtmp (1,2,1) = x(itmp(1),itmp(2)+1,1)
          xtmp (1,2,2) = x(itmp(1),itmp(2)+1,2)
          ytmp (1,2)   = y(itmp(1),itmp(2)+1)
          dytmp(1,2,1)= dy(itmp(1),itmp(2)+1,1)
          dytmp(1,2,2)= dy(itmp(1),itmp(2)+1,2)
          dytmp(1,2,3)= dy(itmp(1),itmp(2)+1,3)

          !upper right corner
          xtmp (2,2,1) = x(itmp(1)+1,itmp(2)+1,1)
          xtmp (2,2,2) = x(itmp(1)+1,itmp(2)+1,2)
          ytmp (2,2)   = y(itmp(1)+1,itmp(2)+1)
          dytmp(2,2,1)= dy(itmp(1)+1,itmp(2)+1,1)
          dytmp(2,2,2)= dy(itmp(1)+1,itmp(2)+1,2)
          dytmp(2,2,3)= dy(itmp(1)+1,itmp(2)+1,3)

          !lower right corner
          xtmp (2,1,1) = x(itmp(1)+1,itmp(2),1)
          xtmp (2,1,2) = x(itmp(1)+1,itmp(2),2)
          ytmp (2,1)   = y(itmp(1)+1,itmp(2))
          dytmp(2,1,1)= dy(itmp(1)+1,itmp(2),1)
          dytmp(2,1,2)= dy(itmp(1)+1,itmp(2),2)
          dytmp(2,1,3)= dy(itmp(1)+1,itmp(2),3)

       else

          !lower left corner
          xtmp(1,1,1)  = x(itmp(1)-1,itmp(2),1)
          xtmp(1,1,2)  = x(itmp(1)-1,itmp(2),2)
          ytmp(1,1)    = y(itmp(1)-1,itmp(2))
          dytmp(1,1,1)= dy(itmp(1)-1,itmp(2),1)
          dytmp(1,1,2)= dy(itmp(1)-1,itmp(2),2)
          dytmp(1,1,3)= dy(itmp(1)-1,itmp(2),3)

          !upper left corner
          xtmp (1,2,1) = x(itmp(1)-1,itmp(2)+1,1)
          xtmp (1,2,2) = x(itmp(1)-1,itmp(2)+1,2)
          ytmp (1,2)   = y(itmp(1)-1,itmp(2)+1)
          dytmp(1,2,1)= dy(itmp(1)-1,itmp(2)+1,1)
          dytmp(1,2,2)= dy(itmp(1)-1,itmp(2)+1,2)
          dytmp(1,2,3)= dy(itmp(1)-1,itmp(2)+1,3)

          !upper right corner
          xtmp (2,2,1) = x(itmp(1),itmp(2)+1,1)
          xtmp (2,2,2) = x(itmp(1),itmp(2)+1,2)
          ytmp (2,2)   = y(itmp(1),itmp(2)+1)
          dytmp(2,2,1)= dy(itmp(1),itmp(2)+1,1)
          dytmp(2,2,2)= dy(itmp(1),itmp(2)+1,2)
          dytmp(2,2,3)= dy(itmp(1),itmp(2)+1,3)

          !lower right corner
          xtmp (2,1,1) = x(itmp(1),itmp(2),1)
          xtmp (2,1,2) = x(itmp(1),itmp(2),2)
          ytmp (2,1)   = y(itmp(1),itmp(2))
          dytmp(2,1,1)= dy(itmp(1),itmp(2),1)
          dytmp(2,1,2)= dy(itmp(1),itmp(2),2)
          dytmp(2,1,3)= dy(itmp(1),itmp(2),3)

       end if

       return

    end if

    if (itmp(2) == jfin) then

       if (x(itmp(1),itmp(2),1) < xint(1)) then

          !lower left corner
          xtmp(1,1,1)  = x(itmp(1),itmp(2)-1,1)
          xtmp(1,1,2)  = x(itmp(1),itmp(2)-1,2)
          ytmp(1,1)    = y(itmp(1),itmp(2)-1)
          dytmp(1,1,1)= dy(itmp(1),itmp(2)-1,1)
          dytmp(1,1,2)= dy(itmp(1),itmp(2)-1,2)
          dytmp(1,1,3)= dy(itmp(1),itmp(2)-1,3)

          !upper left corner
          xtmp (1,2,1) = x(itmp(1),itmp(2),1)
          xtmp (1,2,2) = x(itmp(1),itmp(2),2)
          ytmp (1,2)   = y(itmp(1),itmp(2))
          dytmp(1,2,1)= dy(itmp(1),itmp(2),1)
          dytmp(1,2,2)= dy(itmp(1),itmp(2),2)
          dytmp(1,2,3)= dy(itmp(1),itmp(2),3)

          !upper right corner
          xtmp (2,2,1) = x(itmp(1)+1,itmp(2),1)
          xtmp (2,2,2) = x(itmp(1)+1,itmp(2),2)
          ytmp (2,2)   = y(itmp(1)+1,itmp(2))
          dytmp(2,2,1)= dy(itmp(1)+1,itmp(2),1)
          dytmp(2,2,2)= dy(itmp(1)+1,itmp(2),2)
          dytmp(2,2,3)= dy(itmp(1)+1,itmp(2),3)

          !lower right corner
          xtmp (2,1,1) = x(itmp(1)+1,itmp(2)-1,1)
          xtmp (2,1,2) = x(itmp(1)+1,itmp(2)-1,2)
          ytmp (2,1)   = y(itmp(1)+1,itmp(2)-1)
          dytmp(2,1,1)= dy(itmp(1)+1,itmp(2)-1,1)
          dytmp(2,1,2)= dy(itmp(1)+1,itmp(2)-1,2)
          dytmp(2,1,3)= dy(itmp(1)+1,itmp(2)-1,3)

       else

          !lower left corner
          xtmp(1,1,1)  = x(itmp(1)-1,itmp(2)-1,1)
          xtmp(1,1,2)  = x(itmp(1)-1,itmp(2)-1,2)
          ytmp(1,1)    = y(itmp(1)-1,itmp(2)-1)
          dytmp(1,1,1)= dy(itmp(1)-1,itmp(2)-1,1)
          dytmp(1,1,2)= dy(itmp(1)-1,itmp(2)-1,2)
          dytmp(1,1,3)= dy(itmp(1)-1,itmp(2)-1,3)

          !upper left corner
          xtmp (1,2,1) = x(itmp(1)-1,itmp(2),1)
          xtmp (1,2,2) = x(itmp(1)-1,itmp(2),2)
          ytmp (1,2)   = y(itmp(1)-1,itmp(2))
          dytmp(1,2,1)= dy(itmp(1)-1,itmp(2),1)
          dytmp(1,2,2)= dy(itmp(1)-1,itmp(2),2)
          dytmp(1,2,3)= dy(itmp(1)-1,itmp(2),3)

          !upper right corner
          xtmp (2,2,1) = x(itmp(1),itmp(2),1)
          xtmp (2,2,2) = x(itmp(1),itmp(2),2)
          ytmp (2,2)   = y(itmp(1),itmp(2))
          dytmp(2,2,1)= dy(itmp(1),itmp(2),1)
          dytmp(2,2,2)= dy(itmp(1),itmp(2),2)
          dytmp(2,2,3)= dy(itmp(1),itmp(2),3)

          !lower right corner
          xtmp (2,1,1) = x(itmp(1),itmp(2)-1,1)
          xtmp (2,1,2) = x(itmp(1),itmp(2)-1,2)
          ytmp (2,1)   = y(itmp(1),itmp(2)-1)
          dytmp(2,1,1)= dy(itmp(1),itmp(2)-1,1)
          dytmp(2,1,2)= dy(itmp(1),itmp(2)-1,2)
          dytmp(2,1,3)= dy(itmp(1),itmp(2)-1,3)

       end if

       return

    end if



    !**************************************
    !*       handling of general case
    !**************************************

    if (x(itmp(1),itmp(2),1) < xint(1) .and. &
         x(itmp(1),itmp(2),2) < xint(2)) then
       !the interpolated pt is in the the lower left corner

       !lower left corner
       xtmp(1,1,1)  = x(itmp(1),itmp(2),1)
       xtmp(1,1,2)  = x(itmp(1),itmp(2),2)
       ytmp(1,1)    = y(itmp(1),itmp(2))
       dytmp(1,1,1)= dy(itmp(1),itmp(2),1)
       dytmp(1,1,2)= dy(itmp(1),itmp(2),2)
       dytmp(1,1,3)= dy(itmp(1),itmp(2),3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1),itmp(2)+1,1)
       xtmp (1,2,2) = x(itmp(1),itmp(2)+1,2)
       ytmp (1,2)   = y(itmp(1),itmp(2)+1)
       dytmp(1,2,1)= dy(itmp(1),itmp(2)+1,1)
       dytmp(1,2,2)= dy(itmp(1),itmp(2)+1,2)
       dytmp(1,2,3)= dy(itmp(1),itmp(2)+1,3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1)+1,itmp(2)+1,1)
       xtmp (2,2,2) = x(itmp(1)+1,itmp(2)+1,2)
       ytmp (2,2)   = y(itmp(1)+1,itmp(2)+1)
       dytmp(2,2,1)= dy(itmp(1)+1,itmp(2)+1,1)
       dytmp(2,2,2)= dy(itmp(1)+1,itmp(2)+1,2)
       dytmp(2,2,3)= dy(itmp(1)+1,itmp(2)+1,3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1)+1,itmp(2),1)
       xtmp (2,1,2) = x(itmp(1)+1,itmp(2),2)
       ytmp (2,1)   = y(itmp(1)+1,itmp(2))
       dytmp(2,1,1)= dy(itmp(1)+1,itmp(2),1)
       dytmp(2,1,2)= dy(itmp(1)+1,itmp(2),2)
       dytmp(2,1,3)= dy(itmp(1)+1,itmp(2),3)


    else if (x(itmp(1),itmp(2),1) < xint(1) .and. &
         x(itmp(1),itmp(2),2) > xint(2)) then

       !the interpolated pt is in the the upper left

       !lower left corner
       xtmp(1,1,1)  = x(itmp(1),itmp(2)-1,1)
       xtmp(1,1,2)  = x(itmp(1),itmp(2)-1,2)
       ytmp(1,1)    = y(itmp(1),itmp(2)-1)
       dytmp(1,1,1)= dy(itmp(1),itmp(2)-1,1)
       dytmp(1,1,2)= dy(itmp(1),itmp(2)-1,2)
       dytmp(1,1,3)= dy(itmp(1),itmp(2)-1,3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1),itmp(2),1)
       xtmp (1,2,2) = x(itmp(1),itmp(2),2)
       ytmp (1,2)   = y(itmp(1),itmp(2))
       dytmp(1,2,1)= dy(itmp(1),itmp(2),1)
       dytmp(1,2,2)= dy(itmp(1),itmp(2),2)
       dytmp(1,2,3)= dy(itmp(1),itmp(2),3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1)+1,itmp(2),1)
       xtmp (2,2,2) = x(itmp(1)+1,itmp(2),2)
       ytmp (2,2)   = y(itmp(1)+1,itmp(2))
       dytmp(2,2,1)= dy(itmp(1)+1,itmp(2),1)
       dytmp(2,2,2)= dy(itmp(1)+1,itmp(2),2)
       dytmp(2,2,3)= dy(itmp(1)+1,itmp(2),3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1)+1,itmp(2)-1,1)
       xtmp (2,1,2) = x(itmp(1)+1,itmp(2)-1,2)
       ytmp (2,1)   = y(itmp(1)+1,itmp(2)-1)
       dytmp(2,1,1)= dy(itmp(1)+1,itmp(2)-1,1)
       dytmp(2,1,2)= dy(itmp(1)+1,itmp(2)-1,2)
       dytmp(2,1,3)= dy(itmp(1)+1,itmp(2)-1,3)

       return

    else if (x(itmp(1),itmp(2),1) > xint(1) .and. &
         x(itmp(1),itmp(2),2) < xint(2)) then

       !the interpolated pt is in the the lower right corner

       !lower left corner
       xtmp(1,1,1)  = x(itmp(1)-1,itmp(2),1)
       xtmp(1,1,2)  = x(itmp(1)-1,itmp(2),2)
       ytmp(1,1)    = y(itmp(1)-1,itmp(2))
       dytmp(1,1,1)= dy(itmp(1)-1,itmp(2),1)
       dytmp(1,1,2)= dy(itmp(1)-1,itmp(2),2)
       dytmp(1,1,3)= dy(itmp(1)-1,itmp(2),3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1)-1,itmp(2)+1,1)
       xtmp (1,2,2) = x(itmp(1)-1,itmp(2)+1,2)
       ytmp (1,2)   = y(itmp(1)-1,itmp(2)+1)
       dytmp(1,2,1)= dy(itmp(1)-1,itmp(2)+1,1)
       dytmp(1,2,2)= dy(itmp(1)-1,itmp(2)+1,2)
       dytmp(1,2,3)= dy(itmp(1)-1,itmp(2)+1,3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1),itmp(2)+1,1)
       xtmp (2,2,2) = x(itmp(1),itmp(2)+1,2)
       ytmp (2,2)   = y(itmp(1),itmp(2)+1)
       dytmp(2,2,1)= dy(itmp(1),itmp(2)+1,1)
       dytmp(2,2,2)= dy(itmp(1),itmp(2)+1,2)
       dytmp(2,2,3)= dy(itmp(1),itmp(2)+1,3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1),itmp(2),1)
       xtmp (2,1,2) = x(itmp(1),itmp(2),2)
       ytmp (2,1)   = y(itmp(1),itmp(2))
       dytmp(2,1,1)= dy(itmp(1),itmp(2),1)
       dytmp(2,1,2)= dy(itmp(1),itmp(2),2)
       dytmp(2,1,3)= dy(itmp(1),itmp(2),3)

       return

    else if (x(itmp(1),itmp(2),1) > xint(1) .and. &
         x(itmp(1),itmp(2),2) > xint(2)) then
       !the interpolated pt is in the the upper right corner

       !lower left corner
       xtmp(1,1,1)  = x(itmp(1)-1,itmp(2)-1,1)
       xtmp(1,1,2)  = x(itmp(1)-1,itmp(2)-1,2)
       ytmp(1,1)    = y(itmp(1)-1,itmp(2)-1)
       dytmp(1,1,1)= dy(itmp(1)-1,itmp(2)-1,1)
       dytmp(1,1,2)= dy(itmp(1)-1,itmp(2)-1,2)
       dytmp(1,1,3)= dy(itmp(1)-1,itmp(2)-1,3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1)-1,itmp(2),1)
       xtmp (1,2,2) = x(itmp(1)-1,itmp(2),2)
       ytmp (1,2)   = y(itmp(1)-1,itmp(2))
       dytmp(1,2,1)= dy(itmp(1)-1,itmp(2),1)
       dytmp(1,2,2)= dy(itmp(1)-1,itmp(2),2)
       dytmp(1,2,3)= dy(itmp(1)-1,itmp(2),3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1),itmp(2),1)
       xtmp (2,2,2) = x(itmp(1),itmp(2),2)
       ytmp (2,2)   = y(itmp(1),itmp(2))
       dytmp(2,2,1)= dy(itmp(1),itmp(2),1)
       dytmp(2,2,2)= dy(itmp(1),itmp(2),2)
       dytmp(2,2,3)= dy(itmp(1),itmp(2),3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1),itmp(2)-1,1)
       xtmp (2,1,2) = x(itmp(1),itmp(2)-1,2)
       ytmp (2,1)   = y(itmp(1),itmp(2)-1)
       dytmp(2,1,1)= dy(itmp(1),itmp(2)-1,1)
       dytmp(2,1,2)= dy(itmp(1),itmp(2)-1,2)
       dytmp(2,1,3)= dy(itmp(1),itmp(2)-1,3)

       return

    else if (x(itmp(1),itmp(2),1) == xint(1) .and. &
         x(itmp(1),itmp(2),2) < xint(2)) then

       !the interpolated pt is on the left edge
       ! towards the bottom

       !lower left corner
       xtmp(1,1,1)  = x(itmp(1)-1,itmp(2),1)
       xtmp(1,1,2)  = x(itmp(1)-1,itmp(2),2)
       ytmp(1,1)    = y(itmp(1)-1,itmp(2))
       dytmp(1,1,1)= dy(itmp(1)-1,itmp(2),1)
       dytmp(1,1,2)= dy(itmp(1)-1,itmp(2),2)
       dytmp(1,1,3)= dy(itmp(1)-1,itmp(2),3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1)-1,itmp(2)+1,1)
       xtmp (1,2,2) = x(itmp(1)-1,itmp(2)+1,2)
       ytmp (1,2)   = y(itmp(1)-1,itmp(2)+1)
       dytmp(1,2,1)= dy(itmp(1)-1,itmp(2)+1,1)
       dytmp(1,2,2)= dy(itmp(1)-1,itmp(2)+1,2)
       dytmp(1,2,3)= dy(itmp(1)-1,itmp(2)+1,3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1),itmp(2)+1,1)
       xtmp (2,2,2) = x(itmp(1),itmp(2)+1,2)
       ytmp (2,2)   = y(itmp(1),itmp(2)+1)
       dytmp(2,2,1)= dy(itmp(1),itmp(2)+1,1)
       dytmp(2,2,2)= dy(itmp(1),itmp(2)+1,2)
       dytmp(2,2,3)= dy(itmp(1),itmp(2)+1,3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1),itmp(2),1)
       xtmp (2,1,2) = x(itmp(1),itmp(2),2)
       ytmp (2,1)   = y(itmp(1),itmp(2))
       dytmp(2,1,1)= dy(itmp(1),itmp(2),1)
       dytmp(2,1,2)= dy(itmp(1),itmp(2),2)
       dytmp(2,1,3)= dy(itmp(1),itmp(2),3)

       return

    else if (x(itmp(1),itmp(2),1) == xint(1) .and. &
         x(itmp(1),itmp(2),2) > xint(2)) then
       !the interpolated pt is on the left edge
       ! towards the top

       !lower left corner
       xtmp(1,1,1)  = x(itmp(1)-1,itmp(2)-1,1)
       xtmp(1,1,2)  = x(itmp(1)-1,itmp(2)-1,2)
       ytmp(1,1)    = y(itmp(1)-1,itmp(2)-1)
       dytmp(1,1,1)= dy(itmp(1)-1,itmp(2)-1,1)
       dytmp(1,1,2)= dy(itmp(1)-1,itmp(2)-1,2)
       dytmp(1,1,3)= dy(itmp(1)-1,itmp(2)-1,3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1)-1,itmp(2),1)
       xtmp (1,2,2) = x(itmp(1)-1,itmp(2),2)
       ytmp (1,2)   = y(itmp(1)-1,itmp(2))
       dytmp(1,2,1)= dy(itmp(1)-1,itmp(2),1)
       dytmp(1,2,2)= dy(itmp(1)-1,itmp(2),2)
       dytmp(1,2,3)= dy(itmp(1)-1,itmp(2),3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1),itmp(2),1)
       xtmp (2,2,2) = x(itmp(1),itmp(2),2)
       ytmp (2,2)   = y(itmp(1),itmp(2))
       dytmp(2,2,1)= dy(itmp(1),itmp(2),1)
       dytmp(2,2,2)= dy(itmp(1),itmp(2),2)
       dytmp(2,2,3)= dy(itmp(1),itmp(2),3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1),itmp(2)-1,1)
       xtmp (2,1,2) = x(itmp(1),itmp(2)-1,2)
       ytmp (2,1)   = y(itmp(1),itmp(2)-1)
       dytmp(2,1,1)= dy(itmp(1),itmp(2)-1,1)
       dytmp(2,1,2)= dy(itmp(1),itmp(2)-1,2)
       dytmp(2,1,3)= dy(itmp(1),itmp(2)-1,3)

       return

    else if (x(itmp(1),itmp(2),1) < xint(1) .and. &
         x(itmp(1),itmp(2),2) == xint(2)) then

       !the interpolated pt is on the top edge
       ! towards the left

       !lower left corner
       xtmp(1,1,1)  = x(itmp(1),itmp(2)-1,1)
       xtmp(1,1,2)  = x(itmp(1),itmp(2)-1,2)
       ytmp(1,1)    = y(itmp(1),itmp(2)-1)
       dytmp(1,1,1)= dy(itmp(1),itmp(2)-1,1)
       dytmp(1,1,2)= dy(itmp(1),itmp(2)-1,2)
       dytmp(1,1,3)= dy(itmp(1),itmp(2)-1,3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1),itmp(2),1)
       xtmp (1,2,2) = x(itmp(1),itmp(2),2)
       ytmp (1,2)   = y(itmp(1),itmp(2))
       dytmp(1,2,1)= dy(itmp(1),itmp(2),1)
       dytmp(1,2,2)= dy(itmp(1),itmp(2),2)
       dytmp(1,2,3)= dy(itmp(1),itmp(2),3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1)+1,itmp(2),1)
       xtmp (2,2,2) = x(itmp(1)+1,itmp(2),2)
       ytmp (2,2)   = y(itmp(1)+1,itmp(2))
       dytmp(2,2,1)= dy(itmp(1)+1,itmp(2),1)
       dytmp(2,2,2)= dy(itmp(1)+1,itmp(2),2)
       dytmp(2,2,3)= dy(itmp(1)+1,itmp(2),3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1)+1,itmp(2)-1,1)
       xtmp (2,1,2) = x(itmp(1)+1,itmp(2)-1,2)
       ytmp (2,1)   = y(itmp(1)+1,itmp(2)-1)
       dytmp(2,1,1)= dy(itmp(1)+1,itmp(2)-1,1)
       dytmp(2,1,2)= dy(itmp(1)+1,itmp(2)-1,2)
       dytmp(2,1,3)= dy(itmp(1)+1,itmp(2)-1,3)

       return

    else if (x(itmp(1),itmp(2),1) > xint(1) .and. &
         x(itmp(1),itmp(2),2) == xint(2)) then

       !the interpolated pt is on the top
       ! towards the rigth

       !lower left corner
       xtmp(1,1,1)  = x(itmp(1)-1,itmp(2)-1,1)
       xtmp(1,1,2)  = x(itmp(1)-1,itmp(2)-1,2)
       ytmp(1,1)    = y(itmp(1)-1,itmp(2)-1)
       dytmp(1,1,1)= dy(itmp(1)-1,itmp(2)-1,1)
       dytmp(1,1,2)= dy(itmp(1)-1,itmp(2)-1,2)
       dytmp(1,1,3)= dy(itmp(1)-1,itmp(2)-1,3)

       !upper left corner
       xtmp (1,2,1) = x(itmp(1)-1,itmp(2),1)
       xtmp (1,2,2) = x(itmp(1)-1,itmp(2),2)
       ytmp (1,2)   = y(itmp(1)-1,itmp(2))
       dytmp(1,2,1)= dy(itmp(1)-1,itmp(2),1)
       dytmp(1,2,2)= dy(itmp(1)-1,itmp(2),2)
       dytmp(1,2,3)= dy(itmp(1)-1,itmp(2),3)

       !upper right corner
       xtmp (2,2,1) = x(itmp(1),itmp(2),1)
       xtmp (2,2,2) = x(itmp(1),itmp(2),2)
       ytmp (2,2)   = y(itmp(1),itmp(2))
       dytmp(2,2,1)= dy(itmp(1),itmp(2),1)
       dytmp(2,2,2)= dy(itmp(1),itmp(2),2)
       dytmp(2,2,3)= dy(itmp(1),itmp(2),3)

       !lower right corner
       xtmp (2,1,1) = x(itmp(1),itmp(2)-1,1)
       xtmp (2,1,2) = x(itmp(1),itmp(2)-1,2)
       ytmp (2,1)   = y(itmp(1),itmp(2)-1)
       dytmp(2,1,1)= dy(itmp(1),itmp(2)-1,1)
       dytmp(2,1,2)= dy(itmp(1),itmp(2)-1,2)
       dytmp(2,1,3)= dy(itmp(1),itmp(2)-1,3)

       return
    end if



  end subroutine get_4neighbors

end module interpol
