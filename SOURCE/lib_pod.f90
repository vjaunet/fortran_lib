module lib_pod
  implicit none

  !=================Specification=============================
  !
  !
  ! contains POD routines
  !
  !*  author : Vincent Jaunet
  !*  License: GPL v3.0
  !
  !
  !
  !===========================================================

  integer, private :: i,j,ic,jc,is,imod

  type PODdata
     real(kind=4),   dimension(:)        ,allocatable  ::lambda
     real(kind=4),   dimension(:,:,:,:)  ,allocatable  ::phi
     real(kind=4),   dimension(:,:)      ,allocatable  ::alpha
     real(kind=4),   dimension(:,:,:,:)  ,allocatable  ::upod

   contains
     procedure :: calpod => f_makepod
     procedure :: recons => f_recons
  end type PODdata

  private :: f_gappypod
  interface gappypod
     module procedure :: f_gappypod
  end interface

contains

  !#####################################################################################
  !#
  !#          Temporal POD analysis : Snapshot POD (Sirovich,'87)
  !#
  !#####################################################################################

  subroutine f_makepod(datapod,u)

    !$$
    !$$ This subroutine does POD Decompostion
    !$$
    !$$ JAUNET Vincent, 10/08/2015
    !$$
    !$$
    !$$ ref : - An application of Gappy POD - Exp_ Fluids (2007) 42:79-91
    !$$       - The Karhunene-loeve procedure for gappy data J. Opt. Soc. Am A 12(8):1657-1664
    !$$       - Sirovich 1987
    !$$ needs : LAPACK library
    !$$

    implicit none
    !------------------------------------
    !declarations
    !____________________________________
    class(PODdata)                                           ::datapod

    !input variables
    real(kind=4)         ,dimension(:,:,:,:)                  ::u

    !calculation variables
    !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    integer(kind=4)                                           ::iVmx,jvmx,nc
    integer(kind=4)                                           ::nsetv
    integer(kind=4)                                           ::nmodv

    integer(kind=4)                                           ::k,l
    integer(kind=4)                                           ::iimg,imod
    real(kind=4)                                              ::err,err1,err2
    integer(kind=4)                                           ::cpt,cpt1

    real(kind=4),   dimension(:,:),   allocatable             ::cormat
    integer(kind=4)                                           ::nwork
    integer(kind=4)                                           ::ifail
    real(kind=4),   dimension(:),     allocatable             ::work
    real(kind=4),   dimension(:,:),   allocatable             ::cormat_old
    real(kind=4),   dimension(:),     allocatable             ::lambda_old

    real(kind=4)                                              ::t1,t2
    !-----------------------------------------------------------------------------------

    ivmx = size(u(:,1,1,1))
    jvmx = size(u(1,:,1,1))
    nc   = size(u(1,1,:,1))
    nsetv = size(u(1,1,1,:))
    nmodv = nsetv

    !allocation of the arrays
<<<<<<< HEAD
    allocate(datapod.phi(ivmx,jvmx,nc,nsetv))
    allocate(datapod.alpha(nsetv,nsetv))
    allocate(datapod.lambda(nsetv))
=======
    allocate(datapod%phi(ivmx,jvmx,nc,nsetv))
    allocate(datapod%alpha(nsetv,nsetv))
    allocate(datapod%lambda(nsetv))
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
    allocate(cormat(nsetv,nsetv))

    !initialization of some variables
    nwork       = 64*nsetV

    !___________________
    ! POD procedure

    !calculation of the correlation matrix
    write(06,*)'Computing correlation matrix'
    cormat = 0.d0
    !$OMP PARALLEL DO PRIVATE(l,i,j,ic) SHARED(u,cormat,ivmx,nc,jvmx,nsetv) SCHEDULE(Dynamic)
    do k=1,nsetv
       do l=k,nsetv
          do i=1,ivmx
             do j=1,jvmx
                do ic=1,nc
                   cormat(k,l) = cormat(k,l) +&
                        u(i,j,ic,k)*u(i,j,ic,l)
                end do
             end do
          end do
          cormat(k,l) = cormat(k,l)/(real(ivmx*jvmx*nsetv))
       end do
    end do
    !$OMP END PARALLEL DO

    !----- don't have to fill up cormat entirely
    ! do k=1,nsetv
    !    do l=k,nsetv
    !       cormat(l,k) = cormat(k,l)
    !    end do
    ! end do

    !Solving the POD
    allocate(work(nwork))
    allocate(lambda_old(nsetv))
    allocate(cormat_old(nsetv,nsetv))

<<<<<<< HEAD
    datapod.lambda = 0.d0
    write(06,*)"Entering SSYEV..."
    call SSYEV('V','U',nsetv,cormat,nsetv,datapod.lambda,WORK,nwork,ifail)

    if (ifail .ne. 0) then
       write(06,*)'Error in POD calculation, ifail =',ifail
       write(06,*)'lambda :',datapod.lambda(ifail)
       stop
    end if

    lambda_old = datapod.lambda
    do i=1,nsetv
       datapod.lambda(i) = lambda_old(nsetv-i+1)
    end do
    cpt1 = 0
    do imod=1,nsetv
       if (datapod.lambda(imod) < 0.d0) then
          datapod.lambda(imod) = 1e-20
=======
    datapod%lambda = 0.d0
    write(06,*)"Entering SSYEV..."
    call SSYEV('V','U',nsetv,cormat,nsetv,datapod%lambda,WORK,nwork,ifail)

    if (ifail .ne. 0) then
       write(06,*)'Error in POD calculation, ifail =',ifail
       write(06,*)'lambda :',datapod%lambda(ifail)
       stop
    end if

    lambda_old = datapod%lambda
    do i=1,nsetv
       datapod%lambda(i) = lambda_old(nsetv-i+1)
    end do
    cpt1 = 0
    do imod=1,nsetv
       if (datapod%lambda(imod) < 0.d0) then
          datapod%lambda(imod) = 1e-20
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
          cpt1 = cpt1 + 1
       end if
    end do
    if (cpt1 .ne. 0) then
       print*, cpt1,'   lambda < 0.'
    end if

    cormat_old = cormat
    do i=1,nsetv
       do j=1,nsetv
          !norme at lambda
<<<<<<< HEAD
          cormat(i,j) = cormat_old(i,nsetv-j+1)*sqrt(real(nsetv)*datapod.lambda(j))
=======
          cormat(i,j) = cormat_old(i,nsetv-j+1)*sqrt(real(nsetv)*datapod%lambda(j))
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
       end do
    end do

    deallocate(lambda_old)
    deallocate(cormat_old)
    deallocate(work)

    !calculating the spatial basis functions
    write(06,*)'Computing basis functions'
<<<<<<< HEAD
    datapod.phi = 0.d0
=======
    datapod%phi = 0.d0
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
    !$OMP PARALLEL DO PRIVATE(i,j,ic,iimg) &
    !$OMP& SHARED(u,cormat,datapod,nc,ivmx,jvmx,nsetv) &
    !$OMP& SCHEDULE(Dynamic)
    do imod=1,nsetv
       do i=1,ivmx
          do j=1,jvmx
             do ic=1,nc
                do iimg=1,nsetv

<<<<<<< HEAD
                   datapod.phi(i,j,ic,imod) = datapod.phi(i,j,ic,imod) +&
                        u(i,j,ic,iimg)*cormat(iimg,imod)

                end do
                datapod.phi(i,j,ic,imod) = datapod.phi(i,j,ic,imod)/(real(nsetv)*&
                     datapod.lambda(imod))
=======
                   datapod%phi(i,j,ic,imod) = datapod%phi(i,j,ic,imod) +&
                        u(i,j,ic,iimg)*cormat(iimg,imod)

                end do
                datapod%phi(i,j,ic,imod) = datapod%phi(i,j,ic,imod)/(real(nsetv)*&
                     datapod%lambda(imod))
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !saving alphas
<<<<<<< HEAD
    datapod.alpha=cormat
=======
    datapod%alpha=cormat
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0

    !deallocation of the arrays
    deallocate(cormat)


  end subroutine f_makepod

  !##################################################################
  !#
  !#
  !##################################################################

  subroutine f_recons(datapod,nmodes)
    class(PODdata)                 ::datapod
    integer(kind=4)                ::nmodes

    integer(kind=4)                ::nx,ny,nc,ns
    integer(kind=4)                ::i,j,ic,iimg,imod
    !----------------------------------------------
<<<<<<< HEAD
    nx = size(datapod.phi(:,1,1,1))
    ny = size(datapod.phi(1,:,1,1))
    nc = size(datapod.phi(1,1,:,1))
    ns = size(datapod.phi(1,1,1,:))

    datapod.upod = 0.d0
=======
    nx = size(datapod%phi(:,1,1,1))
    ny = size(datapod%phi(1,:,1,1))
    nc = size(datapod%phi(1,1,:,1))
    ns = size(datapod%phi(1,1,1,:))

    datapod%upod = 0.d0
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0

    !$OMP PARALLEL DO PRIVATE(j,ic,iimg,imod) &
    !$OMP& SHARED(datapod,nc,nx,ny,ns) &
    !$OMP& SCHEDULE(Dynamic)
    do i=1,nx
       do j=1,ny
          do ic=1,nc
             do iimg=1,ns
                do imod=1,nmodes
<<<<<<< HEAD
                   datapod.upod(i,j,ic,iimg) = datapod.upod(i,j,ic,iimg) + &
                        datapod.alpha(iimg,imod)*datapod.phi(i,j,ic,imod)
=======
                   datapod%upod(i,j,ic,iimg) = datapod%upod(i,j,ic,iimg) + &
                        datapod%alpha(iimg,imod)*datapod%phi(i,j,ic,imod)
>>>>>>> f68c88dba04320fe3fd5b93242bdec0455f3b2f0
                end do
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine f_recons

  !##################################################################
  !#
  !#
  !##################################################################

  subroutine f_gappypod(u,w,cptmx,err_min)

    !$$
    !$$ This subroutine corrects the PIV images using POD Decompostion
    !$$ and gives the POd decomposition of this corrected image
    !$$
    !$$ JAUNET Vincent, 10/10/2010
    !$$        11-12-2012 : corrected bug on replacement by mean value
    !$$
    !$$
    !$$ DON'T FORGET to adjust the parameters, depends on your signal
    !$$ err_min : minimum error you want
    !$$ cptmx   : maximum number of loop
    !$$
    !$$ ref : An application of Gappy POD - Exp_ Fluids (2007) 42:79-91
    !$$       The Karhunene-loeve procedure for gappy data J. Opt. Soc. Am A 12(8):1657-1664
    !$$
    !$$ needs : LAPACK library
    !$$
    !$$ Input Variables  :
    !$$     cal_type    : depends on the calculation type wanted 'M' for mean value, 'F' for fluctuating ones
    !$$     cptmx       : integer(kind=8) maximum number of POD loops
    !$$     err_min     : real(kind=8)    minimum error on energy to end Gappy POD
    !$$     iVmx        : integer(kind=8) horizontal domain size
    !$$     jvmx        : integer(kind=8) vertical domain size
    !$$     nsetV       : integer(kind=8) number of snapshot
    !$$     nmodv       : integer(kind=8) number of modes needed out ouf the subroutine
    !$$                                   and used to reconstrcut the field
    !$$     vMask       : integer(kind=8),dimension(iVmx,jVmx,nsetV) contains the erroneous vector location
    !$$ Output Variables :
    !$$     phispatvu1c : real(kind=8),   dimension(iVmx,jVmx,nmodv) spatial modes of the first component
    !$$     phispatvu2c : real(kind=8),   dimension(iVmx,jVmx,nmodv) spatial modes of the first component
    !$$     lambdav     : real(kind=8),   dimension(nsetV) proper values (energie modes), increasing order
    !$$     vU1         : real(kind=8),   dimension(iVmx,jVmx,nsetV) first components of velocity
    !$$     vU2         : real(kind=8),   dimension(iVmx,jVmx,nsetV) second component of velocity


    implicit none
    !------------------------------------
    !declarations
    !____________________________________
    !input variables
    real(kind=4)         ,dimension(:,:,:,:)                  ::u
    real(kind=4)         ,dimension(:,:,:)                    ::w
    real(kind=4)                                              ::err_min
    integer(kind=4)                                           ::cptmx

    !calculation variables
    !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    integer(kind=4)                                           ::iVmx,jvmx,nc
    integer(kind=4)                                           ::nsetv
    integer(kind=4)                                           ::nmodv

    integer(kind=4)                                           ::k,l
    integer(kind=4)                                           ::iimg,imod
    real(kind=4)                                              ::err,err1,err2
    integer(kind=4)                                           ::cpt,cpt1

    real(kind=4),   dimension(:,:,:,:), allocatable             ::upod
    real(kind=4),   dimension(:,:),   allocatable             ::cormat
    real(kind=4),   dimension(:),     allocatable             ::lambda_err
    real(kind=4),   dimension(:,:,:,:), allocatable           ::phi

    integer(kind=4)                                           ::nwork
    integer(kind=4)                                           ::ifail
    real(kind=4),   dimension(:),     allocatable             ::work
    real(kind=4),   dimension(:,:),   allocatable             ::cormat_old
    real(kind=4),   dimension(:),     allocatable             ::lambda_old

    real(kind=4)                                              ::t1,t2

    !output variables
    real(kind=4),   dimension(:) , allocatable                ::lambdav
    !-----------------------------------------------------------------------------------

    print*,"here"

    ivmx = size(u(:,1,1,1))
    jvmx = size(u(1,:,1,1))
    nc   = size(u(1,1,:,1))
    nsetv = size(u(1,1,1,:))
    nmodv = nsetv

    !allocation of the arrays
    allocate(upod(ivmx,jvmx,nc,nsetv))
    allocate(cormat(nsetv,nsetv))
    allocate(lambda_err(nsetv))
    allocate(lambdav(nsetv))
    allocate(phi(ivmx,jvmx,nc,nsetv))

    !initialization of some variables
    nwork       = 64*nsetV
    lambda_err  = 0.d0
    err         = 1e15
    cpt = 1

    !------------------------------------
    !begining of the process
    !Set outliers to 0
    do i=1,ivmx
       do j=1,jvmx
          do ic=1,nc
             do iimg=1,nsetv
                u(i,j,ic,iimg) = u(i,j,ic,iimg)*w(i,j,iimg)
             end do
          end do
       end do
    enddo


    !___________________
    !gappy POD procedure

    do while((abs(err) >= err_min) .and. cpt <= cptmx)

       !calculation of the correlation matrix
       cormat = 0.d0
       do k=1,nsetv
          do l=k,nsetv
             do i=1,ivmx
                do j=1,jvmx
                   do ic=1,nc
                      cormat(k,l) = cormat(k,l) +&
                           u(i,j,ic,k)*u(i,j,ic,l)
                   end do
                end do
             end do
             cormat(k,l) = cormat(k,l)/(real(ivmx*jvmx*nsetv))
          end do
       end do

       ! don't have to be done...
       ! do k=1,nsetv
       !    do l=k,nsetv
       !       cormat(l,k) = cormat(k,l)
       !    end do
       ! end do

       !Solving the POD
       allocate(work(nwork))
       allocate(lambda_old(nsetv))
       allocate(cormat_old(nsetv,nsetv))

       lambdav = 0.d0
       write(06,*)"Entering SSYEV..."
       call SSYEV('V','U',nsetv,cormat,nsetv,lambdav,WORK,nwork,ifail)

       if (ifail .ne. 0) then
          write(06,*)'Error in POD calculation, ifail =',ifail
          write(06,*)'lambda :',lambdav(ifail)
          stop
       end if

       lambda_old = lambdav
       do i=1,nsetv
          lambdav(i) = lambda_old(nsetv-i+1)
       end do

       cpt1 = 0
       do imod=1,nsetv
          if (lambdav(imod) < 0.d0) then
             lambdav(imod) = 1e-20
             cpt1 = cpt1 + 1
          end if
       end do
       if (cpt1 .ne. 0) then
          print*, cpt1,'   lambda < 0.'
       end if

       cormat_old = cormat
       do i=1,nsetv
          do j=1,nsetv
             !norme at lambda
             cormat(i,j) = cormat_old(i,nsetv-j+1)*sqrt(real(nsetv)*lambdav(j))
          end do
       end do

       deallocate(lambda_old)
       deallocate(cormat_old)
       deallocate(work)

       !calculating the spatial basis functions
       phi = 0.d0
       do imod=1,nsetv
          do i=1,ivmx
             do j=1,jvmx
                do ic=1,nc
                   do iimg=1,nsetv

                      phi(i,j,ic,imod) = phi(i,j,ic,imod) +&
                           u(i,j,ic,iimg)*cormat(iimg,imod)

                   end do
                   phi(i,j,ic,imod) = phi(i,j,ic,imod)/(real(nsetv)*lambdav(imod))
                end do
             end do
          end do
       end do

       !repairing the velocity field
       write(06,*)"repairing the velocity"
       upod = 0.d0
       do i=1,ivmx
          do j=1,jvmx
             do ic=1,nc
                do iimg=1,nsetv
                   do imod=1,nmodv
                      upod(i,j,ic,iimg) = upod(i,j,ic,iimg) + &
                           cormat(iimg,imod)*phi(i,j,ic,imod)
                   end do
                end do
             end do
          end do
       end do


       do i=1,ivmx
          do j=1,jvmx
             do ic=1,nc
                do k=1,nsetV
                   u(i,j,ic,k) = u(i,j,ic,iimg)*w(i,j,iimg) +&
                        (1.0 - w(i,j,iimg))*upod(i,j,ic,iimg)
                end do
             end do
          end do
       end do

       !calulating the error
       err  = 0.d0
       err1 = 0.d0
       err2 = 0.d0
       do imod=1,nsetv
          err1 = err1 + (lambdav(imod) - lambda_err(imod))**2
          err2 = err2 + lambdav(imod)**2
       end do
       err = sqrt(err1)/sqrt(err2)
       print*,'error :',err

       lambda_err = lambdav
       cpt = cpt + 1

    end do

    !deallocation of the arrays
    deallocate(upod)
    deallocate(cormat)
    deallocate(lambda_err)
    deallocate(phi)


  end subroutine f_gappypod

end module lib_pod
