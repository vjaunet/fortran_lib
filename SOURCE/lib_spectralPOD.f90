module lib_spectralpod
  !=================Specification=============================
  !
  !
  ! contains POD routines
  !*  author : Vincent Jaunet
  !*  License: GPL v3.0
  !
  !
  !
  !===========================================================

  integer, private :: i,j,ic,jc,is,imod,if
  type SpectralPOD
     real(kind=4),      dimension(:,:)      ,allocatable  ::lambda
     complex(kind=4),   dimension(:,:,:)    ,allocatable  ::phi
     real(kind=4),      dimension(:,:,:)    ,allocatable  ::phiEn
     complex(kind=4),   dimension(:,:,:,:)  ,allocatable  ::upod

   contains
     procedure :: makepod_spectral_scalar
     procedure :: makepod_spectral_vec
     procedure :: makepod_spectral_scalar_f
     procedure :: makepod_spectral_vec_f
     generic   :: make => makepod_spectral_scalar,makepod_spectral_vec, &
          makepod_spectral_scalar_f,makepod_spectral_vec_f
     procedure :: recons => recons_spectral
  end type SpectralPOD

contains
  !#####################################################################################
  !#
  !#          Spectral POD analysis
  !#
  !#               -- Did on radial profiles
  !#               -- Cross-component cross spectrum is taken into acount
  !#
  !#####################################################################################

  !################### Scalar POD with frequency #########################
  subroutine makepod_spectral_scalar_f(datapod,r,u,f,param_psd)
    use lib_spectral
    class(SpectralPOD)                             ::datapod
    class(*)         ,dimension(:,:)               ::u
    real(kind=4)         ,dimension(:)             ::r,f
    type(psd_param)                                ::param_psd
    !----------------------------------------------------------

    call makepod_spectral_scalar(datapod,r,u,param_psd)

    !fill in f
    do if=1,param_psd.nfft
       f(if) = param_psd.fe*&
            dble(if-param_psd.nfft/2-1)/dble(param_psd.nfft)
    end do

    return

  end subroutine makepod_spectral_scalar_f

  !################### Vector POD with frequency #########################
  subroutine makepod_spectral_vec_f(datapod,r,u,f,param_psd)
    use lib_spectral
    class(SpectralPOD)                             ::datapod
    class(*)         ,dimension(:,:,:)             ::u
    real(kind=4)         ,dimension(:)             ::r,f
    type(psd_param)                                ::param_psd
    !----------------------------------------------------------

    call makepod_spectral_vec(datapod,r,u,param_psd)

    !fill in f
    do if=1,param_psd.nfft
       f(if) = param_psd.fe*&
            dble(if-param_psd.nfft/2-1)/dble(param_psd.nfft)
    end do

    return

  end subroutine makepod_spectral_vec_f

  !################### Scalar POD #########################
  subroutine makepod_spectral_scalar(datapod,r,u,param_psd)
    use lib_spectral
    implicit none

    class(SpectralPOD)                             ::datapod
    class(*)         ,dimension(:,:)               ::u
    real(kind=4)         ,dimension(:)             ::r
    type(psd_param)                                ::param_psd

    real(kind=4),    dimension(:,:),   allocatable ::ur4
    complex(kind=8), dimension(:,:),   allocatable ::uc8
    complex(kind=8), dimension(:,:,:), allocatable ::crosspsd
    real(kind=8)   , dimension(:)    , allocatable ::eigenval
    integer(kind=4)                                ::nc,nr,nf

    !for lapack
    INTEGER :: LWORK,INFO
    COMPLEX*16, DIMENSION(:) , allocatable         ::WORK
    REAL*8,     DIMENSION(:),  allocatable         ::RWORK

    logical                                        ::cplx=.false.
    !-------------------------------------------------------------
    !CALL ZHEEV('V','U',nr,crosspsd(:,:,i),nr, eigenval, WORK, LWORK, RWORK, INFO )

    nr=size(u,1)

    select type (u)
    type is (real*4)
       cplx=.false.
       allocate(ur4(size(u,1), size(u,2)))
       ur4 = u
       nf=param_psd.nfft/2+1
    type is (complex)
       cplx=.true.
       allocate(uc8(size(u,1), size(u,2)))
       uc8 = u
       nf=param_psd.nfft
    end select

    allocate(datapod.lambda(nr,nf))
    allocate(eigenval(nr))
    allocate(datapod.phi(nr,nr,nf))
    allocate(datapod.phiEn(nr,nr,nf))

    ! computing the cross_spectral matrix crossPSD
    !   --> XPSD(sqrt(r1)*U(r1),sqrt(r2)*U(r2))
    allocate(crosspsd(nr,nr,nf))
    crosspsd = 0.d0

    !fill-up the cross-component cross-spectral matrix portion
    if (cplx) then !u is complex
       do i=1,nr
          do j=i,nr
             call xpsd(dcmplx(uc8(i,:)*sqrt(r(i))), dcmplx(uc8(j,:)*sqrt(r(j))),&
                  crosspsd(i,j,:), param_psd)
          end do
       end do
    else ! u was real
       do i=1,nr
          do j=i,nr
             call xpsd(dble(ur4(i,:)*sqrt(r(i))), dble(ur4(j,:)*sqrt(r(j))),&
                  crosspsd(i,j,:), param_psd)
          end do
       end do
    end if

    ! making the spectral POD
    !   --> computes eigenvectors of the XPSD matrix at each frequency
    allocate(WORK(2*nr-1))
    allocate(RWORK(3*nr-2))
    LWORK = 2*nr-1
    do if=1,nf
       CALL ZHEEV('V','U',nr,crosspsd(:,:,if),nr, eigenval, WORK, LWORK, RWORK, INFO )
       IF(INFO.ne.0)THEN
          WRITE(*,*)'c_makepod_spectral : ZHEEV failed on exit'
       ENDIF

       ! do i=1,nr
       !    write(10,'(64e15.8,2x)')(real(crosspsd(i,j,if)),j=1,nr)
       ! end do

       !radial POD
       !Gamard et al JFM 2004
       do imod=1,nr
          !sorting lambda/phi in decreasing order
          datapod.lambda(imod,if) = eigenval(nr-imod+1)

          datapod.phi(2:nr,imod,if) = sqrt(datapod.lambda(imod,if))* &
               crosspsd(2:nr,nr-imod+1,if)/sqrt(r(2:nr))
          datapod.phiEn(:,imod,if) = abs(datapod.phi(:,imod,if))**2
       end do
    end do

    ! free some memory
    if (allocated(uc8)) deallocate(uc8)
    if (allocated(ur4)) deallocate(ur4)

  end subroutine makepod_spectral_scalar

  !################# Vector POD #########################
  subroutine makepod_spectral_vec(datapod,r,u,param_psd)
    use lib_spectral
    implicit none

    class(SpectralPOD)                             ::datapod
    class(*)         ,dimension(:,:,:)             ::u
    real(kind=4)         ,dimension(:)             ::r
    type(psd_param)                                ::param_psd

    real(kind=4),    dimension(:,:,:), allocatable ::ur4
    complex(kind=8), dimension(:,:,:), allocatable ::uc8
    complex(kind=8), dimension(:,:,:), allocatable ::crosspsd
    real(kind=8)   , dimension(:)    , allocatable ::eigenval
    integer(kind=4)                                ::nc,nr,nf

    !for lapack
    INTEGER :: LWORK,INFO
    COMPLEX*16, DIMENSION(:) , allocatable         ::WORK
    REAL*8,     DIMENSION(:),  allocatable         ::RWORK

    logical                                        ::cplx=.false.
    !-------------------------------------------------------------

    nr=size(u,1)
    nc=size(u,3)

    select type (u)
    type is (real*4)
       cplx=.false.
       allocate(ur4(size(u,1), size(u,2), size(u,3)))
       ur4 = u
       nf=param_psd.nfft/2+1
    type is (complex*8)
       cplx=.true.
       allocate(uc8(size(u,1), size(u,2), size(u,3)))
       uc8 = u
       nf=param_psd.nfft
    end select

    allocate(datapod.lambda(nc*nr,nf))
    allocate(eigenval(nc*nr))
    allocate(datapod.phi(nc*nr,nc*nr,nf))
    allocate(datapod.phiEn(nc*nr,nc*nr,nf))

    ! computing the cross_spectral matrix crossPSD
    !   --> XPSD(sqrt(r1)*U(r1),sqrt(r2)*U(r2))
    allocate(crosspsd(nc*nr,nc*nr,nf))
    crosspsd = 0.d0

    !fill-up the cross-component cross-spectral matrix portion
    if (cplx) then !u is complex
       do ic=1,nc
          do jc=ic,nc
             do i=1,nr
                do j=i,nr
                   call xpsd(uc8(i,ic,:)*sqrt(r(i)), uc8(j,jc,:)*sqrt(r(j)),&
                        crosspsd((ic-1)*nr + i, (jc-1)*nr + j,:), param_psd)
                end do
             end do
          end do
       end do
    else ! u was real
       do ic=1,nc
          do jc=ic,nc
             do i=1,nr
                do j=i,nr
                   call xpsd(dble(ur4(i,ic,:)*sqrt(r(i))), dble(ur4(j,jc,:)*sqrt(r(j))),&
                        crosspsd((ic-1)*nr + i, (jc-1)*nr + j,:), param_psd)
                end do
             end do
          end do
       end do
    end if

    ! making the spectral POD
    !   --> computes eigenvectors of the XPSD matrix at each frequency
    allocate(WORK(2*nc*nr-1))
    allocate(RWORK(3*nc*nr-2))
    LWORK = 2*nc*nr-1
    do i=1,nf
       CALL ZHEEV('V','U',nr,crosspsd(:,:,i),nr, eigenval, WORK, LWORK, RWORK, INFO )
       IF(INFO.ne.0)THEN
          WRITE(*,*)'c_makepod_spectral : ZHEEV failed on exit'
       ENDIF

       !radial POD
       !Gamard et al JFM 2004
       do imod=1,nc*nr
          !sorting lambda/phi in decreasing order
          datapod.lambda(imod,i) = eigenval(nr-imod+1)
          datapod.phi(:,imod,i) = sqrt(datapod.lambda(imod,i))* &
               crosspsd(:,nr-imod+1,i)/sqrt(r(:))
          datapod.phiEn(:,imod,i) = abs(datapod.phi(:,imod,i))**2
       end do
    end do

    ! free some memory
    if (allocated(uc8)) deallocate(uc8)
    if (allocated(ur4)) deallocate(ur4)

  end subroutine makepod_spectral_vec

  subroutine recons_spectral(datapod)
    implicit none
    class(SpectralPOD)                         ::datapod
    !-------------------------------------------------------


  end subroutine recons_spectral

end module lib_spectralpod
