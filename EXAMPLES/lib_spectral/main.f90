PROGRAM lib_spectral_usage
  use lib_spectral
  implicit none

  !*==================================================================
  !*
  !*
  !*         Spectral library Usage
  !*
  !*
  !*       author : Vincent Jaunet
  !*       date   : 11-06-2014
  !*       license: GPL v3.0
  !*       contact: vincent[dot]jaunet[at]ensma[dot]fr
  !*
  !*
  !======================================================================

  character(len=200)                               ::ifname,arg
  character(len=200)                               ::ofname
  integer                                          ::iarg=1

  integer, parameter                               ::Ns=1e5
  real   , parameter                               ::Fs=200000.
  real(kind=4)    ,dimension(:)  ,    allocatable  ::sig_r4
  real(kind=8)    ,dimension(:)  ,    allocatable  ::sig_r8
  complex(kind=8) ,dimension(:)  ,    allocatable  ::sig_c8

  type(psd_param)                                  ::parampsd

  complex(kind=8) ,dimension(:)  ,    allocatable  ::spectre,cpsd
  real(kind=8)    ,dimension(:)  ,    allocatable  ::f,coh
  real(kind=8)    ,dimension(:)  ,    allocatable  ::tau,R
  integer                                          ::nfft=1024

  real,parameter                                   ::pi=4.*atan(1.0)
  integer                                          ::i,ic
  !=======================================================


  do while (iarg <= command_argument_count())
     call get_command_argument(iarg, arg)

     select case (trim(arg))

     case("-i")
        iarg=iarg+1
        call get_command_argument(iarg, ifname)

     case("-o")
        iarg=iarg+1
        call get_command_argument(iarg, ofname)

     case("-nfft")
        iarg=iarg+1
        call get_command_argument(iarg, arg)
        read(arg,*)nfft

     case default
        write(06,*) "wrong input parameter :", arg
        call print_help()

     end select

     iarg = iarg+1
  end do

  !===================================================
  !* create the data
  !===================================================
  allocate(sig_r4(Ns)) !single precision
  allocate(sig_r8(Ns)) !double precision
  allocate(sig_c8(Ns)) !complex double precision

  do i=1,Ns
     sig_r4(i) = rand() + 10*cos(2*pi*1000*(i-1)/fs)
  end do
  sig_r8 = dble(sig_r4)
  sig_c8 = dcmplx(sig_r8,sig_r8)

  !===================================================
  !* compute the some spectral values
  !===================================================
  parampsd = psd_param(&
       nfft=nfft,&
       overlap=nfft/2,&
       window='A',&
       fe=fs,&
       check_pval=.false. &
       )

  !----- real signals
  allocate(f(nfft/2+1))

  !PSD
  allocate(spectre(nfft/2+1))
  call psd(sig_r4(:),f,spectre(:),parampsd)
  call psd(sig_r8(:),f,spectre(:),parampsd)


  !XPSD : note that I did not implement the r4 version
  allocate(cpsd(nfft/2+1))
  call xpsd(dble(sig_r4(:)),dble(sig_r4(:)),f,cpsd(:),parampsd)
  call xpsd(sig_r8(:),sig_r8(:),f,cpsd(:),parampsd)

  !Coherence : note that I did not implement the r4 version
  allocate(coh(nfft/2+1))
  call mscohere(dble(sig_r4(:)),dble(sig_r4(:)),f,coh(:),parampsd)
  call mscohere(sig_r8(:),sig_r8(:),f,coh(:),parampsd)

  deallocate(spectre,cpsd,coh,f)

  !------ complex signals
  allocate(f(nfft))

  !PSD
  allocate(spectre(nfft))
  call psd(sig_c8(:),f,spectre(:),parampsd)

  !XPSD :
  allocate(cpsd(nfft)) !note the nfft not the nfft/2+1
  call xpsd(sig_c8(:),sig_c8(:),f,cpsd(:),parampsd)

  !Coherence
  allocate(coh(nfft))
  call mscohere(sig_c8(:),sig_c8(:),f,coh(:),parampsd)

  deallocate(spectre,cpsd,coh,f)

  !===================================================
  !* compute the some correlations
  !===================================================
  parampsd = psd_param(&
       nfft=nfft,&
       overlap=nfft/2,&
       window='A',&
       fe=fs,&
       check_pval=.false. &
       )
  allocate(tau(nfft),R(nfft))

  !Auto-corelation using Wiener-Kintchine : r4 not implemented
  call cor(dble(sig_r4(:)),tau,R,parampsd)
  call cor(sig_r8(:),tau,R,parampsd)

  !Cross-corelation using Wiener-Kintchine
  call xcor(sig_r4(:),sig_r4(:),tau,R,parampsd)
  call xcor(sig_r8(:),sig_r8(:),tau,R,parampsd)

  deallocate(tau,R)

  !===================================================
  !* output the results
  !===================================================
  ! open(unit=11, file=trim(ofname), status='replace')
  ! do i=1,nfft/2+1
  !    write(11,'(10(e15.5,2x))')f(i),(abs(spectre(i,ic)),ic=1,nc)
  ! end do
  ! write(11,'(a)')" "
  ! close(11)


contains
  subroutine print_help()

    STOP 'This is the help, please fill me in'

  end subroutine print_help

END PROGRAM LIB_SPECTRAL_USAGE
