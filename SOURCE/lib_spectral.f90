module lib_spectral
  implicit none
  include 'fftw3'
  !*==================================================================
  !*
  !*
  !*         Spectral library for Signal processing
  !*
  !*
  !*       author : Vincent Jaunet
  !*       date   : 11-06-2014
  !*
  !*       contact: v.jaunet@gmail.com
  !*-------------------------------------------------------------------
  !*
  !* Requires: fftw3 to be intalled and linked during compilation
  !* TO DO : Coherence and Bi-Coherence
  !*
  !* Contains :
  !*
  !* fft(s(:),sp(:),psd_param)        , for fft 1d
  !* fft(s(:),f(:),sp(:),psd_param)   , for fft 1d w frequency
  !* ifft(sp(:),s(:),psd_param)       , for inverse fft 1d
  !*
  !* psd(s(:),sp(:),psd_param,check_Pval)      , for 1d PSD welch
  !* psd(s(:),f(:),sp(:),psd_param,check_Pval) , for 1d PSD welch w frequency
  !*
  !* psd(t(:),s(:),f(:),sp(:),psd_param,check_Pval) , for 1d PSD stlotting w frequency
  !
  !* cor(s(:),tau(:),cor(:),psd_param)       , for 1d auto-correlation
  !* cor(t(:),s(:),tau(:),cor(:),psd_param)  , for 1d auto-correlation slotting
  !*
  !* xpsd(s1(:),s2(:),f(:),xsp(:),psd_param)   ,for 1d cross PSD
  !* xcor(s1(:),s2(:),xsp(:),tau(:),psd_param) ,for 1d cross-correlation
  !*
  !* xpsd(t1(:),s1(:),t2(:),s2(:),f(:),xsp(:),psd_param)   ,for 1d X-PSD slotting
  !* xcor(t1(:),s1(:),t2(:),s2(:),tau(:),xcor(:),psd_param) ,for 1d X-cor slotting
  !*
  !*==============================================================================

  integer(kind=8), private ::i,j,k,in,if,ic

  type psd_param
     integer(kind=8)                ::nfft=1024
     integer(kind=8)                ::overlap=512
     character(len=1)               ::window='H'
     !available window type : box = 'B', Hamming = 'H', Hanning 'A'
     real(kind=8)                   ::fe=1.d0      !sampling freq
     real(kind=8)                   ::fmin=0.5d0   !min frequency to be computed (slotting)
     logical                        ::allocated_fft =.false.
     logical                        ::allocated_ifft=.false.
     logical                        ::norm_fft=.true.
     logical                        ::rms_norm=.false.
     logical                        ::check_pval=.false.
     logical                        ::detrend=.false.
     integer(kind=8)                ::plan=0
     integer(kind=8)                ::plan_ifft=0
  end type psd_param

  public  :: fft,ifft,&
       psd,cor,&
       xpsd,xcor

  private :: d_fft_1d,d_fft_1d_f,&
       d_psd_1d, d_psd_1d_f, d_cor_1d,&
       d_xpsd_1d, d_xcor_1d,&
       d_xpsd_1d_f,&
       d_cor_lda, d_psd_lda,&
       d_xcor_lda, d_xpsd_lda,&
       free_fft, rmlintrend,&
       triangle, slottingFuzzy

  interface fft
     module procedure d_fft_1d,d_fft_1d_f,d_fft_2d
  end interface fft

  interface ifft
     module procedure d_ifft_1d,d_ifft_2d
  end interface ifft

  interface psd
     module procedure d_psd_1d, d_psd_1d_f, d_psd_lda
  end interface psd

  interface cor
     module procedure d_cor_1d, d_cor_lda
  end interface cor

  interface xpsd
     module procedure d_xpsd_1d, d_xpsd_1d_f, d_xpsd_lda
  end interface xpsd

  interface xcor
     module procedure d_xcor_1d, d_xcor_lda, d_xcor_2d
  end interface xcor

contains

  !----------------------------------------------------
  !*    destroy plans
  !----------------------------------------------------

  subroutine free_fft(param)
    type(psd_param)                             ::param
    !----------------------------------------------------

    if (param.allocated_fft) then
       !free fftw
       call dfftw_destroy_plan(param.plan)

    end if

  end subroutine free_fft

  subroutine free_ifft(param)
    type(psd_param)                             ::param
    !----------------------------------------------------

    if (param.allocated_ifft) then
       !free ifftw
       call dfftw_destroy_plan(param.plan_ifft)

    end if

  end subroutine free_ifft


  !----------------------------------------------------
  !*    Fourier Transform
  !----------------------------------------------------

  subroutine d_fft_1d(s, sp, param)
    real(kind=8)     ,dimension(:)              ::s
    complex(kind=8)  ,dimension(:)              ::sp
    type(psd_param) ,optional                   ::param

    type(psd_param)                             ::def_param
    integer(kind=8)                             ::plan
    !----------------------------------------------------

    if (size(sp,1) /= size(s)/2+1) then
       STOP 'FFT : size(sp) must be size(s)/2+1.'
    end if

    if (present(param)) then
       def_param = param
    end if

    if (.not.def_param.allocated_fft) then
       !allocate fftw
       call dfftw_plan_dft_r2c_1d(def_param.plan,&
            def_param.nfft,s,sp,FFTW_ESTIMATE)
       !set allocated to true for next call
       def_param.allocated_fft = .true.
    end if

    !compute fftw
    call dfftw_execute(def_param.plan,s,sp)

    !normalization
    if (def_param.norm_fft) then
       sp = sp/def_param.nfft
    end if

    !return parameter values for next call
    if (present(param)) then
       param = def_param
    end if

    return

  end subroutine d_fft_1d

  subroutine d_fft_1d_f(s, f, sp ,param)
    real(kind=8)     ,dimension(:)              ::s
    complex(kind=8)  ,dimension(:)              ::sp
    type(psd_param) ,optional                   ::param
    real(kind=8)  ,dimension(:)                 ::f

    type(psd_param)                             ::def_param
    !----------------------------------------------------


    if (present(param)) then
       def_param = param
    end if

    !call fft
    call d_fft_1d(s, sp, param)

    !fill in F
    do if=1,def_param.nfft/2+1
       f(if) = def_param.fe*&
            dble(if-1)/dble(def_param.nfft)
    end do

    !return parameter values for next call
    if (present(param)) then
       param = def_param
    end if

    return

  end subroutine d_fft_1d_f


  subroutine d_fft_2d(s, sp, param)
    real(kind=8)    ,dimension(:,:)             ::s
    complex(kind=8) ,dimension(:,:)             ::sp
    type(psd_param) ,optional                   ::param

    type(psd_param)                             ::def_param
    integer(kind=8)                             ::nx,ny
    !----------------------------------------------------

    if (size(sp,1) /= size(s,1) .or. &
         size(sp,2) /= size(s,2)) then
       STOP 'size(sp) must be size(s) in d_dfft_2d'
    end if

    nx=size(s,1)
    ny=size(s,2)

    if (present(param)) then
       def_param = param
    end if

    if (.not.def_param.allocated_fft) then
       !allocate fftw
       call dfftw_plan_dft_2d(def_param.plan,&
            nx,ny,sp,sp,FFTW_FORWARD,FFTW_ESTIMATE)
       !set allocated to true for next call
       def_param.allocated_fft = .true.
    end if

    !compute fftw
    !passing S into FFTW
    sp = s
    call dfftw_execute_dft(def_param.plan,sp,sp)

    !normalization
    if (def_param.norm_fft) then
       sp = sp/dble(nx*ny)
    end if

    !return parameter values for next call
    if (present(param)) then
       param = def_param
    end if

    return

  end subroutine d_fft_2d

  !----------------------------------------------------
  !*   END Fourier Transform
  !----------------------------------------------------

  !----------------------------------------------------
  !*   Inverse Fourier Transform
  !----------------------------------------------------

  subroutine d_ifft_1d(sp,s,param)
    real(kind=8)     ,dimension(:)              ::s
    complex(kind=8)  ,dimension(:)              ::sp
    type(psd_param) ,optional                   ::param

    type(psd_param)                             ::def_param
    integer(kind=8)                             ::plan
    !----------------------------------------------------

    if (present(param)) then
       def_param = param
    end if

    if (.not.def_param.allocated_ifft) then
       !allocate fftw
       call dfftw_plan_dft_c2r_2d(def_param.plan_ifft,&
            size(sp,1),sp,s,FFTW_ESTIMATE)
       def_param.allocated_ifft = .true.
    end if

    !compute fftw
    call dfftw_execute(def_param.plan_ifft,sp,s)

    !return parameter values for next call
    if (present(param)) then
       param = def_param
    end if

    return

  end subroutine d_ifft_1d

  subroutine d_ifft_2d(sp,s,param)
    complex(kind=8)  ,dimension(:,:)            ::sp
    real(kind=8)     ,dimension(:,:)            ::s
    type(psd_param)  ,optional                  ::param

    type(psd_param)                             ::def_param
    integer(kind=8)                             ::nx,ny
    !----------------------------------------------------

    if (present(param)) then
       def_param = param
    end if

    nx=size(sp,1)
    ny=size(sp,2)

    if (.not.def_param.allocated_ifft) then
       !allocate fftw
       call dfftw_plan_dft_2d(def_param.plan_ifft,&
            nx,ny,sp,sp,FFTW_BACKWARD,FFTW_ESTIMATE)
       def_param.allocated_ifft = .true.
    end if

    !compute fftw
    call dfftw_execute_dft(def_param.plan_ifft,sp,sp)

    !get real value back
    s = dble(sp)

    !return parameter values for next call
    if (present(param)) then
       param = def_param
    end if

    return

  end subroutine d_ifft_2d

  !----------------------------------------------------
  !*   END Inverse Fourier Transform
  !----------------------------------------------------

  !----------------------------------------------------
  !*    Power spectral density
  !----------------------------------------------------

  subroutine d_psd_1d (s,sp,param)
    real(kind=8)     ,dimension(:)              ::s
    complex(kind=8)  ,dimension(:)              ::sp
    type(psd_param) ,optional                   ::param

    real(kind=8)   ,dimension(:) ,allocatable   ::xk,window
    complex(kind=8),dimension(:) ,allocatable   ::sk
    integer(kind=8)                             ::plan
    integer(kind=8)                             ::nf
    type(psd_param)                             ::def_param

    real(kind=8)                                ::moy,rms,sigmoy,sigrms
    real(kind=8)                                ::powfen2

    integer(kind=8)                             ::is_deb,is_fin
    integer(kind=8)                             ::nn
    !---------------------------------------------------

    nn = size(s,1)

    if (present(param)) then
       def_param = param
    end if

    if (size(sp) .ne. def_param.nfft/2+1) then
       write(06,*) 'sp size invalid : sp(1:nfft/2+1)'
       stop
    end if

    !remove mean and compute input rms
    sigmoy = sum(s(:))/dble(nn)
    sigrms = dsqrt(sum((s(:)-sigmoy)**2)/dble(nn))
    !remove the mean (DC)
    s(:)  = s(:) - sigmoy

    !initialization
    allocate(window(def_param.nfft))
    call get_window(def_param,window,powfen2)

    allocate(xk(def_param.nfft))
    allocate(sk(def_param.nfft/2+1))
    sp = 0.d0
    sk = 0.d0
    ic = 0
    is_deb = 1
    is_fin = def_param.nfft

    !start the psd processing
    do while(is_fin.le.nn)
       !Loop on the diffrent input signals,
       !simultneously acquired
       do i=is_deb,is_fin
          xk(i-ic*def_param.overlap) = s(i)
       end do

       ! !removing trend
       if (def_param.detrend) then
          call rmlintrend(xk,def_param.nfft)
       end if

       !windowing
       xk(:) = xk(:) * window(:)

       !fft
       call d_fft_1d(xk,sk,def_param)

       sp(:) = sp(:) + sk(:)*dconjg(sk)

       !Get new block
       is_deb = is_deb+def_param.overlap
       is_fin = is_fin+def_param.overlap
       ic     = ic + 1
    end do

    !free fftw3
    call free_fft(def_param)

    deallocate(xk,sk,window)

    !account of window energy
    !and one sided spectrum
    sp = 2.d0*sp/powfen2

    !Averaging the blocks
    sp = sp/dble(ic)

    !energy per hertz
    sp = sp/def_param.fe*dble(def_param.nfft)

    !check the parseval theorem if wanted
    if (def_param.check_pval == .true.) then
       rms = 0.d0
       rms = rms + 0.5d0*abs(sp(1))*&
            def_param.fe/dble(def_param.nfft)
       do if=2,def_param.nfft/2
          rms = rms + 1.d0*abs(sp(if))*&
               def_param.fe/dble(def_param.nfft)
       end do
       rms = rms + 0.5d0*abs(sp(def_param.nfft/2+1))*&
            def_param.fe/dble(def_param.nfft)
       write(06,'(a,e15.3,2x,e15.3)')'Parseval rms, sum(psd*df) : ',sigrms,sqrt(rms)
    end if

    !normalize output power
    if (def_param.rms_norm .and. sigrms /= 0.d0) then
       sp = sp/(sigrms**2)
    end if

    !Recovering original signal
    s = s + sigmoy

    return

  end subroutine d_psd_1d

  subroutine d_psd_1d_f (s,f,sp,param)
    real(kind=8)     ,dimension(:)              ::s
    complex(kind=8)  ,dimension(:)              ::sp
    type(psd_param)  ,optional                  ::param
    real(kind=8)     ,dimension(:)              ::f

    type(psd_param)                             ::def_param
    !---------------------------------------------------

    if (present(param)) then
       def_param = param
    end if

    !call psd
    call d_psd_1d(s,sp,def_param)

    !fill in f
    do if=1,def_param.nfft/2+1
       f(if) = def_param.fe*&
            dble(if-1)/dble(def_param.nfft)
    end do

    return

  end subroutine d_psd_1d_f

  !----------------------------------------------------
  !*   END PSD
  !----------------------------------------------------

  !----------------------------------------------------
  !*   Cross-PSD
  !----------------------------------------------------


  subroutine d_xpsd_1d (s1,s2,sp,param)
    real(kind=8)     ,dimension(:)              ::s1,s2
    complex(kind=8)  ,dimension(:)              ::sp
    type(psd_param)  ,optional                  ::param

    real(kind=8)   ,dimension(:) ,allocatable   ::xk,window
    complex(kind=8),dimension(:) ,allocatable   ::sk,sk_tmp
    integer(kind=8)                             ::plan
    integer(kind=8)                             ::nf
    type(psd_param)                             ::def_param

    real(kind=8)                                ::sigmoy2,sigrms2
    real(kind=8)                                ::sigmoy1,sigrms1
    real(kind=8)                                ::powfen2,rms

    integer(kind=8)                             ::is_deb,is_fin,nn
    !---------------------------------------------------

    nn = size(s1,1)

    if (present(param)) then
       def_param = param
    end if

    if (size(sp) .ne. def_param.nfft/2+1) then
       write(06,*) 'sp size invalid : sp(1:nfft/2+1) in d_xpsd_1d'
       stop
    end if

    !remove mean and compute input rms
    sigmoy1 = sum(s1(:))/dble(nn)
    sigrms1 = dsqrt(sum((s1(:)-sigmoy1)**2)/dble(nn))
    !remove the mean (DC)
    s1(:)  = s1(:) - sigmoy1

    sigmoy2 = sum(s2(:))/dble(nn)
    sigrms2 = dsqrt(sum((s2(:)-sigmoy2)**2)/dble(nn))
    !remove the mean (DC)
    s2(:)  = s2(:) - sigmoy2

    !initialization
    allocate(window(def_param.nfft))
    allocate(xk(def_param.nfft))
    allocate(sk(def_param.nfft/2+1))
    allocate(sk_tmp(def_param.nfft/2+1))
    sp = 0.d0
    sk = 0.d0
    sk_tmp = 0.d0
    ic = 0
    is_deb = 1
    is_fin = def_param.nfft

    call get_window(def_param,window,powfen2)

    !start the psd processing
    do while(is_fin.le.nn)

       !extract sample
       do i=is_deb,is_fin
          xk(i-ic*def_param.overlap) = s1(i)
       end do

       !removing trend s1
       !call rmlintrend(xk,def_param.nfft)

       !windowing
       xk(:) = xk(:) * window(:)

       !fft
       call d_fft_1d(xk,sk,def_param)
       sk_tmp = sk

       !extract sample
       do i=is_deb,is_fin
          xk(i-ic*def_param.overlap) = s2(i)
       end do

       !removing trend s2
       !call rmlintrend(xk,def_param.nfft)

       !windowing
       xk(:) = xk(:) * window(:)

       !fft
       call d_fft_1d(xk,sk,def_param)

       !compute inter-spectrum
       sp = sp + sk_tmp*dconjg(sk)

       !Get new block
       is_deb = is_deb+def_param.overlap
       is_fin = is_fin+def_param.overlap
       ic     = ic + 1
    end do

    deallocate(xk)
    deallocate(sk,sk_tmp)
    deallocate(window)

    !free fftw3
    call free_fft(def_param)

    !account of window energy
    !and one sided spectrum
    sp = 2.d0*sp/powfen2

    !Averaging the blocks
    sp = sp/dble(ic)

    !energy per hertz
    sp = sp/def_param.fe*dble(def_param.nfft)

    if (def_param.check_pval) then
       rms = 0.d0
       rms = rms + 0.5d0*abs(sp(1))*&
            def_param.fe/dble(def_param.nfft)
       do if=2,def_param.nfft/2
          rms = rms + 1.d0*abs(sp(if))*&
               def_param.fe/dble(def_param.nfft)
       end do
       rms = rms + 0.5d0*abs(sp(def_param.nfft/2+1))*&
            def_param.fe/dble(def_param.nfft)
       write(06,'(a,e15.3,2x,e15.3)')'Parseval rms**2, sum(xpsd*df) : ', sigrms1*sigrms2,rms
    end if

    !normalize output power
    if (def_param.rms_norm &
         .and. sigrms1 /= 0.d0 &
         .and. sigrms2 /= 0.d0) then
       sp = sp/(sigrms1*sigrms2)
    end if

    return

  end subroutine d_xpsd_1d

  subroutine d_xpsd_1d_f (s1,s2,f,sp,param)
    real(kind=8)     ,dimension(:)              ::s1,s2,f
    complex(kind=8)  ,dimension(:)              ::sp
    type(psd_param)  ,optional                  ::param

    type(psd_param)                             ::def_param
    !---------------------------------------------------

    if (present(param)) then
       def_param = param
    end if

    !call xpsd
    call d_xpsd_1d(s1,s2,sp,def_param)

    !free fftw3
    call free_fft(def_param)

    !fill in f
    do if=1,def_param.nfft/2+1
       f(if) = def_param.fe*&
            dble(if-1)/dble(def_param.nfft)
    end do

    return

  end subroutine d_xpsd_1d_f

  !----------------------------------------------------
  !*  END Cross-PSD
  !----------------------------------------------------

  !----------------------------------------------------
  !*  PSD via Slotting Techniques for LDA
  !----------------------------------------------------

  subroutine d_psd_lda(at,s,f,sp,param)
    !computes the PSD of an unevenly sampled signal
    !using Fuzzy slotting technique and Wiener-Kintchine

    type(psd_param)  ,optional                  ::param
    real(kind=8)     ,dimension(:)              ::at,s
    real(kind=8)     ,dimension(:)              ::f
    complex(kind=8)  ,dimension(:)              ::sp

    real(kind=8)  ,dimension(:)  ,allocatable   ::window
    real(kind=8)                                ::powfen2

    real(kind=8)                                ::rms
    type(psd_param)                             ::def_param
    real(kind=8), dimension(:) ,allocatable     ::tau,xcor
    integer(kind=8)                             ::nn
    !---------------------------------------------------

    nn = size(s)

    if (present(param)) then
       def_param = param
    end if

    !fill in frequencies
    def_param.fe = def_param.fmin*real(def_param.nfft)/2
    do if=1,def_param.nfft/2+1
       f(if) = def_param.fe*&
            dble(if-1)/dble(def_param.nfft)
    end do

    !compute slotting
    allocate(tau(def_param.nfft),xcor(def_param.nfft))
    call d_cor_lda(at,s,tau,xcor,def_param)

    !windowing
    allocate(window(def_param.nfft))
    call get_window(def_param,window,powfen2,.true.)
    xcor(:) = xcor(:) * window(:)

    !compute Fourrier Transform of Correlation
    call d_fft_1d(xcor,sp,def_param)

    !free fft
    call free_fft(def_param)

    !one-sided spectrum
    sp = 2.d0*sp/powfen2

    !Account for window energy
    sp = sp/powfen2

    !energy per hertz
    sp = sp/def_param.fe*dble(def_param.nfft)

    !check the parseval theorem if wanted
    if (def_param.check_pval) then
       rms = 0.d0
       rms = rms + 0.5d0*abs(sp(1))*&
            def_param.fe/dble(def_param.nfft)
       do if=2,def_param.nfft/2
          rms = rms + 1.d0*abs(sp(if))*&
               def_param.fe/dble(def_param.nfft)
       end do
       rms = rms + 0.5d0*abs(sp(def_param.nfft/2+1))*&
            def_param.fe/dble(def_param.nfft)
       write(06,'(a,f15.3,2x,f15.3)')'Parseval rms**2/sum(psd*df) : ',&
            xcor(def_param.nfft/2+1),rms
    end if

    !deallocation of tables
    deallocate(tau,xcor)

    !set fft_parameters for next call
    if (present(param)) then
       param.fe = def_param.fe
    end if

    deallocate(window)

    return

  end subroutine d_psd_lda

  subroutine d_xpsd_lda(at1,s1,at2,s2,f,sp,param)
    !computes the cross-PSD of an unevenly sampled signal
    !using Fuzzy slotting technique and Wiener-Kintchine

    type(psd_param)  ,optional                  ::param
    real(kind=8)     ,dimension(:)              ::at1,s1
    real(kind=8)     ,dimension(:)              ::at2,s2
    real(kind=8)     ,dimension(:)              ::f
    complex(kind=8)  ,dimension(:)              ::sp

    real(kind=8)  ,dimension(:)  ,allocatable   ::window
    real(kind=8)                                ::powfen2
    real(kind=8)                                ::rms
    integer(kind=8)                             ::nn1,nn2

    type(psd_param)                             ::def_param
    real(kind=8), dimension(:) ,allocatable     ::tau,xcor
    !---------------------------------------------------

    nn1 = size(s1)
    nn2 = size(s2)

    if (present(param)) then
       def_param = param
    end if

    !fill in frequencies
    def_param.fe = def_param.fmin*real(def_param.nfft)/2
    do if=1,def_param.nfft/2+1
       f(if) = def_param.fe*&
            dble(if-1)/dble(def_param.nfft)
    end do

    !compute slotting
    allocate(tau(def_param.nfft),xcor(def_param.nfft))
    call d_xcor_lda(at1,s1,at2,s2,tau,xcor,def_param)

    !windowing
    allocate(window(def_param.nfft))
    call get_window(def_param,window,powfen2,.true.)
    xcor(:) = xcor(:) * window(:)

    !compute Fourrier Transform of Correlation
    call d_fft_1d(xcor,sp,def_param)

    !free ftt
    call free_fft(def_param)

    !one-sided spectrum
    sp = 2.d0*sp/powfen2

    !Account for window energy
    sp = sp/powfen2

    !energy per hertz
    sp = sp/def_param.fe*dble(def_param.nfft)

    !check the parseval theorem if wanted
    if (def_param.check_pval) then
       rms = 0.d0
       rms = rms + 0.5d0*abs(sp(1))*&
            def_param.fe/dble(def_param.nfft)
       do if=2,def_param.nfft/2
          rms = rms + 1.d0*abs(sp(if))*&
               def_param.fe/dble(def_param.nfft)
       end do
       rms = rms + 0.5d0*abs(sp(def_param.nfft/2+1))*&
            def_param.fe/dble(def_param.nfft)
       write(06,'(a,f15.3,2x,f15.3)')'Parseval rms**2/sum(psd*df) : ',&
            xcor(def_param.nfft/2+1),rms
    end if

    !set fft_parameters for next call
    if (present(param)) then
       param.plan = def_param.plan
       param.allocated_fft = .true.
       param.fe = def_param.fe
    end if

    !deallocation of tables
    deallocate(tau,xcor)

  end subroutine d_xpsd_lda

  !----------------------------------------------------
  !* END PSD Slotting Techniques for LDA
  !----------------------------------------------------


  !----------------------------------------------------
  !*  Correlations Slotting and Cross-slotting
  !----------------------------------------------------

  subroutine d_cor_lda(at,s,tau,xcor,param)
    !      ..... Calculation of the auto-correlation function
    !      ..... using the slotting correlation process

    !      DECLARATION OF GLOBAL VARIABLES
    !      at         : arrival time of signal samples
    !      s          : velocity samples
    !      tau        : retarded time history of the correlation
    !      cor        : values history of the correlation

    implicit none

    type(psd_param) ,optional                     ::param

    !slotting parameter
    type(psd_param)                               ::def_param
    real(kind=8)                                  ::DTslot,tau_min,tau_max
    integer(kind=8)                               ::nn

    !signal variables
    real(kind=8),dimension(:)                     ::at,s
    real(kind=8)                                  ::Sm,Srms

    !output result
    real(kind=8), dimension(:)                    ::tau,xcor
    !--------------------------------------------------------------

    nn = size(s)

    if (present(param)) then
       def_param = param
    end if

    !Slotting the correlation function
    DTslot  = 1.d0/(def_param.fmin)*2.d0/real(def_param.nfft)
    tau_max = 1.d0/(def_param.fmin)
    tau_min = -tau_max

    !Time delay vector of the correlation
    do k=1,def_param.nfft
       tau(k) = tau_min + DTslot*real(k-1)
    enddo

    !Remove mean parts and norm_fft
    Sm   = sum(s(:))/real(nn)
    Srms = sqrt(sum((s(:)-Sm)**2)/real(nn))
    s(:) = (s(:) - Sm)/Srms

    !Norm_fftd Correlation calculation
    call slottingFuzzy(nn,at,s,&
         nn,at,s,&
         tau,DTslot,xcor,&
         def_param.nfft)

    !Normalization
    xcor = xcor * Srms**2

    !Recover Original signal
    s = s*Srms + Sm

    return
  end subroutine d_cor_lda


  subroutine d_xcor_lda(at1,s1,at2,s2,tau,xcor,param)
    !      ..... Calculation of the cross-correlation function
    !      ..... using the slotting correlation process

    !      DECLARATION OF GLOBAL VARIABLES
    !      at         : arrival time of signal samples
    !      s          : velocity samples
    !      tau        : retarded time history of the correlation
    !      cor        : values history of the correlation

    implicit none

    type(psd_param) ,optional                     ::param

    !slotting parameter
    type(psd_param)                               ::def_param
    real(kind=8)                                  ::DTslot,tau_min,tau_max

    !signal variables
    integer(kind=8)	                          ::nn1,nn2
    real(kind=8),dimension(:)                     ::at1,s1
    real(kind=8),dimension(:)                     ::at2,s2
    real(kind=8)                                  ::S1m,S1rms
    real(kind=8)                                  ::S2m,S2rms

    !output result
    real(kind=8), dimension(:)                    ::tau,xcor
    !--------------------------------------------------------------

    nn1 = size(s1)
    nn2 = size(s2)

    if (present(param)) then
       def_param = param
    end if

    !Slotting the correlation function
    DTslot  = 1.d0/(def_param.fmin)*2.d0/real(def_param.nfft)
    tau_max = 1.d0/(def_param.fmin)
    tau_min = -tau_max

    !Time delay vector of the correlation
    do k=1,def_param.nfft
       tau(k) = tau_min + DTslot*real(k-1)
    enddo

    !Remove mean parts and norm_fft
    S1m   = sum(s1(:))/real(nn1)
    S1rms = sqrt(sum((s1(:)-S1m)**2)/real(nn1))
    s1(:) = (s1(:) - S1m)/S1rms

    S2m   = sum(s2(:))/real(nn2)
    S2rms = sqrt(sum((s2(:)-S2m)**2)/real(nn2))
    s2(:) = (s2(:) - S2m)/S2rms

    !Norm_fftd Correlation calculation
    call slottingFuzzy(nn1,at1,s1,&
         nn2,at2,s2,&
         tau,DTslot,xcor,&
         def_param.nfft)

    !Normalization
    xcor = xcor*S1rms*S2rms

    !Recover Original signals
    s1 = s1*S1rms + S1m
    s2 = s2*S2rms + S2m


    return
  end subroutine d_xcor_lda


  !........................................................................
  !     Subroutine for Fuzzy Slotting with Local Normalisation
  !     M.J. Tummers and D.M. Passchier, 1996, Spectral estimation using
  !     a variable window and the slotting technique with local
  !     normalisation, Meas. Sc. Technol. 10, pp. L4-L7

  SUBROUTINE slottingFuzzy(n1,t1,u1 &
       ,n2,t2,u2 &
       ,dt,dtslot,cu,nfft)

    implicit none
    INTEGER(kind=8) 	                     ::n1,n2
    INTEGER(kind=8)	                     ::p,q,qL,k
    INTEGER(kind=8)	                     ::idx,bug
    INTEGER(kind=8)	                     ::nfft
    REAL(kind=8),DIMENSION(0:n1-1)             ::t1
    REAL(kind=8),DIMENSION(0:n2-1)             ::t2
    REAL(kind=8),DIMENSION(0:n2-1)             ::u1
    REAL(kind=8),DIMENSION(0:n2-1)             ::u2
    REAL(kind=8) 		                     ::lag,w
    REAL(kind=8),DIMENSION(1:nfft)             ::nu1
    REAL(kind=8),DIMENSION(1:nfft)             ::nu2
    REAL(kind=8),DIMENSION(1:nfft)             ::cu
    REAL(kind=8),DIMENSION(1:nfft)             ::dt
    REAL(kind=8)		                     ::dtslot

    integer(kind=8),parameter                  ::un=1,deux=2

    !	Initialisation
    cu(1:nfft) =0.d0
    nu1(1:nfft)=0.d0
    nu2(1:nfft)=0.d0

    qL=0
    DO p=0,n1-1
       !       Beginning of loop on time series of signal 1
       bug=0
       !       qL denotes the lower limit of signal 2 index under which it
       !       is not neccessary to compute the loop since the lag is out of
       !       range.
       q=qL
       DO while (bug.eq.0 .and. q.le.n2-1)
          !       Beginning of loop on time series of signal 2
          lag=t2(q)-t1(p)
          if (abs(lag).ge.dt(nfft)+dtslot) then
             !            lag is out of range
             if (lag.lt.0) then
                !               if lag negative, continue loop incrementing the time
                q=q+1
                qL=q
             else if (lag.gt.0 .or. q.ge.n2-1) then
                !               if lag positive, no need to continue with this
                !            ...time lag, exit loop and increment time signal 1
                bug=1
             endif
          else
             !            lag is in range
             !               slot number for the specific lag
             if (lag.lt.dt(1)) then
                !                lag falls into the lower limit of the first slot,
                !             ...contributes only to the first slot

                call triangle(dt(1)-dtslot,dt(1),lag,w,deux)
                cu(1)  = cu(1) + u1(p)*u2(q)*w
                nu1(1) = nu1(1)+ u1(p)**2*w
                nu2(1) = nu2(1)+ u2(q)**2*w
             else if (lag.gt.dt(nfft)) then
                !                lag falls into the upper limit of the last slot,
                !             ...contributes only to the last slot

                call triangle(dt(nfft),dt(nfft)+dtslot,lag,w,un)
                cu(nfft)  = cu(nfft) + u1(p)*u2(q)*w
                nu1(nfft) = nu1(nfft)+ u1(p)**2*w
                nu2(nfft) = nu2(nfft)+ u2(q)**2*w
             else
                !               in the other cases, the lag contributes to both
                !            ...the k-th and k+1-th slot
                idx = aint(nfft/2 + lag/dtslot + 1)
                if (idx.lt.1 .or. idx.gt.nfft) then
                   stop 'Error in slot index calculation, STOP !'
                else
                   call triangle(dt(idx),dt(idx+1),lag,w,un)
                   cu(idx)  =cu(idx)  + u1(p)*u2(q)*w
                   nu1(idx) =nu1(idx) + u1(p)**2*w
                   nu2(idx) =nu2(idx) + u2(q)**2*w

                   call triangle(dt(idx),dt(idx+1),lag,w,deux)
                   cu(idx+1)  = cu(idx+1) + u1(p)*u2(q)*w
                   nu1(idx+1) = nu1(idx+1) + u1(p)**2*w
                   nu2(idx+1) = nu2(idx+1) + u2(q)**2*w
                endif
             endif
             !               increment next time signal 2
             q=q+1
          endif
          !       End of loop on time series of signal 2
       ENDDO
       !       End of loop on time series of signql 1
    ENDDO

    !	Normalisation
    do p=1,nfft,1
       cu(p)=cu(p)/sqrt(nu1(p)*nu2(p))
    enddo

    RETURN
  END SUBROUTINE slottingFuzzy


  !------------------------------------------------------------------
  !      Subroutine for the calculation of the triangular weighting
  !------------------------------------------------------------------

  SUBROUTINE triangle(tmi,tpi,txi,coefti,op)

    implicit none
    REAL(kind=8)                          ::tmi,tpi,txi
    REAL(kind=8)                          ::coefti
    INTEGER(kind=8)                       ::op

    if (op.eq.1) then
       coefti=txi/(tmi-tpi)+tpi/(tpi-tmi)
    else if (op.eq.2) then
       coefti=txi/(tpi-tmi)+tmi/(tmi-tpi)
    else
       stop 'Error in triangle index'
    endif

    RETURN
  END SUBROUTINE triangle

  !----------------------------------------------------
  !* END Correaltion  Slotting Techniques for LDA
  !----------------------------------------------------


  !----------------------------------------------------
  !*    Correlation Wiener-Kintchine
  !----------------------------------------------------


  subroutine d_cor_1d(s,tau,cor,param)
    real(kind=8)     ,dimension(:)              ::s
    real(kind=8)     ,dimension(:)              ::cor
    real(kind=8)     ,dimension(:)              ::tau
    type(psd_param)  ,optional                  ::param

    complex(kind=8) ,dimension(:) , allocatable ::sp,sp_full
    real(kind=8)    ,dimension(:) , allocatable ::cor_tmp
    type(psd_param)                             ::def_param

    real(kind=8)                                ::rms
    integer(kind=8)                             ::nn
    !---------------------------------------------------

    nn = size(cor)

    if (size(tau) /= size(cor)) then
       STOP 'size(tau) /= size(cor) in d_cor_1d'
    end if

    if (present(param)) then
       def_param = param
    end if

    !fill in tau
    do i=1,nn
       tau(i) = (dble(i)-nn/2-1)/def_param.fe
    end do

    !compute PSD using Welch's method
    allocate(sp(nn/2+1))
    call d_psd_1d(s,sp,def_param)

    !compute inverse FFT of PSD/2.d0 (Wiener-Kintchine)
    allocate(sp_full(nn))
    sp_full(1:nn/2+1) = sp/2.
    sp_full(nn/2+2:nn) = sp(nn/2+2:2:-1)/2.
    call d_ifft_1d(sp_full,cor,def_param)

    !free ifft
    call free_ifft(def_param)

    !swap left an rigth
    allocate(cor_tmp(nn))
    cor_tmp(nn/2+1:nn) = cor(1:nn/2)
    cor_tmp(1:nn/2)    = cor(nn/2+1:nn)

    !norm_fft correlation
    cor = cor_tmp*def_param.fe/real(def_param.nfft)

    deallocate(cor_tmp)

    return

  end subroutine d_cor_1d

  !-----------------------

  subroutine d_xcor_1d(s1,s2,tau,xcor,param)
    real(kind=8)     ,dimension(:)              ::s1,s2
    real(kind=8)     ,dimension(:)              ::xcor
    type(psd_param)  ,optional                  ::param
    real(kind=8)     ,dimension(:)              ::tau

    complex(kind=8) ,dimension(:) , allocatable ::xsp,xsp_full
    real(kind=8)    ,dimension(:) , allocatable ::xcor_tmp
    type(psd_param)                             ::def_param

    real(kind=8)                                ::rms
    integer(kind=8)                             ::nf
    !---------------------------------------------------

    if (size(tau) /= size(xcor)) then
       STOP 'size(tau) /= size(cor) in d_xcor_1d'
    else if (size(s1) .lt. size(tau)) then
       STOP 'size(s1) < size(tau) in d_cor_1d'
    end if

    if (present(param)) then
       def_param = param
    end if

    nf = def_param.nfft

    !fill in tau
    do i=1,nf
       tau(i) = (dble(i)-nf/2-1)/def_param.fe
    end do

    !compute XPSD using Welch's method
    allocate(xsp(nf/2+1))
    call d_xpsd_1d(s1,s2,xsp,def_param)

    !compute inverse FFT of PSD/2.0 (Wiener-Kintchine)
    allocate(xsp_full(nf))
    xsp_full(1:nf/2+1)  = xsp/2.
    xsp_full(nf/2+2:nf) = xsp(nf/2+2:2:-1)/2.
    call d_ifft_1d(xsp_full,xcor,def_param)

    !free ifft
    call free_ifft(def_param)

    !swap left an rigth
    allocate(xcor_tmp(nf))
    xcor_tmp(nf/2+1:nf) = xcor(1:nf/2)
    xcor_tmp(1:nf/2)    = xcor(nf/2+1:nf)
    xcor = xcor_tmp

    !norm_fft the correlation
    xcor =xcor*def_param.fe/real(def_param.nfft)

    deallocate(xcor_tmp)
    deallocate(xsp)

    return

  end subroutine d_xcor_1d

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine d_xcor_2d(s1,s2,tau,xcor,param)
    real(kind=8)     ,dimension(:,:)            ::s1,s2
    real(kind=8)     ,dimension(:,:)            ::xcor
    type(psd_param)  ,optional                  ::param
    real(kind=8)     ,dimension(:,:,:)          ::tau

    complex(kind=8) ,dimension(:,:),allocatable ::xsp1,xsp2,xsp_full
    real(kind=8)    ,dimension(:,:),allocatable ::xcor_tmp_x,xcor_tmp_y
    type(psd_param)                             ::def_param

    real(kind=8)                                ::rms1,rms2
    real(kind=8)                                ::avg1,avg2
    integer(kind=8)                             ::nx,ny
    !---------------------------------------------------

    if (size(tau,1) /= size(xcor,1) .or. &
         size(tau,2) /= size(xcor,2)) then
       STOP 'size(tau) /= size(cor) in d_xcor_2d'
    else if (size(tau,1) /= size(s1,1) .or. &
         size(tau,2) /= size(s1,2)) then
       STOP 'size(tau) /= size(s1) in d_xcor_2d'
    else if (size(tau,1) /= size(s2,1) .or. &
         size(tau,2) /= size(s2,2)) then
       STOP 'size(tau) /= size(s2) in d_xcor_2d'
    end if

    nx = size(s1,1)
    ny = size(s1,2)

    if (present(param)) then
       def_param = param
       ! we handle nfft normalisation in here
       def_param.norm_fft = .false.
    end if

    !fill in tau
    do i=1,nx
       do j=1,ny
          tau(i,j,1) = (dble(i)-nx/2-1)/def_param.fe
          tau(i,j,2) = (dble(j)-ny/2-1)/def_param.fe
       end do
    end do

    !compute interspectrum
    allocate(xsp1(nx,ny))
    allocate(xsp2(nx,ny))
    call d_fft_2d(s1,xsp1,def_param)
    call d_fft_2d(s2,xsp2,def_param)

    !removing the mean value (set zero freq to 0)
    xsp1(1,1) = 0.d0
    xsp2(1,1) = 0.d0

    !compute inverse FFT
    allocate(xsp_full(nx,ny))
    xsp_full  = xsp1*dconjg(xsp2)
    call d_ifft_2d(xsp_full,xcor,def_param)

    !free ffts
    call free_ifft(def_param)
    call free_fft(def_param)

    !swap x
    allocate(xcor_tmp_x(nx/2,ny))
    xcor_tmp_x(:,:)   = xcor(1:nx/2,:)
    xcor(1:nx/2,:)    = xcor(nx/2+1:nx,:)
    xcor(nx/2+1:nx,:) = xcor_tmp_x(:,:)
    deallocate(xcor_tmp_x)

    !swap y
    allocate(xcor_tmp_y(nx,ny/2))
    xcor_tmp_y(:,:)   = xcor(:,1:ny/2)
    xcor(:,1:ny/2)    = xcor(:,ny/2+1:ny)
    xcor(:,ny/2+1:ny) = xcor_tmp_y(:,:)
    deallocate(xcor_tmp_y)


    !norm_fft the correlation
    xcor =xcor*def_param.fe/dble(nx*ny)**2

    !getthe correlation coefficient
    if (def_param.rms_norm) then
       avg1 = sum(s1)/nx/ny
       avg2 = sum(s2)/nx/ny

       rms1 = sum((s1-avg1)**2)/nx/ny
       rms2 = sum((s2-avg2)**2)/nx/ny

       if (rms1 /= 0.d0 .and. rms2 /= 0.d0) then
          xcor = xcor/sqrt(rms1*rms2)
       end if
    end if

    deallocate(xsp1,xsp2,xsp_full)

    return

  end subroutine d_xcor_2d


  !----------------------------------------------------
  !*  End Wienner-Kintchine
  !----------------------------------------------------

  !----------------------------------------------------
  !*  Windowing
  !----------------------------------------------------

  subroutine get_window(type,window,power,Wiener_K)
    type(psd_param)                             ::type
    real(kind=8) ,dimension(:)                  ::window
    complex(kind=8) ,dimension(:),allocatable   ::FTwindow
    real(kind=8)                                ::power
    logical             ,optional               ::Wiener_K
    integer(kind=8)                             ::plan,nn
    real(kind=8)                                ::pi = 4.d0*atan(1.d0)
    !--------------------------------------------------------------

    nn = size(window)
    allocate(FTwindow(nn/2+1))

    select case(type.window)

    case('B')

       window = 1.d0

    case('H')

       do i=1,nn
          window(i) = 0.54d0 - 0.46d0*cos(2.d0*pi*dble(i-1)/dble(nn-1))
       end do

    case('A')

       do i=1,nn
          window(i) = 0.5d0*(1.d0 - cos(2.d0*pi*dble(i-1)/dble(nn-1)))
       end do

    end select

    !allocate FFT
    call dfftw_plan_dft_r2c_1d(&
         plan,nn,window,FTwindow,FFTW_ESTIMATE)

    !execute and scaling FFT
    call dfftw_execute(plan)
    FTwindow = FTwindow/dble(type.nfft)

    !free fft
    call dfftw_destroy_plan(plan)

    !for Wiener_Kintchine method
    if (.not.present(Wiener_K)) then

       !calculate power of the window function
       power = 0.d0
       power = power + 1.0d0*abs(FTwindow(1))**2
       do i=2,nn/2
          power = power + 2.d0*abs(FTwindow(i))**2
       end do
       power = power + 1.0d0*abs(FTwindow(nn/2+1))**2

    else
       if (Wiener_K) then
          !calculate power of the window function
          power = 0.d0
          power = power + 1.0d0*abs(FTwindow(1))
          do i=2,nn/2
             power = power + 2.d0*abs(FTwindow(i))
          end do
          power = power + 1.0d0*abs(FTwindow(nn/2+1))

       end if
    end if

    deallocate(FTwindow)

  end subroutine get_window

  !----------------------------------------------------
  !*  End Windowing
  !----------------------------------------------------

  !----------------------------------------------------
  !*  Remove linear trend
  !----------------------------------------------------

  subroutine rmlintrend(yy,yxdim)
    implicit none

    !calculates linear regression by lest-square-method

    !input varaibles
    integer(kind=8)                                    ::yxdim

    !intern variables
    integer(kind=8)                                    ::i
    real(kind=8),  dimension(:), allocatable           ::xx
    real(kind=8)                                       ::sx,xm
    real(kind=8)                                       ::sxy,ym
    real(kind=8)                                       ::aaa,bbb

    !input-output variables
    real(kind=8),  dimension(yxdim)                    ::yy

    !---------------------------------------------------------------
    allocate(xx(yxdim))
    do i=1,yxdim
       xx(i) = dble(i-1)
    end do

    sx  = 0.d0
    ym  = 0.d0
    xm  = 0.d0
    sxy = 0.d0

    do i=1,yxdim
       xm  = xm  + xx(i)
       ym  = ym  + yy(i)
    end do
    xm = xm/dble(yxdim)
    ym = ym/dble(yxdim)

    do i=1,yxdim
       sx = sx + (xx(i)-xm)**2
       sxy = sxy + (xx(i)-xm)*(yy(i)-ym)
    end do
    sx  = sx/dble(yxdim)
    sxy = sxy/dble(yxdim)

    aaa = Sxy/Sx
    bbb = ym - xm*Sxy/Sx

    !remove tendency
    do i=1,yxdim
       yy(i) = yy(i) - (aaa*dble(i-1) + bbb)
    end do

  end subroutine rmlintrend


  !----------------------------------------------------
  !*  END Remove linear trend
  !----------------------------------------------------


end module lib_spectral
