	  �f  �   k820309              13.0        ۗ'U                                                                                                           
       lib_spectral.f90 LIB_SPECTRAL              D_FFT_1D D_FFT_1D_F D_PSD_1D D_PSD_1D_F D_COR_1D D_XPSD_1D D_XCOR_1D D_XPSD_1D_F D_COR_LDA D_PSD_LDA D_XCOR_LDA D_XPSD_LDA FREE_FFT RMLINTREND TRIANGLE SLOTTINGFUZZY UNWRAP_PHASE_D UNWRAP_PHASE_F I J K IN IF IC                                                        u #D_FFT_1D    #D_FFT_1D_F    #D_FFT_2D    #         @     @X                                               #D_FFT_1D%PRESENT    #D_FFT_1D%SIZE    #S    #SP    #PARAM                                                    PRESENT                                                 SIZE           D@                                                 
               &                                                     D@                                                                &                                                     F @                                    P               #PSD_PARAM    #         @     @X                                                #D_FFT_1D_F%DBLE 	   #D_FFT_1D_F%PRESENT 
   #S    #F    #SP    #PARAM                                               	     DBLE                                            
     PRESENT           D @                                                 
               &                                                     D                                                   
               &                                                     D @                                                                &                                                     F @                                    P               #PSD_PARAM                                                           u #D_IFFT_1D    #D_IFFT_2D                                                           u #D_PSD_1D    #D_PSD_1D_F    #D_PSD_LDA %   #         @     @X                                               #D_PSD_1D%SQRT    #D_PSD_1D%ABS    #D_PSD_1D%DCONJG    #D_PSD_1D%DSQRT    #D_PSD_1D%SUM    #D_PSD_1D%DBLE    #D_PSD_1D%PRESENT    #D_PSD_1D%SIZE    #S    #SP    #PARAM                                                    SQRT                                                 ABS                                                 DCONJG                                                 DSQRT                                                 SUM                                                 DBLE                                                 PRESENT                                                 SIZE           D@                                                 
               &                                                     D@                                                                &                                                      @                                    P               #PSD_PARAM    #         @     @X                                                #D_PSD_1D_F%DBLE    #D_PSD_1D_F%PRESENT     #S !   #F "   #SP #   #PARAM $                                                   DBLE                                                  PRESENT           D @                              !                   
               &                                                     D                                "                   
               &                                                     D @                              #                                  &                                                      @                               $     P               #PSD_PARAM    #         @     @X                             %                   #D_PSD_LDA%ABS &   #D_PSD_LDA%DBLE '   #D_PSD_LDA%PRESENT (   #D_PSD_LDA%SIZE )   #D_PSD_LDA%REAL *   #AT +   #S ,   #F -   #SP .   #PARAM /                                              &     ABS                                            '     DBLE                                            (     PRESENT                                            )     SIZE                                            *     REAL           D @                              +                   
               &                                                     D@                              ,                   
                &                                                     D                                -                   
 !              &                                                     D @                              .                    "              &                                                     F @                               /     P               #PSD_PARAM                                                           u #D_COR_1D 0   #D_COR_LDA 9   #         @     @X                             0                   #D_COR_1D%DBLE 1   #D_COR_1D%PRESENT 2   #D_COR_1D%SIZE 3   #D_COR_1D%REAL 4   #S 5   #TAU 6   #COR 7   #PARAM 8                                              1     DBLE                                            2     PRESENT                                            3     SIZE                                            4     REAL           D @                              5                   
 A              &                                                     D@                              6                   
 C              &                                                     D@                              7                   
 B              &                                                      @                               8     P               #PSD_PARAM    #         @     @X                            9                   #D_COR_LDA%SQRT :   #D_COR_LDA%SUM ;   #D_COR_LDA%PRESENT <   #D_COR_LDA%SIZE =   #D_COR_LDA%REAL >   #AT ?   #S @   #TAU A   #XCOR B   #PARAM C                                              :     SQRT                                            ;     SUM                                            <     PRESENT                                            =     SIZE                                            >     REAL           D @                              ?                   
 /              &                                                     D@                              @                   
 0              &                                                     D @                              A                   
 1              &                                                     D @                              B                   
 2              &                                                      @                               C     P               #PSD_PARAM                                                           u #D_XPSD_1D D   #D_XPSD_1D_F Q   #D_XPSD_LDA Y   #         @     @X                            D                   #D_XPSD_1D%SQRT E   #D_XPSD_1D%ABS F   #D_XPSD_1D%DCONJG G   #D_XPSD_1D%DSQRT H   #D_XPSD_1D%SUM I   #D_XPSD_1D%DBLE J   #D_XPSD_1D%PRESENT K   #D_XPSD_1D%SIZE L   #S1 M   #S2 N   #SP O   #PARAM P                                              E     SQRT                                            F     ABS                                            G     DCONJG                                            H     DSQRT                                            I     SUM                                            J     DBLE                                            K     PRESENT                                            L     SIZE           D@                              M                   
               &                                                     D                                N                   
               &                                                     D@                              O                                  &                                                      @                               P     P               #PSD_PARAM    #         @     @X                             Q                   #D_XPSD_1D_F%DBLE R   #D_XPSD_1D_F%PRESENT S   #S1 T   #S2 U   #F V   #SP W   #PARAM X                                              R     DBLE                                            S     PRESENT           D @                              T                   
               &                                                     D @                              U                   
               &                                                     D                                V                   
               &                                                     D @                              W                                  &                                                      @                               X     P               #PSD_PARAM    #         @     @X                             Y                   #D_XPSD_LDA%ABS Z   #D_XPSD_LDA%DBLE [   #D_XPSD_LDA%PRESENT \   #D_XPSD_LDA%SIZE ]   #D_XPSD_LDA%REAL ^   #AT1 _   #S1 `   #AT2 a   #S2 b   #F c   #SP d   #PARAM e                                              Z     ABS                                            [     DBLE                                            \     PRESENT                                            ]     SIZE                                            ^     REAL           D @                              _                   
 &              &                                                     D@                              `                   
 '              &                                                     D @                              a                   
 (              &                                                     D@                              b                   
 )              &                                                     D                                c                   
 *              &                                                     D @                              d                    +              &                                                     F @                               e     P               #PSD_PARAM                                                           u #D_XCOR_1D f   #D_XCOR_LDA p   #D_XCOR_2D }   #         @     @X                             f                   #D_XCOR_1D%DBLE g   #D_XCOR_1D%PRESENT h   #D_XCOR_1D%SIZE i   #D_XCOR_1D%REAL j   #S1 k   #S2 l   #TAU m   #XCOR n   #PARAM o                                              g     DBLE                                            h     PRESENT                                            i     SIZE                                            j     REAL           D@                              k                   
 G              &                                                     D @                              l                   
 H              &                                                     D@                              m                   
 J              &                                                     D@                              n                   
 I              &                                                      @                               o     P               #PSD_PARAM    #         @     @X                            p                   #D_XCOR_LDA%SQRT q   #D_XCOR_LDA%SUM r   #D_XCOR_LDA%PRESENT s   #D_XCOR_LDA%SIZE t   #D_XCOR_LDA%REAL u   #AT1 v   #S1 w   #AT2 x   #S2 y   #TAU z   #XCOR {   #PARAM |                                              q     SQRT                                            r     SUM                                            s     PRESENT                                            t     SIZE                                            u     REAL           D @                              v                   
 3              &                                                     D@                              w                   
 4              &                                                     D @                              x                   
 5              &                                                     D@                              y                   
 6              &                                                     D @                              z                   
 7              &                                                     D @                              {                   
 8              &                                                      @                               |     P               #PSD_PARAM                                                           u #UNWRAP_PHASE_F ~   #UNWRAP_PHASE_D �   #         @     @X                             ~                   #UNWRAP_PHASE_F%SIGN    #UNWRAP_PHASE_F%ATAN �   #UNWRAP_PHASE_F%ABS �   #UNWRAP_PHASE_F%SIZE �   #PHI �                                                   SIGN                                             �     ATAN                                            �     ABS                                            �     SIZE           D@                              �                   	 \              &                                           #         @     @X                             �                   #UNWRAP_PHASE_D%SIGN �   #UNWRAP_PHASE_D%ATAN �   #UNWRAP_PHASE_D%ABS �   #UNWRAP_PHASE_D%SIZE �   #PHI �                                              �     SIGN                                             �     ATAN                                            �     ABS                                            �     SIZE           D@                              �                   
 [              &                                                                                        �                                                       0                                             �                                                      1                                             �                                                      2                                             �                                                      3                                             �                                                      4                                             �                                                      5                                             �                                                      6                                             �                                                      7                                             �                                                      8                                             �                                       	               9                                             �                                       
               10                                             �                                          ��������                                                     �                                                      1                                             �                                                       0                                             �                                                      1                                             �                                                      2                                             �                                                      4                                             �                                                      8                                             �                                                      16                                             �                                                       32                                             �                                       @               64                                             �                                                       2097152                                             �                                       �               128                                             �                                                      256                                             �                                                      512                                             �                                                      1024                                             �                                                      2048                                             �                                                      4096                                             �                                                       8192                                             �                                        @              16384                                             �                                        �              32768                                             �                                                      65536                                             �                                                      131072                                             �                                                      262144                                             �                                                      524288                                             �                                                      1048576                  @                               'P                    #NFFT �   #OVERLAP �   #WINDOW �   #FE �   #FMIN �   #ALLOCATED_FFT �   #ALLOCATED_IFFT �   #NORM_FFT �   #RMS_NORM �   #CHECK_PVAL �   #DETREND �   #PLAN �   #PLAN_IFFT �               �                              �                                                                                               1024                �                              �                                                                                              512                �                              �                                                                                                        CH                                �                              �              
                                                
                       �?        1.D0                �                              �               
                                                
                       �?        0.5D0                �                               �     (                                                                                                             �                               �     ,                                                                                                             �                               �     0                                                                             ��������                        �                               �     4       	                                                                                                      �                               �     8       
                                                                                                      �                               �     <                                                                                                             �                              �     @                                                                                          0                �                              �     H                                                                                          0    #         @      X                                               #D_FFT_2D%DBLE �   #D_FFT_2D%PRESENT �   #D_FFT_2D%SIZE �   #S �   #SP �   #PARAM �                                              �     DBLE                                            �     PRESENT                                            �     SIZE           @                              �                   
               &                   &                                                     D@                              �                                  &                   &                                                     F @                               �     P               #PSD_PARAM    #         @      X                                               #D_IFFT_1D%PRESENT �   #D_IFFT_1D%SIZE �   #SP �   #S �   #PARAM �                                              �     PRESENT                                            �     SIZE           D@                              �                    	              &                                                     D @                              �                   
               &                                                     F @                               �     P               #PSD_PARAM    #         @      X                                               #D_IFFT_2D%DBLE �   #D_IFFT_2D%PRESENT �   #D_IFFT_2D%SIZE �   #SP �   #S �   #PARAM �                                              �     DBLE                                            �     PRESENT                                            �     SIZE           D@                              �                    
              &                   &                                                     D                                �                   
               &                   &                                                     F @                               �     P               #PSD_PARAM    #         @      X                             }                   #D_XCOR_2D%SQRT �   #D_XCOR_2D%DCONJG �   #D_XCOR_2D%SUM �   #D_XCOR_2D%DBLE �   #D_XCOR_2D%PRESENT �   #D_XCOR_2D%SIZE �   #S1 �   #S2 �   #TAU �   #XCOR �   #PARAM �                                              �     SQRT                                            �     DCONJG                                            �     SUM                                            �     DBLE                                            �     PRESENT                                            �     SIZE           D@                              �                   
 N              &                   &                                                     D@                              �                   
 O              &                   &                                                     D@                              �                   
 Q              &                   &                   &                                                     D@                              �                   
 P              &                   &                                                      @                               �     P               #PSD_PARAM    #         @                                  �                    #PARAM �             D @                               �     P               #PSD_PARAM    #         @                                  �                   #GET_WINDOW%COS �   #GET_WINDOW%ATAN �   #GET_WINDOW%ABS �   #GET_WINDOW%DBLE �   #GET_WINDOW%PRESENT �   #GET_WINDOW%SIZE �   #TYPE �   #WINDOW �   #POWER �   #WIENER_K �                                              �     COS                                             �     ATAN                                            �     ABS                                            �     DBLE                                            �     PRESENT                                            �     SIZE                                            �     P               #PSD_PARAM              D@                              �                   
 W              &                                                     D                                �     
                  @                               �               �   &      fn#fn "   �   �   b   uapp(LIB_SPECTRAL    �  l       gen@FFT      �      D_FFT_1D !   �  @      D_FFT_1D%PRESENT    �  =      D_FFT_1D%SIZE      �   a   D_FFT_1D%S    �  �   a   D_FFT_1D%SP    5  W   a   D_FFT_1D%PARAM    �  �      D_FFT_1D_F     "  =      D_FFT_1D_F%DBLE #   _  @      D_FFT_1D_F%PRESENT    �  �   a   D_FFT_1D_F%S    +  �   a   D_FFT_1D_F%F    �  �   a   D_FFT_1D_F%SP !   C  W   a   D_FFT_1D_F%PARAM    �  ^       gen@IFFT    �  m       gen@PSD    e  �      D_PSD_1D    c	  =      D_PSD_1D%SQRT    �	  <      D_PSD_1D%ABS     �	  ?      D_PSD_1D%DCONJG    
  >      D_PSD_1D%DSQRT    Y
  <      D_PSD_1D%SUM    �
  =      D_PSD_1D%DBLE !   �
  @      D_PSD_1D%PRESENT      =      D_PSD_1D%SIZE    O  �   a   D_PSD_1D%S    �  �   a   D_PSD_1D%SP    g  W   a   D_PSD_1D%PARAM    �  �      D_PSD_1D_F     T  =      D_PSD_1D_F%DBLE #   �  @      D_PSD_1D_F%PRESENT    �  �   a   D_PSD_1D_F%S    ]  �   a   D_PSD_1D_F%F    �  �   a   D_PSD_1D_F%SP !   u  W   a   D_PSD_1D_F%PARAM    �  �      D_PSD_LDA    �  <      D_PSD_LDA%ABS    �  =      D_PSD_LDA%DBLE "     @      D_PSD_LDA%PRESENT    \  =      D_PSD_LDA%SIZE    �  =      D_PSD_LDA%REAL    �  �   a   D_PSD_LDA%AT    b  �   a   D_PSD_LDA%S    �  �   a   D_PSD_LDA%F    z  �   a   D_PSD_LDA%SP       W   a   D_PSD_LDA%PARAM    ]  ]       gen@COR    �  �      D_COR_1D    u  =      D_COR_1D%DBLE !   �  @      D_COR_1D%PRESENT    �  =      D_COR_1D%SIZE    /  =      D_COR_1D%REAL    l  �   a   D_COR_1D%S    �  �   a   D_COR_1D%TAU    �  �   a   D_COR_1D%COR      W   a   D_COR_1D%PARAM    g  �      D_COR_LDA    B  =      D_COR_LDA%SQRT      <      D_COR_LDA%SUM "   �  @      D_COR_LDA%PRESENT    �  =      D_COR_LDA%SIZE    8  =      D_COR_LDA%REAL    u  �   a   D_COR_LDA%AT      �   a   D_COR_LDA%S    �  �   a   D_COR_LDA%TAU      �   a   D_COR_LDA%XCOR     �  W   a   D_COR_LDA%PARAM    �  p       gen@XPSD    l       D_XPSD_1D    {  =      D_XPSD_1D%SQRT    �  <      D_XPSD_1D%ABS !   �  ?      D_XPSD_1D%DCONJG     3  >      D_XPSD_1D%DSQRT    q  <      D_XPSD_1D%SUM    �  =      D_XPSD_1D%DBLE "   �  @      D_XPSD_1D%PRESENT    *   =      D_XPSD_1D%SIZE    g   �   a   D_XPSD_1D%S1    �   �   a   D_XPSD_1D%S2    !  �   a   D_XPSD_1D%SP     "  W   a   D_XPSD_1D%PARAM    b"  �      D_XPSD_1D_F !   #  =      D_XPSD_1D_F%DBLE $   @#  @      D_XPSD_1D_F%PRESENT    �#  �   a   D_XPSD_1D_F%S1    $  �   a   D_XPSD_1D_F%S2    �$  �   a   D_XPSD_1D_F%F    $%  �   a   D_XPSD_1D_F%SP "   �%  W   a   D_XPSD_1D_F%PARAM    &  �      D_XPSD_LDA    �&  <      D_XPSD_LDA%ABS     2'  =      D_XPSD_LDA%DBLE #   o'  @      D_XPSD_LDA%PRESENT     �'  =      D_XPSD_LDA%SIZE     �'  =      D_XPSD_LDA%REAL    )(  �   a   D_XPSD_LDA%AT1    �(  �   a   D_XPSD_LDA%S1    A)  �   a   D_XPSD_LDA%AT2    �)  �   a   D_XPSD_LDA%S2    Y*  �   a   D_XPSD_LDA%F    �*  �   a   D_XPSD_LDA%SP !   q+  W   a   D_XPSD_LDA%PARAM    �+  n       gen@XCOR    6,  �      D_XCOR_1D    �,  =      D_XCOR_1D%DBLE "   <-  @      D_XCOR_1D%PRESENT    |-  =      D_XCOR_1D%SIZE    �-  =      D_XCOR_1D%REAL    �-  �   a   D_XCOR_1D%S1    �.  �   a   D_XCOR_1D%S2    /  �   a   D_XCOR_1D%TAU    �/  �   a   D_XCOR_1D%XCOR     &0  W   a   D_XCOR_1D%PARAM    }0  �      D_XCOR_LDA     p1  =      D_XCOR_LDA%SQRT    �1  <      D_XCOR_LDA%SUM #   �1  @      D_XCOR_LDA%PRESENT     )2  =      D_XCOR_LDA%SIZE     f2  =      D_XCOR_LDA%REAL    �2  �   a   D_XCOR_LDA%AT1    /3  �   a   D_XCOR_LDA%S1    �3  �   a   D_XCOR_LDA%AT2    G4  �   a   D_XCOR_LDA%S2    �4  �   a   D_XCOR_LDA%TAU     _5  �   a   D_XCOR_LDA%XCOR !   �5  W   a   D_XCOR_LDA%PARAM !   B6  h       gen@UNWRAP_PHASE    �6  �      UNWRAP_PHASE_F $   ^7  =      UNWRAP_PHASE_F%SIGN $   �7  =      UNWRAP_PHASE_F%ATAN #   �7  <      UNWRAP_PHASE_F%ABS $   8  =      UNWRAP_PHASE_F%SIZE #   Q8  �   a   UNWRAP_PHASE_F%PHI    �8  �      UNWRAP_PHASE_D $   �9  =      UNWRAP_PHASE_D%SIGN $   �9  =      UNWRAP_PHASE_D%ATAN #   :  <      UNWRAP_PHASE_D%ABS $   G:  =      UNWRAP_PHASE_D%SIZE #   �:  �   a   UNWRAP_PHASE_D%PHI    ;  q       FFTW_R2HC    �;  q       FFTW_HC2R    �;  q       FFTW_DHT    c<  q       FFTW_REDFT00    �<  q       FFTW_REDFT01    E=  q       FFTW_REDFT10    �=  q       FFTW_REDFT11    '>  q       FFTW_RODFT00    �>  q       FFTW_RODFT01    	?  q       FFTW_RODFT10    z?  r       FFTW_RODFT11    �?  p       FFTW_FORWARD    \@  q       FFTW_BACKWARD    �@  q       FFTW_MEASURE #   >A  q       FFTW_DESTROY_INPUT    �A  q       FFTW_UNALIGNED %    B  q       FFTW_CONSERVE_MEMORY     �B  q       FFTW_EXHAUSTIVE $   C  r       FFTW_PRESERVE_INPUT    tC  r       FFTW_PATIENT    �C  r       FFTW_ESTIMATE !   XD  w       FFTW_WISDOM_ONLY &   �D  s       FFTW_ESTIMATE_PATIENT #   BE  s       FFTW_BELIEVE_PCOST !   �E  s       FFTW_NO_DFT_R2HC $   (F  t       FFTW_NO_NONTHREADED "   �F  t       FFTW_NO_BUFFERING $   G  t       FFTW_NO_INDIRECT_OP )   �G  t       FFTW_ALLOW_LARGE_GENERIC $   �G  u       FFTW_NO_RANK_SPLITS %   mH  u       FFTW_NO_VRANK_SPLITS !   �H  u       FFTW_NO_VRECURSE    WI  v       FFTW_NO_SIMD    �I  v       FFTW_NO_SLOW ,   CJ  v       FFTW_NO_FIXED_RADIX_LARGE_N #   �J  w       FFTW_ALLOW_PRUNING    0K  �       PSD_PARAM    .L  �   a   PSD_PARAM%NFFT "   �L  �   a   PSD_PARAM%OVERLAP !   }M  �   a   PSD_PARAM%WINDOW    ;N  �   a   PSD_PARAM%FE    �N  �   a   PSD_PARAM%FMIN (   �O  �   a   PSD_PARAM%ALLOCATED_FFT )   0P  �   a   PSD_PARAM%ALLOCATED_IFFT #   �P  �   a   PSD_PARAM%NORM_FFT #   xQ  �   a   PSD_PARAM%RMS_NORM %   R  �   a   PSD_PARAM%CHECK_PVAL "   �R  �   a   PSD_PARAM%DETREND    dS  �   a   PSD_PARAM%PLAN $   	T  �   a   PSD_PARAM%PLAN_IFFT    �T  �       D_FFT_2D    LU  =      D_FFT_2D%DBLE !   �U  @      D_FFT_2D%PRESENT    �U  =      D_FFT_2D%SIZE    V  �   a   D_FFT_2D%S    �V  �   a   D_FFT_2D%SP    NW  W   a   D_FFT_2D%PARAM    �W  �       D_IFFT_1D "   2X  @      D_IFFT_1D%PRESENT    rX  =      D_IFFT_1D%SIZE    �X  �   a   D_IFFT_1D%SP    ;Y  �   a   D_IFFT_1D%S     �Y  W   a   D_IFFT_1D%PARAM    Z  �       D_IFFT_2D    �Z  =      D_IFFT_2D%DBLE "   �Z  @      D_IFFT_2D%PRESENT    <[  =      D_IFFT_2D%SIZE    y[  �   a   D_IFFT_2D%SP    \  �   a   D_IFFT_2D%S     �\  W   a   D_IFFT_2D%PARAM    ]  �       D_XCOR_2D    
^  =      D_XCOR_2D%SQRT !   G^  ?      D_XCOR_2D%DCONJG    �^  <      D_XCOR_2D%SUM    �^  =      D_XCOR_2D%DBLE "   �^  @      D_XCOR_2D%PRESENT    ?_  =      D_XCOR_2D%SIZE    |_  �   a   D_XCOR_2D%S1     `  �   a   D_XCOR_2D%S2    �`  �   a   D_XCOR_2D%TAU    �a  �   a   D_XCOR_2D%XCOR     $b  W   a   D_XCOR_2D%PARAM    {b  S       FREE_IFFT     �b  W   a   FREE_IFFT%PARAM    %c  �       GET_WINDOW    d  <      GET_WINDOW%COS     Wd  =      GET_WINDOW%ATAN    �d  <      GET_WINDOW%ABS     �d  =      GET_WINDOW%DBLE #   e  @      GET_WINDOW%PRESENT     Me  =      GET_WINDOW%SIZE     �e  W   a   GET_WINDOW%TYPE "   �e  �   a   GET_WINDOW%WINDOW !   mf  @   a   GET_WINDOW%POWER $   �f  @   a   GET_WINDOW%WIENER_K 