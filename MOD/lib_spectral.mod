	  /U  Â   k820309              13.0        ó%jT                                                                                                           
       lib_spectral.f90 LIB_SPECTRAL              D_FFT_1D D_FFT_1D_F D_PSD_1D D_PSD_1D_F D_COR_1D D_XPSD_1D D_XCOR_1D D_XPSD_1D_F D_COR_LDA D_PSD_LDA D_XCOR_LDA D_XPSD_LDA FREE_FFT RMLINTREND TRIANGLE SLOTTINGFUZZY I J K IN IF IC                                                        u #D_FFT_1D    #D_FFT_1D_F    #         @     @X                                               #D_FFT_1D%PRESENT    #D_FFT_1D%SIZE    #S    #SP    #PARAM                                                    PRESENT                                                 SIZE           D@                                                 
               &                                                     D@                                                                &                                                     F @                                    P               #PSD_PARAM    #         @     @X                                                #D_FFT_1D_F%DBLE 	   #D_FFT_1D_F%PRESENT 
   #S    #F    #SP    #PARAM                                               	     DBLE                                            
     PRESENT           D @                                                 
               &                                                     D                                                   
               &                                                     D @                                                                &                                                     F @                                    P               #PSD_PARAM                                                           u #D_IFFT_1D                                                           u #D_PSD_1D    #D_PSD_1D_F    #D_PSD_LDA "   #         @     @X                                               #D_PSD_1D%ABS    #D_PSD_1D%DCONJG    #D_PSD_1D%DSQRT    #D_PSD_1D%SUM    #D_PSD_1D%DBLE    #D_PSD_1D%PRESENT    #D_PSD_1D%SIZE    #S    #SP    #PARAM                                                    ABS                                                 DCONJG                                                 DSQRT                                                 SUM                                                 DBLE                                                 PRESENT                                                 SIZE           D@                                                 
               &                                                     D@                                                  	              &                                                      @                                    P               #PSD_PARAM    #         @     @X                                                #D_PSD_1D_F%DBLE    #D_PSD_1D_F%PRESENT    #S    #F    #SP     #PARAM !                                                   DBLE                                                 PRESENT           D @                                                 
               &                                                     D                                                   
               &                                                     D @                                                                 &                                                      @                               !     P               #PSD_PARAM    #         @     @X                             "                   #D_PSD_LDA%ABS #   #D_PSD_LDA%DBLE $   #D_PSD_LDA%PRESENT %   #D_PSD_LDA%SIZE &   #D_PSD_LDA%REAL '   #AT (   #S )   #F *   #SP +   #PARAM ,                                              #     ABS                                            $     DBLE                                            %     PRESENT                                            &     SIZE                                            '     REAL           D @                              (                   
               &                                                     D@                              )                   
               &                                                     D                                *                   
               &                                                     D @                              +                                  &                                                     F @                               ,     P               #PSD_PARAM                                                           u #D_COR_1D -   #D_COR_LDA 6   #         @     @X                             -                   #D_COR_1D%DBLE .   #D_COR_1D%PRESENT /   #D_COR_1D%SIZE 0   #D_COR_1D%REAL 1   #S 2   #TAU 3   #COR 4   #PARAM 5                                              .     DBLE                                            /     PRESENT                                            0     SIZE                                            1     REAL           D @                              2                   
 =              &                                                     D@                              3                   
 ?              &                                                     D@                              4                   
 >              &                                                      @                               5     P               #PSD_PARAM    #         @     @X                            6                   #D_COR_LDA%SQRT 7   #D_COR_LDA%SUM 8   #D_COR_LDA%PRESENT 9   #D_COR_LDA%SIZE :   #D_COR_LDA%REAL ;   #AT <   #S =   #TAU >   #XCOR ?   #PARAM @                                              7     SQRT                                            8     SUM                                            9     PRESENT                                            :     SIZE                                            ;     REAL           D @                              <                   
 +              &                                                     D@                              =                   
 ,              &                                                     D @                              >                   
 -              &                                                     D @                              ?                   
 .              &                                                      @                               @     P               #PSD_PARAM                                                           u #D_XPSD_1D A   #D_XPSD_1D_F M   #D_XPSD_LDA U   #         @     @X                            A                   #D_XPSD_1D%ABS B   #D_XPSD_1D%DCONJG C   #D_XPSD_1D%DSQRT D   #D_XPSD_1D%SUM E   #D_XPSD_1D%DBLE F   #D_XPSD_1D%PRESENT G   #D_XPSD_1D%SIZE H   #S1 I   #S2 J   #SP K   #PARAM L                                              B     ABS                                            C     DCONJG                                            D     DSQRT                                            E     SUM                                            F     DBLE                                            G     PRESENT                                            H     SIZE           D@                              I                   
               &                                                     D                                J                   
               &                                                     D@                              K                                  &                                                      @                               L     P               #PSD_PARAM    #         @     @X                             M                   #D_XPSD_1D_F%DBLE N   #D_XPSD_1D_F%PRESENT O   #S1 P   #S2 Q   #F R   #SP S   #PARAM T                                              N     DBLE                                            O     PRESENT           D @                              P                   
               &                                                     D @                              Q                   
               &                                                     D                                R                   
               &                                                     D @                              S                                  &                                                      @                               T     P               #PSD_PARAM    #         @     @X                             U                   #D_XPSD_LDA%ABS V   #D_XPSD_LDA%DBLE W   #D_XPSD_LDA%PRESENT X   #D_XPSD_LDA%SIZE Y   #D_XPSD_LDA%REAL Z   #AT1 [   #S1 \   #AT2 ]   #S2 ^   #F _   #SP `   #PARAM a                                              V     ABS                                            W     DBLE                                            X     PRESENT                                            Y     SIZE                                            Z     REAL           D @                              [                   
 "              &                                                     D@                              \                   
 #              &                                                     D @                              ]                   
 $              &                                                     D@                              ^                   
 %              &                                                     D                                _                   
 &              &                                                     D @                              `                    '              &                                                     F @                               a     P               #PSD_PARAM                                                           u #D_XCOR_1D b   #D_XCOR_LDA l   #         @     @X                             b                   #D_XCOR_1D%DBLE c   #D_XCOR_1D%PRESENT d   #D_XCOR_1D%SIZE e   #D_XCOR_1D%REAL f   #S1 g   #S2 h   #TAU i   #XCOR j   #PARAM k                                              c     DBLE                                            d     PRESENT                                            e     SIZE                                            f     REAL           D@                              g                   
 C              &                                                     D @                              h                   
 D              &                                                     D@                              i                   
 F              &                                                     D@                              j                   
 E              &                                                      @                               k     P               #PSD_PARAM    #         @     @X                            l                   #D_XCOR_LDA%SQRT m   #D_XCOR_LDA%SUM n   #D_XCOR_LDA%PRESENT o   #D_XCOR_LDA%SIZE p   #D_XCOR_LDA%REAL q   #AT1 r   #S1 s   #AT2 t   #S2 u   #TAU v   #XCOR w   #PARAM x                                              m     SQRT                                            n     SUM                                            o     PRESENT                                            p     SIZE                                            q     REAL           D @                              r                   
 /              &                                                     D@                              s                   
 0              &                                                     D @                              t                   
 1              &                                                     D@                              u                   
 2              &                                                     D @                              v                   
 3              &                                                     D @                              w                   
 4              &                                                      @                               x     P               #PSD_PARAM                                                 y                                                       0                                             z                                                      1                                             {                                                      2                                             |                                                      3                                             }                                                      4                                             ~                                                      5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                    	               9                                                                                    
               10                                                                                       ÿÿÿÿÿÿÿÿ                                                                                                           1                                                                                                    0                                                                                                   1                                                                                                   2                                                                                                   4                                                                                                   8                                                                                                   16                                                                                                    32                                                                                    @               64                                                                                                    2097152                                                                                                   128                                                                                                   256                                                                                                   512                                                                                                   1024                                                                                                   2048                                                                                                   4096                                                                                                    8192                                                                                     @              16384                                                                                                   32768                                                                                                   65536                                                                                                   131072                                                                                                   262144                                                                                                   524288                                                                                                   1048576                  @                               'P                    #NFFT    #OVERLAP    #WINDOW    #FE     #FMIN ¡   #ALLOCATED_FFT ¢   #ALLOCATED_IFFT £   #NORM_FFT ¤   #RMS_NORM ¥   #CHECK_PVAL ¦   #PLAN §   #PLAN_IFFT ¨                                                                                                                                            1024                                                                                                                                            512                                                                                                                                                      CH                                                                             
                                                
                       ð?        1.D0                                              ¡               
                                                
                       à?        0.5D0                                               ¢     (                                                                                                                                            £     ,                                                                                                                                            ¤     0                                                                                                                                            ¥     4       	                                                                                                                                     ¦     8       
                                                                                                                                    §     @                                                                                          0                                              ¨     H                                                                                          0    #         @      X                                               #D_IFFT_1D%PRESENT ©   #D_IFFT_1D%SIZE ª   #SP «   #S ¬   #PARAM ­                                              ©     PRESENT                                            ª     SIZE           D@                              «                                  &                                                     D @                              ¬                   
               &                                                     F @                               ­     P               #PSD_PARAM    #         @                                  ®                    #PARAM ¯             D @                               ¯     P               #PSD_PARAM    #         @                                  °                   #GET_WINDOW%COS ±   #GET_WINDOW%ATAN ²   #GET_WINDOW%ABS ³   #GET_WINDOW%DBLE ´   #GET_WINDOW%PRESENT µ   #GET_WINDOW%SIZE ¶   #TYPE ·   #WINDOW ¸   #POWER ¹   #WIENER_K º                                              ±     COS                                             ²     ATAN                                            ³     ABS                                            ´     DBLE                                            µ     PRESENT                                            ¶     SIZE                                            ·     P               #PSD_PARAM              D@                              ¸                   
 J              &                                                     D                                ¹     
                  @                               º                   &      fn#fn "   Æ   Å   b   uapp(LIB_SPECTRAL      ^       gen@FFT    é        D_FFT_1D !   t  @      D_FFT_1D%PRESENT    ´  =      D_FFT_1D%SIZE    ñ     a   D_FFT_1D%S    }     a   D_FFT_1D%SP    	  W   a   D_FFT_1D%PARAM    `        D_FFT_1D_F     ö  =      D_FFT_1D_F%DBLE #   3  @      D_FFT_1D_F%PRESENT    s     a   D_FFT_1D_F%S    ÿ     a   D_FFT_1D_F%F         a   D_FFT_1D_F%SP !     W   a   D_FFT_1D_F%PARAM    n  O       gen@IFFT    ½  m       gen@PSD    *  ë      D_PSD_1D    	  <      D_PSD_1D%ABS     Q	  ?      D_PSD_1D%DCONJG    	  >      D_PSD_1D%DSQRT    Î	  <      D_PSD_1D%SUM    

  =      D_PSD_1D%DBLE !   G
  @      D_PSD_1D%PRESENT    
  =      D_PSD_1D%SIZE    Ä
     a   D_PSD_1D%S    P     a   D_PSD_1D%SP    Ü  W   a   D_PSD_1D%PARAM    3        D_PSD_1D_F     É  =      D_PSD_1D_F%DBLE #     @      D_PSD_1D_F%PRESENT    F     a   D_PSD_1D_F%S    Ò     a   D_PSD_1D_F%F    ^     a   D_PSD_1D_F%SP !   ê  W   a   D_PSD_1D_F%PARAM    A  ×      D_PSD_LDA      <      D_PSD_LDA%ABS    T  =      D_PSD_LDA%DBLE "     @      D_PSD_LDA%PRESENT    Ñ  =      D_PSD_LDA%SIZE      =      D_PSD_LDA%REAL    K     a   D_PSD_LDA%AT    ×     a   D_PSD_LDA%S    c     a   D_PSD_LDA%F    ï     a   D_PSD_LDA%SP     {  W   a   D_PSD_LDA%PARAM    Ò  ]       gen@COR    /  »      D_COR_1D    ê  =      D_COR_1D%DBLE !   '  @      D_COR_1D%PRESENT    g  =      D_COR_1D%SIZE    ¤  =      D_COR_1D%REAL    á     a   D_COR_1D%S    m     a   D_COR_1D%TAU    ù     a   D_COR_1D%COR      W   a   D_COR_1D%PARAM    Ü  Û      D_COR_LDA    ·  =      D_COR_LDA%SQRT    ô  <      D_COR_LDA%SUM "   0  @      D_COR_LDA%PRESENT    p  =      D_COR_LDA%SIZE    ­  =      D_COR_LDA%REAL    ê     a   D_COR_LDA%AT    v     a   D_COR_LDA%S         a   D_COR_LDA%TAU         a   D_COR_LDA%XCOR       W   a   D_COR_LDA%PARAM    q  p       gen@XPSD    á  û      D_XPSD_1D    Ü  <      D_XPSD_1D%ABS !     ?      D_XPSD_1D%DCONJG     W  >      D_XPSD_1D%DSQRT      <      D_XPSD_1D%SUM    Ñ  =      D_XPSD_1D%DBLE "     @      D_XPSD_1D%PRESENT    N  =      D_XPSD_1D%SIZE         a   D_XPSD_1D%S1          a   D_XPSD_1D%S2    £      a   D_XPSD_1D%SP     /!  W   a   D_XPSD_1D%PARAM    !  ¡      D_XPSD_1D_F !   '"  =      D_XPSD_1D_F%DBLE $   d"  @      D_XPSD_1D_F%PRESENT    ¤"     a   D_XPSD_1D_F%S1    0#     a   D_XPSD_1D_F%S2    ¼#     a   D_XPSD_1D_F%F    H$     a   D_XPSD_1D_F%SP "   Ô$  W   a   D_XPSD_1D_F%PARAM    +%  ï      D_XPSD_LDA    &  <      D_XPSD_LDA%ABS     V&  =      D_XPSD_LDA%DBLE #   &  @      D_XPSD_LDA%PRESENT     Ó&  =      D_XPSD_LDA%SIZE     '  =      D_XPSD_LDA%REAL    M'     a   D_XPSD_LDA%AT1    Ù'     a   D_XPSD_LDA%S1    e(     a   D_XPSD_LDA%AT2    ñ(     a   D_XPSD_LDA%S2    })     a   D_XPSD_LDA%F    	*     a   D_XPSD_LDA%SP !   *  W   a   D_XPSD_LDA%PARAM    ì*  _       gen@XCOR    K+  É      D_XCOR_1D    ,  =      D_XCOR_1D%DBLE "   Q,  @      D_XCOR_1D%PRESENT    ,  =      D_XCOR_1D%SIZE    Î,  =      D_XCOR_1D%REAL    -     a   D_XCOR_1D%S1    -     a   D_XCOR_1D%S2    #.     a   D_XCOR_1D%TAU    ¯.     a   D_XCOR_1D%XCOR     ;/  W   a   D_XCOR_1D%PARAM    /  ó      D_XCOR_LDA     0  =      D_XCOR_LDA%SQRT    Â0  <      D_XCOR_LDA%SUM #   þ0  @      D_XCOR_LDA%PRESENT     >1  =      D_XCOR_LDA%SIZE     {1  =      D_XCOR_LDA%REAL    ¸1     a   D_XCOR_LDA%AT1    D2     a   D_XCOR_LDA%S1    Ð2     a   D_XCOR_LDA%AT2    \3     a   D_XCOR_LDA%S2    è3     a   D_XCOR_LDA%TAU     t4     a   D_XCOR_LDA%XCOR !    5  W   a   D_XCOR_LDA%PARAM    W5  q       FFTW_R2HC    È5  q       FFTW_HC2R    96  q       FFTW_DHT    ª6  q       FFTW_REDFT00    7  q       FFTW_REDFT01    7  q       FFTW_REDFT10    ý7  q       FFTW_REDFT11    n8  q       FFTW_RODFT00    ß8  q       FFTW_RODFT01    P9  q       FFTW_RODFT10    Á9  r       FFTW_RODFT11    3:  p       FFTW_FORWARD    £:  q       FFTW_BACKWARD    ;  q       FFTW_MEASURE #   ;  q       FFTW_DESTROY_INPUT    ö;  q       FFTW_UNALIGNED %   g<  q       FFTW_CONSERVE_MEMORY     Ø<  q       FFTW_EXHAUSTIVE $   I=  r       FFTW_PRESERVE_INPUT    »=  r       FFTW_PATIENT    ->  r       FFTW_ESTIMATE !   >  w       FFTW_WISDOM_ONLY &   ?  s       FFTW_ESTIMATE_PATIENT #   ?  s       FFTW_BELIEVE_PCOST !   ü?  s       FFTW_NO_DFT_R2HC $   o@  t       FFTW_NO_NONTHREADED "   ã@  t       FFTW_NO_BUFFERING $   WA  t       FFTW_NO_INDIRECT_OP )   ËA  t       FFTW_ALLOW_LARGE_GENERIC $   ?B  u       FFTW_NO_RANK_SPLITS %   ´B  u       FFTW_NO_VRANK_SPLITS !   )C  u       FFTW_NO_VRECURSE    C  v       FFTW_NO_SIMD    D  v       FFTW_NO_SLOW ,   D  v       FFTW_NO_FIXED_RADIX_LARGE_N #    E  w       FFTW_ALLOW_PRUNING    wE  ñ       PSD_PARAM    hF  ¨   a   PSD_PARAM%NFFT "   G  §   a   PSD_PARAM%OVERLAP !   ·G  ¾   a   PSD_PARAM%WINDOW    uH  ¨   a   PSD_PARAM%FE    I  ©   a   PSD_PARAM%FMIN (   ÆI  ¤   a   PSD_PARAM%ALLOCATED_FFT )   jJ  ¤   a   PSD_PARAM%ALLOCATED_IFFT #   K  ¤   a   PSD_PARAM%NORM_FFT #   ²K  ¤   a   PSD_PARAM%RMS_NORM %   VL  ¤   a   PSD_PARAM%CHECK_PVAL    úL  ¥   a   PSD_PARAM%PLAN $   M  ¥   a   PSD_PARAM%PLAN_IFFT    DN         D_IFFT_1D "   ÑN  @      D_IFFT_1D%PRESENT    O  =      D_IFFT_1D%SIZE    NO     a   D_IFFT_1D%SP    ÚO     a   D_IFFT_1D%S     fP  W   a   D_IFFT_1D%PARAM    ½P  S       FREE_IFFT     Q  W   a   FREE_IFFT%PARAM    gQ  ö       GET_WINDOW    ]R  <      GET_WINDOW%COS     R  =      GET_WINDOW%ATAN    ÖR  <      GET_WINDOW%ABS     S  =      GET_WINDOW%DBLE #   OS  @      GET_WINDOW%PRESENT     S  =      GET_WINDOW%SIZE     ÌS  W   a   GET_WINDOW%TYPE "   #T     a   GET_WINDOW%WINDOW !   ¯T  @   a   GET_WINDOW%POWER $   ïT  @   a   GET_WINDOW%WIENER_K 