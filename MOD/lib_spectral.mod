	  �`  �   k820309              13.0        BT                                                                                                           
       lib_spectral.f90 LIB_SPECTRAL              D_FFT_1D D_FFT_1D_F D_PSD_1D D_PSD_1D_F D_COR_1D D_XPSD_1D D_XCOR_1D D_XPSD_1D_F D_COR_LDA D_PSD_LDA D_XCOR_LDA D_XPSD_LDA FREE_FFT RMLINTREND TRIANGLE SLOTTINGFUZZY I J K IN IF IC                                                        u #D_FFT_1D    #D_FFT_1D_F    #D_FFT_2D    #         @     @X                                               #D_FFT_1D%PRESENT    #D_FFT_1D%SIZE    #S    #SP    #PARAM                                                    PRESENT                                                 SIZE           D@                                                 
               &                                                     D@                                                                &                                                     F @                                    P               #PSD_PARAM    #         @     @X                                                #D_FFT_1D_F%DBLE 	   #D_FFT_1D_F%PRESENT 
   #S    #F    #SP    #PARAM                                               	     DBLE                                            
     PRESENT           D @                                                 
               &                                                     D                                                   
               &                                                     D @                                                                &                                                     F @                                    P               #PSD_PARAM                                                           u #D_IFFT_1D    #D_IFFT_2D                                                           u #D_PSD_1D    #D_PSD_1D_F    #D_PSD_LDA $   #         @     @X                                               #D_PSD_1D%ABS    #D_PSD_1D%DCONJG    #D_PSD_1D%DSQRT    #D_PSD_1D%SUM    #D_PSD_1D%DBLE    #D_PSD_1D%PRESENT    #D_PSD_1D%SIZE    #S    #SP    #PARAM                                                    ABS                                                 DCONJG                                                 DSQRT                                                 SUM                                                 DBLE                                                 PRESENT                                                 SIZE           D@                                                 
               &                                                     D@                                                                &                                                      @                                    P               #PSD_PARAM    #         @     @X                                                #D_PSD_1D_F%DBLE    #D_PSD_1D_F%PRESENT    #S     #F !   #SP "   #PARAM #                                                   DBLE                                                 PRESENT           D @                                                  
               &                                                     D                                !                   
               &                                                     D @                              "                                  &                                                      @                               #     P               #PSD_PARAM    #         @     @X                             $                   #D_PSD_LDA%ABS %   #D_PSD_LDA%DBLE &   #D_PSD_LDA%PRESENT '   #D_PSD_LDA%SIZE (   #D_PSD_LDA%REAL )   #AT *   #S +   #F ,   #SP -   #PARAM .                                              %     ABS                                            &     DBLE                                            '     PRESENT                                            (     SIZE                                            )     REAL           D @                              *                   
               &                                                     D@                              +                   
                &                                                     D                                ,                   
 !              &                                                     D @                              -                    "              &                                                     F @                               .     P               #PSD_PARAM                                                           u #D_COR_1D /   #D_COR_LDA 8   #         @     @X                             /                   #D_COR_1D%DBLE 0   #D_COR_1D%PRESENT 1   #D_COR_1D%SIZE 2   #D_COR_1D%REAL 3   #S 4   #TAU 5   #COR 6   #PARAM 7                                              0     DBLE                                            1     PRESENT                                            2     SIZE                                            3     REAL           D @                              4                   
 A              &                                                     D@                              5                   
 C              &                                                     D@                              6                   
 B              &                                                      @                               7     P               #PSD_PARAM    #         @     @X                            8                   #D_COR_LDA%SQRT 9   #D_COR_LDA%SUM :   #D_COR_LDA%PRESENT ;   #D_COR_LDA%SIZE <   #D_COR_LDA%REAL =   #AT >   #S ?   #TAU @   #XCOR A   #PARAM B                                              9     SQRT                                            :     SUM                                            ;     PRESENT                                            <     SIZE                                            =     REAL           D @                              >                   
 /              &                                                     D@                              ?                   
 0              &                                                     D @                              @                   
 1              &                                                     D @                              A                   
 2              &                                                      @                               B     P               #PSD_PARAM                                                           u #D_XPSD_1D C   #D_XPSD_1D_F O   #D_XPSD_LDA W   #         @     @X                            C                   #D_XPSD_1D%ABS D   #D_XPSD_1D%DCONJG E   #D_XPSD_1D%DSQRT F   #D_XPSD_1D%SUM G   #D_XPSD_1D%DBLE H   #D_XPSD_1D%PRESENT I   #D_XPSD_1D%SIZE J   #S1 K   #S2 L   #SP M   #PARAM N                                              D     ABS                                            E     DCONJG                                            F     DSQRT                                            G     SUM                                            H     DBLE                                            I     PRESENT                                            J     SIZE           D@                              K                   
               &                                                     D                                L                   
               &                                                     D@                              M                                  &                                                      @                               N     P               #PSD_PARAM    #         @     @X                             O                   #D_XPSD_1D_F%DBLE P   #D_XPSD_1D_F%PRESENT Q   #S1 R   #S2 S   #F T   #SP U   #PARAM V                                              P     DBLE                                            Q     PRESENT           D @                              R                   
               &                                                     D @                              S                   
               &                                                     D                                T                   
               &                                                     D @                              U                                  &                                                      @                               V     P               #PSD_PARAM    #         @     @X                             W                   #D_XPSD_LDA%ABS X   #D_XPSD_LDA%DBLE Y   #D_XPSD_LDA%PRESENT Z   #D_XPSD_LDA%SIZE [   #D_XPSD_LDA%REAL \   #AT1 ]   #S1 ^   #AT2 _   #S2 `   #F a   #SP b   #PARAM c                                              X     ABS                                            Y     DBLE                                            Z     PRESENT                                            [     SIZE                                            \     REAL           D @                              ]                   
 &              &                                                     D@                              ^                   
 '              &                                                     D @                              _                   
 (              &                                                     D@                              `                   
 )              &                                                     D                                a                   
 *              &                                                     D @                              b                    +              &                                                     F @                               c     P               #PSD_PARAM                                                           u #D_XCOR_1D d   #D_XCOR_LDA n   #D_XCOR_2D {   #         @     @X                             d                   #D_XCOR_1D%DBLE e   #D_XCOR_1D%PRESENT f   #D_XCOR_1D%SIZE g   #D_XCOR_1D%REAL h   #S1 i   #S2 j   #TAU k   #XCOR l   #PARAM m                                              e     DBLE                                            f     PRESENT                                            g     SIZE                                            h     REAL           D@                              i                   
 G              &                                                     D @                              j                   
 H              &                                                     D@                              k                   
 J              &                                                     D@                              l                   
 I              &                                                      @                               m     P               #PSD_PARAM    #         @     @X                            n                   #D_XCOR_LDA%SQRT o   #D_XCOR_LDA%SUM p   #D_XCOR_LDA%PRESENT q   #D_XCOR_LDA%SIZE r   #D_XCOR_LDA%REAL s   #AT1 t   #S1 u   #AT2 v   #S2 w   #TAU x   #XCOR y   #PARAM z                                              o     SQRT                                            p     SUM                                            q     PRESENT                                            r     SIZE                                            s     REAL           D @                              t                   
 3              &                                                     D@                              u                   
 4              &                                                     D @                              v                   
 5              &                                                     D@                              w                   
 6              &                                                     D @                              x                   
 7              &                                                     D @                              y                   
 8              &                                                      @                               z     P               #PSD_PARAM                                                 |                                                       0                                             }                                                      1                                             ~                                                      2                                                                                                   3                                             �                                                      4                                             �                                                      5                                             �                                                      6                                             �                                                      7                                             �                                                      8                                             �                                       	               9                                             �                                       
               10                                             �                                          ��������                                                     �                                                      1                                             �                                                       0                                             �                                                      1                                             �                                                      2                                             �                                                      4                                             �                                                      8                                             �                                                      16                                             �                                                       32                                             �                                       @               64                                             �                                                       2097152                                             �                                       �               128                                             �                                                      256                                             �                                                      512                                             �                                                      1024                                             �                                                      2048                                             �                                                      4096                                             �                                                       8192                                             �                                        @              16384                                             �                                        �              32768                                             �                                                      65536                                             �                                                      131072                                             �                                                      262144                                             �                                                      524288                                             �                                                      1048576                  @                               'P                    #NFFT �   #OVERLAP �   #WINDOW �   #FE �   #FMIN �   #ALLOCATED_FFT �   #ALLOCATED_IFFT �   #NORM_FFT �   #RMS_NORM �   #CHECK_PVAL �   #PLAN �   #PLAN_IFFT �               �                              �                                                                                               1024                �                              �                                                                                              512                �                              �                                                                                                        CH                                �                              �              
                                                
                       �?        1.D0                �                              �               
                                                
                       �?        0.5D0                �                               �     (                                                                                                             �                               �     ,                                                                                                             �                               �     0                                                                                                             �                               �     4       	                                                                                                      �                               �     8       
                                                                                                      �                              �     @                                                                                          0                �                              �     H                                                                                          0    #         @      X                                               #D_FFT_2D%DBLE �   #D_FFT_2D%PRESENT �   #D_FFT_2D%SIZE �   #S �   #SP �   #PARAM �                                              �     DBLE                                            �     PRESENT                                            �     SIZE           @                              �                   
               &                   &                                                     D@                              �                                  &                   &                                                     F @                               �     P               #PSD_PARAM    #         @      X                                               #D_IFFT_1D%PRESENT �   #D_IFFT_1D%SIZE �   #SP �   #S �   #PARAM �                                              �     PRESENT                                            �     SIZE           D@                              �                    	              &                                                     D @                              �                   
               &                                                     F @                               �     P               #PSD_PARAM    #         @      X                                               #D_IFFT_2D%DBLE �   #D_IFFT_2D%PRESENT �   #D_IFFT_2D%SIZE �   #SP �   #S �   #PARAM �                                              �     DBLE                                            �     PRESENT                                            �     SIZE           D@                              �                    
              &                   &                                                     D                                �                   
               &                   &                                                     F @                               �     P               #PSD_PARAM    #         @      X                             {                   #D_XCOR_2D%SQRT �   #D_XCOR_2D%DCONJG �   #D_XCOR_2D%SUM �   #D_XCOR_2D%DBLE �   #D_XCOR_2D%PRESENT �   #D_XCOR_2D%SIZE �   #S1 �   #S2 �   #TAU �   #XCOR �   #PARAM �                                              �     SQRT                                            �     DCONJG                                            �     SUM                                            �     DBLE                                            �     PRESENT                                            �     SIZE           D@                              �                   
 N              &                   &                                                     D@                              �                   
 O              &                   &                                                     D@                              �                   
 Q              &                   &                   &                                                     D@                              �                   
 P              &                   &                                                      @                               �     P               #PSD_PARAM    #         @                                  �                    #PARAM �             D @                               �     P               #PSD_PARAM    #         @                                  �                   #GET_WINDOW%COS �   #GET_WINDOW%ATAN �   #GET_WINDOW%ABS �   #GET_WINDOW%DBLE �   #GET_WINDOW%PRESENT �   #GET_WINDOW%SIZE �   #TYPE �   #WINDOW �   #POWER �   #WIENER_K �                                              �     COS                                             �     ATAN                                            �     ABS                                            �     DBLE                                            �     PRESENT                                            �     SIZE                                            �     P               #PSD_PARAM              D@                              �                   
 W              &                                                     D                                �     
                  @                               �               �   &      fn#fn "   �   �   b   uapp(LIB_SPECTRAL    �  l       gen@FFT    �  �      D_FFT_1D !   �  @      D_FFT_1D%PRESENT    �  =      D_FFT_1D%SIZE    �  �   a   D_FFT_1D%S    �  �   a   D_FFT_1D%SP      W   a   D_FFT_1D%PARAM    n  �      D_FFT_1D_F       =      D_FFT_1D_F%DBLE #   A  @      D_FFT_1D_F%PRESENT    �  �   a   D_FFT_1D_F%S      �   a   D_FFT_1D_F%F    �  �   a   D_FFT_1D_F%SP !   %  W   a   D_FFT_1D_F%PARAM    |  ^       gen@IFFT    �  m       gen@PSD    G  �      D_PSD_1D    2	  <      D_PSD_1D%ABS     n	  ?      D_PSD_1D%DCONJG    �	  >      D_PSD_1D%DSQRT    �	  <      D_PSD_1D%SUM    '
  =      D_PSD_1D%DBLE !   d
  @      D_PSD_1D%PRESENT    �
  =      D_PSD_1D%SIZE    �
  �   a   D_PSD_1D%S    m  �   a   D_PSD_1D%SP    �  W   a   D_PSD_1D%PARAM    P  �      D_PSD_1D_F     �  =      D_PSD_1D_F%DBLE #   #  @      D_PSD_1D_F%PRESENT    c  �   a   D_PSD_1D_F%S    �  �   a   D_PSD_1D_F%F    {  �   a   D_PSD_1D_F%SP !     W   a   D_PSD_1D_F%PARAM    ^  �      D_PSD_LDA    5  <      D_PSD_LDA%ABS    q  =      D_PSD_LDA%DBLE "   �  @      D_PSD_LDA%PRESENT    �  =      D_PSD_LDA%SIZE    +  =      D_PSD_LDA%REAL    h  �   a   D_PSD_LDA%AT    �  �   a   D_PSD_LDA%S    �  �   a   D_PSD_LDA%F      �   a   D_PSD_LDA%SP     �  W   a   D_PSD_LDA%PARAM    �  ]       gen@COR    L  �      D_COR_1D      =      D_COR_1D%DBLE !   D  @      D_COR_1D%PRESENT    �  =      D_COR_1D%SIZE    �  =      D_COR_1D%REAL    �  �   a   D_COR_1D%S    �  �   a   D_COR_1D%TAU      �   a   D_COR_1D%COR    �  W   a   D_COR_1D%PARAM    �  �      D_COR_LDA    �  =      D_COR_LDA%SQRT      <      D_COR_LDA%SUM "   M  @      D_COR_LDA%PRESENT    �  =      D_COR_LDA%SIZE    �  =      D_COR_LDA%REAL      �   a   D_COR_LDA%AT    �  �   a   D_COR_LDA%S      �   a   D_COR_LDA%TAU    �  �   a   D_COR_LDA%XCOR     7  W   a   D_COR_LDA%PARAM    �  p       gen@XPSD    �  �      D_XPSD_1D    �  <      D_XPSD_1D%ABS !   5  ?      D_XPSD_1D%DCONJG     t  >      D_XPSD_1D%DSQRT    �  <      D_XPSD_1D%SUM    �  =      D_XPSD_1D%DBLE "   +  @      D_XPSD_1D%PRESENT    k  =      D_XPSD_1D%SIZE    �  �   a   D_XPSD_1D%S1    4   �   a   D_XPSD_1D%S2    �   �   a   D_XPSD_1D%SP     L!  W   a   D_XPSD_1D%PARAM    �!  �      D_XPSD_1D_F !   D"  =      D_XPSD_1D_F%DBLE $   �"  @      D_XPSD_1D_F%PRESENT    �"  �   a   D_XPSD_1D_F%S1    M#  �   a   D_XPSD_1D_F%S2    �#  �   a   D_XPSD_1D_F%F    e$  �   a   D_XPSD_1D_F%SP "   �$  W   a   D_XPSD_1D_F%PARAM    H%  �      D_XPSD_LDA    7&  <      D_XPSD_LDA%ABS     s&  =      D_XPSD_LDA%DBLE #   �&  @      D_XPSD_LDA%PRESENT     �&  =      D_XPSD_LDA%SIZE     -'  =      D_XPSD_LDA%REAL    j'  �   a   D_XPSD_LDA%AT1    �'  �   a   D_XPSD_LDA%S1    �(  �   a   D_XPSD_LDA%AT2    )  �   a   D_XPSD_LDA%S2    �)  �   a   D_XPSD_LDA%F    &*  �   a   D_XPSD_LDA%SP !   �*  W   a   D_XPSD_LDA%PARAM    	+  n       gen@XCOR    w+  �      D_XCOR_1D    @,  =      D_XCOR_1D%DBLE "   },  @      D_XCOR_1D%PRESENT    �,  =      D_XCOR_1D%SIZE    �,  =      D_XCOR_1D%REAL    7-  �   a   D_XCOR_1D%S1    �-  �   a   D_XCOR_1D%S2    O.  �   a   D_XCOR_1D%TAU    �.  �   a   D_XCOR_1D%XCOR     g/  W   a   D_XCOR_1D%PARAM    �/  �      D_XCOR_LDA     �0  =      D_XCOR_LDA%SQRT    �0  <      D_XCOR_LDA%SUM #   *1  @      D_XCOR_LDA%PRESENT     j1  =      D_XCOR_LDA%SIZE     �1  =      D_XCOR_LDA%REAL    �1  �   a   D_XCOR_LDA%AT1    p2  �   a   D_XCOR_LDA%S1    �2  �   a   D_XCOR_LDA%AT2    �3  �   a   D_XCOR_LDA%S2    4  �   a   D_XCOR_LDA%TAU     �4  �   a   D_XCOR_LDA%XCOR !   ,5  W   a   D_XCOR_LDA%PARAM    �5  q       FFTW_R2HC    �5  q       FFTW_HC2R    e6  q       FFTW_DHT    �6  q       FFTW_REDFT00    G7  q       FFTW_REDFT01    �7  q       FFTW_REDFT10    )8  q       FFTW_REDFT11    �8  q       FFTW_RODFT00    9  q       FFTW_RODFT01    |9  q       FFTW_RODFT10    �9  r       FFTW_RODFT11    _:  p       FFTW_FORWARD    �:  q       FFTW_BACKWARD    @;  q       FFTW_MEASURE #   �;  q       FFTW_DESTROY_INPUT    "<  q       FFTW_UNALIGNED %   �<  q       FFTW_CONSERVE_MEMORY     =  q       FFTW_EXHAUSTIVE $   u=  r       FFTW_PRESERVE_INPUT    �=  r       FFTW_PATIENT    Y>  r       FFTW_ESTIMATE !   �>  w       FFTW_WISDOM_ONLY &   B?  s       FFTW_ESTIMATE_PATIENT #   �?  s       FFTW_BELIEVE_PCOST !   (@  s       FFTW_NO_DFT_R2HC $   �@  t       FFTW_NO_NONTHREADED "   A  t       FFTW_NO_BUFFERING $   �A  t       FFTW_NO_INDIRECT_OP )   �A  t       FFTW_ALLOW_LARGE_GENERIC $   kB  u       FFTW_NO_RANK_SPLITS %   �B  u       FFTW_NO_VRANK_SPLITS !   UC  u       FFTW_NO_VRECURSE    �C  v       FFTW_NO_SIMD    @D  v       FFTW_NO_SLOW ,   �D  v       FFTW_NO_FIXED_RADIX_LARGE_N #   ,E  w       FFTW_ALLOW_PRUNING    �E  �       PSD_PARAM    �F  �   a   PSD_PARAM%NFFT "   <G  �   a   PSD_PARAM%OVERLAP !   �G  �   a   PSD_PARAM%WINDOW    �H  �   a   PSD_PARAM%FE    II  �   a   PSD_PARAM%FMIN (   �I  �   a   PSD_PARAM%ALLOCATED_FFT )   �J  �   a   PSD_PARAM%ALLOCATED_IFFT #   :K  �   a   PSD_PARAM%NORM_FFT #   �K  �   a   PSD_PARAM%RMS_NORM %   �L  �   a   PSD_PARAM%CHECK_PVAL    &M  �   a   PSD_PARAM%PLAN $   �M  �   a   PSD_PARAM%PLAN_IFFT    pN  �       D_FFT_2D    O  =      D_FFT_2D%DBLE !   KO  @      D_FFT_2D%PRESENT    �O  =      D_FFT_2D%SIZE    �O  �   a   D_FFT_2D%S    lP  �   a   D_FFT_2D%SP    Q  W   a   D_FFT_2D%PARAM    gQ  �       D_IFFT_1D "   �Q  @      D_IFFT_1D%PRESENT    4R  =      D_IFFT_1D%SIZE    qR  �   a   D_IFFT_1D%SP    �R  �   a   D_IFFT_1D%S     �S  W   a   D_IFFT_1D%PARAM    �S  �       D_IFFT_2D    �T  =      D_IFFT_2D%DBLE "   �T  @      D_IFFT_2D%PRESENT    �T  =      D_IFFT_2D%SIZE    ;U  �   a   D_IFFT_2D%SP    �U  �   a   D_IFFT_2D%S     �V  W   a   D_IFFT_2D%PARAM    �V  �       D_XCOR_2D    �W  =      D_XCOR_2D%SQRT !   	X  ?      D_XCOR_2D%DCONJG    HX  <      D_XCOR_2D%SUM    �X  =      D_XCOR_2D%DBLE "   �X  @      D_XCOR_2D%PRESENT    Y  =      D_XCOR_2D%SIZE    >Y  �   a   D_XCOR_2D%S1    �Y  �   a   D_XCOR_2D%S2    �Z  �   a   D_XCOR_2D%TAU    B[  �   a   D_XCOR_2D%XCOR     �[  W   a   D_XCOR_2D%PARAM    =\  S       FREE_IFFT     �\  W   a   FREE_IFFT%PARAM    �\  �       GET_WINDOW    �]  <      GET_WINDOW%COS     ^  =      GET_WINDOW%ATAN    V^  <      GET_WINDOW%ABS     �^  =      GET_WINDOW%DBLE #   �^  @      GET_WINDOW%PRESENT     _  =      GET_WINDOW%SIZE     L_  W   a   GET_WINDOW%TYPE "   �_  �   a   GET_WINDOW%WINDOW !   /`  @   a   GET_WINDOW%POWER $   o`  @   a   GET_WINDOW%WIENER_K 