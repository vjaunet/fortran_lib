	  �  M   k820309              13.0        :k�T                                                                                                           
       lib_stat.f90 LIB_STAT              D_AVERAGE_1D_1C D_RMS_1D_1C D_RMS_1D_1C_MOY F_AVERAGE_1D_1C F_RMS_1D_1C F_RMS_1D_1C_MOY I J K IC IN NN NC II NX NY                                                        u #D_AVERAGE_1D_1C    #F_AVERAGE_1D_1C    #         @     @X                                               #D_AVERAGE_1D_1C%SUM    #D_AVERAGE_1D_1C%PRESENT    #D_AVERAGE_1D_1C%SIZE    #VAR    #MOY    #WEIGHT                                                    SUM                                                 PRESENT                                                 SIZE           @                                                 
               &                                                     D                                     
                  @                                                 
               &                                           #         @     @X                                               #F_AVERAGE_1D_1C%SUM 	   #F_AVERAGE_1D_1C%PRESENT 
   #F_AVERAGE_1D_1C%SIZE    #VAR    #MOY    #WEIGHT                                               	     SUM                                            
     PRESENT                                                 SIZE           @                                                 	               &                                                     D                                     	                  @                                                 	               &                                                                                                  u #D_RMS_1D_1C    #D_RMS_1D_1C_MOY    #F_RMS_1D_1C     #F_RMS_1D_1C_MOY (   #         @     @X                                                #D_RMS_1D_1C%SQRT    #D_RMS_1D_1C%SUM    #D_RMS_1D_1C%PRESENT    #D_RMS_1D_1C%SIZE    #VAR    #RMS    #WEIGHT                                                    SQRT                                                 SUM                                                 PRESENT                                                 SIZE           D@                                                 
               &                                                     D @                                   
                  @                                                 
               &                                           #         @     @X                                                #D_RMS_1D_1C_MOY%SQRT    #D_RMS_1D_1C_MOY%SUM    #D_RMS_1D_1C_MOY%PRESENT    #D_RMS_1D_1C_MOY%SIZE    #VAR    #MOY    #RMS    #WEIGHT                                                    SQRT                                                 SUM                                                 PRESENT                                                 SIZE           @                                                 
               &                                                                                          
                 D @                                   
                  @                                                 
               &                                           #         @     @X                                                 #F_RMS_1D_1C%SQRT !   #F_RMS_1D_1C%SUM "   #F_RMS_1D_1C%PRESENT #   #F_RMS_1D_1C%SIZE $   #VAR %   #RMS &   #WEIGHT '                                              !     SQRT                                            "     SUM                                            #     PRESENT                                            $     SIZE           D@                              %                   	 
              &                                                     D @                              &     	                  @                              '                   	               &                                           #         @     @X                            (                   #F_RMS_1D_1C_MOY%SQRT )   #F_RMS_1D_1C_MOY%SUM *   #F_RMS_1D_1C_MOY%PRESENT +   #F_RMS_1D_1C_MOY%SIZE ,   #VAR -   #MOY .   #RMS /   #WEIGHT 0                                              )     SQRT                                            *     SUM                                            +     PRESENT                                            ,     SIZE           @                              -                   	               &                                                                                     .     	                 D @                              /     	                  @                              0                   	               &                                                                                                  u #F_SKEWNESS_1D_1C 1                                                          u #F_FLATNESS_1D_1C 2                                                          u #F_XMOM_1D_1C 3   #         @      X                             1                   #F_SKEWNESS_1D_1C%SUM 4   #F_SKEWNESS_1D_1C%PRESENT 5   #F_SKEWNESS_1D_1C%SIZE 6   #VAR 7   #SKEW 8   #WEIGHT 9                                              4     SUM                                            5     PRESENT                                            6     SIZE           D@                              7                   	               &                                                     D                                8     	                  @                              9                   	               &                                           #         @      X                             2                   #F_FLATNESS_1D_1C%SUM :   #F_FLATNESS_1D_1C%PRESENT ;   #F_FLATNESS_1D_1C%SIZE <   #VAR =   #FLAT >   #WEIGHT ?                                              :     SUM                                            ;     PRESENT                                            <     SIZE           D@                              =                   	               &                                                     D                                >     	                  @                              ?                   	               &                                           #         @      X                             3                   #F_XMOM_1D_1C%SUM @   #F_XMOM_1D_1C%PRESENT A   #F_XMOM_1D_1C%SIZE B   #VAR1 C   #VAR2 D   #XMOM E   #WEIGHT F                                              @     SUM                                            A     PRESENT                                            B     SIZE           D@                              C                   	               &                                                     D @                              D                   	               &                                                     D                                E     	                  @                              F                   	               &                                              �         fn#fn    �   �   b   uapp(LIB_STAT    A  j       gen@AVERAGE     �  �      D_AVERAGE_1D_1C $   a  <      D_AVERAGE_1D_1C%SUM (   �  @      D_AVERAGE_1D_1C%PRESENT %   �  =      D_AVERAGE_1D_1C%SIZE $     �   a   D_AVERAGE_1D_1C%VAR $   �  @   a   D_AVERAGE_1D_1C%MOY '   �  �   a   D_AVERAGE_1D_1C%WEIGHT     r  �      F_AVERAGE_1D_1C $   (  <      F_AVERAGE_1D_1C%SUM (   d  @      F_AVERAGE_1D_1C%PRESENT %   �  =      F_AVERAGE_1D_1C%SIZE $   �  �   a   F_AVERAGE_1D_1C%VAR $   m  @   a   F_AVERAGE_1D_1C%MOY '   �  �   a   F_AVERAGE_1D_1C%WEIGHT    9  �       gen@RMS    �  �      D_RMS_1D_1C !   �  =      D_RMS_1D_1C%SQRT     �  <      D_RMS_1D_1C%SUM $   �  @      D_RMS_1D_1C%PRESENT !   >	  =      D_RMS_1D_1C%SIZE     {	  �   a   D_RMS_1D_1C%VAR     
  @   a   D_RMS_1D_1C%RMS #   G
  �   a   D_RMS_1D_1C%WEIGHT     �
  �      D_RMS_1D_1C_MOY %   �  =      D_RMS_1D_1C_MOY%SQRT $   �  <      D_RMS_1D_1C_MOY%SUM (   %  @      D_RMS_1D_1C_MOY%PRESENT %   e  =      D_RMS_1D_1C_MOY%SIZE $   �  �   a   D_RMS_1D_1C_MOY%VAR $   .  @   a   D_RMS_1D_1C_MOY%MOY $   n  @   a   D_RMS_1D_1C_MOY%RMS '   �  �   a   D_RMS_1D_1C_MOY%WEIGHT    :  �      F_RMS_1D_1C !   �  =      F_RMS_1D_1C%SQRT     7  <      F_RMS_1D_1C%SUM $   s  @      F_RMS_1D_1C%PRESENT !   �  =      F_RMS_1D_1C%SIZE     �  �   a   F_RMS_1D_1C%VAR     |  @   a   F_RMS_1D_1C%RMS #   �  �   a   F_RMS_1D_1C%WEIGHT     H  �      F_RMS_1D_1C_MOY %   !  =      F_RMS_1D_1C_MOY%SQRT $   ^  <      F_RMS_1D_1C_MOY%SUM (   �  @      F_RMS_1D_1C_MOY%PRESENT %   �  =      F_RMS_1D_1C_MOY%SIZE $     �   a   F_RMS_1D_1C_MOY%VAR $   �  @   a   F_RMS_1D_1C_MOY%MOY $   �  @   a   F_RMS_1D_1C_MOY%RMS '   #  �   a   F_RMS_1D_1C_MOY%WEIGHT    �  V       gen@SKEWNESS      V       gen@FLATNESS    [  R       gen@XMOMENT !   �  �       F_SKEWNESS_1D_1C %   g  <      F_SKEWNESS_1D_1C%SUM )   �  @      F_SKEWNESS_1D_1C%PRESENT &   �  =      F_SKEWNESS_1D_1C%SIZE %      �   a   F_SKEWNESS_1D_1C%VAR &   �  @   a   F_SKEWNESS_1D_1C%SKEW (   �  �   a   F_SKEWNESS_1D_1C%WEIGHT !   x  �       F_FLATNESS_1D_1C %   2  <      F_FLATNESS_1D_1C%SUM )   n  @      F_FLATNESS_1D_1C%PRESENT &   �  =      F_FLATNESS_1D_1C%SIZE %   �  �   a   F_FLATNESS_1D_1C%VAR &   w  @   a   F_FLATNESS_1D_1C%FLAT (   �  �   a   F_FLATNESS_1D_1C%WEIGHT    C  �       F_XMOM_1D_1C !   �  <      F_XMOM_1D_1C%SUM %   8  @      F_XMOM_1D_1C%PRESENT "   x  =      F_XMOM_1D_1C%SIZE "   �  �   a   F_XMOM_1D_1C%VAR1 "   A  �   a   F_XMOM_1D_1C%VAR2 "   �  @   a   F_XMOM_1D_1C%XMOM $     �   a   F_XMOM_1D_1C%WEIGHT 