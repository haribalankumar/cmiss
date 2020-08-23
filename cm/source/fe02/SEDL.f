      SUBROUTINE SEDL(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,DL,SE,ERROR,*)

C#### Subroutine: SEDL
C###  Description:
C###    SEDL gets line scale factors from element scale factors.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),nb,
     '  NEELEM(0:NE_R_M,0:NRM),NLL(12,NEM),NNL(0:4,12,NBFM),
     '  NPL(5,0:3,NLM)
      REAL*8 DL(3,NLM),SE(NSM,NBFM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER N,NAE,ne,ni,NI1(5,3),ni2,ni3,NITB,nk,NKTB,nl,nn,
     '  NNK,noelem,nr

      DATA NI1/1,0,0,0,0,
     '         1,2,1,0,0,
     '         1,2,3,1,2/

      CALL ENTERS('SEDL',*9999)
      NITB=NIT(nb)
      NKTB=NKT(0,nb)
C KAT 1/3/00: ni3 needs initializing if not 3d
      ni3=0
      DO nr=1,NRT
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NKTB.GT.1) THEN
            DO NAE=1,NLE(nb)
              nl=NLL(NAE,ne)
              IF(nl.GT.0) THEN
                ni=NPL(1,0,nl)
                IF(IBT(1,ni,nb).EQ.2) THEN
C ***             Hermite
                  ni2=NI1(ni+1,NITB)
                  IF(NITB.EQ.3) ni3=NI1(ni+2,NITB)
                  DO N=1,2
                    nn=NNL(N,NAE,nb)
                    NNK=(nn-1)*NKTB
                    DO nk=2,NKTB
                      IF(NITB.EQ.1.OR.
     '                  (NITB.EQ.2.AND.IDO(nk,nn,ni2,nb).EQ.1).OR.
     '                  (NITB.EQ.3.AND.IDO(nk,nn,ni2,nb).EQ.1
     '                  .AND.IDO(nk,nn,ni3,nb).EQ.1)) THEN
                        DL(N,nl)=SE(nk+NNK,nb,ne)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
            ENDDO

          ENDIF
        ENDDO
      ENDDO

      CALL EXITS('SEDL')
      RETURN
 9999 CALL ERRORS('SEDL',ERROR)
      CALL EXITS('SEDL')
      RETURN 1
      END


