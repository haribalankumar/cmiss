      SUBROUTINE CAPLOT(ISEG,ISPLOT,NEELEM,STRING,ERROR,*)

C#### Subroutine: CAPLOT
C###  Description:
C###    CAPLOT cancels plot PHIGS structure.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISPLOT(NHM,0:NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,ne,nh,noelem,noplot,nr

      CALL ENTERS('CAPLOT',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel plot;s
C###  Description:
C###    Cancel Phigs plot structure.

        OP_STRING(1)=STRING(1:IEND)//';s'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAPLOT',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        iw=9
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO nh=1,NHM
            DO noplot=1,NTPLOT
              IF(ISPLOT(nh,0,noplot).GT.0) THEN
                CALL DELETE_SEGMENT(ISPLOT(nh,0,noplot),ISEG,iw,ERROR,
     '            *9999)
                DO nr=1,NRT
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    IF(ISPLOT(nh,ne,noplot).GT.0) THEN
                      CALL DELETE_SEGMENT(ISPLOT(nh,ne,noplot),ISEG,iw,
     '                  ERROR,*9999)
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
        NTPLOT=0
      ENDIF

      CALL EXITS('CAPLOT')
      RETURN
 9999 CALL ERRORS('CAPLOT',ERROR)
      CALL EXITS('CAPLOT')
      RETURN 1
      END


