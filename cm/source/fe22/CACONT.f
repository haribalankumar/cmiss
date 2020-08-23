      SUBROUTINE CACONT(ISCONO,ISCONT,ISEG,NEELEM,STRING,ERROR,*)

C#### Subroutine: CACONT
C###  Description:
C###    CACONT cancels element contour segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISCONO(NHM,NEM),ISCONT(NHM,NEM,NGRSEGM),
     '  NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,ne,nh,nocont,noelem,nr

      CALL ENTERS('CACONT',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel contours;s
C###  Description:
C###    Cancel contour segment.

        OP_STRING(1)=STRING(1:IEND)//';s'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CACONT',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        iw=2*NJT-3+IMAP
        CALL ACWK(iw,1,ERROR,*9999)
        DO nh=1,NHM
          DO nr=1,NRT
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nocont=1,NTCONT
                IF(ISCONT(nh,ne,nocont).GT.0) THEN
                  CALL DELETE_SEGMENT(ISCONT(nh,ne,nocont),ISEG,iw,
     '              ERROR,*9999)
                ENDIF
              ENDDO
              IF(ISCONO(nh,ne).GT.0) THEN
                CALL DELETE_SEGMENT(ISCONO(nh,ne),ISEG,iw,ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        CALL DAWK(iw,1,ERROR,*9999)
      ENDIF

      CALL EXITS('CACONT')
      RETURN
 9999 CALL ERRORS('CACONT',ERROR)
      CALL EXITS('CACONT')
      RETURN 1
      END


