      SUBROUTINE CASTRA(ISEG,ISSTRA,NEELEM,STRING,ERROR,*)

C#### Subroutine: CASTRA
C###  Description:
C###    CASTRA cancels element strain segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSTRA(NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,ne,noelem,nostra,nr

      CALL ENTERS('CASTRA',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel strain;s
C###  Description:
C###    Cancel strain segment.

        OP_STRING(1)=STRING(1:IEND)//';s'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CASTRA',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        iw=2*NJT-3
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO nr=1,NRT
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nostra=1,NTSTRA
                IF(ISSTRA(ne,nostra).GT.0) THEN
                  CALL DELETE_SEGMENT(ISSTRA(ne,nostra),ISEG,iw,ERROR,
     '              *9999)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)
          NTSTRA=0
        ENDIF
      ENDIF

      CALL EXITS('CASTRA')
      RETURN
 9999 CALL ERRORS('CASTRA',ERROR)
      CALL EXITS('CASTRA')
      RETURN 1
      END


