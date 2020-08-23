      SUBROUTINE CALINE(ISEG,ISLINE,ISLINO,STRING,ERROR,*)

C#### Subroutine: CALINE
C###  Description:
C###    CALINE cancels line segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISLINE(NWM,2*NGRSEGM),ISLINO(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),N3CO,noiw,noline,NTIW
      CHARACTER TYPE*4
      LOGICAL CBBREV

      CALL ENTERS('CALINE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel lines;s
C###  Description:
C###    Cancel line segments.
C###  Parameter:      <(all/last)[all]>
C###    Specify either 'all' lines or just the 'last' line drawn.
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15) //'<(all/last)[all]>'
        OP_STRING(3)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CALINE',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'ALL',1,noco+1,NTCO,N3CO)) THEN
          TYPE='ALL'
        ELSE IF(CBBREV(CO,'LAST',1,noco+1,NTCO,N3CO)) THEN
          TYPE='LAST'
        ELSE
          TYPE='ALL'
        ENDIF

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(TYPE(1:3).EQ.'ALL') THEN
              DO noline=1,NTLINE
                IF(ISLINE(iw,noline).GT.0) THEN
                  CALL DELETE_SEGMENT(ISLINE(iw,noline),ISEG,iw,
     '              ERROR,*9999)
                ENDIF
              ENDDO
            ELSE IF(TYPE(1:4).EQ.'LAST') THEN
              IF(ISLINE(iw,NTLINE).GT.0) THEN
                CALL DELETE_SEGMENT(ISLINE(iw,NTLINE),ISEG,iw,
     '            ERROR,*9999)
              ENDIF
            ENDIF
            IF(ISLINO(iw).GT.0) THEN
              CALL DELETE_SEGMENT(ISLINO(iw),ISEG,iw,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
        IF(TYPE(1:3).EQ.'ALL') THEN
          NTLINE=0
        ELSE IF(TYPE(1:4).EQ.'LAST') THEN
          NTLINE=NTLINE-1
        ENDIF
      ENDIF

      CALL EXITS('CALINE')
      RETURN
 9999 CALL ERRORS('CALINE',ERROR)
      CALL EXITS('CALINE')
      RETURN 1
      END


