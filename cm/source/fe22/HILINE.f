      SUBROUTINE HILINE(ISEG,ISLINE,ISLINO,STRING,ERROR,*)

C#### Subroutine: HILINE
C###  Description:
C###    HILINE hides line segments.

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
      CHARACTER TYPE*8
      LOGICAL CBBREV

      CALL ENTERS('HILINE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide lines
C###  Description:
C###    Hide line segments or numbers.
C###  Parameter:    <(segments/numbers)[segments]>
C###    Option to hide segments or numbers
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.

        OP_STRING(1)=STRING(1:IEND) //' <(segments/numbers)[segments]>'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HILINE',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'SEGMENTS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='SEGMENTS'
        ELSE IF(CBBREV(CO,'NUMBERS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='NUMBERS'
        ELSE
          TYPE='SEGMENTS'
        ENDIF

        IF(TYPE(1:8).EQ.'SEGMENTS') THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO noline=1,NTLINE
                IF(ISLINE(iw,noline).GT.0) THEN
                  IF(ISEG(ISLINE(iw,noline)).EQ.2) THEN
                    CALL VISIB(iw,ISEG,ISLINE(iw,noline),'INVISIBLE',
     '                ERROR,*9999)
                  ENDIF
                ENDIF
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO
        ELSE IF(TYPE(1:7).EQ.'NUMBERS') THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              IF(ISEG(ISLINO(iw)).EQ.2) THEN
                CALL VISIB(iw,ISEG,ISLINO(iw),'INVISIBLE',ERROR,*9999)
              ENDIF
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('HILINE')
      RETURN
 9999 CALL ERRORS('HILINE',ERROR)
      CALL EXITS('HILINE')
      RETURN 1
      END


