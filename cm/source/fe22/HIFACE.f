      SUBROUTINE HIFACE(ISEG,ISFACE,ISFANO,STRING,ERROR,*)

C#### Subroutine: HIFACE
C###  Description:
C###    HIFACE hides face segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISFACE(NWM,NFM),ISFANO(NWM,NFM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),N3CO,nf,noiw,NTIW
      CHARACTER TYPE*8
      LOGICAL CBBREV

      CALL ENTERS('HIFACE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide faces
C###  Parameter:    <(segments/numbers)[segments]>
C###    Hide face segments or numbers
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to hide the
C###    faces on.
C###  Description:
C###    Hide face segments or numbers on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <(segments/numbers)[segments]>'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIFACE',ERROR,*9999)
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
              DO nf=1,NFT
                IF(ISEG(ISFACE(iw,nf)).EQ.2) THEN
                  CALL VISIB(iw,ISEG,ISFACE(iw,nf),'INVISIBLE',ERROR,
     '              *9999)
                ENDIF
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO
        ELSE IF(TYPE(1:7).EQ.'NUMBERS') THEN
          noco=noco+1
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO nf=1,NFT
                IF(ISEG(ISFANO(iw,nf)).EQ.2) THEN
                  CALL VISIB(iw,ISEG,ISFANO(iw,nf),'INVISIBLE',ERROR,
     '              *9999)
                ENDIF
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('HIFACE')
      RETURN
 9999 CALL ERRORS('HIFACE',ERROR)
      CALL EXITS('HIFACE')
      RETURN 1
      END


