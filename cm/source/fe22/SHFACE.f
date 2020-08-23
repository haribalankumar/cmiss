      SUBROUTINE SHFACE(ISEG,ISFACE,ISFANO,STRING,ERROR,*)

C#### Subroutine: SHFACE
C###  Description:
C###    SHFACE shows face segments.

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
      INTEGER IBEG,IEND,iw,IWK(6),nf,noiw,NTIW
      LOGICAL ABBREV

      CALL ENTERS('SHFACE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show faces
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to show
C###    faces on.
C###  Description:
C###    Make the face segments visible on the specified workstation.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM show face numbers
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to show the
C###    face numbers on.
C###  Description:
C###    Make the face segment numbers visible on the specified
C###    workstation.

        OP_STRING(1)=BLANK(1:15) //' numbers'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHFACE',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nf=1,NFT
              IF(ISEG(ISFACE(iw,nf)).EQ.1) THEN
                CALL VISIB(iw,ISEG,ISFACE(iw,nf),'VISIBLE',ERROR,*9999)
              ELSE IF(ISEG(ISFACE(iw,nf)).EQ.0) THEN
                WRITE(OP_STRING,'('' >>Face number '',I4,'
     '            //''' is not defined on '',I1)') nf,iw
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
        IF(ABBREV(CO(noco+1),'NUMBERS',1)) THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO nf=1,NFT
                IF(ISEG(ISFANO(iw,nf)).EQ.1) THEN
                  CALL VISIB(iw,ISEG,ISFANO(iw,nf),'VISIBLE',ERROR,
     '              *9999)
                ELSE IF(ISEG(ISFANO(iw,nf)).EQ.0) THEN
                  WRITE(OP_STRING,'('' >>Face number '',I4,'
     '              //''' is not defined on '',I1)') nf,iw
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('SHFACE')
      RETURN
 9999 CALL ERRORS('SHFACE',ERROR)
      CALL EXITS('SHFACE')
      RETURN 1
      END


