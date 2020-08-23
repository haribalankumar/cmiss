      SUBROUTINE CAFACE(ISEG,ISFACE,ISFANO,STRING,ERROR,*)

C#### Subroutine: CAFACE
C###  Description:
C###    CAFACE cancels face segments.

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

      CALL ENTERS('CAFACE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel faces;s
C###  Description:
C###    Cancel face segment on specified workstations.
C###  Parameter:      <numbers>
C###    Specify the faces which are to be cancelled. The default
C###    is to cancel all faces.
C###    A list of face numbers or the 'all' option is allowed.
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the window number on which to cancel the faces.
C###    The default is to cancel the faces on all windows.

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15) //'<numbers>'
        OP_STRING(3)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAFACE',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        IF(noco.EQ.NTCO) THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO nf=1,NFT
                IF(ISFACE(iw,nf).GT.0) THEN
                  CALL DELETE_SEGMENT(ISFACE(iw,nf),ISEG,iw,ERROR,*9999)
                ENDIF
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF
        IF(ABBREV(CO(noco+1),'NUMBERS',1)) THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO nf=1,NFT
                IF(ISFANO(iw,nf).GT.0) THEN
                  CALL DELETE_SEGMENT(ISFANO(iw,nf),ISEG,iw,ERROR,*9999)
                ENDIF
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('CAFACE')
      RETURN
 9999 CALL ERRORS('CAFACE',ERROR)
      CALL EXITS('CAFACE')
      RETURN 1
      END


