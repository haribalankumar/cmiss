      SUBROUTINE HIELEM(ISEG,ISELNO,ISERR,NEELEM,STRING,ERROR,*)

C#### Subroutine: HIELEM
C###  Description:
C###    HIELEM hides element segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISELNO(NWM,NEM),ISERR(NWM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),N3CO,ne,noelem,noiw,nr,NTIW
      LOGICAL CBBREV,ERR

      CALL ENTERS('HIELEM',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide elements
C###  Description:
C###    Hide element numbers on specified workstations.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify which window to hide the element number on. The
C###    default is to hide the element numbers on all windows.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM hide elements error
C###  Description:
C###    Hide element errors on specified workstations.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify which window to hide the element errors on. The
C###    default is to hide the element errors on all windows.

        OP_STRING(1)=STRING(1:IEND) //' error <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIELEM',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

C CS 22/8/97 new hide element error
        IF(CBBREV(CO,'ERROR',2,noco+1,NTCO,N3CO)) THEN
          ERR=.TRUE.
        ELSE
          ERR=.FALSE.
        ENDIF

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(.NOT.ERR) THEN
                  IF(ISELNO(iw,ne).GT.0) THEN
                    IF(ISEG(ISELNO(iw,ne)).EQ.2) THEN
                      CALL VISIB(iw,ISEG,ISELNO(iw,ne),'INVISIBLE',
     '                  ERROR,*9999)
                    ENDIF
                  ENDIF
                ELSE
                  IF(ISERR(iw,ne).GT.0) THEN
                    IF(ISEG(ISERR(iw,ne)).EQ.2) THEN
                      CALL VISIB(iw,ISEG,ISERR(iw,ne),'INVISIBLE',
     '                  ERROR,*9999)
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIELEM')
      RETURN
 9999 CALL ERRORS('HIELEM',ERROR)
      CALL EXITS('HIELEM')
      RETURN 1
      END


