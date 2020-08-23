      SUBROUTINE HIOBJE(ISEG,ISOBJE,STRING,ERROR,*)

C#### Subroutine: HIOBJE
C###  Description:
C###    HIOBJE hides object segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISOBJE(NWM,NGRSEGM,NGRSEGM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,noobje,nopart,NTIW,NTPART

      CALL ENTERS('HIOBJE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide objects
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to hide the
C###    objects on.
C###  Description:
C###    Hide object segments on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIOBJE',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO noobje=1,NTOBJE
              NTPART=NSOBJE(2,noobje)
              DO nopart=1,NTPART
                IF(ISEG(ISOBJE(iw,noobje,nopart)).EQ.2) THEN
                  CALL VISIB(iw,ISEG,ISOBJE(iw,noobje,nopart),
     '              'INVISIBLE',ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIOBJE')
      RETURN
 9999 CALL ERRORS('HIOBJE',ERROR)
      CALL EXITS('HIOBJE')
      RETURN 1
      END


