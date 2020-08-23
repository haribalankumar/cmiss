      SUBROUTINE HIFIEL(ISEG,ISFIEL,NEELEM,STRING,ERROR,*)

C#### Subroutine: HIFIEL
C###  Description:
C###    HIFIEL hides field segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISFIEL(NWM,NEM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),ne,noelem,noiw,nr,NTIW

      CALL ENTERS('HIFIEL',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide field
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to hide the
C###    field on.
C###  Description:
C###    Hide fields on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIFIEL',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(ISFIEL(iw,ne).GT.0) THEN
                  IF(ISEG(ISFIEL(iw,ne)).EQ.2) THEN
                    CALL VISIB(iw,ISEG,ISFIEL(iw,ne),'INVISIBLE',ERROR,
     '                *9999)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIFIEL')
      RETURN
 9999 CALL ERRORS('HIFIEL',ERROR)
      CALL EXITS('HIFIEL')
      RETURN 1
      END


