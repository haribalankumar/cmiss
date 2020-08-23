      SUBROUTINE HIFIBR(ISEG,ISFIBR,NEELEM,STRING,ERROR,*)

C#### Subroutine: HIFIBR
C###  Description:
C###    HIFIBR hides fibre segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISFIBR(NWM,NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),ne,noelem,nofibr,noiw,nr,NTIW

      CALL ENTERS('HIFIBR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide fibres
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to hide the
C###    fibres on.
C###  Description:
C###    Hide fibres on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIFIBR',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nofibr=1,NTFIBR
              DO nr=1,NRT
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IF(ISEG(ISFIBR(iw,ne,nofibr)).EQ.2) THEN
                    CALL VISIB(iw,ISEG,ISFIBR(iw,ne,nofibr),'INVISIBLE',
     '                ERROR,*9999)
                  ENDIF
                ENDDO
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIFIBR')
      RETURN
 9999 CALL ERRORS('HIFIBR',ERROR)
      CALL EXITS('HIFIBR')
      RETURN 1
      END


