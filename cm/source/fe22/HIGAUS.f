      SUBROUTINE HIGAUS(ISEG,ISGAUS,NBJ,NEELEM,STRING,ERROR,*)

C#### Subroutine: HIGAUS
C###  Description:
C###    HIGAUS hides Gauss segments.


      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISGAUS(NWM,NGM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),ne,ng,noelem,noiw,nr,NTIW

      CALL ENTERS('HIGAUS',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide gauss
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to hide the
C###    gauss on.
C###  Description:
C###    Hide Gauss points on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIGAUS',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                DO ng=1,NGT(NBJ(1,ne))
                  IF(ISEG(ISGAUS(iw,ng,ne)).EQ.2) THEN
                    CALL VISIB(iw,ISEG,ISGAUS(iw,ng,ne),'INVISIBLE',
     '                ERROR,*9999)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIGAUS')
      RETURN
 9999 CALL ERRORS('HIGAUS',ERROR)
      CALL EXITS('HIGAUS')
      RETURN 1
      END


