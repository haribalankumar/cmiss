      SUBROUTINE HINODE(ISEG,ISNONO,NPNODE,STRING,ERROR,*)

C#### Subroutine: HINODE
C###  Description:
C###    HINODE hides node segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISNONO(NWM,NPM),NPNODE(0:NP_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,nonode,np,nr,NTIW

      CALL ENTERS('HINODE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide nodes
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to hide
C###    nodes on.
C###  Description:
C###    Hide nodes on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HINODE',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nr=1,NRT
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                IF(ISNONO(iw,np).GT.0) THEN
                  IF(ISEG(ISNONO(iw,np)).EQ.2) THEN
                    CALL VISIB(iw,ISEG,ISNONO(iw,np),'INVISIBLE',ERROR,
     '                *9999)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HINODE')
      RETURN
 9999 CALL ERRORS('HINODE',ERROR)
      CALL EXITS('HINODE')
      RETURN 1
      END


