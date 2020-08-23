      SUBROUTINE INITVC(NPNODE,NRLIST,NVCNODE,VC,VC_init,ERROR,*)

C#### Subroutine: INITVC
C###  Description:
C###    INITVC copies the VC array into VC_INIT, for use in UPPRESSVORO

C     Alex 4Feb03 added the internal_boundary nodes
      
      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'

!     Parameter List
      INTEGER NRLIST(0:NRM),NVCNODE(2,NP_R_M),NPNODE(0:NP_R_M,0:NRM)
      REAL*8 VC(0:NVCM),VC_INIT(2,NVCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nonode,no_nr,nr,nvc
C      INTEGER nj_pressure
C      REAL*8 A,B,K,XP(NKM,NVM,NJM,NPM)
      LOGICAL ALL_REGIONS
      
      CALL ENTERS('INITVC',*9999)

      CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

      DO no_nr=1,NRLIST(0)
        nr=NRLIST(no_nr)
        DO nonode=1,NPNODE(0,nr)
          IF(NVCNODE(1,nonode).GE.2) THEN
            nvc=nvcnode(2,nonode)
            VC_init(1,nvc)=VC(nvc)
          ENDIF
        ENDDO
      ENDDO


      CALL EXITS('INITVC')
      RETURN
 9999 CALL ERRORS('INITVC',ERROR)
      CALL EXITS('INITVC')
      RETURN 1
      END



