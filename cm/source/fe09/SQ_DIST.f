      REAL*8 FUNCTION SQ_DIST(NODE1,NODE2,XP)

C#### Function: SQ_DIST
C###  Type: REAL*8
C###  Description:
C###    Determines the squared distance between 2 nodes, version
C###    and deriv 1 only of the global CMISS array XP

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NODE1,NODE2
      REAL*8 XP(NKM,NVM,NJM,NPM)
!     Local Variables
      INTEGER nj

      SQ_DIST=0.0d0
      DO nj=1,NJT
        SQ_DIST=SQ_DIST+(XP(1,1,nj,NODE1)-XP(1,1,nj,NODE2))**2
      ENDDO

      RETURN
      END

C FE70 Functions
C ==============

