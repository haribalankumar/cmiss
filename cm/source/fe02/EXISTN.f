      LOGICAL FUNCTION EXISTN(np,NPNODE,XP,XC)

C#### Function: EXISTN
C###  Type: LOGICAL
C###  Description:
C###    EXISTN returns true if x,y world coords lie inside a
C###    node segment.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER np,NPNODE(0:NP_R_M,0:NRM)
      REAL*8 XC(3),XP(NKM,NVM,NJM,NPM)
!     Local Variables
      INTEGER MP,nonode,nj
      REAL*8 DIST

      EXISTN=.FALSE.
      DO nonode=1,NPNODE(0,1)
        MP=NPNODE(nonode,1)
        DIST=0.0d0
        DO nj=1,NJT
          DIST=DIST+(XC(nj)-XP(1,1,nj,MP))**2
        ENDDO
        DIST=DSQRT(DIST)
        IF(DIST.LT.(1.0d0-2.0d0*DBLE(DIAG))) THEN
          np=MP
          EXISTN=.TRUE.
          GO TO 51
        ENDIF
      ENDDO

 51   RETURN
      END


