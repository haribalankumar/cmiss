      SUBROUTINE MAKE_SETS_DELAUNAY(BDRY,NPLIST,N_BDRY,N_IBDRY,XP,ERROR,
     '  *)

C#### Subroutine: MAKE_SETS_DELAUNAY
C###  Description:
C###  MAKE_SETS_DELAUNAY tags each boundary (B) node with the matching
C###  internal boundary (IB) node.  For corner or edge locations, there
C###  will be more than one B node for each IB node.  This routine
C###  assumes that the correct set of B nodes will be closer to the
C###  matching IB node than any other B nodes.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER BDRY(NBOUNDM),NPLIST(0:NPM),N_BDRY,N_IBDRY
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj,no_b,no_ib,np1,np2,n_set
      REAL*8 DIST,MIN_DIST

      CALL ENTERS('MAKE_SETS_DELAUNAY',*9999)

C For each IB node there will be one or more corresponding B node.
      DO no_b=1,N_BDRY
        np1=NPLIST(no_b)
        MIN_DIST=1.d10
        DO no_ib=1,N_IBDRY
          np2=NPLIST(N_BDRY+no_ib)
          DIST=0.d0
          DO nj=1,NJT
            DIST=DIST+(XP(1,1,nj,np1)-XP(1,1,nj,np2))**2.d0
          ENDDO !nj
          DIST=DSQRT(DIST)
          IF(DIST.LT.MIN_DIST)THEN
            MIN_DIST=DIST
            n_set=no_ib
          ENDIF
        ENDDO !IB
        CALL ASSERT(no_b.LE.NBOUNDM,'>>Increase NBOUNDM in genmesh.cmn',
     '    ERROR,*9999)
        BDRY(no_b)=n_set !records the IB node that B is matched with
      ENDDO !no_b

      CALL EXITS('MAKE_SETS_DELAUNAY')
      RETURN
 9999 CALL ERRORS('MAKE_SETS_DELAUNAY',ERROR)
      CALL EXITS('MAKE_SETS_DELAUNAY')
      RETURN 1
      END


