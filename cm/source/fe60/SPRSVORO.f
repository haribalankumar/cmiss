      SUBROUTINE SPRSVORO(ISC_GKK,ISR_GKK,M,N,NFVC,nr,
     '  NVCNODE,nx,WORK_ARRAY,ERROR,*)

C#### Subroutine: SPRSVORO
C###  Description:
C###    Calculates the sparsity pattern for the Voronoi finite volume
C###    structure. This is a relatively simple task as the topology
C###    array NFVC is basically the same as the sparsity structure.

C***  Basically each Voronoi cell is linked to other Voronoi cells, and
C***  this 'connecting' information is stored in NFVC. This 'linking'
C***  is the same as compressed row storage, excepting that the
C***  diagonal, or the link to its own cell is not there.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER ISC_GKK(NISC_GKKM),ISR_GKK(NISR_GKKM),M,N,
     '  NFVC(2,0:NFVCM,NVCM),nr,NVCNODE(2,NP_R_M),
     '  nx
      LOGICAL*1 WORK_ARRAY(N,M)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nvc,nfvl,i,j,NZTOT,cnvc,cnonode

      CALL ENTERS('SPRSVORO',*9999)

C     ..Quick exit if possible
      CALL ASSERT(M.EQ.NVCT.AND.N.EQ.NVCT,'Number of Voronoi cells'//
     '  ' and rows,columns dont match',ERROR,*9999)

C     ..Initialize the work array sparsity pattern
      DO i=1,M
        DO j=1,N
          WORK_ARRAY(j,i)=.FALSE.
        ENDDO
      ENDDO

C     ..Create the work array sparsity pattern
      NZTOT=0
      DO nvc=1,M
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)
          IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
            cnvc=NVCNODE(MAP,cnonode)
            WORK_ARRAY(cnvc,nvc)=.TRUE.
            NZTOT=NZTOT+1
          ENDIF
        ENDDO
        WORK_ARRAY(nvc,nvc)=.TRUE. ! Diagonal - set entry
        NZTOT=NZTOT+1
      ENDDO

C     ..Check that number of non zeroes does not breach array bounds
      CALL ASSERT(NZTOT.LE.NZ_GKK_M,'>>Increase NZ_GKK_M',
     '  ERROR,*9999)

      NZZT(1,nr,nx)=NZTOT
      NZT(1,nx)=NZTOT

C     ..Calculate the sparsity pattern
      CALL CALC_SPARSE(NISC_GKKM,NISR_GKKM,ISC_GKK,ISR_GKK,M,N,
     '  NZTOT,SPARSEGKK(nx),WORK_ARRAY,ERROR,*9999)

      CALL EXITS('SPRSVORO')
      RETURN
 9999 CALL ERRORS('SPRSVORO',ERROR)
      CALL EXITS('SPRSVORO')
      RETURN 1
      END


