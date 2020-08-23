      SUBROUTINE CONSTRUCT(ELEM,COORD,NODE,RESULT)

C#### Subroutine: Construct
C###  Description:
C###    Constructs the intercept of a plane (D) given the
C###    equation Ax + By + Cz + D = 0. The constant D determines whether
C###    the point (x,y,z) is inside (+ve) or outside (-ve) the
C###    tetrahedral.
CC JMB 16-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4)
      REAL*8 COORD(3),NODE(LDNODE,3),RESULT(4)
!     Local Variables
      INTEGER I,INFO,LDA,LWORK
      PARAMETER (LDA=4,LWORK=LDA*8)
      INTEGER IPIV(LDA)
      REAL*8 A(LDA,LDA),WORK(LWORK)

      DO I=1,LDA
        CALL DCOPY(NJT,NODE(ELEM(1,I),1),LDNODE,A(I,1),LDA)
        A(I,LDA)=ONE
      ENDDO !i
      ! Compute the inverse of the A matrix
      CALL DGETRF(LDA,LDA,A,LDA,IPIV,INFO)
c      CALL ASSERT_VORO(INFO .GE. 0, 'DGETRF failed in routine '//
c     '  'CONSTRUCT',ERROR,*9999)
      CALL DGETRI(LDA,A,LDA,IPIV,WORK,LWORK,INFO)
c      CALL ASSERT_VORO(INFO .EQ. 0, 'DGETRI failed in routine '//
c     '  'CONSTRUCT',ERROR,*9999)
      CALL DCOPY(NJT,COORD,1,WORK,1)
      WORK(LDA)=ONE
      CALL DGEMV('T',LDA,LDA,ONE,A,LDA,WORK,1,ZERO,RESULT,1)

      RETURN
      END


