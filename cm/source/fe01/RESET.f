      SUBROUTINE RESET(TRANS)

C#### Subroutine: RESET
C###  Description:
C###    RESET resets the transformation matrix TRANS to the identity
C###    matrix.

      IMPLICIT NONE
!     Parameter List
      REAL*8 TRANS(3,4)
!     Local Variables
      INTEGER i,j

C     CALL ENTERS('RESET',*9999)
      DO i=1,3
        DO j=1,4
          TRANS(i,j)=0.0D0
        ENDDO
        TRANS(i,i)=1.0D0
      ENDDO

C     CALL EXITS('RESET')
      RETURN
      END


