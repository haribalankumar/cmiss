      SUBROUTINE INSERTION(NN,A,ERROR,*)

C#### Subroutine: INSERTION
C###  Description:
C###    Sorts an array A into ascending order, by straight insertion.

C***  Reference: Numerical Recipes in Fortran, Williham H. Press.
CC JMB 16-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER NN,A(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,TMP

      CALL ENTERS('INSERTION',*9999)

      DO J=2,NN
        TMP=A(J)
        DO I=J-1,1,-1
          IF(A(I).LE.TMP) GOTO 10
          A(I+1)=A(I)
        ENDDO !i
        I=0
  10    A(I+1)=TMP
      ENDDO !j

      CALL EXITS('INSERTION')
      RETURN
 9999 CALL ERRORS('INSERTION',ERROR)
      CALL EXITS('INSERTION')
      RETURN 1
      END


