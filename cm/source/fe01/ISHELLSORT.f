      SUBROUTINE ISHELLSORT(N,A)

C#### Subroutine: ISHELLSORT
C###  Description:
C###    Numerical recipes integer shell sort algorithm. Use when N < 50

      IMPLICIT NONE
!     Parameter List
      INTEGER N,A(N)
!     Local Variables
      INTEGER I,J,INC,V

      INC=1
 1    INC=3*INC+1
      IF(INC.LE.N) GOTO 1
 2    CONTINUE
      INC=INC/3
      DO I=INC+1,N
        V=A(I)
        J=I
 3      IF(A(J-INC).GT.V)THEN
          A(J)=A(J-INC)
          J=J-INC
          IF(J.LE.INC) GOTO 4
          GOTO 3
        ENDIF
 4      A(J)=V
      ENDDO
      IF(INC.GT.1) GOTO 2
      RETURN
      END

