      SUBROUTINE FPROTA(COS,SIN,A,B)
C#### Subroutine: FPROTA
C###  Description:
C###    Auxililary routine for CURFIT.
c  subroutine fprota applies a givens rotation to a and b.
c  ..
      IMPLICIT NONE
c  ..scalar arguments..
      REAL*8 COS,SIN,A,B
c ..local scalars..
      REAL*8 STOR1,STOR2
c  ..
      STOR1 = A
      STOR2 = B
      B = COS*STOR2+SIN*STOR1
      A = COS*STOR1-SIN*STOR2
      RETURN
      END

