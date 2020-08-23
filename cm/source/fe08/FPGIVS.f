      SUBROUTINE FPGIVS(PIV,WW,COS,SIN)
C#### Subroutine: FPGIVS
C###  Description:
C###    Auxililary routine for CURFIT.
c  subroutine fpgivs calculates the parameters of a givens
c  transformation .
c  ..
      IMPLICIT NONE
c  ..scalar arguments..
      REAL*8 PIV,WW,COS,SIN
c  ..local scalars..
      REAL*8 DD,ONE,STORE
c  ..function references..
      REAL*8 ABS,SQRT
c  ..
      ONE = 0.1D+01
      STORE = ABS(PIV)
      IF(STORE.GE.WW) DD = STORE*SQRT(ONE+(WW/PIV)**2)
      IF(STORE.LT.WW) DD = WW*SQRT(ONE+(PIV/WW)**2)
      COS = WW/DD
      SIN = PIV/DD
      WW = DD
      RETURN
      END

