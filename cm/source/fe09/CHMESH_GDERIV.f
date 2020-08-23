      REAL*8 FUNCTION CHMESH_GDERIV(THETA)

C#### Function: CHMESH_GDERIV
C###  Type: REAL*8
C###  Description:
C###    CHMESH_GDERIV is one of J Crocombe's (temporary) change
C###    mesh functions.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'chmesh0.cmn'
      REAL*8 THETA
      INTEGER n

      CHMESH_GDERIV=0.0d0
      DO n=1,NCC
        CHMESH_GDERIV=CHMESH_GDERIV+(-COS_COEFFS(n)*DSIN(n*
     '    (THETA+PI*0.5d0)))
      ENDDO
      RETURN
      END


