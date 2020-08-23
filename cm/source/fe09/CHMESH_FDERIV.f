      REAL*8 FUNCTION CHMESH_FDERIV(ZETA)

C#### Function: CHMESH_FDERIV
C###  Type: REAL*8
C###  Description:
C###    CHMESH_FDERIV is one of J Crocombe's (temporary) change
C###    mesh functions.

      IMPLICIT NONE
      INCLUDE 'chmesh0.cmn'
      REAL*8 ZETA
      INTEGER n

      CHMESH_FDERIV=0.0d0
      DO n=2,NPC
        CHMESH_FDERIV = CHMESH_FDERIV+POLY_COEFFS(n)*ZETA**(n-2)
      ENDDO
      RETURN
      END



