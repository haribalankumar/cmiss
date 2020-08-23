      REAL*8 FUNCTION CHMESH_FFUNC(ZETA)

C#### Function: CHMESH_FFUNC
C###  Type: REAL*8
C###  Description:
C###    CHMESH_FFUNC is one of J Crocombe's (temporary) change
C###    mesh functions.

      IMPLICIT NONE
      INCLUDE 'chmesh0.cmn'
      REAL*8 ZETA
      INTEGER n

      CHMESH_FFUNC=0.0d0
      DO n=1,NPC
        CHMESH_FFUNC = CHMESH_FFUNC+POLY_COEFFS(n)*ZETA**(n-1)
      ENDDO
      RETURN
      END


