      REAL*8 FUNCTION CHMESH_GFUNC(THETA)

C#### Function: CHMESH_GFUNC
C###  Type: REAL*8
C###  Description:
C###    CHMESH_GFUNC is one of J Crocombe's (temporary) change
C###    mesh functions.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'chmesh0.cmn'
      REAL*8 THETA
      INTEGER n

      CHMESH_GFUNC=1.0d0
      DO n=1,NCC
        CHMESH_GFUNC=CHMESH_GFUNC+(COS_COEFFS(n)*DCOS(n*(THETA+
     '    PI*0.5d0)))
      ENDDO

      RETURN
      END


