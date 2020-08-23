      REAL*8 FUNCTION PAL2(K,na,XI)

C#### Function: PAL2
C###  Type: REAL*8
C###  Description:
C###    PAL2 evaluates 1D Legendre auxiliary basis functions of
C###    degree na, at Xi, for na>=2

      IMPLICIT NONE
!     Parameter List
      INTEGER K,na
      REAL*8 XI
!     Local Variables
      INTEGER N,NPOWER
      REAL*8 FACT,XIA,XINA
      LOGICAL EVEN

      XIA=2.0d0*XI-1.0d0
      NPOWER=na+1-K
      IF( DABS(XIA).LT.1.0d-08) THEN
        IF(NPOWER.LE.0) THEN
          XINA=1.0d0
        ELSE
          XINA=0.0d0
        ENDIF
      ELSE
        XINA=XIA**NPOWER
      ENDIF
      FACT=1.0d0
      DO N=1,na
        FACT=FACT*N
      ENDDO
      EVEN=.TRUE.
      IF((-1)**na.LT.0) EVEN=.FALSE.
      GO TO (10,20,30),K
 10   IF(EVEN) THEN
        PAL2=(XINA-1.0d0)/FACT
      ELSE
        PAL2=(XINA-XIA)/FACT
      ENDIF
      RETURN
 20   IF(EVEN) THEN
        PAL2=2.0d0*na/FACT*XINA
      ELSE
        PAL2=2.0d0/FACT*(na*XINA-1.0d0)
      ENDIF
      RETURN
 30   PAL2=4.0D0*na*(na-1)/FACT*XINA
      RETURN
      END


