      REAL*8 FUNCTION PSIM(I,K,M,XL)

C#### Function: PSIM
C###  Type: REAL*8
C###  Description:
C###    PSIM evaluates contribution to Simplex basis function at XL.
C**** I is node position index
C**** K is partial derivative index for area coords
C**** M is 1,2,3 for linear,quad or cubic
C**** XL is area coordinate

      IMPLICIT NONE
!     Parameter List
      INTEGER I,K,M
      REAL*8 XL
!     Local Variables
      INTEGER J
      REAL*8 PSIM0,SUM

      PSIM=1.0D0
      DO J=1,I
        PSIM=PSIM*(DBLE(M)*XL-DBLE(J)+1.0D0)/DBLE(J)
      ENDDO
      IF(K.EQ.1) RETURN

      PSIM0=PSIM
      SUM=0.0D0
      DO J=1,I
        SUM=SUM+DBLE(M)/(DBLE(M)*XL-DBLE(J)+1.0D0)
      ENDDO
      PSIM=PSIM*SUM
      IF(K.EQ.2) RETURN

      SUM=0.0D0
      DO J=1,I
        SUM=SUM+(DBLE(M)/(DBLE(M)*XL-DBLE(J)+1.0D0))**2
      ENDDO
      PSIM=PSIM*PSIM/PSIM0-PSIM0*SUM
      RETURN
      END


