      REAL*8 FUNCTION EDP(X1)

C#### Function: EDP
C###  Type: REAL*8
C###  Description:
C###    EDP calculates the elliptic integral of the second kind - E(m).

      IMPLICIT NONE
!     Parameter List
      REAL*8 X1
!     Local Variables
      REAL*8 A1,A2,A3,A4,B1,B2,B3,B4,TERM1,TERM2,X

      A1=0.44325141463d0
      A2=0.06260601220d0
      A3=0.04757383546d0
      A4=0.01736506451d0

      B1=0.24998368310d0
      B2=0.09200180037d0
      B3=0.04069697526d0
      B4=0.00526449639d0

      X=1.0d0-X1
      TERM1=1.0d0+(A1+(A2+(A3+A4*X)*X)*X)*X
      TERM2=(B1+(B2+(B3+B4*X)*X)*X)*X
      EDP=TERM1+TERM2*DLOG(1.0d0/X)
      RETURN
      END


