      REAL*8 FUNCTION KDP(X1)

C#### Function: KDP
C###  Type: REAL*8
C###  Description:
C###    KDP calculates the elliptic integral of the first kind - K(m).

      IMPLICIT NONE
!     Parameter List
      REAL*8 X1
!     Local Variables
      REAL*8 A0,A1,A2,A3,A4,B0,B1,B2,B3,B4,TERM1,TERM2,X

      A0=1.38629436112D0
      A1=0.09666344259D0
      A2=0.03590092383D0
      A3=0.03742563713D0
      A4=0.01451196212D0

      B0=0.5D0
      B1=0.12498593597D0
      B2=0.06880248576D0
      B3=0.03328355346D0
      B4=0.00441787012D0

      X=1.0D0-X1
      TERM1=A0+(A1+(A2+(A3+A4*X)*X)*X)*X
      TERM2=B0+(B1+(B2+(B3+B4*X)*X)*X)*X
      KDP=TERM1+TERM2*DLOG(1.0D0/X)
      RETURN
      END


