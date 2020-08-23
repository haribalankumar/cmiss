      REAL*8 FUNCTION PL3(I,K,XI)

C#### Function: PL3
C###  Type: REAL*8
C###  Description:
C###    PL3 evaluates 1D cubic Lagrange basis funtion at XI.
C**** 08-Jan-1989: The following section added today (ADM)
C****              also corrected 2nd derivs in fn PL2 above
C****              I haven't checked the derivs very carefully!
C**** I is node position index
C**** K is partial derivative index
C**** XI is Xi-coordinate

      IMPLICIT NONE
!     Parameter List
      INTEGER I,K
      REAL*8 XI

      GO TO (10,20,30),K
 10     GO TO (11,12,13,14),I
 11       PL3=0.5D0*(3.0D0*XI-1.0D0)*(3.0D0*XI-2.0D0)*(1.0D0-XI)
          RETURN
 12       PL3=4.5D0*XI*(3.0D0*XI-2.0D0)*(XI-1.0D0)
          RETURN
 13       PL3=4.5D0*XI*(3.0D0*XI-1.0D0)*(1.0D0-XI)
          RETURN
 14       PL3=0.5D0*XI*(3.0D0*XI-1.0D0)*(3.0D0*XI-2.0D0)
          RETURN
 20     GO TO (21,22,23,24),I
 21       PL3=-13.5D0*XI*XI+18.0D0*XI-5.5D0
          RETURN
 22       PL3= 40.5D0*XI*XI-45.0D0*XI+9.0D0
          RETURN
 23       PL3=-40.5D0*XI*XI+36.0D0*XI-4.5D0
          RETURN
 24       PL3= 13.5D0*XI*XI- 9.0D0*XI+1.0D0
          RETURN
 30     GO TO (31,32,33,34),I
 31       PL3=9.0D0*(2.0D0-3.0D0*XI)
          RETURN
 32       PL3=9.0D0*(9.0D0*XI-5.0D0)
          RETURN
 33       PL3=9.0D0*(4.0D0-9.0D0*XI)
          RETURN
 34       PL3=9.0D0*(3.0D0*XI-1.0D0)
          RETURN
C**** 08-Jan-1989: End of addition (ADM)
      END


