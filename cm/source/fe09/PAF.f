      REAL*8 FUNCTION PAF(K,na,XI)

C#### Function: PAF
C###  Type: REAL*8
C###  Description:
C###    PAF evaluates 1D Fourier auxiliary basis (sine) functions at Xi.

C**** na is half the wavenumber of the function
C**** K is the derivative reqd. (K=1 -> value reqd)

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
!     Parameter List
      INTEGER K,na
      REAL*8 XI
!     Local Variables
      REAL*8 PAL0,RNA

      IF(na.EQ.0) THEN        !Use constant basis fn
        PAF=PAL0(K)
      ELSE
        RNA=DBLE(na)
        GO TO (10,20,30),K
 10       PAF=               DSIN(RNA*PI*XI)
          RETURN
 20       PAF=        RNA*PI*DCOS(RNA*PI*XI)
          RETURN
 30       PAF=-RNA*RNA*PI*PI*DSIN(RNA*PI*XI)
          RETURN
      ENDIF

      RETURN
      END


