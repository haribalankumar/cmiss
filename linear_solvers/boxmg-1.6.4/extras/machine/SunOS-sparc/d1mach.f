      REAL*8 FUNCTION D1MACH(I)

C***BEGIN PROLOGUE  D1MACH
C***DATE WRITTEN   790101   (YYMMDD)
C***REVISION DATE  860825   
C                  940420
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns double precision machine dependent constants
C***DESCRIPTION
C
C     D1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          A = D1MACH(I)
C
C     where I=1,...,5.  The (output) value of A above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  Double-Precision Machine Constants
C
C  D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
C  D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C  D1MACH(3) = B**(-T), the smallest relative spacing.
C  D1MACH(4) = B**(1-T), the largest relative spacing.
C  D1MACH(5) = LOG10(B)
C
C***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR
C                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
C                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
C                 PP. 177-188.
C***ROUTINES CALLED (NONE)
C***END PROLOGUE  D1MACH

      IMPLICIT NONE

      INTEGER  I
      REAL*8   DMACH(5)

C -------------------------------------------------------------------------
C     Machine Constants for the SUNs SPARC series
C
      DATA DMACH(1) / 2.2250738585072104D-308  /
      DATA DMACH(2) / 1.7976931348623131D+308  /
      DATA DMACH(3) / 1.1102230246251565D-16   /
      DATA DMACH(4) / 2.2204460492503131D-16   /
      DATA DMACH(5) / 3.0102999566398120D-01   /

C ------------------------------------------------------------------------
C  FIRST EXECUTABLE STATEMENT  D1MACH
C
      IF (I .LT. 1  .OR.  I .GT. 5) THEN
         WRITE(*,*) '***ERROR:  D1MACH -- I out of bounds***'
         STOP
      ENDIF

      D1MACH = DMACH(I)

      RETURN
      END
