      REAL*8 FUNCTION CALC_TIME_FROM_SAMPLE(SAMPLE)

C#### Function: CALC_TIME_FROM_SAMPLE
C###  Description:
C###    CALC_TIME_FROM_SAMPLE calculates the real time corresponding to
C###    a sample number using the transfer frequency. Note that a time
C###    of 0.0 corresponds to a sample number of 1.

      IMPLICIT NONE
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER SAMPLE
!     Local Variables

      CALC_TIME_FROM_SAMPLE=DBLE(SAMPLE-1)/TRSF_FREQUENCY

      RETURN
      END


