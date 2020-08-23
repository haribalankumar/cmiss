      INTEGER FUNCTION CALC_SAMPLE_FROM_TIME(TIME,ERR,ERROR)

C#### Function: CALC_SAMPLE_FROM_TIME
C###  Description:
C###    CALC_SAMPLE_FROM_TIME calculates the sample number corresponding to
C###    a real time using the transfer frequency. Note that a time of 0.0
C###    corresponds to a sample number of 1.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER ERR
      REAL*8 TIME
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER SAMPLE

      SAMPLE=INT(TIME*TRSF_FREQUENCY)+1
      IF(SAMPLE.LT.1) THEN
        ERR=1
        ERROR='>>Zero or negative sample number'
      ELSE IF(SAMPLE.GT.NTSM) THEN
        ERR=1
        ERROR='>>Increase NTSM'
      ENDIF

      CALC_SAMPLE_FROM_TIME=SAMPLE

      RETURN
      END


