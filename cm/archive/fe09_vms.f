C#### Module: FE09_VMS
C###  Description:
C###    System specific real functions (Vax-VMS)

C###  Routine: CM_RANDOM_NUMBER (fn) return random number

      REAL*8 FUNCTION CM_RANDOM_NUMBER(ISEED)

C#### Function: CM_RANDOM_NUMBER
C###  Type: REAL*8
C###  Description:
C###    CM_RANDOM_NUMBER returns REAL*8 random number between 0.0 and 1.0

      IMPLICIT NONE
!     Parameter List
      INTEGER ISEED
!     Local Variables

      CM_RANDOM_NUMBER=DBLE(RAN(ISEED))

      RETURN
      END
