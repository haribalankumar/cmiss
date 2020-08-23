      REAL*8 FUNCTION CM_RANDOM_NUMBER(ISEED)

C#### Function: CM_RANDOM_NUMBER
C###  Type: REAL*8
C###  Description:
C###    CM_RANDOM_NUMBER returns REAL*8 random number between 0.0 and
C###    1.0
C###  Changes:
C###    26jun01.  DB.  Changed to using RAN0 below so that independent
C###      of compiler/os

      IMPLICIT NONE
!     Parameter List
      INTEGER ISEED
!     Functions
C DB.  26Jun01.  Don't use system random
C      REAL*8 RAND
      REAL*8 RAN0

C***  Fortran 90 has
C      CALL RANDOM_NUMBER(CM_RANDOM_NUMBER)
C DB.  26Jun01.  Don't use system random
C      CM_RANDOM_NUMBER=RAND(ISEED)
      CM_RANDOM_NUMBER=RAN0(ISEED)

      RETURN
      END


