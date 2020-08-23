C#######################################################################
      SUBROUTINE AMGTIME(t)
C#######################################################################

      IMPLICIT NONE

      REAL*8  t
      REAL*4  temp(2)
 
      REAL*4   ETIME
      EXTERNAL ETIME

C     The call is of the form
C 
C         t = ETIME(temp)
C
C     upon return temp(1) = "user time" and temp(2)= "system time"
C     The function value is the sum of the two.
C
      t = ETIME(temp)

      RETURN
      END
