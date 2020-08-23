      SUBROUTINE BMG_timer(t)

C =======================================================================
C
C   BMG_timer.f
C
C   -------------------
C   DESCRIPTION:
C   -------------------
C
C     The intrinsic function ETIME is provided as an extension to the
C     standard intrinsic functions of most fortran compilers on Unix
C     systems.  For GNU compilers it is included in libg2c, and for 
C     most vendor compilers it is in libU77.
C 
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/02/15 (JDM)
C               - just a modified version of my earlier timer wrappers.
C
C ==========================================================================

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

