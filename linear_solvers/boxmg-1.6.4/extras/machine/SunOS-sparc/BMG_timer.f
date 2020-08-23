      SUBROUTINE BMG_timer(t)

C =======================================================================
C
C   BMG_timer.f
C
C   -------------------
C   DESCRIPTION:
C   -------------------
C
C     This subroutine uses the intrinsic function SECNDS to obtain
C     the system time as a real number accurate to 0.01 seconds. It is
C     computed as the elapsed time from midnight, minus the START time.
C
C     Sun Solaris:
C
C     There is no SECNDS (or equivalent) intrinsic function provided in 
C     SUN FORTRAN, it is provided as a VAX FORTRAN extension and hence
C     the $(LINK) line of the make file must include the VAX extension 
C     library -lV77 
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

      REAL*4    rZERO
      PARAMETER (rZERO=0)

      REAL*8   t

      REAL*4   SECNDS
      EXTERNAL SECNDS
 
C     The call is of the form
C 
C         t = SECNDS(START) 
C
C     however, for use with boxmg the value of START=0.0

      t = SECNDS(rZERO)

      RETURN
      END
