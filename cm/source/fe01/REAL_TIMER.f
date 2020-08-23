      SUBROUTINE REAL_TIMER(FLAG,REAL_TIME)

C#### Subroutine: REAL_TIMER
C###  Description:
C###     Returns in REAL_TIME(1) the value of time in seconds since
C###     the initial call to the routine.
C###
C###     In the past it returned the time since the Unix epoch,
C###     00:00:00 UTC, January 1,1970 (sort of like Cambodia's year 1,
C###     except it's Unix so it is year 0). However, the mantissa of a
C###     REAL*4 is not big enough to hold this value without rounding,
C###     (at least it hasn't been since mid Feb 1970) so now we use the
C###     routine to return the time since the initial call.

      IMPLICIT NONE
      INCLUDE 'time02.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER FLAG
      REAL REAL_TIME
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR
      LOGICAL INITIAL
      REAL*8 RETURN_TIME,INITIAL_TIME
      CHARACTER ERROR*255
      SAVE INITIAL,INITIAL_TIME
      DATA INITIAL,INITIAL_TIME / .TRUE.,0.0D0 /

      ERR=0
C$OMP CRITICAL(REAL_TIMER_1)
      CALL CTIMER_REAL(RETURN_TIME,FLAG,ERR,CERROR)
      IF(INITIAL) THEN
        INITIAL=.FALSE.
        INITIAL_TIME=RETURN_TIME
      ENDIF
      REAL_TIME = REAL(RETURN_TIME-INITIAL_TIME)
C$OMP END CRITICAL(REAL_TIMER_1)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,OP_STRING(1))
        OP_STRING(1)(CERRLEN:)=' '
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

 9999 RETURN
      END


