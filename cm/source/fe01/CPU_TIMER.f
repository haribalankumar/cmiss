      SUBROUTINE CPU_TIMER(FLAG,CPU_TIME)

C#### Subroutine: CPU_TIMER
C###  Description:
C###    CPU_TIMER returns CPU time in CPU_TIME
C###    IF FLAG=CPU_USER the CPU time used while executing
C###    instructions in the user space of the calling process
C###    is returned in seconds since initialisation.
C###    IF FLAG=CPU_SYSTEM the CPU time used by the system
C###    on behalf of the calling process is returned in seconds
C###    since initialisation.
C###    IF FLAG=CPU_TOTAL the sum of the CPU time used while
C###    executing instructions in the user space of the
C###    calling process and CPU time used by the system
C###    on behalf of the calling process is returned in seconds
C###    since initialisation.
C###    IF FLAG=CPU_TICKS the routine returns the resolution of
C###    the cpu clock in ticks per second.

      IMPLICIT NONE
      INCLUDE 'time02.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER FLAG
      REAL CPU_TIME
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR
      LOGICAL INITIAL
      REAL*8 CURRENT_TIME,INITIAL_TIME
      CHARACTER ERROR*255
      SAVE INITIAL_TIME
      DATA INITIAL,INITIAL_TIME / .TRUE.,0.0d0 /

      ERR=0
C$OMP CRITICAL(CPU_TIMER_1)
      CALL CTIMER_CPU(CURRENT_TIME,FLAG,ERR,CERROR)
      IF(INITIAL) THEN
        INITIAL=.FALSE.
        INITIAL_TIME=CURRENT_TIME
      ENDIF
C$OMP END CRITICAL(CPU_TIMER_1)
      CPU_TIME=REAL(CURRENT_TIME-INITIAL_TIME)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,OP_STRING(1))
        OP_STRING(1)(CERRLEN:)=' '
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

 9999 RETURN
      END


