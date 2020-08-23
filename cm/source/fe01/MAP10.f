      SUBROUTINE MAP10(NTPTS,Y)

C#### Subroutine: MAP10
C###  Description:
C###    MAP10 sets world coords for history display.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'hist00.cmn'
!     Parameter List
      INTEGER NTPTS
      REAL*8 Y(NTPTS)
!     Local Variables
      CHARACTER ERROR*10
      INTEGER n

C     CALL ENTERS('MAP10',*9999)
      IF(DABS(YDMAX).GT.1.0D-6) THEN
        DO n=1,NTPTS
          Y(n)=(Y(n)-(YDMAX+YDMIN)/2.0D0)*2.0D0*FACTOR/(YDMAX-YDMIN)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,*)' map10 y=',y(n)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDDO
      ENDIF

 9999 RETURN
      END


