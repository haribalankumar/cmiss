      INTEGER FUNCTION ILOG10(RDATA)

C#### Function: ILOG10
C###  Type: INTEGER
C###  Description:
C###    ILOG10eRturns the truncated log to base 10 of the real
C###    argument ABS(RDATA).

      IMPLICIT NONE
!     Parameter List
      REAL*8 RDATA
!     Local Variables
      INTEGER i
      REAL*8 ARDATA,R10

C     CALL ENTERS('ILOG10',*9999)
      ARDATA=DABS(RDATA)
      IF(ARDATA.GE.1.0D0) THEN
        R10=1.0D1
        DO 1 i=0,100
          IF(ARDATA.LT.R10) THEN
            ILOG10=i
            GOTO 3
          ENDIF
          R10=R10*1.0D1
    1   CONTINUE
      ELSE
        R10=1.0D-1
        DO 2 i=1,100
          IF(ARDATA.GE.R10) THEN
            ILOG10=-i
            GOTO 3
          ENDIF
          R10=R10*1.0D-1
    2   CONTINUE
      ENDIF
C3    CALL EXITS('ILOG10')
 3    RETURN
      END


