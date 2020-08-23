      SUBROUTINE OPTIME(NTIME_INTERP,NTIME_POINTS,TIME_VALUES,
     '  TIME_VARIABLE_NAMES,ERROR,*)

C#### Subroutine: OPTIME
C###  Description:
C###    OPTIME outputs time variable information
C**** Created by Martin Buist, 15 November 1999.

      IMPLICIT NONE

      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'

!     Parameter List
      INTEGER NTIME_INTERP(NTIMEVARSM),NTIME_POINTS(NTIMEVARSM)
      REAL*8 TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM)
      CHARACTER ERROR*(*),TIME_VARIABLE_NAMES(NTIMEVARSM)*(*)
!     Local Variables
      INTEGER i,IBEG,IEND,j

      CALL ENTERS('OPTIME',*9999)

      DO i=1,NTIMEVARST
        WRITE(OP_STRING,'('' '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        CALL STRING_TRIM(TIME_VARIABLE_NAMES(i),IBEG,IEND)
        WRITE(OP_STRING,'('' Time variable name: '',A)')
     '    TIME_VARIABLE_NAMES(i)(IBEG:IEND)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        WRITE(OP_STRING,'('' Number of time points set: '',I4)')
     '    NTIME_POINTS(i)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C       FORMAT='($,'' Enter the time variable value before '
C    '    //'time point 1 [0.0]: '',E12.4)'
        WRITE(OP_STRING,'('' Value before time point 1: '',E12.4)')
     '    TIME_VALUES(2,0,i)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        DO j=1,NTIME_POINTS(i)
          WRITE(OP_STRING,'('' Time point '',I4,'' is: '',E12.4)')
     '      j,TIME_VALUES(1,j,i)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          WRITE(OP_STRING,'('' Value at time point '',
     '      I4,'' is: '',E12.4)') j,TIME_VALUES(2,j,i)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !j

        WRITE(OP_STRING,'('' Value after time point '',I4,'': '',
     '    E12.4)') NTIME_POINTS(i),TIME_VALUES(2,NTIME_POINTS(i)+1,i)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(NTIME_INTERP(i).EQ.1) THEN
          WRITE(OP_STRING,'('' Interpolation is Constant'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(NTIME_INTERP(i).EQ.2) THEN
          WRITE(OP_STRING,'('' Interpolation is Linear Lagrange'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(NTIME_INTERP(i).EQ.3) THEN
          WRITE(OP_STRING,'('' Interpolation is Quadratic Lagrange'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(NTIME_INTERP(i).EQ.4) THEN
          WRITE(OP_STRING,'('' Interpolation is Cubic Lagrange'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(NTIME_INTERP(i).EQ.5) THEN
          WRITE(OP_STRING,'('' Interpolation is Cubic Hermite'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF !interpolation
      ENDDO !i

      CALL EXITS('OPTIME')
      RETURN
 9999 CALL ERRORS('OPTIME',ERROR)
      CALL EXITS('OPTIME')
      RETURN 1
      END


C LKC 11-NOV-1998 *** archived ***
C      SUBROUTINE OPTRAN(ERROR,*)




