      SUBROUTINE NORMALISE(NUMCMPTS,VECTOR,ERROR,*)

C#### Subroutine: NORMALISE
C###  Description:
C###    NORMALISE divides the components of VECTOR by it's length

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NUMCMPTS
      REAL*8 VECTOR(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni
      REAL*8 VECTOR_LENGTH

      CALL ENTERS('NORMALISE',*9999)

      VECTOR_LENGTH=0.0d0
      DO ni=1,NUMCMPTS
        VECTOR_LENGTH=VECTOR_LENGTH+VECTOR(ni)*VECTOR(ni)
      ENDDO !ni
      VECTOR_LENGTH=DSQRT(VECTOR_LENGTH)
      IF(VECTOR_LENGTH.GT.ZERO_TOL) THEN
        DO ni=1,NUMCMPTS
          VECTOR(ni)=VECTOR(ni)/VECTOR_LENGTH
        ENDDO !ni
      ELSE
        WRITE(OP_STRING,
     '    '('' >>WARNING: Cannot normalise a zero length vector'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('NORMALISE')
      RETURN
 9999 CALL ERRORS('NORMALISE',ERROR)
      CALL EXITS('NORMALISE')
      RETURN 1
      END


