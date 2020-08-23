      SUBROUTINE D_NORMALISE(NITB,D_VECTOR,VECTOR,ERROR,*)

C#### Subroutine: D_NORMALISE
C###  Description:
C###    Given a vector VECTOR and its NITB derivatives D_VECTOR,
C###    D_NORMALISE normalises the VECTOR by dividing the components of
C###    VECTOR by it's length and calculates the derivatives of the
C###    normalised vector.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NITB
      REAL*8 D_VECTOR(3,3),VECTOR(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni,njj
      REAL*8 SUM,VECTOR_LENGTH

      CALL ENTERS('D_NORMALISE',*9999)

      VECTOR_LENGTH=0.0d0
      DO njj=1,3
        VECTOR_LENGTH=VECTOR_LENGTH+VECTOR(njj)*VECTOR(njj)
      ENDDO !njj
      VECTOR_LENGTH=DSQRT(VECTOR_LENGTH)
      IF(VECTOR_LENGTH.GT.ZERO_TOL) THEN
        DO njj=1,3
          VECTOR(njj)=VECTOR(njj)/VECTOR_LENGTH
        ENDDO !njj
        DO ni=1,NITB
          SUM=0.0d0
          DO njj=1,3
            D_VECTOR(njj,ni)=D_VECTOR(njj,ni)/VECTOR_LENGTH
            SUM=SUM+D_VECTOR(njj,ni)*VECTOR(njj)
          ENDDO !njj
          DO njj=1,3
            D_VECTOR(njj,ni)=D_VECTOR(njj,ni)-SUM*VECTOR(njj)
          ENDDO !njj
        ENDDO !ni
      ELSE
        WRITE(OP_STRING,
     '    '('' >>WARNING: Cannot normalise a zero length vector'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('D_NORMALISE')
      RETURN
 9999 CALL ERRORS('D_NORMALISE',ERROR)
      CALL EXITS('D_NORMALISE')
      RETURN 1
      END


