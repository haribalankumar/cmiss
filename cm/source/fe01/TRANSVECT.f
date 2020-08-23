      SUBROUTINE TRANSVECT(TRANSMAT,VECTOR,ERROR,*)

C#### Subroutine: TRANSVECT
C###  Description:
C###    TRANSVECT applies the geometric transformation matrix set in
C###    SETTRANSMAT.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      REAL*8 TRANSMAT(3,4),VECTOR(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj1,nj2
      REAL*8 RESULT(3),SUM

      CALL ENTERS('TRANSVECT',*9999)

      DO nj1=1,NJT
        SUM=0.0d0
        DO nj2=1,NJT
          SUM=SUM+TRANSMAT(nj1,nj2)*VECTOR(nj2)
        ENDDO !nj
        RESULT(nj1)=SUM+TRANSMAT(nj1,NJT+1)
      ENDDO !nj
      DO nj1=1,NJT
        VECTOR(nj1)=RESULT(nj1)
      ENDDO !nj

      CALL EXITS('TRANSVECT')
      RETURN
 9999 CALL ERRORS('TRANSVECT',ERROR)
      CALL EXITS('TRANSVECT')
      RETURN 1
      END


