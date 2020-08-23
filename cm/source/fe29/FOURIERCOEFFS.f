      SUBROUTINE FOURIERCOEFFS(CONSTR,M,MN,U,LDU,B,WORK,ERROR,*)

C#### Subroutine: FOURIERCOEFFS
C###  Description:
C###    FOURIERCOEFFS computes the Fourier coefficients |U'b|.
CC JMB 13-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER LDU, M, MN
      REAL*8 B(*), U(LDU,*), WORK(*)
      CHARACTER CONSTR, ERROR*(*)
!     Local Variables
      INTEGER i
      REAL*8 ONE,ZERO
      PARAMETER (ONE=1.0d0,ZERO=0.0d0)
!     Functions
      LOGICAL LSAME

      CALL ENTERS('FOURIERCOEFFS', *9999)

      ! Initialisation
      CALL DGEMV('T', M, MN, ONE, U, LDU, B, 1, ZERO, WORK(1), 1)

      IF( .NOT.LSAME(CONSTR, 'I') ) THEN
        DO i = 1,MN/2
          CALL DSWAP(1, WORK(i), MN, WORK(MN + 1 - i), MN)
        ENDDO
      ENDIF

      CALL EXITS('FOURIERCOEFFS')
      RETURN
 9999 CALL ERRORS('FOURIERCOEFFS', ERROR)
      CALL EXITS('FOURIERCOEFFS')
      RETURN 1
      END

C---------------------------------------------------------------------
