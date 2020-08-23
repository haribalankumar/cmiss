      SUBROUTINE GSVALUES(CONSTR,MN,SM,LDSM,WORK,ERROR,*)

C#### Subroutine: GSVALUES
C###  Description:
C###    GSVALUES computes the generalised singular values.
CC JMB 13-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER LDSM, MN
      REAL*8 SM(LDSM,*), WORK(*)
      CHARACTER CONSTR, ERROR*(*)
!     Local Variables
      INTEGER i
!     Functions
      LOGICAL LSAME

      CALL ENTERS('GSVALUES', *9999)

      IF( LSAME(CONSTR, 'I') ) THEN
        ! Zero-order Tikhonov or TSVD
        CALL DCOPY(MN, SM(1,1), 1, WORK(1), 1)
      ELSE
        ! First or second-order Tikhonov or TGSVD
        DO i = 1,MN
          WORK(i) = SM(MN + 1 - i,1)/SM(MN + 1 - i,2)
        ENDDO
      ENDIF

      CALL EXITS('GSVALUES')
      RETURN
 9999 CALL ERRORS('GSVALUES',ERROR)
      CALL EXITS('GSVALUES')
      RETURN 1
      END

C---------------------------------------------------------------------
