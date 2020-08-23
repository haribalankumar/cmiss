      SUBROUTINE GETRAN(IUNIT,TRANS)

C#### Subroutine: GETRAN
C###  Description:
C###    GETRAN reads the transformation matrix TRANS from IUNIT.

      IMPLICIT NONE
!     Parameter List
      INTEGER IUNIT
      REAL*8 TRANS(3,4)
!     Local Variables
      INTEGER i,j

C     CALL ENTERS('GETRAN',*9999)
      REWIND( IUNIT )
      READ(IUNIT,'(/)')
      DO i=1,3
        READ(IUNIT,'(15X,4F16.8)') (TRANS(i,j),j=1,4)
      ENDDO

C     CALL EXITS('GETRAN')
      RETURN
      END


