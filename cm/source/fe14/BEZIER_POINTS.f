      SUBROUTINE BEZIER_POINTS(PL,XBEZ,YBEZ,ERROR,*)

C#### Subroutine: BEZIER_POINTS
C###  Description:
C###    BEZIER_POINTS calculates 21 points along Bezier curve given
C###    by XBEZ(i),YBEZ(i),i=1,4.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      REAL*8 PL(3,21),XBEZ(4),YBEZ(4)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,k
      REAL*8 XI,XI2,XI3

      CALL ENTERS('BEZIER_POINTS',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' XBEZ: '',4E12.3)') (XBEZ(i),i=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' YBEZ: '',4E12.3)') (YBEZ(i),i=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO j=0,20
        XI=DBLE(j)/20.0D0
        XI2=XI*XI
        XI3=XI2*XI
        PL(1,j+1)=(1.0D0-3.0D0*XI+3.0D0*XI2-XI3) * XBEZ(1) +
     '            (3.0D0*XI-6.0D0*XI2+3.0D0*XI3) * XBEZ(2) +
     '            (3.0D0*XI2-3.0D0*XI3         ) * XBEZ(3) +
     '             XI3                           * XBEZ(4)
        PL(2,J+1)=(1.0D0-3.0D0*XI+3.0D0*XI2-XI3) * YBEZ(1) +
     '            (3.0D0*XI-6.0D0*XI2+3.0D0*XI3) * YBEZ(2) +
     '            (3.0D0*XI2-3.0D0*XI3         ) * YBEZ(3) +
     '             XI3                           * YBEZ(4)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' J,XI= '',I3,F12.5)') j,XI
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' PL 1,2 = '',2F12.4)')
     '      (PL(k,j+1),k=1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO

      CALL EXITS('BEZIER_POINTS')
      RETURN
 9999 CALL ERRORS('BEZIER_POINTS',ERROR)
      CALL EXITS('BEZIER_POINTS')
      RETURN 1
      END


