      SUBROUTINE NRSPLINTGRAD(XPOS,YPOS,Y2A,N,X,Y,ERROR,*)

C#### Subroutine: NRSPLINE
C###  Description:
C###    Numerical recipes SPLINTGRAD (spline gradient interpolation)

      IMPLICIT NONE
!     Parameter list
      INTEGER N
      REAL*8 X,XPOS(N),Y,YPOS(N),Y2A(N)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER k,khi,klo
      REAL*8 A,B,H

      CALL ENTERS('NRSPLINTGRAD',*9999)

      klo=1
      khi=N
 1    IF(khi-klo.GT.1) THEN
        k=(khi+klo)/2
        IF(XPOS(k).GT.X) THEN
          khi=k
        ELSE
          klo=k
        ENDIF
        GOTO 1
      ENDIF
      H=XPOS(khi)-XPOS(klo)
      IF(H.EQ.0.d0) THEN
        ERROR='Bad XPOS input in nrsplintgrad'
        GOTO 9999
      ENDIF
      A=(XPOS(khi)-X)/H
      B=(X-XPOS(klo))/H
      Y=(YPOS(khi)-YPOS(klo))/H
     '  -((3.d0*A*A-1.d0)*H*Y2A(klo) - (3.d0*B*B-1.d0)*H*Y2A(khi))/6.d0

      CALL EXITS('NRSPLINTGRAD')
      RETURN
 9999 CALL ERRORS('NRSPLINTGRAD',ERROR)
      CALL EXITS('NRSPLINTGRAD')
      RETURN 1
      END


