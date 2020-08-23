      SUBROUTINE NRSPLINT(XPOS,YPOS,Y2A,N,X,Y,ERROR,*)

C#### Subroutine: NRSPLINE
C###  Description:
C###    Numerical recipes SPLINT (spline interpolation)

      IMPLICIT NONE
!     Parameter list
      INTEGER N
      REAL*8 X,XPOS(N),Y,YPOS(N),Y2A(N)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER k,khi,klo
      REAL*8 A,B,H

      CALL ENTERS('NRSPLINT',*9999)

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
        ERROR='Bad XPOS input in nrsplint'
        GOTO 9999
      ENDIF
      A=(XPOS(khi)-X)/H
      B=(X-XPOS(klo))/H
      Y=A*YPOS(klo)+B*YPOS(khi)+
     '  ((A**3-A)*Y2A(klo)+(B**3-B)*Y2A(khi))*(H**2)/6.d0


      CALL EXITS('NRSPLINT')
      RETURN
 9999 CALL ERRORS('NRSPLINT',ERROR)
      CALL EXITS('NRSPLINT')
      RETURN 1
      END


