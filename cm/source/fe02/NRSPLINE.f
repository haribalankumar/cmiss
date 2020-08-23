      SUBROUTINE NRSPLINE(X,Y,N,YP1,YPN,Y2,ERROR,*)

C#### Subroutine: NRSPLINE
C###  Description:
C###    Numerical recipes SPLINE

      IMPLICIT NONE
!     Parameter list
      INTEGER N
      REAL*8 X(N),Y(N),YPN,YP1,Y2(N)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER NMAX
      PARAMETER (NMAX=100)
      INTEGER i,k
      REAL*8 P,QN,SIG,U(NMAX),UN

      CALL ENTERS('NRSPLINE',*9999)

      IF(YP1.GT.0.99d30) THEN
        Y2(1)=0.d0
        U(1)=0.d0
      ELSE
        Y2(1)=-0.5d0
        U(1)=(3.d0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO i=2,N-1
        SIG=(X(i)-X(i-1))/(X(i+1)-X(i-1))
        P=SIG*Y2(i-1)+2.d0
        Y2(i)=(SIG-1.d0)/P
        U(i)=(6.d0*((Y(i+1)-Y(i))/(X(i+1)-X(i))-(Y(i)-Y(i-1))
     '    /(X(i)-X(i-1)))/(X(i+1)-X(i-1))-SIG*U(i-1))/P
      ENDDO
      IF(YPN.GT.0.99d30) THEN
        QN=0.d0
        UN=0.d0
      ELSE
        QN=0.5d0
        UN=(3.d0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.d0)
      DO k=N-1,1,-1
        Y2(k)=Y2(k)*Y2(k+1)+U(k)
      ENDDO

      CALL EXITS('NRSPLINE')
      RETURN
 9999 CALL ERRORS('NRSPLINE',ERROR)
      CALL EXITS('NRSPLINE')
      RETURN 1
      END

