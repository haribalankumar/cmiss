      SUBROUTINE GAUSSLEG(X1,X2,X,W,N,ERROR,*)

C#### Subroutine: GAUSSLEG
C###  Description:
C###    This routine calculates the weights and abscissae for a
C###    Gauss-Legendre quadrature scheme. X1 an X2 are the lower and
C###    upper limits for the integration scheme. N is the order of
C###    the integration scheme. X are the positions and W are the
C###    weights of the gauss points.
C***    This is the routine given in Numerical Recipies.

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INTEGER NGMAX
      PARAMETER (NGMAX = 64)
!     Parameter List
      CHARACTER ERROR*(*)
      INTEGER N
      REAL*8 X1,X2,X(NGMAX),W(NGMAX)
!     Local Variables
      INTEGER i,j,m
      REAL*8 DIFF,DLAMCH,EPS,P1,P2,P3,PP,T1,T2,XL,XM,Z,Z1

      CALL ENTERS('GAUSSLEG',*9999)

      CALL ASSERT(N.LE.NGMAX,'>>ERROR: Max. order of quadrature is 64',
     '  ERROR,*9999)

      EPS=DLAMCH('EPS')
      M=(N+1)/2
      XM=0.5d0*(X2+X1)
      XL=0.5d0*(X2-X1)

      DO i=1,M
        Z=DCOS(DACOS(-1.0d0)*(i-0.25d0)/(N+0.5d0))
 1      CONTINUE
        P1=1.0d0
        P2=0.0d0
        DO j=1,N
          P3=P2
          P2=P1
          P1=((2.0d0*j-1.0d0)*Z*P2-(j-1.0d0)*P3)/J
        ENDDO
        PP=N*(Z*P1-P2)/(Z*Z-1.0d0)
        Z1=Z
        Z=Z1-P1/PP
        IF(DABS(Z-Z1).GT.EPS) GOTO 1
        X(i)=XM-XL*Z
        X(n+1-i)=XM+XL*Z
        W(i)=2.0d0*XL/((1.0d0-Z*Z)*PP*PP)
        W(N+1-i)=W(i)
      ENDDO

      IF(DOP) THEN !Check by integrating y=x+1
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        T1=0.0d0  !Numerical
        T2=0.0d0  !Analytic
        DO i=1,N
          T1=T1+((X(i)+1)*W(i))
        ENDDO
        T2=(X2**2.0d0/2.0d0+X2)-(X1**2.0d0/2.0d0-X1)
        DIFF=DABS(T2-T1)
        WRITE(OP_STRING,'(''Numerical Integration Test Diff: '',F12.6)')
     '    DIFF
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('GAUSSLEG')
      RETURN
 9999 CALL ERRORS('GAUSSLEG',ERROR)
      CALL EXITS('GAUSSLEG')
      RETURN 1
      END


c cpb 4/11/96 Adding logarithmic Gaussian quadrature

