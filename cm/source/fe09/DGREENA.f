      REAL*8 FUNCTION DGREENA(IGREN,CE,XN,RP,ZP,RQ,ZQ)

C#### Function: DGREENA
C###  Type: REAL*8
C###  Description:
C###    DGREENA identifies appropriate normal derivative of complex
C###    Green's function.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'eqt000.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IGREN
      REAL*8 CE(NMM),XN(NJM),ZP,RP,ZQ,RQ
!     Local Variables
      INTEGER ngp
      REAL*8 A,B,DRDNU,EDP,KB1,KB2,KDP,M,POWER,R,THETA

      IF(IGREN.EQ.12) THEN
        !Normal deriv. Yukawa equation,3d cylindrically symmetric
        !Reference: Brebbia, Telles and Wrobel pg 97
        A = RP**2 + RQ**2 + (ZP-ZQ)**2
        B = 2.0D0 * RP * RQ
        M = 2.0D0 * B/(A+B)
        !Evaluate KB1 using quadrature
        KB1=0.0D0
        !Loop over the Gauss points
        DO ngp=1,ng
          THETA=2.0D0*PI*D1(ngp)
          !Evaluate R(P,Q)
          R=DSQRT(A-B*DCOS(THETA))
          !Evaluate the exponent of the exponential term
          POWER=-1.0D0*CE(1)*R
          !Evaluate the normal derivative of R(P,Q)
          !DRDNU = grad(R).nu = drdzq*nuzq + drdrq*nurq
          DRDNU=((ZQ-ZP)*XN(2) + (RQ-RP*DCOS(THETA))*XN(1))/R
          !Sum up the quadrature scheme
          KB1=KB1+W1(ngp)*((POWER-1.0D0)*DEXP(POWER)+1.0D0)/(R*R)*DRDNU
        ENDDO
        !Multiply by the Jacobian
        KB1=KB1*2.0D0*PI

        !Evaluate KB2 using Elliptic integrals
        !Refer Brebbia et. al. pg 97
        KB2=4.0D0/DSQRT(A+B)*(((RP*RP-RQ*RQ+(ZP-ZQ)*(ZP-ZQ))/
     '    (A-B)*EDP(M)-KDP(M))*XN(1)/(2.0D0*RQ) + (ZP-ZQ)/
     '    (A-B)*EDP(M)*XN(2))

        !Evaluate the normal derivative by adding the quadrature and the
        !elliptic integral parts.
        DGREENA=(KB1+KB2)/(4.0D0*PI)
C        IF(DOP) THEN
C          WRITE(IO4,*)'DGREENA= ',DGREENA
C        ENDIF

      ENDIF
      RETURN
      END


