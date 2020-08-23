      REAL*8 FUNCTION GREENA(IGREN,CE,RP,ZP,RQ,ZQ)

C#### Function: GREENA
C###  Type: REAL*8
C###  Description:
C###    GREENA identifies appropriate complex Green's function for
C###    axisymmetric problems.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'eqt000.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IGREN
      REAL*8 CE(NMM),ZP,RP,ZQ,RQ
!     Local Variables
      INTEGER ngp
      REAL*8 A,B,KA1,KA2,KDP,M,POWER,R,THETA

      IF(IGREN.EQ.12) THEN ! Yukawa equation,3d cylindrically symmetric
        !Reference: Brebbia, Telles and Wrobel pg 97
        A = RP**2 + RQ**2 + (ZP-ZQ)**2
        B = 2.0D0 * RP * RQ
        M = 2.0D0 * B/(A+B)

        !Quadrature on first part to find KA1
        KA1=0.0D0
        DO ngp=1,ng
          THETA=2.0D0*PI*D1(ngp)
          R=DSQRT(A-B*DCOS(THETA))
          POWER=-1.0D0*CE(1)*R
          KA1=KA1+W1(ngp)*(DEXP(POWER)-1.0D0)/R
        ENDDO
        KA1=KA1*2.0D0*PI

        !Evaluate KA2 using Elliptic integrals
        KA2=4.0D0*KDP(M)/DSQRT(A+B)

        !Evaluate the value of the Greens function by summing the
        !quadrature and the elliptic integral parts.
        GREENA=(KA1+KA2)/(4.0D0*PI)
C        IF(DOP) THEN
C          WRITE(IO4,*)
C          WRITE(IO4,*)'GREENA= ',GREENA
C        ENDIF
      ENDIF
      RETURN
      END


