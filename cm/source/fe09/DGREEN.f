      REAL*8 FUNCTION DGREEN(IGREN,IMAT,k,l,CE,DRDN,RAD,
     '  RJAC,XN,XR)

C#### Function: DGREEN
C###  Type: REAL*8
C###  Description:
C###    DGREEN calculates the normal derivative of the appropriate
C###    Green's function. RJAC is the Jacobian obtained from any element
C###    subdivision or polar coordinate transformation. If the
C###    generalised Laplace equation is being solved then IMAT
C###    determines whether the conductivity needs to be included in the
C###    DGREEN expression (depends on the coefficient of DGREEN in the
C###    calling routine).  IMAT=0 means no material properties to be
C###    included and IMAT=1 means to include them. This allows different
C###    regions to be coupled correctly. See LE BEM comment in GREEN.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IGREN,IMAT,k,l
      REAL*8 CE(NMM),DRDN,RAD,RJAC,XN(NJM),XR(NJM)
!     Local Variables
      CHARACTER ERROR*10
      REAL*8 C2,C3,C4,K1
      GOTO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150)IGREN
10    CONTINUE !Normal deriv. of Green's function Laplace eqn,2D cart
        DGREEN=-1.0d0/(2.0d0*PI*RAD)*DRDN*RJAC
        GOTO 9999
20    CONTINUE !Normal deriv. of Green's function Laplace eqn,3D cart
        DGREEN=-1.0d0/(4.0d0*PI*RAD*RAD)*DRDN*RJAC
        GOTO 9999
30    CONTINUE !Normal deriv. of Green's function Helmholtz eqn,2D cart
        GOTO 9999
40    CONTINUE !Normal deriv. of Green's function Helmholtz eqn,3D cart
        GOTO 9999
50    CONTINUE !Normal deriv. of Green's function Yukawa eqn,2D cart
        DGREEN=-CE(1)/(2.0d0*PI)*K1(CE(1)*RAD)*DRDN*RJAC
        GOTO 9999
60    CONTINUE !Normal deriv. of Green's function Yukawa eqn,3D cart
        DGREEN=-DEXP(-CE(1)*RAD)*(CE(1)+1.0d0/RAD)/
     '    (4.0d0*PI*RAD)*DRDN*RJAC
        GOTO 9999
70    CONTINUE !Normal deriv. generalised Laplace eqtn, 2d cart
        !Depending on the coefficient of DGREEN the
        !conductivity is included or not
        IF(IMAT.EQ.0) THEN
          DGREEN=-1.0d0/(2.0d0*PI*RAD)*DRDN*RJAC
        ELSE
          DGREEN=-1.0d0/(2.0d0*PI*RAD*CE(1))*DRDN*RJAC
        ENDIF
        GOTO 9999
80    CONTINUE !Normal deriv. generalised Laplace eqtn, 3d cart
        !Depending on the coefficient of DGREEN the
        !conductivity is included or not (if the coefficient is just a
        !potential then it is not needed, but if it is the arc length
        !derivative then it does need to be included).
        IF(IMAT.EQ.0) THEN
          DGREEN=-1.0d0/(4.0d0*PI*RAD*RAD)*DRDN*RJAC
        ELSE
          DGREEN=-1.0d0/(4.0d0*PI*RAD*RAD*CE(1))*DRDN*RJAC
        ENDIF
        GOTO 9999
90    CONTINUE
        !Normal deriv. Poisson eqtn, 2d cartesian, special source term
        !Source term is -div(k.grad(g)), so use Laplace fundamental soln
        DGREEN=-1.0d0/(2.0d0*PI*RAD)*DRDN*RJAC
        GOTO 9999
100   CONTINUE
        !Normal deriv. Poisson eqtn, 2d cartesian, special source term
        !Source term is -div(k.grad(g)), so use Laplace fundamental soln
        DGREEN=-1.0d0/(4.0d0*PI*RAD*RAD)*DRDN*RJAC
        GOTO 9999
110   CONTINUE !Normal deriv. Helmholtz eqtn, COMPLEX valued
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>Error in DGREEN. IGREN=11. Should use '
     '    //'DGREENC'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        GOTO 9999
CC$      call mp_unsetlock()
120    CONTINUE !Modified Helmholtz(Yukawa) eq.,3D cyl. sym. z axis
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>Error in GREEN. IGREN=12. Should use '
     '   //'GREENC'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
        GOTO 9999
130    CONTINUE !Green's function for Linear Elasticity - 2D plane strain
C       CE(1) = Youngs modulus
C       CE(2) = Poissons ratio
        C2=-1.0d0/(4.0d0*PI*(1.0D0-CE(2)))
        C3=1.0d0-2.0d0*CE(2)
        C4=2.0d0
        IF(k.EQ.l) THEN
          DGREEN=C2/RAD*(C3+C4*(XR(k)/RAD)*(XR(k)/RAD))*DRDN
        ELSE
          DGREEN=C2/RAD*(C4*(XR(k)/RAD)*(XR(l)/RAD)*DRDN-
     '      C3*((XR(k)/RAD)*XN(l)-(XR(l)/RAD)*XN(k)))
        ENDIF
        GOTO 9999
140    CONTINUE !Green's function for Linear Elasticity - 2D plane stress
        C2=-(1.0d0+CE(2))/(4.0d0*PI)
        C3=(1.0d0-CE(2))/(1.0d0+CE(2))
        C4=2.0d0
        IF(k.EQ.l) THEN
          DGREEN=C2/RAD*(C3+C4*(XR(k)/RAD)*(XR(k)/RAD))*DRDN
        ELSE
          DGREEN=C2/RAD*(C4*(XR(k)/RAD)*(XR(l)/RAD)*DRDN-
     '      C3*((XR(k)/RAD)*XN(l)-(XR(l)/RAD)*XN(k)))
        ENDIF
        GOTO 9999
150    CONTINUE !Green's function for Linear Elasticity - 3D
        C2=-1.0d0/(8.0d0*PI*(1.0d0-CE(2)))
        C3=1.0d0-2.0d0*CE(2)
        C4=3.0d0
        IF(k.EQ.l) THEN
          DGREEN=C2/(RAD*RAD)*(C3+C4*(XR(k)/RAD)*(XR(k)/RAD))*DRDN
        ELSE
          DGREEN=C2/(RAD*RAD)*(C4*(XR(k)/RAD)*(XR(l)/RAD)*DRDN-
     '      C3*((XR(k)/RAD)*XN(l)-(XR(l)/RAD)*XN(k)))
        ENDIF
        GOTO 9999
 9999 RETURN
      END


