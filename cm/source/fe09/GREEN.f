      REAL*8 FUNCTION GREEN(IGREN,FLAG,k,l,CE,RAD,RJAC,XR)

C#### Function: GREEN
C###  Type: REAL*8
C###  Description:
C###    GREEN identifies appropriate Real Green's function. RJAC is the
C###    Jacobian obtained from any element subdivision or polar
C###    coordinate transformation.  FLAG is a flag to indicate anything
C###    special should be done to the Green's function. RJAC=1 if no
C###    transformation has been used. Linear elasticity formulae taken
C###    from BEER&WATSON, "Intro to finite and boundary element
C###    methods for engineers", p487ish.

C**** k and l are used for Linear Elasticity.  k is the direction of
C**** the point load, l is the direction of displacement.
C**** XR(i) = X(i)-X(i)* ie projection of radius in direction i.
C**** Linear Elasticity:
C**** There are discrepencies wrt the formulae for t*ij.
C**** Essentially, the sign of the dr/dn term wrt the r,i*n,j terms
C**** alternates.  We use the expression found in AJP's notes
C**** - supported by Beskos(1987) (1)
C****                Brebbia and Orszag(1983)
C****                Bannerjee and Butterfield(1979) (1)
C**** - contradicted by Beer and Watson(1992)
C****                   Brebbia and Connor(1987)
C****                   Beskos(1987) (2)
C****                   Bannerjee and Butterfield(1979) (2)
C**** There are no contradictions for u*ij.
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IGREN,FLAG,k,l
      REAL*8 CE(NMM),RAD,RJAC,XR(NJM)
!     Local Variables
      CHARACTER ERROR*10
      REAL*8 C,C1,K0

      GOTO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150)IGREN
10    CONTINUE !Green's function for Laplace equation,2D cartesians
        IF(FLAG.EQ.1) THEN !logarithmic quadrature (non-log part)
          GREEN=-1.0d0/(2.0d0*PI)*DLOG(RAD)*RJAC
        ELSE IF(FLAG.EQ.2) THEN !logarithmic quadrature (log part)
          GREEN=1.0d0/(2.0d0*PI)*RJAC
        ELSE
          GREEN=-1.0d0/(2.0d0*PI)*DLOG(RAD)*RJAC
        ENDIF
        GOTO 9999
20    CONTINUE !Green's function for Laplace equation,3D cartesians
        GREEN=1.0d0/(4.0d0*PI*RAD)*RJAC
        GOTO 9999
30    CONTINUE !Green's function for Helmholtz equation,2D cartesians
        GOTO 9999
40    CONTINUE !Green's function for Helmholtz equation,3D cartesians
        GOTO 9999
50    CONTINUE !Green's function for Yukawa equation,2D cartesians
        GREEN=1.0d0/(2.0d0*PI)*K0(CE(1)*RAD)*RJAC
        GOTO 9999
60    CONTINUE !Green's function for Yukawa equation,3D cartesians
        GREEN=DEXP(-CE(1)*RAD)/(4.0d0*PI*RAD)*RJAC
        GOTO 9999
70    CONTINUE !Green's function for generalised Laplace eqtn, 2d cart
        !Include conductivity in fundamental solution so coupling is
        !handled correctly
        IF(FLAG.EQ.1) THEN !logarithmic quadrature (non-log part)
          GREEN=-1.0d0/(2.0d0*PI*CE(1))*DLOG(RAD)*RJAC
        ELSE IF(FLAG.EQ.2) THEN !logarithmic quadrature (log part)
          GREEN=1.0d0/(2.0d0*PI*CE(1))*RJAC
        ELSE
          GREEN=-1.0d0/(2.0d0*PI*CE(1))*DLOG(RAD)*RJAC
        ENDIF
        GOTO 9999
80    CONTINUE !Green's function for generalised Laplace eqtn, 3d cart
        !Include conductivity in fundamental solution so coupling is
        !handled correctly
        GREEN=1.0d0/(4.0d0*PI*RAD*CE(1))*RJAC
        GOTO 9999
90    CONTINUE
C       !Green's fn for Poisson eqtn, 2d cartesian, special source term
        !Source term is -div(k.grad(g)), so use Laplace fundamental soln
        IF(FLAG.EQ.1) THEN !logarithmic quadrature (non-log part)
          GREEN=-1.0d0/(2.0d0*PI)*DLOG(RAD)*RJAC
        ELSE IF(FLAG.EQ.2) THEN !logarithmic quadrature (log part)
          GREEN=1.0d0/(2.0d0*PI)*RJAC
        ELSE
          GREEN=-1.0d0/(2.0d0*PI)*DLOG(RAD)*RJAC
        ENDIF
        GOTO 9999
100    CONTINUE
C       !Green's fn for Poisson eqtn, 3d cartesian, special source term
        !Source term is -div(k.grad(g)), so use Laplace fundamental soln
        GREEN=1.0d0/(4.0d0*PI*RAD)*RJAC
        GOTO 9999
110    CONTINUE !Green's function for Helmholtz eqtn, COMPLEX valued
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>Error in GREEN. IGREN=11. Should use '
     '   //'GREENC'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
        GOTO 9999
120    CONTINUE !Modified Helmholtz(Yukawa) eq.,3D cyl. sym. z axis
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>Error in GREEN. IGREN=12. Should use '
     '   //'GREENC'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
        GOTO 9999
130    CONTINUE !Green's function for Linear Elasticity - 2D plane strain
C       CE(1) = Youngs modulus
C       CE(2) = Poissons ratio
        C=-(2.0D0*(1.0D0+CE(2)))/(8.0D0*PI*CE(1)*(1.0D0-CE(2)))
        IF(k.EQ.l) THEN
          C1=3.0D0-4.0D0*CE(2)
          GREEN=C*(C1*LOG(RAD)-(XR(k)/RAD)*(XR(k)/RAD))*RJAC
        ELSE
          GREEN=-C*(XR(k)/RAD)*(XR(l)/RAD)*RJAC
        ENDIF
        GOTO 9999
140    CONTINUE !Green's function for Linear Elasticity - 2D plane stress
        C=-(2.0D0*(1.0D0+CE(2))*(1.0D0+CE(2)))/(8.0D0*PI*CE(1))
        IF(k.EQ.l) THEN
          C1=3.0D0-4.0D0*CE(2)/(1.0D0+CE(2))
          GREEN=C*(C1*LOG(RAD)-(XR(k)/RAD)*(XR(k)/RAD))*RJAC
        ELSE
          GREEN=-C*(XR(k)/RAD)*(XR(l)/RAD)*RJAC
        ENDIF
        GOTO 9999
150    CONTINUE !Green's function for Linear Elasticity - 3D
        C=-(2.0D0*(1.0D0+CE(2)))/(16.0D0*PI*CE(1)*(1.0D0-CE(2)))
        IF(k.EQ.l) THEN
          C1=3.0D0-4.0D0*CE(2)
          GREEN=C/RAD*(C1+(XR(k)/RAD)*(XR(k)/RAD))*RJAC
        ELSE
          GREEN=C/RAD*(XR(k)/RAD)*(XR(l)/RAD)*RJAC
        ENDIF
        GOTO 9999
9999  RETURN
      END


