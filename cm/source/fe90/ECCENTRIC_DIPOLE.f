      SUBROUTINE ECCENTRIC_DIPOLE(SIGMA,DIPOLE_ORIENTATION,
     '  DIPOLE_POSITION,POSITION,POTENTIAL,GRAD_POTENTIAL,ERROR,*)

C#### Subroutine: ECCENTRIC_DIPOLE
C###  Description:
C###    ECCENTRIC_DIPOLE evaluates the potential on the surface of
C###    a sphere containing an eccentric dipole.  The result is
C###    returned in POTENTIAL.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
!     Parameter list
      REAL*8 SIGMA,DIPOLE_ORIENTATION(3),
     '  DIPOLE_POSITION(3),POSITION(3),POTENTIAL,
     '  GRAD_POTENTIAL(3)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nj,nj_deriv
      REAL*8 LEAD_FIELD_GRAD
      REAL*8 LEAD_FIELD,p0,r,s,term1

      CALL ENTERS('ECCENTRIC_DIPOLE',*9999)

C GMH 18/4/96 The following formula is taken from Brody, Terry and
C     Ideker 'Eccentric Dipole in a Spherical Medium:  Generalised
C     Expression for Surface Potentials'  IEEE Transactions on
C     Biomedical Engineering, March 1973 p141-143.
C       Initialise constants
      p0=0.0d0
      s=0.0d0
      r=0.0d0
      DO nj=1,3
        p0=p0+(POSITION(nj)-DIPOLE_POSITION(nj))**2
        s=s+(POSITION(nj)*DIPOLE_POSITION(nj))
        r=r+POSITION(nj)**2
      ENDDO !nj
      p0=DSQRT(p0)
      r=DSQRT(r)
      term1=1.0d0/(4.0d0*PI*SIGMA*p0)
C     Sum the lead field components
      POTENTIAL=0.0d0
      DO nj_deriv=1,3
        GRAD_POTENTIAL(nj_deriv)=0.0d0
      ENDDO !nj_deriv
      DO nj=1,3
        DO nj_deriv=1,3
C Add the deriv of the njth field compo wrt nj_deriv to grad_pot(nj_d)
C         dF/dS
          LEAD_FIELD_GRAD=(term1/(r*r)*
     '      (((p0+r-s/r)*POSITION(nj)/r+
     '      (POSITION(nj)*s/r-r*DIPOLE_POSITION(nj))/r)/
     '      ((p0+r-s/r)*(p0+r-s/r))))*DIPOLE_POSITION(nj_deriv)
C         dF/dp0
          LEAD_FIELD_GRAD=LEAD_FIELD_GRAD-
     '      (term1*((6.0d0*(POSITION(nj)-DIPOLE_POSITION(nj))/
     '      (p0*p0*p0))+
     '      (1.0d0/(r*r)*(POSITION(nj)/p0+
     '      ((POSITION(nj)*s/r-r*DIPOLE_POSITION(nj))*
     '      (2.0d0*p0+r-s/r))/
     '      (p0*(p0+r-s/r)*(p0+r-s/r))))))*
     '      ((POSITION(nj_deriv)-DIPOLE_POSITION(nj_deriv))/p0)
          IF(nj.EQ.nj_deriv) THEN
C           dF/dx
            LEAD_FIELD_GRAD=LEAD_FIELD_GRAD+
     '        term1*(2.0d0/(p0*p0)+1.0d0/(r*r)*
     '        (1.0d0+(s/r)/(p0+r-s/r)))
          ENDIF
          GRAD_POTENTIAL(nj_deriv)=GRAD_POTENTIAL(nj_deriv)+
     '      LEAD_FIELD_GRAD*DIPOLE_ORIENTATION(nj)
        ENDDO !nj_deriv
        LEAD_FIELD=term1*(2.0d0*(POSITION(nj)-DIPOLE_POSITION(nj))/
     '    (p0*p0)+1.0d0/(r*r)*(POSITION(nj)+
     '    (POSITION(nj)*s/r-r*DIPOLE_POSITION(nj))/
     '    (p0+r-s/r)))
        POTENTIAL=POTENTIAL+LEAD_FIELD*DIPOLE_ORIENTATION(nj)
      ENDDO !nj

      CALL EXITS('ECCENTRIC_DIPOLE')
      RETURN
 9999 CALL ERRORS('ECCENTRIC_DIPOLE',ERROR)
      CALL EXITS('ECCENTRIC_DIPOLE')
      RETURN 1
      END


