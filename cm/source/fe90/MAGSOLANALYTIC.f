      SUBROUTINE MAGSOLANALYTIC(NPNODE,DIPOLE_CEN,DIPOLE_DIR,H_FIELD,XP,
     '  INT_POINT,MAGDIPOLE,MAGVOLUME,VPOTENTIAL,ERROR,*)

C#### Subroutine: MAGSOLANALYTIC
C###  Description:
C###    <HTML>      
C###    <P>MAGSOLANALYTIC calculates the analytic magnetic field
C###    that results from a dipole source in the x direction
C###    centered on the x axis.</P>
C###    <P>NOTE: That this analytic solution uses a different version
C###    of spherical polar coordinates (theta and phi swapped) than
C###    is used in CMISS. Also uses theta(phi) from Z axis. The
C###    analytic solution comes from Cuffin and Cohen, 1977,
C###    <I>Magnetic fields of a dipole in special volume conductor
C###    shapes</I> IEEETBME 24(4)</P>
C###    <P>Note that the solution in the paper is wrong! but has been
C###    corrected before being implemented here.</P>
C###    </HTML>      

C**** Created by Martin Buist - December 2001

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER NPNODE(0:NP_R_M)
      REAL*8 DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM),H_FIELD(*),INT_POINT(NJM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL MAGDIPOLE,MAGVOLUME,VPOTENTIAL
!     Local Variables
      INTEGER nj,np
      REAL*8 a,costheta,cosphi,DATAN_MOD,dAdPhi,dAdTheta,GAMMA,
     '  JACOBIAN(9),H_phi,H_r,H_theta,Hd_phi,Hd_r,Hd_theta,Hv_phi,Hv_r,
     '  Hv_theta,phi,r,rho,rS,sinphi,sintheta,x,y,z
      LOGICAL PROBLEM

      CALL ENTERS('MAGSOLANALYTIC',*9999)

      PROBLEM=.FALSE.
      CALL ASSERT(NJT.EQ.3,'>>ERROR - analytic solution requires 3D',
     '  ERROR,*9999)

C***  The dipole centre is (0,0,a) in rc and
C***  the dipole vector is (rho,0,0) in rc.
C***  Any other definitions will break the analytic solution
      a=DIPOLE_CEN(3,0,1)
      rho=DIPOLE_DIR(1,0,1)

C***  The rS term is the radius of the sphere.
      np=NPNODE(1)
      rS=0.0d0
      DO nj=1,NJT
        rS=rS+(XP(1,1,nj,np)*XP(1,1,nj,np))
      ENDDO
      rS=DSQRT(rS)

C***  Convert the point of interest into a spherical polar
C***  coordinate, also set up SIN and COS of the two angles.
C***  Note this is done in terms of the true origin and not
C***  from 'a'. The spherical polar system here is different
C***  from the one in the rest of cmiss. Here phi runs in the XY
C***  plane, anticlockwise from the positive X axis and theta runs
C***  from the positive Z axis.
      x=INT_POINT(1)
      y=INT_POINT(2)
      z=INT_POINT(3)

C      x=1.76776695296637d0
C      y=3.06186217847897d0
C      z=3.53553390593274d0

      r=DSQRT((x*x)+(y*y)+(z*z))
      phi=DATAN_MOD(x,y)
      cosphi=DCOS(phi)
      sinphi=DSIN(phi)
      IF(r.GT.ZERO_TOL) THEN
        costheta=z/r
        sintheta=DSQRT((x*x)+(y*y))/r
C        theta=DACOS(z/r)
      ELSE
        PROBLEM=.TRUE.
      ENDIF
C***  Note that the reverse of these is:
C***  x = r cosphi sintheta
C***  y = r sinphi sintheta
C***  z = r costheta

C***  Check all denominator multipliers for zeroes
      IF(DABS(r).LT.ZERO_TOL) PROBLEM=.TRUE.
      IF(DABS(a).LT.ZERO_TOL) PROBLEM=.TRUE.
      IF(DABS(sintheta).LT.ZERO_TOL) PROBLEM=.TRUE.

C***  r*sqrt(GAMMA) is the radius from the field point to the
C     dipole located at z=a. If this radius is called rdash then
C     GAMMA is equal to (rdash**2)/(r**2). This result is easily
C     obtained from the cosine law.
      IF(.NOT.PROBLEM) THEN
        GAMMA=1.0d0-(2.0d0*a*costheta/r)+((a/r)**2.0d0)
        IF(DABS(GAMMA).LT.ZERO_TOL) PROBLEM=.TRUE.
      ENDIF

      IF(.NOT.PROBLEM) THEN
C***    Analytic solution for Hv - field due to the spherical
C       volume conductor only (polar).
        IF(MAGVOLUME) THEN
C***       These are the derivatives dA/dtheta and dA/dphi
          dAdTheta=(-rho*cosphi*rS/(4.0d0*pi*r*r*(GAMMA**1.5d0)))*
     '      (costheta-(a/r)-((GAMMA/(sintheta*sintheta))*
     '      (costheta-(r/a)+(r*DSQRT(GAMMA)/a))))
          dAdPhi=(rho*sinphi*rS/(4.0d0*pi*r*r*DSQRT(GAMMA)))*
     '      (((costheta/sintheta)*(costheta-(r/a)+
     '      (r*DSQRT(GAMMA)/a)))+sintheta)

C***      Multiply through by the Jacobian (d(r,t,p)/d(x,y,z))
          Hv_theta=(1.0d0/(r**2.0d0*sintheta))*dAdPhi
          Hv_phi=(-1.0d0/(r**2.0d0*sintheta))*dAdTheta

          Hv_r=0.0d0
          IF(r.LT.rS) THEN
            Hv_theta=0.0d0
            Hv_phi=0.0d0
          ENDIF
        ELSE
          Hv_r=0.0d0
          Hv_theta=0.0d0
          Hv_phi=0.0d0
        ENDIF

C***    Analytic solution for Hd - field due to a current
C***    dipole only (polar).  (P x R)/4 pi r^3
        IF(MAGDIPOLE) THEN
C***      Note that in the paper this term is missing a (1/r)
          Hd_r=(a*sintheta*sinphi*rho)/(4.0d0*PI*r*r*r*GAMMA**1.5d0)

C***      Note that in the paper this term is missing a (1/r)
          Hd_theta=((rho*sinphi)/(4.0d0*PI*r*r*r*GAMMA**1.5d0))*
     '      ((a*costheta/r)-1.0d0)

C***      Note that in the paper this term is missing a (1/(r sintheta))
          Hd_phi=((rho*cosphi)/(4.0d0*PI*r*r*r*GAMMA**1.5d0*sintheta))*
     '      ((a/r)-costheta)
        ELSE
          Hd_r=0.0d0
          Hd_theta=0.0d0
          Hd_phi=0.0d0
        ENDIF

C***    H=Hd-Hv
        H_r=Hd_r-Hv_r
        H_theta=Hd_theta-Hv_theta
        H_phi=Hd_phi-Hv_phi

C***    Analytic solution for the vector potential field 'A'
        IF(VPOTENTIAL) THEN
C***      Result from the infinite sum
          H_r=((costheta/sintheta)*(costheta-(r/a)+(DSQRT(GAMMA)*r/a)))
     '      +sintheta

C***      Multiplied by the terms outside the sum
          H_r=H_r*(-rho*cosphi*rS/(4.0d0*pi*r*r*DSQRT(GAMMA)))

C***      It is important to note that the solution is not valid
C***      inside the spherical domain, set to zero here.
          IF(r.LT.rS) H_r=0.0d0

C***      Only non-zero component is normal to the sphere
          H_theta=0.0d0
          H_phi=0.0d0
        ENDIF

C***    This is the Jacobian from the spherical->rc transform
        JACOBIAN(1)=cosphi*sintheta     !=dx/dr
        JACOBIAN(2)=r*cosphi*costheta   !=dx/dtheta
        JACOBIAN(3)= -r*sinphi*sintheta !=dx/dphi
        JACOBIAN(4)=sinphi*sintheta     !=dy/dr
        JACOBIAN(5)=r*sinphi*costheta   !=dy/dtheta
        JACOBIAN(6)=r*cosphi*sintheta   !=dy/dphi
        JACOBIAN(7)=costheta            !=dz/dr
        JACOBIAN(8)= -r*sintheta        !=dz/dtheta
        JACOBIAN(9)= 0.0d0              !=dz/dphi

C***    dx=(dx/dr)*dr + (dx/dtheta)*dtheta+ (dx/dphi)*dphi
        H_FIELD(1)=(JACOBIAN(1)*H_r)+(JACOBIAN(2)*H_theta)+
     '    (JACOBIAN(3)*H_phi)
C***    dy=(dy/dr)*dr + (dy/dtheta)*dtheta+ (dy/dphi)*dphi
        H_FIELD(2)=(JACOBIAN(4)*H_r)+(JACOBIAN(5)*H_theta)+
     '    (JACOBIAN(6)*H_phi)
C***    dz=(dz/dr)*dr + (dz/dtheta)*dtheta+ (dz/dphi)*dphi
        H_FIELD(3)=(JACOBIAN(7)*H_r)+(JACOBIAN(8)*H_theta)+
     '    (JACOBIAN(9)*H_phi)

C***    LEAVE - Magnetic field lines from a current
C       source in free space (rc).
C        VD(1)=1.0d0
C        VD(2)=0.0d0
C        VD(3)=0.0d0
C        CD(1)=0.0d0
C        CD(2)=0.0d0
C        CD(3)=1.0d0
C        !TERM1-3 are the 3 vector radius components
C        !TERM4 is the scalar radius
C        !TERM5 is the constant multiplier
C        TERM1=INT_POINT(1)-CD(1)
C        TERM2=INT_POINT(2)-CD(2)
C        TERM3=INT_POINT(3)-CD(3)
C        TERM4=DSQRT(TERM1**2.0d0+TERM2**2.0d0+TERM3**2.0d0)
C        IF(DABS(TERM4).GT.ZERO_TOL) THEN
C          TERM5=1.0d0/(4.0d0*PI*TERM4*TERM4*TERM4)
C          H_FIELD(1)=(VD(2)*TERM3-VD(3)*TERM2)*TERM5
C          H_FIELD(2)=(VD(3)*TERM1-VD(1)*TERM3)*TERM5
C          H_FIELD(3)=(VD(1)*TERM2-VD(2)*TERM1)*TERM5
C        ELSE
C          H_FIELD(1)=0.0d0
C          H_FIELD(2)=0.0d0
C          H_FIELD(3)=0.0d0
C        ENDIF
C       End magnetic field in free space

C***    LEAVE - Potential field from a current source
C       in free space (rc) (into H_FIELD(1)).
C        VD(1)=1.0d0
C        VD(2)=0.0d0
C        VD(3)=0.0d0
C        CD(1)=0.0d0
C        CD(2)=0.0d0
C        CD(3)=1.0d0
C        !TERM1-3 are the 3 vector radius components
C        !TERM4 is the scalar radius
C        !TERM5 is the constant multiplier
C        TERM1=INT_POINT(1)-CD(1)
C        TERM2=INT_POINT(2)-CD(2)
C        TERM3=INT_POINT(3)-CD(3)
C        TERM4=DSQRT(TERM1**2.0d0+TERM2**2.0d0+TERM3**2.0d0)
C        IF(DABS(TERM4).GT.ZERO_TOL) THEN
C          !Assuming unit conductivity here
C          TERM5=1.0d0/(4.0d0*PI*TERM4*TERM4*TERM4)
C          H_FIELD(1)=(VD(1)*TERM1+VD(2)*TERM2+VD(3)*TERM3)*TERM5
C          H_FIELD(2)=0.0d0
C          H_FIELD(3)=0.0d0
C        ELSE
C          H_FIELD(1)=0.0d0
C          H_FIELD(2)=0.0d0
C          H_FIELD(3)=0.0d0
C        ENDIF
C***    End potential in free space

C***    LEAVE - This is the solution in RC for Hd
C        rdash=r*DSQRT(GAMMA) !adjust to radius from 'a'
C        H_FIELD(1)=0.0d0
C        H_FIELD(2)=-rho*(z-a)/(4.0d0*pi*rdash*rdash*rdash)
C        H_FIELD(3)=rho*y/(4.0d0*pi*rdash*rdash*rdash)

C***    LEAVE - This is the solution in RC for Hd with no x,y,z
C        rdash=r*DSQRT(GAMMA) !adjust to radius from 'a'
C        H_FIELD(1)=0.0d0
C        H_FIELD(2)=-rho*(r*costheta-a)/(4.0d0*pi*rdash*rdash*rdash)
C        H_FIELD(3)=rho*r*sinphi*sintheta/(4.0d0*pi*rdash*rdash*rdash)

C        write(*,*) H_FIELD(1),H_FIELD(2),H_FIELD(3)

      ELSE !is a problem
        WRITE(OP_STRING,'('' Warning: zero denominator detected'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        H_FIELD(1)=0.0d0
        H_FIELD(2)=0.0d0
        H_FIELD(3)=0.0d0
      ENDIF !problem

      CALL EXITS('MAGSOLANALYTIC')
      RETURN
9999  CALL ERRORS('MAGSOLANALYTIC',ERROR)
      CALL EXITS('MAGSOLANALYTIC')
      RETURN 1
      END


