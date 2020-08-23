      REAL*8 FUNCTION ANALY_PHI_M(dipole_x,dipole_y,dipole_z,
     '  radius,theta,phi)

C#### Function: ANALY_PHI_M
C###  Type: REAL*8
C###  Description:
C###    <HTML>
C###      Analytic transmembrane potential for a dipole in free space.
C###      The dipole direction is specified and the location
C###      you are evaluating at in spherical polar coordinates.
C###    </HTML>

      REAL*8 dipole_x,dipole_y,dipole_z,radius,theta,phi

      ANALY_PHI_M=(2.D0*radius**3+1)/radius**2
     '  * (
     '  dipole_x*DCOS(theta)*DSIN(phi)
     '  + dipole_y*DSIN(theta)*DSIN(phi)
     '  + dipole_z*DCOS(phi)
     '  )


      RETURN
      END
