      SUBROUTINE MAGSOLCOIL(B,SQUID_1,INT_POINT,
     '  COIL_RADIUS,CURR,COIL_HEIGHT,ERROR,*)

C#### Subroutine: MAGSOLCOIL
C###  Description:
C###    This file calculates the magneticfield at the squid sensor 
C###    positions due to a circular current filament (1 turn coil).
C###  See-Also: EVSOLU, FIRSTELLIPTIC, SECONDELLIPTIC

C**** File created by Huia Burt, November 2004
C**** Copied from code by Leo Cheng (coil.c)
 
      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      REAL*8 SQUID_1(NJM),INT_POINT(NJM)
      REAL*8 COIL_RADIUS, CURR, COIL_HEIGHT
      REAL*8 B(NJM)
      CHARACTER ERROR*(*)

!     Local Variables
      REAL*8 MU0, x_squid, z_squid, a, current
      REAL*8 y,x_cen, z_cen
      REAL*8 sq_x, sq_z, r_final
      REAL*8 al, alt4, be, ga, al2, be2, q, rq, k
      REAL*8 x_input, y_input, z_input
      REAL*8 fk, ek, Ht, Hy, Hr, theta
      REAL*8 Br, Bz, Bx, By

      CALL ENTERS('MAGSOLCOIL',*9999)

C     Note: y and z are switched around from the coil.c code as
C     our squid is positioned facing outwards into -y

C     Define Mu
      MU0 = PI*0.0000004d0

C     Define coil radius (mm)
      a=COIL_RADIUS

C     Define current
      current=CURR
      
C     Define the coordinates of the centre of the dewar, which is
C     where we are locating our dipole
      x_cen=SQUID_1(1)
      z_cen=SQUID_1(3)

C     Define the squid point coordinates
      x_squid=INT_POINT(1)
      z_squid=INT_POINT(3)
 
C     Calculate the distance from the dipole centre to the point
C     of interest (the squid sensor)
      sq_x = (x_squid-x_cen)*(x_squid-x_cen)
      sq_z = (z_squid-z_cen)*(z_squid-z_cen)
      r_final = DSQRT(sq_x + sq_z)

C     Define the height of the sensor above the coil in mm
      y=COIL_HEIGHT

C     Calculate k
      al = r_final/a
      alt4 = al*4.0d0
      be = y/a
      IF(r_final.GT.1.0D-30) THEN
          ga = y/r_final
      ENDIF
      al2 = al*al
      be2 = be*be
      q = (1+al)*(1+al) + be2
      rq = sqrt(q)
      k = sqrt(alt4/q)
 
C     Calculate the magnetic field intensity in the y_dir
      Ht = current/(2.0d0*a*PI*rq)

C     Define inputs for calculating the 1st & 2nd elliptic integrals
      x_input = 0.0d0         
      y_input = 1.0d0 - (k*k) 
      z_input = 1.0d0        

C     Calculate the 1st and 2nd elliptic integrals
      CALL FIRSTELLIPTIC(fk,x_input,y_input,z_input,ERROR,*9999)

      CALL SECONDELLIPTIC(ek,fk,k,x_input,y_input,z_input,ERROR,*9999)
      
C     Calculate the radial magnetic field and the magnetic field in the
C     y direction
      IF (r_final.LE.1.0D-30) THEN
         Hr = 0.0d0
      ELSE
         Hr = (Ht*ga*(ek*(1.0d0 + al2 + be2)/(q-alt4)-fk))
      ENDIF

      Hy = Ht*(ek*(1-al2-be2)/(q-alt4)+fk)
      
C     Calculate the magnetic field
C     Note we multiply by 1000 as H has units T/mm.
      Br = Hr*MU0*1000
      By = Hy*MU0*1000
       
      if (ABS(z_squid-z_cen) .LE. 1.0D-30) then
        theta=0.0d0       ! 0 degrees
      else
        theta=atan(ABS(z_squid-z_cen)/ABS(x_squid-x_cen))
      endif

C     Calculate the magnetic field in the x and z directions
      Bx=Br*cos(theta)
      Bz=Br*sin(theta)

C     Change the direction of the magnetic field depending on where
C     the squid point sits wrt the coil
      IF(x_squid .GT. x_cen) THEN
        Bx = -Bx
      ENDIF
      IF(z_squid .GT. z_cen) THEN
        Bz = -Bz
      ENDIF

C     Fix up the points at x=0
      IF (ABS(x_squid-x_cen) .LE. 1.0D-30) THEN

            Bz = Bx
            Bx = 0.0d0
   
      ENDIF
      
C     Store the magnetic flux results in a vector
      B(1) = Bx
      B(2) = By
      B(3) = Bz
      
      CALL EXITS('MAGSOLCOIL')
      RETURN
9999  CALL ERRORS('MAGSOLCOIL',ERROR)
      CALL EXITS('MAGSOLCOIL')
      RETURN 1
      END


