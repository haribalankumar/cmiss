      SUBROUTINE SECONDELLIPTIC(RESULT,first_elliptic,k,x,y,z,ERROR,*)

C#### Subroutine:  SECONDELLIPTIC
C###  Description:
C###  This file calculates the elliptic integral of the first kind.
C###  See-Also: MAGSOLCOIL

C**** File created by Huia Burt - November 2004
C**** Code copied from Dennison, 1998.

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      REAL*8 first_elliptic, k, x, y, z, RESULT
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER continue
      REAL*8 c1, c2, c3, c4, sigma, power4, mu 
      REAL*8 xn, yn, zn
      REAL*8 xndev, yndev, zndev, epsilon
      REAL*8 xnroot, ynroot, znroot, lambda 
      REAL*8 ea, eb, ec, ed, ef, s1, s2, interim_result
      REAL*8 error_tol

      CALL ENTERS('SECONDELLIPTIC',*9999)
      
      continue = 1
      error_tol = 0.003d0
      c1 = 3.0d0/14.0d0
      c2 = 1.0d0/6.0d0
      c3 = 9.0d0/22.0d0
      c4 = 3.0d0/26.0d0

C     Calculate the elliptic integral of the second kind
      xn = x
      yn = y
      zn = z

      sigma = 0.0d0
      power4 = 1.0d0
      
      DO WHILE (continue .EQ. 1)
         mu = ( xn + yn + (3.0d0*zn))*0.2d0       
         xndev = (mu-xn)/mu
         yndev = (mu-yn)/mu
         zndev = (mu-zn)/mu
         epsilon = max(abs(xndev),abs(yndev),abs(zndev))
         IF (epsilon .LT. error_tol) THEN
            continue = 0
         ELSE
            xnroot = sqrt(xn)
            ynroot = sqrt(yn)
            znroot = sqrt(zn)
            lambda = xnroot*(ynroot+znroot) +ynroot*znroot
            sigma = sigma+power4/(znroot*(zn+lambda))
            power4 = power4*0.25d0
            xn = (xn+lambda)*0.25d0
            yn = (yn+lambda)*0.25d0
            zn = (zn+lambda)*0.25d0
         ENDIF
      ENDDO
      
      ea = xndev*yndev
      eb = zndev*zndev
      ec = ea-eb
      ed = ea-6.0d0*eb
      ef = ed+ec+ec
      s1 = ed*(-c1 + (0.25d0*c3*ed) - (1.5d0*c4*zndev*ef))
      s2 = zndev*(c2*ef + zndev*(-c3*ec+zndev*c4*ea))

      interim_result = 3.0d0*sigma + power4*(1.0d0+s1+s2)/(mu*sqrt(mu))

      RESULT = first_elliptic - ((k*k)/3.0d0)*interim_result 

      CALL EXITS('SECONDELLIPTIC')
      RETURN
9999  CALL ERRORS('SECONDELLIPTIC',ERROR)
      CALL EXITS('SECONDELLIPTIC')
      RETURN 1
      END


