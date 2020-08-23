      SUBROUTINE FIRSTELLIPTIC(RESULT,x,y,z,ERROR,*)

C#### Subroutine:  FIRSTELLIPTIC
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
      REAL*8 x,y,z, RESULT
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER continue
      REAL*8 c1, c2, c3 
      REAL*8 xn, yn, zn, mu 
      REAL*8 xndev, yndev, zndev, epsilon
      REAL*8 xnroot, ynroot, znroot, lambda 
      REAL*8 state_2, state_3, s
      REAL*8 error_tol
      
      CALL ENTERS('FIRSTELLIPTIC',*9999)
      
      continue = 1
      error_tol = 0.003d0
      c1 = 1.0d0/24.0d0 
      c2 = 3.0d0/44.0d0 
      c3 = 1.0d0/14.0d0

C     Calculate the elliptic integral of the first kind

      xn = x
      yn = y
      zn = z

      DO WHILE (continue .EQ. 1)
         mu = (xn+yn+zn)/3.0d0
         xndev = 2.0d0-(mu+xn)/mu
         yndev = 2.0d0-(mu+yn)/mu
         zndev = 2.0d0-(mu+zn)/mu
         epsilon = max(abs(xndev),abs(yndev),abs(zndev))
         IF (epsilon .LT. error_tol) THEN
            continue = 0
         ELSE
            xnroot = sqrt(xn)
            ynroot = sqrt(yn)
            znroot = sqrt(zn)
            lambda = xnroot*(ynroot+znroot) + (ynroot*znroot)
            xn = (xn+lambda)*0.25d0
            yn = (yn+lambda)*0.25d0
            zn = (zn+lambda)*0.25d0
         ENDIF
      ENDDO
      
      state_2 = (xndev*yndev) - (zndev*zndev)
      state_3 = xndev*yndev*zndev
      s = 1.0d0 + (c1*state_2 - 0.1d0 - c2*state_3)*state_2 + c3*state_3
      RESULT = s/sqrt(mu)

      CALL EXITS('FIRSTELLIPTIC')
      RETURN
9999  CALL ERRORS('FIRSTELLIPTIC',ERROR)
      CALL EXITS('FIRSTELLIPTIC')
      RETURN 1
      END


