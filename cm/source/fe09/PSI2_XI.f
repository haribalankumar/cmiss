      REAL*8 FUNCTION PSI2_XI(nb,nu,nn,XI)

C#### Function: PSI2_XI
C###  Type: REAL*8
C###  Description:
C###    PSI2_XI evaluates simplex basis functions with Xi coordinates.

C**** Created by Carey Stevens 12 Aug 1997

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nb,nn,nu
      REAL*8 XI(3)

      IF(NBSC(2,nb).EQ.1) THEN !linear
        IF(nn.EQ.1) THEN
          IF(nu.EQ.1) THEN
            PSI2_XI=1.0d0-XI(1)-XI(2)
          ELSE IF(nu.EQ.2) THEN
            PSI2_XI=-1.0d0
          ELSE IF(nu.EQ.4) THEN
            PSI2_XI=-1.0d0
          ELSE
            PSI2_XI=0.0d0
          ENDIF
        ELSE IF(nn.EQ.2) THEN
          IF(nu.EQ.1) THEN
            PSI2_XI=XI(1)
          ELSE IF(nu.EQ.2) THEN
            PSI2_XI=1.0d0
          ELSE
            PSI2_XI=0.0d0
          ENDIF
        ELSE
          IF(nu.EQ.1) THEN
            PSI2_XI=XI(2)
          ELSE IF(nu.EQ.4) THEN
            PSI2_XI=1.0d0
          ELSE
            PSI2_XI=0.0d0
          ENDIF
        ENDIF
      ELSE IF(NBSC(2,nb).EQ.2) THEN !quadratic
        IF(nn.EQ.1) THEN
          IF(nu.EQ.1) THEN
            PSI2_XI=(1.0d0-XI(1)-XI(2))*(1.0d0-2.0d0*XI(1)-2.0d0*XI(2))
          ELSE IF(nu.EQ.2) THEN
            PSI2_XI=-3.0d0+4.0d0*XI(1)+4.0d0*XI(2)
          ELSE IF(nu.EQ.3) THEN
            PSI2_XI=4.0d0
          ELSE IF(nu.EQ.4) THEN
            PSI2_XI=-3.0d0+4.0d0*XI(1)+4.0d0*XI(2)
          ELSE IF(nu.EQ.5) THEN
            PSI2_XI=4.0d0
          ELSE IF(nu.EQ.6) THEN
            PSI2_XI=4.0d0
          ELSE
            PSI2_XI=0.0d0
          ENDIF
        ELSE IF(nn.EQ.2) THEN
          IF(nu.EQ.1) THEN
            PSI2_XI=4.0d0*XI(1)*(1.0d0-XI(1)-XI(2))
          ELSE IF(nu.EQ.2) THEN
            PSI2_XI=4.0d0-8.0d0*XI(1)-4.0d0*XI(2)
          ELSE IF(nu.EQ.3) THEN
            PSI2_XI=-8.0d0
          ELSE IF(nu.EQ.4) THEN
            PSI2_XI=-4.0d0*XI(1)
          ELSE IF(nu.EQ.6) THEN
            PSI2_XI=-4.0d0
          ELSE
            PSI2_XI=0.0d0
          ENDIF
        ELSE IF(nn.EQ.3) THEN
          IF(nu.EQ.1) THEN
            PSI2_XI=XI(1)*(2.0d0*XI(1)-1.0d0)
          ELSE IF(nu.EQ.2) THEN
            PSI2_XI=4.0d0*XI(1)-1.0d0
          ELSE IF(nu.EQ.3) THEN
            PSI2_XI=4.0d0
          ELSE
            PSI2_XI=0.0d0
          ENDIF
       ELSE IF(nn.EQ.4) THEN
          IF(nu.EQ.1) THEN
            PSI2_XI=4.0d0*XI(2)*(1.0d0-XI(1)-XI(2))
          ELSE IF(nu.EQ.2) THEN
            PSI2_XI=-4.0d0*XI(2)
          ELSE IF(nu.EQ.4) THEN
            PSI2_XI=4.0d0-4.0d0*XI(1)-8.0d0*XI(2)
          ELSE IF(nu.EQ.5) THEN
            PSI2_XI=-8.0d0
          ELSE IF(nu.EQ.6) THEN
            PSI2_XI=-4.0d0
          ELSE
            PSI2_XI=0.0d0
          ENDIF
       ELSE IF(nn.EQ.5) THEN
          IF(nu.EQ.1) THEN
            PSI2_XI=4.0d0*XI(1)*XI(2)
          ELSE IF(nu.EQ.2) THEN
            PSI2_XI=4.0d0*XI(2)
          ELSE IF(nu.EQ.4) THEN
            PSI2_XI=4.0d0*XI(1)
          ELSE IF(nu.EQ.6) THEN
            PSI2_XI=4.0d0
          ELSE
            PSI2_XI=0.0d0
          ENDIF
       ELSE IF(nn.EQ.6) THEN
          IF(nu.EQ.1) THEN
            PSI2_XI=XI(2)*(2.0d0*XI(2)-1.0d0)
          ELSE IF(nu.EQ.4) THEN
            PSI2_XI=4.0d0*XI(2)-1.0d0
          ELSE IF(nu.EQ.5) THEN
            PSI2_XI=4.0d0
          ELSE
            PSI2_XI=0.0d0
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END


