      REAL*8 FUNCTION D2ZX(ICOORD,i,j,k,X)

C#### Function: D2ZX
C###  Type: REAL*8
C###  Description:
C###    Calculates D2Z(I)/DX(J)DX(K) at X, where Z(I) are rectangalar
C###    Cartesian and X(J) are curvilinear coordinates defined by
C###    ICOORD.
C###    ICOORD=1 for rectangular cartesian coordinates
C###           2 for cylindrical polar coordinates
C###           3 for spherical polar coordinates
C###           4 for prolate spheriodal coordinates
C###           5 for oblate spheroidal coordinates

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
!     Parameter List
      INTEGER ICOORD,i,j,k
      REAL*8  X(3)

      GO TO (1000,2000,3000,4000,5000),ICOORD
        D2ZX=0.0D0
        RETURN
 1000   D2ZX=0.0D0
        RETURN
 2000   GO TO (2100,2200,2300),i
 2100     GO TO (2110,2120,2130),j
 2110       GO TO (2111,2112,2113),k
 2111         D2ZX=0.0D0
              RETURN
 2112         D2ZX=-DSIN(X(2))
              RETURN
 2113         D2ZX=0.0D0
              RETURN
 2120       GO TO (2112,2122,2123),k
 2122         D2ZX=-X(1)*DCOS(X(2))
              RETURN
 2123         D2ZX=0.0D0
              RETURN
 2130       GO TO (2113,2123,2133),k
 2133         D2ZX=0.0D0
              RETURN
 2200     GO TO (2210,2220,2230),j
 2210       GO TO (2211,2212,2213),k
 2211         D2ZX=0.0D0
              RETURN
 2212         D2ZX=DCOS(X(2))
              RETURN
 2213         D2ZX=0.0D0
              RETURN
 2220       GO TO (2212,2222,2223),k
 2222         D2ZX=-X(1)*DSIN(X(2))
              RETURN
 2223         D2ZX=0.0D0
              RETURN
 2230       GO TO (2213,2223,2233),k
 2233         D2ZX=0.0D0
              RETURN
 2300         D2ZX=0.0D0
              RETURN
 3000   GO TO (3100,3200,3300),i
 3100     GO TO (3110,3120,3130),j
 3110       GO TO (3111,3112,3113),k
 3111         D2ZX=0.0D0
              RETURN
 3112         D2ZX=-DSIN(X(2))*DCOS(X(3))
              RETURN
 3113         D2ZX=-DCOS(X(2))*DSIN(X(3))
              RETURN
 3120       GO TO (3112,3122,3123),k
 3122         D2ZX=-X(1)*DCOS(X(2))*DCOS(X(3))
              RETURN
 3123         D2ZX=X(1)*DSIN(X(2))*DSIN(X(3))
              RETURN
 3130       GO TO (3113,3123,3133),k
C PJH 29Aug95 3133         D2ZX=X(1)*DCOS(X(2))*DCOS(X(3))
 3133         D2ZX=-X(1)*DCOS(X(2))*DCOS(X(3))
              RETURN
 3200     GO TO (3210,3220,3230),j
 3210       GO TO (3211,3212,3213),k
 3211         D2ZX=0.0D0
              RETURN
 3212         D2ZX=DCOS(X(2))*DCOS(X(3))
              RETURN
 3213         D2ZX=-DSIN(X(2))*DSIN(X(3))
              RETURN
 3220       GO TO (3212,3222,3223),k
 3222         D2ZX=-X(1)*DSIN(X(2))*DCOS(X(3))
              RETURN
 3223         D2ZX=-X(1)*DCOS(X(2))*DSIN(X(3))
              RETURN
 3230       GO TO (3213,3223,3233),k
 3233         D2ZX=-X(1)*DSIN(X(2))*DCOS(X(3))
              RETURN
 3300     GO TO (3310,3320,3330),j
 3310       GO TO (3311,3312,3313),k
 3311         D2ZX=0.0D0
              RETURN
 3312         D2ZX=0.0D0
              RETURN
 3313         D2ZX=DCOS(X(3))
              RETURN
 3320       GO TO (3312,3322,3323),k
 3322         D2ZX=0.0D0
              RETURN
 3323         D2ZX=0.0D0
              RETURN
 3330       GO TO (3313,3323,3333),k
 3333         D2ZX=-X(1)*DSIN(X(3))
              RETURN
 4000   GO TO (4100,4200,4300),i
 4100     GO TO (4110,4120,4130),j
 4110       GO TO (4111,4112,4113),k
 4111         D2ZX=FOCUS*DCOSH(X(1))*DCOS(X(2))
              RETURN
 4112         D2ZX=-FOCUS*DSINH(X(1))*DSIN(X(2))
              RETURN
 4113         D2ZX=0.0D0
              RETURN
 4120       GO TO (4112,4122,4123),k
 4122         D2ZX=-FOCUS*DCOSH(X(1))*DCOS(X(2))
              RETURN
 4123         D2ZX=0.0D0
              RETURN
 4130       GO TO (4113,4123,4133),k
 4133         D2ZX=0.0D0
              RETURN
 4200     GO TO (4210,4220,4230),j
 4210       GO TO (4211,4212,4213),k
 4211         D2ZX=FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
              RETURN
 4212         D2ZX=FOCUS*DCOSH(X(1))*DCOS(X(2))*DCOS(X(3))
              RETURN
 4213         D2ZX=-FOCUS*DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))
              RETURN
 4220       GO TO (4212,4222,4223),k
 4222         D2ZX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
              RETURN
 4223         D2ZX=-FOCUS*DSINH(X(1))*DCOS(X(2))*DSIN(X(3))
              RETURN
 4230       GO TO (4213,4223,4233),k
 4233         D2ZX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
              RETURN
 4300     GO TO (4310,4320,4330),j
 4310       GO TO (4311,4312,4313),k
 4311         D2ZX=FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
              RETURN
 4312         D2ZX=FOCUS*DCOSH(X(1))*DCOS(X(2))*DSIN(X(3))
              RETURN
 4313         D2ZX=FOCUS*DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))
              RETURN
 4320       GO TO (4312,4322,4323),k
 4322         D2ZX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
              RETURN
 4323         D2ZX=FOCUS*DSINH(X(1))*DCOS(X(2))*DCOS(X(3))
              RETURN
 4330       GO TO (4313,4323,4333),k
 4333         D2ZX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
              RETURN
 5000   GO TO (5100,5200,5300),i
 5100         D2ZX=0.0D0
              RETURN
 5200         D2ZX=0.0D0
              RETURN
 5300         D2ZX=0.0D0
              RETURN
      END


