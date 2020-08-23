      REAL*8 FUNCTION DXZ(ICOORD,i,k,X)

C#### Function: DXZ
C###  Type: REAL*8
C###  Description:
C###    DXZ calculates DX(I)/DZ(K) at X, where Z(K) are rectangular
C###    Cartesian and X(I) are curvilinear coordinates defined by
C###    ICOORD.
C###    ICOORD=1 for rectangular cartesian coordinates;
C###           2 for cylindrical polar coordinates;
C###           3 for spherical polar coordinates;
C###           4 for prolate spheriodal coordinates; and
C###           5 for oblate spheroidal coordinates.

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
!     Parameter List
      INTEGER ICOORD,i,k
      REAL*8 RD,X(3)

      GO TO (100,200,300,400,500),ICOORD
        DXZ=0.0D0
        RETURN
  100   DXZ=0.0D0
        IF(i.EQ.k) DXZ=1.0D0
        RETURN
  200   GO TO (210,220,230),i
  210     GO TO (211,212,213),k
  211       DXZ=DCOS(X(2))
            RETURN
  212       DXZ=DSIN(X(2))
            RETURN
  213       DXZ=0.0D0
            RETURN
  220     GO TO (221,222,223),k
  221       DXZ=-DSIN(X(2))/X(1)
            RETURN
  222       DXZ=DCOS(X(2))/X(1)
            RETURN
  223       DXZ=0.0D0
            RETURN
  230     GO TO (231,232,233),k
  231       DXZ=0.0D0
            RETURN
  232       DXZ=0.0D0
            RETURN
  233       DXZ=1.0D0
            RETURN
  300   GO TO (310,320,330),i
  310     GO TO (311,312,313),k
  311       DXZ=DCOS(X(2))*DCOS(X(3))
            RETURN
  312       DXZ=DSIN(X(2))*DCOS(X(3))
            RETURN
  313       DXZ=DSIN(X(3))
            RETURN
  320     GO TO (321,322,323),k
  321       DXZ=-DSIN(X(2))/(X(1)*DCOS(X(3)))
            RETURN
  322       DXZ=DCOS(X(2))/(X(1)*DCOS(X(3)))
            RETURN
  323       DXZ=0.0D0
            RETURN
  330     GO TO (331,332,333),k
  331       DXZ=-DCOS(X(2))*DSIN(X(3))/X(1)
            RETURN
  332       DXZ=-DSIN(X(2))*DSIN(X(3))/X(1)
            RETURN
  333       DXZ=DCOS(X(3))/X(1)
            RETURN
  400   RD=FOCUS*(DCOSH(X(1))*DCOSH(X(1))-DCOS(X(2))*DCOS(X(2)))
        GO TO (410,420,430),i
  410     GO TO (411,412,413),k
  411       DXZ=DSINH(X(1))*DCOS(X(2))/RD
            RETURN
  412       DXZ=DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))/RD
            RETURN
  413       DXZ=DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))/RD
            RETURN
  420     GO TO (421,422,423),k
  421       DXZ=-DCOSH(X(1))*DSIN(X(2))/RD
            RETURN
  422       DXZ=DSINH(X(1))*DCOS(X(2))*DCOS(X(3))/RD
            RETURN
  423       DXZ=DSINH(X(1))*DCOS(X(2))*DSIN(X(3))/RD
            RETURN
  430     GO TO (431,432,433),k
  431       DXZ=0.0D0
            RETURN
  432       DXZ=-DSIN(X(3))/(FOCUS*DSINH(X(1))*DSIN(X(2)))
            RETURN
  433       DXZ=DCOS(X(3))/(FOCUS*DSINH(X(1))*DSIN(X(2)))
            RETURN
  500   GO TO (510,520,530),i
  510     GO TO (511,512,513),k
  511       DXZ=0.0D0
            RETURN
  512       DXZ=0.0D0
            RETURN
  513       DXZ=0.0D0
            RETURN
  520     GO TO (521,522,523),k
  521       DXZ=0.0D0
            RETURN
  522       DXZ=0.0D0
            RETURN
  523       DXZ=0.0D0
            RETURN
  530     GO TO (531,532,533),k
  531       DXZ=0.0D0
            RETURN
  532       DXZ=0.0D0
            RETURN
  533       DXZ=0.0D0
            RETURN
      END


