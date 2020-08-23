      REAL*8 FUNCTION DZX(ICOORD,i,j,X)

C#### Function: DZX
C###  Type: REAL*8
C###  Description:
C###    <HTML>
C###    DZX calculates DZ(i)/DX(j) at X, where Z(i) are rectangular
C###    Cartesian and X(j) are curvilinear coordinates defined by
C###    ICOORD.
C###    <PRE>
C###    ICOORD=1 for rectangular cartesian coordinates;
C###           2 for cylindrical polar coordinates;
C###           3 for spherical polar coordinates;
C###           4 for prolate spheriodal coordinates; and
C###           5 for oblate spheroidal coordinates.
C###    </PRE></HTML>

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
!     Parameter List
      INTEGER ICOORD,i,j
      REAL*8  X(3)

      GO TO (100,200,300,400,500),ICOORD
        DZX=0.0D0
        RETURN
  100   DZX=0.0D0
        IF(i.EQ.j) DZX=1.0D0
        RETURN
  200   GO TO (210,220,230),i
  210     GO TO (211,212,213),j
  211       DZX=DCOS(X(2))
            RETURN
  212       DZX=-X(1)*DSIN(X(2))
            RETURN
  213       DZX=0.0D0
            RETURN
  220     GO TO (221,222,223),j
  221       DZX=DSIN(X(2))
            RETURN
  222       DZX=X(1)*DCOS(X(2))
            RETURN
  223       DZX=0.0D0
            RETURN
  230     GO TO (231,232,233),j
  231       DZX=0.0D0
            RETURN
  232       DZX=0.0D0
            RETURN
  233       DZX=1.0D0
            RETURN
  300   GO TO (310,320,330),i
  310     GO TO (311,312,313),j
  311       DZX=DCOS(X(2))*DCOS(X(3))
            RETURN
  312       DZX=-X(1)*DSIN(X(2))*DCOS(X(3))
            RETURN
  313       DZX=-X(1)*DCOS(X(2))*DSIN(X(3))
            RETURN
  320     GO TO (321,322,323),j
  321       DZX=DSIN(X(2))*DCOS(X(3))
            RETURN
  322       DZX=X(1)*DCOS(X(2))*DCOS(X(3))
            RETURN
  323       DZX=-X(1)*DSIN(X(2))*DSIN(X(3))
            RETURN
  330     GO TO (331,332,333),j
  331       DZX=DSIN(X(3))
            RETURN
  332       DZX=0.0D0
            RETURN
  333       DZX=X(1)*DCOS(X(3))
            RETURN
  400   GO TO (410,420,430),i
  410     GO TO (411,412,413),j
  411       DZX=FOCUS*DSINH(X(1))*DCOS(X(2))
            RETURN
  412       DZX=-FOCUS*DCOSH(X(1))*DSIN(X(2))
            RETURN
  413       DZX=0.0D0
            RETURN
  420     GO TO (421,422,423),j
  421       DZX=FOCUS*DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))
            RETURN
  422       DZX=FOCUS*DSINH(X(1))*DCOS(X(2))*DCOS(X(3))
            RETURN
  423       DZX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
            RETURN
  430     GO TO (431,432,433),j
  431       DZX=FOCUS*DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))
            RETURN
  432       DZX=FOCUS*DSINH(X(1))*DCOS(X(2))*DSIN(X(3))
            RETURN
  433       DZX=FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
            RETURN
  500   GO TO (510,520,530),i
  510     GO TO (511,512,513),j
  511       DZX= FOCUS*DSINH(X(1))*DCOS(X(2))*DCOS(X(3))
            RETURN
  512       DZX=-FOCUS*DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))
            RETURN
  513       DZX=-FOCUS*DCOSH(X(1))*DCOS(X(2))*DSIN(X(3))
            RETURN
  520     GO TO (521,522,523),j
  521       DZX= FOCUS*DCOSH(X(1))*DSIN(X(2))
            RETURN
  522       DZX= FOCUS*DSINH(X(1))*DCOS(X(2))
            RETURN
  523       DZX= 0.0D0
            RETURN
  530     GO TO (531,532,533),j
  531       DZX= FOCUS*DSINH(X(1))*DCOS(X(2))*DSIN(X(3))
            RETURN
  532       DZX=-FOCUS*DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))
            RETURN
  533       DZX= FOCUS*DCOSH(X(1))*DCOS(X(2))*DCOS(X(3))
            RETURN
      END


