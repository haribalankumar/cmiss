      SUBROUTINE ZX(ICOORD,Z,X)

C#### Subroutine: ZX
C###  Description:
C###    ZX performs a coordinate transformation from rectangular
C###    cartesian coordinates at the point with coordinates Z to the
C###    coordinate system identified by ICOORD.
C**** ICOORD=1 for rectangular cartesian coordinates
C****        2 for cylindrical polar coordinates
C****        3 for spherical polar coordinates
C****        4 for prolate spheriodal coordinates
C****        5 for oblate spheroidal coordinates

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ICOORD
      REAL*8 X(3),Z(3)
!     Local Variables
      CHARACTER ERROR*10
      REAL*8 A1,A2,A3,A4,A5,A6,A7,A8,A9,DATAN_MOD,Z2(3)

C     CALL ENTERS('ZX',*9999)
      Z2(1)=Z(1)
      Z2(2)=Z(2)
      Z2(3)=Z(3)
      GO TO (1,2,3,4,5),ICOORD
    1   X(1)=Z2(1)
        X(2)=Z2(2)
        X(3)=Z2(3)
        GO TO 9999
    2   X(1)=DSQRT(Z2(1)**2+Z2(2)**2)
C cpb 27/1/96 Replacing this by a function
C        IF((Z2(1).NE.0.0D0).AND.(Z2(2).NE.0.0D0))THEN
C          X(2)=DATAN2(Z2(2),Z2(1))
C        ELSE IF(Z2(2).EQ.0.0D0.AND.Z2(1).GT.0.0D0)THEN     !AAY 11JUN90
C          X(2)=0.0D0
C        ELSE IF(Z2(2).EQ.0.0D0.AND.Z2(1).LT.0.0D0)THEN     !AAY 11JUN90
C          X(2)=PI
C        ELSE IF(Z2(1).EQ.0.0D0.AND.Z2(2).GT.0.0D0)THEN     !AAY 11JUN90
C          X(2)=PI/2.0D0
C        ELSE IF(Z2(1).EQ.0.0D0.AND.Z2(2).LT.0.0D0)THEN     !AAY 11JUN90
C          X(2)=-PI/2.0D0
C        ENDIF
        X(2)=DATAN_MOD(Z2(1),Z2(2))
        IF(X(2).LT.0.0D0) X(2)=X(2)+2.0D0*PI !ref coords in 0 -> 2*pi
        X(3)=Z2(3)
        GO TO 9999
    3   X(1)=DSQRT(Z2(1)**2+Z2(2)**2+Z2(3)**2)
        IF((Z2(1).NE.0.0D0).OR.(Z2(2).NE.0.0D0))THEN
          X(2)=DATAN2(Z2(2),Z2(1))
        ELSE
          X(2)=0.0D0
        ENDIF
        A1=DSQRT(Z2(1)**2+Z2(2)**2)
        IF((Z2(3).NE.0.0D0).OR.(a1.NE.0.0D0))THEN
          X(3)=DATAN2(Z2(3),A1)
        ELSE
          X(3)=0.0D0
        ENDIF
        GO TO 9999
    4   A1=Z2(1)**2+Z2(2)**2+Z2(3)**2-FOCUS**2
        A2=DSQRT(A1**2+4.0D0*(FOCUS**2)*(Z2(2)**2+Z2(3)**2))
        A3=2.0D0*(FOCUS**2)
        A4=DMAX1((A2+A1)/A3,0.0D0)
        A5=DMAX1((A2-A1)/A3,0.0D0)
        A6=DSQRT(A4)
        A7=DMIN1(DSQRT(A5),1.0D0)
        IF(DABS(A7).LE.1.0D0) THEN
          A8=DASIN(A7)
        ELSE
          A8=0.0D0
          WRITE(OP_STRING,'('' Warning: Put A8=0 in ZX since '','
     '      //'''abs(A8)>1'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF((Z2(3).EQ.0.0D0).OR.(A6.EQ.0.0D0).OR.(A7.EQ.0.0D0)) THEN
          A9=0.0D0
        ELSE
          IF(DABS(A6*A7).GT.0.0D0) THEN
            A9=Z2(3)/(FOCUS*A6*A7)
          ELSE
            A9=0.0D0
            WRITE(OP_STRING,'('' Warning: Put A9=0 in ZX since '','
     '        //'''A6*A7=0'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          IF(A9.GE.1.0D0) then
            A9=PI/2.0D0
          ELSE IF(A9.LE.-1.0D0) THEN
            A9=-PI/2.0D0
          ELSE
            A9=DASIN(A9)
          ENDIF
        ENDIF
        X(1)=DLOG(A6+DSQRT(A4+1.0D0))
        IF(Z2(1).GE.0.0D0) THEN
          X(2)=A8
        ELSE
          X(2)=PI-A8
        ENDIF
        IF(Z2(2).GE.0.0D0) THEN
          X(3)=DMOD(A9+2.0D0*PI,2.0D0*PI)
        ELSE
          X(3)=PI-A9
        ENDIF
 5      CONTINUE
C99     CALL EXITS('ZX')
 9999   RETURN
      END


