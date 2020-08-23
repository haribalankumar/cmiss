      SUBROUTINE XZ_DERIV(ICOORD,nu,XG,XGRC)

C#### Subroutine: XZ_DERIV
C###  Description:
C###    XZ_DERIV performs a coordinate transformation from the
C###    coordinate system identified by ICOORD at the point with
C###    coordinates/derivs XG(nj,nu) to rect cart coordinates/derivs
C###    XGRC(nj,nu).
C**** nu=1 value
C****    2 1st deriv wrt s1
C****    3 2nd deriv wrt s1
C****    4 1st deriv wrt s2
C****    5 2nd deriv wrt s2
C****    6 2nd deriv wrt s1,s2
C****    7 1st deriv wrt s3
C****    8 2nd deriv wrt s3
C****    9 2nd deriv wrt s1,s3
C****   10 2nd deriv wrt s2,s3
C****   11 3rd deriv wrt s1,s2,s3

C**** ICOORD=1 for rectangular cartesian coordinates
C****        2 for cylindrical polar coordinates
C****        3 for spherical polar coordinates
C****        4 for prolate spheriodal coordinates
C****        5 for oblate spheroidal coordinates

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ICOORD,nj,nu
      REAL*8 XG(NJM,NUM),XGRC(NJM,NUM)
!     Local Variables
      REAL*8 ZG_LOCAL(3,11)

C     case icoord of
      GO TO (1,2,3,4,5),ICOORD
    1   IF(NJT.EQ.1) THEN !rect cart
          ZG_LOCAL(1,nu)=XG(1,nu)
        ELSE IF(NJT.EQ.2) THEN
          ZG_LOCAL(1,nu)=XG(1,nu)
          ZG_LOCAL(2,nu)=XG(2,nu)
        ELSE IF(NJT.EQ.3) THEN
          ZG_LOCAL(1,nu)=XG(1,nu)
          ZG_LOCAL(2,nu)=XG(2,nu)
          ZG_LOCAL(3,nu)=XG(3,nu)
        ENDIF
        GO TO 99

    2   GO TO (21,22,23,24,25,26,27,28,29,210,211),nu  !cylin polar
   21     CALL XZ(ICOORD,XG(1,nu),ZG_LOCAL(1,nu))           !value
          GO TO 99
   22     ZG_LOCAL(1,2)=DCOS(XG(2,1))*XG(1,2)-XG(1,1)*
     '      DSIN(XG(2,1))*XG(2,2)                     !d(x)/d(s1)
          ZG_LOCAL(2,2)=DSIN(XG(2,1))*XG(1,2)+XG(1,1)*
     '      DCOS(XG(2,1))*XG(2,2)                     !d(y)/d(s1)
          ZG_LOCAL(3,2)=XG(3,2)                       !d(z)/d(s1)
          GO TO 99
   23     GO TO 99                                    !d2(x)/d(s1)2
   24     ZG_LOCAL(1,4)=DCOS(XG(2,1))*XG(1,4)-XG(1,1)*
     '      DSIN(XG(2,1))*XG(2,4)
!d(x)/d(s2)
          ZG_LOCAL(2,4)=DSIN(XG(2,1))*XG(1,4)+XG(1,1)*
     '      DCOS(XG(2,1))*XG(2,4)
!d(y)/d(s2)
          ZG_LOCAL(3,4)=XG(3,4)                       !d(z)/d(s2)
          GO TO 99
   25     GO TO 99                                    !d2(x)/d(s2)2
   26     ZG_LOCAL(1,6)=DCOS(XG(2,1))*XG(1,6)-XG(1,2)*
     '      DSIN(XG(2,1))*XG(2,4)
     '      -(XG(1,4)*DSIN(XG(2,1))*XG(2,2)+XG(1,1)*DCOS(XG(2,1))*
     '        XG(2,2)*XG(2,4)+XG(1,1)*DSIN(XG(2,1))*XG(2,6))
!d2(x)/d(s1)d(s2)
          ZG_LOCAL(2,6)=DSIN(XG(2,1))*XG(1,6)+XG(1,2)*
     '      DCOS(XG(2,1))*XG(2,4)
     '      +(XG(1,4)*DCOS(XG(2,1))*XG(2,2)-XG(1,1)*DSIN(XG(2,1))*
     '        XG(2,2)*XG(2,4)+XG(1,1)*DCOS(XG(2,1))*XG(2,6))
!d2(y)/d(s1)d(s2)
          ZG_LOCAL(3,6)=XG(3,6)                       !d2(z)/d(s1)d(s2)
          GO TO 99
   27     ZG_LOCAL(1,7)=DCOS(XG(2,1))*XG(1,7)-XG(1,1)*
     '      DSIN(XG(2,1))*XG(2,7)
!d(x)/d(s3)
          ZG_LOCAL(2,7)=DSIN(XG(2,1))*XG(1,7)+XG(1,1)*
     '      DCOS(XG(2,1))*XG(2,7)
!d(y)/d(s3)
          ZG_LOCAL(3,7)=XG(3,7)                       !d(z)/d(s3)
          GO TO 99
   28     GO TO 99                                    !d2(x)/d(s3)2
   29     ZG_LOCAL(1,9)=DCOS(XG(2,1))*XG(1,9)-XG(1,2)*
     '      DSIN(XG(2,1))*XG(2,7)
     '      -(XG(1,7)*DSIN(XG(2,1))*XG(2,2)+XG(1,1)*DCOS(XG(2,1))*
     '        XG(2,2)*XG(2,7)+XG(1,1)*DSIN(XG(2,1))*XG(2,9))
!d2(x)/d(s1)d(s3)
          ZG_LOCAL(2,9)=DSIN(XG(2,1))*XG(1,9)+XG(1,2)*
     '      DCOS(XG(2,1))*XG(2,7)
     '      +(XG(1,7)*DCOS(XG(2,1))*XG(2,2)-XG(1,1)*DSIN(XG(2,1))*
     '        XG(2,2)*XG(2,7)+XG(1,1)*DCOS(XG(2,1))*XG(2,9))
!d2(y)/d(s1)d(s3)
          ZG_LOCAL(3,9)=XG(3,9)                       !d2(z)/d(s1)d(s3)
          GO TO 99
  210     ZG_LOCAL(1,10)=DCOS(XG(2,1))*XG(1,10)-XG(1,4)*
     '      DSIN(XG(2,1))*XG(2,7)
     '      -(XG(1,7)*DSIN(XG(2,1))*XG(2,4)+XG(1,1)*DCOS(XG(2,1))*
     '        XG(2,4)*XG(2,7)+XG(1,1)*DSIN(XG(2,1))*XG(2,10))
!d2(x)/d(s2)d(s3)
          ZG_LOCAL(2,10)=DSIN(XG(2,1))*XG(1,10)+XG(1,4)*
     '      DCOS(XG(2,1))*XG(2,7)
     '      +(XG(1,7)*DCOS(XG(2,1))*XG(2,4)-XG(1,1)*DSIN(XG(2,1))*
     '        XG(2,4)*XG(2,7)+XG(1,1)*DCOS(XG(2,1))*XG(2,10))
!d2(y)/d(s2)d(s3)
          ZG_LOCAL(3,10)=XG(3,10)                     !d2(z)/d(s2)d(s3)
          GO TO 99
  211     ZG_LOCAL(1,11)=                    !d3(x)/d(s1)d(s2)d(s3)
     '      -DCOS(XG(2,1))*XG(2,2)*XG(2,4)*XG(1,7)-
     '       DSIN(XG(2,1))*(XG(2,6)*XG(1,7)+XG(2,4)*XG(1,9))-
     '       DSIN(XG(2,1))*XG(2,2)*XG(1,10)+
     '       DCOS(XG(2,1))*XG(1,11)-
     '       DCOS(XG(2,1))*XG(2,2)*XG(1,4)*XG(2,7)-
     '       DSIN(XG(2,1))*(XG(1,6)*XG(2,7)+XG(1,4)*XG(2,9))+
     '     (-DCOS(XG(2,1))*XG(1,2)+DSIN(XG(2,1))*XG(1,1)*XG(2,2))*
     '        XG(2,4)*XG(2,7)-XG(1,1)*DCOS(XG(2,1))*
     '        (XG(2,6)*XG(2,7)+XG(2,4)*XG(2,9))+
     '     (-DSIN(XG(2,1))*XG(1,2)-XG(1,1)*
     '       DCOS(XG(2,1))*XG(2,2))*XG(2,10)-
     '       XG(1,1)*DSIN(XG(2,1))*XG(2,11)
          ZG_LOCAL(2,11)=                     !d3(y)/d(s1)d(s2)d(s3)
     '      -DSIN(XG(2,1))*XG(2,2)*XG(2,4)*XG(1,7)+
     '       DCOS(XG(2,1))*(XG(2,6)*XG(1,7)+XG(2,4)*XG(1,9))+
     '       DCOS(XG(2,1))*XG(2,2)*XG(1,10)+
     '       DSIN(XG(2,1))*XG(1,11)-
     '       DSIN(XG(2,1))*XG(2,2)*XG(1,4)*XG(2,7)+
     '       DCOS(XG(2,1))*(XG(1,6)*XG(2,7)+XG(1,4)*XG(2,9))+
     '     (-DSIN(XG(2,1))*XG(1,2)-XG(1,1)*DCOS(XG(2,1))*
     '       XG(2,2))*XG(2,4)*XG(2,7)-
     '       XG(1,1)*DSIN(XG(2,1))*(XG(2,6)*XG(2,7)+XG(2,4)*XG(2,9))+
     '     ( DCOS(XG(2,1))*XG(1,2)-XG(1,1)*DSIN(XG(2,1))*XG(2,2))*
     '       XG(2,10)+XG(1,1)*DCOS(XG(2,1))*XG(2,11)
          ZG_LOCAL(3,11)=XG(3,11)             !d3(z)/d(s1)d(s2)d(s3)
          GO TO 99

    3   GO TO (31,31,31,31,31,31),nu  !spherical polar
   31     CALL XZ(ICOORD,XG(1,nu),ZG_LOCAL(1,nu))           !value

    4   GO TO (41,42,43,44,45,46,47,48,49,410,411),nu  !prolate spheroid
   41     CALL XZ(ICOORD,XG(1,nu),ZG_LOCAL(1,nu))           !value
          GO TO 99
   42     ZG_LOCAL(1,nu)=FOCUS*DSINH(XG(1,1))*DCOS(XG(2,1))*XG(1,nu)
     '      -FOCUS*DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(2,nu)
!d(x)/d(s1)
          ZG_LOCAL(2,nu)=FOCUS*DCOSH(XG(1,1))*DSIN(XG(2,1))*
     '      DCOS(XG(3,1))*XG(1,nu)
     '      +FOCUS*DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,nu)
     '      -FOCUS*DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,nu)
!d(y)/d(s1)
          ZG_LOCAL(3,nu)=FOCUS*DCOSH(XG(1,1))*DSIN(XG(2,1))*
     '      DSIN(XG(3,1))*XG(1,nu)
     '      +FOCUS*DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,nu)
     '      +FOCUS*DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,nu)
!d(z)/d(s1)
          GO TO 99
   43     GO TO 99                                    !d2(x)/d(s1)2
   44     ZG_LOCAL(1,nu)=FOCUS*DSINH(XG(1,1))*DCOS(XG(2,1))*XG(1,nu)
     '      -FOCUS*DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(2,nu)
!d(x)/d(s2)
          ZG_LOCAL(2,nu)=FOCUS*DCOSH(XG(1,1))*DSIN(XG(2,1))*
     '      DCOS(XG(3,1))*XG(1,nu)
     '      +FOCUS*DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,nu)
     '      -FOCUS*DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,nu)
!d(y)/d(s2)
          ZG_LOCAL(3,nu)=FOCUS*DCOSH(XG(1,1))*DSIN(XG(2,1))*
     '      DSIN(XG(3,1))*XG(1,nu)
     '      +FOCUS*DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,nu)
     '      +FOCUS*DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,nu)
!d(z)/d(s2)
          GO TO 99
   45     GO TO 99                                    !d2(x)/d(s2)2
   46     ZG_LOCAL(1,6)=FOCUS*                        !d2(x)/d(s1)d(s2)
     '      (DSINH(XG(1,1))*DCOS(XG(2,1))*XG(1,6)+
     '      XG(1,2)*(DCOSH(XG(1,1))*DCOS(XG(2,1))*XG(1,4)-
     '               DSINH(XG(1,1))*DSIN(XG(2,1))*XG(2,4))
     '      -DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(2,6)-
     '      XG(2,2)*(DSINH(XG(1,1))*DSIN(XG(2,1))*XG(1,4)+
     '               DCOSH(XG(1,1))*DCOS(XG(2,1))*XG(2,4)))
          ZG_LOCAL(2,6)=FOCUS*                      !d2(y)/d(s1)d(s2)
     '      (DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,6)+
     '        XG(1,2)*(DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(1,4)+DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,4)-DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,4))+DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,6)+
     '        XG(2,2)*(DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(1,4)-DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,4)-DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,4))-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,6)-
     '        XG(3,2)*(DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(1,4)+DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,4)+DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,4)))
          ZG_LOCAL(3,6)=FOCUS*                      !d2(z)/d(s1)d(s2)
     '      (DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,6)+
     '        XG(1,2)*(DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(1,4)+DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,4)+DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,4))+DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,6)+
     '        XG(2,2)*(DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(1,4)-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,4)+DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,4))+DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,6)+
     '        XG(3,2)*(DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(1,4)+DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,4)-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,4)))
          GO TO 99
   47     ZG_LOCAL(1,nu)=FOCUS*DSINH(XG(1,1))*DCOS(XG(2,1))*XG(1,nu)
     '      -FOCUS*DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(2,nu)
!d(x)/d(s3)
          ZG_LOCAL(2,nu)=FOCUS*DCOSH(XG(1,1))*DSIN(XG(2,1))*
     '      DCOS(XG(3,1))*XG(1,nu)
     '      +FOCUS*DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,nu)
     '      -FOCUS*DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,nu)
!d(y)/d(s3)
          ZG_LOCAL(3,nu)=FOCUS*DCOSH(XG(1,1))*DSIN(XG(2,1))*
     '      DSIN(XG(3,1))*XG(1,nu)
     '      +FOCUS*DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,nu)
     '      +FOCUS*DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,nu)
!d(z)/d(s3)
          GO TO 99
   48     GO TO 99                                    !d2(x)/d(s3)2
   49     ZG_LOCAL(1,9)=FOCUS*                        !d2(x)/d(s1)d(s3)
     '      (DSINH(XG(1,1))*DCOS(XG(2,1))*XG(1,9)+
     '      XG(1,2)*(DCOSH(XG(1,1))*DCOS(XG(2,1))*XG(1,7)-
     '               DSINH(XG(1,1))*DSIN(XG(2,1))*XG(2,7))
     '      -DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(2,9)-
     '      XG(2,2)*(DSINH(XG(1,1))*DSIN(XG(2,1))*XG(1,7)+
     '               DCOSH(XG(1,1))*DCOS(XG(2,1))*XG(2,7)))
          ZG_LOCAL(2,9)=FOCUS*                   !d2(y)/d(s1)d(s3)
     '      (DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,9)+
     '        XG(1,2)*(DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(1,7)+DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,7)-DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,7))+DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,9)+
     '        XG(2,2)*(DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(1,7)-DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,7)-DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,7))-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,9)-
     '        XG(3,2)*(DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(1,7)+DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,7)+DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,7)))
          ZG_LOCAL(3,9)=FOCUS*                   !d2(z)/d(s1)d(s3)
     '      (DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,9)+
     '        XG(1,2)*(DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(1,7)+DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,7)+DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,7))+DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,9)+
     '        XG(2,2)*(DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(1,7)-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,7)+DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,7))+DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,9)+
     '        XG(3,2)*(DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(1,7)+DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,7)-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,7)))
          GO TO 99
  410     ZG_LOCAL(1,10)=FOCUS*                !d2(x)/d(s2)d(s3)
     '      (DSINH(XG(1,1))*DCOS(XG(2,1))*XG(1,10)+
     '      XG(1,4)*(DCOSH(XG(1,1))*DCOS(XG(2,1))*XG(1,7)-
     '               DSINH(XG(1,1))*DSIN(XG(2,1))*XG(2,7))
     '      -DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(2,10)-
     '      XG(2,4)*(DSINH(XG(1,1))*DSIN(XG(2,1))*XG(1,7)+
     '               DCOSH(XG(1,1))*DCOS(XG(2,1))*XG(2,7)))
          ZG_LOCAL(2,10)=FOCUS*                !d2(y)/d(s2)d(s3)
     '      (DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,10)+
     '        XG(1,4)*(DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(1,7)+DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,7)-DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,7))+DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,10)+
     '        XG(2,4)*(DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(1,7)-DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,7)-DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,7))-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,10)-
     '        XG(3,4)*(DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(1,7)+DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,7)+DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,7)))
          ZG_LOCAL(3,10)=FOCUS*                !d2(z)/d(s2)d(s3)
     '      (DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,10)+
     '        XG(1,4)*(DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(1,7)+DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,7)+DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,7))+DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,10)+
     '        XG(2,4)*(DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '          XG(1,7)-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(2,7)+DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,7))+DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(3,10)+
     '        XG(3,4)*(DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '          XG(1,7)+DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '          XG(2,7)-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '          XG(3,7)))
          GO TO 99
  411     ZG_LOCAL(1,11)=FOCUS*                   !d3(x)/d(s1)d(s2)d(s3)
     '      ((DSINH(XG(1,1))*DCOS(XG(2,1))*XG(1,2)-
     '        DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(2,2))*XG(1,4)*XG(1,7)+
     '        DCOSH(XG(1,1))*DCOS(XG(2,1))*
     '          (XG(1,6)*XG(1,7)+XG(1,4)*XG(1,9))+
     '      (-DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(1,2)-
     '        DSINH(XG(1,1))*DCOS(XG(2,1))*XG(2,2))*XG(2,4)*XG(1,7)-
     '        DSINH(XG(1,1))*DSIN(XG(2,1))*
     '          (XG(2,6)*XG(1,7)+XG(2,4)*XG(1,9))+
     '      ( DCOSH(XG(1,1))*DCOS(XG(2,1))*XG(1,2)-
     '        DSINH(XG(1,1))*DSIN(XG(2,1))*XG(2,2))*XG(1,10)+
     '        DSINH(XG(1,1))*DCOS(XG(2,1))*XG(1,11)+
     '      (-DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(1,2)-
     '        DSINH(XG(1,1))*DCOS(XG(2,1))*XG(2,2))*XG(1,4)*XG(2,7)-
     '        DSINH(XG(1,1))*DSIN(XG(2,1))*
     '          (XG(1,6)*XG(2,7)+XG(1,4)*XG(2,9))+
     '      (-DSINH(XG(1,1))*DCOS(XG(2,1))*XG(1,2)+
     '        DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(2,2))*XG(2,4)*XG(2,7)-
     '        DCOSH(XG(1,1))*DCOS(XG(2,1))*
     '          (XG(2,6)*XG(2,7)+XG(2,4)*XG(2,9))+
     '      (-DSINH(XG(1,1))*DSIN(XG(2,1))*XG(1,2)-
     '        DCOSH(XG(1,1))*DCOS(XG(2,1))*XG(2,2))*XG(2,10)-
     '        DCOSH(XG(1,1))*DSIN(XG(2,1))*XG(2,11))
          ZG_LOCAL(2,11)=                         !d3(y)/d(s1)d(s2)d(s3)
     '      FOCUS*
     '        ((DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,2)+
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,2)-
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '            XG(1,4)*XG(1,7)+
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '            (XG(1,6)*XG(1,7)+XG(1,4)*XG(1,9))+
     '        ( DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(1,2)-
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(2,2)-
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '            XG(2,4)*XG(1,7)+
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*(XG(2,6)*
     '            XG(1,7)+XG(2,4)*XG(1,9))+
     '        (-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,2)-
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,2)-
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '            XG(3,4)*XG(1,7)-
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*(XG(3,6)*
     '            XG(1,7)+XG(3,4)*XG(1,9))+
     '        ( DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,2)+
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,2)-
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '            XG(1,10)+
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,11))+
     '     FOCUS*
     '       ((DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(1,2)-
     '         DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(2,2)-
     '         DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '           XG(1,4)*XG(2,7)+
     '         DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '           (XG(1,6)*XG(2,7)+XG(1,4)*XG(2,9))+
     '       (-DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,2)-
     '         DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,2)+
     '         DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '           XG(2,4)*XG(2,7)-
     '         DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '           (XG(2,6)*XG(2,7)+XG(2,4)*XG(2,9))+
     '       (-DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(1,2)+
     '         DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(2,2)-
     '         DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '           XG(3,4)*XG(2,7)-
     '         DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '           (XG(3,6)*XG(2,7)+XG(3,4)*XG(2,9))+
     '       ( DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(1,2)-
     '         DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(2,2)-
     '         DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '           XG(2,10)+
     '         DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,11))+
     '     FOCUS*
     '      ((-DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,2)-
     '         DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,2)-
     '         DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '           XG(1,4)*XG(3,7)-
     '         DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '           (XG(1,6)*XG(3,7)+XG(1,4)*XG(3,9))+
     '       (-DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(1,2)+
     '         DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(2,2)-
     '         DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '           XG(2,4)*XG(3,7)-
     '         DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '           (XG(2,6)*XG(3,7)+XG(2,4)*XG(3,9))+
     '       (-DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,2)-
     '         DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,2)+
     '         DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '           XG(3,4)*XG(3,7)-
     '         DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '           (XG(3,6)*XG(3,7)+XG(3,4)*XG(3,9))+
     '       (-DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,2)-
     '         DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,2)-
     '         DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '           XG(3,10)-
     '         DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,11))
          ZG_LOCAL(3,11)=                         !d3(z)/d(s1)d(s2)d(s3)
     '      FOCUS*
     '        ((DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,2)+
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,2)+
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '            XG(1,4)*XG(1,7)+
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '            (XG(1,6)*XG(1,7)+XG(1,4)*XG(1,9))+
     '        ( DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(1,2)-
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(2,2)+
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '            XG(2,4)*XG(1,7)+
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '            (XG(2,6)*XG(1,7)+XG(2,4)*XG(1,9))+
     '        ( DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,2)+
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,2)-
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '            XG(3,4)*XG(1,7)+
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '            (XG(3,6)*XG(1,7)+XG(3,4)*XG(1,9))+
     '        ( DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,2)+
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,2)+
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '            XG(1,10)+
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,11))+
     '      FOCUS*
     '        ((DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(1,2)-
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(2,2)+
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '            XG(1,4)*XG(2,7)+
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*
     '            (XG(1,6)*XG(2,7)+XG(1,4)*XG(2,9))+
     '        (-DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,2)-
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,2)-
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '            XG(2,4)*XG(2,7)-
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '            (XG(2,6)*XG(2,7)+XG(2,4)*XG(2,9))+
     '        ( DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(1,2)-
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(2,2)-
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '            XG(3,4)*XG(2,7)+
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '            (XG(3,6)*XG(2,7)+XG(3,4)*XG(2,9))+
     '        ( DCOSH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(1,2)-
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(2,2)+
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '            XG(2,10)+
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,11))+
     '      FOCUS*
     '        ((DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,2)+
     '          DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,2)-
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '            XG(1,4)*XG(3,7)+
     '          DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*
     '            (XG(1,6)*XG(3,7)+XG(1,4)*XG(3,9))+
     '        ( DCOSH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(1,2)-
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(2,2)-
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '            XG(2,4)*XG(3,7)+
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*
     '            (XG(2,6)*XG(3,7)+XG(2,4)*XG(3,9))+
     '        (-DCOSH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(1,2)-
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DSIN(XG(3,1))*XG(2,2)-
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,2))*
     '            XG(3,4)*XG(3,7)-
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*
     '            (XG(3,6)*XG(3,7)+XG(3,4)*XG(3,9))+
     '        ( DCOSH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(1,2)+
     '          DSINH(XG(1,1))*DCOS(XG(2,1))*DCOS(XG(3,1))*XG(2,2)-
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DSIN(XG(3,1))*XG(3,2))*
     '            XG(3,10)+
     '          DSINH(XG(1,1))*DSIN(XG(2,1))*DCOS(XG(3,1))*XG(3,11))
          GO TO 99

    5   GO TO (51,51,51,51,51,51),nu  !oblate spheroid
   51     CALL XZ(ICOORD,XG(1,nu),ZG_LOCAL(1,nu))       !value

 99     CONTINUE
C     endcase
      DO nj=1,NJT
        XGRC(nj,nu)=ZG_LOCAL(nj,nu)
      ENDDO

C     CALL EXITS('XZ_DERIV')
      RETURN
      END


