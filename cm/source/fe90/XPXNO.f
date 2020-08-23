      SUBROUTINE XPXNO(NBH,NENP,NP_INTERFACE,np,nr,NW,nx,DSDX,DXDS,XG,
     '  XNO,XP,ERROR,*)

C#### Subroutine: XPXNO
C###  Description:
C###    XPXNO calculates the vector XNO - the direction(s) in
C###    which the normal BIE is differentiated.

C**** XNO(nj,i) is the njth component of the ith direction in
C**** which the equation is differentiated.
C**** Also returns DSDX(ni,nj)=dSi/dXj at the point np.
C**** This subroutine assumes that geometric derivative information
C**** is contained in the XP array.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NP_INTERFACE(0:NPM,0:3),np,nr,NW(NEM,3),nx
      REAL*8 DSDX(3,*),DXDS(3,*),XNO(3,*),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER M(5),NB2,ne,ni,NI2,NITB,nj,nv
      REAL*8 D,SUM1,VLENGTH1,XN_LOCAL(3)
      LOGICAL INTERFACE2

      DATA M/1,2,3,1,2/

      CALL ENTERS('XPXNO',*9999)

      nb2=NBH(NH_LOC(1,nx),1,NENP(np,1,nr))
      NITB=NIT(nb2)

      IF(NITB.EQ.1) THEN
        CALL ASSERT(NUT(nb2).GE.3.,'>>Increase NUM to at least 3',
     '    ERROR,*9999)
      ELSE
        CALL ASSERT(NUT(nb2).GE.4.,'>>Increase NUM to at least 4',
     '    ERROR,*9999)
      ENDIF
C***  Calculate dSi/dXj at the point np and the vector XNO (the
C***  direction(s) in which the normal BIE is differentiated). Need to
C***  firstly find dXj/dSi and then invert.
      nv=1 ! temporary
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        XNO(nj,1)=XP(2,nv,nj,np) !Derivative in the s1 direction
        DXDS(nj,1)=XP(2,nv,nj,np) !dXj/ds1
        IF(NITB.EQ.2) THEN !3D problem
          XNO(nj,2)=XP(3,nv,nj,np) !Derivative in the s2 direction
          DXDS(nj,2)=XP(3,nv,nj,np) !dXj/ds2 .dXj/dn is calculated below
! New AJP 13-1-94
          IF(ktyp92.gt.0.and.ktyp92.NE.2) THEN
! End new
            XNO(nj,3)=XP(4,nv,nj,np) !Derivative in the s1s2 direction
          ENDIF
        ENDIF
      ENDDO
      NI2=2*NITB-1
! New AJP 13-1-94
      IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3.AND.
     '  (KTYP92.EQ.2.OR.KTYP92.EQ.0)) NI2=NI2-1
        !no cross derivative
! End new
      DO ni=1,NI2
        SUM1=0.0d0
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          SUM1=SUM1+XNO(nj,ni)*XNO(nj,ni)
        ENDDO
        VLENGTH1=DSQRT(SUM1)
        IF(VLENGTH1.LE.ZERO_TOL) THEN
          ERROR='>>Zero vector length'
          GOTO 9999
        ENDIF
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          XNO(nj,ni)=XNO(nj,ni)/VLENGTH1
        ENDDO
      ENDDO
      INTERFACE2=.FALSE. !INTERFACE2 is true if node np is on the
C                         interface between regions in a coupled
C                         problem and the region to which
C                         NE2 belongs is NOT the region with
C                         the smallest region number.
      IF((NP_INTERFACE(np,0).GT.1).AND.(NP_INTERFACE(np,1).NE.nr))
     '  INTERFACE2=.TRUE.
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        XG(nj,2)=XNO(nj,1)
        IF(NITB.EQ.2) XG(nj,4)=XNO(nj,2)
      ENDDO
      ne=NENP(np,1,nr)
      CALL NORMAL(ne,nr,NW,XG,XN_LOCAL,INTERFACE2,ERROR,*9999)
      !Find normal vector at np
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        DXDS(nj,NJ_LOC(NJL_GEOM,0,nr))=XN_LOCAL(nj) !dXj/dn
      ENDDO
      IF(KTYP92.EQ.2) THEN !Differentiate in the normal direction
                          !instead of the tangential (2d) or
                          !cross (3d) direction.
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          XNO(nj,2*NITB-1)=XN_LOCAL(nj)
        ENDDO
      ELSE IF(KTYP92.EQ.3) THEN !Use the normal derivative BIE
                              !in place of the conventional BIE.
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          XNO(nj,2*NITB)=XN_LOCAL(nj)
        ENDDO
      ENDIF !End of KTYP92 loop
      !Calculate DSDX(ni,nj)=dSi/dXj
      IF(NITB.EQ.1) THEN !Inverse of 2x2 matrix
        D=DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1)
        CALL ASSERT(DABS(D).GE.RDELTA,'>>Zero determinant?',
     '    ERROR,*9999)
        DSDX(1,1)=DXDS(2,2)/D
        DSDX(1,2)=-DXDS(1,2)/D
        DSDX(2,1)=-DXDS(2,1)/D
        DSDX(2,2)=DXDS(1,1)/D
      ELSE IF(NITB.EQ.2) THEN !Inverse of 3x3 matrix
        D=DXDS(1,1)*(DXDS(2,2)*DXDS(3,3)-DXDS(3,2)*DXDS(2,3))
     '    +DXDS(1,2)*(DXDS(2,3)*DXDS(3,1)-DXDS(3,3)*DXDS(2,1))
     '    +DXDS(1,3)*(DXDS(2,1)*DXDS(3,2)-DXDS(3,1)*DXDS(2,2))
        DO ni=1,3
          DO nj=1,3
            DSDX(ni,nj)=(DXDS(M(nj+1),M(ni+1))*DXDS(M(nj+2),
     '        M(ni+2))-DXDS(M(nj+2),M(ni+1))*DXDS(M(nj+1),
     '        M(ni+2)))/D
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('XPXNO')
      RETURN
 9999 CALL ERRORS('XPXNO',ERROR)
      CALL EXITS('XPXNO')
      RETURN 1
      END





