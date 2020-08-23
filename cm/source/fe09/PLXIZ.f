      REAL*8 FUNCTION PLXIZ(nh,nhx,NPL,nr,nu,DL,XI,ZP)

C#### Function: PLXIZ
C###  Type: REAL*8
C###  Description:
C###    <HTML>
C###    PLXIZ interpolates nodal values ZP(nk,nv,nh,np,nc) in a line
C###    segment definedby NPL.
C###    <PRE>
C###    If  NPL(1,nhx)=1 :linear Lagrange interpolation
C###    "      "       2 :quad       "          "
C###    "      "       3 :cubic      "          "
C###    "      "       4 :cubic  Hermite        "
C###    "      "       5 :Bspline               "
C###    "      "       6 :quadratic Hermite (node 1)
C###    "      "       7 :quadratic Hermite (node 3)
C###    </PRE> </HTML>

C**** Added time interpolation 1/11/90 AAY

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER nh,nhx,NPL(5,0:3),nr,nu
      REAL*8 DL(3),XI,ZP(NKM,NVM,NHM,NPM,NCM)
!     Local Variables
      INTEGER nc,ni,nn,NPL2,nv
      REAL*8 PH2,PH3,PL1,PL2,PL3,XN_LOCAL(4),ZERO_TOLERANCE
      CHARACTER ERROR*10

      PARAMETER(ZERO_TOLERANCE=1.0d-6)

      nc=1 !temporary cpb 23/11/94
      nv=1 !temporary cpb 23/11/94

      ni=NPL(1,0)
      PLXIZ=0.0d0
      NPL2=NPL(1,nhx) !basis type for line - see above

      IF((NPL2.LE.4).OR.(NPL2.EQ.6).OR.(NPL2.EQ.7)) THEN
        XN_LOCAL(1)=ZP(1,nv,nh,NPL(2,1),nc)
        XN_LOCAL(2)=ZP(1,nv,nh,NPL(3,1),nc)
        IF(NPL2.EQ.2) THEN      !quadratic Lagrange
          XN_LOCAL(3)=ZP(1,nv,nh,NPL(4,1),nc)
        ELSE IF(NPL2.EQ.3) THEN !cubic Lagrange
          XN_LOCAL(3)=ZP(1,nv,nh,NPL(4,1),nc)
          XN_LOCAL(4)=ZP(1,nv,nh,NPL(5,1),nc)
        ELSE IF(NPL2.EQ.4) THEN !cubic Hermite
          XN_LOCAL(3)=ZP(NPL(4,1),nv,nh,NPL(2,1),nc) !deriv nk number at 1st node
          XN_LOCAL(4)=ZP(NPL(5,1),nv,nh,NPL(3,1),nc) !deriv nk number at 2nd node
        ELSE IF(NPL2.EQ.6) THEN !quadratic Hermite (node 1)
          XN_LOCAL(3)=ZP(NPL(5,1),nv,nh,NPL(3,1),nc) !deriv nk number at 2nd node
        ELSE IF(NPL2.EQ.7) THEN !quadratic Hermite (node 3)
          XN_LOCAL(3)=ZP(NPL(4,1),nv,nh,NPL(2,1),nc) !deriv nk number at 1st node
        ENDIF
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' PLXIZ: nh='',I1,'' XI='',D12.4)') nh,XI
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' XN_LOCAL(1)='',E12.3,'' XN_LOCAL(2)='','
     '    //'D12.4)')
     '    XN_LOCAL(1),XN_LOCAL(2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(ITYP10(nr).eq.2.and.nhx.EQ.2) THEN
        IF((ni.EQ.1).AND.(XN_LOCAL(1).GE.XN_LOCAL(2)))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(2)-XN_LOCAL(1)).GT.PI))
     '    XN_LOCAL(1)=XN_LOCAL(1)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(1)-XN_LOCAL(2)).GT.PI))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
      ELSE IF(ITYP10(nr).eq.3.and.nhx.EQ.2) THEN
        IF((ni.EQ.1).AND.(XN_LOCAL(1).GE.XN_LOCAL(2)))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(2)-XN_LOCAL(1)).GT.PI))
     '    XN_LOCAL(1)=XN_LOCAL(1)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(1)-XN_LOCAL(2)).GT.PI))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
      ELSE IF(ITYP10(nr).eq.3.and.nhx.EQ.3) THEN
        IF((ni.ne.1).AND.((XN_LOCAL(2)-XN_LOCAL(1)).GT.PI))
     '    XN_LOCAL(1)=XN_LOCAL(1)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(1)-XN_LOCAL(2)).GT.PI))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
      ELSE IF(ITYP10(nr).eq.4.and.nhx.EQ.3) THEN
        IF((ni.EQ.1).AND.(XN_LOCAL(2).GE.XN_LOCAL(1)))
     '    XN_LOCAL(1)=XN_LOCAL(1)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(2)-XN_LOCAL(1)).GT.PI))
     '    XN_LOCAL(1)=XN_LOCAL(1)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(1)-XN_LOCAL(2)).GT.PI))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
      ENDIF
      IF(ITYP10(nr).eq.2.and.nhx.eq.2.and.ni.EQ.2) THEN
        IF(DABS(ZP(1,nv,1,NPL(2,1),nc)).LT.ZERO_TOLERANCE) THEN
          XN_LOCAL(1)=XN_LOCAL(2)
        ELSE IF(DABS(ZP(1,nv,1,NPL(3,1),nc)).LT.ZERO_TOLERANCE) THEN
          XN_LOCAL(2)=XN_LOCAL(1)
        ENDIF
      ELSE IF(ITYP10(nr).eq.4.and.nhx.eq.3.and.ni.EQ.2) THEN
        IF(DABS(ZP(1,nv,2,NPL(2,1),nc)).LT.ZERO_TOLERANCE) THEN
          XN_LOCAL(1)=XN_LOCAL(2)
        ELSE IF(DABS(ZP(1,nv,2,NPL(3,1),nc)).LT.ZERO_TOLERANCE) THEN
          XN_LOCAL(2)=XN_LOCAL(1)
        ENDIF
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,'('' XN_LOCAL(1)='',D12.4,'' XN_LOCAL(2)='','
     '    //'D12.4)')
     '    XN_LOCAL(1),XN_LOCAL(2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(NPL2.EQ.1) THEN      !Linear Lagrange
        DO nn=1,2
          PLXIZ=PLXIZ+PL1(nn,nu,XI)*XN_LOCAL(nn)
        ENDDO
      ELSE IF(NPL2.EQ.2) THEN !Quadratic Lagrange
        DO nn=1,3
          PLXIZ=PLXIZ+PL2(nn,nu,XI)*XN_LOCAL(nn)
        ENDDO
      ELSE IF(NPL2.EQ.3) THEN !Cubic Lagrange
        DO nn=1,4
          PLXIZ=PLXIZ+PL3(nn,nu,XI)*XN_LOCAL(nn)
        ENDDO
      ELSE IF(NPL2.EQ.4) THEN !Cubic Hermite
        DO nn=1,2
          PLXIZ=PLXIZ+PH3(nn,1,nu,XI)*XN_LOCAL(nn)
     '             +PH3(nn,2,nu,XI)*XN_LOCAL(nn+2)*DL(nn)
        ENDDO
      ELSE IF(NPL2.EQ.5) THEN !B-spline
C CS 9OCT97 B-spline no longer supported ?
C        XI2=BKNOT(NPL(12),1)+XI*(BKNOT(NPL(13),1)-BKNOT(NPL(12),1))
C        IF(DABS(XI2-BKNOT(NPL(13),1)).LT.1.0D-06) XI2=XI2-1.0d-04
C        DO nn=1,MBSPL(1)+1
C          PLXIZ=PLXIZ+PSI3(nn,1,1,1,1,XI2)*ZP(1,nv,nh,NPL(5+nn),nc)
C        ENDDO
      ELSE IF(NPL2.EQ.6) THEN !quadratic Hermite (node 1)
        PLXIZ=PLXIZ+PH2(1,1,nu,1,XI)*XN_LOCAL(1)+
     '     PH2(2,1,nu,1,XI)*XN_LOCAL(2)
     '    +PH2(2,2,nu,1,XI)*XN_LOCAL(3)*DL(2)
      ELSE IF(NPL2.EQ.7) THEN !quadratic Hermite (node 3)
        PLXIZ=PLXIZ+PH2(1,1,nu,2,XI)*XN_LOCAL(1)+
     '     PH2(2,1,nu,2,XI)*XN_LOCAL(2)
     '    +PH2(1,2,nu,2,XI)*XN_LOCAL(3)*DL(1)
      ENDIF

 9999 RETURN
      END


