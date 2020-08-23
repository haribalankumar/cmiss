      REAL*8 FUNCTION PLXI(nj,NPL,nr,nu,DL,XI,XP)

C#### Function: PLXI
C###  Type: REAL*8
C###  Description:
C###    <HTML>
C###    PLXI interpolates nodal values XP(nk,nv,nj,np) in a line
C###    segment defined by NPL.
C###    <PRE>
C###    If  NPL(1,nj)=1 :linear Lagrange interpolation
C###    "      "      2 :quad       "          "
C###    "      "      3 :cubic      "          "
C###    "      "      4 :cubic  Hermite        "
C###    "      "      5 :Bspline               "
C###    "      "      6 :quadratic Hermite (node 1)
C###    "      "      7 :quadratic Hermite (node 2)
C###    </PRE> </HTML>

C**** Added time interpolation 1/11/90 AAY


      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER nj,NPL(5,0:3),nr,nu
      REAL*8 DL(3),XI,XP(NKM,NVM,NJM,NPM)
!     Local Variables
      INTEGER ni,nn,NPL2,nv
      REAL*8 PH2,PH3,PL1,PL2,PL3,XN_LOCAL(4),ZERO_TOLERANCE
      CHARACTER ERROR*10

      PARAMETER(ZERO_TOLERANCE=1.0d-6)

      nv=1 ! temporary cpb 22/11/94

      ni=NPL(1,0)
      PLXI=0.0d0
C PJH 16DEC95 Update for drawing lines from nodal field variables
C CS 9OCT97 Not necessary with new NPL structure
C      IF(nj.LE.NJT) THEN !nj is geometric variable
C        NPL2=NPL(1+nj) !basis type for line - see above
C      ELSE !nj is field variable & need to use corresponding geom basis
C        NPL2=NPL(1-NJT,nj) !basis type for line - see above
C      ENDIF
      NPL2=NPL(1,nj) !basis type for line - see above

      IF((NPL2.LE.4).OR.(NPL2.EQ.6).OR.(NPL2.EQ.7)) THEN
        XN_LOCAL(1)=XP(1,nv,nj,NPL(2,1))
        XN_LOCAL(2)=XP(1,nv,nj,NPL(3,1))
        IF(NPL2.EQ.2) THEN      !quadratic Lagrange
          XN_LOCAL(3)=XP(1,nv,nj,NPL(4,1))
        ELSE IF(NPL2.EQ.3) THEN !cubic Lagrange
          XN_LOCAL(3)=XP(1,nv,nj,NPL(4,1))
          XN_LOCAL(4)=XP(1,nv,nj,NPL(5,1))
        ELSE IF(NPL2.EQ.4) THEN !cubic Hermite
          XN_LOCAL(3)=XP(NPL(4,1),nv,nj,NPL(2,1)) !deriv nk number at 1st node
          XN_LOCAL(4)=XP(NPL(5,1),nv,nj,NPL(3,1)) !deriv nk number at 2nd node
        ELSE IF(NPL2.EQ.6) THEN !quadratic Hermite node 1
          XN_LOCAL(3)=XP(NPL(5,1),nv,nj,NPL(3,1)) !deriv nk number at 2nd node
        ELSE IF(NPL2.EQ.7) THEN !quadratic Hermite node 2
          XN_LOCAL(3)=XP(NPL(4,1),nv,nj,NPL(2,1)) !deriv nk number at 1st node
        ENDIF
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' PLXI: nj='',I1,'' XI='',D12.4)') nj,XI
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' XN_LOCAL(1)='',E12.3,'' XN_LOCAL(2)='','
     '    //'D12.4)')
     '    XN_LOCAL(1),XN_LOCAL(2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(ITYP10(nr).eq.2.and.nj.EQ.2) THEN
        IF((ni.EQ.1).AND.(XN_LOCAL(1).GE.XN_LOCAL(2)))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(2)-XN_LOCAL(1)).GT.PI))
     '    XN_LOCAL(1)=XN_LOCAL(1)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(1)-XN_LOCAL(2)).GT.PI))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
      ELSE IF(ITYP10(nr).eq.3.and.nj.EQ.2) THEN
        IF((ni.EQ.1).AND.(XN_LOCAL(1).GE.XN_LOCAL(2)))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(2)-XN_LOCAL(1)).GT.PI))
     '    XN_LOCAL(1)=XN_LOCAL(1)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(1)-XN_LOCAL(2)).GT.PI))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
      ELSE IF(ITYP10(nr).eq.3.and.nj.EQ.3) THEN
        IF((ni.ne.1).AND.((XN_LOCAL(2)-XN_LOCAL(1)).GT.PI))
     '    XN_LOCAL(1)=XN_LOCAL(1)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(1)-XN_LOCAL(2)).GT.PI))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
      ELSE IF(ITYP10(nr).eq.4.and.nj.EQ.3) THEN
        IF((ni.EQ.1).AND.(XN_LOCAL(2).GE.XN_LOCAL(1)))
     '    XN_LOCAL(1)=XN_LOCAL(1)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(2)-XN_LOCAL(1)).GT.PI))
     '    XN_LOCAL(1)=XN_LOCAL(1)+2.0d0*PI
        IF((ni.ne.1).AND.((XN_LOCAL(1)-XN_LOCAL(2)).GT.PI))
     '    XN_LOCAL(2)=XN_LOCAL(2)+2.0d0*PI
      ENDIF
      IF(ITYP10(nr).eq.2.and.nj.eq.2.and.ni.EQ.2) THEN
        IF(DABS(XP(1,nv,1,NPL(2,1))).LT.ZERO_TOLERANCE) THEN
          XN_LOCAL(1)=XN_LOCAL(2)
        ELSE IF(DABS(XP(1,nv,1,NPL(3,1))).LT.ZERO_TOLERANCE) THEN
          XN_LOCAL(2)=XN_LOCAL(1)
        ENDIF
      ELSE IF(ITYP10(nr).eq.4.and.nj.eq.3.and.ni.EQ.2) THEN
        IF(DABS(XP(1,nv,2,NPL(2,1))).LT.ZERO_TOLERANCE) THEN
          XN_LOCAL(1)=XN_LOCAL(2)
        ELSE IF(DABS(XP(1,nv,2,NPL(3,1))).LT.ZERO_TOLERANCE) THEN
          XN_LOCAL(2)=XN_LOCAL(1)
        ENDIF
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,'('' XN_LOCAL(1)='',E12.3,'' XN_LOCAL(2)='','
     '    //'D12.4)')
     '    XN_LOCAL(1),XN_LOCAL(2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(NPL2.EQ.1) THEN      !Linear Lagrange
        DO nn=1,2
          PLXI=PLXI+PL1(nn,nu,XI)*XN_LOCAL(nn)
        ENDDO
      ELSE IF(NPL2.EQ.2) THEN !Quadratic Lagrange
        DO nn=1,3
          PLXI=PLXI+PL2(nn,nu,XI)*XN_LOCAL(nn)
        ENDDO
      ELSE IF(NPL2.EQ.3) THEN !Cubic Lagrange
        DO nn=1,4
          PLXI=PLXI+PL3(nn,nu,XI)*XN_LOCAL(nn)
        ENDDO
      ELSE IF(NPL2.EQ.4) THEN !Cubic Hermite
        DO nn=1,2
          PLXI=PLXI+PH3(nn,1,nu,XI)*XN_LOCAL(nn)
     '             +PH3(nn,2,nu,XI)*XN_LOCAL(nn+2)*DL(nn)
        ENDDO
      ELSE IF(NPL2.EQ.5) THEN !B-spline
C CS 9OCT97 B-spline is no longer supported ?
C        XI2=BKNOT(NPL(12),1)+XI*(BKNOT(NPL(13),1)-BKNOT(NPL(12),1))
C        IF(DABS(XI2-BKNOT(NPL(13),1)).LT.1.0D-06)XI2=XI2-1.0D-04
C        DO nn=1,MBSPL(1)+1
C          PLXI=PLXI+PSI3(nn,1,1,1,1,XI2)*XP(1,nv,nj,NPL(5+nn))
C        ENDDO
      ELSE IF(NPL2.EQ.6) THEN !quadratic Hermite (node 1)
        PLXI=PLXI+PH2(1,1,nu,1,XI)*XN_LOCAL(1)+
     '    PH2(2,1,nu,1,XI)*XN_LOCAL(2)+
     '    PH2(2,2,nu,1,XI)*XN_LOCAL(3)*DL(2)
      ELSE IF(NPL2.EQ.7) THEN !quadratic Hermite (node 2)
        PLXI=PLXI+PH2(1,1,nu,2,XI)*XN_LOCAL(1)+
     '    PH2(2,1,nu,2,XI)*XN_LOCAL(2)+
     '    PH2(1,2,nu,2,XI)*XN_LOCAL(3)*DL(1)
      ENDIF

 9999 RETURN
      END


