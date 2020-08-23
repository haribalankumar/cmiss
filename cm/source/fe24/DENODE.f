      SUBROUTINE DENODE(INDEX,INSTAT,ISEG,ISNONO,iw,NHP,NKH,
     '  np,NPNODE,NYNP,XP,BOTH,DEFORM,FIX,FRAME,CSEG,TYPE,ERROR,*)

C#### Subroutine: DENODE
C###  Description:
C###    DENODE creates node np identified by circle, unless node
C###    already exists.  (If existing node is the one previously
C###    created, NPOLD is set to zero and a new sequence of line
C###    segments begins).  Node creation is terminated when 2nd mouse
C###    button is pressed.
C###    TYPE is 'new' for  creating  new node np, 'old' for recreating
C###    old node np.
C###    Xwc,Ywc are rectangular cartesian world coordinates.
C###    XC(nj) are curvilinear coordinates.
C###    If NJT=3 and BOTH=.true.  both workstation viewports 1 and 2 are
C###    used for defining node position. If NJT=3 and BOTH=.false.
C###    workstation viewport IW is used only.
C###    DEFORM is .true. if new position is a deformed node position
C###    Note: Before requesting locator an inquiry is made to find the
C###    current normalisation transformation number and the required
C###    number is placed above it.

      IMPLICIT NONE
      INCLUDE 'alig00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER INDEX,INSTAT,ISEG(*),ISNONO(NWM,NPM),iw,
     '  NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),np,NPNODE(0:NP_R_M,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER CSEG(*)*(*),ERROR*(*),TYPE*3
      LOGICAL BOTH,DEFORM,FIX(NYM,NIYFIXM,NXM),FRAME
!     Local Variables
      INTEGER i,INDEX_OLD,ISNL,MP,nj,nr,NTDX,nv,nx
      REAL*8 PTS(3,50),THETA,X(3),XC(3),XWC,XWC1,XWC2,YWC,YWC1,YWC2,
     '  Z(3)
      LOGICAL EXISTN

      CALL ENTERS('DENODE',*9999)
      nv=1 ! Temporary MPN 12-Nov-94
      nr=1 !Needs fixing
      nx=1 !Needs fixing
      IF(np.EQ.0) THEN
        XWC1=0.0d0
        YWC1=0.0d0
        XWC2=0.0d0
        YWC2=0.0d0
      ELSE IF(np.GT.0) THEN
        X(1)=XP(1,nv,1,np)
        X(2)=XP(1,nv,2,np)
        X(3)=XP(1,nv,3,np)
        CALL XZ(ITYP10(nr),X,Z)
        IF(NJT.EQ.2) THEN
          XWC1=Z(1)
          YWC1=Z(2)
        ELSE IF(NJT.EQ.3) THEN
          IF(FRAME) THEN
C 30May96 MPN unused?
            CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
c            CALL FBIMGPT(5,REAL(Z(1)),REAL(Z(2)),REAL(Z(3)),
c     '        REAL(XWC1),REAL(YWC1))
c            CALL FBIMGPT(6,REAL(Z(1)),REAL(Z(2)),REAL(Z(3)),
c     '        REAL(XWC2),REAL(YWC2))
          ELSE
            XWC1=Z(1)
            YWC1=Z(3)
            XWC2=Z(2)
            YWC2=Z(3)
          ENDIF
        ENDIF
      ENDIF

      IF(NJT.EQ.2.OR..NOT.BOTH) THEN
        IF(iw.EQ.5.OR.iw.EQ.6) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,*)' Calling ACWK for wkst 1'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL ACWK(1,0,ERROR,*9999)
          PROJEC='NONE' !don't use image transformations in plotting node
        ENDIF
        CALL ACWK(iw,0,ERROR,*9999)
        IF(TYPE(1:3).EQ.'NEW') THEN
          IF(np.EQ.0) THEN
            WRITE(OP_STRING,'('' >>Locate node position on '',I1)')iw
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(np.GT.0) THEN
            WRITE(OP_STRING,'('' >>Locate node '',I4,'' on '',I1)')
     '        np,iw
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(TYPE(1:3).EQ.'OLD') THEN
          WRITE(OP_STRING,'('' >>Relocate node '',I4,'' on '',I1)')
     '      np,iw
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(DEFORM) THEN
          CALL LOCATOR(INSTAT,0.0d0,XWC,0.0d0,YWC,
     '      ERROR,*9999)
        ELSE IF(.NOT.DEFORM) THEN
          IF(iw.LT.5) THEN
            CALL LOCATOR(INSTAT,0.0d0,XWC,0.0d0,YWC,
     '        ERROR,*9999)
          ELSE
C 30May96 MPN unused?
            CALL ASSERT(.FALSE.,'>>ERROR: Not implemented',ERROR,*9999)
C            CALL FBPOINT(iw,INSTAT,XWC,YWC,ERROR,*9999)
          ENDIF
        ENDIF
        IF(iw.EQ.1) THEN
          IF(ALIGNMENT_ON) THEN
            XWC=DBLE(XMIN+REAL(DX_ALIG*NINT((XWC-DBLE(XMIN))/DX_ALIG)))
            YWC=DBLE(YMIN+REAL(DY_ALIG*NINT((YWC-DBLE(YMIN))/DY_ALIG)))
          ENDIF
          Z(1)=XWC
          Z(2)=YWC
          Z(3)=0.0D0
          CALL ZX(ITYP10(nr),Z,XC)
        ELSE IF(iw.EQ.4) THEN
          CALL PSCOORD1(XC,XWC,YWC,ERROR,*9999)
        ELSE IF(iw.GE.5) THEN
          Z(1)=XWC
          Z(2)=YWC
          Z(3)=0.0D0
          CALL ZX(ITYP10(nr),Z,XC)
        ENDIF
        IF(INSTAT.EQ.1) THEN
!PJH 13-feb-92 IF(TYPE(1:3).EQ.'NEW'.AND.EXISTN(MP,XP,XC)) THEN
!           write(*,'('' node exists'')')
!           IF(mp.NE.NPOLD) THEN
!             np=MP
!           ELSE IF(MP.EQ.NPOLD) THEN
!             np=0
!             NPOLD=0
!           ENDIF
!         ELSE
            IF(TYPE(1:3).EQ.'NEW') THEN
              NPT(1)=np
c             np=NPT(1)
            ENDIF
C MAYBE WRITE TO IODI, NOT SURE
            WRITE(OP_STRING,
     '        '('' Node= '',I4,'' Xwc= '',E10.3,'' Ywc= '',E10.3,'
     '        //''' XC: '',2E11.3)') np,XWC,YWC,XC(1),XC(2)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            IF(iw.EQ.4) THEN
              DO nj=2,NJT
                XP(1,nv,nj,np)=XC(nj)
              ENDDO
            ELSE
              DO nj=1,NJT
                XP(1,nv,nj,np)=XC(nj)
              ENDDO
            ENDIF
            IF(NJT.EQ.2) THEN
              CALL SGNODE(INDEX,ISEG,ISNONO(iw,np),iw,
     '          NHP(1,nr,nx),NKH(1,1,1,nr),np,nr,nx,NYNP,
     '          FIX(1,1,nx),XWC,YWC,0.0d0,CSEG,ERROR,*9999)
              IF(iw.EQ.5.OR.iw.EQ.6) THEN
                IF(DOP) THEN
                  WRITE(OP_STRING,*)' Calling sgnode for wkst 1'
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                CALL SGNODE(INDEX,ISEG,ISNONO(1,np),1,
     '            NHP(1,nr,nx),NKH(1,1,1,nr),np,nr,nx,NYNP,
     '            FIX(1,1,nx),XWC,YWC,0.0d0,CSEG,ERROR,*9999)
              ENDIF
            ENDIF
!         ENDIF
        ENDIF
        CALL DAWK(iw,0,ERROR,*9999)
        IF(iw.EQ.5.OR.iw.EQ.6) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,*)' Calling DAWK for wkst 1'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL DAWK(1,0,ERROR,*9999)
        ENDIF
      ELSE IF(NJT.EQ.3.AND.BOTH) THEN
        IF(BOTH.AND.(.NOT.FRAME)) iw=1
        IF(iw.EQ.1.OR.iw.GE.5) THEN
          CALL ACWK(iw,0,ERROR,*9999)
          WRITE(OP_STRING,'('' >>Locate node position on '',I1)') iw
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(iw.EQ.1) THEN
            CALL LOCATOR(INSTAT,0.0D0,XWC1,0.0D0,
     '        YWC1,ERROR,*9999)
          ELSE IF(iw.GE.5) THEN
C 30May96 MPN unused?
            CALL ASSERT(.FALSE.,'>>ERROR: Not implemented',ERROR,*9999)
C            CALL FBPOINT(iw,INSTAT,REAL(XWC1),REAL(YWC1),ERROR,*9999)
          ENDIF

          IF(INSTAT.EQ.1) THEN
            IF(TYPE(1:3).EQ.'NEW') THEN
              NPT(1)=NPT(1)+1
              np=NPT(1)
            ENDIF
            WRITE(OP_STRING,'('' Node='',I4,'' Xwc1= '',D10.3,'
     '        //''' Ywc1= '',D10.3)') np,XWC1,YWC1
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            IF(iw.EQ.1) THEN
              CALL SGNODE(INDEX,ISEG,ISNONO(1,np),1,NHP(1,nr,nx),
     '          NKH(1,1,1,nr),np,nr,nx,NYNP,FIX(1,1,nx),XWC1,0.0d0,
     '          YWC1,CSEG,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,0,ERROR,*9999)
          ELSE
            CALL DAWK(iw,0,ERROR,*9999)
            GO TO 9998
          ENDIF
        ENDIF

        IF(BOTH.AND.(.NOT.FRAME)) iw=2
        IF(iw.EQ.2.OR.iw.GE.5) THEN
          CALL ACWK(iw,0,ERROR,*9999)
          WRITE(OP_STRING,'('' >>Locate node position on '',I1)') iw
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(BOTH) THEN
C ***       Draw guide line in second viewport
            IF(ITYP10(nr).EQ.1) THEN
              NTDX=2
              PTS(1,1)=DBLE(YMIN)
              PTS(1,2)=DBLE(YMAX)
              PTS(2,1)=YWC1
              PTS(2,2)=YWC1
            ELSE IF(ITYP10(nr).GT.1) THEN
              NTDX=50
              DO i=1,50
                THETA=2.0D0*PI*(i-1)/49.0D0
                PTS(1,i)=XWC1*DCOS(THETA)
                PTS(2,i)=XWC1*DSIN(THETA)
              ENDDO
            ENDIF
            IF(iw.EQ.2) THEN
C??? LKC 12-APR-98 ISNL is used before set ???
              ISNL=0
              CALL OPEN_SEGMENT(ISNL,ISEG,iw,'NODE LINE',INDEX,
     '          INDEX_OLD,1,1,CSEG,ERROR,*9999)
              CALL POLYLINE(1,iw,NTDX,PTS,ERROR,*9999)
              CALL CLOSE_SEGMENT(ISNL,iw,ERROR,*9999)
            ELSE IF(iw.GE.5) THEN
C 30May96 MPN unused?
            CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
C              CALL FBAUXPT(iw,INSTAT,REAL(XWC1),REAL(YWC1),REAL(XWC2),
C     '          REAL(YWC2),ERROR,*9999)
            ENDIF
          ENDIF
          IF(iw.EQ.2) THEN
            CALL LOCATOR(INSTAT,0.0D0,XWC2,0.0D0,
     '        YWC2,ERROR,*9999)
            CALL DELETE_SEGMENT(ISNL,ISEG,iw,ERROR,*9999)
          ENDIF
          IF(INSTAT.EQ.1) THEN
            WRITE(OP_STRING,'('' Node='',I4,'' Xwc2= '',D10.3,'
     '        //''' Ywc2= '',D10.3)') np,XWC2,YWC2
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            IF(ITYP10(nr).EQ.1) THEN
              Z(1)=XWC1
              Z(2)=XWC2
              Z(3)=YWC1
            ELSE IF(ITYP10(nr).GT.1) THEN
              Z(1)=XWC2
              Z(2)=YWC2
              Z(3)=YWC1
            ENDIF
            CALL ZX(ITYP10(nr),Z,XC)
            IF(EXISTN(MP,NPNODE,XP,XC)) THEN
              WRITE(OP_STRING,'('' node exists'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL DELETE_SEGMENT(ISNONO(1,np),ISEG,iw,ERROR,*9999)
              ISNONO(1,np)=0
C              NTSG=NTSG-1 !used before set.
              NPT(1)=NPT(1)-1
              np=MP
            ELSE IF(iw.EQ.2) THEN
              CALL SGNODE(INDEX,ISEG,ISNONO(2,np),2,NHP(1,nr,nx),
     '          NKH(1,1,1,nr),np,nr,nx,NYNP,FIX(1,1,nx),0.0d0,XWC2,YWC1,
     '          CSEG,ERROR,*9999)
            ENDIF
          ENDIF
          CALL DAWK(iw,0,ERROR,*9999)
        ENDIF
        IF(FRAME) THEN
C 30May96 MPN unused?
          CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
c          CALL FBWLDPT(iw,REAL(XWC1),REAL(YWC1),REAL(XWC2),REAL(YWC2),
c     '      REAL(Z1),REAL(Z2),REAL(Z3))
C          Z(1)=Z1
C          Z(2)=Z2
C          Z(3)=Z3
        ELSE
          Z(1)=XWC1
          Z(2)=XWC2
          IF(iw.EQ.1) THEN
            Z(3)=YWC1
          ELSE IF(iw.EQ.2) THEN
            Z(3)=YWC2
          ENDIF
        ENDIF
C GMH 8/1/97 Update cmgui link
        CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
        CALL ZX(ITYP10(nr),Z,XC)
        XP(1,nv,1,np)=XC(1)
        XP(1,nv,2,np)=XC(2)
        XP(1,nv,3,np)=XC(3)
        CALL ACWK(3,0,ERROR,*9999)
        CALL SGNODE(INDEX,ISEG,ISNONO(3,np),3,
     '    NHP(1,nr,nx),NKH(1,1,1,nr),np,nr,nx,NYNP,
     '    FIX(1,1,nx),Z(1),Z(2),Z(3),CSEG,ERROR,*9999)
        CALL DAWK(3,0,ERROR,*9999)
      ENDIF

 9998 CALL EXITS('DENODE')
      RETURN
 9999 CALL ERRORS('DENODE',ERROR)
      CALL EXITS('DENODE')
      RETURN 1
      END


