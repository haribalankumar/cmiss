      SUBROUTINE CRHERM(IDO,ISEG,ISLINE,ISLINO,ISL2BE,ISL3BE,
     '  ISN2BE,ISN3BE,iw,NBJ,NEL,NKJ,NLLIST,NPL,NVJL,
     '  DL,XB,XP,ZP,CSEG,ERROR,*)

C#### Subroutine: CRHERM
C###  Description:
C###    CRHERM creates Bezier control points on workstation iw.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fbgr00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),ISEG(*),ISL2BE(NLM),ISL3BE(NLM),
     '  ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),ISN2BE(NLM),ISN3BE(NLM),iw,
     '  NBJ(NJM,NEM),NEL(0:NELM,NLM),NKJ(NJM,NPM),
     '  NLLIST(0:NLM),NPL(5,0:3,NLM),NVJL(4,NJM,NLM)
      REAL*8 DL(3,NLM),XB(2,NJM,NLM),XP(NKM,NVM,NJM,NPM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER I,IBEG,IEND,IFROMC,INDEX,INDEX_CPT,INDEX_TAN,
     '  INDEX_OLD,INSTAT,IPICK,ISEGM,IV1,LIADJ,N2LI,nl,
     '  nline,NOLINE,nolist,np,NP1,NP4,nr,NWIND,nx,INDEX_POLYMARKER,
     '  INDEX_POLYLINE
      REAL X,XBADJ(2),Y,YBADJ(2)
      REAL*8 DXI(2,3),PT(3,2),XBEZ(4),XSLOPE,XWC,YBEZ(4),YSLOPE,YWC,
     '  ZSLOPE

      CALL ENTERS('CRHERM',*9999)
      CALL ACWK(iw,0,ERROR,*9999)
      DO nline=1,NTLINE
        CALL STRING_TRIM(CSEG(ISLINE(iw,nline)),IBEG,IEND)
        IF(ISLINE(iw,nline).GT.0)THEN
          NWIND=IFROMC(CSEG(ISLINE(iw,nline))(IEND:))
          IF(NWIND.EQ.iw) NOLINE=nline
        ENDIF
      ENDDO
      DO 200 nolist=1,NLLIST(0)
        nl=NLLIST(nolist)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,*) 'nl=',nl
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(NPL(1,1,nl).EQ.4) THEN
C ***     Draw polymarkers at Bezier control points and tangent lines
          NP1=NPL(2,1,nl)
          NP4=NPL(3,1,nl)
          IF(NJT.EQ.2) THEN
            XBEZ(1)=XP(1,1,1,NP1)
            YBEZ(1)=XP(1,1,2,NP1)
            XBEZ(4)=XP(1,1,1,NP4)
            YBEZ(4)=XP(1,1,2,NP4)
            XBEZ(2)=XBEZ(1)+XP(NPL(4,1,nl),1,1,NP1)*DL(3,nl)/3.0d0
            YBEZ(2)=YBEZ(1)+XP(NPL(4,1,nl),1,2,NP1)*DL(3,nl)/3.0d0
            XBEZ(3)=XBEZ(4)-XP(NPL(5,1,nl),1,1,NP4)*DL(3,nl)/3.0d0
            YBEZ(3)=YBEZ(4)-XP(NPL(5,1,nl),1,2,NP4)*DL(3,nl)/3.0d0
            XB(1,1,nl)=XBEZ(2)
            XB(1,2,nl)=YBEZ(2)
            XB(2,1,nl)=XBEZ(3)
            XB(2,2,nl)=YBEZ(3)
            IF(RHTRAN) IV1=5
            IF(LFTRAN) IV1=6
          ELSE IF(NJT.EQ.3) THEN
            IF(iw.EQ.1) THEN
              XBEZ(1)=XP(1,1,1,NP1)
              YBEZ(1)=XP(1,1,3,NP1)
              XBEZ(4)=XP(1,1,1,NP4)
              YBEZ(4)=XP(1,1,3,NP4)
              XBEZ(2)=XBEZ(1)+XP(NPL(4,1,nl),1,1,NP1)*DL(3,nl)/3.d0
              YBEZ(2)=YBEZ(1)+XP(NPL(4,1,nl),1,3,NP1)*DL(3,nl)/3.d0
              XBEZ(3)=XBEZ(4)-XP(NPL(5,1,nl),1,1,NP4)*DL(3,nl)/3.d0
              YBEZ(3)=YBEZ(4)-XP(NPL(5,1,nl),1,3,NP4)*DL(3,nl)/3.d0
              XB(1,1,nl)=XBEZ(2)
              XB(1,3,nl)=YBEZ(2)
              XB(2,1,nl)=XBEZ(3)
              XB(2,3,nl)=YBEZ(3)
            ELSE IF(iw.EQ.2) THEN
              XBEZ(1)=XP(1,1,2,NP1)
              YBEZ(1)=XP(1,1,3,NP1)
              XBEZ(4)=XP(1,1,2,NP4)
              YBEZ(4)=XP(1,1,3,NP4)
              XBEZ(2)=XBEZ(1)+XP(NPL(4,1,nl),1,2,NP1)*DL(3,nl)/3.d0
              YBEZ(2)=YBEZ(1)+XP(NPL(4,1,nl),1,3,NP1)*DL(3,nl)/3.d0
              XBEZ(3)=XBEZ(4)-XP(NPL(5,1,nl),1,2,NP4)*DL(3,nl)/3.d0
              YBEZ(3)=YBEZ(4)-XP(NPL(5,1,nl),1,3,NP4)*DL(3,nl)/3.d0
              XB(1,2,nl)=XBEZ(2)
              XB(1,3,nl)=YBEZ(2)
              XB(2,2,nl)=XBEZ(3)
              XB(2,3,nl)=YBEZ(3)
            ELSE IF(iw.EQ.3) THEN
              IF(RHTRAN) THEN
                IV1=5
C                IV2=6
              ELSE
                IV1=6
C                IV2=5
              ENDIF
C old MPN unused?
              CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
c              CALL FBIMGPT(IV1,REAL(XP(1,1,1,NP1)),
c     '          REAL(XP(1,1,2,NP1)),
c     '          REAL(XP(1,1,3,NP1)),XBEZ(1),YBEZ(1))
c              CALL FBIMGPT(IV1,REAL(XP(1,1,1,NP4)),
c     '          REAL(XP(1,1,2,NP4)),
c     '          REAL(XP(1,1,3,NP4)),XBEZ(4),YBEZ(4))
              XB(1,1,nl)=XP(1,1,1,NP1)+XP(NPL(4,1,nl),1,1,NP1)
     '          *DL(3,nl)/3.d0
              XB(1,2,nl)=XP(1,1,2,NP1)+XP(NPL(4,1,nl),1,2,NP1)
     '          *DL(3,nl)/3.d0
              XB(1,3,nl)=XP(1,1,3,NP1)+XP(NPL(4,1,nl),1,3,NP1)
     '          *DL(3,nl)/3.d0
              XB(2,1,nl)=XP(1,1,1,NP4)-XP(NPL(5,1,nl),1,1,NP4)
     '          *DL(3,nl)/3.d0
              XB(2,2,nl)=XP(1,1,2,NP4)-XP(NPL(5,1,nl),1,2,NP4)
     '          *DL(3,nl)/3.d0
              XB(2,3,nl)=XP(1,1,3,NP4)-XP(NPL(5,1,nl),1,3,NP4)
     '          *DL(3,nl)/3.d0
C old MPN unused?
              CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
c              CALL FBIMGPT(IV1,REAL(XB(1,1,nl)),REAL(XB(1,2,nl)),
c     '          REAL(XB(1,3,nl)),XBEZ(2),YBEZ(2))
c              CALL FBIMGPT(IV1,REAL(XB(2,1,nl)),REAL(XB(2,2,nl)),
c     '          REAL(XB(2,3,nl)),XBEZ(3),YBEZ(3))
            ENDIF
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' XBEZ: '',4E12.3)') (XBEZ(I),I=1,4)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' YBEZ: '',4E12.3)') (YBEZ(I),I=1,4)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

          INDEX_CPT=INDEX_POLYMARKER(0,'POINT','SIZE1','BLACK')
          INDEX_TAN=INDEX_POLYLINE(0,'DOTTED','WIDTH1','BLACK')
          CALL OPEN_SEGMENT(ISN2BE(nl),ISEG,iw,'CPT1',INDEX_CPT,
     '      INDEX_OLD,nl,1,CSEG,ERROR,*9999)
          PT(1,2)=XBEZ(2)
          PT(2,2)=YBEZ(2)
          CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
          CALL CLOSE_SEGMENT(ISN2BE(nl),iw,ERROR,*9999)
          CALL DETECT(iw,ISEG,ISN2BE(nl),'DETECTABLE',ERROR,*9999)

          CALL OPEN_SEGMENT(ISL2BE(nl),ISEG,iw,'TAN1',
     '      INDEX_TAN,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
          PT(1,1)=XBEZ(1)
          PT(2,1)=YBEZ(1)
          CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
          CALL CLOSE_SEGMENT(ISL2BE(nl),iw,ERROR,*9999)

          CALL OPEN_SEGMENT(ISN3BE(nl),ISEG,iw,'CPT2',
     '      INDEX_CPT,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
          PT(1,1)=XBEZ(3)
          PT(2,1)=YBEZ(3)
          CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
          CALL CLOSE_SEGMENT(ISN3BE(nl),iw,ERROR,*9999)
          CALL DETECT(iw,ISEG,ISN3BE(nl),'DETECTABLE',ERROR,*9999)

          CALL OPEN_SEGMENT(ISL3BE(nl),ISEG,iw,'TAN2',
     '      INDEX_TAN,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
          PT(1,2)=XBEZ(4)
          PT(2,2)=YBEZ(4)
          CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
          CALL CLOSE_SEGMENT(ISL3BE(nl),iw,ERROR,*9999)
        ENDIF
 200  CONTINUE

      INSTAT=1
      DO WHILE(INSTAT.EQ.1)
C ***   Pick control point segments
        WRITE(OP_STRING,'('' >>Pick control point on '',I1,'':'')') iw
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL PICK(iw,'REQUEST',INSTAT,ISEGM,IPICK,ERROR,*9999)
        IF(INSTAT.EQ.1) THEN
          IF(CSEG(ISEGM)(1:4).EQ.'CPT1') THEN
            nl=IFROMC(CSEG(ISEGM)(53:57))
            np=NPL(2,1,nl)
C ***       Locate new control point
            WRITE(OP_STRING,'('' >>Relocate 1st control point for '','
     '        //'''line '',I4,'' on '',I1,'':'')') nl,iw
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(INSTAT,0.d0,XWC,0.d0,YWC,
     '        ERROR,*9999)
            IF(NJT.EQ.2) THEN
              XB(1,1,nl)=XWC
              XB(1,2,nl)=YWC
              XBEZ(2)=XWC
              YBEZ(2)=YWC
              X=REAL(XP(1,1,1,np))
              Y=REAL(XP(1,1,2,np))
            ELSE IF(NJT.EQ.3) THEN
              IF(iw.EQ.1) THEN
                XB(1,1,nl)=XWC
                XB(1,3,nl)=YWC
                XBEZ(2)=XWC
                YBEZ(2)=YWC
                X=REAL(XP(1,1,1,np))
                Y=REAL(XP(1,1,3,np))
              ELSE IF(iw.EQ.2) THEN
                XB(1,2,nl)=XWC
                XBEZ(2)=XWC
                YBEZ(2)=XB(1,3,nl)
                X=REAL(XP(1,1,2,np))
                Y=REAL(XP(1,1,3,np))
              ELSE IF(iw.EQ.3) THEN
C old MPN unused?
                CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
C               modify control point according to epipolar constraint
C                CALL FBIMGPT(IV2,REAL(XB(1,1,nl)),REAL(XB(1,2,nl)),
C     '            REAL(XB(1,3,nl)),X,Y)
C                IF(AUX) THEN
C                  CALL EPIPOLAR(IV2,X,Y,XWC,YWC)
C                ELSE
C                  CALL EPIPOLAR(IV1,XWC,YWC,X,Y)
C                ENDIF
C                CALL FBWLDPT(IV1,XWC,YWC,X,Y,XB1,XB2,XB3)
C                XB(1,1,nl)=DBLE(XB1)
C                XB(1,2,nl)=DBLE(XB2)
C                XB(1,3,nl)=DBLE(XB3)
C                XBEZ(2)=XWC
C                YBEZ(2)=YWC
c                CALL FBIMGPT(IV1,REAL(XP(1,1,1,np)),
c     '            REAL(XP(1,1,2,np)),REAL(XP(1,1,3,np)),X,Y)
              ENDIF
            ENDIF
            XBEZ(1)=DBLE(X)
            YBEZ(1)=DBLE(Y)

            CALL OPEN_SEGMENT(ISN2BE(nl),ISEG,iw,'CPT1',
     '        INDEX_CPT,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,2)=XBEZ(2)
            PT(2,2)=YBEZ(2)
            CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISN2BE(nl),iw,ERROR,*9999)
            CALL DETECT(iw,ISEG,ISN2BE(nl),'DETECTABLE',ERROR,*9999)

            CALL OPEN_SEGMENT(ISL2BE(nl),ISEG,iw,'TAN1',
     '        INDEX_TAN,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,1)=XBEZ(1)
            PT(2,1)=YBEZ(1)
            CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISL2BE(nl),iw,ERROR,*9999)

            XSLOPE=3.d0*(XB(1,1,nl)-XP(1,1,1,np))/DL(3,nl)
            YSLOPE=3.d0*(XB(1,2,nl)-XP(1,1,2,np))/DL(3,nl)
            ZSLOPE=3.d0*(XB(1,3,nl)-XP(1,1,3,np))/DL(3,nl)
C ***       Find adjacent control point if colinear
            N2LI=LIADJ(1,nl,NPL,NVJL)
            IF(n2li.ne.0.AND.NKJ(1,np).GE.2) THEN
              XBEZ(1)=DBLE(X)
              YBEZ(1)=DBLE(Y)
              XSLOPE=XSLOPE*DL(3,N2LI)
              YSLOPE=YSLOPE*DL(3,N2LI)
              ZSLOPE=ZSLOPE*DL(3,N2LI)
              IF(NJT.EQ.2) THEN
                XB(2,1,N2LI)=XP(1,1,1,np)-XSLOPE/3.d0
                XB(2,2,N2LI)=XP(1,1,2,np)-YSLOPE/3.d0
                XBADJ(1)=REAL(XB(2,1,N2LI))
                YBADJ(1)=REAL(XB(2,2,N2LI))
              ELSE IF(NJT.EQ.3) THEN
                IF(iw.EQ.1) THEN
                  XB(2,1,N2LI)=XP(1,1,1,np)-XSLOPE/3.d0
                  XB(2,3,N2LI)=XP(1,1,3,np)-ZSLOPE/3.d0
                  XBADJ(1)=REAL(XB(2,1,N2LI))
                  YBADJ(1)=REAL(XB(2,3,N2LI))
                ELSE IF(iw.EQ.2) THEN
                  XB(2,2,N2LI)=XP(1,1,2,np)-YSLOPE/3.d0
                  XB(2,3,N2LI)=XP(1,1,3,np)-ZSLOPE/3.d0
                  XBADJ(1)=REAL(XB(2,2,N2LI))
                  YBADJ(1)=REAL(XB(2,3,N2LI))
                ELSE IF(iw.EQ.3) THEN
                  XB(2,1,N2LI)=XP(1,1,1,np)-XSLOPE/3.d0
                  XB(2,2,N2LI)=XP(1,1,2,np)-YSLOPE/3.d0
                  XB(2,3,N2LI)=XP(1,1,3,np)-ZSLOPE/3.d0
C old MPN unused?
                  CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
c                  CALL FBIMGPT(IV1,REAL(XB(2,1,N2LI)),
c     '              REAL(XB(2,2,N2LI)),REAL(XB(2,3,N2LI)),
c     '              XBADJ(1),YBADJ(1))
                ENDIF
              ENDIF
              XBADJ(2)=REAL(XBEZ(1))
              YBADJ(2)=REAL(YBEZ(1))

              CALL OPEN_SEGMENT(ISN3BE(N2LI),ISEG,iw,'CPT2',
     '          INDEX_CPT,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,1)=DBLE(XBADJ(1))
              PT(2,1)=DBLE(YBADJ(1))
              CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISN3BE(N2LI),iw,ERROR,*9999)
              CALL DETECT(iw,ISEG,ISN3BE(N2LI),'DETECTABLE',ERROR,*9999)

              CALL OPEN_SEGMENT(ISL3BE(N2LI),ISEG,iw,'TAN2',
     '          INDEX_TAN,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,2)=DBLE(XBADJ(2))
              PT(2,2)=DBLE(YBADJ(2))
              CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISL3BE(N2LI),iw,ERROR,*9999)

            ENDIF

          ELSE IF(CSEG(ISEGM)(1:4).EQ.'CPT2') THEN
            nl=IFROMC(CSEG(ISEGM)(53:57))
            np=NPL(3,1,nl)
C ***       Locate new control point
            WRITE(OP_STRING,'('' >>Relocate 2nd control point for '','
     '        //'''line '',I4,'' on '',I1,'':'')') nl,iw
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(INSTAT,0.d0,XWC,0.d0,YWC,
     '        ERROR,*9999)
            IF(NJT.EQ.2) THEN
              XB(2,1,nl)=XWC
              XB(2,2,nl)=YWC
              XBEZ(3)=XWC
              YBEZ(3)=YWC
              X=REAL(XP(1,1,1,np))
              Y=REAL(XP(1,1,2,np))
            ELSE IF(NJT.EQ.3) THEN
              IF(iw.EQ.1) THEN
                XB(2,1,nl)=XWC
                XB(2,3,nl)=YWC
                XBEZ(3)=XWC
                YBEZ(3)=YWC
                X=REAL(XP(1,1,1,np))
                Y=REAL(XP(1,1,3,np))
              ELSE IF(iw.EQ.2) THEN
                XB(2,2,nl)=XWC
                XBEZ(3)=XWC
                YBEZ(3)=XB(2,3,nl)
                X=REAL(XP(1,1,2,np))
                Y=REAL(XP(1,1,3,np))
              ELSE IF(iw.EQ.3) THEN
C ***           modify control point according to epipolar constraint
C old MPN unused?
                CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
C                CALL FBIMGPT(IV2,REAL(XB(2,1,nl)),REAL(XB(2,2,nl)),
C     '            REAL(XB(2,3,nl)),X,Y)
C                IF(AUX) THEN
C                  CALL EPIPOLAR(IV2,X,Y,XWC,YWC)
C                ELSE
C                  CALL EPIPOLAR(IV1,XWC,YWC,X,Y)
C                ENDIF
C                CALL FBWLDPT(IV1,XWC,YWC,X,Y,XB(2,1,nl),XB(2,2,nl),
C     '            XB(2,3,nl))
c                XBEZ(3)=XWC
c                YBEZ(3)=YWC
c                CALL FBIMGPT(IV1,REAL(XP(1,1,1,np)),
c     '            REAL(XP(1,1,2,np)),
c     '            REAL(XP(1,1,3,np)),X,Y)
              ENDIF
            ENDIF
            XBEZ(4)=DBLE(X)
            YBEZ(4)=DBLE(Y)

            CALL OPEN_SEGMENT(ISN3BE(nl),ISEG,iw,'CPT2',
     '        INDEX_CPT,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,1)=XBEZ(3)
            PT(2,1)=YBEZ(3)
            CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISN3BE(nl),iw,ERROR,*9999)
            CALL DETECT(iw,ISEG,ISN3BE(nl),'DETECTABLE',ERROR,*9999)

            CALL OPEN_SEGMENT(ISL3BE(nl),ISEG,iw,'TAN2',
     '        INDEX_TAN,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,2)=XBEZ(4)
            PT(2,2)=YBEZ(4)
            CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISL3BE(nl),iw,ERROR,*9999)

            XSLOPE=3.d0*(XP(1,1,1,np)-XB(2,1,nl))/DL(3,nl)
            YSLOPE=3.d0*(XP(1,1,2,np)-XB(2,2,nl))/DL(3,nl)
            ZSLOPE=3.d0*(XP(1,1,3,np)-XB(2,3,nl))/DL(3,nl)
C ***       Find adjacent control point if colinear
            N2LI=LIADJ(2,nl,NPL,NVJL)
            IF(n2li.ne.0.AND.NKJ(1,np).GE.2) THEN
              XBEZ(1)=DBLE(X)
              YBEZ(1)=DBLE(Y)
              XSLOPE=XSLOPE*DL(3,N2LI)
              YSLOPE=YSLOPE*DL(3,N2LI)
              ZSLOPE=ZSLOPE*DL(3,N2LI)
              IF(NJT.EQ.2) THEN
                XB(1,1,N2LI)=XP(1,1,1,np)+XSLOPE/3.d0
                XB(1,2,N2LI)=XP(1,1,2,np)+YSLOPE/3.d0
                XBADJ(2)=REAL(XB(1,1,N2LI))
                YBADJ(2)=REAL(XB(1,2,N2LI))
              ELSE IF(NJT.EQ.3) THEN
                IF(iw.EQ.1) THEN
                  XB(1,1,N2LI)=XP(1,1,1,np)+XSLOPE/3.d0
                  XB(1,3,N2LI)=XP(1,1,3,np)+ZSLOPE/3.d0
                  XBADJ(2)=REAL(XB(1,1,N2LI))
                  YBADJ(2)=REAL(XB(1,3,N2LI))
                ELSE IF(iw.EQ.2) THEN
                  XB(1,2,N2LI)=XP(1,1,2,np)+YSLOPE/3.d0
                  XB(1,3,N2LI)=XP(1,1,3,np)+ZSLOPE/3.d0
                  XBADJ(2)=REAL(XB(1,2,N2LI))
                  YBADJ(2)=REAL(XB(1,3,N2LI))
                ELSE IF(iw.EQ.3) THEN
                  XB(1,1,N2LI)=XP(1,1,1,np)+XSLOPE/3.d0
                  XB(1,2,N2LI)=XP(1,1,2,np)+YSLOPE/3.d0
                  XB(1,3,N2LI)=XP(1,1,3,np)+ZSLOPE/3.d0
C old MPN unused?
                  CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
c                  CALL FBIMGPT(IV1,REAL(XB(1,1,N2LI)),
c     '              REAL(XB(1,2,N2LI)),REAL(XB(1,3,N2LI)),
c     '              XBADJ(2),YBADJ(2))
                ENDIF
              ENDIF
              XBADJ(1)=REAL(XBEZ(4))
              YBADJ(1)=REAL(YBEZ(4))

              CALL OPEN_SEGMENT(ISN2BE(N2LI),ISEG,iw,'CPT1',
     '          INDEX_CPT,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,2)=DBLE(XBADJ(2))
              PT(2,2)=DBLE(YBADJ(2))
              CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISN2BE(N2LI),iw,ERROR,*9999)
              CALL DETECT(iw,ISEG,ISN2BE(N2LI),'DETECTABLE',ERROR,*9999)

              CALL OPEN_SEGMENT(ISL2BE(N2LI),ISEG,iw,'TAN1',
     '          INDEX_TAN,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,1)=DBLE(XBADJ(1))
              PT(2,1)=DBLE(YBADJ(1))
              CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISL2BE(N2LI),iw,ERROR,*9999)
            ENDIF
          ENDIF

C ***     Calculate line length and redraw line segments
          IF(NPL(1,1,nl).EQ.4) THEN
            NP1=NPL(2,1,nl)
            NP4=NPL(3,1,nl)
            IF(NJT.EQ.2) THEN
              XBEZ(2)=DBLE(XB(1,1,nl))
              YBEZ(2)=DBLE(XB(1,2,nl))
              XBEZ(3)=DBLE(XB(2,1,nl))
              YBEZ(3)=DBLE(XB(2,2,nl))
              XBEZ(1)=XP(1,1,1,NP1)
              YBEZ(1)=XP(1,1,2,NP1)
              XBEZ(4)=XP(1,1,1,NP4)
              YBEZ(4)=XP(1,1,2,NP4)
              DXI(1,1)=3.d0*(XBEZ(2)-XBEZ(1))
              DXI(1,2)=3.d0*(YBEZ(2)-YBEZ(1))
              DXI(2,1)=3.d0*(XBEZ(4)-XBEZ(3))
              DXI(2,2)=3.d0*(YBEZ(4)-YBEZ(3))
              CALL DXIDL(IDO,NBJ,NEL(0,nl),NPL(1,0,nl),
     '          DL(1,nl),DXI,XP,ERROR,*9999)
              XP(NPL(4,1,nl),1,1,NP1)=DXI(1,1)/DL(3,nl)
              XP(NPL(4,1,nl),1,2,NP1)=DXI(1,2)/DL(3,nl)
              XP(NPL(5,1,nl),1,1,NP4)=DXI(2,1)/DL(3,nl)
              XP(NPL(5,1,nl),1,2,NP4)=DXI(2,2)/DL(3,nl)
            ELSE IF(NJT.EQ.3) THEN
              IF(iw.EQ.1) THEN
                XBEZ(1)=XP(1,1,1,NP1)
                YBEZ(1)=XP(1,1,3,NP1)
                XBEZ(4)=XP(1,1,1,NP4)
                YBEZ(4)=XP(1,1,3,NP4)
                XBEZ(2)=DBLE(XB(1,1,nl))
                YBEZ(2)=DBLE(XB(1,3,nl))
                XBEZ(3)=DBLE(XB(2,1,nl))
                YBEZ(3)=DBLE(XB(2,3,nl))
                DXI(1,1)=3.d0*(XBEZ(2)-XBEZ(1))
                DXI(1,3)=3.d0*(YBEZ(2)-YBEZ(1))
                DXI(2,1)=3.d0*(XBEZ(4)-XBEZ(3))
                DXI(2,3)=3.d0*(YBEZ(4)-YBEZ(3))
              ELSE IF(iw.EQ.2) THEN
                XBEZ(1)=XP(1,1,2,NP1)
                YBEZ(1)=XP(1,1,3,NP1)
                XBEZ(4)=XP(1,1,2,NP4)
                YBEZ(4)=XP(1,1,3,NP4)
                XBEZ(2)=XB(1,2,nl)
                YBEZ(2)=XB(1,3,nl)
                XBEZ(3)=XB(2,2,nl)
                YBEZ(3)=XB(2,3,nl)
                DXI(1,2)=3.d0*(XBEZ(2)-XBEZ(1))
                DXI(1,3)=3.d0*(YBEZ(2)-YBEZ(1))
                DXI(2,2)=3.d0*(XBEZ(4)-XBEZ(3))
                DXI(2,3)=3.d0*(YBEZ(4)-YBEZ(3))
              ELSE IF(iw.EQ.3) THEN
C old MPN unused?
                CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
c                CALL FBIMGPT(IV1,REAL(XP(1,1,1,NP1)),
c     '            REAL(XP(1,1,2,NP1)),
c     '            REAL(XP(1,1,3,NP1)),REAL(XBEZ(1)),REAL(YBEZ(1)))
c                CALL FBIMGPT(IV1,REAL(XP(1,1,1,NP4)),
c     '            REAL(XP(1,1,2,NP4)),
c     '            REAL(XP(1,1,3,NP4)),REAL(XBEZ(4)),REAL(YBEZ(4)))
c                CALL FBIMGPT(IV1,REAL(XB(1,1,nl)),REAL(XB(1,2,nl)),
c     '            REAL(XB(1,2,nl)),REAL(XBEZ(2)),REAL(YBEZ(2)))
c                CALL FBIMGPT(IV1,REAL(XB(2,1,nl)),REAL(XB(2,2,nl)),
c     '            REAL(XB(2,3,nl)),REAL(XBEZ(3)),REAL(YBEZ(3)))
                DXI(1,1)=3.d0*(XB(1,1,nl)-XP(1,1,1,NP1))
                DXI(1,2)=3.d0*(XB(1,2,nl)-XP(1,1,2,NP1))
                DXI(1,3)=3.d0*(XB(1,3,nl)-XP(1,1,3,NP1))
                DXI(2,1)=3.d0*(XP(1,1,1,NP4)-XB(2,1,nl))
                DXI(2,2)=3.d0*(XP(1,1,2,NP4)-XB(2,2,nl))
                DXI(2,3)=3.d0*(XP(1,1,3,NP4)-XB(2,3,nl))
              ENDIF
              CALL DXIDL(IDO,NBJ,NEL(0,nl),NPL(1,0,nl),
     '          DL(1,nl),DXI,XP,ERROR,*9999)
              IF(iw.EQ.1) THEN
                XP(NPL(4,1,nl),1,1,NP1)=DXI(1,1)/DL(3,nl)
                XP(NPL(4,1,nl),1,3,NP1)=DXI(1,3)/DL(3,nl)
                XP(NPL(5,1,nl),1,1,NP4)=DXI(2,1)/DL(3,nl)
                XP(NPL(5,1,nl),1,3,NP4)=DXI(2,3)/DL(3,nl)
              ELSE IF(iw.EQ.2) THEN
                XP(NPL(4,1,nl),1,2,NP1)=DXI(1,2)/DL(3,nl)
                XP(NPL(4,1,nl),1,3,NP1)=DXI(1,3)/DL(3,nl)
                XP(NPL(5,1,nl),1,2,NP4)=DXI(2,2)/DL(3,nl)
                XP(NPL(5,1,nl),1,3,NP4)=DXI(2,3)/DL(3,nl)
              ELSE IF(iw.EQ.3) THEN
                XP(NPL(4,1,nl),1,1,NP1)=DXI(1,1)/DL(3,nl)
                XP(NPL(4,1,nl),1,2,NP1)=DXI(1,2)/DL(3,nl)
                XP(NPL(4,1,nl),1,3,NP1)=DXI(1,3)/DL(3,nl)
                XP(NPL(5,1,nl),1,1,NP4)=DXI(2,1)/DL(3,nl)
                XP(NPL(5,1,nl),1,2,NP4)=DXI(2,2)/DL(3,nl)
                XP(NPL(5,1,nl),1,3,NP4)=DXI(2,3)/DL(3,nl)
              ENDIF
            ENDIF
          ENDIF
          IF(n2li.ne.0.AND.NKJ(1,np).EQ.2) THEN
            IF(NPL(1,1,N2LI).EQ.4) THEN
              NP1=NPL(2,1,N2LI)
              NP4=NPL(3,1,N2LI)
              IF(NJT.EQ.2) THEN
                XBEZ(1)=XP(1,1,1,NP1)
                YBEZ(1)=XP(1,1,2,NP1)
                XBEZ(4)=XP(1,1,1,NP4)
                YBEZ(4)=XP(1,1,2,NP4)
                XBEZ(2)=DBLE(XB(1,1,N2LI))
                YBEZ(2)=DBLE(XB(1,2,N2LI))
                XBEZ(3)=DBLE(XB(2,1,N2LI))
                YBEZ(3)=DBLE(XB(2,2,N2LI))
                DXI(1,1)=3.d0*(XBEZ(2)-XBEZ(1))
                DXI(1,2)=3.d0*(YBEZ(2)-YBEZ(1))
                DXI(2,1)=3.d0*(XBEZ(4)-XBEZ(3))
                DXI(2,2)=3.d0*(YBEZ(4)-YBEZ(3))
                CALL DXIDL(IDO,NBJ,NEL(0,nl),NPL(1,0,N2LI),
     '            DL(1,N2LI),DXI,XP,ERROR,*9999)
                XP(NPL(4,1,N2LI),1,1,NP1)=DXI(1,1)/DL(3,N2LI)
                XP(NPL(4,1,N2LI),1,2,NP1)=DXI(1,2)/DL(3,N2LI)
                XP(NPL(5,1,N2LI),1,1,NP4)=DXI(2,1)/DL(3,N2LI)
                XP(NPL(5,1,N2LI),1,2,NP4)=DXI(2,2)/DL(3,N2LI)
              ELSE IF(NJT.EQ.3) THEN
                IF(iw.EQ.1) THEN
                  XBEZ(1)=XP(1,1,1,NP1)
                  YBEZ(1)=XP(1,1,3,NP1)
                  XBEZ(4)=XP(1,1,1,NP4)
                  YBEZ(4)=XP(1,1,3,NP4)
                  XBEZ(2)=XB(1,1,N2LI)
                  YBEZ(2)=DBLE(XB(1,3,N2LI))
                  XBEZ(3)=DBLE(XB(2,1,N2LI))
                  YBEZ(3)=DBLE(XB(2,3,N2LI))
                  DXI(1,1)=3.d0*(XBEZ(2)-XBEZ(1))
                  DXI(1,3)=3.d0*(YBEZ(2)-YBEZ(1))
                  DXI(2,1)=3.d0*(XBEZ(4)-XBEZ(3))
                  DXI(2,3)=3.d0*(YBEZ(4)-YBEZ(3))
                ELSE IF(iw.EQ.2) THEN
                  XBEZ(1)=XP(1,1,2,NP1)
                  YBEZ(1)=XP(1,1,3,NP1)
                  XBEZ(4)=XP(1,1,2,NP4)
                  YBEZ(4)=XP(1,1,3,NP4)
                  XBEZ(2)=XB(1,2,N2LI)
                  YBEZ(2)=XB(1,3,N2LI)
                  XBEZ(3)=XB(2,2,N2LI)
                  YBEZ(3)=XB(2,3,N2LI)
                  DXI(1,2)=3.d0*(XBEZ(2)-XBEZ(1))
                  DXI(1,3)=3.d0*(YBEZ(2)-YBEZ(1))
                  DXI(2,2)=3.d0*(XBEZ(4)-XBEZ(3))
                  DXI(2,3)=3.d0*(YBEZ(4)-YBEZ(3))
                ELSE IF(iw.EQ.3) THEN
C old MPN unused?
                  CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
c                  CALL FBIMGPT(IV1,REAL(XP(1,1,1,NP1)),
c     '              REAL(XP(1,1,2,NP1)),
c     '              REAL(XP(1,1,3,NP1)),XBEZ(1),YBEZ(1))
c                  CALL FBIMGPT(IV1,REAL(XP(1,1,1,NP4)),
c     '              REAL(XP(1,1,2,NP4)),
c     '              REAL(XP(1,1,3,NP4)),XBEZ(4),YBEZ(4))
c                  CALL FBIMGPT(IV1,REAL(XB(1,1,N2LI)),
c     '              REAL(XB(1,2,N2LI)),REAL(XB(1,3,N2LI)),
c     '              XBEZ(2),YBEZ(2))
c                  CALL FBIMGPT(IV1,REAL(XB(2,1,N2LI)),
c     '              REAL(XB(2,2,N2LI)),REAL(XB(2,3,N2LI)),
c     '              XBEZ(3),YBEZ(3))
                  DXI(1,1)=3.d0*(XB(1,1,N2LI)-XP(1,1,1,NP1))
                  DXI(1,2)=3.d0*(XB(1,2,N2LI)-XP(1,1,2,NP1))
                  DXI(1,3)=3.d0*(XB(1,3,N2LI)-XP(1,1,3,NP1))
                  DXI(2,1)=3.d0*(XP(1,1,1,NP4)-XB(2,1,N2LI))
                  DXI(2,2)=3.d0*(XP(1,1,2,NP4)-XB(2,2,N2LI))
                  DXI(2,3)=3.d0*(XP(1,1,3,NP4)-XB(2,3,N2LI))
                ENDIF
                CALL DXIDL(IDO,NBJ,NEL(0,nl),NPL(1,0,N2LI),
     '            DL(1,N2LI),DXI,XP,ERROR,*9999)
                IF(iw.EQ.1) THEN
                  XP(NPL(4,1,N2LI),1,1,NP1)=DXI(1,1)/DL(3,N2LI)
                  XP(NPL(4,1,N2LI),1,3,NP1)=DXI(1,3)/DL(3,N2LI)
                  XP(NPL(5,1,N2LI),1,1,NP4)=DXI(2,1)/DL(3,N2LI)
                  XP(NPL(5,1,N2LI),1,3,NP4)=DXI(2,3)/DL(3,N2LI)
                ELSE IF(iw.EQ.2) THEN
                  XP(NPL(4,1,N2LI),1,2,NP1)=DXI(1,2)/DL(3,N2LI)
                  XP(NPL(4,1,N2LI),1,3,NP1)=DXI(1,3)/DL(3,N2LI)
                  XP(NPL(5,1,N2LI),1,2,NP4)=DXI(2,2)/DL(3,N2LI)
                  XP(NPL(5,1,N2LI),1,3,NP4)=DXI(2,3)/DL(3,N2LI)
                ELSE IF(iw.EQ.3) THEN
                  XP(NPL(4,1,N2LI),1,1,NP1)=DXI(1,1)/DL(3,N2LI)
                  XP(NPL(4,1,N2LI),1,2,NP1)=DXI(1,2)/DL(3,N2LI)
                  XP(NPL(4,1,N2LI),1,3,NP1)=DXI(1,3)/DL(3,N2LI)
                  XP(NPL(5,1,N2LI),1,1,NP4)=DXI(2,1)/DL(3,N2LI)
                  XP(NPL(5,1,N2LI),1,2,NP4)=DXI(2,2)/DL(3,N2LI)
                  XP(NPL(5,1,N2LI),1,3,NP4)=DXI(2,3)/DL(3,N2LI)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          nr=1 !temporary
          nx=1 !temporary
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','RED')
          CALL SGLINE(INDEX,ISEG,ISLINE(iw,NOLINE),ISLINO(iw),iw,NLLIST,
     '      NOLINE,NPL,nr,nx,CSEG,'UNDEFORMED',DL,'SOLID',.TRUE.,XP,ZP,
     '      ERROR,*9999)
          IF(RHTRAN.OR.LFTRAN) THEN
C old MPN unused?
            CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
c            CALL FBCLEAR(IV1)
            CALL SGLINE(INDEX,ISEG,ISLINE(IV1,NOLINE),ISLINO(IV1),IV1,
     '        NLLIST,NOLINE,NPL,nr,nx,
     '        CSEG,'UNDEFORMED',DL,'SOLID',.TRUE.,XP,ZP,ERROR,*9999)
          ENDIF
        ENDIF
      ENDDO
      CALL DAWK(iw,0,ERROR,*9999)

      CALL ACWK(iw,1,ERROR,*9999)
      DO 850 nolist=1,NLLIST(0)
        nl=NLLIST(nolist)
        IF(NPL(1,1,nl).EQ.4) THEN
          IF(ISEG(ISN2BE(nl)).GT.0) THEN
            CALL DELETE_SEGMENT(ISN2BE(nl),ISEG,iw,ERROR,*9999)
          ENDIF
          IF(ISEG(ISN3BE(nl)).GT.0) THEN
            CALL DELETE_SEGMENT(ISN3BE(nl),ISEG,iw,ERROR,*9999)
          ENDIF
          IF(ISEG(ISL2BE(nl)).GT.0) THEN
            CALL DELETE_SEGMENT(ISL2BE(nl),ISEG,iw,ERROR,*9999)
          ENDIF
          IF(ISEG(ISL3BE(nl)).GT.0) THEN
            CALL DELETE_SEGMENT(ISL3BE(nl),ISEG,iw,ERROR,*9999)
          ENDIF
        ENDIF
 850  CONTINUE
      CALL DAWK(iw,1,ERROR,*9999)

      CALL EXITS('CRHERM')
      RETURN
 9999 CALL ERRORS('CRHERM',ERROR)
      CALL EXITS('CRHERM')
      RETURN 1
      END

