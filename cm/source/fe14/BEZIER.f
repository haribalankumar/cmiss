      SUBROUTINE BEZIER(index,INDEX_PLIN,ISBEZE,ISEG,ISL2BE_LOC,
     '  ISL3BE_LOC,ISN2BE_LOC,ISN3BE_LOC,iw,NTL,
     '  XPTS,YPTS,CSEG,ERROR,*)

C#### Subroutine: BEZIER
C###  Description:
C###    BEZIER creates Bezier control points on workstation iw.

C**** INDEX is the GKS line type index
C**** INDEX_PLIN is the polyline index
C**** XBEZ(I),YBEZ(I),I=1 & 4 are x,y coords of nodes.
C**** XBEZ(I),YBEZ(I),I=2 & 3 are x,y coords of slope control points.
C**** Note: This routine creates temporary segments which are deleted at
C**** the end and the total segment count NTSG is reduced again.
C GMH 18-4-96 Shouldn't these arrays (ISL3BE etc) be passed in?
C AJP 19-5-96 Probably - for now change the name to reflect the fact that
C the array is local.
C DPN 17-12-97 - modify so these arrays (ISL3BE etc) are passed in

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER NTLM
      PARAMETER(NTLM=20)
      INTEGER index,INDEX_PLIN,ISBEZE,ISEG(*),ISL2BE_LOC(NTLM),
     '  ISL3BE_LOC(NTLM),ISN2BE_LOC(NTLM),ISN3BE_LOC(NTLM),iw,NTL
      REAL*8 XPTS(*),YPTS(*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER i,IFROMC,INDEX_OLD,INSTAT,IPICK,isegm,N2LI,nl
C DPN 17-12-97 - these arrays now passed in
C      INTEGER ISL2BE_LOC(NTLM),ISL3BE_LOC(NTLM),ISN2BE_LOC(NTLM),
C     '  ISN3BE_LOC(NTLM)
      REAL*8 D_XWC,D_YWC
      REAL*8 BEZLEN(NTLM),PL(3,21),PT(3,2),XBEZ(4,NTLM),XSLOPE,
     '  YBEZ(4,NTLM),YSLOPE
      SAVE BEZLEN,XBEZ,YBEZ
C DPN 17-12-97 - these arrays now passed in
C      SAVE ISL2BE_LOC,ISL3BE_LOC,ISN2BE_LOC,ISN3BE_LOC

C      DATA LD1/1/,BLANK/' '/

      CALL ENTERS('BEZIER',*9999)
      CALL ACWK(iw,0,ERROR,*9999)

      CALL ASSERT(NTL.LE.NTLM,
     '  '>> Increase number of lines (NTLM) in BEZIER',ERROR,*9999)
      DO nl=1,NTL
        IF(DOP) THEN
          WRITE(OP_STRING,*) 'nl=',nl
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C ***   Draw polymarkers at Bezier control points and tangent lines
        DO i=1,4
          XBEZ(I,nl)=XPTS(3*(nl-1)+i)
          YBEZ(I,nl)=YPTS(3*(nl-1)+i)
        ENDDO
        BEZLEN(nl)=DSQRT((XBEZ(4,nl)-XBEZ(1,nl))**2+
     '    (YBEZ(4,nl)-YBEZ(1,nl))**2)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' XBEZ: '',4E12.3)') (XBEZ(i,nl),i=1,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' YBEZ: '',4E12.3)') (YBEZ(i,nl),i=1,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        CALL OPEN_SEGMENT(ISN2BE_LOC(nl),ISEG,iw,'CPT1',
     '    index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
        PT(1,2)=XBEZ(2,nl)
        PT(2,2)=YBEZ(2,nl)
        CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISN2BE_LOC(nl),iw,ERROR,*9999)
        CALL DETECT(iw,ISEG,ISN2BE_LOC(nl),'DETECTABLE',ERROR,*9999)

        CALL OPEN_SEGMENT(ISL2BE_LOC(nl),ISEG,iw,'TAN1',
     '    index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
        PT(1,1)=XBEZ(1,nl)
        PT(2,1)=YBEZ(1,nl)
        CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISL2BE_LOC(nl),iw,ERROR,*9999)

        CALL OPEN_SEGMENT(ISN3BE_LOC(nl),ISEG,iw,'CPT2',
     '    index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
        PT(1,1)=XBEZ(3,nl)
        PT(2,1)=YBEZ(3,nl)
        CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISN3BE_LOC(nl),iw,ERROR,*9999)
        CALL DETECT(iw,ISEG,ISN3BE_LOC(nl),'DETECTABLE',ERROR,*9999)

        CALL OPEN_SEGMENT(ISL3BE_LOC(nl),ISEG,iw,'TAN2',
     '    index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
        PT(1,2)=XBEZ(4,nl)
        PT(2,2)=YBEZ(4,nl)
        CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISL3BE_LOC(nl),iw,ERROR,*9999)
      ENDDO

      INSTAT=1
      DO WHILE(INSTAT.EQ.1)
C ***   Pick control point segments
        WRITE(OP_STRING,
     '    '('' >>Pick control point on '',I1,'':'')') iw
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL PICK(iw,'REQUEST',INSTAT,isegm,IPICK,ERROR,*9999)
        IF(INSTAT.EQ.1) THEN
          IF(CSEG(isegm)(1:4).EQ.'CPT1') THEN
            nl=IFROMC(CSEG(isegm)(53:57))
C ***       Locate new control point
            WRITE(OP_STRING,
     '        '('' >>Relocate 1st control point for line '','
     '        //'I4,'' on '',I1,'':'')') nl,iw
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(INSTAT,0.0D0,D_XWC,0.0D0,
     '        D_YWC,ERROR,*9999)
            XBEZ(2,nl)=D_XWC
            YBEZ(2,nl)=D_YWC

            CALL OPEN_SEGMENT(ISN2BE_LOC(nl),ISEG,iw,'CPT1',
     '        index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,2)=XBEZ(2,nl)
            PT(2,2)=YBEZ(2,nl)
            CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISN2BE_LOC(nl),iw,ERROR,*9999)
            CALL DETECT(iw,ISEG,ISN2BE_LOC(nl),'DETECTABLE',ERROR,*9999)

            CALL OPEN_SEGMENT(ISL2BE_LOC(nl),ISEG,iw,'TAN1',
     '        index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,1)=XBEZ(1,nl)
            PT(2,1)=YBEZ(1,nl)
            CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISL2BE_LOC(nl),iw,ERROR,*9999)

C ***       Find adjacent control point if colinear
            IF(nl.ne.1) THEN
              N2LI=nl-1
              XSLOPE=(XBEZ(2,nl)-XBEZ(1,nl))/BEZLEN(nl)
              YSLOPE=(YBEZ(2,nl)-YBEZ(1,nl))/BEZLEN(nl)
              XBEZ(3,N2LI)=XBEZ(1,nl)-XSLOPE*BEZLEN(N2LI)
              YBEZ(3,N2LI)=YBEZ(1,nl)-YSLOPE*BEZLEN(N2LI)

              CALL OPEN_SEGMENT(ISN3BE_LOC(N2LI),ISEG,iw,'CPT2',
     '          index,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,1)=XBEZ(3,N2LI)
              PT(2,1)=YBEZ(3,N2LI)
              CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISN3BE_LOC(N2LI),iw,ERROR,*9999)
              CALL DETECT(iw,ISEG,ISN3BE_LOC(N2LI),'DETECTABLE',
     '          ERROR,*9999)

              CALL OPEN_SEGMENT(ISL3BE_LOC(N2LI),ISEG,iw,'TAN2',
     '          index,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,2)=XBEZ(4,N2LI)
              PT(2,2)=YBEZ(4,N2LI)
              CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISL3BE_LOC(N2LI),iw,ERROR,*9999)
            ENDIF

          ELSE IF(CSEG(isegm)(1:4).EQ.'CPT2') THEN
            nl=IFROMC(CSEG(isegm)(53:57))
C ***       Locate new control point
            WRITE(OP_STRING,
     '        '('' >>Relocate 2nd control point for line '','
     '        //'I4,'' on '',I1,'':'')') nl,iw
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(INSTAT,0.0D0,D_XWC,0.0D0,
     '        D_YWC,ERROR,*9999)
            XBEZ(3,nl)=D_XWC
            YBEZ(3,nl)=D_YWC

            CALL OPEN_SEGMENT(ISN3BE_LOC(nl),ISEG,iw,'CPT2',
     '        index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,1)=XBEZ(3,nl)
            PT(2,1)=YBEZ(3,nl)
            CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISN3BE_LOC(nl),iw,ERROR,*9999)
            CALL DETECT(iw,ISEG,ISN3BE_LOC(nl),'DETECTABLE',ERROR,*9999)

            CALL OPEN_SEGMENT(ISL3BE_LOC(nl),ISEG,iw,'TAN2',
     '        index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,2)=XBEZ(4,nl)
            PT(2,2)=YBEZ(4,nl)
            CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISL3BE_LOC(nl),iw,ERROR,*9999)

C ***       Find adjacent control point if colinear
            IF(nl.ne.NTL) THEN
              N2LI=nl+1
              XSLOPE=(XBEZ(4,nl)-XBEZ(3,nl))/BEZLEN(nl)
              YSLOPE=(YBEZ(4,nl)-YBEZ(3,nl))/BEZLEN(nl)
              XBEZ(2,N2LI)=XBEZ(4,nl)+XSLOPE*BEZLEN(N2LI)
              YBEZ(2,N2LI)=YBEZ(4,nl)+YSLOPE*BEZLEN(N2LI)

              CALL OPEN_SEGMENT(ISN2BE_LOC(N2LI),ISEG,iw,'CPT1',
     '          index,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,2)=XBEZ(2,N2LI)
              PT(2,2)=YBEZ(2,N2LI)
              CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISN2BE_LOC(N2LI),iw,ERROR,*9999)
              CALL DETECT(iw,ISEG,ISN2BE_LOC(N2LI),'DETECTABLE',
     '          ERROR,*9999)

              CALL OPEN_SEGMENT(ISL2BE_LOC(N2LI),ISEG,iw,'TAN1',
     '          index,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,1)=XBEZ(1,N2LI)
              PT(2,1)=YBEZ(1,N2LI)
              CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISL2BE_LOC(N2LI),iw,ERROR,*9999)
            ENDIF
          ENDIF

C ***     Calculate line length and redraw line segments
          CALL OPEN_SEGMENT(ISBEZE,ISEG,iw,'polyline_bezier',index,
     '      INDEX_OLD,INDEX_PLIN,1,CSEG,ERROR,*9999)
          DO nl=1,NTL
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' nl= '',I3)') nl
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL BEZIER_POINTS(PL,XBEZ(1,nl),YBEZ(1,nl),ERROR,*9999)
            CALL POLYLINE(1,iw,21,PL,ERROR,*9999)
          ENDDO
          CALL CLOSE_SEGMENT(ISBEZE,iw,ERROR,*9999)
        ENDIF
      ENDDO
      CALL DAWK(iw,0,ERROR,*9999)

      CALL ACWK(iw,1,ERROR,*9999)
      DO nl=NTL,1,-1
        IF(ISEG(ISL3BE_LOC(nl)).GT.0) THEN
          CALL DELETE_SEGMENT(ISL3BE_LOC(nl),ISEG,iw,ERROR,*9999)
          NTSG=NTSG-1 !reduce total segment count by 1
        ENDIF
        IF(ISEG(ISN3BE_LOC(nl)).GT.0) THEN
          CALL DELETE_SEGMENT(ISN3BE_LOC(nl),ISEG,iw,ERROR,*9999)
          NTSG=NTSG-1 !reduce total segment count by 1
        ENDIF
        IF(ISEG(ISL2BE_LOC(nl)).GT.0) THEN
          CALL DELETE_SEGMENT(ISL2BE_LOC(nl),ISEG,iw,ERROR,*9999)
          NTSG=NTSG-1 !reduce total segment count by 1
        ENDIF
        IF(ISEG(ISN2BE_LOC(nl)).GT.0) THEN
          CALL DELETE_SEGMENT(ISN2BE_LOC(nl),ISEG,iw,ERROR,*9999)
          NTSG=NTSG-1 !reduce total segment count by 1
        ENDIF
        DO i=1,4
          XPTS(3*(nl-1)+i)=XBEZ(i,nl)
          YPTS(3*(nl-1)+i)=YBEZ(i,nl)
        ENDDO
      ENDDO
      CALL DAWK(iw,1,ERROR,*9999)

      CALL EXITS('BEZIER')
      RETURN
 9999 CALL ERRORS('BEZIER',ERROR)
      CALL EXITS('BEZIER')
      RETURN 1
      END


