      SUBROUTINE DROBJE(ISEG,ISOBJE,LD,MXI,NDDL,NDLT,NEELEM,
     '  XID,ZD,ZDD,WD,CSEG,STRING,ERROR,*)

C#### Subroutine: DROBJE
C###  Description:
C###    DROBJE draws graphical objects noobje=1,ntobje.
C**** NSOBJE(1,noobje) is segment number of object
C**** NSOBJE(2,noobje) is number of parts to object
C****   (defined by intersections with other objects)
C**** NSOBJE(3..,noobje) are nd numbers of object 'corners'

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISOBJE(NWM,NGRSEGM,NGRSEGM),LD(NDM),MXI(2,NEM),
     '  NDDL(NEM,NDEM),NDLT(NEM),NEELEM(0:NE_R_M,0:NRM)
      REAL*8 XID(NIM,NDM),ZD(NJM,NDM),ZDD(NJM,NDM),WD(NJM,NDM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER i,IBEG,IBEG5,IDATA(2),ID_SEGMENT,IEND,
     '  IEND5,IFROMC,INDEX,INSTAT,IPICK,ISEGM,ITOT,iw,IWK(6),
     '  n1part,N3CO,naobje,NAPART,nbobje,NBPART,nd,ND1,ND2,ND3,
     '  nda,NDA1,NDA2,ndb,NDB1,NDB2,NDI,ne,NEE,nj,
     '  noelem,noiw,NOOBJE,nopart,nosg,nr,NTIW,NTPART
c     INTEGER ID_DEVICE,ID_WS
      REAL*8 DISTA,DISTB,XID1,XID2,XWC,YWC,ZWC
c     REAL*8 RDATA(1)
      CHARACTER CHAR5*5,CLASS*8,NAME*20
c     CHARACTER SDATA*10
      LOGICAL ABBREV,CBBREV,CONTINUE,FOUND,MOUSE,OBJECT_DEFINED,SEGME

C      DATA LD1/1/

      CALL ENTERS('DROBJE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw object;m/s
C###  Parameter:    <on WS_ID#[1]>
C###    Specify the workstation (GX window) to draw the
C###    objects on.
C###  Description:
C###    Draws graphical objects, "m" option
C###    specifies mouse and "s" specifies the segment
C###    intersecrions
C###

        OP_STRING(1)=STRING(1:IEND)//';m/s'
        OP_STRING(2)=BLANK(1:15)//'<on WS_ID#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw object
C###  Parameter:    <with SEGMENT_ID#[0]>
C###    Specify the segment number of the object
C###  Parameter:    <name=OBJECT_NAME[DATA]>
C###    Specify the name of the object
C###  Description:
C###    Draws graphical objects

C        CHAR5=CFROMI(NTSG,'(I5)')
        WRITE(CHAR5,'(I5)') NTSG
        CALL STRING_TRIM(CHAR5,IBEG5,IEND5)
        OP_STRING(1)=STRING(1:IEND)
     '    //'<with SEGMENT_ID#['//CHAR5(IBEG5:IEND5)//']>'
        OP_STRING(2)=BLANK(1:15)//'<name=OBJECT_NAME[DATA]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw object from DATA_NUMBER to DATA_NUMBER
C###  Parameter:    <name=OBJECT_NAME[DATA]>
C###  Description:
C###    Draws a range of objects specified by DATA_NUMBER range

        OP_STRING(1)=STRING(1:IEND)//' from DATA_NUMBER to DATA_NUMBER'
        OP_STRING(2)=BLANK(1:15)//'<name=OBJECT_NAME[DATA]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DROBJE',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL_OBJE=.TRUE.
        MOUSE=.FALSE.
        SEGME=.FALSE.
        IF(ABBREV(COQU(noco,1),'M',1)) THEN
          MOUSE=.TRUE.
        ELSE IF(ABBREV(COQU(noco,1),'S',1)) THEN
          SEGME=.TRUE.
        ENDIF

        IF(MOUSE.OR.SEGME) THEN
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        ELSE
          IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
            ND1=IFROMC(CO(N3CO+1))
          ELSE
            ID_SEGMENT=0
          ENDIF
          IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
            ND2=IFROMC(CO(N3CO+1))
          ELSE
            ID_SEGMENT=0
          ENDIF
          IF(CBBREV(CO,'NAME',1,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            NAME=CO(N3CO+1)(IBEG:IEND)
          ELSE
            NAME='DATA'
          ENDIF
          IF(CBBREV(CO,'WITH',1,noco+1,NTCO,N3CO)) THEN
            ID_SEGMENT=IFROMC(CO(N3CO+1))
          ELSE
            ID_SEGMENT=NTSG
          ENDIF
        ENDIF

        IF(MOUSE) THEN
          IW=IWK(1)
C CPB 10/9/92 Changed workstation update mode
          CALL ACWK(iw,1,ERROR,*9999)
          DO nosg=1,NTSG
            IF(ISEG(nosg).GT.0.AND.CSEG(nosg)(1:13).EQ.'DRAW polyline')
     '        CALL DETECT(1,ISEG,nosg,'DETECTABLE',ERROR,*9999)
          ENDDO
          WRITE(OP_STRING,'('' >>Pick segment on '',I1,'
     '      //''' (use LH button to pick, middle button to end)'')') IW
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL PICK(iw,'EVENT',INSTAT,ISEGM,IPICK,ERROR,*9999)
          CONTINUE=.TRUE.
          DO WHILE (CONTINUE)
C GMH set class and idata so we dont get a warning
            CLASS='UNUSED'
            IDATA(1)=0
            IDATA(2)=0
c           CALL EVENT(ID_WS,ID_Device,INSTAT,CLASS,IDATA,RDATA,SDATA,
c    '        ERROR,*9999)
            IF(CLASS(1:4).EQ.'PICK') THEN
              ISEGM=IDATA(1)
              IF(DOP) THEN
                WRITE(OP_STRING,*) 'Input class is pick'
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(INSTAT.EQ.1) THEN
                IF(DOP) THEN
                  WRITE(OP_STRING,*) ' ISEGM=',ISEGM,' CSEG=',
     '              CSEG(ISEGM)(1:60)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                OBJECT_DEFINED=.FALSE.
                IF(CSEG(ISEGM)(1:20).EQ.'DRAW polyline_bezier') THEN
                  CALL OBJ_BEZIER(ISEGM,ITOT,ZD,WD,ERROR,*9999)
                  OBJECT_DEFINED=.TRUE.
                ELSE IF(CSEG(ISEGM)(1:17).EQ.'DRAW polyline_box') THEN
                  CALL OBJ_BOX(ISEGM,ITOT,ZD,WD,ERROR,*9999)
                  OBJECT_DEFINED=.TRUE.
                ELSE IF(CSEG(ISEGM)(1:20).EQ.'DRAW polyline_circle')
     '            THEN
                  CALL OBJ_CIRCLE(ISEGM,ITOT,ZD,WD,ERROR,*9999)
                  OBJECT_DEFINED=.TRUE.
                ELSE IF(CSEG(ISEGM)(1:21).EQ.'DRAW polyline_ellipse')
     '            THEN
                  CALL OBJ_ELLIPSE(ISEGM,ERROR,*9999)
                  OBJECT_DEFINED=.TRUE.
                ELSE IF(CSEG(ISEGM)(1:20).EQ.'DRAW polyline_locate')
     '            THEN
                  CALL OBJ_LINE(ISEGM,ITOT,ZD,WD,ERROR,*9999)
                  OBJECT_DEFINED=.TRUE.
                ELSE IF(CSEG(ISEGM)(1:20).EQ.'DRAW polyline_stroke')
     '            THEN
                  CALL OBJ_STROKE(ERROR,*9999)
                  OBJECT_DEFINED=.TRUE.
                ENDIF
                IF(OBJECT_DEFINED) THEN
                  IF(IW.EQ.4.AND.MAP_PROJEC(1:2).EQ.'XI') THEN
                    DO i=1,ITOT
                      XWC=ZD(1,NDT+i)
                      YWC=ZD(2,NDT+i)
                      DO nr=1,NRT
                        DO noelem=1,NEELEM(0,nr)
                          NEE=NEELEM(noelem,nr)
                          XID1=DBLE(MAX_XI)*(XWC+1.0D0)/2.0D0
     '                                     -DBLE(MXI(1,NEE)-1)
                          XID2=DBLE(MAX_XI)*(YWC+1.0D0)/2.0D0
     '                                     -DBLE(MXI(2,NEE)-1)
                          IF(XID1.GE.0.0D0.AND.XID1.LT.1.0D0.AND.
     '                       XID2.GE.0.0D0.AND.XID2.LT.1.0D0) THEN
                            XID(1,NDT+I)=XID1
                            XID(2,NDT+I)=XID2
                            LD(NDT+I)=NEE
                            GO TO 400
                          ENDIF
                        ENDDO
                      ENDDO
 400                  ne=NEE !is element no for current data point
                      MXI1=MXI(1,ne)
                      MXI2=MXI(2,ne)
                      CALL POLYMARKER(1,iw,1,XID(1,NDT+I),ERROR,*9999)
                    ENDDO
                    NDT=NDT+ITOT
                    CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)
                    NXIDEF=2
                  ELSE
                    DO i=1,ITOT
                      ZDD(1,i)=ZD(1,NDT+i)
                      ZDD(2,i)=ZD(2,NDT+i)
                    ENDDO
                    CALL POLYMARKER(1,iw,ITOT,ZDD,ERROR,*9999)
                    NDT=NDT+ITOT
                  ENDIF
                  CALL DELETE_SEGMENT(ISEGM,ISEG,iw,ERROR,*9999)
                ENDIF

              ELSE
c               CALL INPUT_MODE(iw,LD1,'PICK','REQUEST',ERROR,*9999)
                CONTINUE=.FALSE.
              ENDIF
            ENDIF
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)

        ELSE IF(SEGME) THEN
C ***     Calculate object intersections
          DO naobje=1,NTOBJE
C GMH 2/9/95 Unused            NOSEGM=NSOBJE(1,naobje) !is segment number of object naobje
            NAPART=NSOBJE(2,naobje) !is number of parts to object naobje
            nopart=0
            FOUND=.FALSE.
            DO WHILE (nopart.LT.NAPART)
              nopart=nopart+1
              WRITE(OP_STRING,
     '          '('' Object no.='',I3,'' part no.='',I3)')
     '          naobje,nopart
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              NDA1=NSOBJE(2+nopart,naobje)
              IF(FOUND) NDA1=NDA1+2 !to avoid picking same intersection
              NDA2=NSOBJE(3+nopart,naobje)
              DISTA=(ZD(1,NDA1+1)-ZD(1,NDA1))**2
     '             +(ZD(2,NDA1+1)-ZD(2,NDA1))**2
              IF(DOP) THEN
                WRITE(OP_STRING,'('' nda1='',I3,'' nda2='',I3,'
     '            //''' dista='',E12.3)') NDA1,NDA2,DISTA
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              FOUND=.FALSE.
              DO nda=NDA1,NDA2
                DO nbobje=1,NTOBJE
                  IF(nbobje.NE.naobje) THEN
C GMH 2/9/95 Unused                    NOSEGM=NSOBJE(1,nbobje) !is segment number of object nbobje
                    NBPART=NSOBJE(2,nbobje) !is number of parts to object nbobje
                    NDB1=NSOBJE(3,nbobje)
                    NDB2=NSOBJE(3+NBPART,nbobje)
                    IF(DOP) THEN
                      WRITE(IOOP,'('' nbpart='',I3,'' nda='',I3,'
     '                  //''' ndb1='',I3,'' ndb2='',I3)') NBPART,nda,
     '                NDB1,NDB2
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    DO ndb=NDB1,NDB2
                      DISTB=(ZD(1,nda)-ZD(1,ndb))**2
     '                     +(ZD(2,nda)-ZD(2,ndb))**2
                      IF(DABS(DISTB).LT.DABS(DISTA)) THEN
                        WRITE(OP_STRING,
     '                    '('' Intersection found at nda='',I3,'
     '                    //''' ndb='',I3,'' distb='',E12.3)') nda,ndb,
     '                    DISTB
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        NDI=nda !is nd number of object naobje at intersection
                        FOUND=.TRUE.
                        GO TO 20
                      ENDIF
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
 20           IF(FOUND) THEN
                DO n1part=NAPART,nopart,-1 !to shift nd numbers up
                  NSOBJE(3+n1part+1,naobje)=NSOBJE(3+n1part,naobje)
                ENDDO
                NSOBJE(3+nopart,naobje)=NDI
                NSOBJE(2,naobje)=NSOBJE(2,naobje)+1 !updates number of parts
                NAPART=NAPART+1
              ENDIF
            ENDDO
          ENDDO

C ***     Display object parts
          DO nd=1,NDT
            IF(IW.EQ.4.AND.MAP_PROJEC(1:2).EQ.'XI') THEN
            ELSE
              DO nj=1,NJT
                ZDD(nj,nd)=ZD(nj,nd)
              ENDDO
            ENDIF
          ENDDO
          DO noiw=1,NTIW
            IW=IWK(noiw)
            CALL ACWK(iw,1,ERROR,*9999)
            DO NOOBJE=1,NTOBJE
C GMH 2/9/95 Unused              NOSEGM=NSOBJE(1,NOOBJE) !is segment number of complete object
              NTPART=NSOBJE(2,NOOBJE) !is number of parts to object
              DO nopart=1,NTPART
                ND1=NSOBJE(2+nopart,NOOBJE) !is 1st nd number of object part
                ND2=NSOBJE(3+nopart,NOOBJE) !is 2nd nd number of object part
                ND3=INT(0.5*REAL(ND1+ND2))  !is nd number of centre of part
                IF(DOP) THEN
                  WRITE(OP_STRING,
     '              '('' nd1='',I3,'' nd2='',I3,'' nd3='',I3)')
     '              ND1,ND2,ND3
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                XWC=ZD(1,ND3)
                YWC=ZD(2,ND3)
                ZWC=ZD(3,ND3)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Coords='',3E12.3)')
     '              XWC,YWC,ZWC
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                CALL SGOBJE(INDEX,ISEG,ISOBJE(iw,NOOBJE,nopart),iw,ND1,
     '            ND2,NOOBJE,nopart,CSEG,ZDD,ERROR,*9999)
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO

        ELSE !define object from existing data or segment
          IF(ID_SEGMENT.GT.0) THEN
            ISEGM=ID_SEGMENT
            IF(DOP) THEN
               WRITE(OP_STRING,*)
     '           ' ISEGM=',ISEGM,' CSEG=',CSEG(ISEGM)(1:60)
               CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            OBJECT_DEFINED=.FALSE.
            IF(CSEG(ISEGM)(1:20).EQ.'DRAW polyline_bezier') THEN
              CALL OBJ_BEZIER(ISEGM,ITOT,ZD,WD,ERROR,*9999)
              OBJECT_DEFINED=.TRUE.
            ELSE IF(CSEG(ISEGM)(1:17).EQ.'DRAW polyline_box') THEN
              CALL OBJ_BOX(ISEGM,ITOT,ZD,WD,ERROR,*9999)
              OBJECT_DEFINED=.TRUE.
            ELSE IF(CSEG(ISEGM)(1:20).EQ.'DRAW polyline_circle') THEN
              CALL OBJ_CIRCLE(ISEGM,ITOT,ZD,WD,ERROR,*9999)
              OBJECT_DEFINED=.TRUE.
            ELSE IF(CSEG(ISEGM)(1:21).EQ.'DRAW polyline_ellipse') THEN
              CALL OBJ_ELLIPSE(ISEGM,ERROR,*9999)
              OBJECT_DEFINED=.TRUE.
            ELSE IF(CSEG(ISEGM)(1:20).EQ.'DRAW polyline_locate') THEN
              CALL OBJ_LINE(ISEGM,ITOT,ZD,WD,ERROR,*9999)
              OBJECT_DEFINED=.TRUE.
            ELSE IF(CSEG(ISEGM)(1:20).EQ.'DRAW polyline_stroke') THEN
              CALL OBJ_STROKE(ERROR,*9999)
              OBJECT_DEFINED=.TRUE.
            ENDIF
            IF(OBJECT_DEFINED) THEN
              IF(IW.EQ.4.AND.MAP_PROJEC(1:2).EQ.'XI') THEN
                DO i=1,ITOT
                  XWC=ZD(1,NDT+i)
                  YWC=ZD(2,NDT+i)
                  DO nr=1,NRT
                    DO noelem=1,NEELEM(0,nr)
                      NEE=NEELEM(noelem,nr)
                      XID1=DBLE(MAX_XI)*(XWC+1.0D0)/2.0D0-
     '                  DBLE(MXI(1,NEE)-1)
                      XID2=DBLE(MAX_XI)*(YWC+1.0D0)/2.0D0-
     '                  DBLE(MXI(2,NEE)-1)
                      IF(XID1.GE.0.0D0.AND.XID1.LT.1.0D0.AND.
     '                   XID2.GE.0.0D0.AND.XID2.LT.1.0D0) THEN
                        XID(1,NDT+i)=XID1
                        XID(2,NDT+i)=XID2
                        LD(NDT+I)=NEE
                        GO TO 500
                      ENDIF
                    ENDDO
                  ENDDO
 500              CONTINUE
                ENDDO
              ENDIF
              NDT=NDT+ITOT
              CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)
              NXIDEF=2
            ENDIF

          ELSE
            NTOBJE=NTOBJE+1
            OBJECT_NAME(NTOBJE)=NAME
            NSOBJE(1,NTOBJE)=0   !is segment number of object
            NSOBJE(2,NTOBJE)=1   !is number of parts to object
            NSOBJE(3,NTOBJE)=ND1 !is first nd in object
            NSOBJE(4,NTOBJE)=ND2 !is last  nd in object
          ENDIF

        ENDIF
      ENDIF

      CALL EXITS('DROBJE')
      RETURN
 9999 CALL ERRORS('DROBJE',ERROR)
      CALL EXITS('DROBJE')
      RETURN 1
      END


