      SUBROUTINE SGDATA(INDEX,IBT,IDO,INP,ISDANO,ISDAPR,ISDATA,ISDATR,
     '  ISEG,iw,LD,MXI,NBJ,NBH,NDDL,NDLT,NDP,NEELEM,NELIST,
     '  NKHE,NKJE,NPF,NPNE,NRE,NVHE,NVJE,NW,
     '  CE,CG,CGE,CP,CSEG,CURVCORRECT,PG,SE,STATIC,TITLE,TYPE,WD,WDL,
     '  XA,XE,XG,XID,XIDL,XP,ZA,ZD,ZDD,ZDL,ZDLMIN,ZE,ZP,ERROR,*)

C#### Subroutine: SGDATA
C###  Description:
C###    SGDATA creates data segments ISDATA(iw),ISDATR(iw,ne),
C###    ISDANO(iw,ne) and ISDAPR(iw,ne).  Draws lines connecting each
C###    data point ZD with its corresponding point on the ensemble.
C###    The data points are identified by small circles. The points on
C###    the ensemble are identified by a cross corresponding to small
C###    increments in the element coordinates XID.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'colo00.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'four00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'scal00.cmn'
      INCLUDE 'scal01.cmn'
      INCLUDE 'time00.cmn'
      INCLUDE 'trans00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ISDANO(NWM,NEM),ISDAPR(NWM,NEM),ISDATA,
     '  ISDATR(NWM,NEM),ISEG(*),iw,LD(NDM),MXI(2,NEM),
     '  NBJ(NJM,NEM),
     '  NBH(NHM,NCM,NEM),NDDL(NEM,NDEM),NDLT(NEM),NDP(NDM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),WD(NJM,NDM),WDL(NHM,NDEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     '  XIDL(NIM,NDEM),XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM),
     '  ZD(NJM,NDM),ZDD(NJM,NDM),ZDL(NHM,NDEM),ZDLMIN,ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),TITLE*(*),TYPE*(*)
      LOGICAL STATIC
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INDEX_OLD,n,nb,NBO,nc,nd,NDD,nde,ne,ng,
     '  ni,NITB,nj,noelem,nolist,NOPARA,nr,ns,NTNDE,nx
      REAL*8 POINT(3),PXI,T1,THETA1,THETA2,THETA_ND,X(3),
     '  XD(3),XI(3),Z(3),Z2(3),
     '  ZVAL
C      REAL*8 TEXT_POINT(3),XXMAX,XXMIN,YYMAX,YYMIN,ZZMAX,ZZMIN
      CHARACTER CHAR*10,CLABEL*52
      LOGICAL COMPUTE,DRAW,INLIST,XI_DEFINED

      CALL ENTERS('SGDATA',*9999)
      nc=1 !Temporary MPN 12-Nov-94
      nr=1 !may need generalizing
      nx=1 !temporary

      COMPUTE=.FALSE.
      IF(INDEX.EQ.0) COMPUTE=.TRUE.

      IF(TYPE(1:8).EQ.'GEOMETRY'.OR
     '  .TYPE(1:6).EQ.'SHEETS'.OR.TYPE(1:6).EQ.'FIBRES'.OR
     '  .TYPE(1:5).EQ.'FIELD'.OR.TYPE(1:6).EQ.'SERIES'.OR
     '  .TYPE(1:9).EQ.'PARAMETER'.OR.TYPE(1:7).EQ.'HISTORY'.OR
     '  .TYPE(1:5).EQ.'SLICE') THEN
        IF(DOP) THEN
          WRITE(OP_STRING,*)' sgdata/iw=',iw
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)' sgdata/nelist=',
     '      (NELIST(nolist),nolist=1,NELIST(0))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

C CPB 8/4/94 This needs generalising for NJ_LOC
        CLABEL='DATA'
        IF(STATIC) THEN
          CLABEL(43:52)='          '
        ELSE
          WRITE(CLABEL(43:52),'(E10.3)') TIME
        ENDIF
        CALL OPEN_SEGMENT(ISDATA,ISEG,iw,CLABEL,INDEX,INDEX_OLD,
     '    NTDATA,1,CSEG,ERROR,*9999)

C GBS 21/10/94
C Changed to draw individual data points so that colours can be
C computed according to the value.
        IF(TYPE(1:8).EQ.'GEOMETRY'.OR.TYPE(1:5).EQ.'FIELD'.OR
     '    .TYPE(1:6).EQ.'SERIES'.OR.TYPE(1:7).EQ.'HISTORY') THEN
          XI_DEFINED=.FALSE.
          DO nd=1,NDT
            IF(LD(nd).GT.0) XI_DEFINED=.TRUE.
          ENDDO
          NDD=0
C GBS 21/10/94
C Changed to draw individual data points so that colours can be
C computed according to the value.
          DO nd=1,NDT
            DRAW=.TRUE.
            IF(COMPUTE) THEN  !From "rgb=value"
C GBS 19/10/94
C Should extend later to call RANGE if ZMINI=ZMAXI (i.e. not set)
              ZVAL=ZD(NJT+1,nd)
              IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
              IF(COLOUR_WS) THEN
                INDEX=256-NINT((ZVAL-ZMINI)/ZDIFF*215.0D0)
              ELSE
                INDEX=8-NINT((ZVAL-ZMINI)/ZDIFF*3.0D0)
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' nd, index='',2I4)') nd,INDEX
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
            IF(TYPE(1:7).NE.'HISTORY') THEN
              IF(XI_DEFINED) THEN
                IF(STATIC) THEN
                  DRAW=.FALSE.
                  ne=LD(nd)
                  IF(INLIST(ne,NELIST(1),NELIST(0),N)) THEN
                    IF(NXIDEF.LT.3) THEN
                      DRAW=.TRUE.
                    ELSE
                      DRAW=XID(3,nd).GE.XI3MIN.AND.XID(3,nd).LE.XI3MAX
                    ENDIF
                    IF(DRAW) THEN
                      IF(iw.EQ.4.AND.MAP_PROJEC(1:2).EQ.'XI') THEN
                        MXI1=MXI(1,ne) !bottom left coords
                        MXI2=MXI(2,ne) !..for Xi map projection
                        IF(DOP) THEN
                          WRITE(OP_STRING,'('' nd='',I5,'' ne='',I4,'
     '                      //''' mxi1='',I3,'' mxi2='',I3,'' XID:'','
     '                      //'2E12.3)')
     '                      nd,ne,MXI1,MXI2,XID(1,nd),XID(2,nd)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                        PROJEC='XI'
                        DO ni=1,3
                          POINT(ni)=XID(ni,nd)
                        ENDDO
                        CALL POLYMARKER(INDEX,iw,1,POINT,ERROR,*9999)
                      ELSE
                        NDD=NDD+1
C ??? cpb 25/1/95 Uncommenting the zdd lines because it doesn't work
C ??? otherwise. I'm not sure why they were commented in the first
C ??? place

C ??? cpb 25/1/95 This code should be generalised with NJ_LOC

                        DO nj=1,NJT
                          IF(TYPE(1:17).EQ.'GEOMETRY_DEFORMED') THEN
                            Z(nj)=ZD(NJT+nj,nd)
                            ZDD(nj,NDD)=ZD(NJT+nj,nd)
                          ELSE
                            Z(nj)=ZD(nj,nd)
                            ZDD(nj,NDD)=ZD(nj,nd)
                          ENDIF
                        ENDDO
                      ENDIF
                    ENDIF
                  ENDIF

                ELSE IF(.NOT.STATIC) THEN
                  IF(DABS(XID(NXIDEF,nd)-TIME).LT.1.0D-06.
     '                   AND.(NXIDEF.LE.3.OR.(XID(3,nd).GE.XI3MIN.
     '                   AND.XID(3,nd).LE.XI3MAX))) THEN
                    NDD=NDD+1
                    DO nj=1,NJT
                      Z(nj)=ZD(nj,nd)
                      ZDD(nj,NDD)=ZD(nj,nd)
                    ENDDO
                  ENDIF
                ENDIF

              ELSE IF(.NOT.XI_DEFINED) THEN
                NDD=NDD+1
                DO nj=1,NJT
                  Z(nj)=ZD(nj,nd)
                  ZDD(nj,NDD)=ZD(nj,nd)
                ENDDO
              ENDIF

            ELSE !history
              CALL ASSERT(NPPF.GT.0,'>>NPPF is zero',ERROR,*9999)
              T1=XID(NXIDEF,NDT) !take final time to be that of last data point
              IF(MOD(nd,NPPF).EQ.MOD(NDHIST,NPPF)) THEN
                NDD=NDD+1
                IF(DOP) THEN
                  WRITE(OP_STRING,
     '              '('' nd,ld,nbj,nit='',4I4)') nd,LD(nd),
     '              NBJ(NJHIST,LD(nd)),NIT(NBJ(NJHIST,LD(nd)))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                Z(1)=XID(NXIDEF,nd)/T1 !normalised time
c                ZDD(1,NDD)=XID(NXIDEF,nd)/T1 !normalised time
                Z(2)=ZD(NJHIST,nd)
c                ZDD(2,NDD)=ZD(NJHIST,nd)     !position
c                IF(DOP) THEN
c                  WRITE(OP_STRING,*)' zdd(1',NDD,')=',zdd(1,ndd)
c                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c                  WRITE(OP_STRING,*)' zdd(2',NDD,')=',zdd(2,ndd)
c                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c                ENDIF
              ENDIF
            ENDIF
            IF(NJT.EQ.2) THEN
C              IF(iw.NE.10.and.iw.NE.5.and.iw.NE.6) THEN
C                CALL STRING_TRIM(TITLE,IBEG,IEND)
C                XXMIN=DBLE(XMIN)+0.01D0*DBLE(XMAX-XMIN)
C                XXMAX=DBLE(XMAX)-0.01D0*DBLE(XMAX-XMIN)
C                YYMIN=DBLE(YMIN)+0.01D0*DBLE(YMAX-YMIN)
C                YYMAX=DBLE(YMAX)-0.01D0*DBLE(YMAX-YMIN)
C                XXDIF=XXMAX-XXMIN
C                YYDIF=YYMAX-YYMIN
C                TEXT_POINT(1)=XXMAX-0.01D0*XXDIF
C                TEXT_POINT(2)=YYMIN+0.01D0*YYDIF
C CPB   28/8/92 COMMENTING OUT TEXT CALLS BECAUSE OF ERROR WITH POSTSCRIPT
C                CALL TEXT(1,iw,TITLE(IBEG:IEND),TEXT_POINT,ERROR,*9999)
C                CALL POLYMARKER(INDEX,iw,NDD,ZDD,ERROR,*9999)
C                DO nj=1,NJT
C                  POINT(nj)=ZDD(nj,nd)
C                ENDDO
C                CALL POLYMARKER(INDEX,iw,1,POINT,ERROR,*9999)
C              ELSE IF(iw.EQ.5.OR.iw.EQ.6) THEN
CC                CALL POLYMARKER(INDEX,iw,NDD,ZDD,ERROR,*9999)
C                DO nj=1,NJT
C                  POINT(nj)=ZDD(nj,nd)
C                ENDDO
C                CALL POLYMARKER(INDEX,iw,1,POINT,ERROR,*9999)
C              ELSE IF(iw.EQ.10) THEN
C                CALL POLYMARKER(INDEX,iw,NDD,ZDD,ERROR,*9999)
              IF(DRAW) THEN
                DO nj=1,NJT
                  POINT(nj)=ZDD(nj,ndd)
                ENDDO
                CALL POLYMARKER(INDEX,iw,1,POINT,ERROR,*9999)
              ENDIF
            ELSE IF(NJT.EQ.3) THEN
              IF(iw.EQ.4.AND.MAP_PROJEC(1:2).EQ.'XI') THEN
              ELSE
C                XXMIN=DBLE(XMIN)+0.01D0*DBLE(XMAX-XMIN)
C                XXMAX=DBLE(XMAX)-0.01D0*DBLE(XMAX-XMIN)
C                YYMIN=DBLE(YMIN)+0.01D0*DBLE(YMAX-YMIN)
C                YYMAX=DBLE(YMAX)-0.01D0*DBLE(YMAX-YMIN)
C                ZZMIN=DBLE(ZMIN)+0.01D0*DBLE(ZMAX-ZMIN)
C                ZZMAX=DBLE(ZMAX)-0.01D0*DBLE(ZMAX-ZMIN)
C                XXDIF=XXMAX-XXMIN
C                YYDIF=YYMAX-YYMIN
C                TEXT_POINT(1)=XXMAX-0.01D0*XXDIF
C                TEXT_POINT(2)=YYMIN-0.01D0*YYDIF
C                TEXT_POINT(3)=ZZMIN+0.01D0*YYDIF
                IF(iw.EQ.4) THEN
                  PROJEC='RECTANGULAR' !text on iw=4 in rc coords
C                  TEXT_POINT(1)=-0.85D0
C                  TEXT_POINT(2)=0.90D0
c                  DO nd=1,NDD
                    CALL ZX(ITYP10(1),Z,Z) !transform to polar coords
c                  ENDDO
                ENDIF
C CPB   28/8/92 COMMENTING OUT TEXT CALLS BECAUSE OF ERROR WITH POSTSCRIPT
C                IF(iw.NE.3) CALL TEXT(1,iw,TITLE(IBEG:IEND),TEXT_POINT,
C       '          ERROR,*9999)
                PROJEC=MAP_PROJEC    !markers on iw=4
                CALL POLYMARKER(INDEX,iw,1,Z,ERROR,*9999)
              ENDIF
            ENDIF
          ENDDO
C GBS 21/10/94
c          IF(DOP) THEN
c            WRITE(OP_STRING,*) ' NDD,NELIST(0),NELIST=',
c     '        NDD,NELIST(0),(NELIST(NEE),NEE=1,NELIST(0))
c            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c          ENDIF

        ELSE IF(TYPE(1:5).EQ.'FIBRE') THEN
          DO nd=1,NDT
            ne=LD(nd)
            IF(INLIST(ne,NELIST(1),NELIST(0),N)) THEN
              IF(DOP) THEN
                WRITE(OP_STRING,*) ' ne=',ne,' nelist(1)=',NELIST(1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(NIT(NBJ(1,ne)).LT.3) THEN
                CALL FIBRE1(INDEX,IBT,IDO,INP,iw,
     '            NBJ,ne,NKJE,NPF,NPNE,NRE,NVJE,
     '            ZD(NJT+1,nd),0.1D0,SE,XA,
     '            XE,XID(1,nd),XP,ERROR,*9999)
              ELSE IF(XID(3,nd).GE.XI3MIN.AND.XID(3,nd).LE.XI3MAX) THEN
                CALL FIBRE1(INDEX,IBT,IDO,INP,iw,
     '            NBJ,ne,NKJE,NPF,NPNE,NRE,NVJE,
     '            ZD(NJT+1,nd),0.1D0,SE,XA,
     '            XE,XID(1,nd),XP,ERROR,*9999)
              ENDIF
            ENDIF
          ENDDO

        ELSE IF(TYPE(1:6).EQ.'SHEETS') THEN
          IF(ITYP10(1).EQ.2) THEN      !cyl. polar
            nj=2
          ELSE IF(ITYP10(1).EQ.4) THEN !prolate
            nj=3
          ENDIF
          THETA1=SHEET_THETA-0.5D0*SHEET_RANGE
          THETA2=SHEET_THETA+0.5D0*SHEET_RANGE
          IF(DOP) THEN
            WRITE(OP_STRING,'('' SHEET_THETA='',E12.3,'' THETA1='','
     '        //'E12.3,'' THETA2='',E12.3)') SHEET_THETA,THETA1,THETA2
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(.NOT.CALC_SHEET) THEN !sheet angle correct.ns not been made
            DO nd=1,NDT
              CALL ZX(ITYP10(1),ZD(1,nd),XD)
              IF(XD(nj).GE.THETA1.AND.XD(nj).LE.THETA2) THEN
                CALL SHEET2(INDEX,iw,ZD(1,nd),ERROR,*9999)
              ENDIF
            ENDDO
          ELSE IF(CALC_SHEET) THEN !sheet angle correctns have been made
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),NRE(ne),
     '          NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
              nb=NBJ(nj,ne)
              DO nde=1,NDLT(ne)
                nd=NDDL(ne,nde)
                CALL ZX(ITYP10(1),ZD(1,nd),XD)
                XI(1)=XID(1,nd)
                XI(2)=XID(2,nd)
                XI(3)=0.0D0
                THETA_ND=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XI,XE(1,nj))
                IF(THETA_ND.GE.THETA1.AND.THETA_ND.LE.THETA2) THEN
                  CALL ASSERT(.FALSE.,'>> Old code - gone to archive',
     '              ERROR,*9999)
C                  CALL SHEET3(INDEX,IBT,IDO,INP,iw,NAN,NBJ,ne,nr,
C     '              XE,XG,XID(1,nd),
C     '              ZD(1,nd),ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
          ENDIF

        ELSE IF(TYPE(1:9).EQ.'PARAMETER') THEN
          NOPARA=IFROMC(TYPE(10:11))
          IF(DOP) THEN
            WRITE(OP_STRING,'('' parameter number = '',I2)') NOPARA
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            IF(ITYP1(nr,nx).EQ.3) THEN
              CALL CPCG(1,NBJ(1,ne),NPNE(1,1,ne),nr,nx,
     '          CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
            ELSE
              CALL CPCG(NW(ne,1,nx),NBJ(1,ne),NPNE(1,1,ne),nr,nx,
     '          CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Element '',I3)') ne
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                WRITE(OP_STRING,'('' XE(ns,'',I1,''): '',8E11.3)')
     '            nj,(XE(ns,nj),ns=1,NST(NBJ(nj,ne)))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
              WRITE(OP_STRING,'('' CE('',I2,'',ne)= '',E11.3)')
     '          NOPARA,CE(NOPARA,ne)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO ng=1,NGT(NBJ(1,ne))
              CALL XEXG(NBJ(1,ne),ng,NRE(ne),PG,
     '          XE,XG,ERROR,*9999)
              X(1)=XG(1,1)
              X(2)=XG(2,1)+SCALE*CG(NOPARA,ng)
              X(3)=XG(3,1)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' X: '',3E11.3)') X(1),X(2),X(3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL ZZ(X,Z,TRANS)
              CALL POLYMARKER(INDEX,iw,1,Z,ERROR,*9999)
            ENDDO
          ENDDO

        ELSE IF(TYPE(1:5).EQ.'SLICE') THEN
          IF(iw.LE.3) THEN
            DO nd=1,NDT
              IF(Z(iw).GE.XI3MIN.AND.Z(iw).LE.XI3MAX) THEN
                DO nj=1,NJT
                  POINT(nj)=ZD(nj,nd)
                ENDDO
                CALL POLYMARKER(INDEX,1,1,POINT,ERROR,*9999)
              ENDIF
            ENDDO
          ENDIF
        ENDIF

        CALL CLOSE_SEGMENT(ISDATA,iw,ERROR,*9999)

      ELSE IF(TYPE(1:7).EQ.'NUMBERS'.OR.TYPE(1:6).EQ.'VALUES') THEN
        XI_DEFINED=.FALSE.
        DO nd=1,NDT
          IF(LD(nd).GT.0) XI_DEFINED=.TRUE.
        ENDDO
        IF(.NOT.XI_DEFINED) THEN !loop thru all data points
          CALL OPEN_SEGMENT(ISDANO(iw,1),ISEG,iw,'DANO',INDEX,
     '      INDEX_OLD,1,1,CSEG,ERROR,*9999)
          IF(ISEG(ISDANO(iw,1)).GT.0) THEN
            DO nd=1,NDT
              IF(TYPE(1:7).EQ.'NUMBERS') THEN
                WRITE(CHAR,'(I4)') NDP(nd)
              ELSE IF(TYPE(1:6).EQ.'VALUES') THEN
                WRITE(CHAR,'(F4.0)') ZD(NJT+1,nd)
              ENDIF
              CALL STRING_TRIM(CHAR,IBEG,IEND)
              IF(TYPE(1:6).EQ.'VALUES') THEN
                IEND=IEND-1
              ENDIF
              CALL TEXT(1,iw,CHAR(IBEG:IEND),ZD(1,nd),ERROR,*9999)
            ENDDO !nd
          ENDIF
          CALL CLOSE_SEGMENT(ISDANO(iw,1),iw,ERROR,*9999)

        ELSE IF(XI_DEFINED) THEN !loop thru data pts within ea element
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            CALL OPEN_SEGMENT(ISDANO(iw,ne),ISEG,iw,'DANO',INDEX,
     '        INDEX_OLD,ne,1,CSEG,ERROR,*9999)
            IF(ISEG(ISDANO(iw,ne)).GT.0) THEN
              DO nde=1,NDLT(ne)
                nd=NDDL(ne,nde)
                IF(TYPE(1:7).EQ.'NUMBERS') THEN
                  WRITE(CHAR,'(I4)') NDP(nd)
                ELSE IF(TYPE(1:6).EQ.'VALUES') THEN
                  WRITE(CHAR,'(F4.0)') ZD(NJT+1,nd)
                ENDIF
                CALL STRING_TRIM(CHAR,IBEG,IEND)
                IF(TYPE(1:6).EQ.'VALUES') THEN
                  IEND=IEND-1
                ENDIF
C why was this done? PJH 21Feb96
C               CALL ZX(ITYP10(1),ZD(1,nd),Z) !transform to polar coords
C               CALL TEXT(1,iw,CHAR(IBEG:IEND),Z,ERROR,*9999)
                CALL TEXT(1,iw,CHAR(IBEG:IEND),ZD(1,nd),ERROR,*9999)
              ENDDO
            ENDIF
            CALL CLOSE_SEGMENT(ISDANO(iw,ne),iw,ERROR,*9999)
          ENDDO

        ENDIF !xi_defined

      ELSE IF(TYPE(1:11).EQ.'PROJECTIONS') THEN
        DO nolist=1,NELIST(0)
          ne=NELIST(nolist)
          CLABEL(1:4)='DAPR'
          IF(STATIC) THEN
            CLABEL(43:52)='          '
            IF(NXIDEF.EQ.0) NXIDEF=NIT(NBJ(1,ne))
            IF(DOP) THEN
              WRITE(OP_STRING,*) 'nxidef=',nxidef
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            WRITE(CLABEL(43:52),'(E10.3)') TIME
            IF(NXIDEF.EQ.0) NXIDEF=NIT(NBJ(1,ne))+1
          ENDIF
          CALL OPEN_SEGMENT(ISDAPR(iw,ne),ISEG,iw,CLABEL,INDEX,
     '      INDEX_OLD,ne,1,CSEG,ERROR,*9999)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' ne='',I5,'' ndlt(ne)='',I5,'
     '        //''' nddl(ne,nde):'',/(20I5))')
     '        ne,NDLT(ne),(NDDL(ne,nde),nde=1,NDLT(ne))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL ZDZDL(0,NBH(1,1,ne),NBJ(1,ne),NDDL,NDLT(ne),ne,1,WD,WDL,
     '      XID,XIDL,ZD,ZDL,ERROR,*9999)
          IF(.NOT.STATIC) THEN !motion field in ZE
            CALL ZPZE(NBH(1,1,ne),nc,NJ_LOC(NJL_GEOM,0,nr),
     '        NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '        NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '        ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
            ENDIF
          NDD=0
          DO nd=1,NDLT(ne)
            IF(STATIC) THEN
              DRAW=.TRUE.
            ELSE
              DRAW=DABS(XIDL(NXIDEF,nd)-TIME).LT.1.0D-06
            ENDIF
            IF(DRAW) THEN
              IF(NXIDEF.LT.3) THEN
                DRAW=.TRUE.
              ELSE
                DRAW=XID(3,nd).GE.XI3MIN.AND.XID(3,nd).LE.XI3MAX
              ENDIF
              IF(DRAW) THEN
                NDD=NDD+1
                DO nj=1,NJT
                  ZDL(nj,NDD)=ZDL(nj,nd)
                ENDDO
                DO ni=1,NXIDEF
                  XIDL(ni,NDD)=XIDL(ni,nd)
                ENDDO
              ENDIF
            ENDIF
          ENDDO

          CALL LDAT(IBT,IDO,INP,iw,NBH(1,1,ne),NBJ(1,ne),NDD,nr,
     '      STATIC,XE,XIDL,ZDL,ZDLMIN,ZE,ERROR,*9999)

          CALL CLOSE_SEGMENT(ISDAPR(iw,ne),iw,ERROR,*9999)
        ENDDO

      ELSE IF(TYPE(1:5).EQ.'TRACE') THEN
        DO nolist=1,NELIST(0)
          ne=NELIST(nolist)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' ne='',I4)') ne
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CLABEL(1:4)='DATR'
          IF(STATIC) THEN
            CLABEL(43:52)='          '
          ELSE
            WRITE(CLABEL(43:52),'(E10.3)') TIME
          ENDIF
          CALL OPEN_SEGMENT(ISDATR(iw,ne),ISEG,iw,CLABEL,INDEX,
     '      INDEX_OLD,ne,1,CSEG,ERROR,*9999)
          IF(ISEG(ISDATR(iw,ne)).GT.0) THEN
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            CALL ZDZDL(0,NBH(1,1,ne),NBJ(1,ne),NDDL,NDLT(ne),ne,1,WD,
     '        WDL,XID,XIDL,ZD,ZDL,ERROR,*9999)
            NTNDE=0
            ni=NIT(NBJ(1,ne))
            IF(.NOT.STATIC) THEN !motion field in ZE
              NITB=NIT(NBH(NH_LOC(1,nx),1,ne))
              CALL ZPZE(NBH(1,1,ne),nc,NJ_LOC(NJL_GEOM,0,nr),
     '          NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '          NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '          ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
            ELSE
              NITB=NIT(NBJ(1,ne))
            ENDIF
            DO nde=1,NDLT(ne)
              nd=NDDL(ne,nde)
              IF(STATIC) THEN
                DRAW=.TRUE.
              ELSE
                DRAW=DABS(XIDL(NITB,nd)-TIME).LT.1.0D-06
              ENDIF
              IF(DRAW) THEN
              IF(NITB.LT.3) THEN
                DRAW=.TRUE.
              ELSE
                DRAW=XID(3,nd).GE.XI3MIN.AND.XID(3,nd).LE.XI3MAX
              ENDIF
              IF(DRAW) THEN
                NTNDE=NTNDE+1
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' nde='',I4,'' ntnde='',I4,'
     '              //''' XI(ni,nde): '',4E11.3)')
     '              nde,NTNDE,(XIDL(ni,nde),ni=1,NITB)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                DO nj=1,NJT
                  nb=NBJ(nj,ne)
                  X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,XIDL(1,nde),XE(1,nj))
                  Z2(nj)=0.0D0
                  IF(.NOT.STATIC) THEN !motion
                    NBO=NBH(nj,1,ne)
                    Z2(nj)=PXI(IBT(1,1,NBO),IDO(1,1,0,NBO),INP(1,1,NBO),
     '                NBO,1,XIDL(1,nde),ZE(1,nj))
                  ENDIF
                ENDDO
                IF(iw.EQ.4) THEN !leave trace coords as ITYP10(1)
                  DO nj=1,NJT
                    Z(nj)=X(nj)
                  ENDDO
                ELSE !convert to rect. cart.
                  CALL XZ(ITYP10(1),X,Z)
                ENDIF
                DO nj=1,NJT
                  POINT(nj)=Z(nj)
                  IF(KTYP58(nr).EQ.2) THEN
                    POINT(nj)=POINT(nj)+Z2(nj)
                  ENDIF
                ENDDO
                CALL POLYMARKER(INDEX,iw,1,POINT,ERROR,*9999)
              ENDIF
              ENDIF
            ENDDO !nde
          ENDIF !iseg>0
          CALL CLOSE_SEGMENT(ISDATR(iw,ne),iw,ERROR,*9999)
        ENDDO !nolist
      ENDIF !type

      CALL EXITS('SGDATA')
      RETURN
 9999 CALL ERRORS('SGDATA',ERROR)
      CALL EXITS('SGDATA')
      RETURN 1
      END


