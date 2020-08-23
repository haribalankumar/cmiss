      SUBROUTINE CADATA(ISDANO,ISDAPR,ISDATA,ISDATR,ISEG,
     '  LD,NDDL,NDLT,NDP,NEELEM,NELIST,
     '  CSEG,EDD,SQ,STRING,WD,XID,ZD,ERROR,*)

C#### Subroutine: CADATA
C###  Description:
C###    CADATA cancels data segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISDANO(NWM,NEM),ISDAPR(NWM,NEM),ISDATA(NWM,NGRSEGM),
     '  ISDATR(NWM,NEM),ISEG(*),LD(NDM),NDDL(NEM,NDEM),
     '  NDLT(NEM),NDP(NDM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM)
      REAL*8 EDD(NDM),SQ(NDM),WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,
     '  INSTAT,ISBOX,iunused,iw,IWK(6),N1DATA,N2DATA,N3CO,nd,ndd1,
     '  ndd2,NDDTOT,nde,ne,ni,nj,nodata,noelem,noiw,nr,NTIW,
     '  nolist
      REAL*8 XBOX,XWC,YBOX,YWC,ZWC,RFROMC,EDDMIN
      CHARACTER TYPE*11
      LOGICAL ABBREV,CBBREV,CONTINUE,MOUSE,SEGME,ALL_ELEM

      CALL ENTERS('CADATA',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel data;s
C###  Description:
C###    Cancel specified type of data segment on specified workstation.
C###  Parameter:      <(positions/numbers/projections/trace)[position]>
C###    Specify the type of data segment to be cancelled.
C###  Parameter:      <(all/last)[all]>
C###    Cancel either all data segments or just the last segment
C###    to be drawn.
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the window number on which to cancel the data.
C###    The default is to cancel the data on all windows.

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15)
     '    //'<(positions/numbers/projections/trace)[positions]>'
        OP_STRING(3)=BLANK(1:15) //'<(all/last)[all]>'
        OP_STRING(4)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM cancel data;m
C###  Description:
C###    Cancel specified data. This removes data from the
C###    problem.
C###  Parameter:      <on WS#[1]>
C###    Specify the window number on which to cancel the data.
C###    The default is to cancel the data on window 1.
C###  Parameter:      <xi>
C###    Remove the xi calculation of specified data sets.
C###  Parameter:      <in (all/ELEMENT#s)[all]>
C###    Specify which elements to cancel segments in. The default is
C###    to cancel segments in all elements.
C###  Parameter:      <errors <greater VALUE#[0.0]>>
C###    Remove the data point errors for specified data sets.
C###    Specifying 'greater' means only remove those that are greater
C###    than VALUE.

        OP_STRING(1)=STRING(1:IEND) //';m'
        OP_STRING(2)=BLANK(1:15) //'<on WS#[1]>'
        OP_STRING(3)=BLANK(1:15) //'<xi>'
        OP_STRING(4)=BLANK(1:15) //'<in (all/ELEMENT#s)[all]>'
        OP_STRING(5)=BLANK(1:15) //'<errors <greater VALUE#[0.0]>>'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CADATA',ERROR,*9999)
      ELSE
C        nx=1 !temporary

        CALL CHECKQ(' SM',noco,1,CO,COQU,STRING,*1)
        SEGME=.FALSE.
        MOUSE=.FALSE.
        IF(ABBREV(COQU(noco,1),'S',1)) THEN
          SEGME=.TRUE.
        ELSE IF(ABBREV(COQU(noco,1),'M',1)) THEN
          MOUSE=.TRUE.
        ENDIF
        IF(SEGME.OR.MOUSE) THEN
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        ENDIF

        IF(SEGME) THEN
          IF(CBBREV(CO,'POSITIONS',2,noco+1,NTCO,N3CO)) THEN
            TYPE='POSITION'
          ELSE IF(CBBREV(CO,'NUMBERS',2,noco+1,NTCO,N3CO)) THEN
            TYPE='NUMBERS'
          ELSE IF(CBBREV(CO,'PROJECTIONS',2,noco+1,NTCO,N3CO)) THEN
            TYPE='PROJECTIONS'
          ELSE IF(CBBREV(CO,'TRACE',2,noco+1,NTCO,N3CO)) THEN
            TYPE='TRACE'
          ELSE
            TYPE='POSITION'
          ENDIF
          IF(CBBREV(CO,'ALL',1,noco+1,NTCO,N3CO)) THEN
            N1DATA=1
            N2DATA=NRM
            NTDATA=0
          ELSE IF(CBBREV(CO,'LAST',1,noco+1,NTCO,N3CO)) THEN
            N1DATA=NTDATA
            N2DATA=NTDATA
            NTDATA=NTDATA-1
          ELSE
            N1DATA=1
            N2DATA=NRM
            NTDATA=0
          ENDIF

        ELSE IF(MOUSE) THEN
          NELIST(0)=0
          DO nr=1,NRT
            DO noelem=NELIST(0)+1,NELIST(0)+NEELEM(0,nr)
              NELIST(noelem)=NEELEM(noelem,nr)
            ENDDO
            NELIST(0)=NELIST(0)+NEELEM(0,nr)
          ENDDO

        ELSE
          IF(CBBREV(CO,'XI',1,noco+1,NTCO,N3CO)) THEN
            TYPE='XI'
            IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),NEM,NELIST(0),NELIST(1),ERROR,
     '          *9999)
              ALL_ELEM=.FALSE.
            ELSE
              ALL_ELEM=.TRUE.
            ENDIF
          ELSE
            TYPE='DATA'
          ENDIF
        ENDIF

        IF(SEGME) THEN
          IF(TYPE(1:8).EQ.'POSITION') THEN
            DO noiw=1,NTIW
              iw=IWK(noiw)
              IF(IWKS(iw).GT.0) THEN
                CALL ACWK(iw,1,ERROR,*9999)
                DO nodata=N1DATA,N2DATA
                  IF(ISDATA(iw,nodata).GT.0) THEN
                    CALL DELETE_SEGMENT(ISDATA(iw,nodata),ISEG,iw,
     '                ERROR,*9999)
                  ENDIF
                ENDDO
                CALL DAWK(iw,1,ERROR,*9999)
              ENDIF
            ENDDO
          ELSE IF(TYPE(1:7).EQ.'NUMBERS') THEN
            DO noiw=1,NTIW
              iw=IWK(noiw)
              IF(IWKS(iw).GT.0) THEN
                CALL ACWK(iw,1,ERROR,*9999)
                DO nr=1,NRT
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    IF(ISDANO(iw,ne).GT.0) THEN
                      CALL DELETE_SEGMENT(ISDANO(iw,ne),ISEG,iw,ERROR,
     '                  *9999)
                    ENDIF
                  ENDDO
                ENDDO
                CALL DAWK(iw,1,ERROR,*9999)
              ENDIF
            ENDDO
          ELSE IF(TYPE(1:11).EQ.'PROJECTIONS') THEN
            DO noiw=1,NTIW
              iw=IWK(noiw)
              IF(IWKS(iw).GT.0) THEN
                CALL ACWK(iw,1,ERROR,*9999)
                DO nr=1,NRT
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    IF(ISDAPR(iw,ne).GT.0) THEN
                      CALL DELETE_SEGMENT(ISDAPR(iw,ne),ISEG,iw,ERROR,
     '                  *9999)
                    ENDIF
                  ENDDO
                ENDDO
                CALL DAWK(iw,1,ERROR,*9999)
              ENDIF
            ENDDO
          ELSE IF(TYPE(1:5).EQ.'TRACE') THEN
            DO noiw=1,NTIW
              iw=IWK(noiw)
              IF(IWKS(iw).GT.0) THEN
                CALL ACWK(iw,1,ERROR,*9999)
                DO nr=1,NRT
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    IF(ISDATR(iw,ne).GT.0) THEN
                      CALL DELETE_SEGMENT(ISDATR(iw,ne),ISEG,iw,ERROR,
     '                  *9999)
                    ENDIF
                  ENDDO
                ENDDO
                CALL DAWK(iw,1,ERROR,*9999)
              ENDIF
            ENDDO
          ENDIF

        ELSE IF(MOUSE) THEN
          iw=IWK(1)
          WRITE(OP_STRING,'('' >>Cancel data on '',I1)') iw
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL ACWK(iw,0,ERROR,*9999)
          CALL LOCATOR(INSTAT,0.0D0,XWC,0.0D0,YWC,
     '      ERROR,*9999)
C ***     Use choice in event mode, device=4 (mouse buttons), echo off
c          OPTION(1)=' '
c          OPTION(2)=' '
c          OPTION(3)=' '
c         CALL CHOICE('Mouse Buttons',4,1,INS,iw,'EVENT',3,3,NOCH,noco,
c    '      1,CO,OPTION,STRING,0.0,0.0,ERROR,*9999)
          XBOX=0.01D0*DBLE(DIAG) !is half x-dimension of box
          YBOX=0.01D0*DBLE(DIAG) !is half y-dimension of box
          ISBOX=0
          CONTINUE=.TRUE.
          DO WHILE(CONTINUE)
C            BUTTON=0
            CALL LOCATOR(INSTAT,0.0D0,XWC,0.0D0,YWC,
     '        ERROR,*9999)
            IF(ISBOX.GT.0) THEN
              CALL DELETE_SEGMENT(ISBOX,ISEG,iw,ERROR,*9999)
              NTSG=NTSG-1
            ENDIF
            CALL OPEN_SEGMENT(ISBOX,ISEG,iw,'Box',INDEX,iunused,1,1,
     '        CSEG,ERROR,*9999)
C LKC 2-MAY-1998 init ZWC (not used)
            ZWC=0.D0
            CALL DBOX(iw,XBOX,YBOX,XWC,YWC,ZWC,ERROR,*9999)
            CALL CLOSE_SEGMENT(ISBOX,iw,ERROR,*9999)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' XWC='',E12.3,'' YWC='',E12.3)')
     '          XWC,YWC
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
c           CALL EVENT(ID_WS,ID_DEVICE,ID_STATUS,CLASS,IDATA,R4DATA,
c    '        SDATA,ERROR,*9999)
c            IF(DOP) THEN
c              WRITE(OP_STRING,'('' Input_class='',A8,'
c     '          //''' ID_device='',I2)') CLASS,ID_DEVICE
c              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c            ENDIF
c            IF(CLASS(1:6).EQ.'CHOICE')THEN
c              BUTTON=IDATA(1)
c              IF(DOP) THEN
c                 write(op_string,*) 'button=',button
c                 CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c              ENDIF
c              IF(BUTTON.EQ.1) THEN
c                X1=XWC-XBOX !is x-coord of bottom of box
c                X2=XWC+XBOX !is x-coord of top of box
c                Y1=YWC-YBOX !is y-coord of LH edge of box
c                Y2=YWC+YBOX !is y-coord of RH edge of box
c                nd=0
c                DO WHILE(nd.LE.NDT)
c                  nd=nd+1
c                  IF(ZD(1,nd).GE.X1.AND.ZD(1,nd).LE.X2.AND
c     '              .ZD(2,nd).GE.Y1.AND.ZD(2,nd).LE.Y2) THEN
c                    DO nd1=nd,NDT-1
c                      DO nj=1,NJM
c                        ZD(nj,nd1)=ZD(nj,nd1+1)
c                        WD(nj,nd1)=WD(nj,nd1+1)
c                      ENDDO
c                    ENDDO
c                    NDT=NDT-1
c                    WRITE(OP_STRING,
c     '                '('' Data point nd='',I6,'' deleted'')') nd
c                     CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c                    TYPE='GEOMETRY'
c                    ZDLMIN=0.0
c                    CALL SGDATA(1,IBT,IDO,INP,ISDANO,ISDAPR,
c     '                ISDATA(iw,NTDATA),ISDATR,ISEG,iw,LD,MXI,
c     '                NBJ,NBH,NDDL,NDLT,NDP,NEELEM,NELIST,
c     '                NKE,NPF,NPNE,NRE,NVHE,NVHP,NVJE,NVJP,NW,
c     '                CE(1,1,nx),CG,CP(1,1,nx),
c     '                CSEG,PG,SE,.TRUE.,' ',TYPE,VE,WD,WDL,
c     '                XA,XE,XG,XID,XIDL,XP,ZA,ZD,ZDD,ZDL,ZDLMIN,ZE,ZP,
c     '                ERROR,*9999)
c                  ENDIF
c                ENDDO
c              ELSE IF(BUTTON.EQ.2) THEN
cc               CALL INPUT_MODE(iw,1,'LOCATOR','REQUEST',ERROR,*9999)
cc               CALL INPUT_MODE(iw,4,'CHOICE','REQUEST',ERROR,*9999)
c                IF(ISBOX.GT.0) THEN
c                  CALL DELETE_SEGMENT(ISBOX,ISEG,iw,ERROR,*9999)
c                  NTSG=NTSG-1
c                ENDIF
c                CONTINUE=.FALSE.
c              ENDIF
c            ENDIF
          ENDDO
          CALL DAWK(iw,0,ERROR,*9999)

        ELSE
C CPB 13/1/92 added cancel data error greater #
          IF(CBBREV(CO,'ERRORS',1,noco+1,NTCO,N3CO)) THEN
            IF(CBBREV(CO,'GREATER',1,noco+1,NTCO,N3CO)) THEN
              EDDMIN=RFROMC(CO(N3CO+1))
            ELSE
              EDDMIN=0.0D0
            ENDIF
            DO nd=1,NDT
              IF(ABS(EDD(nd)).GE.EDDMIN) THEN
                NDDTOT=0
                DO ndd1=1,NDLT(LD(nd))
                  ndd2=NDDL(LD(nd),ndd1)
                  IF(ndd2.NE.nd) THEN
                    NDDTOT=NDDTOT+1
                    NDDL(LD(nd),NDDTOT)=ndd2
                  ENDIF
                ENDDO
                NDLT(LD(nd))=NDDTOT
C GRC 25-JUN-2001 Do not clear NDP here. This allows us to write data
C files that have been host-mesh warped and still relate their NDP
C numbers back to the source data. Even if a data point is not in an
C element leave NDP as it was -- we want those points to remain
C unchanged in number and position when we warp those data points that
C are in elements. Clearing NDP should never be the default behaviour!
C                NDP(nd)=0
                LD(nd)=0
                EDD(nd)=0.0D0
                SQ(nd)=0.0D0
                DO ni=1,NIM
                  XID(ni,nd)=0.0D0
                ENDDO

C LKC 25-JUL-2001 May as well set the weights of cancelled points to 0
                DO nj=1,NJT
                  WD(nj,nd)=0.0D0
                ENDDO

              ENDIF
            ENDDO
          ELSE
            IF(ALL_ELEM) THEN
              DO nd=1,NDM
                NDP(nd)=0
                LD(nd)=0
                EDD(nd)=0.0D0
                SQ(nd)=0.0D0
                DO ni=1,NIM
                  XID(ni,nd)=0.0D0
                ENDDO
              ENDDO
              DO ne=1,NEM
                NDLT(ne)=0
                DO nde=1,NDEM
                  NDDL(ne,nde)=0
                ENDDO
              ENDDO
              CALC_XI=.FALSE.
C CPB 11/2/92 added cancel data xi in element list
            ELSE
              DO nd=1,NDT
                DO nolist=1,NELIST(0)
                  ne=NELIST(nolist)
                  IF(LD(nd).eq.ne) THEN
                    NDP(nd)=0
                    LD(nd)=0
                    EDD(nd)=0.0D0
                    SQ(nd)=0.0D0
                    DO ni=1,NIM
                      XID(ni,nd)=0.0D0
                    ENDDO
                    NDLT(ne)=0
                    DO nde=1,NDEM
                      NDDL(ne,nde)=0
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
            IF(TYPE(1:4).EQ.'DATA') THEN
              DO nd=1,NDM
                DO nj=1,NJM
                  ZD(nj,nd)=0.0D0
                  WD(nj,nd)=0.0D0
                ENDDO
              ENDDO
              NDT=0
              CALL_DATA_FIBRE=.FALSE.
              CALL_DATA_FIELD=.FALSE.
              CALL_DATA_SHEET=.FALSE.
              CALC_SHEET=.FALSE.
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('CADATA')
      RETURN
 9999 CALL ERRORS('CADATA',ERROR)
      CALL EXITS('CADATA')
      RETURN 1
      END


