      SUBROUTINE CHDATA(IBT,IDO,INP,LD,LN,NBJ,NDDL,
     '  NDLT,NEELEM,NELIST,NKJE,NPF,NPNE,NRLIST,NVJE,
     '  NXI,SE,WD,WDL,XA,XE,XID,XP,ZD,STRING,ERROR,*)

C#### Subroutine: CHDATA
C###  Description:
C###    CHDATA changes data locations.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'dtran00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),
     '  LD(NDM),LN(0:NEM),NBJ(NJM,NEM),NDDL(NEM,NDEM),NDLT(NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),WD(NJM,NDM),WDL(NHM,NDEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,l,n,N3CO,nb,nd,nde,ne,
     '  NEXI(2),nj,nolist,nr,NTRL
      REAL*8 a,ANGLE1,ANGLE2,AXIS(3),AMOUNT,b,CENTRE(3),
     '  DOT_PROD,dxdxi1,dxdxi2,dydxi1,dydxi2,
     '  dzdxi1,dzdxi2,g1(3),g2(3),INCREMENT,
     '  length,MUSP(0:5),POSN,PXI,RAD1,
     '  RAD2,RAD3,RFROMC,
     '  RL1(3),RL2(3),RL3(3),SMU,SUM(3),theta,
     '  TRANS(3,4),WEIGHT,XI3_1,XI3_2,
     '  XPD(3),XX(3),zd1,zd2,ZT(3)
      CHARACTER FILE*100,TYPE*12
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,NEWXI,PROJECT,RAXIS,ROTATE,
     '  RSHIFT,TRANSL

      DATA MUSP/ 0.d0,0.7854d0,0.4363d0,0.3491d0,0.2618d0,0.2618d0 /
      DATA NEXI/8,5/

      CALL ENTERS('CHDATA',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM change<;project><;newxi> data
C###  Parameter:      <translate by DX#,DY#,DZ#>
C###    Translates the data points in the cartesian coordinate system
C###  Parameter:      <scale by SX#,SY#,SZ#>
C###    Scale the coordinates of the data points in cartesian
C###    coordinates
C###  Parameter:      <xform FILENAME>
C###  Parameter:      <rotate by ANGLE axis BX#,BY#,BZ# about AX#,AY#,AZ>
C###    Rotate the data points by an angle around an axis which runs
C###    through a point defined by A
C###  Parameter:      <rotate by TX#,TY#,TZ# about AX#,AY#,AZ#>
C###    Rotate the data by a specified angles TX,TY,TZ in turn about
C###    the respective axes which run though the point AX,AY,AZ.
C###  Parameter:      <centre at CX#[0.0],CY#[0.0],CZ#[0.0]>
C###    Defines the centre of the scaling
C###  Description:
C###    Change the data set coordinates.

        OP_STRING(1)=STRING(1:IEND)//'<;project><;newxi>'
        OP_STRING(2)=BLANK(1:15)//'<translate by DX#,DY#,DZ#>'
        OP_STRING(3)=BLANK(1:15)//'<scale by SX#,SY#,SZ#>'
        OP_STRING(4)=BLANK(1:15)//'<xform FILENAME>'
        OP_STRING(5)=BLANK(1:15)
     '    //'<rotate by TX#,TY#,TZ# about AX#,AY#,AZ#>'
        OP_STRING(6)=BLANK(1:15)//'<axis BX#,BY#,BZ#>'
        OP_STRING(7)=BLANK(1:15)
     '    //'<centre at CX#[0.0],CY#[0.0],CZ#[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change data fibres DEGREES#{degrees} DEGREES#{degrees}
C###  Parameter:      <region (#s/all)[1]>
C###    Defines the region of the fibre adjustment
C###  Parameter:      <element (#s/all)[all]>
C###    Defines the elements to change
C###  Parameter:      <from XI_3#[0.0]{>=0.0,<=1.0}>
C###  Parameter:      <to XI_3#[1.0]{>=0.0,<=1.0}>
C###    Define the Xi_3 limits for which the fibres are to be changed
C###  Description:
C###    Change the fibre data.


        OP_STRING(1)=STRING(1:IEND)//' fibres DEGREES#{degrees}'
     '    //' DEGREES#{degrees}'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<from XI_3#[0.0]{>=0.0,<=1.0}>'
        OP_STRING(5)=BLANK(1:15)//'<to XI_3#[1.0]{>=0.0,<=1.0}>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change data fibres surface_norm
C###  Description:
C###    Changes the fibre angles to correct for difference between
C###    surface normal and rig probe direction

        OP_STRING(1)=STRING(1:IEND)//' surface_norm'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change data coordinates cylinder to rectangular
C###  Description:
C###    Changes the data coordinates from cylindrical-polar
C###    coordinates to rectangular-cartesian

        OP_STRING(1)=STRING(1:IEND)//' coordinates '
     '   //'cylinder to rectangular'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change data fibres wrt_xi1
C###  Description:
C###    Changes the fibre angles to correct for devations
C###    from the fibre angle reference plane with xi1

        OP_STRING(1)=STRING(1:IEND)//' wrt_xi1'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change data fibres radians
C###  Description:
C###    Changes the fibre angles from degrees to radians

        OP_STRING(1)=STRING(1:IEND)//' radians'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change data weights WEIGHT#
C###  Parameter:      <region (#s/all)[1]>
C###    Defines the region of the adjustment of the weights
C###  Parameter:      <element (#s/all)[all]>
C###    Defines the elements to change
C###  Parameter:      <from XI_3#[0.0]{>=0.0,<=1.0}>
C###  Parameter:      <to XI_3#[1.0]{>=0.0,<=1.0}>
C###    Define the Xi_3 limits for which the weights are to be changed
C###  Parameter:      <(normal/stripe/contour)[normal]>
C###  Description:
C###    Change the data set weights.

        OP_STRING(1)=STRING(1:IEND)//' weights WEIGHT#'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<from XI_3#[0.0]{>=0.0,<=1.0}>'
        OP_STRING(5)=BLANK(1:15)//'<to XI_3#[1.0]{>=0.0,<=1.0}>'
        OP_STRING(6)=BLANK(1:15)//'<(normal/stripe/contour)[normal]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change data sheets surface_norm
C###  Description:
C###    Changes sheet angles to correct for difference between
C###    surface normal and measured radial direction

        OP_STRING(1)=STRING(1:IEND)//' surface_norm'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change data xi_increment <INCREMENT[0.1]>
C###  Description:
C###    Changes the xi_1 projection of data points by
C###    a specified increment

        OP_STRING(1)=STRING(1:IEND)//' xi_increment <INCREMENT[0.1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHDATA',ERROR,*9999)
      ELSE

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
        nr=NRLIST(1)

        IF(CBBREV(CO,'FIBRES',2,noco+1,NTCO,N3CO)) THEN
          TYPE='FIBRES'
        ELSE IF(CBBREV(CO,'WEIGHTS',2,noco+1,NTCO,N3CO)) THEN
          TYPE='WEIGHTS'
        ELSE IF(CBBREV(CO,'XI_INCREMENT',2,noco+1,NTCO,N3CO)) THEN
          TYPE='XI_INCREMENT'
        ELSE
          TYPE='VALUES'
        ENDIF

        IF(CBBREV(CO,'COORDINATES',4,noco+1,NTCO,N3CO)) THEN
          DO nd=1,NDT
            zd1=ZD(1,nd)
            zd2=ZD(2,nd)
            ZD(1,nd)=zd1*COS(zd2*PI/180.0d0)
            ZD(2,nd)=zd1*SIN(zd2*PI/180.0d0)
          ENDDO
        ENDIF

        IF(TYPE(1:6).EQ.'VALUES') THEN
          IF(NTCOQU(noco).GT.0) THEN
            PROJECT=ABBREV(COQU(noco,1),'PROJECT',1)
          ELSE
            PROJECT=.FALSE.
          ENDIF
          IF(NTCOQU(noco).GT.1) THEN
            NEWXI=ABBREV(COQU(noco,2),'NEWXI',1)
          ELSE
            NEWXI=.FALSE.
          ENDIF
          CALL RESET(TRANS)
          CALL RESET(DATRAN)

          TRANSL=.FALSE.
          IF(CBBREV(CO,'TRANSLATE',1,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'BY',1)) THEN
              CALL PARSRL(CO(N3CO+2),3,NTRL,ZT,ERROR,*9999)
              TRANSL=.TRUE.
              AMOUNT=0.d0
              DO nj=1,NJT
C PJH 9/7/98 wrong direction for translation   ZT(nj)=-ZT(nj)
                AMOUNT=AMOUNT+ZT(nj)**2
              ENDDO
              AMOUNT=DSQRT(AMOUNT)
              CALL SHIFT(ZT,AMOUNT,TRANS,ERROR,*9999)
              CALL SHIFT(ZT,AMOUNT,DATRAN,ERROR,*9999)
            ENDIF
          ENDIF

          IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
            TRANSL=.TRUE.
            IF(CBBREV(CO,'BY',1,noco+2,noco+3,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTRL,RL1,ERROR,*9999)
            ELSE
              RL1(1)=1.d0
              RL1(2)=1.d0
              RL1(3)=1.d0
            ENDIF
            CALL SCALE1(RL1(1),TRANS)
            CALL SCALE2(RL1(2),TRANS)
            IF(NJT.EQ.3) THEN
              CALL SCALE3(RL1(3),TRANS)
            ENDIF
          ENDIF

          IF(CBBREV(CO,'XFORM',1,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            FILE=CO(N3CO+1)(IBEG:IEND)
            CALL OPENF(IFILE,'DISK',FILE(1:IEND-IBEG+1)//'.TRN','OLD',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
            CALL GETRAN(IFILE,TRANS)
            CLOSE(UNIT=IFILE)
          ENDIF

          ROTATE=.FALSE.
          RAXIS =.FALSE.
          RSHIFT=.FALSE.
          IF(CBBREV(CO,'ROTATE',1,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'BY',1)) THEN
              CALL PARSRL(CO(N3CO+2),3,NTRL,RL1,ERROR,*9999)
              IF(ABBREV(CO(N3CO+3),'ABOUT',1)) THEN
                CALL PARSRL(CO(N3CO+4),3,NTRL,RL2,ERROR,*9999)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,*)' rotate about',RL2(1),RL2(2),
     '              RL2(3)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                RSHIFT=.TRUE.
              ENDIF
              IF(NTCO.GT.N3CO+4.AND.ABBREV(CO(N3CO+5),'AXIS',1)) THEN
                 CALL PARSRL(CO(N3CO+6),3,NTRL,RL3,ERROR,*9999)
                 IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                 call mp_setlock()
                   WRITE(OP_STRING,*)'         axis ',RL3(1),RL3(2),
     '               RL3(3)
                   CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                 call mp_unsetlock()
                 ENDIF
                 RAXIS=.TRUE.
              ENDIF
              ROTATE=.TRUE.
            ENDIF
          ENDIF

          IF(ROTATE) THEN
            IF(RSHIFT) THEN
              AMOUNT=DSQRT(RL2(1)**2.d0+RL2(2)**2.d0+RL2(3)**2.d0)
              CALL SHIFT(RL2,-AMOUNT,TRANS,ERROR,*9999)
              CALL SHIFT(RL2,-AMOUNT,DATRAN,ERROR,*9999)
            ENDIF
            RAD1=RL1(1)*PI/180.d0
            RAD2=RL1(2)*PI/180.d0
            RAD3=RL1(3)*PI/180.d0
            IF(RAXIS) THEN
              CALL ZZ(RL3,AXIS,TRANS)
              CALL TWIST(AXIS,RAD1,TRANS,ERROR,*9999)
              CALL TWIST(AXIS,RAD1,DATRAN,ERROR,*9999)
            ELSE
              AXIS(1)=1.d0
              AXIS(2)=0.d0
              AXIS(3)=0.d0
              CALL TWIST(AXIS,RAD1,TRANS,ERROR,*9999)
              CALL TWIST(AXIS,RAD1,DATRAN,ERROR,*9999)
              AXIS(1)=0.d0
              AXIS(2)=1.d0
              AXIS(3)=0.d0
              CALL TWIST(AXIS,RAD2,TRANS,ERROR,*9999)
              CALL TWIST(AXIS,RAD2,DATRAN,ERROR,*9999)
              AXIS(1)=0.d0
              AXIS(2)=0.d0
              AXIS(3)=1.d0
              CALL TWIST(AXIS,RAD3,TRANS,ERROR,*9999)
              CALL TWIST(AXIS,RAD3,DATRAN,ERROR,*9999)
            ENDIF
            IF(RSHIFT) THEN
              CALL SHIFT(RL2,AMOUNT,TRANS,ERROR,*9999)
              CALL SHIFT(RL2,AMOUNT,DATRAN,ERROR,*9999)
            ENDIF
          ENDIF
          IF(TRANSL.OR.ROTATE) THEN
            CALL TRAINV(DATRAN,DATRAN)
            DO nd=1,NDT
              CALL ZZ(ZD(1,nd),ZD(1,nd),TRANS)
            ENDDO
            NELIST(0)=NET(1)
            DO nolist=1,NELIST(0)
              NELIST(nolist)=nolist
            ENDDO
          ENDIF

          IF(CBBREV(CO,'CENTRE',3,noco+1,NTCO,N3CO)) THEN
            DO nj=1,3
              CENTRE(nj)=0.d0
              SUM(nj)=0.d0
              AXIS(nj)=0.d0
            ENDDO
            IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTRL,CENTRE,ERROR,*9999)
            ENDIF
            DO nd=1,NDT
              DO nj=1,NJT
                SUM(nj)=ZD(nj,nd)+SUM(nj)
               ENDDO
            ENDDO
              AMOUNT=0.d0
            DO nj=1,NJT
              SUM(nj)=SUM(nj)/NDT
              AXIS(nj)=CENTRE(nj)-SUM(nj)
              AMOUNT=AMOUNT+AXIS(nj)**2
            ENDDO
            AMOUNT=DSQRT(AMOUNT)
            CALL RESET(TRANS)
            CALL SHIFT(AXIS,AMOUNT,TRANS,ERROR,*9999)
            XX(1)=0.d0
            XX(2)=0.d0
            XX(3)=0.d0
            DO nd=1,NDT
              DO nj=1,NJT
                XX(nj)=ZD(nj,nd)
              ENDDO
              CALL ZZ(XX,XX,TRANS)
              DO nj=1,NJT
                ZD(nj,nd)=XX(nj)
              ENDDO
            ENDDO
          ENDIF

          IF(NPT(1).GT.0.AND.PROJECT) THEN
            IF(NEWXI) THEN
              ERROR='>> Check code'
              GOTO 9999
C              IF(niot.ne.NJT) THEN
C                CALL NEWXID(IBT,IDO,INP,32,LD,LN,NBJF,
C     '            NDDL,NDLT,NKE,NPF,NPL,NPNE,NRE,NVJE,NVJP,NXI,
C     '            DL,SE,SQ,1.0D-6,XA,XE,XID,XP,ZD,ERROR,*9999)
C              ELSE IF(NIOT.EQ.NJT) THEN
C                DO nd=1,NDT
C                  CALL XCOORD(IBT,IDO,INP,NBJ,ne,NEELEM,
C     '              NKE,NPF,NPNE,NPNODE,NVJE,
C     '              SE,XA,XE,XID(1,nd),XP,ZD(1,nd),ERROR,*9999)
C                  LD(nd)=ne
C                  nb=NBJ(1,ne)
C                  WRITE(OP_STRING,'('' Data point '',I3,'' lies in'','
C     '              //''' element '',I3,'' at Xi coords: '',3E11.3)')
C     '              nd,ne,(XID(ni,nd),ni=1,NIT(nb))
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                ENDDO
C              ENDIF
            ELSE
              DO nd=1,NDT
                !reproject data to redefine xid,ln - need a routine
                !to calculate element projection from curvilinear coords
                !write(*,*)' reprojecting data'
                CALL ZX(ITYP10(1),ZD(1,nd),XPD)
                ZD(4,nd)=XPD(1)
                POSN=DMOD(XPD(3)/(2.d0*PI),1.d0)*NEXI(1)
                XID(1,nd)=1.d0+INT(POSN)-POSN
                SMU=0.d0
                DO n=1,NEXI(2)
                  IF((XPD(2).LT.SMU+MUSP(N)).AND.(XPD(2).GE.SMU)) THEN
                    XID(2,nd)=(XPD(2)-SMU)/MUSP(n)
                    LD(nd)=INT(POSN)+1+NEXI(1)*(NEXI(2)-n)
                  ENDIF
                  SMU=SMU+MUSP(n)
                ENDDO
                IF (XPD(2).GE.SMU) THEN
                  XID(2,nd)=(XPD(2)-SMU)/MUSP(NEXI(2))+1.d0
                  LD(nd)=INT(POSN)+1
                ENDIF
              ENDDO
            ENDIF
            CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)
            DO l=1,LN(0)
              NELIST(l)=LN(l)
            ENDDO
          ENDIF

        ELSE IF(TYPE(1:6).EQ.'FIBRES') THEN
          IF(CBBREV(CO,'SURFACE_NORM',3,noco+1,NTCO,N3CO)) THEN
            CALL FACORR(ZD,ERROR,*9999)
          ELSE IF(CBBREV(CO,'WRT_XI1',3,noco+1,NTCO,N3CO)) THEN
            DO nd=1,NDT
              nr=1
              ne=LD(nd)
              IF(ne.NE.0) THEN
                nb=NBJ(1,ne)
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '            XA(1,1,ne),XE,XP,ERROR,*9999)

                dzdxi1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,2,XID(1,nd),XE(1,3))
                dzdxi2=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,4,XID(1,nd),XE(1,3))
                b=1.0d0
                a=-1.0d0*(dzdxi2/dzdxi1)
                g1(3)=0.0d0

                dydxi1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,2,XID(1,nd),XE(1,2))
                dydxi2=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,4,XID(1,nd),XE(1,2))
                g1(2)=a*dydxi1+b*dydxi2

                dxdxi1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,2,XID(1,nd),XE(1,1))
                dxdxi2=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,4,XID(1,nd),XE(1,1))
                g1(1)=a*dxdxi1+b*dxdxi2

                length=DSQRT(g1(1)**2+g1(2)**2+g1(3)**2)
                g1(1)=g1(1)/length
                g1(2)=g1(2)/length
                g1(3)=g1(3)/length

                dzdxi1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,2,XID(1,nd),XE(1,3))
                g2(3)=dzdxi1
                dydxi1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,2,XID(1,nd),XE(1,2))
                g2(2)=dydxi1
                dxdxi1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,2,XID(1,nd),XE(1,1))
                g2(1)=dxdxi1

                length=DSQRT(g2(1)**2+g2(2)**2+g2(3)**2)
                g2(1)=g2(1)/length
                g2(2)=g2(2)/length
                g2(3)=g2(3)/length

                theta=ACOS(DOT_PROD(g1,g2))
                IF(g2(3).LT.0.0d0) THEN
                  theta=theta-PI
                ENDIF

                ZD(4,nd)=ZD(4,nd)+theta*180.d0/PI
              ENDIF
            ENDDO

          ELSE IF(CBBREV(CO,'RADIANS',3,noco+1,NTCO,N3CO)) THEN
            DO nd=1,NDT
              ZD(4,nd)=ZD(4,nd)*PI/180.0d0
            ENDDO

          ELSE
            ANGLE1=RFROMC(CO(N3CO+1))*PI/180.d0
            ANGLE2=RFROMC(CO(N3CO+2))*PI/180.d0
            IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
              XI3_1=RFROMC(CO(N3CO+1))
            ELSE
              XI3_1=0.d0
            ENDIF
            IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
              XI3_2=RFROMC(CO(N3CO+1))
            ELSE
              XI3_2=1.d0
            ENDIF
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              DO nde=1,NDLT(ne)
                nd=NDDL(ne,nde)
                IF(XID(3,nd).GE.XI3_1.AND.XID(3,nd).LE.XI3_2) THEN
                  IF(ZD(NJ_LOC(NJL_FIBR,1,nr),nd).LT.ANGLE1) THEN
                    ZD(NJ_LOC(NJL_FIBR,1,nr),nd)=
     '                ZD(NJ_LOC(NJL_FIBR,1,nr),nd)+PI
                    WRITE(OP_STRING,'('' Changed nd='',I5,'' in '','
     '                //'''element '',I3,'' at Xi3='',F5.2,'' to '','
     '                //'E12.3)') nd,ne,XID(3,nd),
     '                ZD(NJ_LOC(NJL_FIBR,1,nr),nd)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ELSE IF(ZD(NJ_LOC(NJL_FIBR,1,nr),nd).GT.ANGLE2) THEN
                    ZD(NJ_LOC(NJL_FIBR,1,nr),nd)=
     '                ZD(NJ_LOC(NJL_FIBR,1,nr),nd)-PI
                    WRITE(OP_STRING,'('' Changed nd='',I5,'' in '','
     '                //'''element '',I3,'' at Xi3='',F5.2,'' to '','
     '                //'E12.3)') nd,ne,XID(3,nd),
     '                ZD(NJ_LOC(NJL_FIBR,1,nr),nd)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF

        ELSE IF(TYPE(1:7).EQ.'WEIGHTS') THEN
          WEIGHT=RFROMC(CO(N3CO+1))

          IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
            XI3_1=RFROMC(CO(N3CO+1))
          ELSE
            XI3_1=0.d0
          ENDIF
          IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
            XI3_2=RFROMC(CO(N3CO+1))
          ELSE
            XI3_2=1.d0
          ENDIF
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            DO nde=1,NDLT(ne)
              nd=NDDL(ne,nde)
C MLB 16/4/97 nd0 and nd1 are never set so the if statement
C is not correct.
C              IF(nd.GE.nd0.AND.nd.LE.nd1) THEN
                IF(NIT(NBJ(1,ne)).EQ.3)THEN
                  IF(XID(3,nd).GE.XI3_1.AND.XID(3,nd).LE.XI3_2) THEN
                    DO nj=1,NJM
                      WD(nj,nd)=WEIGHT
                      WDL(nj,nde)=WEIGHT
                    ENDDO
                  ENDIF
                ELSE
                  DO nj=1,NJM
                    WD(nj,nd)=WEIGHT
                    WDL(nj,nde)=WEIGHT
                  ENDDO
                ENDIF !nit=3
C              ENDIF !nd between nd0 & nd1
            ENDDO !nde
          ENDDO !nolist

        ELSE IF(TYPE(1:12).EQ.'XI_INCREMENT') THEN
          CALL PARSRL(CO(N3CO+1),1,NTRL,INCREMENT,ERROR,*9999)
          DO nd=1,NDT
            XID(1,nd)=XID(1,nd)+INCREMENT
            IF(XID(1,nd).GT.1.0d0) THEN   !shift to next element
              XID(1,nd)=XID(1,nd)-1.0d0
              LD(nd)=NXI(1,1,LD(nd))
            ENDIF
          ENDDO

        ENDIF
      ENDIF

      CALL EXITS('CHDATA')
      RETURN
 9999 CALL ERRORS('CHDATA',ERROR)
      CALL EXITS('CHDATA')
      RETURN 1
      END


