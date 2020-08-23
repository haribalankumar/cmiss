      SUBROUTINE SGPROF(INDEX,IBT,IDO,ID_TYPE,INP,ISEG,ISPROF,iw,ixi,LD,
     '  NAN,NBH,NBJ,NELIST,NHE,
     '  NKHE,NKJE,NO_OBJECT,NPF,NPNE,NRE,NVHE,NVJE,NW,
     '  CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RGX,SE,
     '  XA,XE,XG,XID,XIPOS,XP,YG,YMAGN,ZA,ZD,ZE,ZG,ZP,
     '  COORDS,CSEG,CSTYPE,TYPE,ERROR,*)

C#### Subroutine: SGPROF
C###  Description:
C###    SGPROF creates profile segment ISPROF.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),ID_TYPE,INDEX,
     '  INP(NNM,NIM,NBFM),ISEG(*),ISPROF,iw,ixi,LD(NDM),
     '  NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NELIST(0:NEM),NHE(NEM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NO_OBJECT,
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),RGX(NGM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),XIPOS(*),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YMAGN,
     '  ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER COORDS*(*),CSEG(*)*(*),CSTYPE*(*),ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER i,IBEG,IEND,INDEX_OLD,N1XI,nb,nc,NCENTR,
     '  nd,ND1,ND2,ne,ni,NITB,NLEFT,nocurv,nolist,nox,noxi,nr,NRIGHT,
     '  NTCURV,NTX,NTXI,NTXMAX,nx
      PARAMETER(NTXMAX=81)
      REAL*8 DZDX(3,3),EG(3,3),PHI(3),PST(3),PTS(3,NTXMAX),
     '  PXI,R(3,3),RGX2D,RGZ,RGZ2D,RI1,RI2,RI3,RM(3,3),
     '  TC(3,3),TG(3,3),TN(3,3),TNA,
     '  U(3,3),XI(3),XPROF(NTXMAX),YPROF(NTXMAX,3),YPROFM,
     '  Z(3),Z0(3),ZVAL
      CHARACTER CHAR*8,CHAR1*1,LABEL(3)*4,STRESSTYPE*17,TITLE*60

      DATA STRESSTYPE/'Total'/

      CALL ENTERS('SGPROF',*9999)
      nc=1 !Temporary MPN 12-Nov-94
      nr=1 !should be replaced when NRE passed through
      nx=1 !Temporary

      IF(TYPE(1:5).EQ.'FIELD'.OR.TYPE(1:9).EQ.'DEPENDENT') THEN
        WRITE(CHAR1,'(I1)') ID_TYPE
        IF(TYPE(1:5).EQ.'FIELD') THEN
          TITLE=' Field variable #'//CHAR1
        ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
          TITLE=' Dependent variable #'//CHAR1
        ENDIF
        NTCURV=1
        NITB=NIT(NBJ(1,1))
        ND1=NSOBJE(3,NO_OBJECT) !is 1st nd in object
        ND2=NSOBJE(4,NO_OBJECT) !is 2nd nd in object
        nox=0
        DO nd=ND1,ND2
          ne=LD(nd)
          DO ni=1,NITB
            XI(ni)=XID(ni,nd)
          ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,'('' nd='',I4,'' Xi(ni):'',3E12.3)')
     '        nd,(XI(ni),ni=1,NITB)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(TYPE(1:5).EQ.'FIELD') THEN
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            nb=NBJ(NJ_LOC(NJL_GEOM,0,nr)+ID_TYPE,ne)
            ZVAL=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '        XE(1,NJ_LOC(NJL_GEOM,0,nr)+ID_TYPE))
          ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
            nr=1
            CALL ZPZE(NBH(1,1,ne),nc,NHE(ne),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '        ERROR,*9999)
            nb=NBH(ID_TYPE,1,ne)
            ZVAL=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '        ZE(1,ID_TYPE))
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,'('' zval='',E12.3)') ZVAL
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          nox=nox+1
          IF(nox.EQ.1) THEN
            XPROF(nox)=0.0D0
          ELSE IF(nox.GT.1) THEN
            XPROF(nox)=XPROF(nox-1)+DSQRT((ZD(1,nd)-ZD(1,nd-1))**2
     '                                  +(ZD(2,nd)-ZD(2,nd-1))**2)
          ENDIF
          YPROF(nox,1)=ZVAL
        ENDDO
        NTX=nox

      ELSE IF(TYPE(1:10).EQ.'ELASTICITY') THEN
        NITB=NIT(NBH(NH_LOC(1,nx),1,1)) !should use NEELEM(1,nr), or
                                        !NELIST(1) instead of 1 in
                                        !third index here
        IF(iw.EQ.31) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' Stresses/strains are wrt '',A,'' coords'','
     '        //'/'' XI(1)='',F5.3,'', XI(2)='',F5.3,'', XI(3)='',F5.3,'
     '        //''' NELIST(0) ='',I2,/'' NELIST ='',6I3)')
     '        COORDS,(XIPOS(ni),ni=1,3),NELIST(0),
     '        (NELIST(nolist),nolist=1,NELIST(0))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(CSTYPE(1:6).EQ.'Strain') THEN
            NTCURV=NITB
            WRITE(CHAR1,'(I1)') ixi
            IF(COORDS(1:9).EQ.'Principal') THEN
              TITLE='Principal Strain Profiles vs Xi_'//CHAR1
              LABEL(1)='E1'
              LABEL(2)='E2'
              LABEL(3)='E3'
            ELSE
              CALL STRING_TRIM(COORDS,IBEG,IEND)
              TITLE=COORDS(IBEG:IEND)//' Normal Strain Profiles vs Xi_'
     '          //CHAR1
              LABEL(1)='E11'
              LABEL(2)='E22'
              LABEL(3)='E33'
            ENDIF
          ELSE IF(CSTYPE(1:6).EQ.'Cauchy') THEN
            NTCURV=NITB
            WRITE(CHAR1,'(I1)') ixi
            IF(COORDS(1:9).EQ.'Principal') THEN
              TITLE='Principal Cauchy Stress Profiles vs Xi_'//CHAR1
              LABEL(1)='T1'
              LABEL(2)='T2'
              LABEL(3)='T3'
            ELSE
              CALL STRING_TRIM(COORDS,IBEG,IEND)
              TITLE=COORDS(IBEG:IEND)//
     '          ' Normal Cauchy Stress Profiles vs Xi_'//CHAR1
              LABEL(1)='T11'
              LABEL(2)='T22'
              LABEL(3)='T33'
            ENDIF
          ELSE IF(CSTYPE(1:7).EQ.'Nominal') THEN
            NTCURV=NITB
            WRITE(CHAR1,'(I1)') ixi
            IF(COORDS(1:9).EQ.'Principal') THEN
              TITLE='Principal Cauchy Stress Profiles vs Xi_'//CHAR1
              LABEL(1)='T1'
              LABEL(2)='T2'
              LABEL(3)='T3'
            ELSE
              CALL STRING_TRIM(COORDS,IBEG,IEND)
              TITLE=COORDS(IBEG:IEND)//
     '          ' Normal Nominal Stress Profiles vs Xi_'//CHAR1
              LABEL(1)='N11'
              LABEL(2)='N22'
              LABEL(3)='N33'
            ENDIF
          ELSE IF(CSTYPE(1:5).EQ.'Piola') THEN
            NTCURV=NITB
            WRITE(CHAR1,'(I1)') ixi
            IF(COORDS(1:9).EQ.'Principal') THEN
              TITLE='Principal Cauchy Stress Profiles vs Xi_'//CHAR1
              LABEL(1)='T1'
              LABEL(2)='T2'
              LABEL(3)='T3'
            ELSE
              CALL STRING_TRIM(COORDS,IBEG,IEND)
              TITLE=COORDS(IBEG:IEND)//
     '       ' Normal 2nd Piola Kirchhoff Stress Profiles vs Xi_'//CHAR1
              LABEL(1)='P11'
              LABEL(2)='P22'
              LABEL(3)='P33'
            ENDIF
          ENDIF
          NTX=NTXMAX
          NTXI=(NTX-1)/NELIST(0)
          NTX=NTXI*NELIST(0)+1
          nox=0
          XPROF(1)=0.0D0
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            nb=NBH(NH_LOC(1,nx),1,ne)
            NITB=NIT(nb)
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),nc,NHE(ne),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '        ERROR,*9999)
C!!!        14-Feb-89: Not correct for piecewise linear material params
            IF(ITYP1(nr,nx).EQ.5) THEN !finite elasticity
              CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,CE(1,ne),
     '          CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
            ELSE !all other problems
              CALL CPCG(NW(ne,1),nb,NPNE(1,1,ne),nr,nx,CE(1,ne),
     '          CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
            ENDIF
            IF(nolist.EQ.1) THEN
              N1XI=0
            ELSE IF(nolist.GT.1) THEN
              N1XI=1
            ENDIF
            DO noxi=N1XI,NTXI
              nox=nox+1
              XIPOS(ixi)=DBLE(noxi)/DBLE(NTXI)
              CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),NRE(ne),XE,XG,XIPOS,
     '          ERROR,*9999)
              IF(NELIST(0).EQ.1) THEN
                XPROF(nox)=XIPOS(ixi)
              ELSE IF(NELIST(0).GT.1) THEN
                IF(nox.EQ.1) THEN
                  CALL XZ(ITYP10(1),XG(1,1),Z0)
                ELSE IF(nox.GT.1) THEN
                  CALL XZ(ITYP10(1),XG(1,1),Z)
                  XPROF(nox)=XPROF(nox-1)+DSQRT((Z(1)-Z0(1))**2
     '              +(Z(2)-Z0(2))**2+(Z(3)-Z0(3))**2)
                  DO ni=1,NITB
                    Z0(ni)=Z(ni)
                  ENDDO
                ENDIF
              ENDIF
              IF(CSTYPE(1:6).EQ.'Strain') THEN
                IF(COORDS(1:9).EQ.'Principal') THEN
                  CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '              NBJ(1,ne),0,NHE(ne),NPNE(1,1,ne),nr,nx,
     '              DZDX,CE(1,ne),CG,CP,EG,PG,PHI,PST,
     '              R,RGX(1),RI1,RI2,RI3,RM,U,
     '              XE,XG,XIPOS,ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=PST(1)
                  YPROF(nox,2)=PST(2)
                  YPROF(nox,3)=PST(3)
                ELSE
                  CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '              NBJ(1,ne),0,NHE(ne),NPNE(1,1,ne),nr,nx,
     '              DZDX,CE(1,ne),CG,CP,EG,PG,PHI,PST,
     '              R,RGX(1),RI1,RI2,RI3,RM,U,
     '              XE,XG,XIPOS,ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=EG(1,1)
                  YPROF(nox,2)=EG(2,2)
                  YPROF(nox,3)=EG(3,3)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' ne='',I3,''XIPOS('',I2,'')='','
     '                //'F10.4'//','' XPROF('',I2,'')='',F10.4 )')
     '                ne,ixi,XIPOS(ixi),nox,XPROF(nox)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,
     '                '('' YPROF('',I2,'',1..3):'',3F10.4)')
     '                nox,(YPROF(nox,I),I=1,3)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ELSE IF(CSTYPE(1:6).EQ.'Cauchy') THEN
                IF(COORDS(1:9).EQ.'Principal') THEN
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50('Reference',CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=PST(1)
                  YPROF(nox,2)=PST(2)
                  YPROF(nox,3)=PST(3)
                ELSE
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50(COORDS,CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=TC(1,1)
                  YPROF(nox,2)=TC(2,2)
                  YPROF(nox,3)=TC(3,3)
                ENDIF
              ELSE IF(CSTYPE(1:7).EQ.'Nominal') THEN
                IF(COORDS(1:9).EQ.'Principal') THEN
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50('Reference',CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=PST(1)
                  YPROF(nox,2)=PST(2)
                  YPROF(nox,3)=PST(3)
                ELSE
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50(COORDS,CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=TN(1,1)
                  YPROF(nox,2)=TN(2,2)
                  YPROF(nox,3)=TN(3,3)
                ENDIF
              ELSE IF(CSTYPE(1:5).EQ.'Piola') THEN
                IF(COORDS(1:9).EQ.'Principal') THEN
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50('Reference',CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=PST(1)
                  YPROF(nox,2)=PST(2)
                  YPROF(nox,3)=PST(3)
                ELSE
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50(COORDS,CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=TG(1,1)
                  YPROF(nox,2)=TG(2,2)
                  YPROF(nox,3)=TG(3,3)
                ENDIF
              ENDIF
            ENDDO
          ENDDO

        ELSE IF(iw.EQ.32) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' Stresses/strains are wrt '',A,'' coords'','
     '        //'/'' XI(1)='',F5.3,'', XI(2)='',F5.3,'', XI(3)='','
     '        //'F5.3,'' NELIST(0) ='',I2,/'' NELIST ='',6I3)')
     '        COORDS,(XIPOS(ni),ni=1,3),NELIST(0),
     '        (NELIST(nolist),nolist=1,NELIST(0))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(CSTYPE(1:6).EQ.'Strain') THEN
            WRITE(CHAR1,'(I1)') ixi
            IF(COORDS(1:9).EQ.'Principal') THEN
              NTCURV=NITB
              TITLE='Principal Angle Profiles vs Xi_'//CHAR1
              LABEL(1)='Phi1'
              LABEL(2)='Phi2'
              LABEL(3)='Phi3'
            ELSE
              NTCURV=(NITB*(NITB-1))/2
              CALL STRING_TRIM(COORDS,IBEG,IEND)
              TITLE=COORDS(IBEG:IEND)//' Shear Strain Profiles vs Xi_'
     '          //CHAR1
              LABEL(1)='E12'
              LABEL(2)='E13'
              LABEL(3)='E23'
            ENDIF
          ELSE IF(CSTYPE(1:6).EQ.'Cauchy') THEN
            NTCURV=(NITB*(NITB-1))/2
            WRITE(CHAR1,'(I1)') ixi
            IF(COORDS(1:9).EQ.'Principal') THEN
              TITLE='Principal Angle Profiles vs Xi_'//CHAR1
              LABEL(1)='Phi1'
              LABEL(2)='Phi2'
              LABEL(3)='Phi3'
            ELSE
              CALL STRING_TRIM(COORDS,IBEG,IEND)
              TITLE=COORDS(IBEG:IEND)//
     '          ' Shear Cauchy Stress Profiles vs Xi_'//CHAR1
              LABEL(1)='T12'
              LABEL(2)='T13'
              LABEL(3)='T23'
            ENDIF
          ELSE IF(CSTYPE(1:7).EQ.'Nominal') THEN
            NTCURV=(NITB*(NITB-1))/2
            WRITE(CHAR1,'(I1)') ixi
            IF(COORDS(1:9).EQ.'Principal') THEN
              TITLE='Principal Angle Profiles vs Xi_'//CHAR1
              LABEL(1)='Phi1'
              LABEL(2)='Phi2'
              LABEL(3)='Phi3'
            ELSE
              CALL STRING_TRIM(COORDS,IBEG,IEND)
              TITLE=COORDS(IBEG:IEND)//
     '          ' Shear Nominal Stress Profiles vs Xi_'//CHAR1
              LABEL(1)='N12'
              LABEL(2)='N13'
              LABEL(3)='N23'
            ENDIF
          ELSE IF(CSTYPE(1:5).EQ.'Piola') THEN
            NTCURV=(NITB*(NITB-1))/2
            WRITE(CHAR1,'(I1)') ixi
            IF(COORDS(1:9).EQ.'Principal') THEN
              TITLE='Principal Angle Profiles vs Xi_'//CHAR1
              LABEL(1)='Phi1'
              LABEL(2)='Phi2'
              LABEL(3)='Phi3'
            ELSE
              CALL STRING_TRIM(COORDS,IBEG,IEND)
              TITLE=COORDS(IBEG:IEND)//
     '        ' Shear 2nd Piola Kirchhoff Stress Profiles vs Xi_'//CHAR1
              LABEL(1)='P12'
              LABEL(2)='P13'
              LABEL(3)='P23'
            ENDIF
          ENDIF
          NTX=NTXMAX
          NTXI=(NTX-1)/NELIST(0)
          NTX=NTXI*NELIST(0)+1
          nox=0
          XPROF(1)=0.0D0
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            NITB=NIT(NBH(NH_LOC(1,nx),1,ne))
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),nc,NHE(ne),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '        ERROR,*9999)
C!!!        14-Feb-89: Not correct for piecewise linear material params
            IF(ITYP1(nr,nx).EQ.5) THEN !finite elasticity
              CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,CE(1,ne),
     '          CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
            ELSE !all other problems
              CALL CPCG(NW(ne,1),nb,NPNE(1,1,ne),nr,nx,CE(1,ne),
     '          CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
            ENDIF
            IF(nolist.EQ.1) THEN
              N1XI=0
            ELSE IF(nolist.GT.1) THEN
              N1XI=1
            ENDIF
            DO noxi=N1XI,NTXI
              nox=nox+1
              XIPOS(ixi)=DBLE(noxi)/DBLE(NTXI)
              CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),NRE(ne),XE,XG,XIPOS,
     '          ERROR,*9999)
              IF(NELIST(0).EQ.1) THEN
                XPROF(nox)=XIPOS(ixi)
              ELSE IF(NELIST(0).GT.1) THEN
                IF(nox.EQ.1) THEN
                  CALL XZ(ITYP10(1),XG(1,1),Z0)
                ELSE IF(nox.GT.1) THEN
                  CALL XZ(ITYP10(1),XG(1,1),Z)
                  XPROF(nox)=XPROF(nox-1)+DSQRT((Z(1)-Z0(1))**2
     '              +(Z(2)-Z0(2))**2+(Z(3)-Z0(3))**2)
                  DO ni=1,NITB
                    Z0(ni)=Z(ni)
                  ENDDO
                ENDIF
              ENDIF
              IF(CSTYPE(1:6).EQ.'Strain') THEN
                IF(COORDS(1:9).EQ.'Principal') THEN
                  CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '              NBJ(1,ne),0,NHE(ne),NPNE(1,1,ne),nr,nx,
     '              DZDX,CE(1,ne),CG,CP,EG,PG,PHI,PST,
     '              R,RGX(1),RI1,RI2,RI3,RM,U,
     '              XE,XG,XIPOS,ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=PHI(1)
                  YPROF(nox,2)=PHI(2)
                  YPROF(nox,3)=PHI(3)
                ELSE
                  CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '              NBJ(1,ne),0,NHE(ne),NPNE(1,1,ne),nr,nx,
     '              DZDX,CE(1,ne),CG,CP,EG,PG,PHI,PST,
     '              R,RGX(1),RI1,RI2,RI3,RM,U,
     '              XE,XG,XIPOS,ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=EG(1,2)
                  YPROF(nox,2)=EG(1,3)
                  YPROF(nox,3)=EG(2,3)
                ENDIF
              ELSE IF(CSTYPE(1:6).EQ.'Cauchy') THEN
                IF(COORDS(1:9).EQ.'Principal') THEN
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50(COORDS,CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=PHI(1)
                  YPROF(nox,2)=PHI(2)
                  YPROF(nox,3)=PHI(3)
                ELSE
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50(COORDS,CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=TC(1,2)
                  YPROF(nox,2)=TC(1,3)
                  YPROF(nox,3)=TC(2,3)
                ENDIF
              ELSE IF(CSTYPE(1:7).EQ.'Nominal') THEN
                IF(COORDS(1:9).EQ.'Principal') THEN
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50(COORDS,CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=PHI(1)
                  YPROF(nox,2)=PHI(2)
                  YPROF(nox,3)=PHI(3)
                ELSE
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50(COORDS,CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=TN(1,2)
                  YPROF(nox,2)=TN(1,3)
                  YPROF(nox,3)=TN(2,3)
                ENDIF
              ELSE IF(CSTYPE(1:5).EQ.'Piola') THEN
                IF(COORDS(1:9).EQ.'Principal') THEN
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50(COORDS,CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=PHI(1)
                  YPROF(nox,2)=PHI(2)
                  YPROF(nox,3)=PHI(3)
                ELSE
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
                  CALL ZETX50(COORDS,CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '              NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '              NPNE(1,1,ne),nr,ne,nx,
     '              CE(1,ne),CG,CP,
     '              FEXT(1,1,ne),PG,PHI,PST,RGX(1),RGX2D,RGZ,RGZ2D,
     '              RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
                  YPROF(nox,1)=TG(1,2)
                  YPROF(nox,2)=TG(1,3)
                  YPROF(nox,3)=TG(2,3)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      IF(YMAGN.EQ.0.0D0) THEN !find YPROFM=largest value of YPROF
        YPROFM=0.0D0
        DO nox=1,NTX
          DO nocurv=1,NTCURV
            IF(DABS(YPROF(nox,nocurv)).GT.YPROFM) THEN
              YPROFM=DABS(YPROF(nox,nocurv))
            ENDIF
          ENDDO
        ENDDO
      ELSE
        YPROFM=YMAGN
        IF(TYPE(1:10).EQ.'ELASTICITY') THEN
          IF(COORDS(1:9).EQ.'Principal'.AND.iw.EQ.32) YPROFM=90.0D0
        ENDIF
      ENDIF
      DO nox=1,NTX !to normalize x-axis (to 1) and y-axis (to YPROFM)
        XPROF(nox)=XPROF(nox)/XPROF(NTX)
        DO nocurv=1,NTCURV
          YPROF(nox,nocurv)=YPROF(nox,nocurv)/YPROFM
        ENDDO
      ENDDO

      CALL OPEN_SEGMENT(ISPROF,ISEG,iw,'PROF',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)
      PTS(1,1)=0.0D0
      PTS(2,1)=0.0D0
      PTS(1,2)=1.0D0
      PTS(2,2)=0.0D0
      CALL POLYLINE(1,iw,2,PTS,ERROR,*9999)
      PTS(1,1)=0.0D0
      PTS(2,1)=-1.0D0
      PTS(1,2)=0.0D0
      PTS(2,2)=1.0D0
      CALL POLYLINE(1,iw,2,PTS,ERROR,*9999)
      IF(TYPE(1:10).EQ.'ELASTICITY') THEN
        PTS(2,1)=-0.95D0
        PTS(2,2)= 0.95D0
        DO nolist=1,NELIST(0)-1
          nox=nolist*(NTXI+1)+1
          PTS(1,1)=XPROF(nox)
          PTS(1,2)=XPROF(nox)
          CALL POLYLINE(1,iw,2,PTS,ERROR,*9999)
        ENDDO
        PTS(2,1)=0.0D0
        PTS(2,2)=0.03D0
        DO I=1,2*NELIST(0)
          PTS(1,1)=DBLE(I)/DBLE(2*NELIST(0))
          PTS(1,2)=PTS(1,1)
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
          IF(MOD(I,2).EQ.0) THEN
            WRITE(CHAR,'(F4.2)') PTS(1,1)
            CALL STRING_TRIM(CHAR,IBEG,IEND)
            CALL TEXT(INDEX,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF
      PTS(1,2)=0.02D0
      DO I=1,5
        PTS(1,1)=0.0D0
        PTS(2,1)=DBLE(I-1)/2.0D0-1.0D0
        PTS(2,2)=PTS(2,1)
        CALL POLYLINE(1,iw,2,PTS,ERROR,*9999)
        WRITE(CHAR,'(F8.2)') PTS(2,1)*YPROFM
        CALL STRING_TRIM(CHAR,IBEG,IEND)
        PTS(1,1)=-0.01D0
        CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
      ENDDO
      PTS(1,1)=0.5D0
      PTS(2,1)=-1.0D0
      CALL STRING_TRIM(TITLE,IBEG,IEND)
      CALL TEXT(INDEX,iw,TITLE(IBEG:IEND),PTS(1,1),ERROR,*9999)

      NLEFT =NTX/3
      NCENTR=NTX/2
      NRIGHT=2*NCENTR-NLEFT
      DO nox=1,NTX
        PTS(1,nox)=XPROF(nox)
        PTS(2,nox)=YPROF(nox,1)
      ENDDO
      CALL POLYLINE(INDEX,iw,NTX,PTS,ERROR,*9999)
      CALL STRING_TRIM(LABEL(1),IBEG,IEND)
      PTS(1,1)=XPROF(NLEFT)
      PTS(2,1)=YPROF(NLEFT,1)
      CALL TEXT(INDEX,iw,LABEL(1)(IBEG:IEND),PTS(1,1),ERROR,*9999)
      IF(NTCURV.GE.2) THEN
        DO nox=1,NTX
          PTS(1,nox)=XPROF(nox)
          PTS(2,nox)=YPROF(nox,2)
        ENDDO
        CALL POLYLINE(INDEX,iw,NTX,PTS,ERROR,*9999)
        CALL STRING_TRIM(LABEL(2),IBEG,IEND)
        PTS(1,1)=XPROF(NCENTR)
        PTS(2,1)=YPROF(NCENTR,2)
        CALL TEXT(INDEX,iw,LABEL(2)(IBEG:IEND),PTS(1,1),ERROR,*9999)
      ENDIF
      IF(NTCURV.GE.3) THEN
        DO nox=1,NTX
          PTS(1,nox)=XPROF(nox)
          PTS(2,nox)=YPROF(nox,3)
        ENDDO
        CALL POLYLINE(INDEX,iw,NTX,PTS,ERROR,*9999)
        CALL STRING_TRIM(LABEL(3),IBEG,IEND)
        PTS(1,1)=XPROF(NRIGHT)
        PTS(2,1)=YPROF(NRIGHT,3)
        CALL TEXT(INDEX,iw,LABEL(3)(IBEG:IEND),PTS(1,1),ERROR,*9999)
      ENDIF

      CALL CLOSE_SEGMENT(ISPROF,iw,ERROR,*9999)

      CALL EXITS('SGPROF')
      RETURN
 9999 CALL ERRORS('SGPROF',ERROR)
      CALL EXITS('SGPROF')
      RETURN 1
      END


