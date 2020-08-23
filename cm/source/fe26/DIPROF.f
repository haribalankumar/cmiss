      SUBROUTINE DIPROF(IBT,IDO,INP,ISEG,ISPROF,LD,NAN,NBH,NBJ,
     '  NEELEM,NELIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '  NPNODE,NRE,NVHE,NVHP,NVJE,NW,NYNE,NYNP,
     '  CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RGX,SE,
     '  XA,XE,XG,XID,XP,YG,YP,ZA,ZD,ZE,ZG,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: DIPROF
C###  Description:
C###    DIPROF displays profiles of normal, shear and principal
C###    stresses or strains in 3D problems.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISEG(*),ISPROF(2),LD(NDM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),RGX(NGM),
     '  SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),
     '  ZD(NJM,NDM),ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),
     '  STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG2,IBEG3,ID_TYPE,IEND,IEND2,IEND3,IFROMC,INDEX,
     '  INDEX_POLYLINE,IW,IWK(2),IXI,N3CO,NOELEM,NOOBJE,NO_OBJECT,
     '  NOPROF,nr,NTIW,nx
      REAL*8 RFROMC,XIPOS(3),YMAGN
      CHARACTER COORDS*9,CSTYPE*7,OBJECT*20,TYPE*10
      LOGICAL CBBREV

      CALL ENTERS('DIPROF',*9999)

      nr=1 !temporarily
      nx=1 !temporary cpb 22/11/94

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        IF(NJ_LOC(njl_fibr,0,nr).GT.0) THEN  !field variables defined
          TYPE='field'
        ELSE IF(ITYP1(nr,nx).GT.0) THEN !dependent variables defined
          TYPE='dependent'
        ELSE
          TYPE='field'
        ENDIF
        CALL STRING_TRIM(TYPE,IBEG2,IEND2)
        IF(NTOBJE.GT.0) THEN
          OBJECT=OBJECT_NAME(NTOBJE)
        ELSE
          OBJECT=' '
        ENDIF
        CALL STRING_TRIM(OBJECT,IBEG3,IEND3)

C---------------------------------------------------------------------

C#### Command: FEM display profile
C###  Parameter:       <field/dependent ID#[dependent]>
C###   Specify the the ID# to display
C###  Parameter:       <with OBJECT_NAME>
C###   Specify the object name
C###  Parameter:       <rgb=RGB[blue]>
C###    Define colour (eg red,green,blue,cyan)
C###  Description:
C###    Displays a profile.

        OP_STRING(1)=STRING(1:IEND) //' <(field/dependent) ID#['
     '                              //TYPE(1:IEND2)//' 1]>'
        OP_STRING(2)=BLANK(1:15)//'<with OBJECT_NAME['
     '                              //OBJECT(1:IEND3)//']>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[blue]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM display profile strain/cauchy/nominal/piola
C###  Parameter:       <in (ELEMENTS/all)[all]>
C###    Specify which elements to dispaly the profile of.
C###  Parameter:       <at XI_DIRECTION[1]>
C###    Specify the xi location.
C###  Parameter:       <xi_1 VALUE[0.5]>
C###    Specify the xi 1 value
C###  Parameter:       <xi_2 VALUE[0.5]>
C###    Specify the xi 2 value
C###  Parameter:       <xi_3 VALUE[0.0]>
C###    Specify the xi 3 value
C###  Parameter:       <(fibre/principal/reference/wall)[fibre]>
C###  Parameter:       <maximum MAGNITUDE>
C###    Specify a maximum magnitude
C###  Parameter:       <on WS[both]>
C###    Specify a workstation number
C###  Parameter:       <rgb=RGB[blue]>
C###    Define colour (eg red,green,blue,cyan)
C###  Description:
C###    Displays a specifyed component of stress distribution.

        OP_STRING(1)=STRING(1:IEND) //' strain/cauchy/nominal/piola'
        OP_STRING(2)=BLANK(1:15)//'<in (ELEMENTS/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<at XI_DIRECTION[1]>'
        OP_STRING(4)=BLANK(1:15)//'<xi_1 VALUE[0.5]>'
        OP_STRING(5)=BLANK(1:15)//'<xi_2 VALUE[0.5]>'
        OP_STRING(6)=BLANK(1:15)//'<xi_3 VALUE[0.0]>'
        OP_STRING(7)=BLANK(1:15)
     '    //'<(fibre/principal/reference/wall)[fibre]>'
        OP_STRING(8)=BLANK(1:15)//'<maximum MAGNITUDE>'
        OP_STRING(9)=BLANK(1:15)//'<on WS[both]>'
        OP_STRING(10)=BLANK(1:15)//'<rgb=RGB[blue]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','DIPROF',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL ASSERT(NJ_LOC(njl_fibr,0,0).GT.0.OR.
     '    ITYP2(nr,nx).GT.0,
     '    '>>No field or dependent variable defined',ERROR,*9999)
        IF(CBBREV(CO,'FIELD',1,noco+1,NTCO,N3CO)) THEN
          TYPE='FIELD'
          ID_TYPE=IFROMC(CO(N3CO+1))
        ELSE IF(CBBREV(CO,'DEPENDENT',1,noco+1,NTCO,N3CO)) THEN
          TYPE='DEPENDENT'
          ID_TYPE=IFROMC(CO(N3CO+1))
        ELSE IF(CBBREV(CO,'CAUCHY',1,noco+1,noco+3,N3CO)) THEN
          TYPE='ELASTICITY'
          CSTYPE='Cauchy'
        ELSE IF(CBBREV(CO,'NOMINAL',1,noco+1,noco+3,N3CO)) THEN
          TYPE='ELASTICITY'
          CSTYPE='Nominal'
        ELSE IF(CBBREV(CO,'PIOLA',1,noco+1,noco+3,N3CO)) THEN
          TYPE='ELASTICITY'
          CSTYPE='Piola'
        ELSE
          IF(NJ_LOC(NJL_FIBR,0,nr).GT.0) THEN  !field variables defined
            TYPE='FIELD'
          ELSE IF(ITYP2(nr,nx).GT.0) THEN !dependent variables defined
            TYPE='DEPENDENT'
          ELSE
            TYPE='FIELD'
          ENDIF
          ID_TYPE=1
        ENDIF

        IF(CBBREV(CO,'IN',1,noco+1,noco+3,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NEM,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          NELIST(0)=0
          DO nr=1,NRT
            DO NOELEM=NELIST(0)+1,NELIST(0)+NEELEM(0,NR)
              NELIST(NOELEM)=NEELEM(NOELEM,NR)
            ENDDO
            NELIST(0)=NELIST(0)+NEELEM(0,NR)
          ENDDO
        ENDIF

        IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
          IXI=IFROMC(CO(N3CO+1))
        ELSE
          IXI=1
        ENDIF
        IF(IXI.EQ.1) THEN
          IF(CBBREV(CO,'XI_2',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(2)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(2)=0.5D0
          ENDIF
          IF(CBBREV(CO,'XI_3',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(3)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(3)=0.0D0
          ENDIF
        ELSE IF(IXI.EQ.2) THEN
          IF(CBBREV(CO,'XI_1',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(1)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(1)=0.5D0
          ENDIF
          IF(CBBREV(CO,'XI_3',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(3)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(3)=0.0D0
          ENDIF
        ELSE IF(IXI.EQ.3) THEN
          IF(CBBREV(CO,'XI_1',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(1)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(1)=0.5D0
          ENDIF
          IF(CBBREV(CO,'XI_2',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(2)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(2)=0.5D0
          ENDIF
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
        ENDIF

        IF(TYPE(1:10).EQ.'ELASTICITY') THEN
          IF(CBBREV(CO,'REFERENCE',3,noco+1,NTCO,N3CO)) THEN
            COORDS='Reference'
          ELSE IF(CBBREV(CO,'PRINCIPAL',3,noco+1,NTCO,N3CO)) THEN
            COORDS='Principal'
          ELSE IF(CBBREV(CO,'WALL',3,noco+1,NTCO,N3CO)) THEN
            COORDS='Wall'
          ELSE
            COORDS='Fibre'
          ENDIF
          IF(CBBREV(CO,'MAXIMUM',1,noco+1,NTCO,N3CO)) THEN
            YMAGN=RFROMC(CO(N3CO+1))
          ELSE
            YMAGN=0.0D0
          ENDIF
          IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
            NTIW=1
            IWK(1)=IFROMC(CO(N3CO+1))
            IF(IWK(1).NE.31.AND.IWK(1).NE.32) IWK(1)=31
          ELSE
            NTIW=2
            IWK(1)=31
            IWK(2)=32
          ENDIF

        ELSE IF(TYPE(1:5).EQ.'FIELD'.OR.TYPE(1:9).EQ.'DEPENDENT') THEN
          CALL ASSERT(NTOBJE.GT.0,'>>No object defined',ERROR,*9999)
          IF(CBBREV(CO,'OBJECT',1,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            OBJECT=CO(N3CO+1)(IBEG:IEND)
          ELSE
            CALL STRING_TRIM(OBJECT_NAME(NTOBJE),IBEG,IEND)
            OBJECT=OBJECT_NAME(NTOBJE)(IBEG:IEND)
          ENDIF
          CALL STRING_TRIM(OBJECT,IBEG,IEND)
          DO NOOBJE=1,NTOBJE
            IF(OBJECT_NAME(NOOBJE)(IBEG:IEND).EQ.OBJECT(IBEG:IEND)) THEN
              NO_OBJECT=NOOBJE
            ENDIF
          ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,'(1X,A,'' is object #'',I2)')
     '        OBJECT(IBEG:IEND),NO_OBJECT
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          NTIW=1
          IWK(1)=31
        ENDIF

        CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '    NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,
     '    NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)

        DO NOPROF=1,NTIW
          IW=IWK(NOPROF)
          CALL ACWK(IW,0,ERROR,*9999)
          CALL SGPROF(INDEX,IBT,IDO,ID_TYPE,INP,ISEG,ISPROF(NOPROF),IW,
     '      IXI,LD,NAN,NBH,NBJ,NELIST,NHE(1,nx),NKHE,NKJE,NO_OBJECT,
     '      NPF,NPNE,NRE,NVHE,NVJE,NW(1,1,nx),
     '      CE(1,1,nx),CG,CGE(1,1,1,nx),
     '      CP(1,1,nx),CURVCORRECT,FEXT,PG,RGX,SE,XA,
     '      XE,XG,XID,XIPOS,XP,YG,YMAGN,ZA,ZD,ZE,ZG,ZP,COORDS,CSEG,
     '      CSTYPE,TYPE,ERROR,*9999)
          CALL DAWK(IW,0,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('DIPROF')
      RETURN
 9999 CALL ERRORS('DIPROF',ERROR)
      CALL EXITS('DIPROF')
      RETURN 1
      END


