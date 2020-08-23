      SUBROUTINE UPMESH(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,ISFIEL,
     '  ISLINE,ISLINO,ISNONO,MXI,NAN,NBH,NBJ,NEELEM,NELIST,NEP,NENP,
     &  NGAP,NHE,NHP,NKH,NKHE,NKJE,NLATNE,NLL,NLLIST,NORD,NPF,NPL,NPNE,
     &  NPNODE,NQNE,NQNLAT,NQS,NQXI,NRE,NVHE,NVHP,NVJE,NVJP,NW,NXI,NYNE,
     &  NYNP,CE,CURVCORRECT,DL,SE,XA,XAB,XE,XG,XIP,XP,XQ,YG,YP,YQ,ZA,ZE,
     &  ZP,CSEG,STRING,FIX,ERROR,*)

C#### Subroutine: UPMESH
C###  Description:
C###    When the UPDATE MESH POSITION option is used UPMESH updates
C###    the nodal positions of the slave mesh which is embedded
C###    in a host mesh using either the deformed or undeformed geometry
C###    of the host mesh, or the fitted geometry and the xi positions of
C###    those nodes. If updating from a fitted geometry, the nodal
C###    derivatives are also updated from the host mesh. For
C###    spherically based coordinate system a recheck (and if necessary
C###    an adjustment) of the connectivity is calculated to ensure
C###    consistency of angular coordinate increase in coronary elements
C###    with host element xi direction.
C###  The other use of UPMESH is to update mesh geometry (e.g. radii) and
C###  connectivity.      

C###    UPMESH calls UPVIEW to update mesh on workstation viewports.
C###    Calls SEGPAC to repack segments to eliminate zeros in ISEG.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISELNO(NWM,NEM),ISFIBR(NWM,NEM,NGRSEGM),ISFIEL(NWM,NEM),
     '  ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),ISNONO(NWM,NPM),MXI(2,NEM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     &  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NEP(NPM),
     &  NENP(NPM,0:NEPM,0:NRM),NGAP(NIM,NBM),NHE(NEM,NXM),
     &  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     &  NKJE(NKM,NNM,NJM,NEM),NLATNE(NEQM+1),NLL(12,NEM),NLLIST(0:NLM),
     &  NORD(5,NE_R_M),NPF(9,NFM),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),NQNE(NEQM,NQEM),NQNLAT(NEQM*NQEM),
     &  NQS(NEQM),NQXI(0:NIM,NQSCM),NRE(NEM),NRLIST(0:NRM),NVHE(NNM,
     &  NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     &  NVJP(NJM,NPM),NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,NEM),
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),
     &  SE(NSM,NBFM,NEM),XA(NAM,NJM),
     &  XE(NSM,NJM),XG(NJM,NUM),XIP(NIM,NPM),XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),YG(NIYGM,NGM,NEM),
     &  YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),ZA(NAM,NHM,NCM,NEM),
     &  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,IWK(6),N3CO,nb1,nb2,ne1,ne2,ne_start,ni,
     &  nj,njj,nj_volume,noelem,noiw,nonode,np1,np2,nr,nr1,nr2,NTIW,nv,
     &  nx
      REAL*8 AA_RATIO,ACINUS_LENGTH,ACINUS_VOLUME,ld_ratio,PXI,
     &  RATIO_DIAMETER,RFROMC,scale_f,value
      LOGICAL ALL_REGIONS,CALL_UPVIEW,CBBREV,CONTAIN,CUBE_ROOT,DEFORM,
     &  FIT,SCALED,UNSTRAINED,UPDATE_GEOMETRY,UPDATE_POSITION
      CHARACTER RADIUS_SCHEME*255
C     INTEGER njj,nj_field,nj_max,nk,nu,
C     REAL*8 DEF(3,3),DET,OLD(3,11),TEMP(3,3),TRANS(3,3),UNDEF(3,3)

      CALL ENTERS('UPMESH',*9999)

      nx=1 !temporary cpb 22/11/94
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------


C#### Command: FEM update mesh position
C###  Parameter:      <of REGION#[2]>
C###  Specify the region of the mesh to be updated (slave)
C###  Parameter:      <in REGION#[1]>
C###  Specify the region of the mesh used to update (host)
C###  Parameter:    <(undeformed/deformed/fit)[undeformed]>
C###  Specify whether undeformed or deformed or fitted coordinates
C###  are used to update global positions
C###  Parameter:    <contain>
C###    Specify whether only the nodes contained in the host region will
C###    be updated.  If this option is not specified, then if some
C###    nodes in the second region are not contained in the first, the
C###    command will not execute.  Specifying this option forces an
C###    update on all nodes contained within the first region.
C###  Description:
C###    Update global positions XP of a slave mesh from Xi positions
C###    XIP within a host mesh.


        OP_STRING(1)=STRING(1:IEND)//' position'
        OP_STRING(2)=BLANK(1:15)//'<of REGION#[2]> <in REGION#[1]>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<(undeformed/deformed/fit)[undeformed]>'
        OP_STRING(4)=BLANK(1:15)//'<contain>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update mesh geometry 
C###  Parameter:      ->    (radius_field #[1])
C###    Specify that the update will operate on the radius field.
C###  Parameter:      ->    <elements #/all[all]>
C###    Specify the elements (list or group) to update.
C###  Parameter:      ->    <region #[all]>
C###    Specify the region number(s) to update.
C###  Parameter:      <constant/ld_ratio/Rd_S/Rd_H/anatomical/lumped [constant]
C###    Specify which radius scheme to use.
C###  Parameter:      <constant>
C###    Specify that the radius is a constant.
C###  Parameter:           <VALUE>
C###    Specify the constant radius value.
C###  Parameter:      <ld_ratio>
C###    Specify that the radius is calculated using the length to
C###    diameter ratio.
C###  Parameter:           <VALUE>
C###    Specify the length to diameter ratio value.
C###  Parameter:      <Rd_Strahler:>
C###    Specify that the radius is based on Strahler diameter ratio.
C###  Parameter:           <ratio_diameter VALUE>
C###    Specify the Strahler diameter ratio.
C###  Parameter:           <cube_root>
C###    Specify that the Strahler diameter ratio = Strahler branching
C###    ratio ^ (1/3).
C###  Parameter:      <Rd_Horsfield>
C###    Specify that the radius is based on Horsfield diameter ratio.
C###  Parameter:           <ratio_diameter VALUE>
C###    Specify the Horsfield diameter ratio.
C###  Parameter:           <cube_root>
C###    Specify that the Horsfield diameter ratio = Horsfield branching
C###    ratio ^ (1/3).
C###  Parameter:      <anatomical>
C###    Specify that the radius is calulated using anatomical values.
C###  Parameter:           <HBW/Weibel/Horsfield/Arteries/Veins>
C###    Specify the anatomical data to use: HBW refers to respiratory
C###    airway data (by generation) from Haefeli-Bleuer & Weibel
C###    (1988), Weibel refers to conducting airway data (by generation)
C###    from Weibel (1963), Horsfield refers to conducting airway data
C###    (by Horsfield order) from Horsfield et al (1966). Arteries
C###    refers to data from Horsfield et al (1978) for the arteries.
C###    Veins, refers to data for veins from Horsfield et al (1981).
C###  Parameter:      <lumped>
C###    Specify that the update is for a 'lumped parameter model'.
C###  Parameter:           <volume VALUE [82.76]>
C###    Specify the volume of the lumped parameter model.
C###  Parameter:           <length VALUE [5.0]>
C###    Specify the length of the lumped parameter model.
C###  Parameter:           <capillary>
C###    Specifies all capillary parameters that are required for
C###    solution.
C###  Parameter:           <scale_factor [3.0144]>
C###    Specifies the scale factor used to convert from CMISS units to mm. 
C###  Parameter:           <aa_ratio VALUE [1.0]>
C###    Specify the ratio of duct to alveolar x-sectional area.

C###  Description: Updates the specified radius field with
C###    a constant or with values based on anatomical data or diameter ratios.

        OP_STRING(1)=STRING(1:IEND)//' geometry'
        OP_STRING(2)=BLANK(1:15)//'<radius_field #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<elements #/all[all]>'
        OP_STRING(4)=BLANK(1:15)//'<region #/all[all]>'
        OP_STRING(5)=BLANK(1:15)//'<constant/ld_ratio/Rd_S/Rd_H/'
     &   //'anatomical/'//'lumped [constant]>'
        OP_STRING(6)=BLANK(1:15)//'<for constant/ld_ratio:>'
        OP_STRING(7)=BLANK(1:15)//'    <VALUE>'
        OP_STRING(8)=BLANK(1:15)//'<for Rd_Strahler/Rd_Horsfield:>'
        OP_STRING(9)=BLANK(1:15)//'    <ratio_diameter VALUE>'
        OP_STRING(10)=BLANK(1:15)//'    <cube_root>'
        OP_STRING(11)=BLANK(1:15)//'<for anatomical:>'
        OP_STRING(12)=BLANK(1:15)//
     &    '   <sym/HBW/Weibel/Horsfield/comb/Arteries/Veins>'
        OP_STRING(13)=BLANK(1:15)//'<for lumped:>'
        OP_STRING(14)=BLANK(1:15)//'    <volume VALUE [82.76]>'
        OP_STRING(15)=BLANK(1:15)//'    <length VALUE [5.0]>'
        OP_STRING(15)=BLANK(1:15)//'<capillary>'
        OP_STRING(16)=BLANK(1:15)//'       <scale_factor [3.0144]>'
        OP_STRING(17)=BLANK(1:15)//'    <length VALUE   [5.0]>'
        OP_STRING(18)=BLANK(1:15)//'    <aa_ratio VALUE [1.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPMESH',ERROR,*9999)
      ELSE
        CALL_UPVIEW=.TRUE.
        IF(CBBREV(CO,'POSITION',1,noco+1,NTCO,N3CO)) THEN
          UPDATE_POSITION=.TRUE.
          UPDATE_GEOMETRY=.FALSE.
          IF(CBBREV(CO,'OF',2,noco+1,NTCO,N3CO)) THEN
            nr2=IFROMC(CO(N3CO+1))
          ELSE
            nr2=2
          ENDIF
          IF(CBBREV(CO,'IN',2,noco+1,NTCO,N3CO)) THEN
            nr1=IFROMC(CO(N3CO+1))
          ELSE
            nr1=1
          ENDIF
          IF(CBBREV(CO,'DEFORMED',2,noco+1,NTCO,N3CO)) THEN
            DEFORM=.TRUE.
          ELSE
            DEFORM=.FALSE.
          ENDIF
          IF(CBBREV(CO,'FIT',2,noco+1,NTCO,N3CO)) THEN
            FIT=.TRUE.
          ELSE
            FIT=.FALSE.
          ENDIF
          IF(CBBREV(CO,'CONTAIN',2,noco+1,NTCO,N3CO)) THEN
            CONTAIN=.TRUE.
          ELSE
            CONTAIN=.FALSE.
          ENDIF

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' Update position of region '',I1,'
     '        //''' in region '',I1)') nr2,nr1
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ELSE IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
          UPDATE_POSITION=.FALSE.
          UPDATE_GEOMETRY=.FALSE.
          CALL PARSIL(CO(N3CO+1),3,NTIW,IWK,ERROR,*9999)
        ELSE IF(CBBREV(CO,'GEOMETRY',3,noco+1,NTCO,N3CO)) THEN
          UPDATE_POSITION=.FALSE.
          UPDATE_GEOMETRY=.TRUE.
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,
     &      *9999)
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
          nr=NRLIST(1)
          IF(CBBREV(CO,'RADIUS_FIELD',3,noco+1,NTCO,N3CO)) THEN
            njj=IFROMC(CO(N3CO+1)) !field number for radius
C           Check that the field exists
            CALL ASSERT(NJ_LOC(NJL_FIEL,njj,nr).GT.0,
     '        '>>Define field first',ERROR,*9999)
            nj_radius=NJ_LOC(NJL_FIEL,njj,nr)
            ne_start = 1
            IF(CBBREV(CO,'SCALE',3,noco+1,NTCO,N3CO))THEN
              scale_f=RFROMC(CO(N3CO+1))
              SCALED=.TRUE.
            ELSE
              scale_f=1.d0
              SCALED=.FALSE.
            ENDIF
            IF(CBBREV(CO,'CONSTANT',3,noco+1,NTCO,N3CO))THEN
              RADIUS_SCHEME='CONSTANT'
              value=RFROMC(CO(N3CO+1))
            ELSE IF(CBBREV(CO,'LD_RATIO',2,noco+1,NTCO,N3CO))THEN
              RADIUS_SCHEME='LD_RATIO'
              ld_ratio=RFROMC(CO(N3CO+1)) !L:D ratio for radius calculation
            ELSE IF(CBBREV(CO,'RD_STRAHLER',4,noco+1,NTCO,N3CO)) THEN
              RADIUS_SCHEME='RDSTRAHLER'
              IF(CBBREV(CO,'RATIO_DIAMETER',3,noco+1,NTCO,N3CO))THEN
                RATIO_DIAMETER=RFROMC(CO(N3CO+1))
                CUBE_ROOT=.FALSE.
              ELSEIF(CBBREV(CO,'CUBE_ROOT',3,noco+1,NTCO,N3CO))THEN
                CUBE_ROOT=.TRUE.
              ENDIF
            ELSE IF(CBBREV(CO,'RD_HORSFIELD',4,noco+1,NTCO,N3CO)) THEN
              RADIUS_SCHEME='RDHORSFIELD'
              IF(CBBREV(CO,'RATIO_DIAMETER',3,noco+1,NTCO,N3CO))THEN
                RATIO_DIAMETER=RFROMC(CO(N3CO+1))
                CUBE_ROOT=.FALSE.
              ELSEIF(CBBREV(CO,'CUBE_ROOT',3,noco+1,NTCO,N3CO))THEN
                CUBE_ROOT=.TRUE.
              ENDIF
            ELSE IF(CBBREV(CO,'ANATOMICAL',4,noco+1,NTCO,N3CO))
     &          THEN
              IF(CBBREV(CO,'HBW',2,noco+1,NTCO,N3CO)) THEN
                RADIUS_SCHEME='HBW'
             ELSE IF(CBBREV(CO,'WEIBEL',2,noco+1,NTCO,N3CO)) THEN
                RADIUS_SCHEME='WEIBEL'
              ELSE IF(CBBREV(CO,'HORSFIELD',2,noco+1,NTCO,N3CO)) THEN
                RADIUS_SCHEME='HORSFIELD'
              ELSE IF(CBBREV(CO,'COMBINED',2,noco+1,NTCO,N3CO)) THEN
                RADIUS_SCHEME='COMBINED'
              ELSE IF(CBBREV(CO,'ARTERIES',2,noco+1,NTCO,N3CO)) THEN
                RADIUS_SCHEME='ARTERIES'
              ELSE IF(CBBREV(CO,'VEINS',2,noco+1,NTCO,N3CO)) THEN
                RADIUS_SCHEME='VEINS'
              ENDIF
              UNSTRAINED=.FALSE.
              IF(RADIUS_SCHEME(1:8).EQ.'ARTERIES'
     &          .OR.RADIUS_SCHEME(1:8).EQ.'VEINS') THEN
                IF(CBBREV(CO,'UNSTRAINED',5,noco+1,NTCO,N3CO))
     &            UNSTRAINED=.TRUE. !scale anatomical info to unstrained radius values
              ENDIF
            ELSE IF(CBBREV(CO,'CE_FIELD',4,noco+1,NTCO,N3CO))THEN
              RADIUS_SCHEME='CE_FIELD'
            ELSE IF(CBBREV(CO,'LUMPED',3,noco+1,NTCO,N3CO))THEN
              RADIUS_SCHEME='LUMPED'
              IF(CBBREV(CO,'VOLUME',3,noco+1,NTCO,N3CO))THEN
                ACINUS_VOLUME=RFROMC(CO(N3CO+1))/NELIST(0)
              ELSE
                ACINUS_VOLUME=82.76d0
              ENDIF
              IF(CBBREV(CO,'LENGTH',3,noco+1,NTCO,N3CO))THEN
                ACINUS_LENGTH=RFROMC(CO(N3CO+1))
              ELSE
                ACINUS_LENGTH=5.d0
              ENDIF
            ENDIF
            IF(CBBREV(CO,'START_ELEM',2,noco+1,NTCO,N3CO))THEN
              ne_start=IFROMC(CO(N3CO+1)) !'starting' element for assigning radius
            ENDIF
          ELSEIF(CBBREV(CO,'CAPILLARY',3,noco+1,NTCO,N3CO)) THEN
C      KSB 28/06/04 - currently pulmonary capillary model not set up to use
C      fields, diameter values stored in CE, this should be updated.
            IF(CBBREV(CO,'SCALE_FACTOR',3,noco+1,NTCO,N3CO))THEN
              SCALE_FACTOR=RFROMC(CO(N3CO+1))
            ELSEIF(.NOT.CALL_MESH) THEN
              SCALE_FACTOR=3.0144d0 !default value
            ENDIF
            !defines capillary mesh parameters
            CALL CAP_NE(NBJ,NELIST,NENP,NPNE,nr,CE,XP,ERROR,*9999)
            UPDATE_GEOMETRY=.FALSE. !don't call UPMESH_GEOMETRY
            CALL_UPVIEW=.FALSE.
          ENDIF
          IF(CBBREV(CO,'AA_RATIO',2,noco+1,NTCO,N3CO))THEN
            AA_RATIO=RFROMC(CO(N3CO+1))
          ELSE
            AA_RATIO=1.d0
          ENDIF
        ELSE
          UPDATE_POSITION=.FALSE.
          NTIW=2*NJT-3+IMAP
          DO noiw=1,NTIW
            IWK(noiw)=noiw
          ENDDO
        ENDIF
        
        IF(UPDATE_POSITION) THEN
          DO nonode=1,NPNODE(0,nr2) !loop over region 2 nodes
            np2=NPNODE(nonode,nr2)  !region 2 node
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np2,.FALSE.,ERROR,*9999)
            ne1=NEP(np2) !region 1 element containing np2
            IF(.NOT.CONTAIN)THEN
C!!! LKC 6-AUG-2002 Should really check the ne1 is valid.
C MHT 4-MAR-2003 added contain option, so assert only done if option 
C not specified. Then is safe to put in the IF statement below.
              CALL ASSERT(ne1.GT.0,'>> No element associated with node',
     '          ERROR,*9999)
            ENDIF
            IF(ne1.NE.0)THEN !i.e. only possible if CONTAIN is set .TRUE.
              nb1=NBJ(1,ne1) !region 1 basis
              IF(DEFORM) THEN !uses the deformed host mesh geometry
                CALL ZPZE(NBH(1,1,ne1),1,NHE(ne1,nx),NKHE(1,1,1,ne1),
     '            NPF(1,1),NPNE(1,1,ne1),NRE(ne1),NVHE(1,1,1,ne1),
     '            NW(ne1,1,nx),nx,CURVCORRECT(1,1,1,ne1),SE(1,1,ne1),
     '            ZA(1,1,1,ne1),ZE,ZP,ERROR,*9999)
                DO nj=1,NJ_LOC(NJL_GEOM,0,NRE(ne1))
                  nb1=NBJ(nj,ne1) !region 1 basis
                  XP(1,1,nj,np2)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),
     '              INP(1,1,nb1),nb1,1,XIP(1,np2),ZE(1,nj))
                  IF (ITYP10(nr1).EQ.4)THEN !prolate spheroidal coordinates
C ensures the theta coordinates are between zero and 2pi
                    IF(nj.EQ.NJ_LOC(NJL_GEOM,3,nr2)) THEN !theta coordinate
                      DO WHILE (XP(1,1,nj,np2).GT.(2.0d0*PI))
                        XP(1,1,nj,np2)=XP(1,1,nj,np2)-(2.0d0*PI)
                      ENDDO
                    ENDIF !theta coordinate
                  ENDIF !prolate spheroidal coordinates
                ENDDO !nj
C GBS  1-Nov-2000 Adding update from fit (before host mesh is updated)
Cstart
              ELSE IF(FIT) THEN !uses the fitted host mesh geometry
                
                CALL UPSLAVE(IBT,IDO,INP,NBJ,ne1,NEELEM,NKJE,np2,NPF,
     '            NPNE,NRE,nr1,nr2,NVJE,NVJP,SE,XIP,XP,ERROR,*9999)
                
C PM 05-11-02 : Computation in the commented secn moved to UPSLAVE
C               routine. first compute new nodal position
C               CALL XPXE(NBJ(1,ne1),NKJE(1,1,1,ne1),NPF(1,1),
C               '          NPNE(1,1,ne1),nr1,NVJE(1,1,1,ne1),
C               '          SE(1,1,ne1),XA,XE,XP,ERROR,*9999)
C               DO njj=1,NJ_LOC(NJL_FIEL,0,NRE(ne1))
C               nj_field=NJ_LOC(NJL_FIEL,njj,NRE(ne1))
C               nj=NJ_LOC(NJL_GEOM,njj,NRE(ne1)) !corresponding geom var
C               nb1=NBJ(nj,ne1) !region 1 basis
C               XP(1,1,nj,np2)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),
C               '            INP(1,1,nb1),nb1,1,XIP(1,np2),XE(1,
C               nj_field))
C               IF (ITYP10(nr1).EQ.4)THEN !prolate spheroidal coordinates
C ensures the theta coordinates are between zero and 2pi
C               IF(nj.EQ.NJ_LOC(NJL_FIEL,3,nr2)) THEN !theta coordinate
C               DO WHILE (XP(1,1,nj,np2).GT.(2.0d0*PI))
C               XP(1,1,nj,np2)=XP(1,1,nj,np2)-(2.0d0*PI)
C               ENDDO
C                  ENDIF!theta coordinate
C                ENDIF!prolate spheroidal coordinates
C              ENDDO !nj

C then compute new nodal derivatives
C UNDEF is undeformed (dX/dxi)
C DEF   is deformed   (dx/dxi)
C TRANS is transform  (dX/dx)
C              nj_max=NJ_LOC(NJL_GEOM,0,NRE(ne1))
C              CALL ASSERT(nj_max.EQ.NIT(NBJ(1,ne1)),
C     '          'Must have NIT=NJT in host mesh',ERROR,*9999)
C              DO njj=1,nj_max
C                nj_field=NJ_LOC(NJL_FIEL,njj,NRE(ne1))
C                nj=NJ_LOC(NJL_GEOM,njj,NRE(ne1))
C                nb1=NBJ(nj,ne1) !region 1 basis
C                DO ni=1,nj_max
C                  nu=ni*2
C                  IF(nu.EQ.6) nu=7 ! nu=2,4,7
C                  TEMP(njj,ni)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),
C     '              INP(1,1,nb1),nb1,nu,XIP(1,np2),XE(1,nj))
C                  DEF(njj,ni)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),
C     '              INP(1,1,nb1),nb1,nu,XIP(1,np2),XE(1,nj_field))
C                ENDDO !ni
C              ENDDO !njj
C              CALL INVERT(nj_max,TEMP,UNDEF,DET)
C              DO nj=1,nj_max
C                DO njj=1,nj_max
C                  TRANS(nj,njj)=0.0d0
C                  DO nk=1,3
C                    TRANS(nj,njj)=TRANS(nj,njj)+DEF(nj,nk)*UNDEF(nk,njj)
C                  ENDDO
C                ENDDO
C              ENDDO
C              DO nj=1,nj_max
C                DO nk=2,NKJ(nj,np2)+1
C                  OLD(nj,nk)=XP(nk,1,nj,np2)
C                ENDDO
C              ENDDO
C              DO nj=1,nj_max
C                DO nk=2,NKJ(nj,np2)+1
C                  XP(nk,1,nj,np2)=0.0d0
C                  DO njj=1,nj_max
C                    XP(nk,1,nj,np2)=XP(nk,1,nj,np2)
C     '                +TRANS(nj,njj)*OLD(njj,nk)
C                  ENDDO
C                ENDDO
C              ENDDO
Cend
              ELSE !uses the undeformed host mesh geometry
                CALL XPXE(NBJ(1,ne1),NKJE(1,1,1,ne1),NPF(1,1),
     '            NPNE(1,1,ne1),nr1,NVJE(1,1,1,ne1),
     '            SE(1,1,ne1),XA,XE,XP,ERROR,*9999)
                DO nj=1,NJ_LOC(NJL_GEOM,0,NRE(ne1))
                  nb1=NBJ(nj,ne1) !region 1 basis
                  DO nv=1,NVJP(nj,np2)
                    XP(1,nv,nj,np2)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),
     '                INP(1,1,nb1),nb1,1,XIP(1,np2),XE(1,nj))
                    IF (ITYP10(nr1).EQ.4)THEN !prolate spheroidal coordinates
C ensures the theta coordinates are between zero and 2pi
                      IF(nj.EQ.NJ_LOC(NJL_GEOM,3,nr2)) THEN !theta coordinate
                        DO WHILE (XP(1,nv,nj,np2).GT.(2.0d0*PI))
                          XP(1,nv,nj,np2)=XP(1,nv,nj,np2)-(2.0d0*PI)
                        ENDDO
                      ENDIF !theta coordinate
                    ENDIF !prolate spheroidal coordinates
                  ENDDO !nv
                ENDDO !nj
              ENDIF !deformed/undeformed/fit
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C               Critical section is not essential.
CC$            call mp_setlock()
                WRITE(OP_STRING,'('' Node '',I3,'
     '            //'/'' Xi-coords:'',2E12.3,/'' Xj-coords:'',2E12.3)')
     '            np2,(XIP(ni,nonode),ni=1,2),(XP(1,1,nj,np2),nj=1,2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
              ENDIF
            ENDIF !ne1.GT.0
          ENDDO !nonode
          IF (ITYP10(nr1).EQ.4)THEN ! prolate spheroidal coordinates
C this section of code maintains consistencey of theta with elemets with
C the corresponding xi direction of the host mesh by recalculating the
C connectivity.
            DO noelem=1,NEELEM(0,nr2)
              ne2=NEELEM(noelem,nr2)
              nb2=NBJ(1,ne2)
              IF (NIT(nb2).EQ.1)THEN
                np1=NPNE(1,nb2,ne2)
                np2=NPNE(2,nb2,ne2)
                nj=NJ_LOC(NJL_GEOM,3,nr2)
                IF((DABS(XP(1,1,nj,np1)-XP(1,1,nj,np2)))
     '            .LT.(ZERO_TOL*2.0d0)) THEN
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np1,.FALSE.,ERROR,*9999)
                  XP(1,1,nj,np1)=XP(1,1,nj,np1)+(ZERO_TOL*2.0d0) !
C checks theta angles are not equal
                ENDIF
                IF(DABS(XP(1,1,nj,np2)-XP(1,1,nj,np1)).LT.PI) THEN
                  IF(XP(1,1,nj,np1).LT.XP(1,1,nj,np2)) THEN
                    NPNE(1,nb2,ne2)=np2
                    NPNE(2,nb2,ne2)=np1
C swaps locals nodes such that theta is decreasing in xi1 to be
C consistant with cordinate system
                  ELSE
                    NPNE(1,nb2,ne2)=np1
                    NPNE(2,nb2,ne2)=np2
                  ENDIF
                ELSE
                  IF(XP(1,1,nj,np1).LT.XP(1,1,nj,np2)) THEN
                    NPNE(1,nb2,ne2)=np1
                    NPNE(2,nb2,ne2)=np2
C swaps locals nodes such that theta is decreasing in x1 to be
C consistant with cordinate system
                  ELSE
                    NPNE(1,nb2,ne2)=np2
                    NPNE(2,nb2,ne2)=np1
                  ENDIF
                ENDIF ! checks if nodes cross theta equals zero
              ENDIF
            ENDDO
          ELSE IF(ITYP10(nr1).EQ.2)THEN  !cylindrical coords
            DO noelem=1,NEELEM(0,nr2)
              ne2=NEELEM(noelem,nr2)
              nb2=NBJ(1,ne2)
              IF (NIT(nb2).EQ.1)THEN !one dimensional elements
                np1=NPNE(1,nb2,ne2)
                np2=NPNE(2,nb2,ne2)
                nj=NJ_LOC(NJL_GEOM,2,nr2)
                IF((DABS(XP(1,1,nj,np1)-XP(1,1,nj,np2)))
     '            .LT.(ZERO_TOL*2.0d0)) THEN
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np1,.FALSE.,ERROR,*9999)
                  XP(1,1,nj,np1)=XP(1,1,nj,np1)+(ZERO_TOL*2.0d0) !
C checks theat angles are not equal
                ENDIF
                IF(DABS(XP(1,1,nj,np1)-XP(1,1,nj,np2)).LT.PI) THEN
                  IF(XP(1,1,nj,np1).LT.XP(1,1,nj,np2)) THEN
                    NPNE(1,nb2,ne2)=np1
                    NPNE(2,nb2,ne2)=np2
C swaps locals nodes such that theta is decreasing in x1 to be
C consistant with cordinate system
                  ELSE
                    NPNE(1,nb2,ne2)=np2
                    NPNE(2,nb2,ne2)=np1
                  ENDIF
                ELSE
                  IF(XP(1,1,nj,np2).LT.XP(1,1,nj,np1)) THEN
                    NPNE(1,nb2,ne2)=np2
                    NPNE(2,nb2,ne2)=np1
C swaps locals nodes such that theta is decreasing in x1 to be
C consistant with cordinate system
                  ELSE
                    NPNE(1,nb2,ne2)=np1
                    NPNE(2,nb2,ne2)=np2
                  ENDIF
                ENDIF ! checks if nodes cross theta equals zero
              ENDIF !NIT(nb2).EQ.1-one dimensional elements
            ENDDO !noelem
          ENDIF !cylindrical coords
        ELSE IF(UPDATE_GEOMETRY)THEN
          CALL UPMESH_GEOMETRY(NBJ,NEELEM(0,nr),NELIST,ne_start,NORD,
     &      NPNE,NVJE,NXI,AA_RATIO,ACINUS_LENGTH,ACINUS_VOLUME,CE,
     &      ld_ratio,RATIO_DIAMETER,scale_f,value,XAB,XP,CUBE_ROOT,
     &      UNSTRAINED,RADIUS_SCHEME,SCALED,ERROR,*9999)
        ELSE IF(CALL_UPVIEW) THEN !postion
          CALL UPVIEW(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,ISFIEL,ISLINE,
     '      ISLINO,ISNONO,IWK,MXI,NAN,NBH,NBJ,NEELEM,NGAP,NHE(1,nx),
     '      NHP(1,0,nx),NKH,NKHE,NKJE,NLATNE,NLL,NLLIST,NPF,NPL,NPNE,
     '      NPNODE,NQNE,NQNLAT,NQS,NQXI,NRE,NTIW,NVHE,NVHP,NVJE,
     &      NW(1,1,nx),nx,NYNE,NYNP,CURVCORRECT,DL,SE,XA,XE,XG,XP,XQ,YG,
     '      YP(1,1,nx),YQ,ZA,ZE,ZP,CSEG,FIX(1,1,nx),ERROR,*9999)
C         CALL SEGPAC
        ENDIF !update_position
      ENDIF

      CALL EXITS('UPMESH')
      RETURN
 9999 CALL ERRORS('UPMESH',ERROR)
      CALL EXITS('UPMESH')
      RETURN 1
      END


