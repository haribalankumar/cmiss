      SUBROUTINE OPSTRE(IBT,IDO,INP,IPOINTTYP,IXI,NAN,NBH,NBJ,
     '  NDDL,NDLT,NEELEM,NELIST,
     '  NGLIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,NPNODE,nr,NTPOIN,NVHE,
     '  NVHP,NVJE,NW,nx,nxc,NYNE,NYNP,
     '  CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RG,SE,WG,XA,XE,
     '  XG,XID,XIG,XIPOS,XP,YG,YP,ZA,ZE,ZG,ZP,COORDS,STRESSTYPE,
     '  EXTREMA_ONLY,FULL,STRAINENERGY,ERROR,*)

C#### Subroutine: OPSTRE
C###  Description:
C###    OPSTRE outputs stress tensors & principal stresses at Gauss pts
C###    or at equally spaced points along an arbitrary xi coordinate
C###    line.

C**** PRSTMAX is principle stress maximum.
C**** PRSTMIN is principle stress minimum (absolute value).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grow00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'stra00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IPOINTTYP,
     '  IXI,NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NDDL(NEM,NDEM),NDLT(NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NGLIST(0:NGM),NHE(NEM),NHP(NPM),
     '  NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),
     '  nr,NTPOIN,NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),nx,nxc,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XIG(NIM,NGM,NBM),XID(NIM,NDM),
     '  XIPOS(3),XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER COORDS*(*),ERROR*(*),STRESSTYPE*(*)
      LOGICAL EXTREMA_ONLY,FULL,STRAINENERGY
!     Local Variables
      INTEGER cpcounter,IBEG,IEND,ig,INTWORK(1),ISEG(1),nb,
     '  nd,nde,ne,ng,NGLIST_TEMP(0:NGM),nolist,nopoin,i,j
      REAL*8 ENERGY_ELEMENT,SE_GAUSS,se_dev,se_vol,
     '  ENERGY_TOTAL,PHI(3),PST(3),REALWORK(1),
     '  RG2D,RGZ,RGZ2D,RM(3,3),TC(3,3),TG(3,3),TN(3,3),
     '  TNA,DZDX(3,3),EG(3,3),R(3,3),RI1,RI2,RI3,U(3,3),
     '  C(3,3),detC,Ed(3,3),Ev(3,3)
!     Functions
      REAL*8 DET

      CHARACTER COORDINATES*9,CSEG(1)*(1),STRING*(MXCH)
      LOGICAL END

      CALL ENTERS('OPSTRE',*9999)

      DO i=1,3
        DO j=1,3
          DZDX(i,j)=0.0d0
          EG(i,j)=0.0d0
          U(i,j)=0.0d0
          C(i,j)=0.0d0
          Ed(i,j)=0.0d0
          Ev(i,j)=0.0d0
        ENDDO
      ENDDO

      COORDINATES=' ' !previously not set MLB 4-Dec-97

      WRITE(OP_STRING,'(/'' Class '',I1,'' (nx='',I1,'
     '  //''') Region '',I1,'':'')') nxc,nx,nr
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '  nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)

      IF(ITYP1(nr,nx).EQ.5) THEN !finite elasticity
        IF(TABLE.AND..NOT.EXTREMA_ONLY) THEN
          CALL STRING_TRIM(COORDS,IBEG,IEND)
          FORMAT='('' '//COORDS(IBEG:IEND)//' stresses: '//
     '      'np,((TC(mi,ni),mi=1,ni),ni=1,NIT),'//
     '      '(PST(ni),ni=1,NIT),(Phi(ni),ni=1,NIT):'')'
          WRITE(OP_STRING,FORMAT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='('' ne='',30(I5,'',''),:/(4X,30(I5,'','')))'
          WRITE(OP_STRING,FORMAT) (NELIST(nolist),nolist=1,NELIST(0))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      ENERGY_TOTAL=0.0d0
      DO nolist=1,NELIST(0)
        ne=NELIST(nolist)
        nb=NBH(NH_LOC(1,nx),1,ne)
        ENERGY_ELEMENT=0.0d0
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '    NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '    CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '    ERROR,*9999)

        IF(ITYP1(nr,nx).EQ.5) THEN !finite elasticity
          CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
     '      CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
        ELSE !all other problem types
          CALL CPCG(NW(ne,1),NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
     '      CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
          IF(KTYP60.EQ.1) THEN !Growth law: so get density from YG(5)
            DO ng=1,NGT(NBJ(1,ne))
              CG(5,ng)=YG(5,ng,ne) !density
              CG(1,ng)=CG(1,ng)*CG(5,ng)**2 !Young's mod from density
            ENDDO !ng
          ELSE IF(KTYP60.EQ.2) THEN !Growth law: get density from YG(5)
            DO ng=1,NGT(NBJ(1,ne))
              CG(5,ng)=YG(5,ng,ne) !density
              CG(1,ng)=CG(1,ng)*(CG(5,ng)/GROW1)**2 !Y's mod frm density
            ENDDO !ng
          ENDIF
        ENDIF

        IF(ITYP1(nr,nx).EQ.4) THEN !linear elasticity
          CALL OPST40(nb,NBH(1,1,ne),NBJ(1,ne),ne,NHE(ne),
     '      NW(ne,1),nx,CG,ENERGY_ELEMENT,PG,'WRITES',
     '      WG,XE,XG,YG(1,1,ne),ZE,ZG,.TRUE.,
     '      COORDINATES,ERROR,*9999)
          ENERGY_TOTAL=ENERGY_TOTAL+ENERGY_ELEMENT

        ELSE IF(ITYP1(nr,nx).EQ.5) THEN !finite elasticity
          IF(IXI.EQ.0) THEN ! all xi dirs

C MPN 5Aug2014: adding stress output at data pts
            IF(IPOINTTYP.EQ.2) THEN !stress at data points
              CALL ASSERT(.NOT.STRAINENERGY,
     '          '>>Strain energy at data pts not implemented',
     '          ERROR,*9999)

              IF(.NOT.TABLE.AND..NOT.EXTREMA_ONLY.AND.
     '          NDLT(ne).GE.1) THEN
                CALL STRING_TRIM(COORDS,IBEG,IEND)
                FORMAT='(/'' Element'',I3,'' Data point stresses'
     '            //' wrt '//COORDS(IBEG:IEND)//' coords:'')'
                WRITE(OP_STRING,FORMAT) ne
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF

              DO nde=1,NDLT(ne)
                nd=NDDL(ne,nde)

            IF(KTYP54(nr).EQ.3) THEN !Gauss point stresses (grid coupling)
C             NGLIST_TEMP holds NGLIST values as NGLIST gets reset to
C             1,6 after
C             UPGRID function has been executed
              DO cpcounter=0,NGLIST(0)
                NGLIST_TEMP(cpcounter)=NGLIST(cpcounter)
              ENDDO !cpcounter
C NOTE: The arrays to store the green strain have been hardcoded to be RCQS
C       The indices of RCQS are also hardcoded. If these numbers and array change
C       a way of parsing perl variables into CMISS code must be coded
C
              CO(1)='FEM'
              CO(2)='UPDATE'
              CO(3)='GRID'
              CO(4)='GREEN_STRAIN'
              CO(5)='no_ze_calc'
              CO(6)='COMPONENT'
              CO(7)='1,2,3,4,5,6'
              CO(8)='RCQS'
              CO(9)='3,6,8,4,5,7'
              CO(10)='ALL_VARIANTS'
              CO(11)='ELEMENTS'
              WRITE(CO(12),'(I5)') ne
              CO(13)='REGION'
              WRITE(CO(14),'(I5)') nr
              CO(15)='ATXI'
              CO(16)='XI_1'
              WRITE(CO(17),'(F19.15)') XID(1,nd)
              CO(18)='XI_2'
              WRITE(CO(19),'(F19.15)') XID(2,nd)
              CO(20)='XI_3'
              WRITE(CO(21),'(F19.15)') XID(3,nd)
              NTCO=21
              STRING='FEM update grid green_strain no_ze_calc '
     '          //'comp 1,2,3,4,5,6 RCQS 3,6,8,4,5,7 ALL_VARIANTS '
     '          //'element '//CO(12)//' region '//CO(14)
     '          //' atxi xi_1='//CO(17)//' xi_2='//CO(19)
     '          //' xi_3='//CO(21)
              CALL UPGRID(IBT,IDO,INP,%VAL(ICQS_SPATIAL_PTR),
     '          %VAL(IRCQS_SPATIAL_PTR),NAN,%VAL(NAQ_PTR),
     '          %VAL(NBH_PTR),%VAL(NBJ_PTR),%VAL(NEELEM_PTR),
     '          %VAL(NELIST_PTR),%VAL(NENQ_PTR),%VAL(NGAP_PTR),
     '          %VAL(NHE_PTR),%VAL(NHP_PTR),%VAL(NKH_PTR),
     '          %VAL(NKHE_PTR),%VAL(NKJE_PTR),%VAL(NLL_PTR),
     '          %VAL(NLQ_PTR),%VAL(NPF_PTR),%VAL(NPL_PTR),NPNE,
     '          %VAL(NPNODE_PTR),%VAL(NQGP_PTR),%VAL(NQLIST_PTR),
     '          %VAL(NQNE_PTR),%VAL(NQXI_PTR),%VAL(NQS_PTR),
     '          %VAL(NRLIST_PTR),%VAL(NVHE_PTR),
     '          %VAL(NVHP_PTR),%VAL(NVJE_PTR),%VAL(NW_PTR),
     '          %VAL(NWQ_PTR),%VAL(NXLIST_PTR),
     '          %VAL(NXQ_PTR),%VAL(NYNE_PTR),%VAL(NYNP_PTR),
     '          %VAL(AQ_PTR),%VAL(CE_PTR),%VAL(CP_PTR),%VAL(CQ_PTR),
     '          %VAL(CURVCORRECT_PTR),%VAL(DL_PTR),%VAL(DNUDXQ_PTR),
     '          %VAL(DXDXIQ_PTR),%VAL(DXDXIQ2_PTR),
     '          %VAL(FEXT_PTR),%VAL(GCHQ_PTR),
     '          %VAL(GUQ_PTR),PG,%VAL(PROPQ_PTR),
     '          %VAL(RCQS_SPATIAL_PTR),SE,%VAL(XA_PTR),
     '          XE,XG,%VAL(XIQ_PTR),
     '          %VAL(XP_PTR),%VAL(XQ_PTR),%VAL(YG_PTR),%VAL(YP_PTR),
     '          %VAL(YQ_PTR),%VAL(YQS_PTR),%VAL(ZA_PTR),
     '          ZE,ZG,%VAL(ZP_PTR),
     '          STRING,ERROR,*9999)
C NOTE: Class is hardcoded to be class 2. If class changes
C       code must be written to identify class
              CO(1)='FEM'
              CO(2)='SOLVE'
              CO(3)='CLASS'
              CO(4)='2'
              CO(5)='RESTART'
              CO(6)='TO'
              CO(7)='0'
              NTCO=7
              STRING='FEM solve class 2 restart to 0'
              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '          ERROR,*9999)
C             Calling UPGAUS through FEM.f to update gauss array with
C             stresses calculated via grid coupling
C NOTE:       Stresses at grid points currently YQS array 1,2,3,4,5,6
              CO(1)='FEM'
              CO(2)='UPDATE'
              CO(3)='GAUSS'
              CO(4)='GRIDVARS'
              CO(5)='YQS'
              CO(6)='2,3,4,5,6,7'
              CO(7)='YG'
              CO(8)='1,2,3,4,5,6'
              CO(9)='INCLUDE'
              CO(10)='ELEMENT'
              WRITE(CO(11),'(I5)') ne
              CO(12)='REGION'
              WRITE(CO(13),'(I5)') nr
              NTCO=13
              STRING='FEM update gauss gridvars yqs 2,3,4,5,6,7 '
     '          //'yg 1,2,3,4,5,6 include element '//CO(11)
     '          //' region '//CO(13)

C!!!!!! MAY NEED TO UPDATE THIS TO ONLY SET YG FOR FIRST GRID POINT??
C!!!!!! (default is to do it for all grid points in the element, but 
C!!!!!!  we only need the first here).

              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '          ERROR,*9999)
C             copy NGLIST_TEMP back to NGLIST
              DO cpcounter=0,NGLIST(0)
                NGLIST(cpcounter)=NGLIST_TEMP(cpcounter)
              ENDDO !cpcounter
            ENDIF  !KTYP54(nr).EQ.3 Gauss point stresses (grid coupling)

                ng=1 !a hack for RG, FEXT, YG below!!?
                CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,
     '            IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne),
     '            NPNE(1,1,ne),nr,ne,nx,
     '            CE(1,ne),CG,CP,FEXT(1,ng,ne),PG,PHI,PST,RG(ng),RG2D,
     '            RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XID(1,nd),
     '            YG(1,ng,ne),ZE,ZG,ERROR,*9999)
C                 (For cellml might only need YG(*,1,ne) above?)
                IF(.NOT.EXTREMA_ONLY) THEN
                  CALL OPTG50(COORDS,IPOINTTYP,STRESSTYPE,
     '              NBH(1,1,ne),nd,nr,nx,
     '              PHI,PST,RM,TC,TG,TN,TNA,
     '              XG,XID(1,nd),ZG,FULL,ERROR,*9999)
                ENDIF
              ENDDO !nde data point in element loop

            ELSE !IPOINTTYP.EQ.1  i.e. stress at Gauss pts

              IF(.NOT.TABLE.AND..NOT.EXTREMA_ONLY) THEN
                CALL STRING_TRIM(COORDS,IBEG,IEND)
                FORMAT='(/'' Element'',I3,'' Gauss point stresses'
     '            //' wrt '//COORDS(IBEG:IEND)//' coords:'')'
                WRITE(OP_STRING,FORMAT) ne
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
C news VJ 13Jan2004.
            IF(KTYP54(nr).EQ.3) THEN !Gauss point stresses (grid coupling)
C             NGLIST_TEMP holds NGLIST values as NGLIST gets reset to 1,6 after
C             UPGRID function has been executed
              DO cpcounter=0,NGLIST(0)
                NGLIST_TEMP(cpcounter)=NGLIST(cpcounter)
              ENDDO !cpcounter
C news VJ 17Jan2004: Calling UPGRID to update grid points with green
C                    strain components
C       E11
C NOTE: The arrays to store the green strain have been hardcoded to be RCQS
C       The indices of RCQS are also hardcoded. If these numbers and array change
C       a way of parsing perl variables into CMISS code must be coded
C
              CO(1)='FEM'
              CO(2)='UPDATE'
              CO(3)='GRID'
              CO(4)='GREEN_STRAIN'
              CO(5)='no_ze_calc'
              CO(6)='COMPONENT'
              CO(7)='1,2,3,4,5,6'
              CO(8)='RCQS'
              CO(9)='3,6,8,4,5,7'
              CO(10)='ALL_VARIANTS'
              CO(11)='ELEMENTS'
              WRITE(CO(12),'(I5)') ne
              CO(13)='REGION'
              WRITE(CO(14),'(I5)') nr
              NTCO=14
C MPN 5Aug2014 - this looks like a bug? STRING did not match CO etc
              STRING='FEM update grid green_strain no_ze_calc '
     '          //'comp 1,2,3,4,5,6 RCQS 3,6,8,4,5,7 ALL_VARIANTS '
     '          //'element '//CO(12)//' region '//CO(14)
C old              STRING='FEM update grid green_strain no_ze_calc 
C old     '          comp 1,2,3,4,5,6 RCQS 3,6,8,4,5,6 ALL_VARIANTS ELEM '
              CALL UPGRID(IBT,IDO,INP,%VAL(ICQS_SPATIAL_PTR),
     '          %VAL(IRCQS_SPATIAL_PTR),NAN,%VAL(NAQ_PTR),
     '          %VAL(NBH_PTR),%VAL(NBJ_PTR),%VAL(NEELEM_PTR),
     '          %VAL(NELIST_PTR),%VAL(NENQ_PTR),%VAL(NGAP_PTR),
     '          %VAL(NHE_PTR),%VAL(NHP_PTR),%VAL(NKH_PTR),
     '          %VAL(NKHE_PTR),%VAL(NKJE_PTR),%VAL(NLL_PTR),
     '          %VAL(NLQ_PTR),%VAL(NPF_PTR),%VAL(NPL_PTR),NPNE,
     '          %VAL(NPNODE_PTR),%VAL(NQGP_PTR),%VAL(NQLIST_PTR),
     '          %VAL(NQNE_PTR),%VAL(NQXI_PTR),%VAL(NQS_PTR),
     '          %VAL(NRLIST_PTR),%VAL(NVHE_PTR),
     '          %VAL(NVHP_PTR),%VAL(NVJE_PTR),%VAL(NW_PTR),
     '          %VAL(NWQ_PTR),%VAL(NXLIST_PTR),
     '          %VAL(NXQ_PTR),%VAL(NYNE_PTR),%VAL(NYNP_PTR),
     '          %VAL(AQ_PTR),%VAL(CE_PTR),%VAL(CP_PTR),%VAL(CQ_PTR),
     '          %VAL(CURVCORRECT_PTR),%VAL(DL_PTR),%VAL(DNUDXQ_PTR),
     '          %VAL(DXDXIQ_PTR),%VAL(DXDXIQ2_PTR),
     '          %VAL(FEXT_PTR),%VAL(GCHQ_PTR),
     '          %VAL(GUQ_PTR),PG,%VAL(PROPQ_PTR),
     '          %VAL(RCQS_SPATIAL_PTR),SE,%VAL(XA_PTR),
     '          XE,XG,%VAL(XIQ_PTR),
     '          %VAL(XP_PTR),%VAL(XQ_PTR),%VAL(YG_PTR),%VAL(YP_PTR),
     '          %VAL(YQ_PTR),%VAL(YQS_PTR),%VAL(ZA_PTR),
     '          ZE,ZG,%VAL(ZP_PTR),
     '          STRING,ERROR,*9999)
C news 24Jan2004: Calling solve to evaluate cell through FEM.f
C NOTE: Class is hardcoded to be class 2. If class changes
C       code must be written to identify class
              CO(1)='FEM'
              CO(2)='SOLVE'
              CO(3)='CLASS'
              CO(4)='2'
              CO(5)='RESTART'
              CO(6)='TO'
              CO(7)='0'
              NTCO=7
              STRING='FEM solve class 2 restart to 0'
              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '          ERROR,*9999)
C             Calling UPGAUS through FEM.f to update gauss array with
C             stresses calculated via grid coupling
C NOTE:       Stresses at grid points currently YQS array 1,2,3,4,5,6
              CO(1)='FEM'
              CO(2)='UPDATE'
              CO(3)='GAUSS'
              CO(4)='GRIDVARS'
              CO(5)='YQS'
              CO(6)='2,3,4,5,6,7'
              CO(7)='YG'
              CO(8)='1,2,3,4,5,6'
              CO(9)='INCLUDE'
              CO(10)='ELEMENT'
              WRITE(CO(11),'(I5)') ne
              CO(12)='REGION'
              WRITE(CO(13),'(I5)') nr
              NTCO=13
C MPN 5Aug2014 - this looks like a bug? STRING did not match CO etc
              STRING='FEM update gauss gridvars yqs 2,3,4,5,6,7 '
     '          //'yg 1,2,3,4,5,6 include element '//CO(11)
     '          //' region '//CO(13)
C old              STRING='FEM update gauss gridvars yqs 2 yg 1'
              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '          ERROR,*9999)
C newe VJ
C             copy NGLIST_TEMP back to NGLIST
              DO cpcounter=0,NGLIST(0)
                NGLIST(cpcounter)=NGLIST_TEMP(cpcounter)
              ENDDO !cpcounter
            ENDIF !KTYP54(nr).EQ.3 Gauss point stresses (grid coupling)

CC newe VJ 13Jan2004

            DO ig=1,NGLIST(0)
              ng=NGLIST(ig)
              CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,
     '          IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),ng,NHE(ne),
     '          NPNE(1,1,ne),nr,ne,nx,
     '          CE(1,ne),CG,CP,FEXT(1,ng,ne),PG,PHI,PST,RG(ng),RG2D,
     '          RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XIG(1,ng,nb),
     '          YG(1,ng,ne),ZE,ZG,ERROR,*9999)
              IF(.NOT.EXTREMA_ONLY) THEN
                CALL OPTG50(COORDS,IPOINTTYP,STRESSTYPE,
     '            NBH(1,1,ne),ng,nr,nx,
     '            PHI,PST,RM,TC,TG,TN,TNA,
     '            XG,XIG(1,ng,nb),ZG,FULL,ERROR,*9999)
              ENDIF
C news GR 08August2008
              IF(STRAINENERGY) THEN
                ! calculate Green's strain
                CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '            NBJ(1,ne),ng,NHE(ne),NPNE(1,1,ne),nr,nx,
     '            DZDX,CE(1,ne),CG,CP,EG,PG,PHI,PST,
     '            R,RG(ng),RI1,RI2,RI3,RM,U,
     '            XE,XG,XIG(1,ng,nb),ZE,ZG,ERROR,*9999)

                ! calculate the deviatoric component of Green's strain
                ! C = 2*E+I, Cd = det(C)**(-1/3)*C, Ed = (1/2)*(Cd-I)
                ! volumetric component is EG-Ed
                DO i=1,3
                  DO j=1,3
                    C(i,j)=2*EG(i,j)
                    IF(i.EQ.j) C(i,j)=C(i,j)+1
                  ENDDO
                ENDDO
                detC=DET(C)**(-1.0D0/3.0D0)
      !PRINT *,'OPSTRE detC**(-1/3) ',detC
                DO i=1,3
                  DO j=1,3
                    Ed(i,j)=0.5d0*(detC*C(i,j))
                    IF(i.EQ.j) Ed(i,j)=Ed(i,j)-0.5d0
                    Ev(i,j)=EG(i,j)-Ed(i,j)
                  ENDDO
                ENDDO
!       PRINT *,'OPSTRE EG ',((EG(i,j),i=1,3),j=1,3)
!       PRINT *,'OPSTRE Ed ',((Ed(i,j),i=1,3),j=1,3)
!       PRINT *,'OPSTRE Ev ',((Ev(i,j),i=1,3),j=1,3)

                FORMAT='(7X,''Strain energy density:'')'
                WRITE(OP_STRING,FORMAT)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
 
                CALL STRAIN_ENERGY(nr,ne,ng,Ed,se_dev,ERROR,*9999)
                CALL STRAIN_ENERGY(nr,ne,ng,Ev,se_vol,ERROR,*9999)
                CALL STRAIN_ENERGY(nr,ne,ng,EG,SE_GAUSS,ERROR,*9999)

                FORMAT='(7X,''Deviatoric:'',D12.4,'' Volumetric:'''
     &            //',D12.4,'' Total:'',D12.4 )'
                WRITE(OP_STRING,FORMAT) se_dev,se_vol,SE_GAUSS
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

!               PRINT *,'OPSTRE dev vol tot ', se_dev, se_vol, SE_GAUSS

                ! Integrate over undeformed volume since we are using
                ! 2PK stress & Green strain
                se_dev=se_dev*RG(ng)*WG(ng,nb)
                se_vol=se_vol*RG(ng)*WG(ng,nb)
                SE_GAUSS=SE_GAUSS*RG(ng)*WG(ng,nb)
                ENERGY_ELEMENT=ENERGY_ELEMENT+SE_GAUSS
                FORMAT='(7X,''Strain energy:'')'
                WRITE(OP_STRING,FORMAT)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                FORMAT='(7X,''Deviatoric:'',D12.4,'' Volumetric:'''
     &            //',D12.4,'' Total:'',D12.4 )'
                WRITE(OP_STRING,FORMAT) se_dev,se_vol,SE_GAUSS
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF !STRAINENERGY
C newe GR 08August08
            ENDDO !ig
C news GR 08August2008
            IF(STRAINENERGY) THEN
              ENERGY_TOTAL=ENERGY_TOTAL+ENERGY_ELEMENT
            ENDIF
C newe GR 08August08

            ENDIF !IPOINTTYP.EQ.2  (stress at data pts)

          ELSE IF(IXI.GE.1 .AND. IXI.LE.3) THEN
            IF(.NOT.TABLE.AND..NOT.EXTREMA_ONLY) THEN
              FORMAT='(/'' Element'',I3,'' stress solutions:'')'
              WRITE(OP_STRING,FORMAT) ne
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            DO nopoin=1,NTPOIN
              XIPOS(IXI)=DBLE(nopoin-1)/DBLE(NTPOIN-1)
C!!! WARNING: FEXT and YG are not passed correctly: should find nearest
C!!!          Gauss pt to Xi point and pass eg. YG(1,ng,ne).
              CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,
     '          IBT,IDO,INP,NAN,NBH(1,1,ne),
     '          NBJ(1,ne),0,NHE(ne),NPNE(1,1,ne),nr,ne,nx,
     '          CE(1,ne),CG,CP,FEXT(1,1,ne),PG,PHI,PST,RG(1),RG2D,
     '          RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XIPOS,
     '          YG(1,1,ne),ZE,ZG,ERROR,*9999)
              IF(.NOT.EXTREMA_ONLY) THEN
                CALL OPTG50(COORDS,IPOINTTYP,STRESSTYPE,
     '            NBH(1,1,ne),nopoin,nr,nx,
     '            PHI,PST,RM,TC,TG,TN,TNA,XG,XIPOS,ZG,FULL,ERROR,*9999)
              ENDIF
            ENDDO !nopoin
          ENDIF
        ELSEIF(ITYP1(nr,nx).EQ.9.AND.ITYP2(nr,nx).EQ.1) THEN
          !BEM and linear elasticity
          CALL OPST40(nb,NBH(1,1,ne),NBJ(1,ne),ne,NHE(ne),
     '      NW(ne,1),nx,CG,ENERGY_ELEMENT,PG,'WRITES',
     '      WG,XE,XG,YG(1,1,ne),ZE,ZG,.TRUE.,
     '      COORDINATES,ERROR,*9999)
          ! GR ENERGY_ELEMENT is the sum of the strain energy densities at each
          ! gauss point. It has not been integrated over the element volume, so
          ! is probably not correct for some problem types.
          ENERGY_TOTAL=ENERGY_TOTAL+ENERGY_ELEMENT
        ENDIF

      ENDDO !nolist (ne)

      IF((ITYP1(nr,nx).EQ.4 !linear elasticity
     '  .OR.ITYP1(nr,nx).EQ.9.AND.ITYP2(nr,nx).EQ.1) !BEM+lin elast
     '  .OR.STRAINENERGY.EQV..TRUE.) THEN !finite elasticity strain energy
        WRITE(OP_STRING,'(/'' Total strain energy '//
     '    'for listed elements = '','//
     '    'D12.4)') ENERGY_TOTAL
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(FULL) THEN
        WRITE(OP_STRING,'(/'' Min princ stress for listed points = '','
     '    //'D12.4)') PRSTMIN
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Max princ stress for listed points = '','
     '    //'D12.4)') PRSTMAX
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('OPSTRE')
      RETURN
 9999 CALL ERRORS('OPSTRE',ERROR)
      CALL EXITS('OPSTRE')
      RETURN 1
      END


