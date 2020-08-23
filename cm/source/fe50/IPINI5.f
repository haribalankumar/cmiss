      SUBROUTINE IPINI5(IBT,IDO,INP,NAN,NBH,NBHF,NEELEM,NEL,NENP,NFF,
     '  NHE,NHP,NKEF,NKH,NKHE,NLL,NNF,NNL,NPF,NPL,NPLIST,NPNE,NPNODE,
     '  nr,NVHE,NVHP,NW,nx,NYNE,NYNP,CE,CG,CGE,CP,DF,DL,PG,SE,WG,
     '  XE,XG,XP,YG,YP,ZA,ZE,ZG,ZP,FIX,NOFIX,UPDATE_REF,ERROR,*)

C#### Subroutine: IPINI5
C###  Description:
C###    Inputs initial conditions and boundary conditions for finite
C###    elasticity problems.

C**** Global node var.s are carried in ZP(nk,nv,nh,np,nc) as follows:
C****   nh=1
C****      :         nodal coordinates of deformed state
C****     NJP(np) ! RGB NJP replaced by NJ_LOC 13/11/97
C****     NJP(np)+1  hyd.press. or lateral ext. (if reqd as a nodal var)
C**** Auxiliary element variables (e.g. press. or lat. ext.) are carried
C****   in ZA(na,nh,nc,ne).

C**** If KTYP57(nr)>1 pressure b.c.'s are applied to Xi(3)=0,1 faces of
C****   elements with NW(ne,1)=2,3, respectively, or both if NW(ne,1)=4.
C****   Pressure b.c.'s are applied using auxiliary variables with
C****   basis fn indices -1 and -2 for Xi3=1 and Xi3=1 faces respec.
C****   When setting up auxiliary basis for pressure, choose 5 aux
C****   params and 'Pressure' aux basis. Then use the -1/-2 indices
C****   to signify which params are for pressure b.c's. To apply
C****   pressure b.c's fix the appropriate parameters in each
C****   element with the desired increments during define initial.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='IPINI5')

      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'mach00.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM), 
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBHF(NHM,NCM,NFM),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),
     '  NHE(NEM),NHP(NPM),NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NLL(12,NEM),NNF(0:17,6,NBFM),
     '  NNL(0:4,12,NBFM),NPF(9,NFM),NPL(5,0:3,NLM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),
     '  CGE(NMM,NGM,NEM),CP(NMM,NPM),DF(NFM),DL(3,NLM),
     '  PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM),NOFIX,UPDATE_REF
!     Local Variables
      INTEGER i,INFO,INTWORK(1),ISEG(1),
     '  na,nb,NBP,ne,nh_pressure,nk,
     '  noelem,nonode,NOQUES,np,nv,temp_i,temp_j
      REAL*8 AZL(3,3),DW(6),HZERO,REALWORK(1)
      CHARACTER CSEG(1)*(1),STRING*(MXCH)
      LOGICAL END,FILEIP,FOUNDna1,FOUNDna2
      EXTERNAL CELLML_DUMMY_ROUTINE

      CALL ENTERS(ROUTINENAME,*9999)

      FILEIP=.FALSE.
      NOQUES=0

      IF(KTYP51(nr).GE.3) THEN !3D,membrane,string or shell
        FORMAT='('' Specify option for Xi(3) faces [1]:'''//
     '    '/''   (1) No pressure boundary conditions'''//
     '    '/''   (2) Boundary pressure increments entered'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP57(nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IONE,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP57(nr)=IDATA(1)

        IF((KTYP57(nr).EQ.2).AND.(KTYP51(nr).EQ.3)) THEN !pressure increments entered
C         Check that the pressure basis function contains parameters for
C         pressure boundary conditions (NBP is press basis for elem 1)
          NBP=NBH(NH_LOC(NH_LOC(0,nx),nx),1,NEELEM(1,nr))
          FOUNDna1=.FALSE.
          FOUNDna2=.FALSE.
          DO na=1,NAT(NBP)
C           Check that pressure varies in Xi(3) dirn only for
C           element based pressure interpolation
            IF(NST(NBP).EQ.0.AND.  !elem based press interp
     '        (NAN(1,na,NBP).NE.0.OR.NAN(2,na,NBP).NE.0)) THEN
              ERROR='>>Pressure must vary only with Xi(3) '
     '          //'for elem based pressure interpolation'
              GOTO 9999
            ENDIF
C           check for element based hyd press interp accounting for
C           pressure boundary condition parameters
            IF(NAN(3,na,NBP).EQ.-1) THEN
              FOUNDna1=.TRUE.
            ELSE IF(NAN(3,na,NBP).EQ.-2) THEN
              FOUNDna2=.TRUE.
            ENDIF
          ENDDO !na
C old MPN 7Apr97: now handled below
C          IF(KTYP51(nr).EQ.3.AND. !3D
C     '      .NOT.(FOUNDna1.AND.FOUNDna2).OR.
C     '      KTYP51(nr).GT.3.AND.  !membrane, string or shell
C     '      .NOT.(FOUNDna1.OR.FOUNDna2)) THEN
C            ERROR='>>Pressure basis fn must have params for press bcs'
C            GO TO 9999
C          ENDIF

C          IF(KTYP51(nr).EQ.3) THEN !3D
C           Prompt for eqn type for undetermined element pressure vars
            FORMAT='('' Determine free hydrostatic pressure vars '
     '        //'using Xi3 face constraints to enforce [1]:'''
     '        //'/''   (1) incompressibility'''
     '        //'/''   (2) continuous normal Cauchy stress'''
     '        //'/$,''    '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP5A(nr)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IONE,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP5A(nr)=IDATA(1)
C new MPN 7Apr97
            IF(KTYP5A(nr).EQ.2.AND. !match norm stress on Xi3 faces
     '        .NOT.(FOUNDna1.AND.FOUNDna2)) THEN
              ERROR='>>ERROR: To match norm stress, press basis fn '
     '          //'must have params for press bcs'
              GO TO 9999
            ENDIF
C          ELSE IF(KTYP51(nr).GT.3) THEN  !membrane, string or shell
C            IF(.NOT.(FOUNDna1.OR.FOUNDna2)) THEN
C              ERROR='>>Pressure basis fn must have params for press bcs'
C              GO TO 9999
C            ENDIF
C end new
C          ENDIF !KTYP51(nr).EQ.3/.GT.3
        ENDIF !(KTYP57(nr).EQ.2).AND.(KTYP51(nr).EQ.3)
      ELSE
        KTYP57(nr)=0
      ENDIF

      CALL IPINIT_ELAS(NAN,NBH,NBHF,NEELEM,NEL,NENP,NFF,
     '  NHE,NHP,NKEF,NKH,NKHE,NLL,NNF,NNL,NPF,NPL,NPLIST,NPNE,
     '  NPNODE,nr,NVHE,NVHP,NW,nx,NYNE,NYNP,DF,DL,PG,SE,WG,XP,
     '  YP,ZA,ZP,FIX,NOFIX,ERROR,*9999)
      
      IF(IOTYPE.NE.3 !new conditions
     '  .AND.(KTYP5.NE.3 !not restart from previous solution
     '  .OR.UPDATE_REF)) THEN
        ! correct and/or save the current solution.
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,
     '    NPNODE,nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        HZERO=0.0d0
        IF(KTYP5.EQ.1.AND. ! zero initial conditions
     '    KTYP52(nr).GE.2.AND.KTYP52(nr).NE.6.AND.
     '    KTYP51(nr).NE.4.AND.KTYP51(nr).NE.5.AND.
     '    .NOT.(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3)) THEN
C old - JHC 23-nov-2004 this block is moved to inside the elem and node loop          
C          IF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.1)) THEN !if using grid coupling and not grid at gauss scheme
CC news JHC 22-NOV-2004: Calling UPGRID to update grid points with green
CC                     strain components
C
CC NOTE: The arrays to store the green strain have been hardcoded to be RCQS
CC       The indices of RCQS are also hardcoded. If these numbers and array change
CC       a way of parsing perl variables into CMISS code must be coded
C
C            STRING='fem update grid green_strain no_ze_calc 
C     '              comp 1,2,3,4,5,6 RCQS 3,6,8,4,5,6 
C     '              ALL_VARIANTS ELEM '
C            CO(1)='FEM'
C            CO(2)='UPDATE'
C            CO(3)='GRID'
C            CO(4)='GREEN_STRAIN'
C            CO(5)='no_ze_calc'
C            CO(6)='COMPONENT'
C            CO(7)='1,2,3,4,5,6'
C            CO(8)='RCQS'
C            CO(9)='3,6,8,4,5,7'
C            CO(10)='ALL_VARIANTS'
C            CO(11)='ELEMENTS'
C            WRITE(CO(12),'(I5)') ne
C            CO(13)='region'
C            WRITE(CO(14),'(I5)') nr
C            NTCO=14
C            CALL UPGRID(IBT,IDO,INP,%VAL(ICQS_SPATIAL_PTR),
C     '        %VAL(IRCQS_SPATIAL_PTR),NAN,%VAL(NAQ_PTR),
C     '        %VAL(NBH_PTR),%VAL(NBJ_PTR),%VAL(NEELEM_PTR),
C     '        %VAL(NELIST_PTR),%VAL(NENQ_PTR),%VAL(NGAP_PTR),
C     '        %VAL(NHE_PTR),%VAL(NHP_PTR),%VAL(NKH_PTR),
C     '        %VAL(NKHE_PTR),%VAL(NKJE_PTR),%VAL(NLL_PTR),%VAL(NLQ_PTR),
C     '        %VAL(NPF_PTR),%VAL(NPL_PTR),NPNE,%VAL(NPNODE_PTR),
C     '        %VAL(NQGP_PTR),%VAL(NQLIST_PTR),%VAL(NQNE_PTR),
C     '        %VAL(NQXI_PTR),%VAL(NQS_PTR),%VAL(NRLIST_PTR),
C     '        %VAL(NVHE_PTR),%VAL(NVHP_PTR),%VAL(NVJE_PTR),%VAL(NW_PTR),
C     '        %VAL(NWQ_PTR),%VAL(NXLIST_PTR),%VAL(NXQ_PTR),
C     '        %VAL(NYNE_PTR),%VAL(NYNP_PTR),%VAL(AQ_PTR),%VAL(CE_PTR),
C     '        %VAL(CP_PTR),%VAL(CQ_PTR),%VAL(CURVCORRECT_PTR),
C     '        %VAL(DL_PTR),%VAL(DNUDXQ_PTR),%VAL(DXDXIQ_PTR),
C     '        %VAL(DXDXIQ2_PTR),%VAL(FEXT_PTR),%VAL(GCHQ_PTR),
C     '        %VAL(GUQ_PTR),PG,%VAL(PROPQ_PTR),
C     '        %VAL(RCQS_SPATIAL_PTR),SE,%VAL(XA_PTR),
C     '        XE,XG,%VAL(XIQ_PTR),
C     '        %VAL(XP_PTR),%VAL(XQ_PTR),%VAL(YG_PTR),%VAL(YP_PTR),
C     '        %VAL(YQ_PTR),%VAL(YQS_PTR),%VAL(ZA_PTR),
C     '        ZE,ZG,%VAL(ZP_PTR),
C     '        STRING,ERROR,*9999)
CC Calling solve to evaluate cell through FEM.f
CC NOTE: Class is hardcoded to be class 2. If class changes
CC       code must be written to identify class
C            STRING='FEM solve class 2 restart to 0'
C            CO(1)='FEM'
C            CO(2)='SOLVE'
C            CO(3)='CLASS'
C            CO(4)='2'
C            CO(5)='RESTART'
C            CO(6)='TO'
C            CO(7)='0'
C            NTCO=7
C            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)
CC       Calling UPGAUS through FEM.f to update gauss array with
CC       stresses calculated via grid coupling
CC NOTE:       Stresses at grid points currently YQS array 2,3,4,5,6,7
C            STRING='FEM update gauss gridvars yqs 2 yg 1'
C            CO(1)='FEM'
C            CO(2)='UPDATE'
C            CO(3)='GAUSS'
C            CO(4)='GRIDVARS'
C            CO(5)='YQS'
C            CO(6)='2,3,4,5,6,7'
C            CO(7)='YG'
C            CO(8)='1,2,3,4,5,6'
C            CO(9)='INCLUDE'
C            CO(10)='ELEMENT'
C            WRITE(CO(11),'(I5)') ne
C            CO(12)='region'
C            WRITE(CO(13),'(I5)') nr
C            NTCO=13
C            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)        
C            HZERO=-(YG(1,1)+YG(4,1)+YG(6,1))/3.0d0
CC Adding code to provide more efficient gauss point
CC stress with grid coupling functionality.
C          ELSEIF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.2)) THEN
CC Added ZGTG53ATGRID call to update current grid/gauss point with green strain values
C
C            CALL ZGTG53GRIDFROMGAUS(AZL,%VAL(ICQS_SPATIAL_PTR),
C     '        %VAL(IRCQS_SPATIAL_PTR),ne,ng,%VAL(NQLIST_PTR),
C     '        %VAL(NQNE_PTR),%VAL(RCQS_SPATIAL_PTR),ERROR,*9999)
CC evaluate cellml file for 1 grid/gauss point
C            CALL ZGTG53EVALCELL(%VAL(CELL_ICQS_VALUE_PTR),
C     &        %VAL(CELL_RCQS_VALUE_PTR),
C     &        %VAL(ICQS_SPATIAL_PTR),
C     &        %VAL(IICQS_SPATIAL_PTR),
C     &        %VAL(IRCQS_SPATIAL_PTR),
C     &        ne,ng,%VAL(NQNE_PTR),
C     &        %VAL(RCQS_SPATIAL_PTR),
C     &        CELLML_DUMMY_ROUTINE,
C     &        %VAL(YQS_PTR),ERROR,*9999)
CC NOTE:   Stresses at grid points currently YQS array 2,3,4,5,6,7
CC         update stresses from YQS to YG for 1 gauss/grid point
C            CALL ZGTG53GAUSFROMGRID(ne,ng,%VAL(NGLIST_PTR),
C     &        %VAL(NQLIST_PTR),%VAL(NQNE_PTR),%VAL(YG_PTR),
C     &        %VAL(YQS_PTR),ERROR,*9999)
C
C            HZERO=-(YG(1,1)+YG(4,1)+YG(6,1))/3.0d0
CC newe JHC
CC         incomp (+fluid), not membrane or string and not const vol constraint
C
C          ELSEIF(KTYP55(nr).EQ.1) THEN      !princ strain invariants
CC news VJ 10Dec2003 Added YG to ENERGY param list
C            nb=NBH(NH_LOC(1,nx),1,NEELEM(1,nr))
C            CALL CPCG(1,nb,NPNE(1,1,NEELEM(1,nr)),nr,nx,
C     '      CE(1,NEELEM(1,nr)),CG,CGE(1,1,NEELEM(1,nr)),
C     '      CP,PG,ERROR,*9999)
C            CALL ENERGY(nr,CG(1,1),DW,3.0d0,3.0d0,1.0d0,0.0d0,0.0d0,
C     '        0.0d0,YG(1,1,NEELEM(1,nr)),ERROR,*9999)
C            HZERO=-(DW(1)+2.0d0*DW(2))
CC TVK 12/01/2000 Set logical var to calculate correct equilm press
C            IF(KTYP52(nr).EQ.4) THEN
C              EQUIM_PRESSURE=.TRUE.
C            ENDIF
CC TVK 02/12/1999 Correct pressure at equilibrium
CC            IF(KTYP52(nr).EQ.4) THEN
CC              IF(ITYP10(nr).EQ.1) THEN
CC                NITB=NIT(nb)
CC                DO i=1,3
CC                  DO j=1,3
CC                    AXU(i,j)=0.0d0
CC                  ENDDO
CC                  AXU(i,i)=1.0d0
CC                ENDDO
CC                CALL ZGMG(NBH(NH_LOC(1,nx),1,NEELEM(1,nr)),nr,AZ,AZL,
CC     '            AZU,ZG,ERROR,*9999)
CC                IF(KTYP53(nr).GT.1) THEN
CC                  EG13=AZL(1,3)/2.0d0
CC                  EG12=AZL(1,2)/2.0d0
CC                  RK1=(AZL(1,1)-1.0d0)/2.0d0
CC                  RK2=EG13*EG13+EG12*EG12
CC                  BG(1,1)=3.0d0-AZL(1,1)
CC                  BG(2,2)=3.0d0-AZL(2,2)
CC                  BG(3,3)=3.0d0-AZL(3,3)
CC                  BG(2,1)=   -AZL(2,1)
CC                  BG(3,1)=   -AZL(3,1)
CC                  BG(3,2)=   -AZL(3,2)
CC                  CALL ENERGY(nr,CG(1,1),DW,3.0d0,3.0d0,1.0d0,0.0d0,
CC     '              0.0d0,0.0d0,YG(1,1,NEELEM(1,nr)),ERROR,*9999)
CC                  IF(KTYP52(nr).EQ.4) THEN !compressible + fluid
CC TVK no compliance on I3 term  W3=RI3*DW(3)*ZG(NH_LOC(0,nx),1)
CC                    W3=1.0d0*DW(3)
CC                  ENDIF
CC                  DO mj=1,NITB
CC                    DO nj=1,mj
CC                      TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj)
CC     '                  + W3*AZU(mj,nj))
CC                    ENDDO
CC                  ENDDO
CC         Note: DW(4) & DW(5) are zero for isotropic case
CC                  TG(3,1)=TG(3,1)+EG13*DW(5)
CC                  TG(2,1)=TG(2,1)+EG12*DW(5)
CC                  TG(1,1)=TG(1,1)+AXU(1,1)*DW(4)
CC                  TG(1,2)=TG(2,1)
CC                  TG(1,3)=TG(3,1)
CC                  TG(2,3)=TG(3,2)
CC                ENDIF
CC              ENDIF
CC              HZERO=TG(1,1)/AXU(1,1)
CC            ENDIF
C
C          ELSE IF(KTYP55(nr).EQ.2) THEN !princ extensions
CC news VJ 10Dec2003 Added YG to ENERGY param list
C            CALL ENERGY(nr,CG(1,1),DW,1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,
C     '        0.0d0,YG(1,1,NEELEM(1,nr)),ERROR,*9999)
CC!!!        30Jun88: Check this expression for W a function of L1..L3 ??
C            HZERO=-(DW(1)+2.0d0*DW(2))
C          ELSE IF(KTYP55(nr).EQ.3) THEN !fibre strains
CC news VJ 10Dec2003 Added YG to ENERGY param list
C            CALL ENERGY(nr,CG(1,1),DW,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
C     '        0.0d0,YG(1,1,NEELEM(1,nr)),ERROR,*9999)
CC!!!        30Jun88: Put in a consistency check here
C!            HZERO=-DW(1)/2.0d0
C            HZERO=0.d0
C          ENDIF
C          IF(DOP) THEN
CC KAT 14May01: Can't branch out of critical section.
CC              Critical section is not essential.
CCC$          call mp_setlock()
C            WRITE(OP_STRING,'('' HZERO ='',D12.5)') HZERO
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            WRITE(OP_STRING,
C     '        '('' DW(i) ='',6(D10.2,2X))') (DW(i),i=1,6)
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CCC$          call mp_unsetlock()
C          ENDIF
C end old
          IF(KTYP52(nr).EQ.5) THEN !incomp+inext
            nh_pressure=NH_LOC(NH_LOC(0,nx)-1,nx)
          ELSE !not incomp+inext
            nh_pressure=NH_LOC(NH_LOC(0,nx),nx)
          ENDIF
          IF(NNT(NBH(nh_pressure,1,NEELEM(1,nr))).NE.0) THEN
C NEWS JHC 25-NOV-2004 Check if grid coupling is used, and if so
C warn user to re-define init condition for appropriate calculation 
C of initial hydrostatic pressure
            IF(KTYP54(nr).EQ.3)THEN
              WRITE(OP_STRING,'('' >>Warning: You will need to '
     '          //'define initial condition again after defining '
     '          //'material laws for CellML.'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
C NEWE  
            DO nonode=1,NPNODE(0,nr)
              ne=NEELEM(1,nr)
              IF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.1)) THEN !if using grid coupling and not grid at gauss scheme
C news JHC 22-NOV-2004: Calling UPGRID to update grid points with green
C                     strain components

C NOTE: The arrays to store the green strain have been hardcoded to be RCQS
C       The indices of RCQS are also hardcoded. If these numbers and array change
C       a way of parsing perl variables into CMISS code must be coded

                STRING='fem update grid green_strain no_ze_calc 
     '                  comp 1,2,3,4,5,6 RCQS 3,6,8,4,5,6 
     '                  ALL_VARIANTS ELEM '
                CALL UPGRID(IBT,IDO,INP,%VAL(ICQS_SPATIAL_PTR),
     '            %VAL(IRCQS_SPATIAL_PTR),NAN,%VAL(NAQ_PTR),
     '            %VAL(NBH_PTR),%VAL(NBJ_PTR),%VAL(NEELEM_PTR),
     '            %VAL(NELIST_PTR),%VAL(NENQ_PTR),%VAL(NGAP_PTR),
     '            %VAL(NHE_PTR),%VAL(NHP_PTR),%VAL(NKH_PTR),
     '            %VAL(NKHE_PTR),%VAL(NKJE_PTR),%VAL(NLL_PTR),
     '            %VAL(NLQ_PTR),%VAL(NPF_PTR),%VAL(NPL_PTR),NPNE,
     '            %VAL(NPNODE_PTR),%VAL(NQGP_PTR),%VAL(NQLIST_PTR),
     '            %VAL(NQNE_PTR),%VAL(NQXI_PTR),%VAL(NQS_PTR),
     '            %VAL(NRLIST_PTR),%VAL(NVHE_PTR),%VAL(NVHP_PTR),
     '            %VAL(NVJE_PTR),%VAL(NW_PTR),%VAL(NWQ_PTR),
     '            %VAL(NXLIST_PTR),%VAL(NXQ_PTR),%VAL(NYNE_PTR),
     '            %VAL(NYNP_PTR),%VAL(AQ_PTR),%VAL(CE_PTR),%VAL(CP_PTR),
     '            %VAL(CQ_PTR),%VAL(CURVCORRECT_PTR),%VAL(DL_PTR),
     '            %VAL(DNUDXQ_PTR),%VAL(DXDXIQ_PTR),%VAL(DXDXIQ2_PTR),
     '            %VAL(FEXT_PTR),%VAL(GCHQ_PTR),%VAL(GUQ_PTR),PG,
     '            %VAL(PROPQ_PTR),%VAL(RCQS_SPATIAL_PTR),SE,
     '            %VAL(XA_PTR),XE,XG,%VAL(XIQ_PTR),
     '            %VAL(XP_PTR),%VAL(XQ_PTR),%VAL(YG_PTR),%VAL(YP_PTR),
     '            %VAL(YQ_PTR),%VAL(YQS_PTR),%VAL(ZA_PTR),
     '            ZE,ZG,%VAL(ZP_PTR),
     '            STRING,ERROR,*9999)
C Calling solve to evaluate cell through FEM.f
C NOTE: Class is hardcoded to be class 2. If class changes
C       code must be written to identify class
                STRING='FEM solve class 2 restart to 0'
                CO(1)='FEM'
                CO(2)='SOLVE'
                CO(3)='CLASS'
                CO(4)='2'
                CO(5)='RESTART'
                CO(6)='TO'
                CO(7)='0'
                NTCO=7
                CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,
     '            *9999)
C       Calling UPGAUS through FEM.f to update gauss array with
C       stresses calculated via grid coupling
C NOTE:       Stresses at grid points currently YQS array 2,3,4,5,6,7
                STRING='FEM update gauss gridvars yqs 2 yg 1'
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
                CO(12)='region'
                WRITE(CO(13),'(I5)') nr
                NTCO=13
                CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,
     '            *9999)        
                HZERO=-(YG(1,1,ne)+YG(4,1,ne)+YG(6,1,ne))/6.0d0
C Adding code to provide more efficient gauss point
C stress with grid coupling functionality.
              ELSEIF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.2)) THEN
C Added ZGTG53ATGRID call to update current grid/gauss point with green strain values
                DO temp_i=1,3
                  DO temp_j=1,3
                    IF(temp_i.EQ.temp_j)THEN
                      AZL(temp_i,temp_j)=1.0d0
                    ELSE
                      AZL(temp_i,temp_j)=0.0d0
                    ENDIF
                  ENDDO
                ENDDO
                CALL ZGTG53GRIDFROMGAUS(AZL,%VAL(ICQS_SPATIAL_PTR),
     '            %VAL(IRCQS_SPATIAL_PTR),ne,1,%VAL(NQLIST_PTR),
     '            %VAL(NQNE_PTR),%VAL(RCQS_SPATIAL_PTR),ERROR,*9999)
C evaluate cellml file for 1 grid/gauss point
                CALL ZGTG53EVALCELL(%VAL(CELL_ICQS_VALUE_PTR),
     &            %VAL(CELL_RCQS_VALUE_PTR),
     &            %VAL(ICQS_SPATIAL_PTR),
     &            %VAL(IICQS_SPATIAL_PTR),
     &            %VAL(IRCQS_SPATIAL_PTR),
     &            ne,1,%VAL(NQNE_PTR),
     &            %VAL(RCQS_SPATIAL_PTR),
     &            %VAL(YQS_PTR),ERROR,*9999)
C NOTE:   Stresses at grid points currently YQS array 2,3,4,5,6,7
C         update stresses from YQS to YG for 1 gauss/grid point
                CALL ZGTG53GAUSFROMGRID(ne,1,%VAL(NGLIST_PTR),
     &            %VAL(NQLIST_PTR),%VAL(NQNE_PTR),%VAL(YG_PTR),
     &            %VAL(YQS_PTR),ERROR,*9999)

                HZERO=-(YG(1,1,ne)+YG(4,1,ne)+YG(6,1,ne))/6.0d0
C newe JHC
C         incomp (+fluid), not membrane or string and not const vol constraint
              ELSE
C
C news JHC 10-03-2005 call CPCG for principle strain invariant, 
C principle extensions and fibre strains
                nb=NBH(NH_LOC(1,nx),1,ne)
                CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,
     '            CE(1,ne),CG,CGE(1,1,ne),
     '            CP,PG,ERROR,*9999)
C newe
                IF(KTYP55(nr).EQ.1) THEN      !princ strain invariants
C news VJ 10Dec2003 Added YG to ENERGY param list

C news JHC 10-03-2005 move this block outside if-statement so that CPCG can be 
C used for principle strain invariant, principle extensions and fibre strains
C                nb=NBH(NH_LOC(1,nx),1,ne)
C                CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,
C     '            CE(1,ne),CG,CGE(1,1,ne),
C     '            CP,PG,ERROR,*9999)
C newe
                  CALL ENERGY(nr,CG(1,1),DW,3.0d0,3.0d0,1.0d0,0.0d0,
     '              0.0d0,
     '              0.0d0,YG(1,1,ne),ERROR,*9999)
                  HZERO=-(DW(1)+2.0d0*DW(2))
C TVK 12/01/2000 Set logical var to calculate correct equilm press
                  IF(KTYP52(nr).EQ.4) THEN
                    EQUIM_PRESSURE=.TRUE.
                  ENDIF
                ELSE IF(KTYP55(nr).EQ.2) THEN !princ extensions
C news VJ 10Dec2003 Added YG to ENERGY param list
                  CALL ENERGY(nr,CG(1,1),DW,1.0d0,1.0d0,1.0d0,0.0d0,
     '              0.0d0,0.0d0,YG(1,1,ne),ERROR,*9999)
C!!!        30Jun88: Check this expression for W a function of L1..L3 ??
                  HZERO=-(DW(1)+2.0d0*DW(2))
                ELSE IF(KTYP55(nr).EQ.3) THEN !fibre strains
C news VJ 10Dec2003 Added YG to ENERGY param list
                  CALL ENERGY(nr,CG(1,1),DW,0.0d0,0.0d0,0.0d0,0.0d0,
     '              0.0d0,0.0d0,YG(1,1,ne),ERROR,*9999)
C!!!        30Jun88: Put in a consistency check here
!            HZERO=-DW(1)/2.0d0
                  HZERO=0.d0
                ENDIF
              ENDIF

              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
                WRITE(OP_STRING,'('' HZERO ='',D12.5)') HZERO
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' DW(i) ='',6(D10.2,2X))') (DW(i),i=1,6)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$               call mp_unsetlock()
              ENDIF
              np=NPNODE(nonode,nr)
              DO nv=1,NVHP(nh_pressure,np,1)
                ZP(1,nv,nh_pressure,np,1)=HZERO
                DO nk=2,NKH(nh_pressure,np,1)
                  ZP(nk,nv,nh_pressure,np,1)=0.0d0
                ENDDO !nk
              ENDDO !nv
            ENDDO !nonode
          ELSE
C NEWS JHC 25-NOV-2004 Check if grid coupling is used, and if so
C warn user to re-define init condition for appropriate calculation 
C of initial hydrostatic pressure
            IF(KTYP54(nr).EQ.3)THEN
              WRITE(OP_STRING,'('' >>Warning: You will need to '
     '          //'define initial condition again after defining '
     '          //'material laws for CellML.'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
C NEWE
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.1)) THEN !if using grid coupling and not grid at gauss scheme
C news JHC 22-NOV-2004: Calling UPGRID to update grid points with green
C                     strain components

C NOTE: The arrays to store the green strain have been hardcoded to be RCQS
C       The indices of RCQS are also hardcoded. If these numbers and array change
C       a way of parsing perl variables into CMISS code must be coded

                STRING='fem update grid green_strain no_ze_calc 
     '                  comp 1,2,3,4,5,6 RCQS 3,6,8,4,5,6 
     '                  ALL_VARIANTS ELEM '
                CALL UPGRID(IBT,IDO,INP,%VAL(ICQS_SPATIAL_PTR),
     '            %VAL(IRCQS_SPATIAL_PTR),NAN,%VAL(NAQ_PTR),
     '            %VAL(NBH_PTR),%VAL(NBJ_PTR),%VAL(NEELEM_PTR),
     '            %VAL(NELIST_PTR),%VAL(NENQ_PTR),%VAL(NGAP_PTR),
     '            %VAL(NHE_PTR),%VAL(NHP_PTR),%VAL(NKH_PTR),
     '            %VAL(NKHE_PTR),%VAL(NKJE_PTR),%VAL(NLL_PTR),
     '            %VAL(NLQ_PTR),%VAL(NPF_PTR),%VAL(NPL_PTR),NPNE,
     '            %VAL(NPNODE_PTR),%VAL(NQGP_PTR),%VAL(NQLIST_PTR),
     '            %VAL(NQNE_PTR),%VAL(NQXI_PTR),%VAL(NQS_PTR),
     '            %VAL(NRLIST_PTR),%VAL(NVHE_PTR),%VAL(NVHP_PTR),
     '            %VAL(NVJE_PTR),%VAL(NW_PTR),%VAL(NWQ_PTR),
     '            %VAL(NXLIST_PTR),%VAL(NXQ_PTR),%VAL(NYNE_PTR),
     '            %VAL(NYNP_PTR),%VAL(AQ_PTR),%VAL(CE_PTR),%VAL(CP_PTR),
     '            %VAL(CQ_PTR),%VAL(CURVCORRECT_PTR),%VAL(DL_PTR),
     '            %VAL(DNUDXQ_PTR),%VAL(DXDXIQ_PTR),%VAL(DXDXIQ2_PTR),
     '            %VAL(FEXT_PTR),%VAL(GCHQ_PTR),%VAL(GUQ_PTR),PG,
     '            %VAL(PROPQ_PTR),%VAL(RCQS_SPATIAL_PTR),SE,
     '            %VAL(XA_PTR),XE,XG,%VAL(XIQ_PTR),
     '            %VAL(XP_PTR),%VAL(XQ_PTR),%VAL(YG_PTR),%VAL(YP_PTR),
     '            %VAL(YQ_PTR),%VAL(YQS_PTR),%VAL(ZA_PTR),
     '            ZE,ZG,%VAL(ZP_PTR),
     '            STRING,ERROR,*9999)
C Calling solve to evaluate cell through FEM.f
C NOTE: Class is hardcoded to be class 2. If class changes
C       code must be written to identify class
                STRING='FEM solve class 2 restart to 0'
                CO(1)='FEM'
                CO(2)='SOLVE'
                CO(3)='CLASS'
                CO(4)='2'
                CO(5)='RESTART'
                CO(6)='TO'
                CO(7)='0'
                NTCO=7
                CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,
     '            *9999)
C       Calling UPGAUS through FEM.f to update gauss array with
C       stresses calculated via grid coupling
C NOTE:       Stresses at grid points currently YQS array 2,3,4,5,6,7
                STRING='FEM update gauss gridvars yqs 2 yg 1'
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
                CO(12)='region'
                WRITE(CO(13),'(I5)') nr
                NTCO=13
                CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,
     '            *9999)        
                HZERO=-(YG(1,1,ne)+YG(4,1,ne)+YG(6,1,ne))/6.0d0
C Adding code to provide more efficient gauss point
C stress with grid coupling functionality.
              ELSEIF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.2)) THEN
C Added ZGTG53ATGRID call to update current grid/gauss point with green strain values
                DO temp_i=1,3
                  DO temp_j=1,3
                    IF(temp_i.EQ.temp_j)THEN
                      AZL(temp_i,temp_j)=1.0d0
                    ELSE
                      AZL(temp_i,temp_j)=0.0d0
                    ENDIF
                  ENDDO
                ENDDO
                CALL ZGTG53GRIDFROMGAUS(AZL,%VAL(ICQS_SPATIAL_PTR),
     '            %VAL(IRCQS_SPATIAL_PTR),ne,1,%VAL(NQLIST_PTR),
     '            %VAL(NQNE_PTR),%VAL(RCQS_SPATIAL_PTR),ERROR,*9999)
C evaluate cellml file for 1 grid/gauss point
                CALL ZGTG53EVALCELL(%VAL(CELL_ICQS_VALUE_PTR),
     &            %VAL(CELL_RCQS_VALUE_PTR),
     &            %VAL(ICQS_SPATIAL_PTR),
     &            %VAL(IICQS_SPATIAL_PTR),
     &            %VAL(IRCQS_SPATIAL_PTR),
     &            ne,1,%VAL(NQNE_PTR),
     &            %VAL(RCQS_SPATIAL_PTR),
     &            %VAL(YQS_PTR),ERROR,*9999)
C NOTE:   Stresses at grid points currently YQS array 2,3,4,5,6,7
C         update stresses from YQS to YG for 1 gauss/grid point
                CALL ZGTG53GAUSFROMGRID(ne,1,%VAL(NGLIST_PTR),
     &            %VAL(NQLIST_PTR),%VAL(NQNE_PTR),%VAL(YG_PTR),
     &            %VAL(YQS_PTR),ERROR,*9999)

                HZERO=-(YG(1,1,ne)+YG(4,1,ne)+YG(6,1,ne))/6.0d0
C newe JHC
C         incomp (+fluid), not membrane or string and not const vol constraint
              ELSE
C news JHC 10-03-2005 call CPCG for principle strain invariant, 
C principle extensions and fibre strains
                nb=NBH(NH_LOC(1,nx),1,ne)
                CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,
     '            CE(1,ne),CG,CGE(1,1,ne),
     '            CP,PG,ERROR,*9999)
C newe
                IF(KTYP55(nr).EQ.1) THEN      !princ strain invariants
C news VJ 10Dec2003 Added YG to ENERGY param list

C news JHC 10-03-2005 move this block outside if-statement so that CPCG can be 
C used for principle strain invariant, principle extensions and fibre strains
C                nb=NBH(NH_LOC(1,nx),1,ne)
C                CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,
C     '            CE(1,ne),CG,CGE(1,1,ne),
C     '            CP,PG,ERROR,*9999)
C newe
                  CALL ENERGY(nr,CG(1,1),DW,3.0d0,3.0d0,1.0d0,0.0d0,
     '              0.0d0,
     '              0.0d0,YG(1,1,NEELEM(1,nr)),ERROR,*9999)
                  HZERO=-(DW(1)+2.0d0*DW(2))
C TVK 12/01/2000 Set logical var to calculate correct equilm press
                  IF(KTYP52(nr).EQ.4) THEN
                    EQUIM_PRESSURE=.TRUE.
                  ENDIF
                ELSE IF(KTYP55(nr).EQ.2) THEN !princ extensions
C news VJ 10Dec2003 Added YG to ENERGY param list
                  CALL ENERGY(nr,CG(1,1),DW,1.0d0,1.0d0,1.0d0,0.0d0,
     '              0.0d0,0.0d0,YG(1,1,ne),ERROR,*9999)
C!!!        30Jun88: Check this expression for W a function of L1..L3 ??
                  HZERO=-(DW(1)+2.0d0*DW(2))
                ELSE IF(KTYP55(nr).EQ.3) THEN !fibre strains
C news VJ 10Dec2003 Added YG to ENERGY param list
                  CALL ENERGY(nr,CG(1,1),DW,0.0d0,0.0d0,0.0d0,0.0d0,
     '              0.0d0,0.0d0,YG(1,1,ne),ERROR,*9999)
C!!!        30Jun88: Put in a consistency check here
!            HZERO=-DW(1)/2.0d0
                  HZERO=0.d0
                ENDIF
              ENDIF
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
                WRITE(OP_STRING,'('' HZERO ='',D12.5)') HZERO
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' DW(i) ='',6(D10.2,2X))') (DW(i),i=1,6)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$               call mp_unsetlock()
              ENDIF
              
              DO na=1,NAT(NBH(nh_pressure,1,ne))
                ZA(na,nh_pressure,1,ne)=0.0d0
              ENDDO !na
              ZA(1,nh_pressure,1,ne)=HZERO
            ENDDO !noelem
          ENDIF

          ! Update current solution
          CALL ZPYP(1,NBH,NEELEM,NHE,NHP,NKH,
     '      NPNODE,nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        ENDIF

        IF(KTYP5.NE.3) THEN !not restart from previous solution
          ! Save initial solution
          CALL ZPYP(3,NBH,NEELEM,NHE,NHP,NKH,
     '      NPNODE,nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        ENDIF

        IF(UPDATE_REF) THEN !update reference state stored in YP(ny,10)
C         Copy current solution into YP(ny,10) to act as reference
C         state (used for cavity problems in ZERE55)
          CALL ZPYP(10,NBH,NEELEM,NHE,NHP,NKH,
     '      NPNODE,nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        ENDIF
      ENDIF ! correct and/or save the current solution.
        
      IF(IOTYPE.NE.3 !new conditions
     '  .AND..NOT.UPDATE_REF.AND.
     '  ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
        WRITE(OP_STRING,
     '    '('' >>WARNING: Set up reference state for '//
     '    'cavity region using "up solu cav"'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C LC 25/2/97 archived section :
C       MPN 15-Sep-95: not used anymore
C       MPN 16-Aug-94: KTYP57(nr)=5 now used for entering cavity pressures

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


