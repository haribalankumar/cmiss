      SUBROUTINE ZERE50(IBT,IDO,INP,NAN,NBH,NBJ,NBJF,ne,NFF,NGAP,
     '  NHE,NKEF,NNF,NPNE,nr,NRE,NW,nx,NXI,
     '  CE,CG,CP,FEXT,PG,RE,RGX,SE,WG,XE,XG,YG,ZE,ZEA,ZG,ERROR,*)

C#### Subroutine: ZERE50
C###  Description:
C###    ZERE50 calculates element residual RE from current dependent
C###    variable array ZE.

C**** X variables refer to orthog curvilinear coords in reference state.
C**** Z     "       "         "        "         "      deformed    "
C**** Material Theta-coordinates (reference for deformation) coincide
C****   with Xj-coords in ref state.
C**** Material Xi-coords are the finite element mesh coordinates.
C**** Material Nu-coordinates (reference for stresses): are orthogonal
C****   and (Nu1,Nu2) lie in the (Xi1-Xi2) plane such that Nu(1)
C****   is aligned with the 'fibres' to which material aeolotropy
C****   is referred; The undeformed base vectors are defined such that
C****   the undeformed metric tensors wrt the Nu are delta(i,j).
C****
C**** ITYP2(nr,nx) is 1..15 for problem type
C**** ITYP4(nr,nx) is 1..4: fem/direct bem/indirect bem/orthog colloc.
C**** ITYP5(nr,nx) is 1..5: static/time integration/modal analysis
C****                       /Fourier analysis/buckling analysis
C**** ITYP6(nr,nx) is 1,2:  linear/nonlinear problem
C**** KTYP5  is 1..3: initial solution zero/read in/restarted
C**** KTYP7  is 1..3: eqn parameters constant wrt time/user defined
C****                 /read from file
C**** KTYP8  is 1..6: geom/fibre/field/potential/Fourier/opt.n fitting
C**** ITYP9(nr,nx) is 1..3: solution by full Newton/modified Newton
C****                 /BFGS inverse/Conjugate gradient
C**** KTYP10 is 1,2 : solution with no search/linear search
C**** KTYP12 is 1,0 : fitting with/without constraints
C**** KTYP13 is 1 if pressure read from file (PRESS.VSAERO)
C**** KTYP16 is 1,2 : lowest/highest eigenvalue
C**** KTYP17 is number of eigenvalue pairs
C**** KTYP19 is number of starting vectors
C**** KTYP22 is 1..3: time integ.n algorithm linear/quadratic/cubic
C**** KTYP23 is 1,2:  time step fixed/calculated to control error
C**** KTYP25 is 1..3: b.c. in form of impulse/step/sine wave
C**** KTYP26 is 1..2: opt of material params/geometric params
C**** KTYP27 is 1..5: type of minimization objective function
C**** KTYP31 is 1,2:  activation model forwards/backwards
C**** KTYP43 is 1..3: linear elastic material isotropic/trans.isotr.
C****                 /orthotropic
C**** KTYP51(nr) is 1..6: plane stress/plane strain/3D/membrane
C****                 /string/shell
C**** KTYP52(nr) is 1..5: compress/incomp/incomp+fluid/comp+fluid/incomp+inext/comp+fluid for lung
C**** KTYP53(nr) is 1..3: isotropic/aeleotropic/aeleo + active fibres
C**** KTYP54(nr) is 1..3: hyperelastic/Cauchy-elasticity/creep
C**** KTYP55(nr) is 1..3: strain invariants/ext ratios/fibre strains
C**** KTYP56(nr) is 1..3: polynomial/special function/user defined
C**** KTYP57(nr) is type of pressure loading applied to elements
C**** KTYP58(nr) is 1,2: conventional/isochoric element
C**** KTYP59(nr)  is elastance/Hill-type/fading-memory formulation
C**** KYTP59S(nr) is type of active stress added to the 2nd PKST
C**** KTYP59S1(nr) is the row component the active stress gets added to
C**** KTYP59S2(nr) is the column component the active stress gets added to
C**** ASCOMP(nr) is the name of the active stress component within CELLML
C**** KTYP5I(nr) is 0..1: No inertia/inertia

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'acti00.cmn'
      INCLUDE 'aero00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NBJF(NJM,NFM),ne,NFF(6),
     '  NGAP(NIM,NBM),NHE,NKEF(0:4,16,6,NBFM),NNF(0:17,6,NBFM),
     '  NPNE(NNM,NBFM,NEM),nr,NRE(NEM),NW,nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),FEXT(NIFEXTM,NGM),
     '  PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),RGX(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM),YG(NIYGM,NGM),
     &  ZE(NSM,NHM),ZEA(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER adj_dirn,iface,IFE,IXF,i,INTWORK(1),ISEG(1),j,JP,k,
     '  mi,mj,na,nb,nb_pressure,NBE,NBFF,NBP,
     '  neadj,nf,ng,nh,nh_pressure,nh1,nh2,nh3,nhx,ni,NITB,nj,njj2,
     '  nk,nkbf,nn,nnbf,ns,nse,nsf,NSP(-2:2),NU1(0:3)
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 ACC(NHM),AG,AGE,AGE1,AGE2,AGE3,AXU(3,3),AZ,AZL(3,3),
     '  AZU(3,3),BG(3,3),CHTOFF(3,3,3),CLZ,CMZ,Comp_resid,CSLZ,CSMZ,
     '  D11,D12,D13,D21,D22,D23,D31,D32,D33,DBM(3,3,3),delsqP,
     '  DET_F_NU,DLA1,DLA2,DLA3,DMA1,DMA2,DMA3,DTA1,DTA2,DTA3,DW(6),
     '  DXIX(3,3),DXIZN(3,3),DXNZN(3,3),dZ1_dNu1,dZ1_dNu2,dZ2_dNu1,
     '  dZ2_dNu2,DZDX(3,3),DZNXI(3,3),DZNXN(3,3),E1,E2,EG(3,3),EG12,
     '  EG13,factor,G1,GXL(3,3),GXU(3,3),Inext_resid,Incomp_resid,
     '  PF(2),PGA1,PGA2,PGA3,PGG,PGX,PPG(64,4),PPGG(4),REALWORK(1),
     '  RF(32,6),RI1,RI2,RI3,RWG,SLZ,
     '  SMZ,SUM,SUM1,TG(3,3),TNA,Volume,W3,XI(3),ZGA(NHM,NUM),
     '  ZG_temp(NHM,NUM),DET,ACTIVE_STRESS
      CHARACTER STRESSTYPE*17,TYPE*9,CSEG(1)*(1),STRING*(MXCH)
      LOGICAL ADJ_XI3_ELEM,ELEMPRESS,END,NODEPRESS,SAMEDEPBASIS
      EXTERNAL CELLML_DUMMY_ROUTINE,USER_CELL1,USER_CELL2,USER_CELL3
      EXTERNAL USER_CELL4,USER_CELL5
      DATA STRESSTYPE/' '/
      DATA NU1/1,2,4,7/

      CALL ENTERS('ZERE50',*9999)
      NITB=NIT(NBJ(1))

C *** Test whether to use unrolled loops
      SAMEDEPBASIS=.TRUE.
      nh1=NH_LOC(1,nx)
      DO nhx=2,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        IF(NBH(nh).NE.NBH(nh1)) SAMEDEPBASIS=.FALSE.
      ENDDO !nhx
      IF(ITYP10(nr).EQ.1.AND.KTYP51(nr).EQ.3.AND.KTYP53(nr).LE.3.
     '  AND.KTYP52(nr).NE.4.AND. !exclude compressible + fluid
     '  NJ_LOC(NJL_GEOM,0,nr).EQ.3.AND.NIT(NBH(NH_LOC(1,nx))).EQ.3) THEN
        TYPE='RC3D'
      ELSE IF(ITYP10(nr).EQ.4.AND.KTYP51(nr).EQ.3.AND.
     '    KTYP53(nr).LE.3.AND.SAMEDEPBASIS.AND..NOT.DOP.AND.
     '    NJ_LOC(NJL_GEOM,0,nr).EQ.3.AND.
     &    NIT(NBH(NH_LOC(1,nx))).EQ.3) THEN
        TYPE='PROLATE'
      ELSE IF(NJT.EQ.2.AND.KTYP51(nr).EQ.4) THEN
        TYPE='BIAXIAL'
      ELSE
        TYPE=' '
      ENDIF

      DO nhx=1,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        DO ns=1,NST(NBH(nh))+NAT(NBH(nh))
          RE(ns,nh)=0.0d0
        ENDDO !ns
      ENDDO !nhx (nh)

      IF(KTYP57(nr).GT.1) THEN !Boundary pressure increments entered
        NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for pressure vars
C       Determine if hydrostatic pressure field is nodally interpolated
        NODEPRESS=.FALSE.
        IF(NBC(NBP).GT.0) THEN
          NODEPRESS=.TRUE.
          DO ni=1,NIT(NBP)
C CS 14/9/200 allow sectors as well
            CALL ASSERT(IBT(1,ni,NBP).EQ.1.OR.IBT(1,ni,NBP).EQ.5.OR.
     '        IBT(1,ni,NBP).EQ.6,
     '        '>> Nodal hyd press interpolation must be trilinear',
     &        ERROR,*9999)
          ENDDO !ni
        ENDIF !NBC(NBP).GT.0
C       Determine if hydrostatic press is interpolated with elem vars
C       and put Xi3 face pressures (aux vars in ZE) into PF array
        ELEMPRESS=.FALSE.
        NSP(-1)=0
        NSP(-2)=0
        NSP(1)=0
        NSP(2)=0
        DO na=1,NAT(NBP)
          IF(NAN(3,na,NBP).EQ.0) THEN
C           Pick up param assoc with const pressure term
            NSP(1)=na
            ELEMPRESS=.TRUE.
          ELSE IF(NAN(3,na,NBP).EQ.1.OR.NAN(3,na,NBP).EQ.3) THEN
C           Pick up param assoc with linear or cubic press term
            NSP(2)=na
            ELEMPRESS=.TRUE.
          ELSE IF(NAN(3,na,NBP).EQ.-1) THEN
C           Pick up param assoc with Xi3=0 face pressure bc
            PF(1)=ZE(NST(NBP)+na,NH_LOC(0,nx))
            NSP(-1)=na
          ELSE IF(NAN(3,na,NBP).EQ.-2) THEN
C           Pick up param assoc with Xi3=1 face pressure bc
            PF(2)=ZE(NST(NBP)+na,NH_LOC(0,nx))
            NSP(-2)=na
          ELSE
            ELEMPRESS=.TRUE.
          ENDIF !NAN
        ENDDO !na
      ENDIF !KTYP57(nr).GT.1

C new MPN 26Mar97
      IF(IWRIT4(nr,nx).GE.1.AND.DOP) THEN
C KAT 14May01: mp_setlock not OPENMP.
C              Critical section is not essential.
CC$      call mp_setlock()
        CALL CPU_TIMER(CPU_USER,TIME_START)
CC$      call mp_unsetlock()
      ENDIF !IWRIT4(nr,nx).GE.1.AND.DOP

C news VJ 31Jan2004: Found UPGRID takes up too much time to interpolate from grid to gauss
C If using grid at gauss scheme, gauss loop further below does necessary calcs also done in UPGRID
C avoid this update of grid properties if using grid at gauss

      IF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.1)) THEN !if using grid coupling and not grid at gauss scheme
C news VJ 17Jan2004: Calling UPGRID to update grid points with green
C                    strain components

C NOTE: The arrays to store the green strain have been hardcoded to be RCQS
C       The indices of RCQS are also hardcoded. If these numbers and array change
C       a way of parsing perl variables into CMISS code must be coded

        STRING='fem update grid green_strain no_ze_calc 
     '          comp 1,2,3,4,5,6 RCQS 3,6,8,4,5,6 ALL_VARIANTS ELEM '
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
        CO(13)='region'
        WRITE(CO(14),'(I5)') nr
        NTCO=14
        CALL UPGRID(IBT,IDO,INP,%VAL(ICQS_SPATIAL_PTR),
     '    %VAL(IRCQS_SPATIAL_PTR),NAN,%VAL(NAQ_PTR),
     '    %VAL(NBH_PTR),%VAL(NBJ_PTR),%VAL(NEELEM_PTR),
     '    %VAL(NELIST_PTR),%VAL(NENQ_PTR),%VAL(NGAP_PTR),
     '    %VAL(NHE_PTR),%VAL(NHP_PTR),%VAL(NKH_PTR),
     '    %VAL(NKHE_PTR),%VAL(NKJE_PTR),%VAL(NLL_PTR),%VAL(NLQ_PTR),
     '    %VAL(NPF_PTR),%VAL(NPL_PTR),NPNE,%VAL(NPNODE_PTR),
     '    %VAL(NQGP_PTR),%VAL(NQLIST_PTR),%VAL(NQNE_PTR),
     '    %VAL(NQXI_PTR),%VAL(NQS_PTR),%VAL(NRLIST_PTR),
     '    %VAL(NVHE_PTR),%VAL(NVHP_PTR),%VAL(NVJE_PTR),%VAL(NW_PTR),
     '    %VAL(NWQ_PTR),%VAL(NXLIST_PTR),%VAL(NXQ_PTR),%VAL(NYNE_PTR),
     '    %VAL(NYNP_PTR),%VAL(AQ_PTR),%VAL(CE_PTR),%VAL(CP_PTR),
     '    %VAL(CQ_PTR),%VAL(CURVCORRECT_PTR),%VAL(DL_PTR),
     '    %VAL(DNUDXQ_PTR),%VAL(DXDXIQ_PTR),%VAL(DXDXIQ2_PTR),
     '    %VAL(FEXT_PTR),%VAL(GCHQ_PTR),
     '    %VAL(GUQ_PTR),PG,%VAL(PROPQ_PTR),
     '    %VAL(RCQS_SPATIAL_PTR),SE,%VAL(XA_PTR),
     '    XE,XG,%VAL(XIQ_PTR),
     '    %VAL(XP_PTR),%VAL(XQ_PTR),%VAL(YG_PTR),%VAL(YP_PTR),
     '    %VAL(YQ_PTR),%VAL(YQS_PTR),%VAL(ZA_PTR),
     '    ZE,ZG,%VAL(ZP_PTR),
     '    STRING,ERROR,*9999)
C news 24Jan2004: Calling solve to evaluate cell through FEM.f
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
        CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)
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
        CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)        

C newe VJ
      ENDIF !gauss stress with grid coupling with non grid at gauss scheme

C newe VJ 13Jan2004

C *** Main Gauss point loop
      DO 50 ng=1,NGT(NBH(NH_LOC(1,nx)))
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(/'' Gauss pt '',I3)') NG
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' ------------''/)')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF !DOP

C       Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
C!!     18Sep88: Eventually, replace following 6 stmts with
C!!      call xgmg(KTYP53(nr),...) where xgmg is modified so that DXIX
C!!      are derivs of Xi wrt Xj if 1st arg=1 or wrt Nu if >1
        IF(KTYP53(nr).EQ.1) THEN
C         stresses are referred to Xj in constitutive law.
          JP=0
        ELSE IF(KTYP53(nr).GE.2) THEN
C         stresses referred to Nu in constitutive law.
          JP=1
        ENDIF !KTYP53(nr)
C       Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C       derivs (DXIX) of Xi wrt Xj (JP=0) or Nu (JP=1) coords.
        CALL XGMG(JP,NIT(NBJ(1)),NBJ(1),nr,DXIX,GXL,GXU,RGX(ng),XG,
     '    ERROR,*9999)
C       Calculate the Jacobian for integration wrt undef coords:
        RWG=RGX(ng)*WG(ng,NBH(NH_LOC(1,nx)))
        IF(JTYP4.EQ.2) RWG=RWG*2.0d0*PI*XG(1,1)    !cyl symm about x
        IF(JTYP4.EQ.3) RWG=RWG*2.0d0*PI*XG(2,1)    !cyl symm about y
        IF(JTYP4.EQ.4) RWG=RWG*4.0d0*PI*XG(1,1)**2 !spherical symm
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' Integrating element resids wrt '
     '      //'undeformed coords: RWG='',D12.4)') RWG
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF !DOP

        IF(KTYP53(nr).EQ.3) THEN !Active stress component included
C         Get derivs of Xi wrt deformed Nu coords, DXIZN
          CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '      DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,'Fibre',ERROR,*9999)

C         Calculate derivs of deformed Nu wrt undeformed Nu (DZNXN)
C         (note that DXIX is derivs of Xi wrt undef Nu from XGMG above)
          DO ni=1,NITB
            DO mi=1,NITB
              SUM1=0.0d0
              DO k=1,NITB
                SUM1=SUM1+DZNXI(ni,k)*DXIX(k,mi)
              ENDDO !k
              DZNXN(ni,mi)=SUM1
            ENDDO !mi
          ENDDO !ni
C new VYW/MPN 12Apr2011
C DXNZN,DET_F_NU not needed in ZGTG5A 	
C          CALL INVERT(NITB,DZNXN,DXNZN,DET_F_NU)
        ENDIF !KTYP53(nr).EQ.3 (active stresses)

C JWF 2/5/03  Find ZGA values from ZEA for inertia problems
        IF (KTYP5I(nr).EQ.1) THEN ! inertia
          CALL ZEZG(1,NBH,ng,NHE,nx,DXIX,PG,ZEA,ZGA,ERROR,*9999) 
        ENDIF
C       Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
        IF(KTYP56(nr).EQ.6) THEN      !linear viscous relation
C         ..call to ZEZG with JP=3 computes 2nd derivs also
          CALL ZEZG(3,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
        ELSE! IF(KTYP56(nr).NE.6) THEN !all others
          CALL ZEZG(1,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
        ENDIF !ktyp56
C KFA 12/08/2004 need to check determinant of deformation tensor to
C   guard against inverting. A warning is printed becuase at the
C   moment the stress calculations in ZGTG53 does not multiply the 
C   hydrostatic pressure contribution by volume (it explicitly assumes
C   that J==1) - (see Bonet and Wood, 1997, pp. 127) This would seem
C   to be a poor assumption when the volume is negative as the pressure
C   is acting in the wrong way.

        IF(KTYP51(nr).EQ.3) THEN
          CALL DEFMGRADRC(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,DZDX,
     '      PG,XG,XI,ZE,ZG_temp,ERROR,*9999)
            Volume = DET(DZDX)
          IF(Volume .LT. 0.0D0) THEN
            WRITE(OP_STRING,'('' >>Warning: Volume at ng='',I5,'
     '        //''' ne='',I5,'' less than zero'')') ng,ne
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' DET(DZDX)='',D12.4)') Volume
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            !Volume = 0 ! GR ensure volume isn't -ve, probably not useful since
                        ! if this happens the solution has probably already failed.
          ENDIF
C           WRITE(OP_STRING,'('' DET(DZDX)='',D12.4)') Volume
C           CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

C       Calculate deformed metric tensors wrt Nu (AZL,AZU)
        CALL ZGMG(NBH(NH_LOC(1,nx)),nr,AZ,AZL,AZU,ZG,ERROR,*9999)

C news VJ 31Jan2004: Adding code to provide more efficient gauss point
C stress with grid coupling functionality.
        IF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.2)) THEN
C VJ 1Feb2004: 
C Added ZGTG53ATGRID call to update current grid/gauss point with green strain values
          CALL ZGTG53GRIDFROMGAUS(AZL,%VAL(ICQS_SPATIAL_PTR),
     '      %VAL(IRCQS_SPATIAL_PTR),ne,ng,%VAL(NQLIST_PTR),
     '      %VAL(NQNE_PTR),%VAL(RCQS_SPATIAL_PTR),ERROR,*9999)
C evaluate cellml file for 1 grid/gauss point
          CALL ZGTG53EVALCELL(%VAL(CELL_ICQS_VALUE_PTR),
     &      %VAL(CELL_RCQS_VALUE_PTR),
     &      %VAL(ICQS_SPATIAL_PTR),
     &      %VAL(IICQS_SPATIAL_PTR),
     &      %VAL(IRCQS_SPATIAL_PTR),
     &      ne,ng,%VAL(NQNE_PTR),
     &      %VAL(RCQS_SPATIAL_PTR),
     &      %VAL(YQS_PTR),ERROR,*9999)
C NOTE:   Stresses at grid points currently YQS array 2,3,4,5,6,7
C         update stresses from YQS to YG for 1 gauss/grid point
          CALL ZGTG53GAUSFROMGRID(ne,ng,%VAL(NGLIST_PTR),
     &      %VAL(NQLIST_PTR),%VAL(NQNE_PTR),%VAL(YG_PTR),%VAL(YQS_PTR),
     &      ERROR,*9999)

        ENDIF
C newe VJ 2Feb2004
C       Get contravariant cpts of 2nd Piola-Kirchhoff stress
C       tensor (TG) wrt undeformed Nu coordinates
        IF(KTYP51(nr).EQ.1) THEN      !plane stress
          CALL ZGTG51(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,
     '      CG(1,ng),RI1,RI2,RI3,TG,YG(1,ng),ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.2) THEN !plane strain
          CALL ZGTG52(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,
     '      CG(1,ng),RI1,RI2,RI3,TG,YG(1,ng),ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.3) THEN !3D
          CALL ZGTG53(STRESSTYPE,NBH(NH_LOC(1,nx)),
     '      nr,nx,AXU,AZ,AZL,AZU,
     '      CG(1,ng),EG,RI1,RI2,RI3,TG,XG,YG(1,ng),ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.4) THEN !membrane
          CALL ZGTG54(NBH(NH_LOC(1,nx)),nr,AXU,AZ,AZL,AZU,
     &      CG(1,ng),EG,RI1,RI2,RI3,TG,YG(1,ng),ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.5) THEN !string
          CALL ZGTG55(nr,AZL,CG(1,ng),EG,TG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.6) THEN !shell
        ENDIF !KTYP51(nr)

        FEXT(1,ng)=DSQRT(AZL(1,1))
        IF(KTYP53(nr).EQ.3) THEN !Active stress component included
C
C     OR 15-08-06
C
C     Changes to determine the proper active stress value depending 
C     on the definitions in IPACTI. KTYP59==3 got introduced. It allows
C     the user to specify a particular CellML variable to be added
C     to a specified component of the stress tensor
C
          IF (KTYP59(nr).EQ.3) THEN !Active stress component defined
                                ! within CellML
            CALL EVALASC(ne,ng,%VAL(NQNE_PTR),%VAL(YQS_PTR)
     &           ,%VAL(RCQS_SPATIAL_PTR),%VAL(ICQS_SPATIAL_PTR)
     &           ,ASC_ARRAYNAME(nr) ,ASC_CELLVARINDEX(nr),ACTIVE_STRESS
     &           ,ERROR,*9999)
          ELSE
            ACTIVE_STRESS = YG(1,ng)
          ENDIF
C     09-Dec-1989: NOTE: Don't have to define a separate face array
          CALL ZGTG5A(NBH(NH_LOC(1,nx)),nr,FEXT(1,ng),DXNZN,DZNXN,
     &         DET_F_NU,TG,TNA,ACTIVE_STRESS,ERROR,*9999)
        ENDIF !KTYP53(nr).EQ.3 (active stresses)

C ***   Main element residual
        IF(TYPE(1:7).EQ.'BIAXIAL') THEN !biaxial membrane testing

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' >>Using unrolled loops'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP
          dZ1_dNu1=ZG(1,NU1(1))
          dZ1_dNu2=ZG(1,NU1(2))
          dZ2_dNu1=ZG(2,NU1(1))
          dZ2_dNu2=ZG(2,NU1(2))
          nh1=NH_LOC(1,nx)
          nh2=NH_LOC(2,nx)
          nb=NBH(nh1)
          DO ns=1,NST(nb)+NAT(nb)
            PGA1 = PG(ns,2,ng,nb)*DXIX(1,1) +
     '             PG(ns,4,ng,nb)*DXIX(2,1)
            PGA2 = PG(ns,2,ng,nb)*DXIX(1,2) +
     '             PG(ns,4,ng,nb)*DXIX(2,2)
            AGE1 = (TG(1,1)*dZ1_dNu1+TG(1,2)*dZ1_dNu2)*PGA1 +
     '             (TG(2,1)*dZ1_dNu1+TG(2,2)*dZ1_dNu2)*PGA2
            AGE2 = (TG(1,1)*dZ2_dNu1+TG(1,2)*dZ2_dNu2)*PGA1 +
     '             (TG(2,1)*dZ2_dNu1+TG(2,2)*dZ2_dNu2)*PGA2
            RE(ns,nh1)=RE(ns,nh1)+AGE1*RWG*SE(ns,nb,ne)
            RE(ns,nh2)=RE(ns,nh2)+AGE2*RWG*SE(ns,nb,ne)
          ENDDO !ns

        ELSE IF(TYPE(1:4).EQ.'RC3D') THEN     !3D rect.cart.

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' >>Using unrolled loops'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr) !loop over geom dep vars
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh) !basis function for dep var
            DO ns=1,NST(nb)+NAT(nb) !element variables

C             PPGG(1) = PG(ns,1,ng,nb)
              PPGG(2) = PG(ns,NU1(1),ng,nb)*DXIX(1,1) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,1) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,1)
              PPGG(3) = PG(ns,NU1(1),ng,nb)*DXIX(1,2) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,2) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,2)
              PPGG(4) = PG(ns,NU1(1),ng,nb)*DXIX(1,3) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,3) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,3)

              AGE=(TG(1,1)*ZG(nhx,2)+TG(1,2)*ZG(nhx,4)+TG(1,3)*
     '          ZG(nhx,7))*PPGG(2)
     '          +(TG(2,1)*ZG(nhx,2)+TG(2,2)*ZG(nhx,4)+TG(2,3)*ZG(nhx,7))
     '          *PPGG(3)
     '          +(TG(3,1)*ZG(nhx,2)+TG(3,2)*ZG(nhx,4)+TG(3,3)*ZG(nhx,7))
     '          *PPGG(4)
C             Main residual
              RE(ns,nh)=RE(ns,nh)+AGE*RWG*SE(ns,nb,ne)
            ENDDO !ns
          ENDDO !nhx

        ELSE IF(TYPE(1:7).EQ.'PROLATE') THEN  !3D prolate spheroidal

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' >>Using unrolled loops'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP
          SLZ=SINH(ZG(1,1))
          SMZ=SIN (ZG(2,1))
          CLZ=DSQRT(1.0d0+SLZ*SLZ)
          CMZ=DSQRT(1.0d0-SMZ*SMZ)
          CSLZ=CLZ/SLZ
          CSMZ=CMZ/SMZ
          G1=SLZ*SLZ+SMZ*SMZ
          E1=CLZ*SLZ/G1
          E2=CMZ*SMZ/G1

          nb=NBH(NH_LOC(1,nx))
          DO ns=1,NST(nb)+NAT(nb)
            PPG(ns,1) = PG(ns,1,ng,nb)
            PPG(ns,2) = PG(ns,NU1(1),ng,nb)*DXIX(1,1) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,1) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,1)
            PPG(ns,3) = PG(ns,NU1(1),ng,nb)*DXIX(1,2) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,2) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,2)
            PPG(ns,4) = PG(ns,NU1(1),ng,nb)*DXIX(1,3) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,3) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,3)
          ENDDO !ns

          D11=ZG(1,NU1(1))
          D12=ZG(1,NU1(2))
          D13=ZG(1,NU1(3))
          D21=ZG(2,NU1(1))
          D22=ZG(2,NU1(2))
          D23=ZG(2,NU1(3))
          D31=ZG(3,NU1(1))
          D32=ZG(3,NU1(2))
          D33=ZG(3,NU1(3))

          DLA1=ZG(1,NU1(1))
          DLA2=ZG(1,NU1(2))
          DLA3=ZG(1,NU1(3))
          DMA1=ZG(2,NU1(1))
          DMA2=ZG(2,NU1(2))
          DMA3=ZG(2,NU1(3))
          DTA1=ZG(3,NU1(1))
          DTA2=ZG(3,NU1(2))
          DTA3=ZG(3,NU1(3))

          nh1=NH_LOC(1,nx)
          nh2=NH_LOC(2,nx)
          nh3=NH_LOC(3,nx)
          DO ns=1,NST(nb)+NAT(nb)
            PGG =PPG(ns,1)
            PGA1=PPG(ns,2)
            PGA2=PPG(ns,3)
            PGA3=PPG(ns,4)

            AGE1 = TG(1,1)*(D11*(PGA1-(E1*DLA1+E2*DMA1)*PGG) +
     '        D21*(E1*DMA1-E2*DLA1)*PGG+D31*E1*SMZ*SMZ*DTA1*PGG) +
     '        TG(2,1)*(D11*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '        D21*(E1*DMA2-E2*DLA2)*PGG+D31*E1*SMZ*SMZ*DTA2*PGG) +
     '        TG(3,1)*(D11*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '        D21*(E1*DMA3-E2*DLA3)*PGG+D31*E1*SMZ*SMZ*DTA3*PGG) +
     '        TG(1,2)*(D12*(PGA1-(E1*DLA1+E2*DMA1)*PGG) +
     '        D22*(E1*DMA1-E2*DLA1)*PGG+D32*E1*SMZ*SMZ*DTA1*PGG) +
     '        TG(2,2)*(D12*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '        D22*(E1*DMA2-E2*DLA2)*PGG+D32*E1*SMZ*SMZ*DTA2*PGG) +
     '        TG(3,2)*(D12*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '        D22*(E1*DMA3-E2*DLA3)*PGG+D32*E1*SMZ*SMZ*DTA3*PGG) +
     '        TG(1,3)*(D13*(PGA1-(E1*DLA1+E2*DMA1)*PGG) +
     '        D23*(E1*DMA1-E2*DLA1)*PGG+D33*E1*SMZ*SMZ*DTA1*PGG) +
     '        TG(2,3)*(D13*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '        D23*(E1*DMA2-E2*DLA2)*PGG+D33*E1*SMZ*SMZ*DTA2*PGG) +
     '        TG(3,3)*(D13*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '        D23*(E1*DMA3-E2*DLA3)*PGG+D33*E1*SMZ*SMZ*DTA3*PGG)

            AGE2 = TG(1,1)*(D11*(E2*DLA1-E1*DMA1)*PGG +
     '        D21*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D31*E2*SLZ*SLZ*DTA1*PGG)+
     '        TG(2,1)*(D11*(E2*DLA2-E1*DMA2)*PGG +
     '        D21*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D31*E2*SLZ*SLZ*DTA2*PGG)+
     '        TG(3,1)*(D11*(E2*DLA3-E1*DMA3)*PGG +
     '        D21*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D31*E2*SLZ*SLZ*DTA3*PGG)+
     '        TG(1,2)*(D12*(E2*DLA1-E1*DMA1)*PGG +
     '        D22*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D32*E2*SLZ*SLZ*DTA1*PGG)+
     '        TG(2,2)*(D12*(E2*DLA2-E1*DMA2)*PGG +
     '        D22*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D32*E2*SLZ*SLZ*DTA2*PGG)+
     '        TG(3,2)*(D12*(E2*DLA3-E1*DMA3)*PGG +
     '        D22*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D32*E2*SLZ*SLZ*DTA3*PGG)+
     '        TG(1,3)*(D13*(E2*DLA1-E1*DMA1)*PGG +
     '        D23*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D33*E2*SLZ*SLZ*DTA1*PGG)+
     '        TG(2,3)*(D13*(E2*DLA2-E1*DMA2)*PGG +
     '        D23*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D33*E2*SLZ*SLZ*DTA2*PGG)+
     '        TG(3,3)*(D13*(E2*DLA3-E1*DMA3)*PGG +
     '        D23*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D33*E2*SLZ*SLZ*DTA3*PGG)

            AGE3 = TG(1,1)*(-(D11*CSLZ+D21*CSMZ)*DTA1*PGG +
     '        D31*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '        TG(2,1)*(-(D11*CSLZ+D21*CSMZ)*DTA2*PGG +
     '        D31*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '        TG(3,1)*(-(D11*CSLZ+D21*CSMZ)*DTA3*PGG +
     '        D31*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG)) +
     '        TG(1,2)*(-(D12*CSLZ+D22*CSMZ)*DTA1*PGG +
     '        D32*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '        TG(2,2)*(-(D12*CSLZ+D22*CSMZ)*DTA2*PGG +
     '        D32*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '        TG(3,2)*(-(D12*CSLZ+D22*CSMZ)*DTA3*PGG +
     '        D32*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG)) +
     '        TG(1,3)*(-(D13*CSLZ+D23*CSMZ)*DTA1*PGG +
     '        D33*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '        TG(2,3)*(-(D13*CSLZ+D23*CSMZ)*DTA2*PGG +
     '        D33*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '        TG(3,3)*(-(D13*CSLZ+D23*CSMZ)*DTA3*PGG +
     '        D33*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG))

            RE(ns,nh1)=RE(ns,nh1)+AGE1*RWG*SE(ns,nb,ne)
            RE(ns,nh2)=RE(ns,nh2)+AGE2*RWG*SE(ns,nb,ne)
            RE(ns,nh3)=RE(ns,nh3)+AGE3*RWG*SE(ns,nb,ne)

          ENDDO !ns

        ELSE !all other cases

          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh)
            DO ns=1,NST(nb)+NAT(nb)
              PPG(ns,1)=PG(ns,1,ng,nb)
              DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj2,nr)
                PPG(ns,1+nj)=PGX(nb,nj,ns,DXIX,PG(1,1,ng,nb))
              ENDDO !njj2
            ENDDO !ns
            DO ns=1,NST(nb)+NAT(nb)
              AGE=AG(nb,nhx,nr,ns,PPG,TG,ZG)
              RE(ns,nh)=RE(ns,nh)+AGE*RWG*SE(ns,nb,ne)
C MPN 28Feb97 IF(DOP) THEN !too much DOP
              IF(DOP.AND.IWRIT4(nr,nx).GE.2) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' >>Residual calcs:'
     '            //' nh='',I2,'' ns='',I2,'' nb='',I2,'
     '            //''' AG='',D10.3,'' RE='',D10.3/,'
     '            //''' PG(ns,1,ng,nb):'',D10.3)')
     '            nh,ns,nb,AGE,RE(ns,nh),PG(ns,1,ng,nb)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF !DOP
            ENDDO !ns
          ENDDO !nhx

        ENDIF !TYPE
        
C JWF 2/5/03  Update acceleration (ACC) vector with inertia terms
        IF (KTYP5I(nr).EQ.1) THEN ! inertia
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr) !loop over geom dep vars
            nh=NH_LOC(nhx,nx)
            ACC(nh)=ZGA(nh,1)
         ENDDO                  !nhx        
        ELSE ! no inertia
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr) !loop over geom dep vars
            nh=NH_LOC(nhx,nx)
            ACC(nh)=0.0d0
          ENDDO !nhx 
        ENDIF  !KTYP5I      

C ***   Gravity
        IF(TYPE(1:4).EQ.'RC3D') THEN     !3D rect.cart.
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr) !loop over geom dep vars
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh) !basis function for dep var
            DO ns=1,NST(nb)+NAT(nb) !element variables
              AGE=CG(IL_density,ng)*(b_inc(nh)+acc(nh))*PG(ns,1,ng,nb)
C             Main residual
              RE(ns,nh)=RE(ns,nh)+AGE*RWG*SE(ns,nb,ne)
            ENDDO !ns
          ENDDO !nhx
        ENDIF

C ***   Incompressibility (+ fluid) constraint(s) in 3D case
C KAT .AND. has higher precedence than .OR..  Is incomp+* always 3D?
        IF(KTYP51(nr).EQ.3.AND.KTYP52(nr).EQ.2.         !3D + incomp
     '                      OR.KTYP52(nr).EQ.3.         !incomp+fluid
     '                      OR.KTYP52(nr).EQ.5) THEN    !incomp+inext
          IF(KTYP52(nr).EQ.2) THEN !standard incomp constraint
            NBP=NBH(NH_LOC(NH_LOC(0,nx),nx))   !basis fn for press var
C            Incomp_resid=DSQRT(RI3)-1.0d0
C KFA 12/08/2004 Changing the incompressible constraint to only
C   have one minima at 1.0 rather than two at +/- 1.0
            Incomp_resid=Volume-1.0d0

          ELSE IF(KTYP52(nr).EQ.3) THEN !incomp + fluid constraint
            NBP=NBH(NH_LOC(NH_LOC(0,nx),nx))   !basis fn for press var
C           Calculate deformed Christoffel symbols wrt undef Nu coords
C           NOTE: ZG needs derivs wrt Nu not Xi !
C KAT 1Nov00: X3G should not be used in TOFFEL as it is not set up.
            CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF,DBM,AZU,ZG,%VAL(0),
     '        .FALSE.,ERROR,*9999)
C            CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF,DBM,AZU,ZG,X3G,
C     '        .FALSE.,ERROR,*9999)
C           Calc del-squared(p) wrt deformed coords, where p is the
C           hydrostatic pressure that varies with Xi(3) only.
            SUM=0.0d0
            DO j=1,NITB
              DO k=1,NITB
                SUM=SUM+CHTOFF(3,j,k)*AZU(j,k)
              ENDDO !k
            ENDDO !j
            delsqP=0.0d0
            DO ns=1,NST(NBP)+NAT(NBP)
              delsqP=delsqP+ZE(ns,NH_LOC(0,nx))*
     '          (PG(ns,8,ng,NBP)*AZU(3,3)-SUM*PG(ns,7,ng,NBP))
            ENDDO !ns
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
C new MPN 18Mar98: added DOP
              WRITE(OP_STRING,'('' fluid cond='',D12.4,'' DT='',D12.4)')
     '          CG(IL_fluid_conductivity,ng),DT
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C end new
              WRITE(OP_STRING,'('' del-squared(pressure)='',D12.4)')
     '          delsqP
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF !DOP
C ***       Calculate residuals associated with Darcy's Law
C new MPN 18Mar98: fixed error
            Incomp_resid=CG(IL_fluid_conductivity,ng)*DT*delsqP
     '        +DSQRT(RI3)-1.d0
C old            Incomp_resid=CG(IL_fluid_conductivity,ng)*DT*delsqP
C old     '        +(DSQRT(RI3)-1.d0)/DSQRT(RI3)

          ELSE IF(KTYP52(nr).EQ.5) THEN !incomp + inextensible
            NBP=NBH(NH_LOC(NH_LOC(0,nx)-1,nx)) !basis fn for press var
            IF(DOP) THEN
              WRITE(OP_STRING,'('' RI3='',D12.4)') RI3
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' FEXT(1,ng)='',D12.4)') FEXT(1,ng)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
C!!! Do we know that this is 3D?
            Incomp_resid=Volume-1.0d0
            Inext_resid =DABS(FEXT(1,ng)-1.0d0)
          ENDIF !incomp/incomp+fluid/incomp+inext

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' SQRT(RI3)='',D12.4)') DSQRT(RI3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Incomp_resid='',D12.4)') Incomp_resid
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP

C CS 25/8/99 slight fix for PJH's changes for RD
C         Galerkin incompressibility (+ fluid or + inext) residuals
C          DO ns=1,NST(NBP)+NAT(NBP)
C PJH 3Mar99   RE(ns,NH_LOC(NH_LOC(0,nx),nx))=
C     '        RE(ns,NH_LOC(NH_LOC(0,nx),nx))+
C            RE(ns,NH_LOC(NH_LOC(4,nx),nx))=   !nhx=4 is pressure var
C     '        RE(ns,NH_LOC(NH_LOC(4,nx),nx))+
C     '        Incomp_resid*PG(ns,1,ng,NBP)*RWG
C            IF(KTYP52(nr).EQ.5) THEN !incomp + inextensible
C              RE(ns,NH_LOC(NH_LOC(5,nx),nx))=   !nhx=5 is tension var
C     '          RE(ns,NH_LOC(NH_LOC(5,nx),nx))+
C     '          Inext_resid*PG(ns,1,ng,NBP)*RWG
C            ENDIF
C          ENDDO !ns
C
C         Galerkin incompressibility (+ fluid or + inext) residuals
          DO ns=1,NST(NBP)+NAT(NBP)
            IF(.NOT.KTYP52(nr).EQ.5) THEN
              RE(ns,NH_LOC(NH_LOC(0,nx),nx))=
     '          RE(ns,NH_LOC(NH_LOC(0,nx),nx))+
     '          Incomp_resid*PG(ns,1,ng,NBP)*RWG
            ELSE !incomp + inextensible
              RE(ns,NH_LOC(NH_LOC(4,nx),nx))= !nhx=4 is pressure var
     '          RE(ns,NH_LOC(NH_LOC(4,nx),nx))+
     '          Incomp_resid*PG(ns,1,ng,NBP)*RWG
              RE(ns,NH_LOC(NH_LOC(5,nx),nx))= !nhx=5 is tension var
     '          RE(ns,NH_LOC(NH_LOC(5,nx),nx))+
     '          Inext_resid*PG(ns,1,ng,NBP)*RWG
            ENDIF
          ENDDO !ns

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' NBP='',I2)') NBP
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' PG(1..,1,ng,NBP): '',5D10.3,'
     '        //'/(19X,5D10.3))')
     '        (PG(ns,1,ng,NBP),ns=1,NST(NBP)+NAT(NBP))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP

C old MPN 17Mar97: use standard Galerkin resids (above)
C          IF(NODEPRESS.AND.ELEMPRESS) THEN
CC           node and element based hydrostatic pressure interp
CC           NOTE: pressure basis function is not used to
CC                 weight the residual here (non-Galerkin)
C            DO ns=1,NST(NBP)+NAT(NBP)
C              RE(ns,NH_LOC(NH_LOC(0,nx),nx))=
C     '          RE(ns,NH_LOC(NH_LOC(0,nx),nx))+Incomp_resid*RWG
C            ENDDO !ns
C          ELSE IF(NODEPRESS) THEN
CC           Node based hydrostatic pressure interpolation.
CC           Weight with Xi3=0 face basis (ie no Xi3 variation for wght)
C            IFE=5
C            nf=NFF(IFE) !face number for Xi3=0 face in current elem
C            NBFF=NBJF(3,nf) !basis fn for third geom var on face nf
CC                         !!!WARNING: Need to store pressure face bases
CC                         !!!in NPF
C            ng_face=MOD(ng-1,NGAP(1,NBP)*NGAP(2,NBP))+1 !Xi3=0 face GP#
C            DO ns=1,NST(NBP)
CC             Determine Xi3=0 face vertex for element vertex ns
C              nsf=MOD(ns-1,NST(NBFF))+1
C              RE(ns,NH_LOC(NH_LOC(0,nx),nx))=
C     '          RE(ns,NH_LOC(NH_LOC(0,nx),nx))+
C     '          Incomp_resid*PG(nsf,1,ng_face,NBFF)*RWG
C            ENDDO !ns
C            IF(DOP) THEN
CC$            call mp_setlock()
C              WRITE(OP_STRING,'('' NBP='',I2,'' nbff='',I2,'
C     '          //''' ng_face='',I5)') NBP,NBFF,ng_face
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              WRITE(OP_STRING,'('' PG(1..,1,ng_face,nbff): '','
C     '          //'5D10.3,/(25X,5D10.3))')
C     '          (PG(MOD(ns-1,NST(NBFF))+1,1,ng_face,NBFF),
C     '          ns=1,NST(NBP))
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
C            ENDIF !DOP
C          ELSE
CC           element based hyd. pressure interp
CC           (standard Galerkin)
C            DO ns=1,NST(NBP)+NAT(NBP)
C              RE(ns,NH_LOC(NH_LOC(0,nx),nx))=
C     '          RE(ns,NH_LOC(NH_LOC(0,nx),nx))+
C     '          Incomp_resid*PG(ns,1,ng,NBP)*RWG
C            ENDDO !ns
C            IF(DOP) THEN
CC$            call mp_setlock()
C              WRITE(OP_STRING,'('' NBP='',I2)') NBP
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              WRITE(OP_STRING,'('' PG(1..,1,ng,NBP): '',5D10.3,'
C     '          //'/(19X,5D10.3))')
C     '          (PG(ns,1,ng,NBP),ns=1,NST(NBP)+NAT(NBP))
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
C            ENDIF !DOP
C          ENDIF !NODEPRESS/ELEMPRESS

C ***   Compressibility + fluid in 3D case
        ELSE IF(KTYP51(nr).EQ.3.AND.KTYP52(nr).EQ.4) THEN !3D + compressible
C TVK 12/01/1999 Correct equilibrium pressure term
          IF(EQUIM_PRESSURE) THEN
            IF(ITYP10(nr).EQ.1) THEN
              NITB=NIT(nb)
              DO i=1,3
                DO j=1,3
                  AXU(i,j)=0.0d0
                ENDDO
                AXU(i,i)=1.0d0
              ENDDO
              CALL ZGMG(NBH(NH_LOC(1,nx)),nr,AZ,AZL,
     '          AZU,ZG,ERROR,*9999)
              IF(KTYP53(nr).GT.1) THEN
                EG13=AZL(1,3)/2.0d0
                EG12=AZL(1,2)/2.0d0
C                RK1=(AZL(1,1)-1.0d0)/2.0d0
C                RK2=EG13*EG13+EG12*EG12
                BG(1,1)=3.0d0-AZL(1,1)
                BG(2,2)=3.0d0-AZL(2,2)
                BG(3,3)=3.0d0-AZL(3,3)
                BG(2,1)=   -AZL(2,1)
                BG(3,1)=   -AZL(3,1)
                BG(3,2)=   -AZL(3,2)
                CALL ENERGY(nr,CG(1,1),DW,3.0d0,3.0d0,1.0d0,0.0d0,
     '            0.0d0,0.0d0,YG(1,ng),ERROR,*9999)
                IF(KTYP52(nr).EQ.4) THEN !compressible + fluid
C TVK no compliance on I3 term  W3=RI3*DW(3)*ZG(NH_LOC(0,nx),1)
                  W3=1.0d0*DW(3)
                ENDIF
                DO mj=1,NITB
                  DO nj=1,mj
                    TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj)
     '                + W3*AZU(mj,nj))
                  ENDDO
                ENDDO
C         Note: DW(4) & DW(5) are zero for isotropic case
                TG(3,1)=TG(3,1)+EG13*DW(5)
                TG(2,1)=TG(2,1)+EG12*DW(5)
                TG(1,1)=TG(1,1)+AXU(1,1)*DW(4)
                TG(1,2)=TG(2,1)
                TG(1,3)=TG(3,1)
                TG(2,3)=TG(3,2)
              ENDIF
            ENDIF
            TG0=TG(1,1)/AXU(1,1)
            EQUIM_PRESSURE=.FALSE.
          ENDIF


          nh_pressure=NH_LOC(NH_LOC(0,nx),nx)
          nb_pressure=NBH(nh_pressure)   !basis fn for press var

C TVK 07/01/2000 Fixing residual term for compression of tissue
C          Comp_resid=(CG(IL_compliance,ng)*ZG(nh_pressure,1)
C     '               -(DSQRT(RI3)-1.d0))*RWG
          Comp_resid=(ZG(nh_pressure,1)-TG0-CG(IL_compliance,ng)*
     '      (1.d0-DSQRT(RI3)))*RWG

          IF(DOP) THEN
            WRITE(OP_STRING,'('' RI3='',D12.4)') RI3
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ZG(nh_pressure,1)='',D12.4)')
     '        ZG(nh_pressure,1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Comp_resid='',D12.4)') Comp_resid
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          DO ns=1,NST(nb_pressure)+NAT(nb_pressure)
            RE(ns,nh_pressure)=RE(ns,nh_pressure)+Comp_resid
          ENDDO

C ***   External press loads in case of 2D or 3D membrane
        ELSE IF(KTYP51(nr).EQ.4.AND.NW.GT.1) THEN !membrane
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
            nh=NH_LOC(nhm,nx)
            nb=NBH(nh)
            SUM=0.0d0
            DO ni=1,NIT(nb)
              SUM=SUM+AZU(ni,3)*ZG(nhx,NU1(ni))
            ENDDO
            DO ns=1,NST(nb)+NAT(nb)
              RE(ns,nh)=RE(ns,nh)+PF(1)*SUM*PG(ns,1,ng,nb)*RWG
            ENDDO
          ENDDO

C ***   External press loads in case of 2D string
        ELSE IF(KTYP51(nr).EQ.5) THEN !2D string
          IF(NRT.EQ.2) THEN !special case for coupled sail-flow problem
            PF(1)=PRESS_DIFF_AERO(ng,ne-NET(1))
          ENDIF
C         For x residual use ZG(2,2)=dy/dNu
          nb=NBH(NH_LOC(1,nx))
          DO ns=1,NST(nb)+NAT(nb)
            RE(ns,1)=RE(ns,1)+PF(1)*ZG(2,2)*PG(ns,1,ng,nb)*RWG
          ENDDO
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' ne='',I4,'' ng='',I3,'
     '        //''' PF(1)='',D12.3,'' ZG(2,2)='',D12.3,'
     '        //''' RE(ns,1):'',6D12.3)')
     '        ne,ng,PF(1),ZG(2,2),(RE(ns,1),ns=1,NST(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
C         For y residual use ZG(1,2)=dx/dNu. Note -ve.
          nb=NBH(NH_LOC(2,nx))
          DO ns=1,NST(nb)+NAT(nb)
            RE(ns,2)=RE(ns,2)-PF(1)*ZG(1,2)*PG(ns,1,ng,nb)*RWG
          ENDDO
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' ne='',I4,'' ng='',I3,'
     '        //''' PF(1)='',D12.3,'' ZG(1,2)='',D12.3,'
     '        //''' RE(ns,2):'',6D12.3)')
     '        ne,ng,PF(1),ZG(1,2),(RE(ns,2),ns=1,NST(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDIF !special cases

 50   CONTINUE !end of ng loop
C new MPN 26Mar97
      IF(IWRIT4(nr,nx).GE.1.AND.DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        WRITE(OP_STRING,
     '    '(/'' CPU time for all Gauss pt calcs:           '','
     '    //'D11.4)') ELAPSED_TIME
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF !IWRIT4(nr,nx).GE.1.AND.DOP

C *** External pressure loads in 3D case
      IF(KTYP51(nr).EQ.3) THEN !3D
        NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for hydr. pressure
        IF(KTYP52(nr).GE.2.AND.KTYP57(nr).GT.1) THEN !incomp; p entered
C         Boundary pressure constraint equations
C new MPN 17Mar97: If high order model with nodal pressure interp
C                  only apply specific pressure bcs if KTYP5A(nr)=2
C                  If low order model with element pressure interp
C                  apply specific pressure bcs if there are 2 free
C                  hyd pressure params [ie NSP(1) and NSP(2) non-zero]
          IF(NODEPRESS) THEN
C           node based hydrostatic pressure interp
            IF(KTYP5A(nr).EQ.2) THEN !match norm stress on Xi3 faces
C             Use explicit boundary press constraint eqns to determine
C             free hydrostatic pressure variables (these replace the
C             incomp constraints above).
              CALL PFRE_NP(IBT,IDO,INP,NAN,NBH,NBJ,NBJF,NFF,NGAP,NHE,
     '          NKEF,NNF,NPNE(1,1,ne),nr,ne,NRE,NW,nx,NXI(-NIM,0,ne),
     '          CE,CP,FEXT,PF,PG,RE(1,NH_LOC(NH_LOC(0,nx),nx)),
     '          WG,XE,XG,YG,ZE,ZG,ERROR,*9999)
            ENDIF !KTYP5A(nr).EQ.2
          ELSE IF(ELEMPRESS) THEN
C           element based hydrostatic pressure interp
            IF(NSP(1).NE.0.AND.NSP(2).NE.0) THEN !na #s for press vars
C             Determine hydrostatic pressure vars using face press bcs
C             instead of incompressibility constraints
              CALL PFRE_NE(IBT,IDO,INP,NAN,NBH,NBJ,NFF,NGAP,NHE,
     '          NPNE(1,1,ne),nr,ne,NRE,NSP,NW,nx,NXI(-NIM,0,ne),
     '          CE,CP,FEXT,PF,PG,RE(1,NH_LOC(NH_LOC(0,nx),nx)),
     '          XE,XG,YG,ZE,ZG,ERROR,*9999)
            ENDIF !NSP(1).NE.0.AND.NSP(2).NE.0
          ENDIF !NODEPRESS/ELEMPRESS
        ENDIF !KTYP52(nr).GE.2.AND.KTYP57(nr).GT.1

C       Contribution of external pressure loads to stress equilibrium
C       equation residuals at nodes
        DO iface=1,2

C new MPN 16APR97: to prevent access violations
          IF(iface.EQ.1) THEN
            adj_dirn=-3
          ELSE
            adj_dirn=3
          ENDIF
          neadj=NXI(adj_dirn,1,ne)
          ADJ_XI3_ELEM=.FALSE.
          IF(neadj.NE.0) THEN
            IF(nr.EQ.NRE(neadj)) ADJ_XI3_ELEM=.TRUE.
          ENDIF
          IF(iface.EQ.1.AND..NOT.ADJ_XI3_ELEM
     '      .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext press bc on Xi3=0 face
     '      iface.EQ.2.AND..NOT.ADJ_XI3_ELEM
     '      .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !ext pres bc on Xi3=1 face
            IF(DABS(PF(iface)).GT.1.0D-10) THEN !non-zero press bc
              IXF=iface-1
              IFE=iface+4
              nf=NFF(IFE)
              IF(nf.NE.0) THEN
                CALL PFRF(IBT,IDO,INP,IXF,NAN,NBH,NBJF(1,nf),NGAP,
     '            nr,nx,PF(iface),PG,RF,WG,ZE,ZG,ERROR,*9999)
                DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                  nh=NH_LOC(nhx,nx)
                  nsf=0
                  NBFF=NBJF(nhx,nf)
                  NBE=NBH(nh)
                  DO nnbf=1,NNT(NBFF)
                    nn=NNF(1+nnbf,IFE,NBE)
                    DO nkbf=1,NKT(nnbf,NBFF)
                      nk=NKEF(nkbf,nnbf,IFE,NBE)
                      nse=nk+(NN-1)*NKT(nn,NBE)
                      nsf=nsf+1
                      RE(nse,nh)=RE(nse,nh)-RF(nsf,nh)*SE(nse,NBE,ne)
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' >>ZERE50 diagnostic op:'
     '                    //' RF('',I2,'','',I1,''):'',D10.3,'
     '                    //''' RE('',I2,'','',I1,''):'',D10.3)')
     '                    nsf,nh,RF(nsf,nh),
     '                    nse,nh,RE(nse,nh)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF !DOP
                    ENDDO !nkbf (nk)
                  ENDDO !nnbf (nn)
                ENDDO !nhx
              ENDIF !nf.NE.0
            ENDIF !PF(iface)>0
          ENDIF !ext press bc on face
        ENDDO !iface
      ENDIF !KTYP51(nr).EQ.3

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
C new MPN 26Mar97
        IF(IWRIT4(nr,nx).GE.1) THEN
          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-ELAPSED_TIME
          WRITE(OP_STRING,
     '      '(/'' CPU time for pressure bc calcs:           '','
     '      //'D11.4)') ELAPSED_TIME
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        WRITE(OP_STRING,'(/'' Residuals for element'',I4,'':'')') ne
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          WRITE(OP_STRING,'('' RE(ns,'',I1,''): '',8D11.3,'
     '      //'/(11X,8D11.3))')
     '      nh,(RE(ns,nh),ns=1,NST(NBH(nh))+NAT(NBH(nh)))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZERE50')
      RETURN
 9999 CONTINUE
C KAT 14May01: mp_setlock not OPENMP.
C              Critical section is not essential.
C      IF(IWRIT4(nr,nx).GE.1.AND.DOP) THEN
CC$      call mp_setlock()
CC$      call mp_unsetlock()
C      ENDIF !IWRIT4(nr,nx).GE.1.AND.DOP
      CALL ERRORS('ZERE50',ERROR)
      CALL EXITS('ZERE50')
      RETURN 1
      END


