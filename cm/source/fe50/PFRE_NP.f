      SUBROUTINE PFRE_NP(IBT,IDO,INP,NAN,NBH,NBJ,NBJF,NFF,NGAP,NHE,NKEF,
     '  NNF,NPNE,nr,ne,NRE,NW,nx,NXI,
     '  CE,CP,FEXT,PF,PG,RE,WG,XE,XW,YG,ZE,ZW,ERROR,*)

C#### Subroutine: PFRE_NP
C###  Description:
C###    PFRE_NP evaluates contribution to element residuals RE(ns,4)
C###    (associated with the hydrostatic pressure variable) for
C###    incompressible materials, arising from the stress constraint
C###    due to pressure bcs.

C**** MPN 11-Jan-95: This routine is used when the hydrostatic pressure
C****                is interpolated using nodal based variables.

      IMPLICIT NONE
      INCLUDE 'acti00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ptr00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NBJF(NJM,NFM),NFF(6),
     '  NGAP(NIM,NBM),NHE,NKEF(0:4,16,6,NBFM),NNF(0:17,6,NBFM),
     '  NPNE(NNM,NBFM),nr,ne,NRE(NEM),NW,nx,NXI(-NIM:NIM,0:NEIM)
      REAL*8 CE(NMM),CP(NMM,NPM),FEXT(NIFEXTM,NGM),PF(2),
     &  PG(NSM,NUM,NGM,NBM),RE(NSM),WG(NGM,NBM),XE(NSM,NJM),
     &  XW(NJM,NUM),YG(NIYGM,NGM),ZE(NSM,NHM),ZW(NHM,NUM),ACTIVE_STRESS
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER adj_dirn,iface,IFE,k,mi,mj,mz,NBFF,NBP,NCW,neadj,nf,
     '  ng,ng1,ng2,NGI1,NGI2,nhx,ni,NITB,nj,nk,nkbf,
     '  nn,nnbf,nse,nsf,nz
      PARAMETER (NCW=35) !CW must be dimen.d the same size as CE array
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CW(NCW),
     '  D(5,5),DETERM,DET_F_NU,
     '  DNUREFDNU(3,3),DNUREFDZ(3,3),DXIXN(3,3),DXIZN(3,3),DXNZN(3,3),
     '  DZDNU(3,3),DZDNUREF(3,3),DZNXI(3,3),DZNXN(3,3),
     '  EG(3,3),GXL(3,3),GXU(3,3),GZ,GZL(3,3),GZU(3,3),
     '  RGX_local,RI1,RI2,RI3,RWG,SUM,
     '  TCNU(3,3),TCNUREF33,TW(3,3),TWA,XI(3)
      LOGICAL ADJ_XI3_ELEM
      CHARACTER STRESSTYPE*17
      DATA D/5*0.0d0,-0.288675134594813d0,0.288675134594813d0,3*0.0d0,
     '       -0.387298334620741d0,0.0d0,0.387298334620741d0,2*0.0d0,
     '       -0.430568155797026d0,    -0.169990521792428d0,
     '        0.169990521792428d0,     0.430568155797026d0,  0.0d0,
     '       -0.453089922969332d0,    -0.269234655052841d0,  0.0d0,
     '        0.269234655052841d0,     0.453089922969332d0/
      DATA STRESSTYPE/' '/

      CALL ENTERS('PFRE_NP',*9999)

      CALL ASSERT(NCW.EQ.NMM,'>>Dimension of CW array (NCW)'
     '  //' must equal dimension of CE array (NMM)',ERROR,*9999)
      NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis function for pressure vars
      NITB=NIT(NBP)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>>PFRE_NP diagnostic op:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF !DOP

      DO iface=1,2
C new MPN 16APR97: to prevent access violations
        IF(iface.EQ.1) THEN
          adj_dirn=-3
        ELSE
          adj_dirn=3
        ENDIF
        neadj=NXI(adj_dirn,1)
        ADJ_XI3_ELEM=.FALSE.
        IF(neadj.NE.0) THEN
          IF(nr.EQ.NRE(neadj)) ADJ_XI3_ELEM=.TRUE.
        ENDIF
        IF(iface.EQ.1.AND..NOT.ADJ_XI3_ELEM
     '    .AND.(NW.EQ.2.OR.NW.EQ.4).OR.   !ext press bc on Xi3=0 face
     '     iface.EQ.2.AND..NOT.ADJ_XI3_ELEM
     '    .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !ext press bc on Xi3=1 face
C old
C        IF(iface.EQ.1.AND.(NXI(-3,1).EQ.0.OR.nr.NE.NRE(NXI(-3,1))).AND.
C     '    (NW.EQ.2.OR.NW.EQ.4).OR. !ext press bc applied on Xi3=0 face
C     '     iface.EQ.2.AND.(NXI( 3,1).EQ.0.OR.nr.NE.NRE(NXI( 3,1))).AND.
C     '    (NW.EQ.3.OR.NW.EQ.4)) THEN !ext press bc applied on Xi3=1 face

          IFE=iface+4
          nf=NFF(IFE)
          IF(nf.NE.0) THEN
            XI(3)=DBLE(iface-1)
            NBFF=NBJF(3,nf) !basis fn for third geometric variable
C                           !!!WARNING: Need to store pressure face bases
C                           !!!in NPF
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(/'' Xi(3)='',I1,'' face (nf='',I3,'').'
     '          //' PF ='',D12.4)') iface-1,nf,PF(iface)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Pressure face basis nbff='',I2)')
     '          nbff
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF !DOP

C MLB 16/4/97 initialising nsf
            nsf=0
C new MPN 17Mar97: only initialise resids on ext face, since incomp
C                  resids have been set up already.
C           Initialise face pressure residuals
            DO nnbf=1,NNT(NBFF)
              nn=NNF(1+nnbf,IFE,NBP)
              DO nkbf=1,NKT(nnbf,NBP)
                nk=NKEF(nkbf,nnbf,IFE,NBP)
                nsf=nsf+1
                nse=nk+(nn-1)*NKT(nn,NBP)
                RE(nse)=0.0d0
              ENDDO !nkbf (nk)
            ENDDO !nnbf (nn)

C           Loop over gauss pts on face
            ng=0
            NGI1=NGAP(1,NBP) !Gauss pts in Xi1 for pressure basis
            NGI2=NGAP(2,NBP) !Gauss pts in Xi2 for pressure basis
            DO ng2=1,NGI2
              DO ng1=1,NGI1
                ng=ng+1
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,
     '              '('' >>>PFRE_NP diag op at Gauss pt '',I2)') NG
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF !DOP
                XI(1)=0.5d0+D(ng1,NGI1)
                XI(2)=0.5d0+D(ng2,NGI2)
C               Calculate stress at XI (ie position
C               of current Gauss pt on face).
C               Interpolate material parameters at XI
                CALL CPXI(1,IBT,IDO,INP,NPNE,nr,nx,CE,CP,CW,XI,
     '            ERROR,*9999)
C               Interpolate midwall geom vars XW and derivs wrt Xi
                CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XW,XI,ERROR,*9999)
C               Calculate undef metric tensors wrt Xi (GXL,GXU) and
C               derivs of Xi wrt undef Nu coords, DXIXN, (JP=1).
                CALL XGMG(1,NITB,NBJ(1),nr,DXIXN,GXL,GXU,RGX_local,XW,
     '            ERROR,*9999)
C               Interpolate dep var.s ZW and derivs wrt Xi (JP=0)
                CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXN,
     '            ZE,ZW,XI,ERROR,*9999)
C               Calculate deformed metric tensors wrt Xi (GZL,GZU)
                CALL ZGMG(NBH(NH_LOC(1,nx)),nr,GZ,GZL,GZU,ZW,ERROR,
     '            *9999)
C               Calc deformed 2D Jacobian for the face integration
                RWG=DSQRT(GZ*GZU(3,3))*WG(ng,NBH(NH_LOC(1,nx)))
                IF(JTYP4.EQ.2) RWG=RWG*2.0d0*PI*ZW(1,1) !cyl sym in x
                IF(JTYP4.EQ.3) RWG=RWG*2.0d0*PI*ZW(2,1) !cyl sym in y
                IF(JTYP4.EQ.4) RWG=RWG*4.0d0*PI*ZW(1,1)**2 !sph sym
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' Integrating pressure resids wrt '
     '              //'deformed area: RWG='',D12.4)') RWG
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF !DOP
C               Get derivs of Xi wrt deformed Nu coords, DXIZN
C               and inverse, DZNXI
                CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '            DXIZN,DZNXI,PG,XE,XW,XI,ZE,ZW,'Fibre',ERROR,*9999)
C               Calc derivs of deformed Nu wrt undeformed Nu (DZNXN)
C               and inverse (DXNZN)
                DO ni=1,NITB
                  DO mi=1,NITB
                    SUM=0.0d0
                    DO k=1,NITB
                      SUM=SUM+DZNXI(ni,k)*DXIXN(k,mi)
                    ENDDO !k
                    DZNXN(ni,mi)=SUM
                  ENDDO !mi
                ENDDO !ni
                CALL INVERT(NITB,DZNXN,DXNZN,DET_F_NU)
C MPN 28Feb97 IF(DOP) THEN
                IF(DOP.AND.IWRIT4(nr,nx).GE.2) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                    WRITE(OP_STRING,'('' DZNXN('',I1,'',nj)   : '','
     '                //'3D12.4)')nhx,
     '                (DZNXN(nhx,nj),nj=1,NJ_LOC(NJL_GEOM,0,nr))
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDDO !nhx
CC$                call mp_unsetlock()
                ENDIF !DOP

C               Interpolate dep var.s ZW and derivs wrt Nu (JP=1)
                CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXN,
     '            ZE,ZW,XI,ERROR,*9999)
C               Calculate deformed metric tensors wrt Nu (AZL,AZU)
                CALL ZGMG(NBH(NH_LOC(1,nx)),nr,AZ,AZL,AZU,ZW,ERROR,
     '            *9999)

C               Get contravariant cpts of 2nd Piola-Kirchhoff stress
C               tensor (TW) wrt undeformed Nu coordinates
                IF(KTYP51(nr).EQ.1) THEN !plane stress
                  CALL ZGTG51(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,
     '              CW,RI1,RI2,RI3,TW,YG(1,ng),ZW,ERROR,*9999)
                ELSE IF(KTYP51(nr).EQ.2) THEN !plane strain
                  CALL ZGTG52(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,
     '              CW,RI1,RI2,RI3,TW,YG(1,ng),ZW,ERROR,*9999)
                ELSE IF(KTYP51(nr).EQ.3) THEN !3D
                  CALL ZGTG53(STRESSTYPE,NBH(NH_LOC(1,nx)),
     '              nr,nx,AXU,AZ,AZL,AZU,
     '              CW,EG,RI1,RI2,RI3,TW,XW,YG(1,ng),ZW,ERROR,*9999)
                ELSE IF(KTYP51(nr).EQ.4) THEN !membrane
                  CALL ZGTG54(NBH(NH_LOC(1,nx)),nr,AXU,AZ,AZL,AZU,
     &              CW,EG,RI1,RI2,RI3,TW,YG(1,ng),ERROR,*9999)
                ELSE IF(KTYP51(nr).EQ.5) THEN !string
                  CALL ZGTG55(nr,AZL,CW,EG,TW,ERROR,*9999)
                ENDIF

                IF(KTYP53(nr).EQ.3) THEN !Active stress cmpt included
C!!!              Pick elem Gauss pt nearest to current face GP.
C!!!              To be strictly correct need to either store FEXT
C!!!              at face basis Gauss pts.
C                NGI3=NGAP(3,NBP) !Gauss pts in Xi3 for press basis
C                ng_elem=ng+(iface-1)*(NGI3-1)*NGI2*NGI1
C
C     OR 15-08-06
C     
C     Changes to determine the proper active stress value depending on
C     the definitions in IPACTI. KTYP59==3 got introduced. It allows the
C     user to specify a particular CellML variable to be added to a
C     specified component of the stress tensor
C     
                  IF (KTYP59(nr).EQ.3) THEN !Active stress component defined
                                ! within CellML
                    CALL EVALASC(ne,ng,%VAL(NQNE_PTR),%VAL(YQS_PTR)
     &                   ,%VAL(RCQS_SPATIAL_PTR),%VAL(ICQS_SPATIAL_PTR)
     &                   ,ASC_ARRAYNAME(nr) ,ASC_CELLVARINDEX(nr)
     &                   ,ACTIVE_STRESS,ERROR,*9999)
                  ELSE
                    ACTIVE_STRESS = YG(1,ng)
                  ENDIF
                  CALL ZGTG5A(NBH(NH_LOC(1,nx)),nr,FEXT(1,ng),DXNZN
     &                 ,DZNXN,DET_F_NU,TW,TWA,ACTIVE_STRESS,ERROR,*9999)
C                  CALL ZGTG5A(NBH(NH_LOC(1,nx)),nr,FEXT(1,ng),
C     '                 DXNZN,DZNXN,DET_F_NU,TW,TWA,YG(1,ng),ERROR,*9999)
                ENDIF           !KTYP53(nr).EQ.3 (active stresses)

C new MPN 17Mar97: correcting pressure bc transformations
C               Get Physical Cauchy stress tensor wrt deformed
C               nu-material coordinates (TCNU) from TW
                DO mz=1,3
                  DO nz=1,3
                    SUM=0.0d0
                    DO mj=1,NITB
                      DO nj=1,NITB
                        SUM=SUM+DZNXN(mz,mj)*TW(mj,nj)*DZNXN(nz,nj)
                      ENDDO !nj
                    ENDDO !mj
                    TCNU(mz,nz)=SUM/DSQRT(RI3)
                  ENDDO !nz
                ENDDO !mz
C MPN 28Feb97 IF(DOP) THEN
                IF(DOP.AND.IWRIT4(nr,nx).GE.2) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' TCNU:'',10X,3D12.4,'
     '              //'/(16X,3D12.4))')
     '              ((TCNU(mz,nz),nz=1,3),mz=1,3)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF !DOP
C               Calc def anatomical fibre vects wrt rc coords at XI
                CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,
     '            nr,nx,DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
     '            PG,XE,XW,XI,ZE,ZW,.TRUE.,ERROR,*9999)
C               Calc deformed fibre ref vectors wrt rc coords at XI
                CALL FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH,0,
     '            NHE,NITB,nr,nx,
     '            DZDNUREF(1,1),DZDNUREF(1,2),DZDNUREF(1,3),
     '            PG,XI,ZE,ZW,ERROR,*9999)
                CALL INVERT(NITB,DZDNUREF,DNUREFDZ,DETERM)
                mz=3 !only need derivs wrt third wall coord
                DO nz=1,3
                  SUM=0.0d0
                  DO mi=1,3
                    SUM=SUM+DNUREFDZ(mz,mi)*DZDNU(mi,nz)
                  ENDDO !mi
                  DNUREFDNU(mz,nz)=SUM
                ENDDO !nz
C               Compute cmpts of physical Cauchy stress wrt deformed
C               fibre reference coords using wall coord derivs
C               wrt nu coords
                TCNUREF33=0.0d0
                DO mi=1,3
                  DO ni=1,3
                    TCNUREF33=TCNUREF33+DNUREFDNU(3,mi)*TCNU(mi,ni)*
     '                DNUREFDNU(3,ni)
                  ENDDO !ni
                ENDDO !mi
C old MPN 17Mar97: incorrect pressure bc transformations
C                  DO mz=1,NITB
C                    DO nz=1,NITB
C                      SUM=0.0d0
C                      DO mix=1,NITB
C                        DO nix=1,NITB
C                          SUM=SUM+DZDX(mz,mix)*TW(mix,nix)*DZDX(nz,nix)
C                        ENDDO !nix
C                      ENDDO !mix
C                      TC(mz,nz)=SUM/DSQRT(RI3)
C                    ENDDO !nz
C                  ENDDO !mz
C                  IF(DOP) THEN
CC$                  call mp_setlock()
C                    WRITE(OP_STRING,'('' TC:'',12X,3D12.4,/'//
C     '                '(16X,3D12.4))')
C     '                ((TC(mz,nz),nz=1,3),mz=1,3)
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
C                  ENDIF !DOP
C                  IF(JTYP9.GE.2) THEN !imbric (+ sheet) angles defined
CC                   Calc def anatomical fibre vects wrt rc coords at XI
C                    CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,
C     '                nr,nx,DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
C     '                PG,XE,XW,XI,ZE,ZW,ERROR,*9999)
C                    CALL INVERT(NITB,DZDNU,DNUDZ,DETERM)
CC                   Compute cmpts of Cauchy stress tensor wrt rc coords
CC                   by rotating def material coord system into rc coords
C                    DO mi=1,3
C                      DO ni=1,3
C                        SUM=0.0d0
C                        DO mj=1,3
C                          DO nj=1,3
C                            SUM=SUM+DZDNU(mi,mj)*TC(mj,nj)*DZDNU(ni,nj)
C                          ENDDO !nj
C                        ENDDO !mj
C                        TCRC(mi,ni)=SUM
C                      ENDDO !ni
C                    ENDDO !mi
C                    IF(DOP) THEN
CC$                    call mp_setlock()
C                      WRITE(OP_STRING,'('' TCRC:'',12X,3D12.4,'
C     '                  //'/(18X,3D12.4))')
C     '                  ((TCRC(mi,ni),ni=1,3),mi=1,3)
C                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
C                    ENDIF !DOP
CC                   Calc deformed fibre ref vectors wrt rc coords at XI
C                    CALL FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH,0,
C     '                NHE,NITB,nr,nx,
C     '                DZDNUREF(1,1),DZDNUREF(1,2),DZDNUREF(1,3),
C     '                PG,XI,ZE,ZW,ERROR,*9999)
C                    CALL INVERT(NITB,DZDNUREF,DNUREFDZ,DETERM)
CC                   Compute cmpts of Cauchy stress tensor wrt deformed
CC                   fibre reference coords by rotating rc coord system
CC                   into deformed fibre ref coords
CC                   NOTE: only need TCNUREF(3,3)
C                    mi=3
C                    ni=3
C                    SUM=0.0d0
C                    DO mj=1,3
C                      DO nj=1,3
C                        SUM=SUM+DNUREFDZ(mi,mj)*TCRC(mj,nj)*
C     '                    DNUREFDZ(ni,nj)
C                      ENDDO !nj
C                    ENDDO !mj
C                    TCNUREF(mi,ni)=SUM
Cc                    IF(DOP) THEN
CcC$                    call mp_setlock()
Cc                      WRITE(OP_STRING,'('' TCNUREF:'',12X,3D12.4,'
Cc    '                   //'/(21X,3D12.4))')
Cc    '                   ((TCNUREF(mi,ni),ni=1,3),mi=1,3)
Cc                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CcC$                    call mp_unsetlock()
Cc                    ENDIF !DOP
C                  ELSE !at most fibres defined
C                    TCNUREF33 is aligned with TC(3,3) for
C                    zero imbric/sheet angles
C                    TCNUREF(3,3)=TC(3,3)
C                  ENDIF !JTYP9.GE.2
C                  TCNUREF33=TCNUREF(3,3)
C end old
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' TCNUREF33='',D12.4)') TCNUREF33
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF !DOP

C               Loop over face nodes and derivs for pressure basis
                nsf=0
C               NOTE: nnbf must loop over nn's for face basis fn of third
C               geometric variable since pressure basis faces aren't
C               stored in NPF.
                DO nnbf=1,NNT(NBFF)
                  nn=NNF(1+nnbf,IFE,NBP)
                  DO nkbf=1,NKT(nnbf,NBP)
                    nk=NKEF(nkbf,nnbf,IFE,NBP)
                    nsf=nsf+1
                    nse=nk+(nn-1)*NKT(nn,NBP)
C MPN 28Feb97     IF(DOP) THEN
                    IF(DOP.AND.IWRIT4(nr,nx).GE.2) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'(/'' nn='',I2,'' nk='',I2,'
     '                  //''' nsf='',I2,'' nse='',I2)') nn,nk,nsf,nse
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF !DOP

                    RE(nse)=RE(nse)+
     '                (TCNUREF33+PF(iface))*PG(nsf,1,ng,NBP)*RWG

C MPN 28Feb97     IF(DOP) THEN
                    IF(DOP.AND.IWRIT4(nr,nx).GE.2) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'('' wg='',D12.4,'' PG(nsf..)='','
     '                  //'D12.4,'' RE('',I2,'')='',D12.4)')
     '                  WG(ng,NBH(NH_LOC(1,nx))),PG(nsf,1,ng,NBP),
     '                  nse,RE(nse)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF !DOP
                  ENDDO !nkbf (nk)
                ENDDO !nnbf (nn)
              ENDDO !ng1
            ENDDO !ng2
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(/'' Pressure residuals for face'','
     '          //'I2,'':'')') iface
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' RE(nse): '',8D11.3,'
     '          //'/(10X,8D11.3))') (RE(nse),nse=1,NST(NBP)+NAT(NBP))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF !DOP
          ENDIF !nf.NE.0
        ENDIF !ext press bc
      ENDDO !iface


      CALL EXITS('PFRE_NP')
      RETURN
 9999 CALL ERRORS('PFRE_NP',ERROR)
      CALL EXITS('PFRE_NP')
      RETURN 1
      END


