      SUBROUTINE D_ZERE50(PARAMTYPE,IBT,IDO,INP,NAN,NBH,NBJ,NBJF,
     '  NGAP,ne,NFF,NHE,NKEF,NMNO,NNF,NPNE,nr,NRE,NW,nx,NXI,
     '  CE,CG,CP,D_RE,D_RI3,D_TG,D_ZG,ES,FEXT,PG,RGX,SE,WG,XE,XG,
     '  YG,ZE,D_ZE,ZG,ZG1,ERROR,*)

C#### Subroutine: D_ZERE50
C###  Description:
C###    D_ZERE50 calculates derivatives of element residual D_RE from
C###    current dependent variable array ZE.

C**** The following code was copied from ZERE50 and altered and should
C**** be kept in synch with ZERE50.

      IMPLICIT NONE
      INCLUDE 'acti00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NBJF(NJM,NFM),
     '  NGAP(NIM,NBM),ne,NFF(6),NHE,NKEF(0:4,16,6,NBFM),
     '  NMNO(1:2,0:NOPM),NNF(0:17,6,NBFM),
     '  NPNE(NNM,NBFM,NEM),nr,NRE(NEM),NW,nx,NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),D_RE(NSM,NHM,NOPM),
     '  D_RI3(NHM*NSM),
     '  D_TG(3,3,NHM*NSM),D_ZG(NHM,NUM,NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),FEXT(NIFEXTM,NGM),
     '  PG(NSM,NUM,NGM,NBM),RGX(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM),YG(NIYGM,NGM),
     '  ZE(NSM,NHM),D_ZE(NSM,NHM),ZG(NHM,NUM),ZG1(NHM,NUM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER adj_dirn,i,iface,IFE,IXF,j,k,mi,na,nb,NB1,NBE,
     '  NBFF,NBP,neadj,nf,ng,nh,nh1,nhs,nhs1,NHSGMAX,NHSTART,nhx,nhx1,
     '  ni,NITB,nj,njj,njj1,njj2,nk,nkbf,nn,nnbf,noopti,
     '  ns,ns1,nsa,nse,nsf,NSP(-2:2),NU1(0:3)
      REAL*8 AXU(3,3),
     '  AZ,AZ1,AZL(3,3),AZL1(3,3),AZU(3,3),AZU1(3,3),
     '  CHTOFF(3,3,3),CHTOFF1(3,3,3),CLZ,CMZ,CSLZ,CSMZ,
     '  D11,D12,D13,D21,D22,D23,D31,D32,D33,
     '  D_AG,D_AGE,D_AGE1,D_AGE2,D_AGE3,D_AZ,
     '  D_AZL(3,3),D_AZU(3,3),DBM(3,3,3),
     '  D_CHTOFF(3,3,3),D_delsqP,
     '  D_EG(3,3),DELTA_ZE,DET_F_NU,D_Incomp_resid,
     '  DLA1,DLA2,DLA3,DMA1,DMA2,DMA3,D_RF(32,6),
     '  D_SUM,D_SUM_MEMB(3,40),
     '  DTA1,DTA2,DTA3,DXIX(3,3),
     '  DXIZN(3,3),DXNZN(3,3),
     '  DZNXI(3,3),DZNXN(3,3),
     '  E1,E2,EG(3,3),G1,GXL(3,3),GXU(3,3),
     '  PF(2),PGA1,PGA2,PGA3,PGG,PGX,PPG(64,4),PPGG(4),
     '  RI1,RI2,RI3,RWG,
     '  SLZ,SMZ,SUM,SUM1,TG(3,3),TNA,XI(3),ACTIVE_STRESS,
     &     DUMMY_AZL(3,3),DUMMY_AZU(3,3)
      CHARACTER CHAR1*3,CHAR2*3,STRESSTYPE*17,TYPE*9
      LOGICAL ADJ_XI3_ELEM,ELEMPRESS,NODEPRESS,SAMEDEPBASIS
      DATA DELTA_ZE/1.0D-8/
      DATA NU1/1,2,4,7/
C     OR 23-08-06 : Initialize arrays
      DATA DUMMY_AZL/9 * 0.0d1/
      DATA DUMMY_AZU/9 * 0.0d1/

      CALL ENTERS('D_ZERE50',*9999)
      NITB=NIT(NBJ(1))

C *** Test whether to use unrolled loops
      SAMEDEPBASIS=.TRUE.
      nh1=NH_LOC(1,nx)
      DO nhx=2,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        IF(NBH(nh).NE.NBH(nh1)) SAMEDEPBASIS=.FALSE.
      ENDDO !nhx (nh)
      IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN
        IF(ITYP10(nr).EQ.1.AND.KTYP51(nr).EQ.3.AND.KTYP53(nr).LE.3.AND.
     '    NJ_LOC(NJL_GEOM,0,nr).EQ.3.AND.NIT(NBH(NH_LOC(1,nx))).EQ.3)
     '    THEN
           TYPE='RC3D'
        ELSE IF(ITYP10(nr).EQ.4.AND.KTYP51(nr).EQ.3.AND.
     '      KTYP53(nr).LE.3.AND.SAMEDEPBASIS.AND..NOT.DOP.AND.
     '      NJ_LOC(NJL_GEOM,0,nr).EQ.3.AND.
     '      NIT(NBH(NH_LOC(1,nx))).EQ.3) THEN
          TYPE='PROLATE'
        ELSE
          TYPE=' '
        ENDIF
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          DO ns=1,NST(NBH(nh))+NAT(NBH(nh))
            DO noopti=1,NTOPTI
              D_RE(ns,nh,noopti)=0.0d0
            ENDDO !noopti
          ENDDO !ns
        ENDDO !nhx (nh)
        CALL ASSERT(NTOPTI.LE.NHM*NSM,'>>Increase dimension of D_RE',
     '    ERROR,*9999)

      ELSE IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN
        TYPE=' '
        nhs=0
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          DO ns=1,NST(NBH(nh))+NAT(NBH(nh))
            D_ZE(ns,nhx)=0.0d0
            nhs=nhs+1
            nhs1=0
            DO nhx1=1,NH_LOC(0,nx)
              nh1=NH_LOC(nhx1,nx)
              DO ns1=1,NST(NBH(nh1))+NAT(NBH(nh1))
                nhs1=nhs1+1
                ES(nhs,nhs1)=0.0d0
              ENDDO !ns1
            ENDDO !nhx1
          ENDDO !ns
        ENDDO !nhx
        IF(KTYP51(nr).EQ.4.AND.NW.GT.1) THEN !membrane + press bc(s)
          WRITE(CHAR1,'(I3)') NJ_LOC(NJL_GEOM,0,nr)
          WRITE(CHAR2,'(I3)') nhs1
          CALL ASSERT(NJ_LOC(NJL_GEOM,0,nr).LE.3.AND.nhs1.LE.40,
     '      '>>Increase dimension of D_SUM_MEMB to be ('
     '      //CHAR1(1:3)//','//CHAR2(1:3)//')',ERROR,*9999)
        ENDIF!membrane + press bc(s)
      ENDIF !PARAMTYPE

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
     '        ERROR,*9999)
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

      CALL ASSERT(KTYP51(nr).EQ.3.OR.KTYP51(nr).EQ.4,
     '  '>>Analytic derivs only implemented for '
     '  //'3D and membrane problems',ERROR,*9999)
      CALL ASSERT(KTYP53(nr).GT.1,
     '  '>>Analytic derivs only implemented for '
     '  //'stresses ref.d to nu in constit. law',ERROR,*9999)

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
C       Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C       derivs (DXIX) of Xi wrt Nu coords.
        CALL XGMG(1,NIT(NBJ(1)),NBJ(1),nr,DXIX,GXL,GXU,
     '    RGX(ng),XG,ERROR,*9999)
C       Calculate the Jacobian for integration wrt undef coords:
        RWG=RGX(ng)*WG(ng,NBH(NH_LOC(1,nx)))
        IF(JTYP4.EQ.2) RWG=RWG*2.0d0*PI*XG(1,1)    !cyl symm about x
        IF(JTYP4.EQ.3) RWG=RWG*2.0d0*PI*XG(2,1)    !cyl symm about y
        IF(JTYP4.EQ.4) RWG=RWG*4.0d0*PI*XG(1,1)**2 !spherical symm

        IF(KTYP53(nr).EQ.3) THEN !Active stress component included
C old MPN 18Mar97: DXIXN is already calc'ed in XGMG above (DXIX)
CC         Get derivs of Xi wrt undeformed Nu (body/fibre) coords,DXIXN
C          CALL DXIDXM(NBJ(1),nr,DXIXN,DXNXI,GXL,GXU,XG,'Fibre',ERROR,*9999)
C         Get derivs of Xi wrt deformed Nu coords, DXIZN
          CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '      DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,'Fibre',ERROR,*9999)
C         Calculate derivs of deformed Nu wrt undeformed Nu (DZNXN)
          DO ni=1,NITB
            DO mi=1,NITB
              SUM1=0.0d0
              DO k=1,NITB
                SUM1=SUM1+DZNXI(ni,k)*DXIX(k,mi)
              ENDDO !k
              DZNXN(ni,mi)=SUM1
            ENDDO !mi
          ENDDO !ni
          CALL INVERT(NITB,DZNXN,DXNZN,DET_F_NU)
        ENDIF !KTYP53(nr).EQ.3 (active stresses)

C       Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
        CALL ZEZG(1,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
C       Calculate deformed metric tensors wrt Xj or Nu (AZL,AZU)
        CALL ZGMG(NBH(NH_LOC(1,nx)),nr,AZ,AZL,AZU,ZG,ERROR,*9999)

        IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN
C         Calc deriv's of contravariant cpts of 2nd Piola-Kirchhoff
C         stress tensor (D_TW - undeformed Nu coordinates) wrt
C         each of the material parameters
C         Note: AZ,AZL,AZU,EG are all independent of material params
C         so the calculation of D_AZ,D_AZL,D_AZU is not
C         needed here.
          IF(KTYP51(nr).EQ.3) THEN
            DO noopti=1,NTOPTI
C     OR 22-08-06 Apparently the values for arguments D_AZL and D_AZU
C     are not important, but %VAL(0) causes a warning, hence I replaced
C     %VAL(0) with DUMMY_AZL and DUMMY_AZU respectively
C             CALL D_ZGTG53(PARAMTYPE,NBH(NH_LOC(1,nx)),
C     '          NMNO(1,noopti),nr,nx,
C     '          AXU,AZ,AZL,AZU,CG(1,ng),0.d0,%VAL(0),%VAL(0),D_EG,
C     '          D_RI3(noopti),D_TG(1,1,noopti),D_ZG(1,1,noopti),EG,
C     '          ZG,ERROR,*9999)
              CALL D_ZGTG53(PARAMTYPE,NBH(NH_LOC(1,nx)), NMNO(1,noopti)
     &             ,nr,nx, AXU,AZ,AZL,AZU,CG(1,ng),0.d0,DUMMY_AZL
     &             ,DUMMY_AZU,D_EG, D_RI3(noopti),D_TG(1,1,noopti)
     &             ,D_ZG(1,1,noopti),EG, ZG,ERROR,*9999)
            ENDDO !noopti
          ELSE IF(KTYP51(nr).EQ.4) THEN
            DO noopti=1,NTOPTI
              CALL D_ZGTG54(PARAMTYPE,NMNO(1,noopti),nr,
     '          AXU,AZ,AZL,AZU,CG(1,ng),D_AZ,D_AZL,D_AZU,D_EG,
     '          D_RI3(noopti),D_TG(1,1,noopti),EG,YG(1,ng),ERROR,*9999)
            ENDDO !noopti
          ENDIF !KTYP51(nr)

        ELSE IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN
C         Calc contravariant cpts of 2nd Piola-Kirchhoff stress
C         tensor (TG) wrt undeformed Nu coordinates
          IF(KTYP51(nr).EQ.3) THEN
            CALL ZGTG53(STRESSTYPE,NBH(NH_LOC(1,nx)),
     '        nr,nx,AXU,AZ,AZL,AZU,
     '        CG(1,ng),EG,RI1,RI2,RI3,TG,XG,YG(1,ng),ZG,ERROR,*9999)
          ELSE IF(KTYP51(nr).EQ.4) THEN
            CALL ZGTG54(NBH(NH_LOC(1,nx)),nr,AXU,AZ,AZL,AZU,
     '        CG(1,ng),EG,RI1,RI2,RI3,TG,YG,ERROR,*9999)
          ENDIF !KTYP51(nr)

C         Calc deriv's of contravariant cpts of 2nd Piola-Kirchhoff
C         stress tensor (D_TG - undeformed Nu coordinates) and deriv
C         of third strain invariant (D_RI3) wrt each of the geometric
C         variables
          nhs1=0
          DO nhx1=1,NH_LOC(0,nx)
            nh1=NH_LOC(nhx1,nx)
            NB1=NBH(nh1)
            DO ns1=1,NST(NB1)+NAT(NB1)
              nhs1=nhs1+1
C!!!          The deriv of the deformed Jacobian D_RWG could be
C!!!          calc'ed here and stored in an array D_RWG(nhs1)
C!!!          (same size as D_RI3), that needs to be set up and passed
C!!!          to this routine (instead of being declared above).
C!!!          Instead a finite difference approx is computed below.
C             Calc D_ZG, deriv of def coords (ZG + Nu derivs)
C             wrt elem coords (ZE)
C             NOTE: code up to the next ! is similar to ZEZG with JP=1
              DO ni=0,NIT(NB1)
                IF(NI.EQ.0) THEN
                  D_ZG(nhx1,NU1(ni),nhs1)=PG(ns1,NU1(ni),ng,NB1)
                ELSE  !derivs of basis fns wrt Nu coords
                  D_ZG(nhx1,NU1(ni),nhs1)=PG(ns1,2,ng,NB1)*DXIX(1,ni) +
     '                                   PG(ns1,4,ng,NB1)*DXIX(2,ni) +
     '                                   PG(ns1,7,ng,NB1)*DXIX(3,ni)
                ENDIF !NI.EQ.0
              ENDDO !ni
C             !end of ZEZG stuff
C             Calculate derivs of deformed metric tensors, D_AZL and
C             D_AZU (Nu coords) and deriv of the determinant of AZL,
C             D_AZ wrt the current geom. var. using finite
C             difference approximations.
C             Perturb deformed element coords
              ZE(ns1,nhx1)=ZE(ns1,nhx1)+DELTA_ZE
C             Interpolate perturbed dep. var.s ZG1 and derivs wrt Nu
              CALL ZEZG(1,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG1,ERROR,*9999)
C             Reset deformed element coords
              ZE(ns1,nhx1)=ZE(ns1,nhx1)-DELTA_ZE
C             Calculate deformed metric tensors for perturbed coords
              CALL ZGMG(NBH(NH_LOC(1,nx)),nr,AZ1,AZL1,AZU1,ZG1,
     '          ERROR,*9999)
C             Finite diff approx to derivs
              DO I=1,3
                DO J=1,3
                  D_AZL(i,j)=(AZL1(i,j)-AZL(i,j))/DELTA_ZE
                  D_AZU(i,j)=(AZU1(i,j)-AZU(i,j))/DELTA_ZE
                ENDDO !J
              ENDDO !I
              D_AZ=(AZ1-AZ)/DELTA_ZE
C             Calc derivs of TG wrt current geom var analytically
              IF(KTYP51(nr).EQ.3) THEN
                CALL D_ZGTG53(PARAMTYPE,NBH(NH_LOC(1,nx)),0,nr,nx,
     '            AXU,AZ,AZL,AZU,CG(1,ng),D_AZ,D_AZL,D_AZU,D_EG,
     '            D_RI3(nhs1),D_TG(1,1,nhs1),D_ZG(1,1,nhs1),EG,ZG,
     '            ERROR,*9999)
              ELSE IF(KTYP51(nr).EQ.4) THEN
                CALL D_ZGTG54(PARAMTYPE,0,nr,
     '            AXU,AZ,AZL,AZU,CG(1,ng),D_AZ,D_AZL,D_AZU,D_EG,
     '            D_RI3(nhs1),D_TG(1,1,nhs1),EG,YG(1,ng),ERROR,*9999)
                IF(NW.GT.1) THEN !membrane+press bc(s)
C                 calculate this quantity here for use later in Gauss pt
C                 loop (saves recalculation of D_AZU later)
                  DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                    nh=NH_LOC(nhx,nx)
                    nb=NBH(nh)
                    D_SUM_MEMB(nh,nhs1)=0.0d0
                    DO ni=1,NIT(nb)
                      D_SUM_MEMB(nhx,nhs1)=D_SUM_MEMB(nhx,nhs1)
     '                  +D_AZU(ni,3)*ZG(nhx,NU1(ni))
     '                  +AZU(ni,3)*D_ZG(nhx,NU1(ni),nhs1)
                    ENDDO !ni
                  ENDDO !nhx
                ENDIF !membrane+press bc(s)
              ENDIF !KTYP51(nr)
            ENDDO !ns1
          ENDDO !nhx1
        ENDIF !PARAMTYPE

C ***   Active stress component 
        IF(KTYP53(nr).EQ.3) THEN 
C     09-Dec-1989: NOTE: Don't have to define a separate face array
C     CALL ZGTG5A(NBH(NH_LOC(1,nx)),nr,FEXT(1,ng),DXNZN,DZNXN,
C     '      DET_F_NU,TG,TNA,YG(1,ng),ERROR,*9999)
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
          CALL ZGTG5A(NBH(NH_LOC(1,nx)),nr,FEXT(1,ng),DXNZN,DZNXN,
     &         DET_F_NU,TG,TNA,ACTIVE_STRESS,ERROR,*9999)          
        ENDIF !KTYP53(nr)

C ***   Derivative of main element residual
        IF(TYPE(1:4).EQ.'RC3D') THEN     !3D rect.cart.
C         Calculate derivs wrt material params using unrolled loops
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' >>Using unrolled loops'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh)
            DO ns=1,NST(nb)+NAT(nb)
              PPGG(2) = PG(ns,NU1(1),ng,nb)*DXIX(1,1) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,1) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,1)
              PPGG(3) = PG(ns,NU1(1),ng,nb)*DXIX(1,2) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,2) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,2)
              PPGG(4) = PG(ns,NU1(1),ng,nb)*DXIX(1,3) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,3) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,3)
              DO noopti=1,NTOPTI
                D_AGE=(D_TG(1,1,noopti)*ZG(nhx,2)+
     '                 D_TG(1,2,noopti)*ZG(nhx,4)+
     '                 D_TG(1,3,noopti)*ZG(nhx,7))*PPGG(2)
     '               +(D_TG(2,1,noopti)*ZG(nhx,2)+
     '                 D_TG(2,2,noopti)*ZG(nhx,4)+
     '                 D_TG(2,3,noopti)*ZG(nhx,7))*PPGG(3)
     '               +(D_TG(3,1,noopti)*ZG(nhx,2)+
     '                 D_TG(3,2,noopti)*ZG(nhx,4)+
     '                 D_TG(3,3,noopti)*ZG(nhx,7))*PPGG(4)
                D_RE(ns,nh,noopti)=D_RE(ns,nh,noopti)+
     '            D_AGE*RWG*SE(ns,nb,ne)
              ENDDO !noopti
            ENDDO !ns
          ENDDO !nhx (nh)

        ELSE IF(TYPE(1:7).EQ.'PROLATE') THEN  !3D prolate spheroidal
C         Calculate derivs wrt material params using unrolled loops
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
          CLZ=SQRT(1.0d0+SLZ*SLZ)
          CMZ=SQRT(1.0d0-SMZ*SMZ)
          CSLZ=CLZ/SLZ
          CSMZ=CMZ/SMZ
          G1=SLZ*SLZ+SMZ*SMZ
          E1=CLZ*SLZ/G1
          E2=CMZ*SMZ/G1

          nb=NBH(NH_LOC(1,nx))
          DO ns=1,NST(nb)+NAT(nb)
            PPG(ns,1) = PG(ns,1,ng,NBH(NH_LOC(1,nx)))
            PPG(ns,2) = PG(ns,NU1(1),ng,NBH(NH_LOC(1,nx)))*DXIX(1,1) +
     '                  PG(ns,NU1(2),ng,NBH(NH_LOC(1,nx)))*DXIX(2,1) +
     '                  PG(ns,NU1(3),ng,NBH(NH_LOC(1,nx)))*DXIX(3,1)
            PPG(ns,3) = PG(ns,NU1(1),ng,NBH(NH_LOC(1,nx)))*DXIX(1,2) +
     '                  PG(ns,NU1(2),ng,NBH(NH_LOC(1,nx)))*DXIX(2,2) +
     '                  PG(ns,NU1(3),ng,NBH(NH_LOC(1,nx)))*DXIX(3,2)
            PPG(ns,4) = PG(ns,NU1(1),ng,NBH(NH_LOC(1,nx)))*DXIX(1,3) +
     '                  PG(ns,NU1(2),ng,NBH(NH_LOC(1,nx)))*DXIX(2,3) +
     '                  PG(ns,NU1(3),ng,NBH(NH_LOC(1,nx)))*DXIX(3,3)
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

          DO ns=1,NST(nb)+NAT(nb)
            PGG =PPG(ns,1)
            PGA1=PPG(ns,2)
            PGA2=PPG(ns,3)
            PGA3=PPG(ns,4)

            DO noopti=1,NTOPTI
              D_AGE1=D_TG(1,1,noopti)*(D11*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+
     '          D21*(E1*DMA1-E2*DLA1)*PGG+D31*E1*SMZ*SMZ*DTA1*PGG) +
     '          D_TG(2,1,noopti)*(D11*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '          D21*(E1*DMA2-E2*DLA2)*PGG+D31*E1*SMZ*SMZ*DTA2*PGG) +
     '          D_TG(3,1,noopti)*(D11*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '          D21*(E1*DMA3-E2*DLA3)*PGG+D31*E1*SMZ*SMZ*DTA3*PGG) +
     '          D_TG(1,2,noopti)*(D12*(PGA1-(E1*DLA1+E2*DMA1)*PGG) +
     '          D22*(E1*DMA1-E2*DLA1)*PGG+D32*E1*SMZ*SMZ*DTA1*PGG) +
     '          D_TG(2,2,noopti)*(D12*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '          D22*(E1*DMA2-E2*DLA2)*PGG+D32*E1*SMZ*SMZ*DTA2*PGG) +
     '          D_TG(3,2,noopti)*(D12*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '          D22*(E1*DMA3-E2*DLA3)*PGG+D32*E1*SMZ*SMZ*DTA3*PGG) +
     '          D_TG(1,3,noopti)*(D13*(PGA1-(E1*DLA1+E2*DMA1)*PGG) +
     '          D23*(E1*DMA1-E2*DLA1)*PGG+D33*E1*SMZ*SMZ*DTA1*PGG) +
     '          D_TG(2,3,noopti)*(D13*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '          D23*(E1*DMA2-E2*DLA2)*PGG+D33*E1*SMZ*SMZ*DTA2*PGG) +
     '          D_TG(3,3,noopti)*(D13*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '          D23*(E1*DMA3-E2*DLA3)*PGG+D33*E1*SMZ*SMZ*DTA3*PGG)

              D_AGE2 = D_TG(1,1,noopti)*(D11*(E2*DLA1-E1*DMA1)*PGG +
     '          D21*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D31*E2*SLZ*SLZ*DTA1*
     '            PGG)+
     '          D_TG(2,1,noopti)*(D11*(E2*DLA2-E1*DMA2)*PGG +
     '          D21*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D31*E2*SLZ*SLZ*DTA2*
     '            PGG)+
     '          D_TG(3,1,noopti)*(D11*(E2*DLA3-E1*DMA3)*PGG +
     '          D21*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D31*E2*SLZ*SLZ*DTA3*
     '            PGG)+
     '          D_TG(1,2,noopti)*(D12*(E2*DLA1-E1*DMA1)*PGG +
     '          D22*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D32*E2*SLZ*SLZ*DTA1*
     '            PGG)+
     '          D_TG(2,2,noopti)*(D12*(E2*DLA2-E1*DMA2)*PGG +
     '          D22*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D32*E2*SLZ*SLZ*DTA2*
     '            PGG)+
     '          D_TG(3,2,noopti)*(D12*(E2*DLA3-E1*DMA3)*PGG +
     '          D22*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D32*E2*SLZ*SLZ*DTA3*
     '            PGG)+
     '          D_TG(1,3,noopti)*(D13*(E2*DLA1-E1*DMA1)*PGG +
     '          D23*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D33*E2*SLZ*SLZ*DTA1*
     '            PGG)+
     '          D_TG(2,3,noopti)*(D13*(E2*DLA2-E1*DMA2)*PGG +
     '          D23*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D33*E2*SLZ*SLZ*DTA2*
     '            PGG)+
     '          D_TG(3,3,noopti)*(D13*(E2*DLA3-E1*DMA3)*PGG +
     '          D23*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D33*E2*SLZ*SLZ*DTA3*
     '            PGG)

              D_AGE3 = D_TG(1,1,noopti)*(-(D11*CSLZ+D21*CSMZ)*DTA1*PGG +
     '          D31*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '          D_TG(2,1,noopti)*(-(D11*CSLZ+D21*CSMZ)*DTA2*PGG +
     '          D31*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '          D_TG(3,1,noopti)*(-(D11*CSLZ+D21*CSMZ)*DTA3*PGG +
     '          D31*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG)) +
     '          D_TG(1,2,noopti)*(-(D12*CSLZ+D22*CSMZ)*DTA1*PGG +
     '          D32*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '          D_TG(2,2,noopti)*(-(D12*CSLZ+D22*CSMZ)*DTA2*PGG +
     '          D32*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '          D_TG(3,2,noopti)*(-(D12*CSLZ+D22*CSMZ)*DTA3*PGG +
     '          D32*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG)) +
     '          D_TG(1,3,noopti)*(-(D13*CSLZ+D23*CSMZ)*DTA1*PGG +
     '          D33*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '          D_TG(2,3,noopti)*(-(D13*CSLZ+D23*CSMZ)*DTA2*PGG +
     '          D33*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '          D_TG(3,3,noopti)*(-(D13*CSLZ+D23*CSMZ)*DTA3*PGG +
     '          D33*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG))

              D_RE(ns,1,noopti)=D_RE(ns,1,noopti)
     '          +D_AGE1*RWG*SE(ns,NBH(NH_LOC(1,nx)),ne)
              D_RE(ns,2,noopti)=D_RE(ns,2,noopti)
     '          +D_AGE2*RWG*SE(ns,NBH(NH_LOC(2,nx)),ne)
              D_RE(ns,3,noopti)=D_RE(ns,3,noopti)
     '          +D_AGE3*RWG*SE(ns,NBH(NH_LOC(3,nx)),ne)

            ENDDO !noopti
          ENDDO !ns

        ELSE
C         Calculate derivs wrt mat or geometric params for general case
          nhs=0
          DO njj1=1,NJ_LOC(NJL_GEOM,0,nr)
            nh=NJ_LOC(NJL_GEOM,njj1,nr)
            nb=NBH(nh)
            DO ns=1,NST(nb)+NAT(nb)
              PPG(ns,1)=PG(ns,1,ng,nb)
              DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj2,nr)
                PPG(ns,1+nj)=PGX(nb,nj,ns,DXIX,PG(1,1,ng,nb))
              ENDDO !njj2
            ENDDO !ns

            IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN
              DO ns=1,NST(nb)+NAT(nb)
                DO noopti=1,NTOPTI
                  D_AGE=D_AG(PARAMTYPE,nb,nhx,nr,ns,D_TG(1,1,noopti),
     '              D_ZG(1,1,noopti),PPG,TG,ZG)
                  D_RE(ns,nh,noopti)=D_RE(ns,nh,noopti)+
     '              D_AGE*RWG*SE(ns,nb,ne)
                ENDDO !noopti
              ENDDO !ns

            ELSE IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN
              DO ns=1,NST(nb)+NAT(nb)
                nhs=nhs+1
                nhs1=0
                DO nhx1=1,NH_LOC(0,nx)
                  nh1=NH_LOC(nhx1,nx)
                  NB1=NBH(nh1)
                  DO ns1=1,NST(NB1)+NAT(NB1)
                    nhs1=nhs1+1
                    D_AGE=D_AG(PARAMTYPE,nb,nhx,nr,ns,D_TG(1,1,nhs1),
     '                D_ZG(1,1,nhs1),PPG,TG,ZG)
                    ES(nhs,nhs1)=ES(nhs,nhs1)+D_AGE*RWG*SE(ns,nb,ne)
                  ENDDO !ns1
                ENDDO !nhx1
              ENDDO !ns
            ENDIF !PARAMTYPE
          ENDDO !njj1
          NHSGMAX=nhs   !max row # for geom resids
        ENDIF !TYPE

C ***   Incompressibilty (+fluid or +inext) constraints
        IF(KTYP51(nr).EQ.3.AND.KTYP52(nr).EQ.2       !3D & incompressible
     '                     .OR.KTYP52(nr).EQ.3       !3D & incomp + fluid
     '                     .OR.KTYP52(nr).EQ.5) THEN !3D & incomp + inext
          NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for press vars
          IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN
            DO ns=1,NST(NBP)+NAT(NBP)
              DO noopti=1,NTOPTI
                D_RE(ns,NH_LOC(NH_LOC(0,nx),nx),noopti)=0.0d0
              ENDDO !noopti
            ENDDO !ns
          ELSE IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN
            IF(KTYP52(nr).EQ.3) THEN !Incompressible + fluid
C             Calc deformed Christoffel symbols wrt undef Nu coords
C             NOTE: ZG needs derivs wrt Nu not Xi !
C KAT 1Nov00: X3G should not be used in TOFFEL as it is not set up.
              CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF,DBM,AZU,ZG,
     '          %VAL(0),.FALSE.,ERROR,*9999)
C              CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF,DBM,AZU,ZG,X3G,
C     '          .FALSE.,ERROR,*9999)
            ENDIF !KTYP52(nr)=3
            DO ns=1,NST(NBP)+NAT(NBP)
              nhs=nhs+1
              nhs1=0
              DO nhx1=1,NH_LOC(0,nx)
                nh1=NH_LOC(nhx1,nx)
                NB1=NBH(nh1)
                DO ns1=1,NST(NB1)+NAT(NB1)
                  nhs1=nhs1+1
                  IF(KTYP52(nr).EQ.2) THEN !standard incomp constraint
                    D_Incomp_resid=D_RI3(nhs1)/(2.0d0*DSQRT(RI3))
                  ELSE IF(KTYP52(nr).EQ.3) THEN !Incompressible + fluid
C                   Calculate deformed Christoffel symbols wrt undef
C                   Nu coords at perturbed geom coords
C                   NOTE: ZG1 needs derivs wrt Nu not Xi!
C KAT 1Nov00: X3G should not be used in TOFFEL as it is not set up.
                    CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF1,DBM,AZU,
     '                ZG1,%VAL(0),.FALSE.,ERROR,*9999)
C                    CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF1,DBM,AZU,
C     '                ZG1,X3G,.FALSE.,ERROR,*9999)
                    DO I=1,3
                      DO J=1,3
                        D_AZU(i,j)=(AZU1(i,j)-AZU(i,j))/DELTA_ZE
                        DO K=1,3
                          D_CHTOFF(i,j,k)=(CHTOFF1(i,j,k)-CHTOFF(i,j,k))
     '                      /DELTA_ZE
                        ENDDO !K
                      ENDDO !J
                    ENDDO !I
C                   Calc deriv wrt current geom var of del-squared(p)
C                   (wrt def coords), where p is the hyd press
C                   that varies with Xi3 only.
                    SUM=0.0d0
                    D_SUM=0.0d0
                    DO J=1,NITB
                      DO K=1,NITB
                        SUM=SUM+CHTOFF(3,j,k)*AZU(j,k)
                        D_SUM=D_SUM+D_CHTOFF(3,j,k)*AZU(j,k)
     '                             +CHTOFF(3,j,k)*D_AZU(j,k)
                      ENDDO !K
                    ENDDO !J
                    D_delsqP=0.0d0
                    DO nsa=1,NST(NBP)+NAT(NBP)
                      D_delsqP=D_delsqP+ZE(nsa,NH_LOC(0,nx))*
     '                  (PG(nsa,8,ng,NBP)*D_AZU(3,3)-
     '                  D_SUM*PG(nsa,7,ng,NBP))
                    ENDDO !nsa
C                   Check if taking deriv wrt the current ZE aux var
                    IF(nh1.EQ.NH_LOC(NH_LOC(0,nx),nx)) THEN
                      D_delsqP=D_delsqP+
     '                  (PG(ns1,8,ng,NBP)*AZU(3,3)-SUM*PG(ns1,7,ng,NBP))
                    ENDIF !nh1=pressure var
C                   Calc derivs of resids associated with Darcy's Law
C new MPN 18Mar98: fixed error
                    D_Incomp_resid=CG(IL_fluid_conductivity,ng)*DT*
     '                D_delsqP+D_RI3(nhs1)/(2.0d0*DSQRT(RI3))
C old                    D_Incomp_resid=CG(IL_fluid_conductivity,ng)*DT*
C old     '                D_delsqP+D_RI3(nhs1)/(2.0d0*DSQRT(RI3*RI3*RI3))
                  ELSE IF(KTYP52(nr).EQ.5) THEN !Incompressible + inextensible
                  ENDIF !KTYP52(nr)

C                 Derivs of Galerkin incompressibility (+ fluid) resids
                  ES(nhs,nhs1)=ES(nhs,nhs1)+
     '              D_Incomp_resid*RWG*PG(ns,1,ng,NBP)

C old MPN 17Mar97: use standard Galerkin resids (above)
C                  IF(NODEPRESS.AND.ELEMPRESS) THEN
CC                   node and element based hydrostatic pressure interp
CC                   NOTE: pressure basis function is not used to
CC                         weight the residual here (non-Galerkin)
C                    ES(nhs,nhs1)=ES(nhs,nhs1)+D_Incomp_resid*RWG
C                  ELSE IF(NODEPRESS) THEN
CC                   Node based hydrostatic pressure interpolation.
CC                   Weight with Xi3=0 face basis (no Xi3 var'n for wght)
C                    IFE=5
C                    nf=NFF(IFE) !face # for Xi3=0 face in current elem
C                    NBFF=NBJF(3,nf) !bas fn for 3rd geom var on face nf
CC                             !!!WARNING: Need to store press face bases
CC                             !!!in NPF
CC                   Determine Xi3=0 face GP#
C                    ng_face=MOD(ng-1,NGAP(1,NBP)*NGAP(2,NBP))+1
CC                   Determine Xi3 face vertex for element vertex ns
C                    nsf=MOD(ns-1,NST(NBFF))+1
C                    ES(nhs,nhs1)=ES(nhs,nhs1)+D_Incomp_resid*RWG*
C     '                PG(nsf,1,ng_face,NBFF)
C                  ELSE
CC                   element based hyd. pressure interp
CC                   (standard Galerkin)
C                    ES(nhs,nhs1)=ES(nhs,nhs1)+D_Incomp_resid*RWG*
C     '                PG(ns,1,ng,NBP)
C                  ENDIF !NODEPRESS/ELEMPRESS
                ENDDO !ns1
              ENDDO !nhx1
            ENDDO !ns
          ENDIF !PARAMTYPE

        ELSE IF(KTYP51(nr).EQ.4.AND.NW.GT.1) THEN !membrane+press bc(s)
          IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN
            nhs=0
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              nb=NBH(nj)
C             NOTE: deriv of SUM, D_SUM_MEMB(nj,nhs1) was calc'ed above
              DO ns=1,NST(nb)+NAT(nb)
                nhs=nhs+1
                nhs1=0
                DO nhx1=1,NH_LOC(0,nx)
                  nh1=NH_LOC(nhx1,nx)
                  NB1=NBH(nh1)
                  DO ns1=1,NST(NB1)+NAT(NB1)
                    nhs1=nhs1+1
                    ES(nhs,nhs1)=ES(nhs,nhs1)
     '                +PF(1)*D_SUM_MEMB(nj,nhs1)*PG(ns,1,ng,nb)*RWG
                  ENDDO !ns1
                ENDDO !nhx1
              ENDDO !ns
            ENDDO !njj
          ENDIF !PARAMTYPE

        ENDIF !KTYP51(nr)...

 50   CONTINUE !end of ng loop

C *** External pressure loads in 3D case
      IF(KTYP51(nr).EQ.3) THEN !3D
        NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for hydr. pressure
C       Boundary pressure constraint equations
        IF(   KTYP52(nr).EQ.2       !incomp
     '    .OR.KTYP52(nr).EQ.3       !incomp + fluid
     '    .OR.KTYP52(nr).EQ.5       !incomp + inext
     '   .AND.KTYP57(nr).GT.1) THEN !pressure entered
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
              ERROR='>>Derivatives not implemented'
              GOTO 9999
            ENDIF !KTYP5A(nr).EQ.2
          ELSE IF(ELEMPRESS) THEN
C           element based hydrostatic pressure interp
            IF(NSP(1).NE.0.AND.NSP(2).NE.0) THEN
C             Determine hydrostatic pressure vars using face press bcs
C             instead of incompressibility constraints
              CALL D_PFRE_NE(PARAMTYPE,IBT,IDO,INP,
     '          NAN,NBH,NBJ,NFF,NGAP,NHE,NHSGMAX,NMNO,
     '          NPNE(1,1,ne),nr,ne,NRE,NSP,NW,nx,NXI(-NIM,0,ne),
     '          CE,CP,D_RE,D_RI3,D_TG,D_ZE,D_ZG,FEXT,PG,
     '          XE,XG,YG,ZE,ZG,ZG1,ERROR,*9999)
            ENDIF !NSP(1).NE.0.AND.NSP(2).NE.0
          ENDIF !NODEPRESS/ELEMPRESS
        ENDIF !KTYP52(nr)=2,3,5 & KTYP57(nr)>1

C       Contribution of external pressure loads to stress equilibrium
C       equation residuals at nodes.
C       Note that this needs to be done only for the derivatives wrt
C       geometric parameters as these contributions are independent
C       of the material parameters.
        IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN
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
     '        .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext press bc on Xi3=0 face
     '        iface.EQ.2.AND..NOT.ADJ_XI3_ELEM
     '        .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !ext pres bc on Xi3=1 face
C old
C            IF(iface.EQ.1.AND.(NXI(-3,1,ne).EQ.0.OR.nr.NE.
C     '        NRE(NXI(-3,1,ne)))
C     '        .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext press bc on Xi3=0 face
C     '         iface.EQ.2.AND.(NXI( 3,1,ne).EQ.0.OR.nr.NE.
C     '        NRE(NXI(3,1,ne)))
C     '        .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !ext pres bc on Xi3=1 fac

              IF(DABS(PF(iface)).GT.ZERO_TOL) THEN
                IXF=iface-1
                IFE=iface+4
                nf=NFF(IFE)
                IF(nf.NE.0) THEN
                  nhs1=0
                  DO nhx1=1,NH_LOC(0,nx)
                    nh1=NH_LOC(nhx1,nx)
                    NB1=NBH(nh1)
                    DO ns1=1,NST(NB1)+NAT(NB1)
                      nhs1=nhs1+1
C                     Analytic derivs of face press resids wrt current
C                     geom variable
                      CALL D_PFRF(IBT,IDO,INP,IXF,NAN,NBH,NBJ,
     '                  NBJF(1,nf),NGAP,NHE,nhx1,nr,ns1,nx,D_RF,D_ZE,
     '                  D_ZG(1,1,nhs1),PF(iface),PG,WG,XE,XG,ZE,ZG,
     '                  ZG1,ERROR,*9999)
                      NHSTART=0
                      DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                        nh=NH_LOC(nhx,nx)
                        nsf=0
                        NBFF=NBJF(nhx,nf)
                        NBE=NBH(nh)
                        DO nnbf=1,NNT(NBFF)
                          nn=NNF(1+nnbf,IFE,NBE)
                          DO nkbf=1,NKT(nnbf,NBFF)
                            nk=NKEF(nkbf,nnbf,IFE,NBE)
                            nse=nk+(nn-1)*NKT(nn,NBE)
                            nhs=NHSTART+nk+(nn-1)*NKT(nn,NBE)
                            nsf=nsf+1
                            ES(nhs,nhs1)=ES(nhs,nhs1)
     '                        -D_RF(nsf,nh)*SE(nse,NBE,ne)
                          ENDDO !nkbf (nk)
                        ENDDO !nnbf (nn)
                        DO nn=1,NNT(NBE)
                          NHSTART=NHSTART+NKT(nn,NBE)
                        ENDDO !nn
                      ENDDO !nhx
                    ENDDO !ns1
                  ENDDO !nhx1
                ENDIF !nf.NE.0
              ENDIF !PF(iface)>0
            ENDIF !ext press bc on face
          ENDDO !iface
        ENDIF !PARAMTYPE
      ENDIF !KTYP51(nr).EQ.3

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>D_ZERE50 diagnostic output'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN
          WRITE(OP_STRING,'(/('' D_RE('',I2,'','',I1,'
     '      //''',noopti): '',5(1X,D12.4),/(20X,5(1X,D12.4))))')
     '      ((nh,(ns,D_RE(ns,nh,noopti),noopti=1,NTOPTI),
     '      ns=1,NST(NBH(NH_LOC(nhx,nx)))
     '      +NAT(NBH(NH_LOC(nhx,nx)))),nhx=1,NH_LOC(0,nx))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN
          nhs=0
          DO nhx=1,NH_LOC(0,nx)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh)
            DO ns=1,NST(nb)+NAT(nb)
              nhs=nhs+1
              WRITE(OP_STRING,'(/'' ES('',I3,'
     '          //''',nhs1): '',5(1X,D12.4),/(15X,5(1X,D12.4)))')
     '          nhs,((ES(nhs,ns1+(nh1-1)*(NST(NBH(NH_LOC(nhx1,nx)))
     '          +NAT(NBH(NH_LOC(nhx1,nx))))),
     '          ns1=1,NST(NBH(NH_LOC(nhx1,nx)))+
     '          NAT(NBH(NH_LOC(nhx1,nx)))),nhx1=1,NH_LOC(0,nx))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !ns
          ENDDO !nhx
        ENDIF !PARAMTYPE
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('D_ZERE50')
      RETURN
 9999 CALL ERRORS('D_ZERE50',ERROR)
      CALL EXITS('D_ZERE50')
      RETURN 1
      END


