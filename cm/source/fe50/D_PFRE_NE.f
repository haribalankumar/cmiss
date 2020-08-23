      SUBROUTINE D_PFRE_NE(PARAMTYPE,IBT,IDO,INP,NAN,NBH,NBJ,
     '  NFF,NGAP,NHE,nhs,NMNO,NPNE,nr,ne,NRE,NSP,NW,nx,NXI,
     '  CE,CP,D_RE,D_RI3,D_TW,D_ZE,D_ZW,FEXT,PG,XE,XW,YG,ZE,ZW,ZW1,
     '  ERROR,*)

C#### Subroutine: D_PFRE_NE
C###  Description:
C###    D_PFRE_NE evaluates contribution to derivatives of element
C###    residuals wrt material parameters D_RE(na,4) (associated with
C###    the hydrostatic pressure var) for incompressible materials,
C###    arising from the stress constraintdue to pressure bcs.
C###    OR
C###    Evaluates contribution to derivatives of element residuals wrt
C###    geometric variables ES (associated with the hydrostatic
C###    pressure variable) for incompressible materials, arising from
C###    the stress constraint due to pressure b.c's.

C**** The following code is similar to PFRE_NE and should be kept in
C**** synch with PFRE_NE.

      IMPLICIT NONE
      INCLUDE 'acti00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'ptr00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NFF(6),NGAP(NIM,NBM),NHE,
     '  nhs,NMNO(1:2,0:NOPM),NPNE(NNM,NBFM),nr,ne,NRE(NEM),
     '  NSP(-2:2),NW,nx,NXI(-NIM:NIM,0:NEIM)
      REAL*8 CE(NMM),CP(NMM,NPM),D_RE(NSM,NHM,NOPM),D_RI3(NHM*NSM),
     '  D_TW(3,3,NHM*NSM),D_ZE(NSM,NHM),D_ZW(NHM,NUM,NHM*NSM),
     '  FEXT(NIFEXTM,NGM),
     '  PG(NSM,NUM,NGM,NBM),XE(NSM,NJM),XW(NJM,NUM),YG(NIYGM,NGM),
     '  ZE(NSM,NHM),ZW(NHM,NUM),ZW1(NHM,NUM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER adj_dirn,i,iface,IFE,j,k,mi,mj,mz,
     '  NB1,NBP,NCW,neadj,ng_elem,ngi1,ngi2,ngi3,nh1,nhs1,nhx1,ni,
     '  NITB,nj,noopti,ns,ns1,nz
      PARAMETER (NCW=35) !CW must be dimen.d the same size as CE array
      REAL*8 AXU(3,3),AZ,AZ1,AZL(3,3),AZL1(3,3),
     '  AZU(3,3),AZU1(3,3),CW(NCW),D_AZ,D_AZL(3,3),D_AZU(3,3),
     '  D_EG(3,3),DELTA_ZE,DETERM,DET_F_NU,
     '  DNUREFDNU(3,3),DNUREFDZ(3,3),D_TCNU(3,3),
     '  D_TCNUREF33,DXIXN(3,3),DXIZN(3,3),
     '  DXNZN(3,3),DZDNU(3,3),
     '  DZDNUREF(1,3),DZNXI(3,3),DZNXN(3,3),
     '  EG(3,3),GXL(3,3),GXU(3,3),
     '  RI1,RI2,RI3,RWX,SUM,TW(3,3),TWA,XI(3),ACTIVE_STRESS,
     &     DUMMY_AZL(3,3),DUMMY_AZU(3,3)
      CHARACTER CHAR1*1,STRESSTYPE*17
      LOGICAL ADJ_XI3_ELEM
      DATA DELTA_ZE/1.D-8/
C     OR 23-08-06 : Initialize arrays
      DATA DUMMY_AZL/9 * 0.0d1/
      DATA DUMMY_AZU/9 * 0.0d1/

      CALL ENTERS('D_PFRE_NE',*9999)
      CALL ASSERT(NCW.EQ.NMM,'>>Dimension of CW array (NCW)'
     '  //' must equal dimension of CE array (NMM)',ERROR,*9999)
      NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis function for pressure vars
      NITB=NIT(NBP)
      DO mi=1,NITB
        XI(MI)=0.5d0
      ENDDO !mi

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
     '    .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext press bc on Xi3=0 face
     '     iface.EQ.2.AND..NOT.ADJ_XI3_ELEM
     '    .AND.(NW.EQ.3.OR.NW.EQ.4).OR. !ext press bc on Xi3=1 face
     '    KTYP5A(nr).EQ.2) THEN !match norm stress on Xi3 faces
C old
C        IF(iface.EQ.1.AND.(NXI(-3,1).EQ.0.OR.nr.NE.NRE(NXI(-3,1)))
C     '    .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext press bc on Xi3=0 face
C     '     iface.EQ.2.AND.(NXI( 3,1).EQ.0.OR.nr.NE.NRE(NXI(3,1)))
C     '    .AND.(NW.EQ.3.OR.NW.EQ.4).OR. !ext press bc on Xi3=1 face
C     '    KTYP5A(nr).EQ.2) THEN !match norm stress on Xi3 faces

          IFE=iface+4
          XI(3)=DBLE(iface-1)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' Xi(3)='',I1,'' face (nf='',I3,'
     '        //''')'')') iface-1,NFF(IFE)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP
C new MPN 15-Apr-96: nodal interpolation of material params
          CALL CPXI(1,IBT,IDO,INP,NPNE,nr,nx,CE,CP,CW,XI,ERROR,*9999)
C old
C          DO il=1,ILT(1,nr,nx)
C            IF(ILP(il,1,nr,nx).EQ.3) THEN
C              CW(il)=0.0d0
C              DO nnbf=1,NNT(NBFF)
C                CW(il)=CW(il)+CP(IL,NPNE(nnbf,NBJ(1)))
C              ENDDO
C              CW(il)=CW(il)/DBLE(NNT(nbff))
C            ELSE IF(ILP(il,1,nr,nx).NE.3) THEN
C              CW(il)=CE(il)
C            ENDIF
C          ENDDO !il

C         Interpolate midwall geometric vars XW and derivs wrt Xi
          CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XW,XI,ERROR,*9999)
C         Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C         derivs (DXIXN) of Xi wrt Nu  coords.
          CALL XGMG(1,NITB,NBJ(1),nr,DXIXN,GXL,GXU,RWX,XW,ERROR,*9999)

C new MPN 4-May-96: new way of handling sheets
C         Get derivs of Xi wrt def Nu coords, DXIZN, and inv, DZNXI
          CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '      DXIZN,DZNXI,PG,XE,XW,XI,ZE,ZW,'Fibre',ERROR,*9999)
C         Calculate derivs of deformed Nu wrt undeformed Nu (DZNXN)
C         and inverse (DXNZN)
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
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(CHAR1,'(I1)') NITB
            WRITE(OP_STRING,'(''  DZNXN:'','//CHAR1(1:1)//'D12.4,'
     '        //'''  DXNZN:'','//CHAR1(1:1)//'D12.4,'//'/(8X,'
     '        //CHAR1(1:1)//'D12.4,8X,'//CHAR1(1:1)//'D12.4))')
     '        ((DZNXN(mi,ni),ni=1,NITB),(DXNZN(mi,ni),ni=1,NITB),
     '        mi=1,NITB)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP
C old
CC         Get derivs of Xi wrt undef Nu (body/fibre) coords,DXIXN
C          CALL DXIDXM(NBJ(1),nr,DXIXN,DXNXI,GXL,GXU,XW,'Fibre',ERROR,*9999)
CC         Interpolate dependent var.s ZG and derivs wrt Xi
C          CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
C     '      DXIX,ZE,ZW,XI,ERROR,*9999)
C          CALL ZGMG(NBH(NH_LOC(1,nx)),nr,GZ,GZL,GZU,ZW,ERROR,*9999)
C          CALL DXIDXM(NBJ(1),nr,DXIZN,DZNXI,GZL,GZU,XW,'Fibre',ERROR,*9999)
CC         Calculate derivs of deformed Nu wrt undeformed Nu (DZNXN)
CC         and inverse (DXNZN)
C          DO ni=1,NITB
C            DO mi=1,NITB
C              SUM1=0.0d0
C              SUM2=0.0d0
C              DO k=1,NITB
C                SUM1=SUM1+DZNXI(ni,k)*DXIXN(k,mi)
C                SUM2=SUM2+DXNXI(ni,k)*DXIZN(k,mi)
C              ENDDO !k
C              DZNXN(ni,mi)=SUM1
C              DXNZN(ni,mi)=SUM2
C            ENDDO !mi
C          ENDDO !ni
C end old

C         Interpolate dependent var.s ZW and derivs wrt Nu
          CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '      DXIXN,ZE,ZW,XI,ERROR,*9999)
C         Calculate deformed metric tensors wrt Nu (AZL,AZU)
          CALL ZGMG(NBH(NH_LOC(1,nx)),nr,AZ,AZL,AZU,ZW,ERROR,*9999)

          CALL ASSERT(KTYP51(nr).EQ.3.OR.KTYP51(nr).EQ.4,
     '      '>>Analytic derivatives only'
     '      //' implemented for 3D and membrane problems',ERROR,*9999)

C news VJ 10Dec2003 calc ng_elem
          NGI1=NGAP(1,NBP)
          NGI2=NGAP(2,NBP)
          NGI3=NGAP(3,NBP)
          ng_elem=(NGI1+1)/2 + ((NGI2-1)/2)*NGI1
     '      + (iface-1)*(NGI3-1)
     '      *NGI2*NGI1 !need integer division            
C newe VJ 10Dec2003

          IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN
C           Calc deriv's of contravariant cpts of 2nd Piola-Kirchhoff
C           stress tensor (D_TW - undeformed Nu coordinates) wrt
C           each of the material parameters
C           Note: AZ,AZL,AZU,EG are all independent of material params
C           so the calculation of D_AZ,D_AZL,D_AZU is not
C           needed here.
            IF(KTYP51(nr).EQ.3) THEN
              RI3=AZ !this saves calling ZGTG53 for RI3
              DO noopti=1,NTOPTI
C     OR 22-08-06 Apparently the values for arguments D_AZL and D_AZU
C     are not important, but %VAL(0) causes a warning, hence I replaced
C     %VAL(0) with DUMMY_AZL and DUMMY_AZU respectively
C                CALL D_ZGTG53(PARAMTYPE,NBH(NH_LOC(1,nx)),
C     '            NMNO(1,noopti),nr,nx,
C     '            AXU,AZ,AZL,AZU,CW,0.d0,%VAL(0),%VAL(0),D_EG,
C     '            D_RI3(noopti),D_TW(1,1,noopti),D_ZW(1,1,noopti),
C     '            EG,ZW,ERROR,*9999)
                CALL D_ZGTG53(PARAMTYPE,NBH(NH_LOC(1,nx)),NMNO(1,noopti)
     &               ,nr,nx, AXU,AZ,AZL,AZU,CW,0.d0,DUMMY_AZL,DUMMY_AZU
     &               ,D_EG, D_RI3(noopti),D_TW(1,1 ,noopti),D_ZW(1,1
     &               ,noopti), EG,ZW,ERROR,*9999)
              ENDDO !noopti
            ELSE IF(KTYP51(nr).EQ.4) THEN
              RI3=1.0d0 !this saves calling ZGTG54 for RI3
              DO noopti=1,NTOPTI
                CALL D_ZGTG54(PARAMTYPE,NMNO(1,noopti),nr,
     '            AXU,AZ,AZL,AZU,CW,D_AZ,D_AZL,D_AZU,D_EG,
     '            D_RI3(noopti),D_TW(1,1,noopti),EG,YG(1,ng_elem),
     '            ERROR,*9999)
              ENDDO !noopti
            ENDIF !KTYP51(nr)
          ELSE IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN
C           Calc contravariant cpts of 2nd Piola-Kirchhoff stress
C           tensor (TW) wrt undeformed Nu coordinates
            IF(KTYP51(nr).EQ.3) THEN
              WRITE(STRESSTYPE,'('''')')
              CALL ZGTG53(STRESSTYPE,NBH(NH_LOC(1,nx)),
     '          nr,nx,AXU,AZ,AZL,AZU,
     '          CW,EG,RI1,RI2,RI3,TW,XW,YG(1,ng_elem),ZW,ERROR,*9999)
            ELSE IF(KTYP51(nr).EQ.4) THEN
              CALL ZGTG54(NBH(NH_LOC(1,nx)),nr,AXU,AZ,AZL,AZU,
     '          CW,EG,RI1,RI2,RI3,TW,YG(1,ng_elem),ERROR,*9999)
            ENDIF !KTYP51(nr)
C           Calc deriv's of contravariant cpts of 2nd Piola-Kirchhoff
C           stress tensor (D_TW - undeformed Nu coordinates) and deriv
C           of third strain invariant (D_RI3) wrt each of the geometric
C           variables
            nhs1=0
            DO nhx1=1,NH_LOC(0,nx)
              nh1=NH_LOC(nhx1,nx)
              NB1=NBH(nh1)
              DO ns1=1,NST(NB1)+NAT(NB1)
                nhs1=nhs1+1
C               Calc D_ZW, deriv of def coords (ZW) wrt ele coords (ZE)
                D_ZE(ns1,nhx1)=1.0d0 !diff'ing wrt ns1,nh1 elem coord
                CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '            DXIXN,D_ZE,D_ZW(1,1,nhs1),XI,ERROR,*9999)
                D_ZE(ns1,nhx1)=0.0d0
C               Calculate derivs of deformed metric tensors, D_AZL and
C               D_AZU (Nu coords) and deriv of the det of AZL,
C               D_AZ wrt the current geom var using finite diff
C               approximations
C               Perturb deformed element coords
                ZE(ns1,nhx1)=ZE(ns1,nhx1)+DELTA_ZE
C               Interpolate perturbed dep vars ZW1 and derivs wrt Nu
                CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '            DXIXN,ZE,ZW1,XI,ERROR,*9999)
C               Reset deformed element coords
                ZE(ns1,nhx1)=ZE(ns1,nhx1)-DELTA_ZE
C               Calc deformed metric tensors for perturbed coords
                CALL ZGMG(NBH(NH_LOC(1,nx)),nr,AZ1,AZL1,AZU1,ZW1,
     '            ERROR,*9999)
C               Finite diff approx to derivs
                DO i=1,3
                  DO j=1,3
                    D_AZL(i,j)=(AZL1(i,j)-AZL(i,j))/DELTA_ZE
                    D_AZU(i,j)=(AZU1(i,j)-AZU(i,j))/DELTA_ZE
                  ENDDO !j
                ENDDO !i
                D_AZ=(AZ1-AZ)/DELTA_ZE
C               Calc derivs of TW wrt current geom var analytically
                IF(KTYP51(nr).EQ.3) THEN
                  CALL D_ZGTG53(PARAMTYPE,NBH(NH_LOC(1,nx)),0,nr,nx,
     '              AXU,AZ,AZL,AZU,CW,D_AZ,D_AZL,D_AZU,D_EG,
     '              D_RI3(nhs1),D_TW(1,1,nhs1),D_ZW(1,1,nhs1),EG,ZW,
     '              ERROR,*9999)
                ELSE IF(KTYP51(nr).EQ.4) THEN
                  CALL D_ZGTG54(PARAMTYPE,0,nr,
     '              AXU,AZ,AZL,AZU,CW,D_AZ,D_AZL,D_AZU,D_EG,
     '              D_RI3(nhs1),D_TW(1,1,nhs1),EG,YG(1,ng_elem),
     '              ERROR,*9999)
                ENDIF !KTYP51(nr)
              ENDDO !ns1
            ENDDO !nhx1
          ENDIF !paramtype=material/geometric

          IF(KTYP53(nr).EQ.3) THEN !Active stress component included
C!!!  Pick Gauss pt nearest to centre of current face in
C!!!        current element. To be strictly correct need to either
C!!!        store FEXT at face basis Gauss pts or somehow interpolate
C!!!        element Gauss pt values to the central point on the face
C VJ 10Dec2003 commented out as ng_elem calculated before ZGTG53 called
C            NGI1=NGAP(1,NBP)
C            NGI2=NGAP(2,NBP)
C            NGI3=NGAP(3,NBP)
C            ng_elem=(NGI1+1)/2 + ((NGI2-1)/2)*NGI1
C     '        + (iface-1)*(NGI3-1)*NGI2*NGI1 !need integer division
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
              CALL EVALASC(ne,ng_elem,%VAL(NQNE_PTR),%VAL(YQS_PTR)
     &             ,%VAL(RCQS_SPATIAL_PTR),%VAL(ICQS_SPATIAL_PTR)
     &             ,ASC_ARRAYNAME(nr) ,ASC_CELLVARINDEX(nr)
     &             ,ACTIVE_STRESS,ERROR,*9999)
            ELSE
              ACTIVE_STRESS = YG(1,ng_elem)
            ENDIF
C     09-Dec-1989: NOTE: Don't have to define a separate face array
            
            CALL ZGTG5A(NBH(NH_LOC(1,nx)),nr,FEXT(1,ng_elem),
     '        DXNZN,DZNXN,DET_F_NU,TW,TWA,ACTIVE_STRESS,ERROR,*9999)
C!!!        NOTE: The code below assumes that the active fibre stress
C!!!        component is independent of the
C!!!        passive material/geometrical parameters
          ENDIF !KTYP53(nr).EQ.3 (active stresses)

          IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN
C new MPN 17Mar97: correcting pressure bc transformations
C           Compute def anatomical fibre vects wrt rc coords at XI
            CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '        DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
     '        PG,XE,XW,XI,ZE,ZW,.TRUE.,ERROR,*9999)
C           Compute deformed fibre ref vectors wrt rc coords at XI
            CALL FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH,0,
     '        NHE,NITB,nr,nx,DZDNUREF(1,1),DZDNUREF(1,2),
     '        DZDNUREF(1,3),PG,XI,ZE,ZW,ERROR,*9999)
            CALL INVERT(NITB,DZDNUREF,DNUREFDZ,DETERM)
            mz=3 !only need derivs wrt 3rd wall coord (Xi1/2 face norm)
            DO nz=1,3
              SUM=0.0d0
              DO mi=1,3
                SUM=SUM+DNUREFDZ(mz,mi)*DZDNU(mi,nz)
              ENDDO !mi
              DNUREFDNU(mz,nz)=SUM
            ENDDO !nz

            DO noopti=1,NTOPTI
C             Get derivs of Physical Cauchy stress tensor wrt deformed
C             nu-material coordinates (D_TCNU) from D_TW
              DO mz=1,3
                DO nz=1,3
                  SUM=0.0d0
                  DO mj=1,NITB
                    DO nj=1,NITB
                      SUM=SUM+DZNXN(mz,mj)*D_TW(mj,nj,noopti)*
     '                  DZNXN(nz,nj)
                    ENDDO !nj
                  ENDDO !mj
                  D_TCNU(mz,nz)=SUM/DSQRT(RI3)
                ENDDO !nz
              ENDDO !mz
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
C new MPN 18Mar98: fixed format
                WRITE(OP_STRING,'('' D_TCNU:'',8X,3D12.4,'
     '            //'/(16X,3D12.4))')
     '            ((D_TCNU(mz,nz),nz=1,3),mz=1,3)
C old                WRITE(OP_STRING,'('' D_TCNU:'',12X,3D12.4,'
C old     '            //'/(16X,3D12.4))')
C old     '            ((D_TCNU(mz,nz),nz=1,3),mz=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF !DOP
C             Compute cmpts of physical Cauchy stress wrt deformed
C             fibre reference coords using wall coord derivs
C             wrt nu coords
              D_TCNUREF33=0.0d0
              DO mi=1,3
                DO ni=1,3
                  D_TCNUREF33=D_TCNUREF33+DNUREFDNU(3,mi)*D_TCNU(mi,ni)*
     '              DNUREFDNU(3,ni)
                ENDDO !ni
              ENDDO !mi
C old MPN 17Mar97: incorrect pressure bc transformations
CC             Get deriv of physical Cauchy stress tensor D_TC wrt
CC             deformed nu-material coordinates from D_TW
C              DO mz=1,3
C                DO nz=1,3
C                  SUM=0.0d0
C                  DO mj=1,NITB
C                    DO nj=1,NITB
C                      SUM=SUM+DZNXN(mz,mj)*D_TW(mj,nj,noopti)
C     '                  *DZNXN(nz,nj)
C                    ENDDO !nj
C                  ENDDO !mj
C                  D_TC(mz,nz)=SUM/DSQRT(RI3)
C                ENDDO !nz
C              ENDDO !mz
C              IF(DOP) THEN
CC$              call mp_setlock()
C                WRITE(OP_STRING,'('' D_TC:'',12X,3D12.4,/'//
C     '            '(16X,3D12.4))')
C     '            ((D_TC(mz,nz),nz=1,3),mz=1,3)
C                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
C              ENDIF !DOP
C              IF(JTYP9.GE.2) THEN !imbric (+ sheet) angles defined
CC               Compute def anatomical fibre vects wrt rc coords at XI
C                CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
C     '            DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
C     '            PG,XE,XW,XI,ZE,ZW,ERROR,*9999)
C                CALL INVERT(NITB,DZDNU,DNUDZ,DETERM)
CC               Compute deriv of Cauchy stress cmpts wrt rc coords
CC               by rotating def material coord system into rc coords
C                DO mi=1,3
C                  DO ni=1,3
C                    SUM=0.0d0
C                    DO mj=1,3
C                      DO nj=1,3
C                        SUM=SUM+DZDNU(mi,mj)*D_TC(mj,nj)*DZDNU(ni,nj)
C                      ENDDO !nj
C                    ENDDO !mj
C                    D_TCRC(mi,ni)=SUM
C                  ENDDO !ni
C                ENDDO !mi
C                IF(DOP) THEN
CC$                call mp_setlock()
C                  WRITE(OP_STRING,'('' D_TCRC:'',12X,3D12.4,'
C     '              //'/(18X,3D12.4))') ((D_TCRC(mi,ni),ni=1,3),mi=1,3)
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
C                ENDIF !DOP
CC               Compute deformed fibre ref vectors wrt rc coords at XI
C                CALL FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH,0,
C     '            NHE,NITB,nr,nx,DZDNUREF(1,1),DZDNUREF(1,2),
C     '            DZDNUREF(1,3),PG,XI,ZE,ZW,ERROR,*9999)
C                CALL INVERT(NITB,DZDNUREF,DNUREFDZ,DETERM)
CC               Compute deriv of Cauchy stress cmpts wrt deformed
CC               fibre reference coords by rotating rc coord system
CC               into deformed fibre ref coords
CC               NOTE: only need D_TCNUREF(3,3)
C                mi=3
C                ni=3
C                SUM=0.0d0
C                DO mj=1,3
C                  DO nj=1,3
C                    SUM=SUM+DNUREFDZ(mi,mj)*D_TCRC(mj,nj)
C     '                *DNUREFDZ(ni,nj)
C                  ENDDO !nj
C                ENDDO !mj
C                D_TCNUREF(mi,ni)=SUM
Cc                IF(DOP) THEN
CcC$                call mp_setlock()
Cc                  WRITE(OP_STRING,'('' D_TCNUREF:'',12X,3D12.4,'
Cc     '              //'/(21X,3D12.4))')
Cc     '              ((D_TCNUREF(mi,ni),ni=1,3),mi=1,3)
Cc                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CcC$                call mp_unsetlock()
Cc                ENDIF !DOP
C              ELSE !at most fibres defined
CC!!! WARNING: When the definition of the sheet angle in routine
CC!!!          MAT_VEC_ROTATE (FE02.F) changes from GAMA-PI/2 to
CC!!!          GAMA, the following 3 lines must be replaced
CC!!!          with the 3 below them, since D_TC is defined
CC!!!          wrt the nu-material coordinates
CC               D_TCNUREF33 is aligned with D_TC(2,2) for
CC               zero imbric/sheet angles
C                D_TCNUREF(3,3)=D_TC(2,2)
CC new           D_TCNUREF33 is aligned with D_TC(3,3) for
CC new           zero imbric/sheet angles
CC new           D_TCNUREF(3,3)=D_TC(3,3)
C              ENDIF !JTYP9.GE.2
C              D_TCNUREF33=D_TCNUREF(3,3)
CC older
CC             Get deriv of Physical Cauchy stress D_TC(3,3)
CC             wrt deformed Nu_3 from D_TW
CC              D_TC(3,3)=0.0d0
CC              DO nj=1,NITB
CC                DO mj=1,NITB
CC                  D_TC(3,3)=D_TC(3,3)+DZNXN(3,nj)*DZNXN(3,mj)
CC     '              *D_TW(nj,mj,noopti)
CC                ENDDO !mj
CC              ENDDO !nj
CC              D_TC(3,3)=D_TC(3,3)/DSQRT(RI3)
C end old
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' D_TCNUREF33='',D12.4)') D_TCNUREF33
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF !DOP

              ns=NSP(iface)
C old
C             IF(iface.EQ.1.AND.(NXI(-3,1).EQ.0.OR.nr.NE.NRE(NXI(-3,1)))
C     '         .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext p bc on Xi3=0 face
C     '          iface.EQ.2.AND.(NXI( 3,1).EQ.0.OR.nr.NE.NRE(NXI( 3,1)))
C     '         .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !ext p bc on Xi3=1 face
C               match norm stress for press bcs
C new MPN 16APR97: to prevent access violations
              IF(iface.EQ.1.AND..NOT.ADJ_XI3_ELEM
     '          .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext pres bc on Xi3=0 face
     '          iface.EQ.2.AND..NOT.ADJ_XI3_ELEM
     '          .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !press bc on Xi3=1 face
                D_RE(ns,NH_LOC(NH_LOC(0,nx),nx),noopti)=D_TCNUREF33
              ELSE !no external pressure bc applied on current face
                IF(KTYP5A(nr).EQ.2) THEN !match norm stress on Xi3 faces
                  IF(iface.EQ.1) THEN !Xi3=0 face
                    D_RE(ns,NH_LOC(NH_LOC(0,nx),nx),noopti)= D_TCNUREF33
                  ELSE IF(iface.EQ.2) THEN !Xi3=1 face
                    D_RE(ns,NH_LOC(NH_LOC(0,nx),nx),noopti)=-D_TCNUREF33
                  ENDIF !iface
                ENDIF !KTYP5A(nr)=2
              ENDIF !ext press bc

            ENDDO !noopti

          ELSE IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN

            CALL ASSERT(.FALSE.,'>>ERROR: Old code, needs updating '//
     '        'for new material vector calcs',ERROR,*9999)

            IF(NW.EQ.3) THEN
              nhs=nhs+2
            ELSE
              nhs=nhs+1
            ENDIF
            nhs1=0
            DO nhx1=1,NH_LOC(0,nx)
              nh1=NH_LOC(nhx1,nx)
              NB1=NBH(nh1)
              DO ns1=1,NST(NB1)+NAT(NB1)
                nhs1=nhs1+1
C               Perturb deformed element coords
                ZE(ns1,nhx1)=ZE(ns1,nhx1)+DELTA_ZE
C               Interpolate perturbed dep var.s ZW1 and derivs wrt Xi
                CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '            DXIXN,ZE,ZW1,XI,ERROR,*9999)
C               Reset deformed element coords
                ZE(ns1,nhx1)=ZE(ns1,nhx1)-DELTA_ZE

C old MPN 17Mar97: incorrect pressure bc transformations
C                CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
C     '            DXIZN1,DZNXI1,PG,XE,XW,XI,ZE,ZW1,'Fibre',ERROR,*9999)
CC               Calculate derivs of deformed Nu wrt undeformed Nu
CC               for perturbed elem coords
C                DO ni=1,NITB
C                  DO mi=1,NITB
C                    SUM1=0.0d0
C                    SUM2=0.0d0
C                    DO k=1,NITB
C                      SUM1=SUM1+DZNXI1(ni,k)*DXIXN(k,mi)
C                      SUM2=SUM2+DXNXI(ni,k)*DXIZN1(k,mi)
C                    ENDDO !k
C                    DZNXN1(ni,mi)=SUM1
C                  ENDDO !mi
C                ENDDO !ni
CC               Finite diff approx to derivatives
C                DO i=1,3
C                  DO j=1,3
C                    D_DZNXN(i,j)=(DZNXN1(i,j)-DZNXN(i,j))/DELTA_ZE
C                  ENDDO !j
C                ENDDO !i
C end old

C!!! needs active stuff here also

C!!! needs new mat vec calcs here

C!!! THE FOLLOWING CODE NEEDS UPDATING

CC               Calc deriv of Cauchy Stress wrt current geometric var.
C                D_TC(3,3)=0.0d0
C                DO nj=1,NITB
C                  DO mj=1,NITB
C                    D_TC(3,3)=D_TC(3,3)
C     '                +D_DZNXN(3,nj)*DZNXN(3,mj)*TW(nj,mj)/DSQRT(RI3)
C     '                +DZNXN(3,nj)*D_DZNXN(3,mj)*TW(nj,mj)/DSQRT(RI3)
C     '                +DZNXN(3,nj)*DZNXN(3,mj)*D_TW(nj,mj,nhs1)
C     '                /DSQRT(RI3)
C     '                -DZNXN(3,nj)*DZNXN(3,mj)*TW(nj,mj)*D_RI3(nhs1)
C     '                /(2.0d0*DSQRT(RI3*RI3*RI3))
C                  ENDDO !mj
C                ENDDO !nj
C
C                IF(iface.EQ.1.AND..NOT.ADJ_XI3_ELEM
C     '            .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext pres bc on Xi3=0 face
C     '            iface.EQ.2.AND..NOT.ADJ_XI3_ELEM
C     '            .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !press bc on Xi3=1 face
CC                 match norm stress for press bcs
C                  ES(nhs,nhs1)= D_TC(3,3)
C                ELSE !no external pressure bc applied on current face
C                  IF(KTYP5A(nr).EQ.2) THEN !match norm stress on Xi3 faces
C                    IF(iface.EQ.1) THEN !Xi3=0 face
C                      ES(nhs,nhs1)= D_TC(3,3)
C                    ELSE IF(iface.EQ.2) THEN !Xi3=1 face
C                      ES(nhs,nhs1)=-D_TC(3,3)
C                    ENDIF !iface
C                  ENDIF !KTYP5A(nr)=2
C                ENDIF !ext press bc
              ENDDO !ns1
            ENDDO !nhx1
          ENDIF !paramtype=material/geometric
        ENDIF !ext press bc or KTYP5A(nr)=2
      ENDDO !iface

      CALL EXITS('D_PFRE_NE')
      RETURN
 9999 CALL ERRORS('D_PFRE_NE',ERROR)
      CALL EXITS('D_PFRE_NE')
      RETURN 1
      END


