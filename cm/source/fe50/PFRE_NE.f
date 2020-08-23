      SUBROUTINE PFRE_NE(IBT,IDO,INP,NAN,NBH,NBJ,NFF,NGAP,NHE,
     '  NPNE,nr,ne,NRE,NSP,NW,nx,NXI,
     '  CE,CP,FEXT,PF,PG,RE,XE,XW,YG,ZE,ZW,ERROR,*)

C#### Subroutine: PFRE_NE
C###  Description:
C###    PFRE_NE evaluates contribution to element residuals RE(na,4)
C###    (associated with the hydrostatic pressure variable) for
C###    incompressible materials, arising from the stress constraint
C###    due to pressure bcs.

C**** Note: XW,ZW,TW,CW etc are element geometric, dependent, stress &
C****       material parameters etc interpolated at the centre of the
C****       pressure-loaded wall.
C**** 30APR89: NSP(i), i=1..2 are the hydrostatic press params coupled
C**** with the pressure contraint on the inside and the outside of the
C**** element, respectively.
C**** MPN 11-Jan-95: This routine is used when the hydrostatic pressure
C****                is interpolated using element based variables.

      IMPLICIT NONE
      INCLUDE 'acti00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ptr00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NFF(6),NGAP(NIM,NBM),NHE,
     '  NPNE(NNM,NBFM),nr,ne,NRE(NEM),NSP(-2:2),NW,
     '  nx,NXI(-NIM:NIM,0:NEIM)
      REAL*8 CE(NMM),CP(NMM,NPM),FEXT(NIFEXTM,NGM),
     '  PF(2),PG(NSM,NUM,NGM,NBM),RE(NSM),
     '  XE(NSM,NJM),XW(NJM,NUM),YG(NIYGM,NGM),ZE(NSM,NHM),ZW(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER adj_dirn,iface,IFE,k,mi,mj,mz,NBP,NCW,neadj,nf,
     '  ng_elem,NGI1,NGI2,NGI3,ni,NITB,nj,ns,nz
      PARAMETER (NCW=35) !CW must be dimen.d the same size as CE array
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CW(NCW),DETERM,DET_F_NU,
     '  DNUREFDNU(3,3),DNUREFDZ(3,3),DXIXN(3,3),DXIZN(3,3),DXNZN(3,3),
     '  DZDNU(3,3),DZDNUREF(3,3),DZNXI(3,3),DZNXN(3,3),
     '  EG(3,3),GXL(3,3),GXU(3,3),RI1,RI2,RI3,RWX,SUM,
     '  TCNU(3,3),TCNUREF33,TW(3,3),TWA,XI(3),ACTIVE_STRESS
      CHARACTER CHAR1*1,STRESSTYPE*17
      LOGICAL ADJ_XI3_ELEM
      DATA STRESSTYPE/' '/

      CALL ENTERS('PFRE_NE',*9999)
      CALL ASSERT(NCW.EQ.NMM,'>>Dimension of CW array (NCW)'
     '  //' must equal dimension of CE array (NMM)',ERROR,*9999)
      CALL ASSERT(KTYP53(nr).GT.1,'stresses must be referred to '
     '  //'Nu coords',ERROR,*9999)
      NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis function for pressure vars
      NITB=NIT(NBP)
      DO mi=1,NITB
        XI(mi)=0.50d0
      ENDDO !mi

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>>PFRE_NE  diagnostic op'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

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
          nf=NFF(IFE)
          XI(3)=DBLE(iface-1)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' Xi(3)='',I1,'' face (nf='',I3,'').'
     '        //' PF ='',D12.4)') iface-1,nf,PF(iface)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP

          CALL CPXI(1,IBT,IDO,INP,NPNE,nr,nx,CE,CP,CW,XI,ERROR,*9999)
C old
C          DO il=1,ILT(1,nr,nx)
C            IF(ILP(il,1,nr,nx).NE.1.OR.ILP(il,1,nr,nx).NE.2) THEN
CC             constant spatially or defined by elements
C              CW(il)=CE(il)
C            ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
C              CW(il)=0.0d0
C              DO nnbf=1,NNT(NBFF)
C                CW(il)=CW(il)+CP(il,NPNE(nnbf,NBJ(1)))
C              ENDDO
C              CW(il)=CW(il)/DBLE(NNT(NBFF))
C            ELSE IF(ILP(il,1,nr,nx).EQ.4) THEN !defined by Gauss points
C              CALL ASSERT(.FALSE.,' >>> Gauss pt mat param variation '
C     '          //'not implemented',ERROR,*9999)
C            ENDIF
C          ENDDO !il

C         Interpolate midwall geometric vars XW and derivs wrt Xi
          CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XW,XI,ERROR,*9999)
C         Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C         derivs of Xi wrt undef Nu coords, DXIXN, (JP=1).
          CALL XGMG(1,NITB,NBJ(1),nr,DXIXN,GXL,GXU,RWX,XW,ERROR,*9999)

C new MPN 4-May-96: new way of handling sheets
C         Get derivs of Xi wrt def Nu coords, DXIZN, and inv, DZNXI
          CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '      DXIZN,DZNXI,PG,XE,XW,XI,ZE,ZW,'Fibre',ERROR,*9999)
C old
CC         Interpolate dependent var.s ZG and derivs wrt Xi
C          CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
C     '      DXIX,ZE,ZW,XI,ERROR,*9999)
C          CALL ZGMG(NBH(NH_LOC(1,nx)),nr,GZ,GZL,GZU,ZW,ERROR,*9999)
C          CALL DXIDXM(NBH(NH_LOC(1,nx)),nr,DXIZN,DZNXI,GZL,GZU,XW,
C     '      'Fibre',ERROR,*9999)
C end old
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

C         Interpolate dependent var.s ZW and derivs wrt Nu (JP=1)
          CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '      DXIXN,ZE,ZW,XI,ERROR,*9999)
C         Calculate deformed metric tensors wrt Nu (AZL,AZU)
          CALL ZGMG(NBH(NH_LOC(1,nx)),nr,AZ,AZL,AZU,ZW,ERROR,*9999)

C         Get contravariant cpts of 2nd Piola-Kirchhoff stress
C         tensor (TW) wrt undeformed Nu coordinates
C news VJ 10Dec2003 calc ng_elem
          NGI1=NGAP(1,NBP)
          NGI2=NGAP(2,NBP)
          NGI3=NGAP(3,NBP)
          ng_elem=(NGI1+1)/2 + ((NGI2-1)/2)*NGI1
     '      + (iface-1)*(NGI3-1)*NGI2*NGI1 !need integer division            
C newe VJ 10Dec2003
          IF(KTYP51(nr).EQ.1) THEN
            CALL ZGTG51(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,CW,
     '        RI1,RI2,RI3,TW,YG(1,ng_elem),ZW,ERROR,*9999)
          ELSE IF(KTYP51(nr).EQ.2) THEN
            CALL ZGTG52(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,CW,
     '        RI1,RI2,RI3,TW,YG(1,ng_elem),ZW,ERROR,*9999)
          ELSE IF(KTYP51(nr).EQ.3) THEN
            CALL ZGTG53(STRESSTYPE,NBH(NH_LOC(1,nx)),
     '        nr,nx,AXU,AZ,AZL,AZU,CW,EG,
     '        RI1,RI2,RI3,TW,XW,YG(1,ng_elem),ZW,ERROR,*9999)
          ELSE IF(KTYP51(nr).EQ.4) THEN
            CALL ZGTG54(NBH(NH_LOC(1,nx)),nr,AXU,AZ,AZL,AZU,CW,EG,
     '        RI1,RI2,RI3,TW,YG(1,ng_elem),ERROR,*9999)
          ELSE IF(KTYP51(nr).EQ.5) THEN !string
            CALL ZGTG55(nr,AZL,CW,EG,TW,ERROR,*9999)
          ENDIF

          IF(KTYP53(nr).EQ.3) THEN !Active stress component included
C!!!        Pick Gauss pt nearest to centre of current face in
C!!!        current element. To be strictly correct need to either
C!!!        store FEXT at face basis Gauss pts or somehow interpolate
C!!!        element Gauss pt values to the central point on the face
C VJ 10Dec2003 commented out ng_elem calc as already done in 
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
            CALL ZGTG5A(NBH(NH_LOC(1,nx)),nr,FEXT(1,ng_elem), DXNZN
     &           ,DZNXN,DET_F_NU,TW,TWA,ACTIVE_STRESS,ERROR,*9999)
C            CALL ZGTG5A(NBH(NH_LOC(1,nx)),nr,FEXT(1,ng_elem), DXNZN
C     &           ,DZNXN,DET_F_NU,TW,TWA,YG(1,ng_elem),ERROR,*9999)
          ENDIF                 !KTYP53(nr).EQ.3 (active stresses)

C new MPN 17Mar97: correcting pressure bc transformations
C         Get Physical Cauchy stress tensor wrt deformed
C         nu-material coordinates (TCNU) from TW
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
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
C new MPN 18Mar98: fixed format
            WRITE(OP_STRING,'('' TCNU:'',10X,3D12.4,/(16X,3D12.4))')
     '        ((TCNU(mz,nz),nz=1,3),mz=1,3)
C old            WRITE(OP_STRING,'('' TCNU:'',12X,3D12.4,/(16X,3D12.4))')
C old     '        ((TCNU(mz,nz),nz=1,3),mz=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP
C         Compute def anatomical fibre vects wrt rc coords at XI
          CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '      DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
     '      PG,XE,XW,XI,ZE,ZW,.TRUE.,ERROR,*9999)
C         Compute deformed fibre ref vectors wrt rc coords at XI
          CALL FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH,0,
     '      NHE,NITB,nr,nx,DZDNUREF(1,1),DZDNUREF(1,2),
     '      DZDNUREF(1,3),PG,XI,ZE,ZW,ERROR,*9999)
          CALL INVERT(NITB,DZDNUREF,DNUREFDZ,DETERM)
          mz=3 !only need derivs wrt third wall coord (Xi1/2 face norm)
          DO nz=1,3
            SUM=0.0d0
            DO mi=1,3
              SUM=SUM+DNUREFDZ(mz,mi)*DZDNU(mi,nz)
            ENDDO !mi
            DNUREFDNU(mz,nz)=SUM
          ENDDO !nz
C         Compute cmpts of physical Cauchy stress wrt deformed
C         fibre reference coords using wall coord derivs
C         wrt nu coords
          TCNUREF33=0.0d0
          DO mi=1,3
            DO ni=1,3
              TCNUREF33=TCNUREF33+DNUREFDNU(3,mi)*TCNU(mi,ni)*
     '          DNUREFDNU(3,ni)
            ENDDO !ni
          ENDDO !mi
C old MPN 17Mar97: incorrect pressure bc transformations
CC         Get Physical Cauchy stress tensor TC wrt deformed
CC         nu-material coordinates from TW
C          DO mz=1,3
C            DO nz=1,3
C              SUM=0.0d0
C              DO mj=1,NITB
C                DO nj=1,NITB
C                  SUM=SUM+DZNXN(mz,mj)*TW(mj,nj)*DZNXN(nz,nj)
C                ENDDO !nj
C              ENDDO !mj
C              TC(mz,nz)=SUM/DSQRT(RI3)
C            ENDDO !nz
C          ENDDO !mz
C          IF(DOP) THEN
CC$          call mp_setlock()
C              WRITE(OP_STRING,'('' TC:'',12X,3D12.4,/(16X,3D12.4))')
C     '        ((TC(mz,nz),nz=1,3),mz=1,3)
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
C          ENDIF !DOP
C          IF(JTYP9.GE.2) THEN !imbric (+ sheet) angles defined
CC           Compute def anatomical fibre vects wrt rc coords at XI
C            CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
C     '        DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
C     '        PG,XE,XW,XI,ZE,ZW,ERROR,*9999)
CC           Compute cmpts of Cauchy stress tensor wrt rc coords
CC           by rotating def material coord system into rc coords
C            DO mi=1,3
C              DO ni=1,3
C                SUM=0.0d0
C                DO mj=1,3
C                  DO nj=1,3
C                    SUM=SUM+DZDNU(mi,mj)*TC(mj,nj)*DZDNU(ni,nj)
C                  ENDDO !nj
C                ENDDO !mj
C                TCRC(mi,ni)=SUM
C              ENDDO !ni
C            ENDDO !mi
C            IF(DOP) THEN
CC$            call mp_setlock()
C              WRITE(OP_STRING,'('' TCRC:'',12X,3D12.4,'
C     '          //'/(18X,3D12.4))') ((TCRC(mi,ni),ni=1,3),mi=1,3)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
C            ENDIF !DOP
CC           Compute deformed fibre ref vectors wrt rc coords at XI
C            CALL FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH,0,
C     '        NHE,NITB,nr,nx,DZDNUREF(1,1),DZDNUREF(1,2),
C     '        DZDNUREF(1,3),PG,XI,ZE,ZW,ERROR,*9999)
C            CALL INVERT(NITB,DZDNUREF,DNUREFDZ,DETERM)
CC           Compute cmpts of Cauchy stress tensor wrt deformed
CC           fibre reference coords by rotating rc coord system
CC           into deformed fibre ref coords
CC           NOTE: only need TCNUREF(3,3)
C            mi=3
C            ni=3
C            SUM=0.0d0
C            DO mj=1,3
C              DO nj=1,3
C                SUM=SUM+DNUREFDZ(mi,mj)*TCRC(mj,nj)*DNUREFDZ(ni,nj)
C              ENDDO !nj
C            ENDDO !mj
C            TCNUREF(mi,ni)=SUM
Cc            IF(DOP) THEN
CcC$            call mp_setlock()
Cc              WRITE(OP_STRING,'('' TCNUREF:'',12X,3D12.4,'
Cc     '          //'/(21X,3D12.4))')
Cc     '        ((TCNUREF(mi,ni),ni=1,3),mi=1,3)
Cc              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CcC$            call mp_unsetlock()
Cc            ENDIF !DOP
C          ELSE !at most fibres defined
C            TCNUREF33 is aligned with TC(3,3) for
C            zero imbric/sheet angles
C            TCNUREF(3,3)=TC(3,3)
C          ENDIF !JTYP9.GE.2
C          TCNUREF33=TCNUREF(3,3)
CC older
CCC         Get Physical Cauchy stress TC33 wrt deformed Nu_3 from TW
CC          TC33=0.0d0
CC          DO nj=1,NITB
CC            DO mj=1,NITB
CC              TC33=TC33+DZNXN(3,nj)*DZNXN(3,mj)*TW(nj,mj)
CC            ENDDO !mj
CC          ENDDO !nj
CC          TC33=TC33/DSQRT(RI3)
C end old
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' TCNUREF33='',D12.4)') TCNUREF33
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP

          ns=NSP(iface)
C old
C          IF(iface.EQ.1.AND.(NXI(-3,1).EQ.0.OR.nr.NE.NRE(NXI(-3,1)))
C     '      .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext. press bc on Xi3=0 face
C     '       iface.EQ.2.AND.(NXI( 3,1).EQ.0.OR.nr.NE.NRE(NXI( 3,1)))
C     '      .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !ext pres bc on Xi3=1 face
C new MPN 16APR97: to prevent access violations
          IF(iface.EQ.1.AND..NOT.ADJ_XI3_ELEM
     '      .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext pres bc on Xi3=0 face
     '      iface.EQ.2.AND..NOT.ADJ_XI3_ELEM
     '      .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !press bc on Xi3=1 face
            RE(ns)=TCNUREF33+PF(iface) !match norm. stress for press bcs
          ELSE !no external pressure bc applied on current face
            IF(KTYP5A(nr).EQ.2) THEN !match norm stress on Xi3 faces
C             Replace the incomp constraint with an explicit
C             normal Cauchy stress continuity constraint
              IF(iface.EQ.1) THEN !Xi3=0 face
                RE(ns)= TCNUREF33
              ELSE IF(iface.EQ.2) THEN !Xi3=1 face
                RE(ns)=-TCNUREF33
              ENDIF !iface
            ENDIF !KTYP5A(nr)=2
          ENDIF !ext press bc

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' RE('',I2,'')= '',D12.3)') ns,RE(ns)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP
        ENDIF !ext press bc or KTYP5A(nr)=2
      ENDDO !iface

      CALL EXITS('PFRE_NE')
      RETURN
 9999 CALL ERRORS('PFRE_NE',ERROR)
      CALL EXITS('PFRE_NE')
      RETURN 1
      END


