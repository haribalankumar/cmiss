      SUBROUTINE ZETX50_1(COORDS,CSTYPE,STRESSTYPE,IBT,IDO,INP,
     '  NAN,NBH,NBJ,ng,NHE,NPNE,nr,ne,nx,CE,CG,CP,FEXT,PG,PHI,PST,
     '  RGX,RGX2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XI,YG,ZE,ZG,ERROR,*)

C#### Subroutine: ZETX50
C###  Description:
C###    ZETX50 calculates 2nd Piola-Kirchhoff, Cauchy and  Nominal
C###    stresses with respect to 'Reference' or 'Fibre' coords
C###    (as specified by COORDS) at position XI in current element
C###    if ng=0 else at Gauss point ng.  STRESSTYPE ('Total', 'Passive'
C###    or 'Active') denotes the components of stress to be computed.

C**** Since base vectors of theta coords are not orthonormal, cpts
C**** of Cauchy and Nominal stress are converted to 'physical' values.
C**** AZ,AZL,AZU  are deformed metric tensors wrt undeformed coords
C**** RI1,RI2,RI3 are principal invariants of AZL
C**** AXU  are contravariant cpts of undeformed metric tensor
C**** XG   are undeformed theta coords and derivs wrt Xi
C**** ZG   are deformed theta coords and derivs wrt undeformed coords
C**** TG   are tensor cpts of 2nd Piola-Kirchhoff stresses
C**** TN   are physical cpts of Nominal stresses
C**** TNA  is the physical (Cauchy) active stress component
C**** TC   are physical cpts of Cauchy stress
C**** PST  are principal stresses
C**** RM   is the modal matrix whose cols are the eigenvectors
C****      associated with PST
C**** PHI  are the Euler angles of the principal stresses
C**** AZLZ are deformed Eulerian metrics (wrt deformed theta coords)

      IMPLICIT NONE
      INCLUDE 'acti00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ptr00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ng,NHE,NPNE(NNM,NBFM),
     '  nr,ne,nx
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),FEXT(NIFEXTM),
     '  PG(NSM,NUM,NGM,NBM),PHI(3),PST(3),RGX,RGX2D,RGZ,RGZ2D,RM(3,3),
     '  TC(3,3),TG(3,3),TN(3,3),TNA,XE(NSM,NJM),XG(NJM,NUM),
     '  XI(3),YG(NIYGM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER COORDS*(*),CSTYPE*(*),ERROR*(*),STRESSTYPE*(*)
!    Local Variables
      INTEGER i,il,j,k,mi,mj,mjj,mix,mz,
     '  nb,NCW,nhx,ni,NITB,nj,njj,NU1(0:3),nix,nz
      PARAMETER (NCW=35) !CW must be dimen.d the same size as CE array
      REAL*8 AA,AXU(3,3),AZ,AZL(3,3),AZLZ(3,3),AZU(3,3),CW(NCW), DETERM
     &     ,DET_DZDX,DXDZ(3,3),DXIXJ(3,3),DXIXN(3,3),DXIZN(3,3), DXJXN(3
     &     ,3),DXNXI(3,3),DXNXJ(3,3),DZDX(3,3),DZNXI(3,3),EG(3,3), G1,G3
     &     ,GXL(3,3),GXU(3,3),GZ,GZL(3,3),GZU(3,3), RC,RI1,RI2,RI3,RR
     &     ,RWX,SLX,SMX,SUM,TGX(3,3),TOL,VALTMP, ACTIVE_STRESS
      DATA NU1/1,2,4,7/

      CALL ENTERS('ZETX50',*9999)
      CALL ASSERT(COORDS.NE.'Wall','>>ERROR: Not implemented '
     '  //'for wall coordinates (yet).',ERROR,*9999)
      CALL ASSERT(NCW.EQ.NMM,'>>Dimension of CW array (NCW)'
     '  //' must equal dimension of CE array (NMM)',ERROR,*9999)
      nb=NBH(NH_LOC(1,nx))
      NITB=NIT(nb)

      DO i=1,3
        DO j=1,3
          GZU(i,j)=0.0d0
        ENDDO
      ENDDO

      IF(ng.EQ.0) THEN
C       Interpolate material parameters at XI

C Hari commented Jun17 2014
C       CALL CPXI(1,IBT,IDO,INP,NPNE,nr,nx,CE,CP,CW,XI,ERROR,*9999)
      ELSE
C       Put Gauss pt params into CW array
        DO il=1,ILT(1,nr,nx)
          CW(il)=CG(il,ng)
        ENDDO !il
      ENDIF

C-stress referred to Xj--------------------------------------------------

      IF(KTYP53(nr).EQ.1) THEN !stress ref'ed to Xj in constitutive law

C-stress referred to Nu--------------------------------------------------

      ELSE IF(KTYP53(nr).GT.1) THEN !stress ref'ed to Nu in constit law
        IF(ng.EQ.0) THEN
C ***     Interpolate midwall geometric var.s XG and derivs wrt Xi
          CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
        ELSE
C ***     Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
        ENDIF
C ***   Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C ***   derivatives of Xi wrt Xj (reference) coords, DXIXJ (IP=0)
        CALL XGMG(0,NITB,NBJ(1),nr,DXIXJ,GXL,GXU,RWX,XG,
     '    ERROR,*9999)
C ***   Calculate 2D Jacobian wrt undef coords for face integrals
C MPN 28Feb97
        RGX2D=RWX*DSQRT(GXU(3,3))
c        GX2D=GXL(1,1)*GXL(2,2)-GXL(1,2)*GXL(2,1)
c        RGX2D=DSQRT(GX2D*GXU(3,3))
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' RGX='',D12.4,'' RGX2D='',D12.4)')
     '      RGX,RGX2D
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
C ***   Get derivs of Xi wrt undeformed Nu (body/fibre) coords,DXIXN
        CALL DXIDXM(NITB,nr,DXIXN,DETERM,XG,'Fibre',ERROR,*9999)
        IF(COORDS.EQ.'Fibre'.OR.COORDS.EQ.'Principal'
     '     .OR.KTYP53(nr).EQ.3) THEN
          IF(ng.EQ.0) THEN
C ***       Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
            CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXJ,ZE,ZG,XI,
     '        ERROR,*9999)
          ELSE
C ***       Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
            CALL ZEZG(0,NBH,ng,NHE,nx,DXIXJ,PG,ZE,ZG,ERROR,*9999)
          ENDIF
C ***     Calculate deformed metric tensors wrt Xi (GZL,GZU)
          CALL ZGMG(nb,nr,GZ,GZL,GZU,ZG,ERROR,*9999)
          RGZ=DSQRT(GZ)
C ***     Calculate 2D Jacobian wrt def coords for face integrals
C MPN 28Feb97
          RGZ2D=DSQRT(GZ*GZU(3,3))
c          GZ2D=GZL(1,1)*GZL(2,2)-GZL(1,2)*GZL(2,1)
c          RGZ2D=DSQRT(GZ2D*GZU(3,3))
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' RGZ='',D12.4,'' RGZ2D='',D12.4)')
     '        RGZ,RGZ2D
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
C new MPN 4-May-96: new way of handling sheets
C         Get derivs of Xi wrt deformed Nu coords, DXIZN
          CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '      DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,'Fibre',ERROR,*9999)
C old
CC ***     Get derivs of Xi wrt deformed Nu coords, DXIZN
C          CALL DXIDXM(NBJ(1),nr,DXIZN,DZNXI,GZL,GZU,XG,'Fibre',ERROR,*9999)
C end old
C ***     Calculate derivs of deformed Nu wrt undeformed Nu (DZDX)
          DO ni=1,NITB
            DO mi=1,NITB
              SUM=0.0d0
              DO k=1,NITB
                SUM=SUM+DZNXI(ni,k)*DXIXN(k,mi)
              ENDDO !k
              DZDX(ni,mi)=SUM
            ENDDO !mi
          ENDDO !ni
          CALL INVERT(NITB,DZDX,DXDZ,DET_DZDX)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
              WRITE(OP_STRING,'('' DZDX('',I1,'',nj)   : '',3D12.4)')
     '          nhx,(DZDX(nhx,nj),nj=1,NJ_LOC(NJL_GEOM,0,nr))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nhx
CC$          call mp_unsetlock()
          ENDIF
        ENDIF !COORDS=Fibre or Princ
        IF(ng.EQ.0) THEN
C ***     Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
          CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXN,ZE,ZG,XI,
     '      ERROR,*9999)
        ELSE
C ***     Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
          IF(KTYP56(nr).EQ.6) THEN      !linear viscous relation
C           ..call to ZEZG with JP=3 computes 2nd derivs also
            CALL ZEZG(3,NBH,ng,NHE,nx,DXIXN,PG,ZE,ZG,ERROR,*9999)
          ELSE IF(KTYP56(nr).NE.6) THEN !all others
            CALL ZEZG(1,NBH,ng,NHE,nx,DXIXN,PG,ZE,ZG,ERROR,*9999)
          ENDIF !ktyp56
        ENDIF
C ***   Calculate deformed metric tensors wrt Nu (AZL,AZU)
        CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)
C ***   Get contravariant cpts of 2nd Piola-Kirchhoff stress
C ***   tensor (TG) wrt undeformed Nu coordinates
        IF(KTYP51(nr).EQ.1) THEN !plane stress
          CALL ZGTG51(nb,nr,nx,AXU,AZ,AZL,AZU,CW,
     '      RI1,RI2,RI3,TG,YG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.2) THEN !plane strain
          CALL ZGTG52(nb,nr,nx,AXU,AZ,AZL,AZU,CW,
     '      RI1,RI2,RI3,TG,YG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.3) THEN !3D
          CALL ZGTG53(STRESSTYPE,nb,nr,nx,AXU,AZ,AZL,AZU,CW,EG,
     '      RI1,RI2,RI3,TG,XG,YG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.4) THEN !membrane
          CALL ZGTG54(nb,nr,AXU,AZ,AZL,AZU,CW,EG,
     '      RI1,RI2,RI3,TG,YG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.5) THEN !string
          CALL ZGTG55(nr,AZL,CW,EG,TG,ERROR,*9999)
        ENDIF
C MPN 24May2000: option for total, or just passive or active stress cmpts
        IF(STRESSTYPE(1:6).EQ.'Active') THEN
C         Set TG to zero since just active components selected
          DO i=1,3
            DO j=1,3
              TG(i,j)=0.0d0
            ENDDO !j
          ENDDO !i
        ENDIF !STRESSTYPE = active
        IF(KTYP53(nr).EQ.3.AND.
     '    (STRESSTYPE(1:5).EQ.'Total'.OR.
     '     STRESSTYPE(1:14).EQ.'Total_no_hydro'.OR.
     '     STRESSTYPE(1:6).EQ.'Active')) THEN
C ***     Add active fibre stress component to TG
          FEXT(1)=DSQRT(AZL(1,1))
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
            ACTIVE_STRESS = YG(1)
          ENDIF
C          CALL ZGTG5A(nb,nr,FEXT,DXDZ,DZDX,DET_DZDX,TG,TNA,YG,ERROR,
C     '      *9999)   
          CALL ZGTG5A(nb,nr,FEXT,DXDZ,DZDX,DET_DZDX,TG,TNA,ACTIVE_STRESS
     &         ,ERROR,*9999)

        ENDIF !KTYP53=3 (active stresses) and (STRESSTYPE = total/active)
C newe MPN 24May2000

        IF(COORDS.EQ.'Reference') THEN
          IF(ng.EQ.0) THEN
C ***       Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
            CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXJ,ZE,ZG,XI,
     '        ERROR,*9999)
          ELSE
C ***       Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
            CALL ZEZG(1,NBH,ng,NHE,nx,DXIXJ,PG,ZE,ZG,ERROR,*9999)
          ENDIF
C ***     Put partial derivs of deformed theta wrt undef Xj into DZDX
          CALL DLZJDX(1,nb,nr,DZDX,XG,ZG,ERROR,*9999)
C ***     Calculate deformed metric tensors wrt Xj (AZL,AZU)
          CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)
C ***     Calculate derivs of undeformed Xj wrt undeformed Nu (DXJXN)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            DO mi=1,NITB
              SUM=0.0d0
              DO k=1,NITB
                SUM=SUM+DXIXN(k,mi)*XG(nj,NU1(k))
              ENDDO
              DXJXN(nj,mi)=SUM
            ENDDO
          ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,'('' DXJXN:'',10X,3D12.4,/(17X,3D12.4))')
     '        ((DXJXN(mz,nz),nz=1,3),mz=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !DOP
C ***     Transform TG to (undeformed) Xj coords
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            DO mjj=1,NJ_LOC(NJL_GEOM,0,nr)
              mj=NJ_LOC(NJL_GEOM,mjj,nr)
              SUM=0.0d0
              DO ni=1,NITB
                DO mi=1,NITB
                  SUM=SUM+DXJXN(nj,ni)*DXJXN(mj,mi)*TG(ni,mi)
                ENDDO !mi
              ENDDO !ni
              TGX(nj,mj)=SUM
            ENDDO !mjj (mj)
          ENDDO !njj (nj)
          DO ni=1,NITB
            DO mi=1,NITB
              TG(ni,mi)=TGX(ni,mi)
            ENDDO !mi
          ENDDO !ni
          IF(DOP) THEN
            WRITE(OP_STRING,'('' TG:'',13X,3D12.4,/(17X,3D12.4))')
     '        ((TG(mz,nz),nz=1,3),mz=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !DOP
        ENDIF !COORDS=Refer
      ENDIF !KTYP53(nr) (reference axes for stresses)

C-----------------------------------------------------------------------

      IF(KTYP51(nr).NE.5) THEN !plane stress,strain/3D/membrane/shell
        DO mz=1,NITB
          DO nz=1,NITB
            SUM=0.0d0
            DO mix=1,NITB
              DO nix=1,NITB
                SUM=SUM+DZDX(mz,mix)*TG(mix,nix)*DZDX(nz,nix)
              ENDDO !nix
            ENDDO !mix
            TC(mz,nz)=SUM/DSQRT(RI3)
          ENDDO !nz
        ENDDO !mz
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' TC:'',12X,3D12.4,/(16X,3D12.4))')
     '      ((TC(mz,nz),nz=1,3),mz=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF !DOP
        DO nix=1,NITB
          DO nz=1,NITB
            SUM=0.0d0
            DO mix=1,NITB
              SUM=SUM+TG(nix,mix)*DZDX(nz,mix)
            ENDDO !mix
            TN(nix,nz)=SUM
          ENDDO !nz
        ENDDO !nix
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' TN:'',12X,3D12.4,/(16X,3D12.4))')
     '      ((TN(mz,nz),nz=1,3),mz=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF !DOP

        IF(COORDS.EQ.'Reference') THEN
C ***     Compute undeformed metric tensor wrt theta coordinates AXU
C ***     and Eulerian metric tensor AZLZ (i.e. wrt deformed theta)
          DO i=1,3
            DO j=1,3
              AXU (i,j)=0.0d0
              AZLZ(i,j)=0.0d0
            ENDDO !j
            AXU (i,i)=1.0d0
            AZLZ(i,i)=1.0d0
          ENDDO !i
          IF(NITB.EQ.2) THEN
            IF(ITYP10(nr).EQ.2) THEN
              RR=XG(1,1)**2
              AXU(2,2)=1.0d0/RR
              AZLZ(2,2)=ZG(1,1)**2
            ELSE IF(ITYP10(nr).EQ.3) THEN
              RR=XG(1,1)**2
              AXU(2,2)=1.0d0/RR
              AXU(3,3)=1.0d0
              RR=ZG(1,1)**2
              AZLZ(3,3)=1.0d0
              AZLZ(2,2)=RR
            ELSE IF(ITYP10(nr).EQ.4) THEN
              AA=FOCUS*FOCUS
              SLX=DSINH(XG(1,1))
              SMX=DSIN (XG(2,1))
              G1=AA*(SLX*SLX+SMX*SMX)
              G3=AA* SLX*SLX*SMX*SMX
              AXU(1,1)=1.0d0/G1
              AXU(2,2)=1.0d0/G1
              AXU(3,3)=1.0d0
              SLX=DSINH(ZG(1,1))
              SMX=DSIN (ZG(2,1))
              G1=AA*(SLX*SLX+SMX*SMX)
              G3=AA* SLX*SLX*SMX*SMX
              AZLZ(1,1)=G1
              AZLZ(2,2)=G1
              AZLZ(3,3)=1.0d0
            ELSE IF(ITYP10(nr).EQ.5) THEN
            ENDIF
          ELSE IF(NITB.EQ.3) THEN
            IF(ITYP10(nr).EQ.2) THEN
              RR=XG(1,1)**2
              AXU(2,2)=1.0d0/RR
              AZLZ(2,2)=ZG(1,1)**2
            ELSE IF(ITYP10(nr).EQ.3) THEN
              RR=XG(1,1)**2
              RC=RR*DCOS(XG(3,1))**2
              AXU(2,2)=1.0d0/RC
              AXU(3,3)=1.0d0/RR
              RR=ZG(1,1)**2
              AZLZ(3,3)=RR
              AZLZ(2,2)=RR*DCOS(ZG(3,1))**2
            ELSE IF(ITYP10(nr).EQ.4) THEN
              AA=FOCUS*FOCUS
              SLX=DSINH(XG(1,1))
              SMX=DSIN (XG(2,1))
              G1=AA*(SLX*SLX+SMX*SMX)
              G3=AA* SLX*SLX*SMX*SMX
              AXU(1,1)=1.0d0/G1
              AXU(2,2)=1.0d0/G1
              AXU(3,3)=1.0d0/G3
              SLX=DSINH(ZG(1,1))
              SMX=DSIN (ZG(2,1))
              G1=AA*(SLX*SLX+SMX*SMX)
              G3=AA* SLX*SLX*SMX*SMX
              AZLZ(1,1)=G1
              AZLZ(2,2)=G1
              AZLZ(3,3)=G3
            ELSE IF(ITYP10(nr).EQ.5) THEN
            ENDIF
          ENDIF
C ***     Compute physical cpts of Cauchy & Nominal stresses
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            DO mjj=1,NJ_LOC(NJL_GEOM,0,nr)
              mj=NJ_LOC(NJL_GEOM,mjj,nr)
              TC(nj,mj)=TC(nj,mj)*DSQRT(AZLZ(nj,nj)*AZLZ(mj,mj))
              TN(nj,mj)=TN(nj,mj)*DSQRT(AZLZ(mj,mj)/ AXU(nj,nj))
            ENDDO !mjj (mj)
          ENDDO !njj (nj)
        ENDIF !COORDS=Refer
C ***   Compute principal Stresses, Euler angles and eigenvectors
        IF(CSTYPE.EQ.'Piola') THEN
          DO ni=1,NITB
            DO mi=1,NITB
              TGX(ni,mi)=TG(ni,mi)
            ENDDO !mi
          ENDDO !ni

C          CALL F02ABF(TGX,3,NITB,PST,RM,3,WK1_LOCAL,IFAIL)
C          DO ni=1,3
C            DO mi=1,3
C              RM(ni,mi)=TGX(ni,mi)
C            ENDDO
C          ENDDO
C MLB 19/3/97
C This may not give evectors as accurately as NAG
C          CALL DSYEV('V','L',NITB,RM,3,PST,WK1_LOCAL,10,IFAIL)

          CALL ESOLVE(.TRUE.,NITB,TGX,3,PST,RM,ERROR,*9999)

        ELSE IF(CSTYPE.EQ.'Nominal') THEN

C Hari deleted Jun 17, 2014

        ELSE IF(CSTYPE.EQ.'Cauchy') THEN
! Hari commented Oct 28, 2014
!         CALL ESOLVE(.TRUE.,NITB,TC,3,PST,RM,ERROR,*9999)

        ENDIF
C news MPN 30-Jun-94
C ***   Order PST from max->min and change cols of RM accordingly
        IF(PST(2).GT.PST(1)) THEN
          VALTMP=PST(1)
          PST(1)=PST(2)
          PST(2)=VALTMP
          DO ni=1,NITB
            VALTMP=RM(ni,1)
            RM(ni,1)=RM(ni,2)
            RM(ni,2)=VALTMP
          ENDDO
        ENDIF
        IF(NITB.EQ.3) THEN
          DO mi=1,2
            IF(PST(3).GT.PST(MI)) THEN
              VALTMP=PST(MI)
              PST(MI)=PST(3)
              PST(3)=VALTMP
              DO ni=1,NITB
                VALTMP=RM(ni,mi)
                RM(ni,mi)=RM(ni,3)
                RM(ni,3)=VALTMP
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C newe
        IF(NITB.EQ.2) THEN
          IF(DABS(RM(1,1)).LE.1.0d0) THEN
            PHI(1)=DACOS(RM(1,1))
          ELSE
            PHI(1)=0.0d0
          ENDIF
          IF(DABS(RM(1,2)).LE.1.0d0) THEN
            PHI(2)=DACOS(RM(1,2))
          ELSE
            PHI(2)=0.0d0
          ENDIF
        ELSE IF(NITB.EQ.3) THEN
          TOL=1.0d-08
          IF(DABS(RM(1,1)).GT.TOL) THEN
            PHI(1)=DATAN2(RM(2,1),RM(1,1))
          ELSE IF(RM(2,1)*RM(1,1).GT.0.0d0) THEN
            PHI(1)=90.0d0
          ELSE
            PHI(1)=-90.0d0
          ENDIF
          IF(DABS(RM(3,1)).LE.1.0d0) THEN
            PHI(2)=DASIN(RM(3,1))
          ELSE IF(RM(3,1).GT.1.0d0) THEN
            PHI(2)=90.0d0
          ELSE
            PHI(2)=-90.0d0
          ENDIF
           IF(DABS(DCOS(PHI(1))).GT.TOL .AND.
     '       DABS(DCOS(PHI(1))).GE.DABS(RM(3,3))) THEN
             PHI(3)=DACOS(RM(3,3)/DCOS(PHI(1)))
          ELSE
            PHI(3)=0.0d0
          ENDIF
        ENDIF
        DO ni=1,NITB
          PHI(ni)=PHI(ni)*180.0d0/PI
          IF(PHI(ni).GT.90.0d0) PHI(ni)=PHI(ni)-180.0d0
          IF(PHI(ni).LT.-90.0d0) PHI(ni)=PHI(ni)+180.0d0
C         Find Max/Min principal stresses
          IF(PST(ni).GT.PRSTMAX) PRSTMAX=PST(ni)
          IF(PST(ni).LT.PRSTMIN) PRSTMIN=PST(ni)
        ENDDO
      ENDIF !KTYP51(nr)=plane stress,strain/3D/membrane/shell

      CALL EXITS('ZETX50')
      RETURN
 9999 CALL ERRORS('ZETX50',ERROR)
      CALL EXITS('ZETX50')
      RETURN 1
      END


