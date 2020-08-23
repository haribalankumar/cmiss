      SUBROUTINE ZEEX50(COORDS,IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,
     '  NPNE,nr,nx,DZDX,CE,CG,CP,EG,PG,PHI,PST,R,RGX,RI1,RI2,RI3,
     '  RM,U,XE,XG,XI,ZE,ZG,ERROR,*)

C#### Subroutine: ZEEX50
C###  Description:
C###    ZEEX50 calculates normal, shear and principal Greens strains
C###    EG with respect to 'Reference', 'Fibre' or 'Wall' coords
C###    (as specified by COORDS) at position XI in current element.

C**** Since base vectors of theta coords are not orthonormal, cpts
C**** of Greens strain are converted to 'physical' values.
C**** AZ,AZL,AZU  are deformed metric tensors wrt (undeformed) COORDS
C**** RI1,RI2,RI3 are principal invariants of AZL
C**** AXU  are contravariant cpts of undeformed metric tensor
C**** D    is the diagonal matrix of principal stretches
C**** DZDX are components of the deformation gradient tensor
C**** XG   are undeformed theta coords and derivs wrt Xi
C**** ZG   are deformed theta coords and derivs wrt undeformed COORDS
C**** EG   are physical cpts of Green's strain
C**** EXR  are extension ratios wrt COORDS
C**** PHI  are the Euler angles wrt COORDS of the principal extensions
C**** PST  are principal strains
C**** R    is the orthogonal rotation tensor
C**** RM   is the modal matrix whose cols are the eigenvectors
C**** U    is the right stretch tensor
C**** UI   is the inverse of U

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ng,NHE,
     '  NPNE(NNM,NBFM),nr,nx
      REAL*8 DZDX(3,3),CE(NMM),CG(NMM,NGM),CP(NMM,NPM),
     '  EG(3,3),PG(NSM,NUM,NGM,NBM),PHI(3),PST(3),
     '  R(3,3),RGX,RI1,RI2,RI3,RM(3,3),U(3,3),
     '  XE(NSM,NJM),XG(NJM,NUM),XI(NIM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER COORDS*(*),ERROR*(*)
!     Local Variables
      INTEGER i,IFAIL,il,j,k,m,mi,mj,mjj,mix,
     '  n,nb,NCW,ni,NITB,nj,njj,nix,nz
      PARAMETER (NCW=35) !CW must be dimen.d the same size as CE array
      REAL*8 AA,AXU(3,3),AZ,AZL(3,3),AZL_tmp(3,3),AZLZ(3,3),AZU(3,3),
     '  CW(NCW),D(3,3),DXIXJ(3,3),DXIXN(3,3),DXIZN(3,3),
     '  DETERM,DZDX_tmp(3,3),
     '  DZNXI(3,3),Fgrowth(3,3),G1,G3,GXL(3,3),GXU(3,3),RC,RR,SLX,SMX,
     '  SUM,SUM1,SUM2,TOL,UI(3,3),VALTMP,WK1_LOCAL(10)

      DATA TOL /1.0d-12/

      CALL ENTERS('ZEEX50',*9999)
      CALL ASSERT(NCW.EQ.NMM,'>>Dimension of CW array (NCW)'
     '  //' must equal dimension of CE array (NMM)',ERROR,*9999)
      nb=NBH(NH_LOC(1,nx))
      NITB=NIT(nb)

      DO i=1,3
        DO j=1,3
          AZL(i,j)=0.0d0
        ENDDO
      ENDDO

      IF(ng.EQ.0) THEN
C ***   Interpolate Xi-coord geometric var.s XG and derivs wrt Xi
        CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
      ELSE
C ***   Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
      ENDIF
C *** Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C *** ..derivatives of Xi wrt Xj (reference) coords, DXIXJ (IP=0)
      CALL XGMG(0,NITB,NBJ(1),nr,DXIXJ,GXL,GXU,RGX,XG,ERROR,*9999)

      IF(COORDS(1:5).EQ.'Refer'.OR.
     '  (COORDS(1:5).EQ.'Princ'.AND.NITB.EQ.3)) THEN
        IF(ng.EQ.0) THEN
C ***     Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
          CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '      DXIXJ,ZE,ZG,XI,ERROR,*9999)
        ELSE
C ***     Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
          CALL ZEZG(1,NBH,ng,NHE,nx,DXIXJ,PG,ZE,ZG,ERROR,*9999)
        ENDIF
C ***   Put partial derivs of deformed theta wrt undef Xj into DZDX
        CALL DLZJDX(1,nb,nr,DZDX,XG,ZG,ERROR,*9999)
C ***   Calculate deformed metric tensors wrt Xj (AZL,AZU)
        CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)

      ELSE IF(COORDS(1:5).EQ.'Fibre'.OR.COORDS(1:4).EQ.'Wall'.OR.
     '  (COORDS(1:5).EQ.'Princ'.AND.NITB.EQ.2)) THEN
C new VYW/MPN 2010-03-16
C ***   Get derivs of Xi wrt undeformed Wall/Nu(body/fibre) coords (DXIXN)
        CALL DXIDXM(NITB,nr,DXIXN,DETERM,XG,COORDS,ERROR,*9999)
C ***   Get derivs of Xi wrt deformed Wall/Nu coords (DXIZN)
        CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '    DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,COORDS,ERROR,*9999)
C old VYW/MPN 2010-03-16
CC news MPN 31Jul97
C        IF(COORDS(1:4).EQ.'Wall') THEN
CC ***     Get derivs of Xi wrt undeformed Wall (cardiac) coords (DXIXN)
C          CALL DXIDXM(NITB,nr,DXIXN,DETERM,XG,'Wall',ERROR,*9999)
CC ***     Get derivs of Xi wrt deformed Wall coords (DXIXN)
C          CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
C     '      DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,'Wall',ERROR,*9999)
C        ELSE
CC newe
CC ***     Get derivs of Xi wrt undeformed Nu (body/fibre) coords (DXIXN)
C          CALL DXIDXM(NITB,nr,DXIXN,DETERM,XG,'Fibre',ERROR,*9999)
CC ***     Get derivs of Xi wrt deformed Nu coords (DXIZN)
C          CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
C     '      DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,'Fibre',ERROR,*9999)
C        ENDIF
C old VYW/MPN
C ***   Calc deformation gradients (DZDX) wrt Nu/Wall coords
        DO ni=1,NITB
          DO mi=1,NITB
            SUM=0.0D0
            DO k=1,NITB
              SUM=SUM+DZNXI(ni,k)*DXIXN(k,mi)
            ENDDO
            DZDX(ni,mi)=SUM
          ENDDO
        ENDDO
        IF(ng.EQ.0) THEN
C ***     Interpolate dependent var.s ZG and derivs wrt Nu/Wall (JP=1)
          CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '      DXIXN,ZE,ZG,XI,ERROR,*9999)
        ELSE
C ***     Interpolate dependent var.s ZG and derivs wrt Nu/Wall (JP=1)
          IF(KTYP56(nr).EQ.6) THEN      !linear viscous relation
            CALL ZEZG(3,NBH,ng,NHE,nx,DXIXN,PG,ZE,ZG,ERROR,*9999)
          ELSE! IF(KTYP56(nr).NE.6) THEN !all others
            CALL ZEZG(1,NBH,ng,NHE,nx,DXIXN,PG,ZE,ZG,ERROR,*9999)
          ENDIF
        ENDIF !ng
C ***   Calculate deformed metric tensors wrt Nu/Wall (AZL,AZU)
        CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)

C old
C        IF(COORDS(1:4).EQ.'Wall') THEN
CC         Calc def anatomical fibre vects wrt rc coords at XI
C          CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
C     '      DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),PG,XE,XG,XI,ZE,ZG,
C     '      ERROR,*9999)
C          CALL INVERT(NITB,DZDNU,DNUDZ,DETERM)
CC         Calc deformed wall vectors wrt rc coords at XI
C          CALL WALL_VEC_DEF(IBT,IDO,INP,NAN,NBH,0,NHE,NITB,nr,nx,
C     '      DZDW(1,1),DZDW(1,2),DZDW(1,3),PG,XI,ZE,ZG,ERROR,*9999)
CC         Calc derivs of deformed wall vectors wrt nu coords
C          DO mz=1,3
C            DO nz=1,3
C              SUM=0.0d0
C              DO mi=1,3
C                SUM=SUM+DNUDZ(mz,mi)*DZDW(mi,nz)
C              ENDDO !mi
C              DNUDW(mz,nz)=SUM
C            ENDDO !nz
C          ENDDO !mz
C          CALL INVERT(NITB,DNUDW,DWDNU,DETERM)
CC         Compute metric and deformation tensors wrt wall coords
CC         (transform from fibre cords to wall coords using DWDNU)
C          DO i=1,3
C            DO j=1,3
C              AZL_tmp(i,j)=AZL(i,j)
C              DZDX_tmp(i,j)=DZDX(i,j)
C            ENDDO !j
C          ENDDO !i
C          DO i=1,3
C            DO j=1,3
C              AZL(i,j)=0.0d0
C              DZDX(i,j)=0.0d0
C              DO m=1,3
C                DO n=1,3
C                  AZL(i,j)=AZL(i,j)+
C     '              AZL_tmp(m,n)*DNUDW(m,i)*DNUDW(n,j)
C                  DZDX(i,j)=DZDX(i,j)+
C     '              DWDNU(i,m)*DZDX_tmp(m,n)*DNUDW(n,j)
CC!!! WARNING: DZDX needs checking!!!
C                ENDDO !n
C              ENDDO !m
C            ENDDO !j
C          ENDDO !i
CC         Recompute AZU,AZ from transformed AZL
C          CALL INVERT(NITB,AZL,AZU,AZ)
C          IF(DOP) THEN
CC$          call mp_setlock()
C            WRITE(OP_STRING,'(''  Wall AZL:'',3D12.4,'
C     '        //'''  Wall AZU:'',3D12.4,/(11X,3D12.4,11X,3D12.4))')
C     '        ((AZL(mi,ni),ni=1,NITB),
C     '        (AZU(mi,ni),ni=1,NITB),mi=1,NITB)
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
C          ENDIF !DOP
C        ENDIF !coords='Wall'
      ENDIF

      DO i=1,3
        DO j=1,3
          AXU(i,j) =0.0d0
          AZLZ(i,j)=0.0d0
          D(i,j)   =0.0d0
        ENDDO !j
        AXU(i,i) =1.0d0
        AZLZ(i,i)=1.0d0
      ENDDO !i

      IF(NITB.EQ.1) THEN
      ELSE IF(NITB.EQ.2) THEN
        IF(COORDS(1:5).EQ.'Refer') THEN
          IF(ITYP10(nr).EQ.1) THEN
C ***       Rectangular Cartesian Coordinates
            RI1= AZL(1,1)+AZL(2,2)
            RI2= AZ
          ELSE IF(ITYP10(nr).EQ.2) THEN
C ***       Cylindrical Polar Coordinates
            AZLZ(2,2)=ZG(1,1)**2
            RR=XG(1,1)**2
            RI1= AZL(1,1)+AZL(2,2)/RR
            RI2=AZ/RR
            AXU(2,2)=1.0D0/RR
          ELSE IF(ITYP10(nr).EQ.3) THEN
C ***       Spherical Polar Coordinates
            AZLZ(2,2)=ZG(1,1)**2
            RR=XG(1,1)**2
            RI1= AZL(1,1)+AZL(2,2)/RR
            RI2=AZ/RR
            AXU(2,2)=1.0D0/RR
          ELSE IF(ITYP10(nr).EQ.4) THEN
C ***       Prolate Spheroidal Coordinates
            AA=FOCUS*FOCUS
            SLX=DSINH(ZG(1,1))
            SMX=DSIN (ZG(2,1))
            G1=AA*(SLX*SLX+SMX*SMX)
            G3=AA* SLX*SLX*SMX*SMX
            AZLZ(1,1)=G1
            AZLZ(2,2)=G1
            SLX=DSINH(XG(1,1))
            SMX=DSIN (XG(2,1))
            G1=AA*(SLX*SLX+SMX*SMX)
            G3=AA* SLX*SLX*SMX*SMX
            RI1= (AZL(1,1)+AZL(2,2))/G1
            RI2=AZ/(G1*G1)
            AXU(1,1)=1.0D0/G1
            AXU(2,2)=1.0D0/G1
          ELSE IF(ITYP10(nr).EQ.5) THEN
C ***       Oblate Spheroidal Coordinates
          ENDIF
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            DO mjj=1,NJ_LOC(NJL_GEOM,0,nr)
              mj=NJ_LOC(NJL_GEOM,mjj,nr)
              DZDX(nj,mj)=DZDX(nj,mj)*DSQRT(AZLZ(nj,nj)*AXU(mj,mj))
            ENDDO
          ENDDO
        ELSE IF(COORDS(1:5).EQ.'Fibre'.OR.COORDS(1:4).EQ.'Wall'.OR.
     '    COORDS(1:5).EQ.'Princ') THEN
          RI1= AZL(1,1)+AZL(2,2)
          RI2= AZ
        ENDIF

      ELSE IF(NITB.EQ.3) THEN
        IF(COORDS(1:5).EQ.'Refer'.OR.
     '    COORDS(1:5).EQ.'Princ') THEN
          IF(ITYP10(nr).EQ.1) THEN
C ***       Rectangular Cartesian Coordinates
            RI3= AZ
            RI1= AZL(1,1)+AZL(2,2)+AZL(3,3)
            RI2=(AZU(1,1)+AZU(2,2)+AZU(3,3))*RI3
          ELSE IF(ITYP10(nr).EQ.2) THEN
C ***       Cylindrical Polar Coordinates
            AZLZ(2,2)=ZG(1,1)**2
            RR=XG(1,1)**2
            RI3=AZ/RR
            RI1= AZL(1,1)+AZL(2,2)/RR+AZL(3,3)
            RI2=(AZU(1,1)+AZU(2,2)*RR+AZU(3,3))*RI3
            AXU(2,2)=1.0D0/RR
          ELSE IF(ITYP10(nr).EQ.3) THEN
C ***       Spherical Polar Coordinates
            AZLZ(3,3)=ZG(1,1)**2
            AZLZ(2,2)=ZG(1,1)**2*DCOS(ZG(3,1))**2
            RR=XG(1,1)**2
            RC=RR*DCOS(XG(3,1))**2
            RI3=AZ/(RR*RC)
            RI1= AZL(1,1)+AZL(2,2)/RC+AZL(3,3)/RR
            RI2=(AZU(1,1)+AZU(2,2)*RC+AZU(3,3)*RR)*RI3
            AXU(2,2)=1.0D0/RC
            AXU(3,3)=1.0D0/RR
          ELSE IF(ITYP10(nr).EQ.4) THEN
C ***       Prolate Spheroidal Coordinates
            AA=FOCUS*FOCUS
            SLX=DSINH(ZG(1,1))
            SMX=DSIN (ZG(2,1))
            G1=AA*(SLX*SLX+SMX*SMX)
            G3=AA* SLX*SLX*SMX*SMX
            AZLZ(1,1)=G1
            AZLZ(2,2)=G1
            AZLZ(3,3)=G3
            SLX=DSINH(XG(1,1))
            SMX=DSIN (XG(2,1))
            G1=AA*(SLX*SLX+SMX*SMX)
            G3=AA* SLX*SLX*SMX*SMX
            RI3=AZ/(G1*G1*G3)
            RI1= (AZL(1,1)+AZL(2,2))/G1+AZL(3,3)/G3
            RI2=((AZU(1,1)+AZU(2,2))*G1+AZU(3,3)*G3)*RI3
            AXU(1,1)=1.0D0/G1
            AXU(2,2)=1.0D0/G1
            AXU(3,3)=1.0D0/G3
          ELSE IF(ITYP10(nr).EQ.5) THEN
C ***       Oblate Spheroidal Coordinates
          ENDIF
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            DO mjj=1,NJ_LOC(NJL_GEOM,0,nr)
              mj=NJ_LOC(NJL_GEOM,mjj,nr)
              DZDX(nj,mj)=DZDX(nj,mj)*DSQRT(AZLZ(nj,nj)*AXU(mj,mj))
            ENDDO
          ENDDO
        ELSE IF(COORDS(1:5).EQ.'Fibre'.OR.COORDS(1:4).EQ.'Wall') THEN
          RI3= AZ
          SUM1=0.0D0
          SUM2=0.0D0
          DO ni=1,NITB
            SUM1=SUM1+AZL(ni,ni)
            SUM2=SUM2+AZU(ni,ni)
          ENDDO
          RI1=SUM1
          RI2=SUM2*RI3
        ENDIF
      ENDIF

C!!! NOTE can only use resid strains for pole zero law until init extn
C!!!      mat params are set up for other problem types
      IF(KTYP55(nr).EQ.3.AND.KTYP56(nr).EQ.3) THEN !pole zero law

        IF(ng.EQ.0) THEN
C         Interpolate material parameters at XI
          CALL CPXI(1,IBT,IDO,INP,NPNE,nr,nx,CE,CP,CW,XI,ERROR,*9999)
        ELSE
C         Put Gauss pt params into CW array
          DO il=1,ILT(1,nr,nx)
            CW(il)=CG(il,ng)
          ENDDO !il
        ENDIF

C       Calc growth defm tens for resid strain and copy AZL to AZL_tmp
        DO i=1,3
          DO j=1,3
            Fgrowth(i,j)=0.0d0 !growth defm tensor for resid strains
            AZL_tmp(i,j)=AZL(i,j)
            DZDX_tmp(i,j)=DZDX(i,j)
          ENDDO !j
          Fgrowth(i,i)=CW(27+i) !NOTE init extns CG(28),CG(29),CG(30)
        ENDDO !i
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' Fgrowth: '',3D12.4,/(10X,3D12.4))')
     '      ((Fgrowth(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF !DOP
C       Apply growth defm to deformed covariant metric tensor AZL
        DO i=1,3
          DO j=1,3
            AZL(i,j)=0.0d0
            DZDX(i,j)=0.0d0
            DO m=1,3
              DO n=1,3
                AZL(i,j)=AZL(i,j)+
     '            Fgrowth(i,m)*AZL_tmp(m,n)*Fgrowth(n,j)
              ENDDO !n
              DZDX(i,j)=DZDX(i,j)+DZDX_tmp(i,m)*Fgrowth(m,j)
            ENDDO !m
          ENDDO !j
        ENDDO !i
C       Recompute AZU,AZ from transformed AZL
        CALL INVERT(NITB,AZL,AZU,AZ)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(''  "Grown" AZL:'',3D12.4,'
     '      //'''  "Grown" AZU:'',3D12.4,/(14X,3D12.4,14X,3D12.4))')
     '      ((AZL(mi,ni),ni=1,NITB),
     '      (AZU(mi,ni),ni=1,NITB),mi=1,NITB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF !DOP
      ENDIF !pole zero law

      IF(NITB.EQ.1) THEN
        EG(1,1)=0.5d0*(AZL(1,1)-1.0d0)
      ELSE IF(NITB.GE.2) THEN
        DO mi=1,NITB
          DO ni=1,NITB
C ***       Calculate physical shear cpts of Green's strain
            EG(mi,ni)=0.5d0*AZL(mi,ni)*DSQRT(AXU(mi,mi)*AXU(ni,ni))
          ENDDO
C ***     Calculate diagonal physical cpts of Green's strain
          EG(mi,mi)=0.5d0*(AZL(mi,mi)*AXU(mi,mi)-1.0d0)
        ENDDO !mi

        IFAIL=0
C        CALL F02ABF(EG,3,NITB,PST,RM,3,WK1_LOCAL,IFAIL)
        DO ni=1,3
          DO mi=1,3
            RM(ni,mi)=EG(ni,mi)
          ENDDO
        ENDDO
C MLB 19/3/97
C This may not give evectors as accurately as NAG
        CALL DSYEV('V','L',NITB,RM,3,PST,WK1_LOCAL,10,IFAIL)
!news MPN 30-Jun-94
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
!newe
        DO ni=1,NITB
          D(ni,ni)=DSQRT(2.0D0*PST(ni)+1.0D0)
        ENDDO
        DO nix=1,NITB
          DO mix=1,NITB
            SUM1=0.0D0
            SUM2=0.0D0
            DO ni=1,NITB
              SUM1=SUM1+RM(nix,ni)*D(ni,ni)*RM(mix,ni)
              SUM2=SUM2+RM(nix,ni)/D(ni,ni)*RM(mix,ni)
            ENDDO
            U (nix,mix)=SUM1
            UI(nix,mix)=SUM2
          ENDDO
        ENDDO
        DO nz=1,NITB
          DO nix=1,NITB
            SUM=0.0D0
            DO mix=1,NITB
              SUM=SUM+DZDX(nz,mix)*UI(mix,nix)
            ENDDO
            R(nz,nix)=SUM
          ENDDO
        ENDDO
      ENDIF

      IF(NITB.EQ.1) THEN
      ELSE IF(NITB.EQ.2) THEN
        IF(DABS(RM(1,1)).GT.TOL.OR.DABS(RM(2,1)).GT.TOL) THEN
          PHI(1)=DATAN2(RM(2,1),RM(1,1))
        ELSE
          PHI(1)=0.0D0
        ENDIF
        IF(DABS(RM(1,2)).GT.TOL.OR.DABS(RM(2,2)).GT.TOL) THEN
          PHI(2)=DATAN2(RM(2,2),RM(1,2))
        ELSE
          PHI(2)=0.0D0
        ENDIF
      ELSE IF(NITB.EQ.3) THEN
        IF(DABS(RM(1,1)).GT.TOL.OR.DABS(RM(2,1)).GT.TOL) THEN
          PHI(1)=DATAN2(RM(2,1),RM(1,1))
        ELSE IF(RM(2,1)*RM(1,1).GT.0.0D0) THEN
          PHI(1)=90.0D0
        ELSE
          PHI(1)=-90.0D0
        ENDIF
        IF(DABS(RM(3,1)).LE.1.0D0) THEN
          PHI(2)=DASIN(RM(3,1))
        ELSE IF(RM(3,1).GT.1.0D0) THEN
          PHI(2)=90.0D0
        ELSE
          PHI(2)=-90.0D0
        ENDIF
         IF(DABS(DCOS(PHI(1))).GT.TOL .AND.
     '     DABS(DCOS(PHI(1))).GE.DABS(RM(3,3))) THEN
           PHI(3)=DACOS(RM(3,3)/DCOS(PHI(1)))
        ELSE
          PHI(3)=0.0D0
        ENDIF
      ENDIF

      DO ni=1,NITB
        PHI(ni)=PHI(ni)*180.0D0/PI
        IF(PHI(ni).GT.90.0D0) THEN
          PHI(ni)=PHI(ni)-180.0D0
        ELSE IF(PHI(ni).LT.-90.0D0) THEN
          PHI(ni)=PHI(ni)+180.0D0
        ENDIF
C       Find Max/Min principal strains
        IF(PST(ni).GT.PRSTMAX) PRSTMAX=PST(ni)
        IF(PST(ni).LT.PRSTMIN) PRSTMIN=PST(ni)
      ENDDO

      CALL EXITS('ZEEX50')
      RETURN
 9999 CALL ERRORS('ZEEX50',ERROR)
      CALL EXITS('ZEEX50')
      RETURN 1
      END


