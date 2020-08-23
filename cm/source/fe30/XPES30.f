      SUBROUTINE XPES30(IBT,IDO,INP,NBH,NBJ,ne,NHE,NORD1,NPNE,nr,
     '  nx,CE,CG,CGE,CP,ED,EM,ER,ES,PG,RG,SE,WG,XE,XG,YG,ZE,ZG,
     '  UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*)

C#### Subroutine: XPES30
C###  Description:
C###    XPES30 calculates element matrices for linear static or
C###    time-dependent equations.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'equation00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     &  NBH(NHM,NCM),NBJ(NJM),ne,NHE,NORD1,NPNE(NNM,NBFM),nr,nx
      REAL*8 CE(NMM),CG(NMM,NGM),CGE(NMM,NGM),CP(NMM,NPM),
     &  ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     &  ES(NHM*NSM,NHM*NSM),PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     &  WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM),YG(NIYGM,NGM),ZE(NSM,NHM),
     &  ZG(NHM,NUM)
      CHARACTER ERROR*(*)
      LOGICAL UPDATE_MATRIX,UPDATE_VECTOR
!     Local Variables
      INTEGER CGFLAG,k,mh,mhs,mi,ms,nb,nbn,nc,ng,
     &  NGTB,nh,nhs,NHST,ni,NITB,nj,NJ_LOOP,
     &  no_tip,ns,NSTB,NSTBM,NSTBN,NU1(0:3),r,s,t
      REAL*8 AREA,DETERM,DXNXI(3,3),DSDXI,DXIXN(3,3),DXIX(3,3),
     &  GL(3,3),GU(3,3),
     &  PA,PB,PGM,PGMS,PGMSI(3),PGMX,PGN,PGNS,PGNSI(3),
     &  PGNX,PGX,PSI1,RK2,rm,RWG,S2,S3,S4,SCALE,SOURCE,SS,SUM,SUM1,
     &  SUM2,V,VELOC,DIFF,gam,Qc,cb,Qb,Vol_acinii
      CHARACTER COMMENT*64
      
      DATA NU1/1,2,4,7/

      CALL ENTERS('XPES30',*9999)
C GBS 10-NOV-1994
      nc=1 !Temporary

      IF(UPDATE_MATRIX.OR.UPDATE_VECTOR) THEN
        nb=NBH(NH_LOC(1,nx),nc)
        NGTB=NGT(nb)
        NITB=NIT(nb)
        NSTB=NST(nb)
        NHST=NHE*NSTB

        CALL CPCG(1,nb,NPNE,nr,nx,CE,CG,CGE,CP,PG,ERROR,*9999)

        DO ms=1,NHST
          IF(UPDATE_MATRIX) THEN
            DO ns=1,NHST
              EM(ms,ns)=0.0d0
              ED(ms,ns)=0.0d0
              ES(ms,ns)=0.0d0
            ENDDO !ns
          ENDIF !update
          IF(UPDATE_VECTOR) ER(ms)=0.0d0
        ENDDO !ms

C       Gauss point loop

        DO ng=1,NGTB
          
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
          
          !static div(kgrad(u))=f or time-dep advection-diffusion
          IF (((ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4).AND.
     '      ITYP2(nr,nx).EQ.5).OR.
     '      (ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.3))THEN
C           DXIX is wrt Nu-coords (PJH 6Mar96)
            CALL XGMG(1,NITB,nb,nr,DXIX,GL,GU,RG(ng),XG,ERROR,*9999)
          ELSE
C           DXIX is wrt Xj-coords
            CALL XGMG(0,NITB,nb,nr,DXIX,GL,GU,RG(ng),XG,ERROR,*9999)
          ENDIF
          IF(JTYP4.EQ.1)THEN !unsymmetric
            RWG=RG(ng)*WG(ng,nb)
          ELSE IF(JTYP4.EQ.2)THEN !cylindrical symmetry about x-axis
            RWG=RG(ng)*WG(ng,nb)*2.d0*PI*XG(2,1)
          ELSE IF(JTYP4.EQ.3)THEN !cylindrical symmetry about y-axis
            RWG=RG(ng)*WG(ng,nb)*2.d0*PI*XG(1,1)
          ELSE IF(JTYP4.EQ.4)THEN !spherical symmetry
            RWG=RG(ng)*WG(ng,nb)
          ENDIF
          IF(ITYP2(nr,nx).EQ.3.AND.ITYP3(nr,nx).EQ.2.AND.
     '      NJ_LOC(NJL_FIBR,0,nr).GT.0)THEN
            CALL DXIDXM(NITB,nr,DXIXN,DETERM,XG,'Fibre',ERROR,*9999)
            CALL INVERT(NITB,DXIXN,DXNXI,DETERM)
          ENDIF
          
          IF(USE_EQUATION_OBJECT)THEN
            IF((ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.6)
     &        .OR.(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.3))THEN
            
              WRITE(COMMENT,'(''Element:'',I6,'' Gauss Point:'',I2)') 
     &          ne,ng
              CALL SET_MATHS_ALL_UNSOLVED(ERROR,*9999)
              CALL UPDATE_EQUATION_INPUTS(nx,CG(1,ng),XG,DXIX,COMMENT,
     &          ERROR,*9999)
     
            ENDIF
          ENDIF
          
C MHT 28-5-98 Implementing loops over nh
          mhs=0
          DO mh=1,NHE
C         Loop over element rows
            NSTBM=NST(NBH(NH_LOC(mh,nx),nc))
            DO ms=1,NSTBM
              mhs=mhs+1
              nhs=0
              PGM=PG(ms,1,ng,nb)
              IF(UPDATE_MATRIX)THEN
C             Loop over element cols
                DO nh=1,NHE
                  NSTBN=NST(NBH(NH_LOC(nh,nx),nc))
                  DO ns=1,NSTBN
                    nhs=nhs+1
                    PGN=PG(ns,1,ng,nb)
                    IF(ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4)THEN
C                ********** Static or quasi-static solution *********
                      IF(ITYP2(nr,nx).EQ.3)THEN
C                     Laplace equation
                        IF(ITYP3(nr,nx).EQ.2.AND.
     '                    NJ_LOC(njl_fibr,0,nr).GT.0)THEN
C                       Generalised Laplace with anisotropic cond.ty
                          DO k=1,NITB
                            PGMSI(k)=PG(ms,NU1(k),ng,nb)
                            SUM1=0.0d0
                            DO s=1,NITB
                              SUM1=SUM1+PG(ns,NU1(s),ng,nb)*DXIXN(s,k)
                            ENDDO !s
                            PGNSI(k)=SUM1*CG(k,ng)
                          ENDDO !k
                          SUM=0.0d0
                          DO r=1,NITB
                            SUM1=0.0d0
                            DO k=1,NITB
                              SUM1=SUM1+PGNSI(k)*DXNXI(k,r)
                            ENDDO !k
                            SUM2=0.0d0
                            DO t=1,NITB
                              SUM2=SUM2+PGMSI(t)*GU(r,t)
                            ENDDO !t
                            SUM=SUM+SUM1*SUM2
                          ENDDO !r
                          ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                        ELSE
C                       Laplace or generalised Laplace with isotropic
C                       conductivity
                          DO ni=1,NITB
                            PGMSI(ni)=PG(ms,NU1(ni),ng,nb)
                            PGNSI(ni)=PG(ns,NU1(ni),ng,nb)
                          ENDDO
                          SUM =0.0d0
                          DO mi=1,NITB
                            DO ni=1,NITB
                              SUM=SUM+PGMSI(mi)*PGNSI(ni)*GU(mi,ni)
                            ENDDO
                          ENDDO
                          IF(ITYP3(nr,nx).EQ.2)THEN
                            ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG*CG(1,ng)
                          ELSE
                            ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                          ENDIF
                        ENDIF
                      ELSE IF(ITYP2(nr,nx).EQ.4)THEN
C                     Helmholtz equation
                        DO ni=1,NITB
                          PGMSI(ni)=PG(ms,NU1(ni),ng,nb)
                          PGNSI(ni)=PG(ns,NU1(ni),ng,nb)
                        ENDDO
                        SUM=0.0d0
                        DO mi=1,NITB
                          DO ni=1,NITB
                            SUM=SUM+PGMSI(mi)*PGNSI(ni)*GU(mi,ni)
                          ENDDO
                        ENDDO
                        RK2=CG(1,ng)**2 !is k^2
                        ES(mhs,nhs)=ES(mhs,nhs)
     '                    +(SUM-RK2*PG(ms,1,ng,nb)*PG(ns,1,ng,nb))*RWG
                      ELSE IF(ITYP2(nr,nx).EQ.5)THEN
C                     div(kgrad(u))=f
                        SUM=0.0d0
                        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                          PGMX=PGX(NBJ(nj),nj,ms,DXIX,PG(1,1,ng,nb))
                          PGNX=PGX(NBJ(nj),nj,ns,DXIX,PG(1,1,ng,nb))
                          SUM=SUM+CG(1+nj,ng)*PGMX*PGNX !CG(2..4)=cndty
                        ENDDO
                        ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                      ELSE IF(ITYP2(nr,nx).EQ.6)THEN
C                     2nd order elliptic
C DMAL 04 MARCH 2004  Replacing old code with a more general 
C 2nd order linear elliptic equation
C <OLDCODE>
C                        DO ni=1,NITB
C                          PGMSI(ni)=PG(ms,NU1(ni),ng,nb)
C                          PGNSI(ni)=PG(ns,NU1(ni),ng,nb)
C                        ENDDO
C                        SUM=0.0d0
C                        DO mi=1,NITB
C                          DO ni=1,NITB
C                            SUM=SUM+PGMSI(mi)*PGNSI(ni)*GU(mi,ni)
C                          ENDDO
C                        ENDDO
C                        ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
C <\OLDCODE>
                        IF(mh.EQ.nh) THEN
                          CGFLAG=1
                        ELSE
                          CGFLAG=0
                        ENDIF
                        S2=0.0d0
                        S3=0.0d0
                        S4=0.0d0
                        IF((NIT(nb).EQ.1).AND.(NJT.GT.1)) THEN
                          NJ_LOOP=1
                        ELSE
                          NJ_LOOP=NJ_LOC(NJL_GEOM,0,nr)
                        ENDIF
                        DO nj=1,NJ_LOOP
                          PGMX=PGX(NBJ(nj),nj,ms,DXIX,PG(1,1,ng,nb))
                          PGNX=PGX(NBJ(nj),nj,ns,DXIX,PG(1,1,ng,nb))
                          S2=S2+CG((mh-1)*(1+3*NJT)-1+3*nj,ng)*
     &                          PGMX*PGNX*CGFLAG
                          S3=S3+CG((mh-1)*(1+3*NJT)+3*nj,ng)*
     &                          PGM*PGNX*CGFLAG
                          S4=S4+CG((mh-1)*(1+3*NJT)+1+3*nj,ng)*
     &                          PGM*PGN*CGFLAG
                        ENDDO
C                       advection term & diffusion term
                        ES(mhs,nhs)=ES(mhs,nhs)-S2*RWG+S3*RWG+S4*RWG

                      ELSE IF(ITYP2(nr,nx).EQ.7)THEN
C                     Biharmonic equation
                      ELSE IF(ITYP2(nr,nx).EQ.8)THEN
C                     Fluid interface
                        DO ni=1,NITB !- Laplace equation
                          PGMSI(ni)=PG(ms,NU1(ni),ng,nb)
                          PGNSI(ni)=PG(ns,NU1(ni),ng,nb)
                        ENDDO
                        SUM=0.0d0
                        DO mi=1,NITB
                          DO ni=1,NITB
                            SUM=SUM+PGMSI(mi)*PGNSI(ni)*GU(mi,ni)
                          ENDDO
                        ENDDO
                        ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                      ELSE IF(ITYP2(nr,nx).EQ.9)THEN
C                     Oxygen transport
                        IF(ITYP3(nr,nx).EQ.1)THEN !multi-field
                          IF(KTYP15.EQ.1)THEN !one field only
                            DO ni=1,NITB
                              PGMSI(ni)=PG(ms,NU1(ni),ng,nb)
                              PGNSI(ni)=PG(ns,NU1(ni),ng,nb)
                            ENDDO
                            SUM=0.0d0
                            DO mi=1,NITB
                              DO ni=1,NITB
                                SUM=SUM+PGMSI(mi)*PGNSI(ni)*GU(mi,ni)
                              ENDDO
                            ENDDO
C                           add diffusive term & exchange term
                            ES(mhs,nhs)=ES(mhs,nhs)
     '                        +CG(1,ng)*CG(3,ng)*SUM*RWG
     '                        +CG(1,ng)*CG(9,ng)*PGM*PGN*RWG
                          ENDIF
                        ELSE IF(ITYP3(nr,nx).EQ.2)THEN
C                       Glucose-oxygen transport
C                      S3=0.0d0
C                      S4=0.0d0
C                      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C                        PGMX=PGX(NBJ(nj),nj,ms,DXIX,PG(1,1,ng,nb))
C                        PGNX=PGX(NBJ(nj),nj,ns,DXIX,PG(1,1,ng,nb))
C                        S3=S3 + CG(3,ng)*PGMX*PGNX
C                        S4=S4 + CG(4,ng)*PGMX*PGNX
C                      ENDDO
C                      SUM3=SUM3 + S3*RWG !diffusive term for oxygen
C                      SUM4=SUM4 + S4*RWG !diffusive term for glucose
                        ENDIF !ityp3
                      ENDIF !ityp2

                    ELSE IF(ITYP5(nr,nx).EQ.2)THEN
C                 *********** Time integration **************
                      IF(ITYP2(nr,nx).EQ.3)THEN
C                     Advection-diffusion
                        IF(mh.EQ.nh) THEN
                          CGFLAG=1
                        ELSE
                          CGFLAG=0
                        ENDIF
                        ED(mhs,nhs)=ED(mhs,nhs)+CGFLAG*
     &                    CG((mh-1)*(3+3*NJT)+3,ng)*PGM*PGN*RWG
                        S2=0.0d0
                        S3=0.0d0
                        S4=0.0d0
                        IF((NIT(nb).EQ.1).AND.(NJT.GT.1)) THEN
                          NJ_LOOP=1
                        ELSE
                          NJ_LOOP=NJ_LOC(NJL_GEOM,0,nr)
                        ENDIF
                        DO nj=1,NJ_LOOP
                          PGMX=PGX(NBJ(nj),nj,ms,DXIX,PG(1,1,ng,nb))
                          PGNX=PGX(NBJ(nj),nj,ns,DXIX,PG(1,1,ng,nb))
                          S2=S2+CG((mh-1)*(3+3*NJT)+1+3*nj,ng)*PGMX*
     &                      PGNX*CGFLAG
                          S3=S3+CG((mh-1)*(3+3*NJT)+2+3*nj,ng)*
     &                      CG((mh-1)*(3+3*NJT)+3,ng)*PGM*PGNX*CGFLAG
                          S4=S4+CG((mh-1)*(3+3*NJT)+3+3*nj,ng)*
     &                      PGM*PGN*CGFLAG
                        ENDDO
C                       advection term & diffusion term
                        ES(mhs,nhs)=ES(mhs,nhs)+S2*RWG+S3*RWG+S4*RWG
     &                    -CG((mh-1)*(3+3*NJT)+2,ng)*PGM*PGN*RWG*CGFLAG
                      ELSE IF(ITYP2(nr,nx).EQ.5)THEN
C                     VELOC not set
                        VELOC=0.d0
                        IF(ITYP3(nr,nx).EQ.1)THEN
C                       Fluid in an elastic tube
                          DSDXI=DSQRT(XG(1,2)**2+XG(2,2)**2)
                          PGMS=PG(ms,2,ng,nb)/DSDXI
                          PGNS=PG(ns,2,ng,nb)/DSDXI
C                         transient term
                          ED(mhs,nhs)=ED(mhs,nhs)+CG(2,ng)*PGM*PGN*RWG
C                         advection term & diffusion term
                          ES(mhs,nhs)=ES(mhs,nhs)
     '                      +CG(2,ng)*VELOC*PGM*PGNS
     '                      +CG(1,ng)*CG(3,ng)*PGMS*PGNS
                        ELSE IF(ITYP3(nr,nx).EQ.2)THEN
C                       Lung gas transport
                        ELSE IF(ITYP3(nr,nx).EQ.3)THEN
C                       General Navier-Stokes' equations
                        ELSE IF(ITYP3(nr,nx).EQ.4)THEN
C                       Stokes Flow equations
                          nbm=NBH(NH_LOC(mh,nx),nc)
                          nbn=NBH(NH_LOC(nh,nx),nc)
                          !continuity equation
                          IF (mh.EQ.1) THEN
                            IF (nh.GE.2) then !div.u
                              PGM=PG(ms,1,ng,nbm)
                              PGNX=PGX(NBJ(nh-1),nh-1,ns,DXIX,
     &                            PG(1,1,ng,nbn))
                              ES(mhs,nhs)=ES(mhs,nhs)+
     &                          PGM*PGNX*RWG
                            ENDIF
                          !momentum equations
                          ELSEIF (mh.GT.1) THEN
                            !dpdx
                            IF (nh.EQ.1) THEN
                              PGM=PG(ms,1,ng,nbm)
                              PGNX=PGX(NBJ(mh-1),mh-1,ns,DXIX,
     &                            PG(1,1,ng,nbn))
                              ES(mhs,nhs)=ES(mhs,nhs)+
     &                          PGM*PGNX*RWG/CG(2,ng)
                            !d2udx2
                            ELSEIF ((nh.GT.1).AND.(mh.eq.nh)) THEN 
                              PGM=PG(ms,1,ng,nbm)
                              PGN=PG(ns,1,ng,nbn)
                              ED(mhs,nhs)=ED(mhs,nhs)+PGM*PGN*RWG
                              SUM=0.0d0
                              DO nj=1,NJ_LOC(NJL_GEOM,0,nr) !add d2udx2
                                PGMX=PGX(NBJ(nj),nj,ms,DXIX,
     &                            PG(1,1,ng,nbm))
                                PGNX=PGX(NBJ(nj),nj,ns,DXIX,
     &                            PG(1,1,ng,nbn))
                                SUM=SUM+CG(nj+2,ng)*PGMX*PGNX
                              ENDDO
                              ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                            ENDIF
                          ENDIF
                        ELSE IF(ITYP3(nr,nx).EQ.5)THEN
C                       Bidoamin Stokes Flow equations
                          NJT=NJ_LOC(NJL_GEOM,0,nr)
                          nbm=NBH(NH_LOC(mh,nx),nc)
                          nbn=NBH(NH_LOC(nh,nx),nc)
                          !ext. continuity equation
                          IF(mh.EQ.1)THEN
                            IF(nh.EQ.1)THEN !Pme*pe
                              PGM=PG(ms,1,ng,nbm)
                              PGN=PG(ns,1,ng,nbn)
                              ES(mhs,nhs)=ES(mhs,nhs)+
     &                          PGM*PGN*RWG*CG(3,ng)
                            ENDIF
                            IF(((nh.GE.2).AND.(nh.LE.(NJT+1))))THEN !div.ue
                              nj=nh-1
                              PGM=PG(ms,1,ng,nbm)
                              PGNX=PGX(NBJ(nj),nj,ns,DXIX,
     &                            PG(1,1,ng,nbn))
                              ES(mhs,nhs)=ES(mhs,nhs)+
     &                          PGM*PGNX*RWG
                            ENDIF
                            IF(nh.EQ.(NJT+2))THEN !-Pmi*pi
                              PGM=PG(ms,1,ng,nbm)
                              PGN=PG(ns,1,ng,nbn)
                              ES(mhs,nhs)=ES(mhs,nhs)-
     &                          PGM*PGN*RWG*CG(4,ng)
                            ENDIF
                          !int. continuity equation
                          ELSEIF(mh.EQ.(NJT+2))THEN
                            IF(nh.EQ.1)THEN !-Pme*pe
                              PGM=PG(ms,1,ng,nbm)
                              PGN=PG(ns,1,ng,nbn)
                              ES(mhs,nhs)=ES(mhs,nhs)-
     &                          PGM*PGN*RWG*CG(3,ng)
                            ENDIF
                            IF(((nh.GE.(NJT+3)).AND.
     &                        (nh.LE.(2*NJT+2))))THEN !div.ui
                              nj=nh-NJT-2
                              PGM=PG(ms,1,ng,nbm)
                              PGNX=PGX(NBJ(nj),nj,ns,DXIX,
     &                            PG(1,1,ng,nbn))
                              ES(mhs,nhs)=ES(mhs,nhs)+
     &                          PGM*PGNX*RWG
                            ENDIF
                            IF(nh.EQ.(NJT+2))THEN !Pmi*pi
                              PGM=PG(ms,1,ng,nbm)
                              PGN=PG(ns,1,ng,nbn)
                              ES(mhs,nhs)=ES(mhs,nhs)+
     &                          PGM*PGN*RWG*CG(4,ng)
                            ENDIF
                           !momentum equations
                          ELSEIF((mh.GE.2).AND.(mh.LE.(NJT+1)))THEN
                            !dpdx
                            IF(nh.EQ.1)THEN
                              nj=mh-1
                              PGM=PG(ms,1,ng,nbm)
                              PGNX=PGX(NBJ(nj),nj,ns,DXIX,
     &                            PG(1,1,ng,nbn))
                              ES(mhs,nhs)=ES(mhs,nhs)+
     &                          PGM*PGNX*RWG/CG(5,ng)
                            !d2uedx2
                            ELSEIF((nh.GE.2).AND.(mh.LE.(NJT+1)).AND.
     &                        (mh.eq.nh))THEN 
                              PGM=PG(ms,1,ng,nbm)
                              PGN=PG(ns,1,ng,nbn)
                              ED(mhs,nhs)=ED(mhs,nhs)+PGM*PGN*RWG
                              SUM=0.0d0
                              DO nj=1,NJ_LOC(NJL_GEOM,0,nr) !add d2uedx2
                                PGMX=PGX(NBJ(nj),nj,ms,DXIX,
     &                            PG(1,1,ng,nbm))
                                PGNX=PGX(NBJ(nj),nj,ns,DXIX,
     &                            PG(1,1,ng,nbn))
                                SUM=SUM+CG(2*nj+5,ng)*PGMX*PGNX
                              ENDDO
                              ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                            ENDIF
                          ELSEIF((mh.GE.(NJT+3)).AND.
     &                      (mh.LE.(2*NJT+2)))THEN
                            !dpdx
                            IF(nh.EQ.(NJT+2))THEN
                              nj=mh-NJT-2
                              PGM=PG(ms,1,ng,nbm)
                              PGNX=PGX(NBJ(nj),nj,ns,DXIX,
     &                            PG(1,1,ng,nbn))
                              ES(mhs,nhs)=ES(mhs,nhs)+
     &                          PGM*PGNX*RWG/CG(6,ng)
                            !d2uidx2
                            ELSEIF((nh.GE.(NJT+3)).AND.
     &                        (mh.LE.(2*NJT+2)).AND.
     &                        (mh.eq.nh))THEN 
                              PGM=PG(ms,1,ng,nbm)
                              PGN=PG(ns,1,ng,nbn)
                              ED(mhs,nhs)=ED(mhs,nhs)+PGM*PGN*RWG
                              SUM=0.0d0
                              DO nj=1,NJ_LOC(NJL_GEOM,0,nr) !add d2uidx2
                                PGMX=PGX(NBJ(nj),nj,ms,DXIX,
     &                            PG(1,1,ng,nbm))
                                PGNX=PGX(NBJ(nj),nj,ns,DXIX,
     &                            PG(1,1,ng,nbn))
                                SUM=SUM+CG(2*nj+6,ng)*PGMX*PGNX
                              ENDDO
                              ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                            ENDIF
                          ENDIF
                        ENDIF !ityp3
                      ELSE IF(ITYP2(nr,nx).EQ.6)THEN
C                     Bioheat transfer equation
                        ED(mhs,nhs)=ED(mhs,nhs)+CG(2,ng)*PGM*PGN*RWG
                        SUM=0.0d0
                        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                          PGMX=PGX(NBJ(nj),nj,ms,DXIX,PG(1,1,ng,nb))
                          PGNX=PGX(NBJ(nj),nj,ns,DXIX,PG(1,1,ng,nb))
                          SUM=SUM+PGMX*PGNX
                        ENDDO
                        ES(mhs,nhs)=ES(mhs,nhs)+CG(3,ng)*SUM*RWG
     '                    +CG(5,ng)*PGM*PGN*RWG
                      ELSE IF(ITYP2(nr,nx).EQ.8)THEN
                        IF(ITYP3(nr,nx).EQ.2)THEN
C                       Cable equation (1D)
                          ED(mhs,nhs)=ED(mhs,nhs)+CG(1,ng)*PGM*PGN*RWG
                          PGMX=PGX(NBJ(nj),nj,ms,DXIX,PG(1,1,ng,nb))
                          PGNX=PGX(NBJ(nj),nj,ns,DXIX,PG(1,1,ng,nb))
                          ES(mhs,nhs)=ES(mhs,nhs)+CG(2,ng)*PGMX*PGNX*RWG
     '                      -CG(3,ng)*PGM*PGN*RWG
                        ELSE IF(ITYP3(nr,nx).EQ.3)THEN
C                       HH/Noble-dif equations
                          ED(mhs,nhs)=ED(mhs,nhs)+CG(2,ng)*PGM*PGN*RWG
                          SUM=0.0d0
                          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                            PGMX=PGX(NBJ(nj),nj,ms,DXIX,PG(1,1,ng,nb))
                            PGNX=PGX(NBJ(nj),nj,ns,DXIX,PG(1,1,ng,nb))
                            SUM=SUM+CG(1,ng)*PGMX*PGNX
                          ENDDO
                          ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                        ENDIF
                      ELSE IF(ITYP2(nr,nx).EQ.9)THEN
                        IF(ITYP3(nr,nx).EQ.1)THEN
C                       Multi-field oxygen transport
                        ELSE IF(ITYP3(nr,nx).EQ.2)THEN
C                       Glucose-oxygen transport
                          IF(mh.EQ.1.AND.nh.EQ.1)THEN
C                         transient term for oxygen
                            ED(mhs,nhs)=ED(mhs,nhs)+CG(1,ng)*PGM*PGN*RWG
                          ELSE IF(mh.EQ.2.AND.nh.EQ.2)THEN
C                         transient term for glucose
                            ED(mhs,nhs)=ED(mhs,nhs)+CG(2,ng)*PGM*PGN*RWG
                          ENDIF
                          S3=0.0d0
                          S4=0.0d0
                          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                            PGMX=PGX(NBJ(nj),nj,ms,DXIX,PG(1,1,ng,nb))
                            PGNX=PGX(NBJ(nj),nj,ns,DXIX,PG(1,1,ng,nb))
                            S3=S3+CG(3,ng)*PGMX*PGNX
                            S4=S4+CG(4,ng)*PGMX*PGNX
                          ENDDO
                          IF(mh.EQ.1.AND.nh.EQ.1)THEN
C                         diffusive term for oxygen
                            ES(mhs,nhs)=ES(mhs,nhs)+S3*RWG
                          ELSE IF(mh.EQ.2.AND.nh.EQ.2)THEN
C                         diffusive term for glucose
                            ES(mhs,nhs)=ES(mhs,nhs)+S4*RWG
                          ENDIF
C                         implicit nonlinear terms to go here
C                         if couple in at new time step
                        ENDIF
C.....................Pulmonary Transport
                      ELSE IF(ITYP2(nr,nx).EQ.11)THEN
                        SUM=0.0d0
                        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                          SUM=SUM+XG(nj,2)**2.0d0
                        ENDDO
                        DSDXI=DSQRT(SUM)
                        PGMS=PG(ms,2,ng,nb)/DSDXI !deriv of psi wrt s
                        PGNS=PG(ns,2,ng,nb)/DSDXI !deriv of psi wrt s
                        SCALE=1.d0 !area scaling for symmetric tree
                        IF(ITYP3(nr,nx).EQ.1)THEN
C.........................1-D advection and binary diffusion
                          IF(NORD1.EQ.1.AND.mhs.EQ.1) SCALE=2.d0
                            DIFF=PULMAT(1)
                          ED(mhs,nhs)=ED(mhs,nhs)+PI*XG(nj_radius,1)**2
     '                        *PGM*PGN*RWG*SCALE
                          IF(REMOVE_a)THEN
                            ES(mhs,nhs)=ES(mhs,nhs)+DIFF*PGMS*PI
     '                        *XG(nj_radius,1)**2*PGNS
     &                        *RWG*SCALE !*XG(nj_alveoli,1)
                          ELSE
                          ES(mhs,nhs)=ES(mhs,nhs)+DIFF*PGMS*PI
     '                        *XG(nj_radius,1)**2*PGNS
     &                        *RWG*SCALE *XG(nj_alveoli,1)
                          ENDIF
                          
                          EM(mhs,nhs)=EM(mhs,nhs)+XG(nj_flow,1)
     &                      *INLET_FLOW(nx)*PGM*PGNS*RWG*SCALE

                        ELSE IF(ITYP3(nr,nx).EQ.2)THEN
C...................... Water vapour concentration and temperature
                          IF(NORD1.EQ.1.AND.(mhs.EQ.1.OR.mhs.EQ.3))
     &                      SCALE=2.d0
                          rm=XG(nj_radius,1)
                          AREA=PI*rm**2.d0
                          VELOC=XG(nj_flow,1)*INLET_FLOW(nx)/AREA

                          IF(INLET_FLOW(nx).GE.0.d0)THEN !inspiration
                            gam=MAX(2.d0,WHT(1)*XG(nj_coeff,1))
                            pb =MAX(2.d0,WHT(2)*XG(nj_coeff,1))
                          ELSE
                             gam=MAX(2.d0,WHT(3)*XG(nj_coeff,1))
                             pb =MAX(2.d0,WHT(4)*XG(nj_coeff,1))
                          ENDIF
                          pa=MAX(2.d0,XG(nj_coeff,1)) !velocity coefficient
                          
                          IF(mh.EQ.1.AND.nh.EQ.1)THEN !heat
                             ED(mhs,nhs)=ED(mhs,nhs)+AREA*PULMAT(3)*
     '                            PULMAT(4)*PGM*PGN*RWG*SCALE
                             ES(mhs,nhs)=ES(mhs,nhs)+(XG(nj_flow,1)
     &                            *INLET_FLOW(nx)/AREA*PULMAT(3)
     &                            *PULMAT(4)*(PB+PA+4.d0)/(PB+PA+2.d0)
     &                            *PGM*PGNS+PULMAT(5)*PGMS*PGNS)*RWG
     &                            *AREA*SCALE+(2.d0*PULMAT(5)*(PB+2.d0)
     &                            /rm**2.d0*PGM*PGN)*RWG*AREA*SCALE
                             EM(mhs,nhs)=EM(mhs,nhs)+(-XG(nj_flow,1)
     &                            *INLET_FLOW(nx)/AREA*PULMAT(3)
     &                            *PULMAT(4)*2.d0*(PA+2.d0)/(PB*(PA+PB
     &                            +2.d0))*PGM*PGNS-PULMAT(5)*2.d0/PB
     &                            *PGMS
     &                            *PGNS)*RWG*AREA*SCALE+(2.d0*PULMAT(5)
     &                            *(PB+2.d0)/rm**2.d0*PGM*PGN)*RWG*AREA
     &                            *SCALE
                             
                          ELSEIF(mh.EQ.1.AND.nh.EQ.2)THEN !water in heat
                             IF(ne.GT.ne_trachea)THEN !i.e. not et tube
                                !1/A.int( 2.PI.r.DH.h.del(C)/del(r))
                                ES(mhs,nhs)=ES(mhs,nhs)+2.d0*(PULMAT(1)
     '                               *PULMAT(6)*(PB+2.d0)*gam/(PB*(gam
     '                               +1.d0))*PGM*PGN)*2.d0*PI*RWG/1.d9
                                EM(mhs,nhs)=EM(mhs,nhs)+2.d0*(PULMAT(1)
     '                               *PULMAT(6)*(PB+2.d0)*gam/(PB*(gam
     '                               +1.d0))*PGM*PGN)*2.d0*PI*RWG/1.d9
                             ENDIF
                             
                          ELSEIF(mh.EQ.2.AND.nh.EQ.2)THEN !water
                             ED(mhs,nhs)=ED(mhs,nhs)+PGM*PGN*RWG*AREA
     '                            *SCALE
                             ES(mhs,nhs)=ES(mhs,nhs)+(XG(nj_flow,1)
     &                            *INLET_FLOW(nx)/AREA*(PA+gam+4.d0)/(PA
     &                            +gam+2.d0)*PGM*PGNS+PULMAT(1)*PGMS
     &                            *PGNS)*RWG*AREA*SCALE
                             EM(mhs,nhs)=EM(mhs,nhs)+(-XG(nj_flow,1)
     &                            *INLET_FLOW(nx)/AREA*2.d0*(PA+2.d0)
     &                            /(gam*(PA+gam+2.d0))*PGM*PGNS
     &                            -PULMAT(1)*2.d0/gam*PGMS*PGNS)*RWG
     &                            *AREA*SCALE
                             IF(ne.GE.ne_trachea)THEN
                                ES(mhs,nhs)=ES(mhs,nhs)+2.d0*PULMAT(1)*
     '                               (gam+2.d0)/rm**2.d0*PGM*PGN*RWG
     &                               *AREA*SCALE
                                EM(mhs,nhs)=EM(mhs,nhs)+2.d0*PULMAT(1)*
     '                               (gam+2.d0)/rm**2.d0*PGM*PGN*RWG
     &                               *AREA*SCALE
                             ENDIF !ne.GE.ne_trachea
                          ENDIF !mh.EQ.1.AND.nh.EQ.1

                        ELSEIF(ITYP3(nr,nx).EQ.5)THEN !MCT
                          IF(NORD1.EQ.1.AND.mhs.EQ.1) SCALE=2.d0
                          ED(mhs,nhs)=ED(mhs,nhs)+PI*XG(nj_radius,1)**2
     &                      *PGM*PGN*RWG*SCALE
                          EM(mhs,nhs)=EM(mhs,nhs)+PI*XG(nj_radius,1)**2
     &                      *XG(nj_flow,1)*PGM*PGNS*RWG*SCALE

                        ELSEIF(ITYP3(nr,nx).EQ.5)THEN !1D lung Navier-Stokes
                          IF(NORD1.EQ.1.AND.mhs.EQ.1) SCALE=2.d0
                          ED(mhs,nhs)=ED(mhs,nhs)+PGM*PGN*RWG*SCALE
                          IF(mh.EQ.1.AND.nh.EQ.1)THEN ! area
c                            ES(mhs,nhs)=ES(mhs,nhs)+(A*PGM*PGNS+K*PGMS
c     &                        *PGNS)*RWG*SCALE
                          ELSEIF(mh.EQ.2.AND.nh.EQ.2)THEN ! flow
                          ELSEIF(mh.EQ.1.AND.nh.EQ.2)THEN ! 
                          ENDIF

                        ENDIF !ityp3
                      ENDIF !ityp2

                    ELSE IF(ITYP5(nr,nx).EQ.3)THEN
C                 ********* Modal analysis ************
                      IF(ITYP2(nr,nx).EQ.3)THEN
C                     Laplace equation
                        DO ni=1,NITB
                          PGMSI(ni)=PG(ms,NU1(ni),ng,nb)
                          PGNSI(ni)=PG(ns,NU1(ni),ng,nb)
                        ENDDO
                        SUM=0.0d0
                        DO mi=1,NITB
                          DO ni=1,NITB
                            SUM=SUM+PGMSI(mi)*PGNSI(ni)*GU(mi,ni)
                          ENDDO
                        ENDDO
                        ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                      ELSE IF(ITYP2(nr,nx).EQ.4)THEN
C                     Helmholtz equation
                        DO ni=1,NITB
                          PGMSI(ni)=PG(ms,NU1(ni),ng,nb)
                          PGNSI(ni)=PG(ns,NU1(ni),ng,nb)
                        ENDDO
                        SUM=0.0d0
                        DO mi=1,NITB
                          DO ni=1,NITB
                            SUM=SUM+PGMSI(mi)*PGNSI(ni)*GU(mi,ni)
                          ENDDO
                        ENDDO
                        RK2=CG(1,ng)**2 !is k^2
                        ES(mhs,nhs)=ES(mhs,nhs)
     '                    +(SUM-RK2*PG(ms,1,ng,nb)*PG(ns,1,ng,nb))*RWG
                      ELSE IF(ITYP2(nr,nx).EQ.5)THEN
C                     div(kgrad(u))=f
                        SUM=0.0d0
                        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                          PGMX=PGX(NBJ(nj),nj,ms,DXIX,PG(1,1,ng,nb))
                          PGNX=PGX(NBJ(nj),nj,ns,DXIX,PG(1,1,ng,nb))
                          SUM=SUM+CG(1+nj,ng)*PGMX*PGNX !CG(2..4)=condty
                        ENDDO
                        ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                      ELSE IF(ITYP2(nr,nx).EQ.6)THEN
C                     2nd order elliptic equation
                        DO ni=1,NITB
                          PGMSI(ni)=PG(ms,NU1(ni),ng,nb)
                          PGNSI(ni)=PG(ns,NU1(ni),ng,nb)
                        ENDDO
                        SUM=0.0d0
                        DO mi=1,NITB
                          DO ni=1,NITB
                            SUM=SUM+PGMSI(mi)*PGNSI(ni)*GU(mi,ni)
                          ENDDO
                        ENDDO
                        ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG
                      ELSE IF(ITYP2(nr,nx).EQ.7)THEN
C                     Biharmonic equation
                      ELSE IF(ITYP2(nr,nx).EQ.8)THEN
C                     Vocal tract equations
                        EM(mhs,nhs)=EM(mhs,nhs)+PGM*PGN*RWG/CG(1,ng)
                        PGMX=PGX(nb,1,ms,DXIX,PG(1,1,ng,nb))
                        PGNX=PGX(nb,1,ns,DXIX,PG(1,1,ng,nb))
                        ES(mhs,nhs)=ES(mhs,nhs)+PGMX*PGNX*RWG*
     '                    CG(2,ng)**2/CG(1,ng)
                      ENDIF !ityp2
                    ENDIF !ityp5
                  ENDDO !ns
                ENDDO !nh
              ENDIF !update_matrix

C***        RHS terms
              IF(UPDATE_VECTOR)THEN
                IF(ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4)THEN
C               ****** Static or quasi-static solution ******
                  IF(ITYP2(nr,nx).EQ.3)THEN
C                 Laplace equation
                    ER(mhs)=0.0d0
                  ELSE IF(ITYP2(nr,nx).EQ.4)THEN
C                 Helmholtz equation
                    ER(mhs)=0.0d0
                  ELSE IF(ITYP2(nr,nx).EQ.5)THEN
C                 div(kgrad(u))=f
                    IF(ILP(1,1,nr,nx).EQ.5)THEN
                      SOURCE=YG(1,ng) !domain source f
! PJH 26Jan96     IF(KTYP32.EQ.2) THEN !bidomain
!                   SOURCE=YG(1,ng) !coupling to intracellular terms
                    ELSE
                      SOURCE=CG(1,ng) !domain source f
                    ENDIF
                    ER(mhs)=ER(mhs)+SOURCE*PGM*RWG
                  ELSE IF(ITYP2(nr,nx).EQ.6)THEN
C                 2nd order elliptic equation
                    ER(mhs)=ER(mhs)-CG((mh-1)*(1+3*NJT)+1,ng)*PGM*RWG
                  ELSE IF(ITYP2(nr,nx).EQ.7)THEN
C                 Biharmonic equation
                    ER(mhs)=0.0d0
                  ELSE IF(ITYP2(nr,nx).EQ.8)THEN
C                 Fluid interface
                    ER(mhs)=0.0d0
                  ELSE IF(ITYP2(nr,nx).EQ.9)THEN
C                 Oxygen transport
                    IF(ITYP3(nr,nx).EQ.1)THEN
C                   Multi-field oxygen transport
                      IF(KTYP15.EQ.1)THEN
C                     Single-field only
                        ER(mhs)=ER(mhs)+(CG(1,ng)*CG(9,ng)*CG(10,ng)
     '                    -CG(8,ng))*PGM*RWG
                      ENDIF
                    ELSE IF(ITYP3(nr,nx).EQ.2)THEN
C                   Glucose-oxygen transport
                    ENDIF
                  ENDIF

                ELSE IF(ITYP5(nr,nx).EQ.2)THEN
C               ***** Time integration *****
                  IF(ITYP2(nr,nx).EQ.3)THEN
C                 Advection-diffusion equation
                    ER(mhs)=ER(mhs)+CG((mh-1)*(3+3*NJT)+1,ng)*PGM*RWG
                  ELSE IF(ITYP2(nr,nx).EQ.5)THEN
                    IF(ITYP3(nr,nx).EQ.4)THEN
                    ! Stokes' Flow
                      IF (mh.EQ.1) THEN
                        nb=NBH(NH_LOC(mh,nx),nc)
                        PGM=PG(ms,1,ng,nb)
                        ER(mhs)=ER(mhs)+CG(1,ng)*PGM*RWG
                      ENDIF
                    ENDIF
                    IF(ITYP3(nr,nx).EQ.5)THEN
                    ! Bidomain Stokes' Flow
                      IF (mh.EQ.1) THEN
                        nb=NBH(NH_LOC(mh,nx),nc)
                        PGM=PG(ms,1,ng,nb)
                        ER(mhs)=ER(mhs)+CG(1,ng)*PGM*RWG
                      ELSEIF (mh.EQ.(NJT+2)) THEN
                        nb=NBH(NH_LOC(mh,nx),nc)
                        PGM=PG(ms,1,ng,nb)
                        ER(mhs)=ER(mhs)+CG(2,ng)*PGM*RWG
                      ENDIF
                    ENDIF
                  ELSE IF(ITYP2(nr,nx).EQ.6)THEN
C                 Bioheat transfer equation
                    ER(mhs)=ER(mhs)+(CG(5,ng)*CG(4,ng)+CG(1,ng))*PGM*RWG
                  ELSE IF(ITYP2(nr,nx).EQ.8)THEN
                    IF(ITYP3(nr,nx).EQ.2)THEN
C                   Cable equation (1D)
                      CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,
     '                  *9999)
                      V=ZG(1,1)
                      ER(mhs)=ER(mhs)+CG(3,ng)*V*V*(V/(CG(4,ng)
     '                  *CG(5,ng))-(1.d0/CG(4,ng)+1.d0/CG(5,ng)))
                    ENDIF
                  ELSE IF(ITYP2(nr,nx).EQ.9)THEN
C                 Oxygen transport
                    IF(ITYP3(nr,nx).EQ.1)THEN !multi-field oxy transport
                      IF(KTYP15.EQ.1)THEN !single-field only
                        ER(mhs)=ER(mhs)+(CG(1,ng)*CG(9,ng)*CG(10,ng)+
     '                    CG(8,ng))*PGM*RWG
                      ENDIF
                    ELSE IF(ITYP3(nr,nx).EQ.2)THEN
C                   Glucose-oxygen transport
                      CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,
     '                  *9999)
                      ER(mhs)=0.d0 !Explicit nonlinear terms to go here
                    ENDIF
                  ELSE IF(ITYP2(nr,nx).EQ.11)THEN !pulmonary transport
                    IF(ITYP3(nr,nx).EQ.1)THEN
                      ER(mhs)=0.d0
c
                      IF(PULMAT(2).NE.0.0d0)THEN
                        Qc=110.0d3
                        Qb=Qc/3.65d6 !need to not hard code this but breaks example91.com
                        ER(mhs)=ER(mhs)+Qb*PULMAT(2)*PI*(XG(nj_radius,1)
     &                    **2.0d0)*PGM*RWG*SCALE
c                        copy_PGN=PGN ! Note global if used
                      ENDIF
c
c                      write(*,*)'ER(mhs)',ER(mhs)
                      
                    ELSE IF(ITYP3(nr,nx).EQ.2)THEN
                      ER(mhs)=0.d0
c                      write(*,*)'2'
                    ELSE IF(ITYP3(nr,nx).EQ.5)THEN
                      ER(mhs)=ER(mhs)+2.d0*PI*XG(nj_radius,1)
     &                  *XG(nj_source,1)*PGM*RWG
                    ENDIF
                  ENDIF !ityp2

                ELSE IF(ITYP5(nr,nx).EQ.3)THEN
C***            ****** Modal analysis ******
                  ER(mhs)=0.d0
                ENDIF

                IF(ITYP2(nr,nx).EQ.5.AND.KTYP60.GT.0)THEN
C               div(kgrad(u))=f with growing tips
                  SOURCE=0.0d0
C PJH 1-JUN-91  Note!!!! use of ns for nn here is temporary -
                  DO no_tip=1,NP_GROWING_TIP(0,ne)
                    SOURCE=SOURCE+FLOW_GROWING_TIP(no_tip,ne)*
     '                PSI1(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),1,1,1,
     '                ns,XI_GROWING_TIP)
                  ENDDO
                  ER(mhs)=ER(mhs)+SOURCE
                ENDIF
              ENDIF !update_vector
            ENDDO !ms
          ENDDO !mh
        ENDDO !ng

C   Scale factor adjustment and mass lumping
        mhs=0
        DO mh=1,NHE
          DO ms=1,NSTB
            mhs=mhs+1
            nhs=0
            IF(UPDATE_MATRIX)THEN
              DO nh=1,NHE
                DO ns=1,NSTB
                  nhs=nhs+1
                  SS=SE(ms,nb)*SE(ns,nb)
                  EM(mhs,nhs)=EM(mhs,nhs)*SS
                  ED(mhs,nhs)=ED(mhs,nhs)*SS
                  ES(mhs,nhs)=ES(mhs,nhs)*SS
                  IF(DOP)THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,'('' mhs='',I2,'' nhs='',I2,'
     '                //''' ED(mhs,nhs)='',E10.3,'' ES(mhs,nhs)='','
     '                //'E10.3)') mhs,nhs,ED(mhs,nhs),ES(mhs,nhs)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF
                ENDDO !ns
              ENDDO !nh
            ENDIF
            IF(UPDATE_VECTOR) ER(mhs)=ER(mhs)*SE(ms,nb)
          ENDDO !ms
        ENDDO !mh

        IF(UPDATE_MATRIX)THEN
          IF(LUMP)THEN
            DO mhs=1,NHST
              SUM=0.0d0
              DO nhs=1,NHST
                SUM=SUM+EM(mhs,nhs)
                EM(mhs,nhs)=0.0d0
              ENDDO
              EM(mhs,mhs)=SUM
            ENDDO
          ENDIF !lump
        ENDIF !update_matrix
      ENDIF !update

      CALL EXITS('XPES30')
      RETURN
 9999 CALL ERRORS('XPES30',ERROR)
      CALL EXITS('XPES30')
      RETURN 1
      END



