      SUBROUTINE PLATE(NBH,NBJ,NHE,nr,NW,nx,
     '  CE,CG,ER,ES,PG,RG,SE,WG,XE,XG,ERROR,*)

C#### Subroutine: PLATE
C###  Description:
C###    PLATE calculates the element load vector ER, the stiffness
C###    matrix ES and, if ITYP5(nr,nx) >1, the mass matrix EM and the
C###    damping matrix ED for the Kirchhoff plate element (NW=6):
C###    lines normal to undeformed middle surface remain straight and
C###    normal to deformed surface (transverse shearis zero and
C###    governing equation is biharmonic in transverse displacements).
C###    Plate element allows in-plane deformation. Geometry of plate
C###    is assumed to be bilinear.
C###    The 1st in-plane displacement is in the direction of 1st Xi(1)
C###    edge of element;
C###    The 2nd in-plane displacement is in the direction normal to
C###    Xi(1) edge of element.

C**** Note: SUM1 accumulates the in-plane (mh=1,2) stiffness matrix
C****                         or flexural (mh=3)       "        "
C****       SUM2      "          geometric             "        "
C****       SUM3      "          mass matrix.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),NHE,nr,NW,nx
      REAL*8 CE(NMM),CG(NMM,NGM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     '  WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,mb,mh,mhg,mhs,mhsg,mk,mn,ms,MSTOT,mvar,
     '  nb,ng,nh,nhg,nhs,nhsg,nhx,ni,nj,nk,nm,nn,ns,NST_12,
     '  NU2(3,3),nvar
      REAL*8 A(3),AA(16,9),ASUM,B(3),BSUM,C(3),CHTOFF(3,3,3),CSUM,
     '  D(3),DBM(3,3,3),DENS,DXIX(3,3),
     '  EP(24,24),FLEX,GL(3,3),GU(3,3),
     '  PM,PMX,PMY,PN,PNX,PNY,POIS,RT(3,3),RWG,
     '  SUM,SUM1,SUM2,SUM3,THETA,THIC,X3G(4,3),YMOD,YY
      DATA NU2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('PLATE',*9999)
      IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3.AND.NBH(NH_LOC(3,nx)).GT.0) THEN !3D & bending included
        NST_12=4*NKT(0,NBH(NH_LOC(3,nx))) !is NST for nh=1,2
      ELSE                              !2D or no bending
        NST_12=4*NKT(0,NBH(NH_LOC(1,nx))) !is NST for nh=1,2
      ENDIF
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' NST_12 = '',I3)') NST_12
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      !Note that XGMG is called with ip=1 to enable DXIX to be calc.d
      !wrt fibre coords.
      DO ng=1,NGT(NBJ(1))
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
        CALL XGMG(1,NIT(NBJ(1)),NBJ(1),nr,DXIX,GL,GU,
     '    RG(ng),XG,ERROR,*9999)
      ENDDO

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' RG(ng): '',9E10.3)') (RG(ng),ng=1,9)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' CE(nm): '',14F6.2)') (CE(nm),nm=1,NMM)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' NVE='',I1,'' NHE='',I1)') NVE(NW),NHE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO ni=1,NIT(NBJ(1))
          WRITE(OP_STRING,
     '    '('' DXIX('',I1,'',j) :'',2E12.3)') ni,(DXIX(ni,j),j=1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      !Compute direction cosines assuming plate has bilinear geometry
      ASUM=0.0d0
      BSUM=0.0d0
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        A(nj)=XE(2,nj)-XE(1,nj)
        B(nj)=XE(3,nj)-XE(1,nj)
        ASUM=ASUM+A(nj)*A(nj)
        BSUM=BSUM+B(nj)*B(nj)
      ENDDO
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        A(nj)=A(nj)/DSQRT(ASUM) !is unit vector along Xi(1) side of plate
        B(nj)=B(nj)/DSQRT(BSUM) !is unit vector along Xi(2) side of plate
      ENDDO
      IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
        THETA=DATAN2(A(2),A(1))
        RT(1,1)= DCOS(THETA)
        RT(1,2)= DSIN(THETA)
        RT(2,1)=-DSIN(THETA)
        RT(2,2)= DCOS(THETA)
      ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
        CALL CROSS(A,B,C)
        CSUM=C(1)*C(1)+C(2)*C(2)+C(3)*C(3)
        DO nj=1,3
          C(nj)=C(nj)/DSQRT(CSUM) !is unit vector normal to plate
        ENDDO
        CALL CROSS(C,A,D) !D is unit vector normal to A and C (2nd in-plane)
        DO NJ=1,3
          RT(1,nj)=A(nj)
          RT(2,nj)=D(nj)
          RT(3,nj)=C(nj)
        ENDDO
      ENDIF
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,3
          nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
          WRITE(OP_STRING,
     '      '('' RT('',I1,'',nhg): '',3F7.3)') nh,(RT(nh,nhg),nhg=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      IF(NHE.EQ.3) THEN !Transverse behaviour included
        nb=NBH(NH_LOC(3,nx))
        X3G(1,1)=0.0d0
        X3G(2,1)=0.0d0
        X3G(3,1)=0.0d0
        X3G(4,1)=0.0d0
        X3G(1,2)=0.0d0
        X3G(2,2)=0.0d0
        X3G(3,2)=0.0d0
        X3G(4,2)=0.0d0
        X3G(1,3)=0.0d0
        X3G(2,3)=0.0d0
        X3G(3,3)=0.0d0
        X3G(4,3)=0.0d0
        DO ng=1,NGT(nb)
          CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF,DBM,GU,XG,X3G,
     '      .FALSE.,ERROR,*9999)
          DO ns=1,NST(nb)
            SUM=0.d0
            DO i=1,2
              DO j=1,2
                SUM=SUM
     '            +(PG(ns,NU2(i,j),ng,nb)-CHTOFF(1,i,j)*PG(ns,2,ng,nb)
     '            -CHTOFF(2,i,j)*PG(ns,4,ng,nb))*GU(i,j)
              ENDDO
            ENDDO
            AA(ns,ng)=SUM !is Laplacian of psi(ns)
          ENDDO
        ENDDO
      ENDIF

      mhs=0
      DO mvar=1,NHE
        mh=NHV(mvar,NW)
        mb=NBH(mh)
        DO ms=1,NST(mb)
          mhs=mhs+1
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '      '('' mvar='',I1,'' mh='',I1,'' ms='',I2)') mvar,mh,ms
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(mvar.EQ.3) THEN !Apply transverse loading
            SUM=0.d0
            DO ng=1,NGT(mb)
              RWG=RG(ng)*WG(ng,mb)
              SUM=SUM+CG(NMP(NW)+NGP(NW)+1,ng)*PG(ms,1,ng,mb)*RWG
            ENDDO
            SUM=SUM*SE(ms,mb)
            MSTOT=0
            DO mhg=1,NHE
              mhsg=ms+MSTOT
              MSTOT=MSTOT+NST_12
              ER(mhsg)=ER(mhsg)+RT(mh,mhg)*SUM
            ENDDO
          ENDIF

          nhs=0
          DO nvar=1,NHE
            nhx=NHV(nvar,NW)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh)
            DO ns=1,NST(nb)
              nhs=nhs+1
              SUM1=0.0d0
              SUM2=0.0d0
              SUM3=0.0d0
              DO ng=1,NGT(nb)
                RWG =RG(ng)*WG(ng,nb)
                THIC=CG(NMP(NW)+1,ng)                  !plate thickness
                YMOD=CG(1,ng)                          !Young's moulus
                POIS=CG(2,ng)                          !Poisson's ratio
                FLEX=YMOD*THIC**3/(12.0d0*(1.0d0-POIS**2)) !flexural rigidity
                IF(ITYP5(nr,nx).EQ.1) THEN
                  DENS=0.d0
                ELSE
                  DENS=CG(ILT(NW,nr,nx),ng)*CG(NMP(NW)+1,ng) !mass density
                ENDIF
                YY=YMOD/(1.d0-POIS*POIS)
C GBS 19/09/95 - TEMP not used
c                IF(KTYP43.GT.0) THEN !thermal strain
c                  IF(ILP(4,6,nr,nx).LE.3) THEN
c                    TEMP=CG(4,ng)
c                  ELSE IF(ILP(4,6,nr,nx).EQ.4) THEN
c                    TEMP=YG(1,ng)                      !temperature
c                  ENDIF
c                ENDIF
                PM =PG(ms,1,ng,mb)
                PN =PG(ns,1,ng,nb)
                PMX=PG(ms,2,ng,mb)*DXIX(1,1)+PG(ms,4,ng,mb)*DXIX(2,1)
                PMY=PG(ms,2,ng,mb)*DXIX(1,2)+PG(ms,4,ng,mb)*DXIX(2,2)
                PNX=PG(ns,2,ng,nb)*DXIX(1,1)+PG(ns,4,ng,nb)*DXIX(2,1)
                PNY=PG(ns,2,ng,nb)*DXIX(1,2)+PG(ns,4,ng,nb)*DXIX(2,2)
                IF(mh.EQ.1.AND.nh.EQ.1) THEN !in-plane
                  SUM1=SUM1+THIC*YY*(PMX*PNX+(1.d0-POIS)/2.d0*PMY*PNY)
     '              *RWG
                  SUM2=0.d0
                  SUM3=SUM3+DENS*PM*PN*RWG
                ELSE IF(mh.EQ.1.AND.nh.EQ.2) THEN !in-plane
                  SUM1=SUM1+THIC*YY*(POIS*PMX*PNY+(1.d0-POIS)/2.d0*PMY
     '              *PNX)*RWG
                  SUM2=0.d0
                  SUM3=0.d0
                ELSE IF(mh.EQ.2.AND.nh.EQ.1) THEN !in-plane
                  SUM1=SUM1+THIC*YY*(POIS*PMY*PNX+(1.d0-POIS)/2.d0*PMX
     '              *PNY)*RWG
                  SUM2=0.d0
                  SUM3=0.d0
                ELSE IF(mh.EQ.2.AND.nh.EQ.2) THEN !in-plane
                  SUM1=SUM1+THIC*YY*(PMY*PNY+(1.d0-POIS)/2.d0*PMX
     '              *PNX)*RWG
                  SUM2=0.d0
                  SUM3=SUM3+DENS*PM*PN*RWG
                ELSE IF(mh.EQ.3.AND.nh.EQ.3) THEN !bending
                  SUM1=SUM1+FLEX*RWG*AA(ms,ng)*AA(ns,ng)
                  SUM2=0.d0
c???? PJH         SUM2=SUM2+PG(ms,2,ng,mb)*PG(ns,2,ng,nb)*RWG
                  SUM3=SUM3+DENS*PM*PN*RWG
                ENDIF
              ENDDO
              SUM1=SUM1*SE(ms,mb)*SE(ns,nb)
              SUM2=SUM2*SE(ms,mb)*SE(ns,nb)
              SUM3=SUM3*SE(ms,mb)*SE(ns,nb)
c             IF(DOP) THEN
c               WRITE(IO4,'('' mvar='',I1,'' mh='',I1,'' ms='',I2,
c    '            '' nvar='',I1,'' nh='',I1,'' ns='',I2)')
c    '            mvar,mh,ms,nvar,nh,ns
c               WRITE(IO4,'('' sum1='',E12.3,'' sum2='',E12.3,'
c    '            //''' sum3='',E12.3)') SUM1,SUM2,SUM3
c             ENDIF

              IF(ITYP5(nr,nx).EQ.1) THEN !static
                EP(mhs,nhs)=SUM1
              ELSE IF(ITYP5(nr,nx).EQ.2.OR.ITYP5(nr,nx).EQ.3) THEN !time or modal
c cpb 1/5/95 Replacing Fourier analysis with Quasi-static analysis
c              ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Fourier
              ENDIF

c             MSTOT=0
c             DO mhg=1,NHE
c               mhsg=ms+MSTOT
c               MSTOT=MSTOT+NST_12
c               NSTOT=0
c               DO nhg=1,NHE
c                 nhsg=ns+NSTOT
c                 NSTOT=NSTOT+NST_12
c                 IF(DOP) THEN
c                   WRITE(IO4,'('' mhg='',I2,'' mhsg='',I2,'' nhg='',I2,
c    '                '' nhsg='',I2)') mhg,mhsg,nhg,nhsg
c                 ENDIF
c                 IF(ITYP5(nr,nx).EQ.1) THEN !static
c                   ES(mhsg,nhsg)=ES(mhsg,nhsg)
c    '                +RT(mh,mhg)*(    SUM1 + SUM2         )*RT(nh,nhg)
c                 ELSE IF(ITYP5(nr,nx).EQ.2.OR.ITYP5(nr,nx).EQ.3) THEN !time or modal
c                   ES(mhsg,nhsg)=ES(mhsg,nhsg)
c    '                +RT(mh,mhg)*(SUM1+CG(ILT(NW,nr,nx),ng)*SUM2)*RT(nh,nhg)
c                   EM(mhsg,nhsg)=EM(mhsg,nhsg)
c    '                +RT(mh,mhg)*(       SUM3             )*RT(nh,nhg)
c                 ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Fourier
c                   ES(mhsg,nhsg)=ES(mhsg,nhsg)
c    '               +RT(mh,mhg)*(        SUM1             )*RT(nh,nhg)
c                   EM(mhsg,nhsg)=EM(mhsg,nhsg)
c    '               +RT(mh,mhg)*(        SUM2             )*RT(nh,nhg)
c                 ENDIF
c               ENDDO
c             ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO mhs=1,24
          WRITE(OP_STRING,'('' EP('',I2,'',nhs):'',24F6.2)')
     '      mhs,(EP(mhs,nhs),nhs=1,24)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      DO mhg=1,3
       ms=0
        DO mn=1,4
        DO mk=1,NKT(0,NBH(mhg))
         ms=ms+1
          mhsg=1+4*(ms-1)+16*(mhg-1)
          DO nhg=1,3
            ns=0
            DO nn=1,4
            DO nk=1,NKT(0,NBH(nhg))
              ns=ns+1
              nhsg=1+4*(ns-1)+16*(nhg-1)
              ES(mhsg,nhsg)=ES(mhsg,nhsg)
     '          +RT(1,mhg)*EP(  mn,  nn)*RT(1,nhg)
     '          +RT(2,mhg)*EP(4+mn,  nn)*RT(1,nhg)
     '          +RT(3,mhg)*EP(8+ms,  nn)*RT(1,nhg)
     '          +RT(1,mhg)*EP(  mn,4+nn)*RT(2,nhg)
     '          +RT(2,mhg)*EP(4+mn,4+nn)*RT(2,nhg)
     '          +RT(3,mhg)*EP(8+ms,4+nn)*RT(2,nhg)
     '          +RT(1,mhg)*EP(  mn,8+ns)*RT(3,nhg)
     '          +RT(2,mhg)*EP(4+mn,8+ns)*RT(3,nhg)
     '          +RT(3,mhg)*EP(8+ms,8+ns)*RT(3,nhg)
            ENDDO
            ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDDO

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO mhsg=1,48
          WRITE(OP_STRING,'('' ES('',I2,'',nhs):'',48F6.2)')
     '      mhsg,(ES(mhsg,nhsg),nhsg=1,48)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('PLATE')
      RETURN
 9999 CALL ERRORS('PLATE',ERROR)
      CALL EXITS('PLATE')
      RETURN 1
      END


