      SUBROUTINE BEAM(NBH,NBJ,NHE,nw,
     '  CE,CG,EM,ER,ES,PG,RG,SE,WG,XE,XG,ERROR,*)

C#### Subroutine: BEAM
C###  Description:
C###    BEAM calculates the element load vector ER, the stiffness
C###    matrix ES and, if ITYP5(nr,nx)>1, the mass matrix EM and the
C###    damping matrix ED for the beam element (nw=3).

C**** Kirchhoff beam: lines normal to undeformed middle surface remain
C****   straight and normal to deformed surface (transverse shear zero
C****   and governing equation biharmonic in transverse displacements).
C**** Mindlin beam  : lines normal to undeformed middle surface remain
C****   straight but not normal to deformed surface (transverse shear
C****   included with rotation variables in 2nd order equation).
C**** Discrete Kirchoff: shear deformation ignored in energy;
C****   rotations = slope of transverse displacement at discrete points.
C**** All beam elements allow axial deformation and include bending
C****   stiffness in two directions when element is in 3D space.
C**** Mindlin beam elements also have torsional stiffness and allow
C****   axial twist.
C****
C**** KTYP45 = 1 beam cross-section is rectangular solid
C****   "    = 2 beam cross-section is rectangular tube
C****   "    = 3 beam cross-section is ellipsoidal solid
C****   "    = 4 beam cross-section is ellipsoidal tube
C****
C**** NMP(3)=number of material params=2 or 4 if thermal strain included
C**** NGP(3)=number of geometric params=2(solid) or 3(tube)
C**** NLP(3)=number of load params=4(njt=2) or 6(njt=3)
C****
C**** Element degrees of freedom:
C**** mh=1 is axial displacement and (optionally) axial rotation
C**** mh=2 is y-dir.n displacement and rotation in x,y-plane
C**** mh=3 is z-dir.n displacement and rotation in x,z-plane (if njt=3)
C****
C**** Element material parameters:
C**** CG(1,ng) = Young's modulus (GPa)
C**** CG(2,ng) = Poisson's ratio
C**** CG(3,ng) = Depth of beam (m)
C**** CG(4,ng) = Width of beam (m)
C**** CG(5,ng) = Wall thickness of beam (m) (if tubular x-section only)
C**** CG(NGP(3)+1,ng) = axial load
C**** CG(NGP(3)+2,ng) = axial torque
C**** CG(NGP(3)+3,ng) = transverse load (y-direction)
C**** CG(NGP(3)+4,ng) = x,y-plane moment
C**** CG(NGP(3)+5,ng) = transverse load (z-direction)
C**** CG(NGP(3)+6,ng) = x,z-plane moment
C**** CG(ILT(3,nr,nx),ng) is density (kg/m^3)
C****
C**** SUM1 accumulates the axial or flexural stiffness matrix
C**** SUM2     "           geometric            "        "
C**** SUM3     "           mass matrix.
C**** RT(nh,nhg) rotates variables nh=1,nhe in the element coord system
C****   to equivalent variables nhg in the global system.
C**** Direction cosines of beam axis are calculated from beam geometry.
C**** Those for princ bending axis are stored in
C**** CE(ILT(nw,nr,nx)+1..,ne).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp40.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),NHE,nw
      REAL*8 CE(NMM),CG(NMM,NGM),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     '  WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mb,mh,mhg,mhs,mhsg,mk,mn,ms,MSTOT,mvar,
     '  nb,ng,nh,nhg,nhs,nhsg,nhx,nj,nk,NLOAD,nm,nn,nr,ns,NSTOT,nvar,nx
      REAL*8 AREA,DENS,DEPTH,DX(3),DXIX(3,3),EAIZZ,GAKYY,GL(3,3),
     '  GU(3,3),POIS,RIYY,RIZZ,RT(6,6),RWG,SMOD,SUM,
     '  SUM1,SUM2,SUM3,THICK,WIDTH,XX,XX2,YMOD

      CALL ENTERS('BEAM',*9999)
      nr=1   !Temporary
      nx=1   !Temporary

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' NVE='',I1,'' NHE='',I1)')
     '    NVE(nw),NHE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' CG(nm,1): '',14F6.2)') (CG(nm,1),nm=1,10)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

C *** Calculate rotation matrix RT
      SUM=0.D0
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        nb=NBJ(nj)
        DX(nj)=XE(1+(NNT(nb)-1)*NKT(0,nb),nj)-XE(1,nj)
        SUM=SUM+DX(nj)*DX(nj)
      ENDDO
      XX=DSQRT(SUM) !is length of element
      XX2=SUM      !is length of element squared
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Element length = '',E11.4)') XX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        RT(1,nj)=DX(nj)/XX !dir.n cosines of beam axis
      ENDDO
      IF(NJ_LOC(NJL_GEOM,0,nr).EQ.1) THEN
        RT(1,1)=1.D0
        RT(1,2)=0.D0
        RT(2,1)=0.D0
        RT(2,2)=1.D0
        IF(NHE.EQ.1) THEN
          RT(2,1)=1.D0
        ELSE IF(NHE.EQ.2) THEN
          RT(3,1)=1.D0
          RT(4,2)=1.D0
          RT(3,2)=0.D0
          RT(4,1)=0.D0
        ENDIF
      ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
        RT(2,1)=-RT(1,2)
        RT(2,2)= RT(1,1)
      ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
        DO nj=1,3
          RT(2,nj)=CE(ILT(nw,nr,nx)+nj)
          RT(3,nj)=CE(ILT(nw,nr,nx)+nj+3)
        ENDDO
      ENDIF
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Beam axis   direction cosines:'',3F7.4)')
     '    (RT(1,nj),nj=1,NJ_LOC(NJL_GEOM,0,nr))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' 1st transv. axis dir. cosines:'',3F7.4)')
     '    (RT(2,nj),nj=1,NJ_LOC(NJL_GEOM,0,nr))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' 2nd transv. axis dir. cosines:'',3F7.4)')
     '    (RT(3,nj),nj=1,NJ_LOC(NJL_GEOM,0,nr))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      DO ng=1,NGT(NBJ(1)) !To determine RG(ng)
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
        CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXIX,GL,GU,
     '    RG(ng),XG,ERROR,*9999)
      ENDDO
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' RG(ng): '',9E10.3)')
     '    (RG(ng),ng=1,NGT(NBJ(1)))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      NLOAD=NMP(nw)+NGP(nw) !is total # of material & geometric params
      mhs=0
      DO mvar=1,NVE(nw) !loops over local variables
        mh=NHV(mvar,nw) !is global variable corr to local variable mvar
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '       '(/'' Local variable mvar='',I2,'
     '       //'''(is global variable mh='',I2,'')'')') mvar,mh
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        mb=NBH(mh)
       ms=0
        DO mn=1,NNT(mb)
          DO mk=1,NKT(mn,mb)
           ms=ms+1
            mhs=mhs+1

C **        Calculate RHS loading term
            SUM=0.D0
            DO ng=1,NGT(mb)
              RWG=RG(ng)*WG(ng,mb)
              SUM=SUM+CG(NLOAD+mh,ng)*PG(ms,1,ng,mb)*RWG
            ENDDO
            SUM=SUM*SE(ms,mb)
            MSTOT=0
            DO mhg=1,NHE
              mhsg=ms+MSTOT
              MSTOT=MSTOT+NST(NBH(NH_LOC(1,nx)))
              ER(mhsg)=ER(mhsg)+RT(mh,mhg)*SUM
            ENDDO

            nhs=0
            DO nvar=1,NVE(nw)
              nhx=NHV(nvar,nw)
              nh=NH_LOC(nhx,nx)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,
     '             '(/'' Local variable nvar='',I2,''(is global '
     '             //'variable nh='',I2,'')'')') nvar,nh
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              nb=NBH(nh)
              ns=0
              DO nn=1,NNT(nb)
                DO nk=1,NKT(nn,nb)
                  ns=ns+1
                  nhs=nhs+1
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,
     '                 '('' mn='',I2,'' mk='',I2,'' mb='',I2,'
     '                //''' nn='',I2,'' nk='',I2,'' nb='',I2)')
     '                mn,mk,mb,nn,nk,nb
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF
                  SUM1=0.D0
                  SUM2=0.D0
                  SUM3=0.D0
                  DO ng=1,NGT(nb)
                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'('' Gauss point ng='',I2)') ng
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF
                    RWG=RG(ng)*WG(ng,nb)
                    YMOD =CG(1,ng) !is Young's modulus (GPa)
                    POIS =CG(2,ng) !is Poisson's ratio
                    DEPTH=CG(3,ng) !is depth of beam
                    WIDTH=CG(4,ng) !is width of beam

                    IF(KTYP45.EQ.1) THEN      !rectangular solid
                      AREA=WIDTH*DEPTH
                      RIZZ=WIDTH*DEPTH**3/12.D0 !is 2nd mom of area about z-axis
                      RIYY=DEPTH*WIDTH**3/12.D0 !is 2nd mom of area about y-axis
                    ELSE IF(KTYP45.EQ.2) THEN !rectangular tube
                      THICK=CG(5,ng)           !is wall thickness
                      AREA=WIDTH*DEPTH-(WIDTH-2.0D0*THICK)*
     '                  (DEPTH-2.0D0*THICK)
                      RIZZ=WIDTH*DEPTH**3/12.D0 !is 2nd mom of area about z-axis
     '                  -(WIDTH-2.0D0*THICK)*(DEPTH-2.0D0*THICK)**3/
     '                  12.D0
                      RIYY=DEPTH*WIDTH**3/12.D0 !is 2nd mom of area about y-axis
     '                  -(DEPTH-2.0D0*THICK)*(WIDTH-2.0D0*THICK)**3/
     '                  12.D0
                    ELSE IF(KTYP45.EQ.3) THEN !ellipsoidal solid
                      AREA=PI*WIDTH*DEPTH/4.D0
                      RIZZ=PI*WIDTH*DEPTH**3/64.D0 !is 2nd mom about z-axis
                      RIYY=PI*DEPTH*WIDTH**3/64.D0 !is 2nd mom about y-axis
                    ELSE IF(KTYP45.EQ.4) THEN !ellipsoidal tube
                      THICK=CG(5,ng)           !is wall thickness
                      AREA=PI*(WIDTH*DEPTH-(WIDTH-2.0D0*THICK)*
     '                  (DEPTH-2.0D0*THICK))/4.D0
                      RIZZ=PI*WIDTH*DEPTH**3/64.D0 !is 2nd mom about z-axis
     '                  -PI*(WIDTH-2.0D0*THICK)*(DEPTH-2.0D0*THICK)**3/
     '                  64.D0
                      RIYY=PI*DEPTH*WIDTH**3/64.D0 !is 2nd mom about y-axis
     '                  -PI*(DEPTH-2.0D0*THICK)*(WIDTH-2.0D0*THICK)**3/
     '                  64.D0
                    ENDIF

                    IF(ITYP5(nr,nx).EQ.1) THEN      !static analysis
                      DENS=0.D0
                    ELSE IF(ITYP5(nr,nx).GT.1) THEN !time dependent
                      DENS=CG(ILT(3,nr,nx),ng)*1.0D-3*AREA !MN/m
                    ENDIF

                    SMOD=YMOD/(2.D0*(1.D0+POIS)) !is shear modulus (GPa)
                    EAIZZ=YMOD*AREA*RIZZ  !is bending rigidity (z-axis)
c                    GAKZZ=SMOD*AREA*RIZZ  !is bending rigidity (z-axis)
                    GAKYY=SMOD*AREA*RIYY  !is bending rigidity (y-axis)

                    IF(mh.EQ.1.AND.nh.EQ.1) THEN
                      IF(mk.EQ.1.AND.nk.EQ.1) THEN
C **                  axial stiffness (assume linear basis)
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          WRITE(OP_STRING,
     '                      '('' *** Axial stiffness (AREA='','
     '                      //'E12.4,'')'')') AREA
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF
                        IF(mn.eq.nn) THEN !prod of PG derivs is +1
                          SUM1=SUM1+YMOD*AREA*RWG/XX2
                        ELSE              !prod of PG derivs is -1
                          SUM1=SUM1-YMOD*AREA*RWG/XX2
                        ENDIF
                        SUM2=0.D0
                        SUM3=SUM3+DENS*PG(ms,1,ng,mb)*PG(ns,1,ng,nb)
     '                                *RWG

                      ELSE IF(mk.EQ.2.AND.nk.EQ.2) THEN
C **                  torsional stiffness
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          WRITE(OP_STRING,
     '                      '('' *** Torsional stiffness'')')
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF

                      ELSE
C **                  Axial-Torsional coupling
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          WRITE(OP_STRING,
     '                      '('' *** Axial-Torsional coupling'')')
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF
                      ENDIF

                    ELSE IF(mh.EQ.2.AND.nh.EQ.2) THEN
C **                  bending stiffness about z-axis
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' *** Bending stiffness '
     '                    //'about z-axis *** (EAIZZ='',E12.4,'')'')')
     '                    EAIZZ
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF
                      SUM1=SUM1+EAIZZ*PG(ms,3,ng,mb)*PG(ns,3,ng,nb)
     '                  *RWG/XX2
c                     IF(mk.EQ.1.AND.nk.EQ.1) THEN      !(displ. term)
c                       SUM1=SUM1
c    '                      +GAKZZ*PG(ms,2,ng,mb)*PG(ns,2,ng,nb)*RWG/XX2
c                       SUM3=SUM3
c                     ELSE IF(mk.EQ.2.AND.nk.EQ.2) THEN !(rotation term)
c                       SUM1=SUM1+YMOD*CG(NMP(nw)+4,ng)*PG(ms,2,ng,mb)
c    '                    *PG(ns,2,ng,nb)*RWG/XX2
c    '                    +GAKZZ*PG(ms,1,ng,MB)*PG(ns,2,ng,nb)*RWG
c                       SUM3=SUM3
c                     ELSE IF(mk.EQ.1.AND.nk.EQ.2) THEN !(off-diag term)
c                       SUM1=SUM1-GAKZZ*PG(ms,2,ng,mb)*PG(ns,1,ng,nb)*RWG/XX
c                       SUM3=SUM3
c                     ELSE IF(mk.EQ.2.AND.nk.EQ.1) THEN !(off-diagonal term)
c                       SUM1=SUM1-GAKZZ*PG(ms,1,ng,mb)*PG(ns,2,ng,nb)*RWG/XX
c                       SUM3=SUM3
c                     ENDIF

                    ELSE IF(mh.EQ.3.AND.nh.EQ.3) THEN
C **                  bending stiffness about y-axis
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' *** Bending stiffness '
     '                    //'about y-axis *** (GAKYY='',E12.4,'')'')')
     '                    GAKYY
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF
                      IF(mk.EQ.1.AND.nk.EQ.1) THEN      !(displ. term)
                        SUM1=SUM1
     '                      +GAKYY*PG(ms,2,ng,mb)*PG(ns,2,ng,nb)*RWG/XX2
                        SUM3=SUM3
                      ELSE IF(mk.EQ.2.AND.nk.EQ.2) THEN !(rotation term)
                        SUM1=SUM1+YMOD*CG(NMP(nw)+4,ng)*PG(ms,2,ng,mb)
     '                    *PG(ns,2,ng,nb)*RWG/XX2
     '                    +GAKYY*PG(ms,1,ng,mb)*PG(ns,2,ng,nb)*RWG
                        SUM3=SUM3
                      ELSE IF(mk.EQ.1.AND.nk.EQ.2) THEN !(off-diag term)
                        SUM1=SUM1
     '                      -GAKYY*PG(ms,2,ng,mb)*PG(ns,1,ng,nb)*RWG/XX
                        SUM3=SUM3
                      ELSE IF(mk.EQ.2.AND.nk.EQ.1) THEN !(off-diag term)
                        SUM1=SUM1
     '                      -GAKYY*PG(ms,1,ng,mb)*PG(ns,2,ng,nb)*RWG/XX
                        SUM3=SUM3
                      ENDIF
                    ENDIF
                  ENDDO
                  SUM1=SUM1*SE(ms,mb)*SE(ns,nb)
                  SUM2=SUM2*SE(ms,mb)*SE(ns,nb)
                  SUM3=SUM3*SE(ms,mb)*SE(ns,nb)
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,
     '                '('' SUM1='',E11.3,'' SUM2='',E11.3,'
     '                //''' SUM3='',E11.3)') SUM1,SUM2,SUM3
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF
                  MSTOT=0
                  DO mhg=1,NHE
                    mhsg=ms+MSTOT
                    MSTOT=MSTOT+NST(NBH(mhg))
                    NSTOT=0
                    DO nhg=1,NHE
                      nhsg=ns+NSTOT
                      NSTOT=NSTOT+NST(NBH(nhg))
                      IF(ITYP5(nr,nx).EQ.1) THEN !static
                        ES(mhsg,nhsg)=ES(mhsg,nhsg)+RT(mh,mhg)*SUM1
     '                    *RT(nh,nhg)
                      ELSE IF(ITYP5(nr,nx).EQ.2.OR.ITYP5(nr,nx).EQ.3)
     '                  THEN
                        !time int. or modal
                        ES(mhsg,nhsg)=ES(mhsg,nhsg)+RT(mh,mhg)*
     '                    (SUM1+CG(NLOAD+1,ng)*SUM2)*RT(nh,nhg)
                        EM(mhsg,nhsg)=EM(mhsg,nhsg)+RT(mh,mhg)*SUM3
     '                    *RT(nh,nhg)
c cpb 1/5/95 Replacing Fourier analysis with Quasi-static analysis
c                      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Fourier
                      ELSE IF(ITYP5(nr,nx).EQ.6) THEN !buckling
                        ES(mhsg,nhsg)=ES(mhsg,nhsg)+RT(mh,mhg)*SUM1
     '                    *RT(nh,nhg)
                        EM(mhsg,nhsg)=EM(mhsg,nhsg)+RT(mh,mhg)*SUM2
     '                    *RT(nh,nhg)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('BEAM')
      RETURN
 9999 CALL ERRORS('BEAM',ERROR)
      CALL EXITS('BEAM')
      RETURN 1
      END


