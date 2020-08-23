      SUBROUTINE ZERE40(NBH,NBJ,ne,NHE,nw,nx,CE,CG,PG,RE,
     '  SE,WG,XE,XG,ZE,ZG,ERROR,*)

C#### Subroutine: ZERE40
C###  Description:
C###    ZERE40 calculates element residual RE from current dependent
C###    variable array ZE.

C**** X variables refer to orthog curvilinear coords in reference state.
C**** Z     "       "         "        "         "      deformed    "
C**** Material Theta-coordinates (reference for deformation) coincide
C****   with Xj-coords in ref state.
C**** Material Xi-coords are the finite element mesh coordinates.
C**** Material Nu-coordinates (reference for stresses): are orthogonal
C****   and  if 3D problem (Nu2,Nu3) lie in the (Xi1-Xi2) plane  s.t.
C****   Nu3 is aligned with the 'fibres' to which material aeolotropy
C****   is referred;
C****   else if 2D problem (Nu1,Nu2) lie in the (Xi1-Xi2) plane and
C****   Nu3 is normal to this plane.
C****   Also: The base vectors are defined such that the undeformed
C****   metric tensors wrt the Nu & Theta systems are equal.
C****   Fibre direction field must be defined for Nu-coords.
C**** AXL,AXU,AZL & AZU are metric tensor cpts wrt Nu-coords.
C**** GXL,GXU are undeformed   "      "     "   "  Xi-coords.
C**** For membranes CG(4,ng) is transverse load in kPa if isotropic
C****  "     "      CG(6,ng) is transverse load in kPa if orthotropic
C**** For beams     CG(5,ng) is transverse load in kN/m.
C**** For battens   CG(4,ng) is axial load in kN/m.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'press00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),ne,NHE,nw,nx
      REAL*8 CE(NMM),CG(NMM,NGM),PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),
     '  SE(NSM,NBFM,NEM),WG(NGM,NBM),XE(NSM,NJM),
     '  XG(NJM,NUM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ng,nh,nhx,ni,nj,nr,ns,NU1(0:3)
      REAL*8 AXL(3,3),AXIAL_STRAIN,AZ,AZL(3,3),DS0,DS1,DS2,
     '  DSDXI,
     '  DX1,DX2,DXIX(3,3),DXINU(3,3),EAU,EAUV,EG(3,3),EIV,EIV1,EIV2,
     '  FX,FY,FYX,FZ,GXL(3,3),GXU(3,3),PGX,PLOAD,PPG(9,2,16),PV,
     '  RGX,RGZ,RL,RLENGTH,RWG,
     '  STRAIN,SZ,TENSION,TG(3,3),XN_LOCAL(3)
      DATA NU1/1,2,4,7/

      CALL ENTERS('ZERE40',*9999)
      nr=1   !Temporary
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/'' ZERE40 *** element '',I3,'' ***'')')
     '    NE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' nw='',I3)') nw
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      DO nhx=1,NHE
        nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
        DO ns=1,NSM
          RE(ns,nh)=0.0d0
        ENDDO
      ENDDO
      NB  =NBH(Nh_LOC(1,nx))

      IF(nw.EQ.1) THEN
C ***   Cable or truss elements
C ***   Calculate tension (in kN)
        CALL ZETG41(CE,STRAIN,SZ,TENSION,XE,ZE,ERROR,*9999)

C ***   Calculate residuals from truss elasticity
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          RE(1,nj)=TENSION*(ZE(1,nj)-ZE(2,nj))/SZ
          RE(2,nj)=-RE(1,nj)
        ENDDO

      ELSE IF(nw.EQ.2) THEN
C ***   Batten elements
        DO ng=1,NGT(NBH(NH_LOC(2,nx)))
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
          CALL XGMG(0,1,NBJ(1),nr,DXIX,GXL,GXU,RGX,XG,
     '      ERROR,*9999)
          RWG=RGX*WG(ng,nb)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' RGX='',E11.3,'' RWG='',E11.3)')
     '        RGX,RWG
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

          CALL ZETG42(NBH,ng,CG(1,ng),DS0,DS1,DSDXI,EIV1,EIV2,PG,PV,
     '      RL,XE,ZE,ERROR,*9999)

C ***     Calculate residuals from length constraint
          nb=NBH(NH_LOC(1,nx))
          DO ns=1,NST(nb)
            RE(ns,1)=RE(ns,1)+RL*PG(ns,2,ng,nb)*DS1*RWG
          ENDDO

C ***     Calculate residuals from flexural rigidity
c         DS2=DS0*DS0
          DS2=DS1*DS1
          nb=NBH(NH_LOC(2,nx))
          DO ns=1,NST(nb)
            RE(ns,2)=RE(ns,2)+(EIV1*PG(ns,3,ng,nb)*DS2+
     '        (EIV2-PV)*PG(ns,2,ng,nb)*DS1)*RWG
          ENDDO

        ENDDO

      ELSE IF(nw.EQ.3) THEN
C ***   Beam elements
        DO ng=1,NGT(NBH(NH_LOC(1,nx)))

          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
          CALL XGMG(0,1,NBJ(1),nr,DXIX,GXL,GXU,RGX,XG,
     '      ERROR,*9999)
          RWG=RGX*WG(ng,nb)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' RGX='',E11.3,'' RWG='',E11.3)')
     '        RGX,RWG
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          DX1=DXIX(1,1)
          DX2=DX1*DX1

          CALL ZETG43(NBH,ng,AXIAL_STRAIN,CG(1,ng),DX2,EAU,EAUV,EIV,
     '      PG,XE,ZE,ERROR,*9999)

C ***     Calculate residuals from beam axial stiffness
          nb=NBH(NH_LOC(1,nx))
          DO ns=1,NST(nb)
            RE(ns,1)=RE(ns,1)+(EAU*PG(ns,2,ng,nb)*DX1
     '        -CG(4,ng)*PG(ns,1,ng,nb))*RWG
          ENDDO

C ***     Calculate residuals from beam flexural rigidity
          nb=NBH(NH_LOC(2,nx))
          DO ns=1,NST(nb)
            RE(ns,2)=RE(ns,2)+(EIV*PG(ns,3,ng,nb)*DX2
     '        +EAUV*PG(ns,2,ng,nb)*DX1-CG(5,ng)*PG(ns,1,ng,nb))*RWG
          ENDDO

        ENDDO

      ELSE IF(nw.EQ.4) THEN
C ***   Link elements
C ***   Calculate forces in X,Y & Z directions (in kN)
        CALL ZETG44(CE,FX,FY,FYX,FZ,ZE,ERROR,*9999)

C ***   Calculate residuals
        RE(1,1)=FX
        RE(1,2)=FY
        RE(2,2)=FYX
        RE(1,3)=FZ
        RE(2,1)=-FX
        RE(3,2)=-FY
        RE(4,2)=-FYX
        RE(2,3)=-FZ

      ELSE IF(nw.EQ.5) THEN
C ***   Membrane elements
        DO ng=1,NGT(NBH(NH_LOC(1,nx)))

C ***     Calculate DXINU: derivs of Xi wrt Nu (fibre) coords
C ***     (unless no fibres defined; then stresses are wrt Xi coords)
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
          CALL XGMG(1,2,NBJ(1),nr,DXINU,GXL,GXU,RGX,XG,
     '      ERROR,*9999)
          DO nb=1,NBFT
            IF(nb.EQ.NBH(NH_LOC(1,nx)).OR.nb.EQ.NBH(NH_LOC(2,nx))
     '        .OR.nb.EQ.NBH(NH_LOC(3,nx))) THEN
              DO ns=1,NST(nb)
                DO ni=1,2
                  PPG(nb,ni,ns)=PGX(nb,ni,ns,DXINU,PG(1,1,ng,nb))
                ENDDO
              ENDDO
            ENDIF
          ENDDO

C ***     Calc XG & ZG wrt Nu coords and metric tensors wrt Nu coords
          CALL ZEAZ45(NBJ,ng,nr,AXL,AZ,AZL,PG,PPG,
     '      XE,XG,ZE,ZG,ERROR,*9999)

          RGZ=DSQRT(AZ)
          nb=NBH(NH_LOC(1,nx))
          RWG=RGX*WG(ng,nb)
          IF(ng.EQ.1.AND.DOP) THEN
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              nb=NBJ(nj)
              WRITE(OP_STRING,'('' XE(ns,'',I1,''): '',9E12.4)')
     '          nj,(XE(ns,nj),ns=1,NST(nb))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
            WRITE(OP_STRING,'('' DXINU:'',4E11.3)')
     '        DXINU(1,1),DXINU(1,2),DXINU(2,1),DXINU(2,2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nhx=1,NHE
              nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
              nb=NBH(nh)
              WRITE(OP_STRING,'('' ZE(ns,'',I1,''): '',9E12.4)')
     '          nhx,(ZE(ns,nhx),ns=1,NST(nb))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
            DO nhx=1,NHE
              nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
              WRITE(OP_STRING,'('' XG('',I1,'',NU1(ni)): '',8E15.8)')
     '          nh,(XG(nh,NU1(ni)),ni=0,NIT(NBH(nh)))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
            DO nhx=1,NHE
              nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
              WRITE(OP_STRING,'('' ZG('',I1,'',NU1(ni)): '',8E15.8)')
     '          nh,(ZG(nh,NU1(ni)),ni=0,NIT(NBH(nh)))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
            WRITE(OP_STRING,'('' RGX='',E11.3,'' RGZ='',E11.3)')
     '        RGX,RGZ
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C ***     Calculate stress tensor wrt Nu coords (in kN/m)
          CALL AZTG45(nw,AXL,AZL,CG(1,ng),EG,TG,ERROR,*9999)

C ***     Calculate residuals from membrane elasticity
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            nb=NBH(nj)
            DO ns=1,NST(nb)
C             RE(ns,nj)=RE(ns,nj)+RWG*(PPG(nb,1,ns)*
C    '          (TG(1,1)*ZG(nj,NU1(1))+TG(1,2)*ZG(nj,NU1(2)))
C    '          +PPG(nb,2,ns)*
C    '          (TG(2,1)*ZG(nj,NU1(1))+TG(2,2)*ZG(nj,NU1(2))))
              RE(ns,nj)=RE(ns,nj)+RWG*(PPG(nb,1,ns)*
     '          (TG(1,1)*XG(nj,NU1(1))+TG(1,2)*XG(nj,NU1(2)))
     '          +PPG(nb,2,ns)*
     '          (TG(2,1)*XG(nj,NU1(1))+TG(2,2)*XG(nj,NU1(2))))
            ENDDO
          ENDDO

C ***     Calculate residuals from pressure on membrane
          PLOAD=PRESS(ng,ne)
          CALL NORM40(ZG,XN_LOCAL,ERROR,*9999)
          RLENGTH=DSQRT(XN_LOCAL(1)**2+XN_LOCAL(2)**2+XN_LOCAL(3)**2)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            nb=NBH(nj)
            XN_LOCAL(nj)=-XN_LOCAL(nj)/RLENGTH
            DO ns=1,NST(nb)
              RE(ns,nj)=RE(ns,nj)-PLOAD*XN_LOCAL(nj)*PG(ns,1,ng,nb)*RWG
            ENDDO
          ENDDO

        ENDDO
      ENDIF

      DO nhx=1,NHE
        nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
        nb=NBH(nh)
        DO ns=1,NST(nb)
          RE(ns,nh)=RE(ns,nh)*SE(ns,nb,ne)
        ENDDO
      ENDDO

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,NHE
          nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
          nb=NBH(nh)
          WRITE(OP_STRING,'('' RE(ns,'',I1,''): '',8E15.7)')
     '      nh,(RE(ns,nh),ns=1,NST(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZERE40')
      RETURN
 9999 CALL ERRORS('ZERE40',ERROR)
      CALL EXITS('ZERE40')
      RETURN 1
      END


