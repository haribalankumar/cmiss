      SUBROUTINE ZERE30(NBH,NBJ,NHE,nr,nx,
     '  CG,PG,RE,WG,XE,XG,YG,ZE,ZG,ERROR,*)

C#### Subroutine: ZERE30
C###  Description:
C###    ZERE30 calculates element residual RE for non-linear problems
C###    from current element dependent variable array ZE.

C**** (ityp5=5 ityp2=3)
C****   Elliptic eikonal equation for wavefront path determination.
C**** (ityp5=1 ityp2=5 ityp3=1)
C****   Poisson equation with radiation boundary condition.
C**** (ityp5=1 ityp2=8 ityp3=1)
C****   Nonlinear static Prandtl boundary layer equations:
C****     ZG(nh=1..nje,1) = velocities (m/s)
C****     XG(nj=nje+1..nje+nje,1) = free stream velocity (m/s)
C****     CG(1,ng) is Kinematic viscosity (m^2/s)
C****     CG(2,ng) is Fluid density (kg/m^3)
C**** (ityp5=1 ityp2=9 ityp3=1)
C****   Nonlinear static multi-field oxygen diffusion problem:
C****      ktyp15 =  number of fields
C****   Dependent variables are:
C****     ZG(1,1) = oxygen concentration in muscle (nmols/mm^3/ks)
C****     ZG(2,1) = oxygen concentration in blood  (nmols/mm^3/ks)
C****   CG(1,ng) is Oxy. m. solubility (nmol/kPa/mm^3)
C****      2        Oxy. b. solubility (nmol/kPa/mm^3)
C****      3        Oxy. m. diffusivity (mm^2/ks)
C****      4        Oxy. b. diffusivity (mm^2/ks)
C****      5        Myoglobin diffusivity (mm^2/ks)
C****      6        Max myoglobin   conc. (nmol/mm^3)
C****      7        Max haemoglobin conc. (nmol/mm^3)
C****      8        Max mitochon. rate (nmol/ks/mm^3)
C****      9        Exchange coefficient (1/ks)
C****     10        Arterial blood oxygen p.p. (kPa)
C****     11        p_50 for myoglobin (kPa)
C****     12        p_50 for haemoglobin (kPa)
C****     13        n for haemoglobin (kPa)
C****     14        p_50 for mitochondria (kPa)
C****     15        Capillary blood velocity   (mm/ks)
C****     16        Relative volume of vascular bed
C****   Conc_MB_Max is max conc. of myoglobin   in muscle (nmols/mm^3)
C****   Conc_HB_Max is max conc. of haemoglobin in blood  (nmols/mm^3)
C****   P50_MB & RN_MB are p_50 & n for myoglobin
C****   P50_HB & RN_HB are p_50 & n for haemoglobin
C****   P50_MC & RN_MC are p_50 & n for mitochondria
C****   F is myoglobin    binding function f; dF_dP is d(f)/d(p)
C****   H is haemoglobin  binding function h; dH_dP is d(h)/d(p)
C****   G is mitochondria binding function g
C**** (ityp5=1 ityp2=9 ityp3=2)
C****   Nonlinear static glucose-oxygen diffusion problem:
C****   Dependent variables are:
C****     ZG(1,1) = oxygen  concentration (nmols/mm^3/ks)
C****     ZG(2,1) = glucose concentration (nmols/mm^3/ks)
C****   CG(1,ng) is Oxygen  solubility (nmol/kPa/mm^3)
C****      2        Glucose solubility (nmol/kPa/mm^3)
C****      3        Oxygen  diffusivity (mm^2/ks)
C****      4        Glucose diffusivity (mm^2/ks)
C****      5        Oxygen  Michaelis p.pressure (kPa)
C****      6        Glucose Michaelis p.pressure (kPa)
C****      7        Vmax (nmol/ks/mm^3)
C****      8        Oxygen reaction weight (nu)

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),NHE,nr,nx
      REAL*8 CG(NMM,NGM),PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),
     '  WG(NGM,NBM),XE(NSM,NJM),
     '  XG(NJM,NUM),YG(NIYGM,NGM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER JP,mi,nb,nbh_e,ng,nh,nhx,ni,NITB,niyg,nj,njj,NJTR,ns,NSTB,
     '  NU1(0:3),NU2(3,3)
      INTEGER*4 NULL
      REAL*8 AA,ADV,CONC_HB_MAX,CONC_MB_MAX,CONT0(2),
     '  dB_dX(3),dF_dP,dF_dX(3),dGL_dX(3),dH_dP,dH_dX(3),
     '  DIFF,DIFF_B,DIFF_M,DIFF_MB,DIFFM1(3),DIFFMODM(3),
     '  dM_dX(3),DNORM,dO2_dX(3),dPSI_dX(3),dT_dX(3),
     '  dU_dX(3,3),dV_dX(3,3),DXIX(3,3),DZXI(3),
     '  EXCH_COEFF,G,G11,G33,GL(3,3),GU(3,3),GUxDZ(3),MGAL,MPET,
     '  P50_HB,P50_MB,P50_MC,PGL,PGN,PGX,PO2,PO2_B,PO2_M,RAD,
     '  RATE_MC_MAX,RATIO,RATIO_N,REL_VOL,RES,RG,RN_HB,RN_MB,RN_MC,RWG,
     '  SLX,SMX,SOLU_B,SOLU_M,SUM,SUM1,SUM2,SUM3,SUM4,TEMP,TEMP2,U(3),
     '  V(3),VELOC,WGHTM1(3)
      LOGICAL MATERIALPET,MODDIFF,PETROV
!     External Functions
      REAL*8 DDOT

      PARAMETER (NULL=0)
      PARAMETER (RN_MB=1.d0,RN_MC=1.d0)

      DATA NU1/1,2,4,7/
      DATA NU2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('ZERE30',*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' NHE='',I2,'' NBH(1)='',I2)') NHE,NBH(1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      DO nhx=1,NHE
        nh=nh_loc(nhx,nx)
        DO ns=1,NST(NBH(nh))+NAT(NBH(nh))
          RE(ns,nh)=0.0d0
        ENDDO
      ENDDO

      NITB=NIT(NBJ(1))
      NJTR=NJ_LOC(NJL_GEOM,0,nr)

      IF(ITYP5(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.3) THEN !eikonal equation
        MATERIALPET=ITYP15(nr,nx).EQ.1
        PETROV=MATERIALPET.OR.ITYP15(nr,nx).EQ.2
C       The last field may contain a diffusion modifier that
C       approaches zero at singular initiation points.
        njj=NJ_LOC(NJL_FIEL,0,nr)
        MODDIFF=njj.EQ.3
        IF(MODDIFF) THEN
          nj=NJ_LOC(NJL_FIEL,njj,nr)
          MODDIFF=NBJ(nj).NE.0
        ENDIF
      ENDIF !eikonal equation

      DO ng=1,NGT(NBH(NH_LOC(1,nx)))

        IF(ITYP5(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.3) THEN !eikonal equation
          IF(PETROV) THEN !Petrov-Galerkin
            njj=1
C           This field should contain a C0 continuous variable that is
C           1 over most of the domain but approaches zero at Dirichlet
C           boundary conditions.
            CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).GE.njj,
     '        '>>Boundary condition field not defined',ERROR,*9999)
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            nb=NBJ(nj)
            IF(nb.NE.0) THEN
C             Variable
              CONT0(1)=DDOT(NST(nb)+NAT(nb),PG(1,1,ng,nb),1,XE(1,nj),1)
C             Multiplier for unit derivatives
              IF(NKT(0,nb).EQ.1) THEN !linear element
                CONT0(2)=1.0d0
              ELSE !cubic element
                CONT0(2)=1.0d0/3.0d0
              ENDIF
            ELSE
              CONT0(1)=1.0d0
              CONT0(2)=1.0d0
            ENDIF
            JP=2 !need 2nd derivs in ZG
          ELSE
            JP=0 !calculate only 1st xi derivs in ZEZG
          ENDIF
C***      Get information contained in YG
          niyg=1
          RWG=YG(niyg,ng)
C         Set up coupling tensor in xi coords GU from YG.
          DO mi=1,NITB
            DO ni=1,mi
              niyg=niyg+1
              GU(mi,ni)=YG(niyg,ng)
              IF(mi.NE.ni) THEN
                GU(ni,mi)=GU(mi,ni)
              ENDIF
            ENDDO !ni
          ENDDO !mi
          IF(PETROV) THEN
C           Multiplier of first derivatives for diffusion term.
            DO ni=1,NITB
              niyg=niyg+1
              DIFFM1(ni)=YG(niyg,ng)
            ENDDO !ni
          ENDIF !PETROV
          IF(MODDIFF) THEN
C           Multiplier of first derivatives for diffusion modifier.
            DO ni=1,NITB
              niyg=niyg+1
              DIFFMODM(ni)=YG(niyg,ng)
            ENDDO !ni
          ENDIF !MODDIFF
        ELSE
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
C         DXIX is wrt Xj-coords
          CALL XGMG(0,NITB,NBJ(1),nr,DXIX,GL,GU,RG,XG,ERROR,*9999)
          IF(ITYP10(1).EQ.4) THEN
            AA=FOCUS*FOCUS
            SLX=DSINH(XG(1,1)) !sinh(lamda)
            SMX=DSIN (XG(2,1)) !sin(mu)
            G11=1.0d0/(AA*(SLX*SLX+SMX*SMX))
            G33=1.0d0/(AA* SLX*SLX*SMX*SMX )
          ENDIF
          JP=1 !calculate only 1st Xj derivs in ZEZG
          RWG=RG*WG(ng,NBH(NH_LOC(1,nx)))
        ENDIF
        CALL ZEZG(JP,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
C       Jacobian * Gauss weight
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '      '('' ng='',I2,'' RG='',E11.4,'' WG='',E11.4,'' XG:'','
     '      //'8E11.4)') ng,RG,WG(ng,1),XG(1,1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' RWG='',E11.4,'' ZG:'',4E11.4)')
     '      RWG,(ZG(1,NU1(nj)),nj=0,NJTR)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(JTYP4.EQ.2) THEN
          IF(ITYP10(1).EQ.1) THEN
            RAD=XG(2,1)
          ELSE IF(ITYP10(1).EQ.2.OR.ITYP10(1).EQ.3) THEN
          ELSE IF(ITYP10(1).EQ.4) THEN
          ELSE IF(ITYP10(1).EQ.5) THEN
          ENDIF
          RWG=RWG*2.0d0*PI*RAD
        ELSE IF(JTYP4.EQ.3) THEN
          IF(ITYP10(1).EQ.1) THEN
            RAD=XG(1,1)
          ELSE IF(ITYP10(1).EQ.2.OR.ITYP10(1).EQ.3) THEN
          ELSE IF(ITYP10(1).EQ.4) THEN
          ELSE IF(ITYP10(1).EQ.5) THEN
          ENDIF
          RWG=RWG*2.0d0*PI*RAD
        ELSE IF(JTYP4.EQ.4) THEN
          RWG=RWG*4.0d0*PI*XG(1,1)**2
        ENDIF

        IF(ITYP5(nr,nx).EQ.5) THEN !Wavefront path analysis
          IF(ITYP2(nr,nx).EQ.3) THEN
C*** Elliptic eikonal equation for wavefront path determination.
C***   CG(1,ng) is the time constant source term,
C***   CG(2,ng) is the dimensionless coefficient of advection,
C***   CG(3..,ng) are coupling coefficients.
C***   FACTOR is the continuation parameter.
            nhx=1
            nh=NH_LOC(nhx,nx)
            nbh_e=NBH(nh)
            NSTB=NST(nbh_e)

C           Copy first derivatives of ZG wrt xi into grad(time).
            DO ni=1,NITB
              DZXI(ni)=ZG(1,NU1(ni))
            ENDDO
C           Find temporary vector of coupling and grad(time), GUxDZ(mi).
            DO mi=1,NITB
              SUM=0.0d0
              DO ni=1,NITB
                SUM=SUM+DZXI(ni)*GU(ni,mi)
              ENDDO
              GUxDZ(mi)=SUM
            ENDDO
C           Coupling norm of first derivatives.
            DNORM=0.0d0
            DO ni=1,NITB
              DNORM=DNORM+GUxDZ(ni)*DZXI(ni)
            ENDDO
C***        Advection and Source Terms:
            VELOC=CG(2,ng)*DSQRT(DNORM)
            ADV=FACTOR*(VELOC-CG(1,ng))*RWG
C            VELOC=DSQRT(CG(2,ng)**2*DNORM+CG(1,ng)**2)
C            ADV=FACTOR*(VELOC-1.4142136d0*CG(1,ng))*RWG
C           Diffusion modification near singularities.
            IF(MODDIFF) THEN
              SUM=0.0d0
              DO ni=1,NITB
                SUM=SUM+DIFFMODM(ni)*DZXI(ni)
              ENDDO
              ADV=ADV-SUM
            ENDIF

            IF(PETROV) THEN !Petrov-Galerkin

C***          Diffusion Term:
              DIFF=0.0d0
              DO ni=1,NITB
                DO mi=1,NITB
                  DIFF=DIFF+GU(mi,ni)*ZG(1,NU2(mi,ni))
                ENDDO
              ENDDO
              DIFF=DIFF*RWG
              DO ni=1,NITB
                DIFF=DIFF+DIFFM1(ni)*DZXI(ni)
              ENDDO

C***          Residual
              RES=(ADV-DIFF)

              CALL PETROV_PREP(NITB,CG(1,ng),CONT0,DNORM,DZXI,GU,GUxDZ,
     '          %VAL(NULL),%VAL(NULL),%VAL(NULL),
     '          MGAL,MPET,MATERIALPET,.FALSE.,ERROR,*9999)

C***          Make additions to the residual integrals.
C             Part of Galerkin weights.
              CALL DAXPY(NSTB,MGAL*ADV,PG(1,1,ng,nbh_e),1,RE(1,nh),1)
C             Petrov derivatives and part of Galerkin weights.
              IF(MATERIALPET) THEN
C***            For material coord based derivative weights:
                TEMP=MPET*RES+MGAL*RWG
                DO ni=1,NITB
                  WGHTM1(ni)=TEMP*GUxDZ(ni)
                ENDDO
              ELSE
C***            For element local coord based derivative weights:
                TEMP=MPET*RES
                TEMP2=MGAL*RWG
                DO ni=1,NITB
                  WGHTM1(ni)=TEMP*DZXI(ni)+TEMP2*GUxDZ(ni)
                ENDDO
              ENDIF !materialpet
              DO ni=1,NITB
                CALL DAXPY(NSTB,WGHTM1(ni),
     '            PG(1,NU1(ni),ng,nbh_e),1,RE(1,nh),1)
              ENDDO

            ELSE !Galerkin

C             Set up temporary vector U(ni) to include everything for
C             diffusion except weighting function derivative wrt xi.
              DO ni=1,NITB
                U(ni)=GUxDZ(ni)*RWG
              ENDDO

C             Add diffusion to residual
              DO ni=1,NITB
                CALL DAXPY(NSTB,U(ni),PG(1,NU1(ni),ng,nbh_e),1,
     '            RE(1,nh),1)
              ENDDO
C             Add advection and source term to residual
              CALL DAXPY(NSTB,ADV,PG(1,1,ng,nbh_e),1,RE(1,nh),1)

            ENDIF
          ELSE
            ERROR='>>Only elliptic eikonal equation is implemented'
            GO TO 9999
          ENDIF

        ELSE IF(ITYP2(nr,nx).EQ.5.AND.ITYP3(nr,nx).EQ.1) THEN
C ***     Poisson equation with radiation bdry condition
          nb=NBH(NH_LOC(1,nx))
          DO nj=1,NJTR
            dT_dX(nj)=ZG(1,NU1(nj)) !are cpts of temperature grad
          ENDDO
          DO ns=1,NST(nb)
            PGN=PG(ns,1,ng,nb)
            DO nj=1,NJTR!are cpts of grad(psi)
              dPSI_dX(nj)=PGX(nb,nj,ns,DXIX,PG(1,1,ng,nb))
            ENDDO
            IF(ITYP10(1).EQ.1) THEN !rectangular cartesian coordinates
              SUM1=0.0d0
              DO nj=1,NJTR
                SUM1=SUM1+CG(1+nj,ng)*dT_dX(nj)*dPSI_dX(nj)
              ENDDO
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' SUM1='',E11.3)') SUM1
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ELSE
              ERROR='>>Only rectangular Cartesian coordinates '
     '          //'are implemented yet'
              GO TO 9999
            ENDIF
            RE(ns,1)=RE(ns,1)+SUM1*RWG
          ENDDO !ns

        ELSE IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.1) THEN
C ***     Prandtl boundary layer equations
          DO nhx=1,NHE
            nh=nh_loc(nhx,nx)
            V(nh)=ZG(nh,1)     !are bdry layer velocs
            U(nh)=XG(NJTR+nh,1) !are free stream velocs
          ENDDO
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' Bdry layer velocs:'',3E12.3)')
     '        (V(nh_loc(nhx,nx)),nhx=1,NHE)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Free stream velocs:'',3E12.3)')
     '        (U(nh_loc(nhx,nx)),nhx=1,NHE)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

C         Calculate gradient vector of velocs
          DO nj=1,NJTR
            DO nhx=1,NHE
              nh=nh_loc(nhx,nx)
              dV_dX(nh,nj)=ZG(nh,NU1(nj))
              dU_dX(nh,nj)=XG(NJTR+nh,NU1(nj))
            ENDDO
          ENDDO

C         Calculate element residuals for momentum equations
          DO nhx=1,NJTR-1
            nh=nh_loc(nhx,nx)
            nb=NBH(nh)
            DO ns=1,NST(nb)
              PGN=PG(ns,1,ng,nb)
              DO nj=1,NJTR !are cpts of grad(psi)
                dPSI_dX(nj)=PGX(nb,nj,ns,DXIX,PG(1,1,ng,nb))
              ENDDO

              IF(ITYP10(1).EQ.1) THEN !rectangular cartesian coordinates
                SUM1=0.0d0
                SUM2=0.0d0
                DO nj=1,NJTR
                  SUM1=SUM1+V(nj)*dV_dX(nh,nj)       !convective term
                  SUM2=SUM2+dV_dX(nh,nj)*dPSI_dX(nj) !viscous term
                ENDDO
                SUM1=SUM1*PGN
                SUM2=SUM2*CG(1,ng)      !CG(1) is kinematic viscosity
                SUM3=U(1)*dU_dX(1,nh)*PGN        !Prandtl press term
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,
     '              '('' SUM1='',E11.3,'' SUM2='',E11.3,'
     '              //''' SUM3='',E11.3)') SUM1,SUM2,SUM3
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
              ELSE
                ERROR='>>Other coordinate systems not implemented yet'
                GO TO 9999
              ENDIF
              RE(ns,nh)=RE(ns,nh)+(SUM1+SUM2-SUM3)*RWG
            ENDDO !ns
          ENDDO !nh

C         Calculate element residuals for conservation of mass equation
          nb=NBH(NJTR) !nb for v velocity
          DO ns=1,NST(nb)
            PGN=PG(ns,1,ng,nb)

            IF(ITYP10(1).EQ.1) THEN !rectangular cartesian coordinates
              SUM1=0.0d0
              DO nj=1,NJTR
                SUM1=SUM1+dV_dX(nj,nj) !is divergence of velocity
              ENDDO
              SUM1=SUM1*PGN
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' SUM1='',E11.3)') SUM1
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ELSE
              ERROR='>>Other coordinate systems not implemented yet'
              GO TO 9999
            ENDIF
            RE(ns,NJTR)=
     '        RE(ns,NJTR)+SUM1*RWG
          ENDDO

        ELSE IF(ITYP2(nr,nx).EQ.9.AND.ITYP3(nr,nx).EQ.1) THEN
C ***     Oxygen transport equation (single or multi-field)
          Solu_M =CG( 1,ng)
          Solu_B =CG( 2,ng)
          Diff_M =CG( 3,ng)
          Diff_B =CG( 4,ng)
          Diff_MB=CG( 5,ng)
          Conc_MB_Max=CG( 6,ng)
          Conc_HB_Max=CG( 7,ng)
          Rate_MC_Max=CG( 8,ng)
          Exch_coeff =CG( 9,ng)
          P50_MB =CG(11,ng)
          P50_HB =CG(12,ng)
          RN_HB  =CG(13,ng)
          P50_MC =CG(14,ng)
          Veloc  =CG(15,ng)
          Rel_vol=CG(16,ng)
          PO2_M =ZG(1,1) !is oxygen partial pressure in muscle
          IF(po2_m.le.0.0d0) PO2_M=1.0d0 !to avoid -ve pressure
          IF(KTYP15.EQ.1) THEN !one field
            PO2_B=CG(10,ng)
          ELSE                 !two field
            PO2_B=ZG(2,1) !is oxygen partial pressure in blood
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' PO2_M='',E12.3,'' PO2_B='',E12.3,'
     '        //''' P50_MB='',E12.3,'' P50_MC='',E12.3)')
     '        PO2_M,PO2_B,P50_MB,P50_MC
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

          !Calculate gradient vector of oxygen partial pressure in muscle
          DO nj=1,NJTR
            dM_dX(nj)=ZG(1,NU1(nj))
          ENDDO

          IF(KTYP15.GT.1) THEN !multi-field problem
            !Calculate gradient vector of oxygen partial pressure in blood
            DO nj=1,NJTR
              dB_dX(nj)=ZG(2,NU1(nj))
            ENDDO
          ENDIF

          !Calculate gradient vector of myoglobin binding function f
          IF(DABS(PO2_M).gt.1.d-5) THEN !pO2 is gt zero
            RATIO=PO2_M/P50_MB
            RATIO_N=RATIO**RN_MB
            dF_dP=Conc_MB_Max*(RN_MB*RATIO_N/(1.d0+RATIO_N)**2)/PO2_M
          ELSE                         !pO2 is zero
            dF_dP=0.0d0
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' dF_dP='',E12.3)') dF_dP
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          DO nj=1,NJTR
            dF_dX(nj)=dF_dP*ZG(1,NU1(nj)) !is cpt of grad(f)
          ENDDO

          IF(KTYP15.GT.1) THEN !multi-field problem
            !Calculate gradient vector of haemoglobin binding function h
            IF(DABS(PO2_B).gt.1.d-5) THEN !pO2 is gt zero
              RATIO=PO2_B/P50_HB
              RATIO_N=RATIO**RN_HB
              dH_dP=Conc_HB_Max*(RN_HB*RATIO_N/(1.d0+RATIO_N)**2)
     '          /PO2_B
            ELSE                         !pO2 is zero
              dH_dP=0.0d0
            ENDIF
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' dH_dP='',E12.3)') dH_dP
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            DO nj=1,NJTR
              dH_dX(nj)=dH_dP*ZG(2,NU1(nj)) !is cpt of grad(h)
            ENDDO
          ENDIF

          !Calculate mitochondria binding function g
          RATIO=PO2_M/P50_MC
          RATIO_N=RATIO**RN_MC
          IF(DABS(PO2_M).gt.1.d-5) THEN !pO2 is gt zero
            G=Rate_MC_Max*RATIO_N/(1.d0+RATIO_N)
          ELSE                         !pO2 is zero
            G=0.0d0
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' G='',E12.3)') G
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

          !Calculate element residuals for muscle equation
          nb=NBH(NH_LOC(1,nx))
          DO ns=1,NST(nb)
            PGN=PG(ns,1,ng,nb)
            DO nj=1,NJTR !are cpts of grad(psi)
              dPSI_dX(nj)=PGX(nb,nj,ns,DXIX,PG(1,1,ng,nb))
            ENDDO

            IF(ITYP10(1).EQ.1) THEN !rectangular cartesian coordinates
              SUM1=0.0d0
              SUM2=0.0d0
              DO nj=1,NJTR
                SUM1=SUM1+ Solu_M*Diff_M*dM_dX(nj)*dPSI_dX(nj) !is oxy diffusion
                SUM2=SUM2+ Diff_MB*dF_dX(nj)*dPSI_dX(nj)       !is myo diffusion
              ENDDO
              SUM3=G*PGN                               !is metabolic sink
              SUM4=Solu_M*Exch_coeff*(PO2_M-PO2_B)*PGN !is exchange term

            ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal coordinates
              SUM1=Solu_M*Diff_M*(G11*(dM_dX(1)*dPSI_dX(1)
     '                           +dM_dX(2)*dPSI_dX(2))
     '                           +G33*dM_dX(3)*dPSI_dX(3))
              SUM2=Diff_MB*(G11*(dF_dX(1)*dPSI_dX(1)
     '                          +dF_dX(2)*dPSI_dX(2))
     '                     +G33* dF_dX(3)*dPSI_dX(3))
              SUM3=G*PGN                               !is metabolic sink
              SUM4=Solu_M*Exch_coeff*(PO2_M-PO2_B)*PGN !is exchange term

            ELSE
              ERROR='>>Other coordinate systems not implemented yet'
              GO TO 9999
            ENDIF
            RE(ns,1)=RE(ns,1)+(1.d0-Rel_vol)*(SUM1+SUM2+SUM3+SUM4)*RWG
          ENDDO

          IF(KTYP15.GT.1) THEN !multi-field
            !Calculate element residuals for blood equation
            !Note assumption that blood velocity is isotropic
            nb=NBH(NH_LOC(2,nx))
            DO ns=1,NST(nb)
              PGN=PG(ns,1,ng,nb)
              DO nj=1,NJTR
                dPSI_dX(nj)=PGX(nb,nj,ns,DXIX,PG(1,1,ng,nb)) !is cpt of grad(psi)
              ENDDO

              IF(ITYP10(1).EQ.1) THEN !rectangular cartesian coordinates
                SUM1=0.0d0
                SUM2=0.0d0
                SUM3=0.0d0
                DO nj=1,NJTR
                  SUM1=SUM1+ Solu_B*Diff_B *dB_dX(nj)*dPSI_dX(nj) !is oxy diff.n
                  SUM2=SUM2+ Solu_B*Veloc*dB_dX(nj)*PGN         !is oxy advec.n
                  SUM3=SUM3+ Veloc*dH_dX(nj)*PGN                !is haemo advec.
                ENDDO
                SUM4=Solu_B*Exch_coeff*(PO2_B-PO2_M)*PGN        !is exchange
              ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal coordinates
                SUM1=Solu_B*Diff_B*(G11*(dB_dX(1)*dPSI_dX(1)
     '                             +dB_dX(2)*dPSI_dX(2))
     '                             +G33*dB_dX(3)*dPSI_dX(3))
                SUM2=Solu_B*Veloc*(G11*(dB_dX(1)+dB_dX(2))
     '            +G33*dB_dX(3))*PGN
                SUM3=Veloc*(G11*(dH_dX(1)+dH_dX(2))+G33*dH_dX(3))*PGN
                SUM4=Solu_B*Exch_coeff*(PO2_B-PO2_M)*PGN     !is exchange term
              ELSE
                ERROR='>>Other coordinate systems not implemented yet'
                GO TO 9999
              ENDIF
              RE(ns,2)=RE(ns,2)+Rel_vol*(SUM1+SUM2+SUM3+SUM4)*RWG
            ENDDO
          ENDIF

        ELSE IF(ITYP2(nr,nx).EQ.9.AND.ITYP3(nr,nx).EQ.2) THEN
C ***     Glucose-oxygen transport equation (two-field)
          PO2=ZG(1,1) !is oxygen  partial pressure
          PGL=ZG(2,1) !is glucose partial pressure
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' PO2='',E12.3,'' PGL='',E12.3)')
     '        PO2,PGL
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

C         Calculate gradient vector of oxygen partial pressure
          DO nj=1,NJTR
            dO2_dX(nj)=ZG(1,NU1(nj))
          ENDDO

C         Calculate gradient vector of glucose partial pressure
          DO nj=1,NJTR
            dGL_dX(nj)=ZG(2,NU1(nj))
          ENDDO

C         Calculate element residuals for oxygen equation
          nb=NBH(NH_LOC(1,nx))
          DO ns=1,NST(nb)
            PGN=PG(ns,1,ng,nb)
            DO nj=1,NJTR !are cpts of grad(psi)
              dPSI_dX(nj)=PGX(nb,nj,ns,DXIX,PG(1,1,ng,nb))
            ENDDO

            IF(ITYP10(1).EQ.1) THEN !rectangular cartesian coordinates
              SUM1=0.0d0
              SUM2=0.0d0
              DO nj=1,NJTR
                SUM1=SUM1+CG(3,ng)*dO2_dX(nj)*dPSI_dX(nj)!oxy diff.n
              ENDDO
              SUM4=0.d0                                  !reaction term
            ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal coordinates
              SUM1=CG(3,ng)*(G11*(dO2_dX(1)*dPSI_dX(1)
     '                           +dO2_dX(2)*dPSI_dX(2))
     '                      +G33* dO2_dX(3)*dPSI_dX(3))
              SUM4=0.d0                                  !reaction term
            ELSE
              ERROR='>>Other coordinate systems not implemented yet'
              GO TO 9999
            ENDIF
            RE(ns,1)=RE(ns,1)+(SUM1+SUM2+SUM3+SUM4)*RWG
          ENDDO

C         Calculate element residuals for glucose equation
          nb=NBH(NH_LOC(2,nx))
          DO ns=1,NST(nb)
            PGN=PG(ns,1,ng,nb)
            DO nj=1,NJTR !are cpts of grad(psi)
              dPSI_dX(nj)=PGX(nb,nj,ns,DXIX,PG(1,1,ng,nb))
            ENDDO

            IF(ITYP10(1).EQ.1) THEN !rectangular cartesian coordinates
              SUM1=0.0d0
              SUM2=0.0d0
              SUM3=0.0d0
              DO nj=1,NJTR
                SUM1=SUM1+CG(6,ng)*dGL_dX(nj)*dPSI_dX(nj)!glucose diff.
              ENDDO
              SUM4=0.d0                                  !reaction term
            ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal coordinates
              SUM2=CG(6,ng)*(G11*(dGL_dX(1)*dPSI_dX(1)
     '                           +dGL_dX(2)*dPSI_dX(2))
     '                      +G33* dGL_dX(3)*dPSI_dX(3))
              SUM4=0.d0                                   !is reaction term
            ELSE
              ERROR='>>Other coordinate systems not implemented yet'
              GO TO 9999
            ENDIF
            RE(ns,2)=RE(ns,2)+(SUM1+SUM2+SUM3+SUM4)*RWG
          ENDDO

        ELSE
          ERROR='>>This nonlinear equation not implemented yet'
          GO TO 9999
        ENDIF
      ENDDO !ng

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,NHE
          nh=nh_loc(nhx,nx)
          WRITE(OP_STRING,
     '      '(/'' RE(ns,'',I1,''): '',10E12.3,(/10E12.3))')
     '      nh,(RE(ns,nh),ns=1,NST(NBH(nh)))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZERE30')
      RETURN
 9999 CALL ERRORS('ZERE30',ERROR)
      CALL EXITS('ZERE30')
      RETURN 1
      END



