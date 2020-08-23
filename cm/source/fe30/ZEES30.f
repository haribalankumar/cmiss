      SUBROUTINE ZEES30(NBH,NBJ,NHE,nr,nx,
     '  CG,DDIFF,ES,PG,WG,XE,XG,YG,ZE,ZG,ERROR,*)

C#### Subroutine: ZEES30
C###  Description:
C###    ZEES30 calculates derivatives of element residual from
C###    current element dependent variable array ZE to form
C###    the tangent stiffness matrix ES for non-linear problems .

C**** (ityp5=5 ityp2=3)
C****   Elliptic eikonal equation for wavefront path determination.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),NHE,nr,nx
      REAL*8 CG(NMM,NGM),DDIFF(NSM),ES(NHM*NSM,NHM*NSM),
     '  PG(NSM,NUM,NGM,NBM),WG(NGM,NBM),XE(NSM,NJM),
     '  XG(NJM,NUM),YG(NIYGM,NGM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mb,mh,mi,mhs,mhx,ms,MSTB,nb,ng,nh,nhs,nhx,
     '  ni,NITB,niyg,nj,njj,ns,NSTNAT,NU1(0:3),NU2(3,3)
      REAL*8 ADVM(3),CONT0(2),DADV,DDEN,DDENM(3),DGDIFF(3),DIFFM1(3),
     '  DIFFM2(3,3),DIFFMODM(3),DNORM,DPGXI(3),DRES,DXIX(3,3),
     '  DZXI(3),GL(3,3),GU(3,3),GUxDPG(3),GUxDZ(3),GU2DZ(3),MDDEN1,
     '  MDDEN2,MDPET,MGAL,MPET,PETM(3),PGNSI,RAD,RES,RG,
     '  RWG,SUM,TEMP,VELOC,WGHTM0,WGHTM1(3)
      LOGICAL MATERIALPET,MODDIFF,PETROV
!     External Functions
      REAL*8 DDOT

      DATA NU1/1,2,4,7/
      DATA NU2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('ZEES30',*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' NHE='',I2,'' NBH(1)='',I2)') NHE,NBH(1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      nhs=0
      DO nhx=1,NHE !Melge uses NHE instead of NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        nb=NBH(nh)
        DO ns=1,NST(nb)+NAT(nb)
          nhs=nhs+1
          mhs=0
          DO mhx=1,NHE !Melge uses NHE instead of NH_LOC(0,nx)
            mh=NH_LOC(mhx,nx)
            mb=NBH(mh)
            DO ms=1,NST(mb)+NAT(mb)
              mhs=mhs+1
              ES(mhs,nhs)=0.0d0
            ENDDO !ms
          ENDDO !mhx
        ENDDO !ns
      ENDDO !nhx

      NITB=NIT(NBJ(1))

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

        IF(ITYP5(nr,nx).EQ.5) THEN !wavefront problem
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
          RWG=RG*WG(ng,NBH(NH_LOC(1,nx)))
        ENDIF
        CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
C       Jacobian * Gauss weight
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' ng='',I2,'' RG='',E11.4,'' WG='',E11.4)')
     '      ng,RG,WG(ng,1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' RWG='',E11.4,'' ZG:'',4E11.4)')
     '      RWG,(ZG(1,NU1(ni)),ni=0,NITB)
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
            nb=NBH(nh)
            MSTB=NST(nb)
            NSTNAT=MSTB+NAT(nb)

C           Copy first derivatives of ZG wrt xi into grad(time).
            DO ni=1,NITB
              DZXI(ni)=ZG(1,NU1(ni))
            ENDDO !ni
C           Find temporary vector of coupling and grad(time), GUxDZ(mi).
            DO ni=1,NITB
              SUM=0.0d0
              DO mi=1,NITB
                SUM=SUM+DZXI(mi)*GU(mi,ni)
              ENDDO !mi
              GUxDZ(ni)=SUM
            ENDDO !ni
C           Coupling norm of first derivatives.
            DNORM=0.0d0
            DO ni=1,NITB
              DNORM=DNORM+GUxDZ(ni)*DZXI(ni)
            ENDDO
C           For advection term
            VELOC=DSQRT(DNORM)
C            VELOC=DSQRT(DNORM+(CG(1,ng)/CG(2,ng))**2)
C           Multiply GU by weights RWG for diffusion.
C           Multiply GUxDZ by factors for advection.
            IF(VELOC.NE.0.0d0) THEN
              TEMP=FACTOR*CG(2,ng)/VELOC*RWG
            ELSE
              TEMP=0.0d0
            ENDIF
            DO ni=1,NITB
              DO mi=1,NITB
                DIFFM2(mi,ni)=GU(mi,ni)*RWG
              ENDDO !mi
              ADVM(ni)=GUxDZ(ni)*TEMP
            ENDDO !ni
C           Store diffusion modifier with advection as it multiplies
C           weighting function not its derivative.
            IF(MODDIFF) THEN
              DO ni=1,NITB
                ADVM(ni)=ADVM(ni)-DIFFMODM(ni)
              ENDDO
            ENDIF

            IF(PETROV) THEN !Petrov-Galerkin

C             Diffusion term for Petrov weights,
C             including Jacobian and Gauss weights
              DO ns=1,NSTNAT
                DDIFF(ns)=0.0d0
              ENDDO !ns
              DO ni=1,NITB
                CALL DAXPY(NSTNAT,DIFFM1(ni),
     '            PG(1,NU1(ni),ng,nb),1,DDIFF,1)
                DO mi=1,NITB
                  CALL DAXPY(NSTNAT,DIFFM2(mi,ni),
     '              PG(1,NU2(mi,ni),ng,nb),1,DDIFF,1)
                ENDDO
              ENDDO !ni

C***          Residual
              RES=-FACTOR*CG(1,ng)*RWG
     '          -DDOT(NSTNAT,DDIFF,1,ZE(1,nh),1)
              DO ni=1,NITB
                RES=RES+ADVM(ni)*DZXI(ni)
              ENDDO
C              RES=FACTOR*(CG(2,ng)*VELOC-1.4142136d0*CG(1,ng))*RWG
C     '          -DDOT(NSTNAT,DDIFF,1,ZE(1,nh),1)
C              IF(MODDIFF) THEN
C                DO ni=1,NITB
C                  RES=RES-DIFFMODM(ni)*DZXI(ni)
C                ENDDO
C              ENDIF

C             Calculate multipliers for Galerkin and derivative weights.
              CALL PETROV_PREP(NITB,CG(1,ng),CONT0,DNORM,DZXI,GU,GUxDZ,
     '          GU2DZ,MDDEN1,MDDEN2,MGAL,MPET,
     '          MATERIALPET,.TRUE.,ERROR,*9999)

C***          Apply weights and assemble into element stiffness matrix.

              MDPET=RES*MPET
              IF(MATERIALPET) THEN !deriv based on material direction
                DO ni=1,NITB
                  PETM(ni)=MPET*GUxDZ(ni)
                  DDENM(ni)=MDDEN1*GUxDZ(ni)+MDDEN2*GU2DZ(ni)
                ENDDO
              ELSE !deriv based on element xi direction
                DO ni=1,NITB
                  PETM(ni)=MPET*DZXI(ni)
                  DDENM(ni)=MDDEN1*DZXI(ni)+MDDEN2*GUxDZ(ni)
                ENDDO
              ENDIF !MATERIALPET
C             Loop over basis functions
              DO ns=1,NSTNAT
C***            Derivatives of products of residual and weights.
                DO mi=1,NITB
                  DGDIFF(mi)=0.0d0
                ENDDO !mi
                DADV=0.0d0
                DDEN=0.0d0
                DO ni=1,NITB
                  DPGXI(ni)=PG(ns,NU1(ni),ng,nb)
C                 Set up temporary vector DGDIFF(ni) to include
C                 everything for diffusion except weighting function
C                 derivative wrt xi.
                  DO mi=1,NITB
                    DGDIFF(mi)=DGDIFF(mi)+DIFFM2(mi,ni)*DPGXI(ni)
                  ENDDO !mi
C                 Set up DADV to include everything for advection except
C                 weighting function.
                  DADV=DADV+ADVM(ni)*DPGXI(ni)
C                 Petrov weight term
                  DDEN=DDEN+DDENM(ni)*DPGXI(ni)
                ENDDO !ni
                DRES=DADV-DDIFF(ns)
C               Galerkin weight and part of residual.
                WGHTM0=MGAL*DADV
C               Multipliers of weights (0th deriv).
                CALL DAXPY(MSTB,WGHTM0,
     '            PG(1,1,ng,nb),1,ES(1,ns),1)
C               Petrov derivative weights and deriv part of Galerkin.
                DO ni=1,NITB
                  GUxDPG(ni)=0.0d0
                  DO mi=1,NITB
                    GUxDPG(ni)=GUxDPG(ni)+DPGXI(mi)*GU(mi,ni)
                  ENDDO !mi
                  WGHTM1(ni)=MGAL*RWG*GUxDPG(ni)+ !Galerkin diffusion
     '              PETM(ni)*DRES      !Petrov weight * res deriv
                ENDDO !ni
                IF(MATERIALPET) THEN
                  DO ni=1,NITB
                    WGHTM1(ni)=WGHTM1(ni)+
     '                MDPET*(GUxDPG(ni)+DDEN*GUxDZ(ni)) !res * w deriv
                  ENDDO !ni
                ELSE
                  DO ni=1,NITB
                    WGHTM1(ni)=WGHTM1(ni)+
     '                MDPET*(DPGXI(ni)+DDEN*DZXI(ni)) !res * w deriv
                  ENDDO !ni
                ENDIF !MATERIALPET
                DO ni=1,NITB
                  CALL DAXPY(MSTB,WGHTM1(ni),
     '              PG(1,NU1(ni),ng,nb),1,ES(1,ns),1)
                ENDDO !ni
              ENDDO !ns

            ELSE !Galerkin weights
              DO ns=1,NSTNAT
C               Set up temporary vector DGDIFF(ni) to include everything for
C               diffusion except weighting function derivative wrt xi.
                DO mi=1,NITB
                  DGDIFF(mi)=0.0d0
                ENDDO !mi
                DO ni=1,NITB
                  PGNSI=PG(ns,NU1(ni),ng,nb)
                  DO mi=1,NITB
                    DGDIFF(mi)=DGDIFF(mi)+DIFFM2(mi,ni)*PGNSI
                  ENDDO !mi
                ENDDO !ni
C               Add diffusion to stiffness matrix.
                DO mi=1,NITB
                  CALL DAXPY(MSTB,DGDIFF(mi),PG(1,NU1(mi),ng,nb),1,
     '              ES(1,ns),1)
                ENDDO !mi
C               Set up DADV to include everything for advection except
C               weighting function.
                DADV=0.0d0
                DO ni=1,NITB
                  DADV=DADV+ADVM(ni)*PG(ns,NU1(ni),ng,nb)
                ENDDO !ni
C               Add advection to stiffness matrix.
                CALL DAXPY(MSTB,DADV,PG(1,1,ng,nb),1,ES(1,ns),1)
              ENDDO !ns

            ENDIF !Petrov / Galerkin

          ELSE
            ERROR='>>Only elliptic eikonal equation is implemented'
            GO TO 9999
          ENDIF

        ELSE
          ERROR=' Analytic derivs not implemented'
          GO TO 9999
        ENDIF
      ENDDO

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        mhs=0
        DO mhx=1,NHE !Melge uses NHE instead of NH_LOC(0,nx)
          mh=NH_LOC(mhx,nx)
          mb=NBH(mh)
          DO ms=1,NST(mb)+NAT(mb)
            mhs=mhs+1
            nhs=0
            DO nhx=1,NHE !Melge uses NHE instead of NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
              nb=NBH(nh)
              NSTNAT=NST(nb)+NAT(nb)
              WRITE(OP_STRING,
     '          '(/'' ES('',I1,'',nhs): '',10E12.3,(/10E12.3))')
     '          nhs,(ES(mhs,nhs+ns),ns=1,NSTNAT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              nhs=nhs+NSTNAT
            ENDDO !mhx
          ENDDO !ns
        ENDDO !nhx
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZEES30')
      RETURN
 9999 CALL ERRORS('ZEES30',ERROR)
      CALL EXITS('ZEES30')
      RETURN 1
      END



