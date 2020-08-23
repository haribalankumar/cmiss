      SUBROUTINE ZFFS30(IDOXFT,JDOXFT,NBHF,NBJF,NEAXT,NHF,NIEF,
     '  NORMALSIGN1,nr,nx,CG,FS,PG,XDF,YGF,ZDF,ERROR,*)

C#### Subroutine: ZFFS30
C###  Description:
C###    ZFFS30 calculates face stiffness matrices FS from integrals of
C###    surfaces (dimension NJ_LOC(NJL_GEOM,0,nr)-1).
C***  Assumes nc=1.
C***  (Upwind terms will probably only be used for this matrix.)

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'petg00.inc'
!     Parameter List
      INTEGER IDOXFT(NHM),JDOXFT,NBHF(NHM),NBJF(NJM),NEAXT,NHF,
     '  NIEF(0:2),NORMALSIGN1,nr,nx
      REAL*8 CG(NMM,NGM),FS(NSFM,2,NSFM,2,2,NHM),PG(NSM,NUM,NGM,NBM),
     '  XDF(NSFM,2,2,NJM),
     '  YGF(NIYGFM,NGFM),ZDF(NSFM,2,2,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER idoxf,jdoxf,mi,mi_f,ms,MSTB,nb,nbh_f,neax,NEAFST,niyg,ng,
     '  nhx,nh,ni,ni_f,NITE,NITF,nj,njj,NJTR,ns,NSTNAT,NU1(3)
      REAL*8 CONT0(2),DDEN,DDENM(0:2),DENOM1,DENOM2,
     '  DIFF,DINFL(0:2,0:2),DNORM,DPET(0:2),DPGN(2),DRES,DZN,DZNM(3,2),
     '  DZXI(3),GU(3,3),GUxDPG,GUxDZ(3),GU2DZ(3),MDDEN1,MDDEN2,MDIN1,
     '  MDIN2,MDIN3,MDPET,MGAL,MNCPT,MPET,MRESI,PETM(0:2),PPET,PGN,RES,
     '  RESI,RESM(0:2,2),RWG(2),SUM,TEMP,U(3),V(3)
      LOGICAL EXTERN,FLUX,INFLOW,MATERIALPET,PETROV
!     External Functions
      REAL*8 DDOT

      DATA NU1/2,4,7/

      CALL ENTERS('ZFFS30',*9999)

      nbh_f=NBHF(NH_LOC(1,nx))
      NITF=NIT(nbh_f) !number of xi dirns in the face
      NITE=NITF+1
      NJTR=NJ_LOC(NJL_GEOM,0,nr)
      CALL ASSERT(NJTR.EQ.NITE,'>>NITF!=NJTR-1',ERROR,*9999)
      MATERIALPET=ITYP15(nr,nx).EQ.1
      PETROV=MATERIALPET.OR.ITYP15(nr,nx).EQ.2
      EXTERN=NEAXT.EQ.1 !domain boundary face
C!!!  need to check order of stabilizing derivs
      FLUX=PETROV.OR..TRUE.

      DO nhx=1,NHF
        nh=NH_LOC(nhx,nx)
        nb=NBHF(nh)
        MSTB=NST(nb)
        NSTNAT=MSTB+NAT(nb)
        IF(FLUX) THEN
          NEAFST=NEAXT
        ELSE
          NEAFST=1
        ENDIF
        DO neax=1,NEAFST
          DO jdoxf=1,JDOXFT
            DO ns=1,NSTNAT
              DO idoxf=1,IDOXFT(nhx)
                DO ms=1,MSTB
                  FS(ms,idoxf,ns,jdoxf,neax,nhx)=0.0d0
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL ASSERT(JTYP10.LT.2,
     '  '>>Only conventional interpolation is implemented',ERROR,*9999)

      DO ng=1,NGT(nbh_f) !assume same gauss points for all dep. vars.

        IF(PETROV) THEN !Petrov-Galerkin
          njj=1
C         This field should contain a C0 continuous variable that is
C         1 over most of the domain but approaches zero at Dirichlet
C         boundary conditions.
          CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).GE.njj,
     '      '>>Boundary condition field not defined',ERROR,*9999)
          nj=NJ_LOC(NJL_FIEL,njj,nr)
          nb=NBJF(nj)
          IF(nb.NE.0) THEN
C           Variable
            CONT0(1)=
     '        DDOT(NST(nb)+NAT(nb),PG(1,1,ng,nb),1,XDF(1,1,1,nj),1)
C           Multiplier for unit derivatives
            IF(NKT(0,nb).EQ.1) THEN !linear element
              CONT0(2)=1.0d0
            ELSE !cubic element
              CONT0(2)=1.0d0/3.0d0
            ENDIF
          ELSE
            CONT0(1)=1.0d0
            CONT0(2)=1.0d0
          ENDIF
          IF(.NOT.EXTERN) THEN
            CALL ASSERT(CONT0(1).EQ.0.0d0,
     '        '>>Continuity field not zero at discontinuity',
     '        ERROR,*9999)
          ENDIF
        ENDIF

        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' ng='',I2)') ng
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

C***    Get information contained in YG
        niyg=0
        IF(PETROV) THEN
          IF(NEAXT.EQ.1) THEN !boundary face
            niyg=niyg+1
            RWG(1)=YGF(niyg,ng)
C***        Set up coupling tensor in xi coords GU from YG.
            DO mi=1,NITE
              DO ni=1,mi
                niyg=niyg+1
                GU(mi,ni)=YGF(niyg,ng)
                IF(mi.NE.ni) THEN
                  GU(ni,mi)=GU(mi,ni)
                ENDIF
              ENDDO !ni
            ENDDO !mi
            niyg=niyg+1
            MNCPT=YGF(niyg,ng) !for inflow b.c.
C***        Evaluate residual multiplier (including Gauss weights)
            DO ni_f=0,NITF
              RESM(ni_f,1)=NORMALSIGN1*RWG(1)*GU(NIEF(ni_f),NIEF(0))
            ENDDO
          ELSE !internal face
            DO neax=1,NEAXT
              niyg=niyg+1
              RWG(neax)=YGF(niyg,ng)
C***          Evaluate residual multipliers (including Gauss weights)
              DO ni=1,NITE
                niyg=niyg+1
                DZNM(ni,neax)=YGF(niyg,ng)
              ENDDO !ni_f
              DO ni_f=0,NITF
                RESM(ni_f,neax)=
     '            NORMALSIGN1*RWG(neax)*DZNM(NIEF(ni_f),neax)
              ENDDO !ni_f
            ENDDO !neax
          ENDIF !NEAXT
        ELSE IF(FLUX) THEN !and not PETROV
          DO neax=1,NEAXT
C***        Evaluate residual multipliers (including Gauss weights)
            DO ni=1,NITE
              niyg=niyg+1
              DZNM(ni,neax)=YGF(niyg,ng)
            ENDDO !ni
            DO ni_f=0,NITF
              RESM(ni_f,neax)=FACTOR*NORMALSIGN1*DZNM(NIEF(ni_f),neax)
            ENDDO !ni_f
          ENDDO !neax
        ELSE !higher derivative discontinuity
          niyg=niyg+1
          RWG(1)=YGF(niyg,ng)
        ENDIF !PETROV

C*** Elliptic eikonal equation for wavefront path determination.
C***   CG(1,ng) is the time constant source term,
C***   CG(2,ng) is the dimensionless coefficient of advection, co
C***   CG(3..,ng) are coupling coefficients.
        IF(FLUX) THEN
          nhx=1
          IF(EXTERN) THEN
C           First derivative of dependent var wrt out of face xi.
            DZXI(NIEF(0))=
     '        DDOT(NSTNAT,PG(1,1,ng,nbh_f),1,ZDF(1,2,1,nhx),1)
C           First derivatives wrt in face xi.
            DO ni_f=1,NITF
              DZXI(NIEF(ni_f))=
     '          DDOT(NSTNAT,PG(1,NU1(ni_f),ng,nbh_f),1,ZDF(1,1,1,nhx),1)
            ENDDO !ni
          ENDIF !EXTERN
        ENDIF !FLUX

        IF(PETROV) THEN !Flux Term.

          IF(EXTERN) THEN !derivative weights
C           Find temporary vector of coupling and grad(time), GUxDZ(ni).
            DO ni=1,NITE
              SUM=0.0d0
              DO mi=1,NITE
                SUM=SUM+DZXI(mi)*GU(mi,ni)
              ENDDO
              GUxDZ(ni)=SUM
            ENDDO
C           Coupling norm of first derivatives.
            DNORM=0.0d0
            DO ni=1,NITE
              DNORM=DNORM+GUxDZ(ni)*DZXI(ni)
            ENDDO
            DZN=GUxDZ(NIEF(0))
            RES=NORMALSIGN1*DZN*RWG(1)
C           Calculate multipliers for Galerkin and derivative weights.
            CALL PETROV_PREP(NITE,CG(1,ng),CONT0,DNORM,DZXI,GU,GUxDZ,
     '        GU2DZ,MDDEN1,MDDEN2,MGAL,MPET,
     '        MATERIALPET,.TRUE.,ERROR,*9999)
          ELSE !only Galerkin weights
            MPET=0.0d0
            MGAL=1.0d0
          ENDIF

          IF(MPET.NE.0.0d0) THEN
            MDPET=RES*MPET
            IF(MATERIALPET) THEN !deriv based on material direction
              DO ni_f=0,NITF
                ni=NIEF(ni_f)
                PETM(ni_f)=MPET*GUxDZ(ni)
                DDENM(ni_f)=MDDEN1*GUxDZ(ni)+MDDEN2*GU2DZ(ni)
              ENDDO
            ELSE !deriv based on element xi direction
              DO ni_f=0,NITF
                ni=NIEF(ni_f)
                PETM(ni_f)=MPET*DZXI(ni)
                DDENM(ni_f)=MDDEN1*DZXI(ni)+MDDEN2*GUxDZ(ni)
              ENDDO
            ENDIF !MATERIALPET
            INFLOW=EXTERN.AND.(NORMALSIGN1.GT.0.NEQV.DZN.GT.0.0d0)
            IF(INFLOW) THEN
C             Apply inflow boundary condition
              TEMP=CG(2,ng)*FACTOR*DSQRT(MNCPT) !Peclet number
              PPET=PETRATIO*TEMP
              PPET=TEMP*(1.0d0+PPET/(PETLIMIT*PPET+1.0d0))
C              MRESI=CONT0(1)*PPET*RWG(1)/2
              MRESI=PPET*RWG(1)/2
              TEMP=DNORM+PETSMOOTH*CG(1,ng)*CG(1,ng)/(CG(2,ng)*CG(2,ng))
              DIFF=TEMP-MNCPT*DZN*DZN
              DENOM1=DSQRT(TEMP)
              DENOM2=DSQRT(DIFF)
              RESI=MRESI*(DENOM1-DENOM2)
              TEMP=MNCPT*DZN
              DO ni=1,NITE
                U(ni)=(GUxDZ(ni)-TEMP*GU(ni,NIEF(0)))/DENOM2
                V(ni)=GUxDZ(ni)/DENOM1-U(ni)
              ENDDO !ni
              TEMP=RESI/DENOM1
              MDIN1=RESI/DENOM2
              MDIN2=TEMP-MDIN1
              MDIN3=TEMP/(DENOM1*DENOM1)
              DO ni_f=0,NITF
                ni=NIEF(ni_f)
                DO mi_f=0,NITF
                  mi=NIEF(mi_f)
                  DINFL(mi_f,ni_f)=
     '              MDIN1*(U(mi)*U(ni)
     '                     +MNCPT*GU(mi,NIEF(0))*GU(ni,NIEF(0)))
     '              +MDIN2*GU(mi,ni)
     '              -MDIN3*GUxDZ(mi)*GUxDZ(ni)
     '              +MRESI*V(mi)*V(ni)
                ENDDO !mi
              ENDDO !ni
C              MDIN1=MRESI*MNCPT
C              DO ni_f=0,NITF
C                ni=NIEF(ni_f)
C                DO mi_f=0,NITF
C                  mi=NIEF(mi_f)
C                  DINFL(mi_f,ni_f)=
C     '              MDIN1*GU(mi,NIEF(0))*GU(ni,NIEF(0))
C                ENDDO !mi
C              ENDDO !ni
            ENDIF !inflow
          ENDIF !MPET!=0

C***      Apply weights and assemble into face stiffness matrix.

C         Loop over dep var dofs
          DO neax=1,NEAXT
            DO jdoxf=1,JDOXFT
              DO ns=1,NSTNAT
C***            Derivatives of products of residual and weights.
C               Assemble weight * res deriv + res * weight deriv into
C               stiffness matrix.
                IF(jdoxf.EQ.1) THEN
C                 Residual deriv wrt in face derivs of dep var
                  DRES=0.0d0
                  DO ni_f=1,NITF
                    DPGN(ni_f)=PG(ns,NU1(ni_f),ng,nbh_f)
                    DRES=DRES+RESM(ni_f,neax)*DPGN(ni_f)
                  ENDDO
                ELSE
C                 Residual deriv wrt out of face derivs of dep var
                  PGN=PG(ns,1,ng,nbh_f)
                  DRES=RESM(0,neax)*PGN
                ENDIF
!C               Galerkin weights.
!                DGAL=MGAL*DRES
!                CALL DAXPY(MSTB,DGAL,
!     '            PG(1,1,ng,nbh_f),1,FS(1,1,ns,jdoxf,neax,nhx),1)
C               Petrov derivative weights:
                IF(MPET.NE.0.0d0) THEN
C                 Weight deriv * residual.
                  IF(jdoxf.EQ.1) THEN !weight deriv wrt in face derivs
                    DDEN=0.0d0
                    DO ni_f=1,NITF
                      DDEN=DDEN+DDENM(ni_f)*DPGN(ni_f)
                    ENDDO
                  ELSE !weight deriv wrt out of face derivs
                    DDEN=DDENM(0)*PGN
                  ENDIF !jdoxf
                  IF(MATERIALPET) THEN !material deriv
                    DO ni_f=0,NITF
                      ni=NIEF(ni_f)
                      IF(jdoxf.EQ.1) THEN !weight deriv wrt in face derivs
                        SUM=0.0d0
                        DO mi_f=1,NITF
                          SUM=SUM+DPGN(mi_f)*GU(ni,NIEF(mi_f))
                        ENDDO !mi
                        GUxDPG=SUM
                      ELSE !weight deriv wrt out of face derivs
                        GUxDPG=PGN*GU(ni,NIEF(0))
                      ENDIF !jdoxf
                      DPET(ni_f)=MDPET*(GUxDPG+DDEN*GUxDZ(NIEF(ni_f)))
                    ENDDO
                  ELSE !xi deriv
                    IF(jdoxf.EQ.1) THEN !weight deriv wrt in face derivs
                      DPET(0)=MDPET*DDEN*DZXI(NIEF(0))
                      DO ni_f=1,NITF
                        DPET(ni_f)=
     '                    MDPET*(DPGN(ni_f)+DDEN*DZXI(NIEF(ni_f)))
                      ENDDO
                    ELSE !weight deriv wrt out of face derivs
                      DPET(0)=MDPET*(PGN+DDEN*DZXI(NIEF(0)))
                      DO ni_f=1,NITF
                        DPET(ni_f)=MDPET*DDEN*DZXI(NIEF(ni_f))
                      ENDDO
                    ENDIF !jdoxf
                  ENDIF !material/xi deriv
                  IF(INFLOW) THEN
                    DO ni_f=0,NITF
                      IF(jdoxf.EQ.1) THEN !weight deriv wrt in face derivs
                        DO mi_f=1,NITF
                          DPET(ni_f)=DPET(ni_f)+
     '                      DPGN(mi_f)*DINFL(mi_f,ni_f)
                        ENDDO !mi
                      ELSE !weight deriv wrt out of face derivs
                        DPET(ni_f)=DPET(ni_f)+PGN*DINFL(0,ni_f)
                      ENDIF !jdoxf
                    ENDDO !ni
                  ENDIF !inflow
C                 Out of face deriv weights.
                  CALL DAXPY(MSTB,PETM(0)*DRES+DPET(0),
     '              PG(1,1,ng,nbh_f),1,FS(1,2,ns,jdoxf,neax,nhx),1)
C                 In face deriv weights.
                  DO ni_f=1,NITF
                    CALL DAXPY(MSTB,PETM(ni_f)*DRES+DPET(ni_f),
     '                PG(1,NU1(ni_f),ng,nbh_f),1,
     '                FS(1,1,ns,jdoxf,neax,nhx),1)
                  ENDDO
                ENDIF !MPET.NE.0.0d0
              ENDDO !ns
            ENDDO !jdoxf
          ENDDO !neax

        ELSE IF(FLUX) THEN !first deriv discont
          IF(EXTERN) THEN !exterior face
C           Evaluate normal derivative
            RES=0.0d0
            DO ni_f=0,NITF
              RES=RES+RESM(ni_f,1)*DZXI(NIEF(ni_f))
            ENDDO
            INFLOW=RES.LT.0d0
          ELSE
            INFLOW=.TRUE. !apply discontinuity term
          ENDIF
C         Loop over dep var dofs
          IF(INFLOW) THEN
            DO neax=1,NEAXT
              DO jdoxf=1,JDOXFT
                DO ns=1,NSTNAT
C***              Derivatives of products of residual and weights.
C                 Assemble weight * res deriv into stiffness matrix.
                  IF(jdoxf.EQ.1) THEN
C                   Residual deriv wrt in face derivs of dep var
                    DRES=0.0d0
                    DO ni_f=1,NITF
                      DPGN(ni_f)=PG(ns,NU1(ni_f),ng,nbh_f)
                      DRES=DRES+RESM(ni_f,neax)*DPGN(ni_f)
                    ENDDO
                  ELSE
C                   Residual deriv wrt out of face derivs of dep var
                    PGN=PG(ns,1,ng,nbh_f)
                    DRES=RESM(0,neax)*PGN
                  ENDIF
                  CALL DAXPY(MSTB,DRES,
     '              PG(1,1,ng,nbh_f),1,FS(1,1,ns,jdoxf,neax,nhx),1)
                ENDDO !ns
              ENDDO !jdoxf
            ENDDO !neax
          ENDIF

CC           Loop over dep var dofs
C          DO neax=1,NEAXT
C            DO ns=1,NSTNAT
CC***          Derivatives of products of residual and weights.
CC             Assemble weight * res deriv into stiffness matrix.
CC             Residual deriv wrt out of face derivs of dep var
C              PGN=PG(ns,1,ng,nbh_f)
C              DRES=FACTOR*RESM(0,neax)*PGN
CC             Only derivative weights wrt out of face xi.
C              DO idoxf=1,IDOXFT(1)
C                CALL DAXPY(MSTB,RESM(0,idoxf)*DRES,
C     '            PG(1,1,ng,nbh_f),1,FS(1,idoxf,ns,2,neax,nhx),1)
C              ENDDO !idoxf
C            ENDDO !ns
C          ENDDO !neax

        ELSE !higher derivs
          DO nhx=1,NHF
            nh=NH_LOC(nhx,nx)
            nb=NBHF(nh)
            MSTB=NST(nb)
            NSTNAT=MSTB+NAT(nb)
C           Tensor product of basis functions multiplied by weights RWG.
            DO ns=1,NSTNAT
              CALL DAXPY(MSTB,PG(ns,1,ng,nb)*RWG(1),PG(1,1,ng,nb),1,
     '          FS(1,1,ns,1,1,nhx),1)
            ENDDO !ns
          ENDDO !nhx

        ENDIF !Flux / Derivative Discontinuity

      ENDDO !ng

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,NHF
          nh=NH_LOC(nhx,nx)
          nb=NBHF(nh)
          MSTB=NST(nb)
          NSTNAT=MSTB+NAT(nb)
          IF(FLUX) THEN
            NEAFST=NEAXT
          ELSE
            NEAFST=1
          ENDIF
          DO neax=1,NEAFST
            DO jdoxf=1,JDOXFT
              WRITE(OP_STRING,
     '          '('' FS(ms,'',I1,'',ns,'',I1,'','','//
     '          'I1,'','',I1,''):'')') idoxf,jdoxf,neax,nhx
              WRITE(OP_STRING,'('' '')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO ms=1,NSTNAT
                DO idoxf=1,IDOXFT(nhx)
                  WRITE(OP_STRING,
     '              '('' FS('',I2,'','',I1,'',ns,'',I1,'','','//
     '              'I1,'','',I1,''): '',8E12.4,(/12X,8E12.4))')
     '              ms,idoxf,jdoxf,neax,nhx,
     '              (FS(ms,idoxf,ns,jdoxf,neax,nhx),ns=1,NSTNAT)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZFFS30')
      RETURN
 9999 CALL ERRORS('ZFFS30',ERROR)
      CALL EXITS('ZFFS30')
      RETURN 1
      END



