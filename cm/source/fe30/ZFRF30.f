      SUBROUTINE ZFRF30(IDOXFT,NBHF,NBJF,NEAXT,NHF,NIEF,NORMALSIGN1,nr,
     '  nx,CG,PG,RDF,XDF,YGF,ZDF,ERROR,*)

C#### Subroutine: ZFRF30
C###  Description:
C###    ZFRF30 calculates residual components RDF from integrals
C###    of surfaces (dimension NJ_LOC(NJL_GEOM,0,nr)-1)
C###    using current element dependent variable array ZE.
C***  For flux integration, the normal used crosses the face in the
C***  direction that the cross face xi increases.
C***  Assumes nc=1.
C***  (Upwind terms will probably only be used for this matrix.)

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'petg00.inc'
!     Parameter List
      INTEGER IDOXFT(NHM),NBHF(NHM),NBJF(NJM),NEAXT,NHF,NIEF(0:2),
     '  NORMALSIGN1,nr,nx
      REAL*8 CG(NMM,NGM),PG(NSM,NUM,NGM,NBM),RDF(NSFM,2,NHM),
     '  XDF(NSFM,2,2,NJM),YGF(NIYGFM,NGFM),ZDF(NSFM,2,2,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER idoxf,mi,nb,nbh_f,neax,ng,nhx,nh,ni,ni_f,NITE,NITF,niyg,
     '  nj,njj,NJTR,ns,NSTB,NSTNAT,NU1(3)
      INTEGER*4 NULL
      REAL*8 CONT0(2),DENOM1,DENOM2,DIFF,DNORM,DZN(2),DZNM(3,2),
     '  DZXI(3,2),GU(3,3),GUxDZ(3),MGAL,MNCPT,MPET,MRESI,PETM(0:2),PPET,
     '  RESI,RWG(2),SUM,TEMP,TEMP1,TEMP2,ZGN
      LOGICAL EXTERN,FLUX,MATERIALPET,PETROV
!     External Functions
      REAL*8 DDOT

      PARAMETER (NULL=0)

      DATA NU1/2,4,7/

      CALL ENTERS('ZFRF30',*9999)

      nbh_f=NBHF(NH_LOC(1,nx))
      NSTNAT=NST(nbh_f)+NAT(nbh_f)
      NITF=NIT(nbh_f) !number of xi dirns in the face
      NITE=NITF+1
      NJTR=NJ_LOC(NJL_GEOM,0,nr)
      CALL ASSERT(NJTR.EQ.NITE,'>>NITF!=NJTR-1',ERROR,*9999)
      MATERIALPET=ITYP15(nr,nx).EQ.1
      PETROV=MATERIALPET.OR.ITYP15(nr,nx).EQ.2
      EXTERN=PETROV.AND.IDOXFT(1).EQ.2 !domain boundary flux
C!!!  need to check order of stabilizing derivs
      FLUX=PETROV.OR..TRUE.

      DO nhx=1,NHF
        nh=NH_LOC(nhx,nx)
        nb=NBHF(nh)
        NSTB=NST(nb)
        DO idoxf=1,IDOXFT(nhx)
          DO ns=1,NSTB
            RDF(ns,idoxf,nhx)=0.0d0
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
C***        Evaluate residual multiplier (not including Gauss weights)
            IF(.NOT.EXTERN) THEN !I don't think this'll happen
              DO ni=1,NITE
                DZNM(ni,1)=GU(ni,NIEF(0))
              ENDDO !ni
            ENDIF !EXTERN
          ELSE !internal face
            DO neax=1,NEAXT
              niyg=niyg+1
              RWG(neax)=YGF(niyg,ng)
C***          Evaluate residual multipliers (not including Gauss weights)
              DO ni=1,NITE
                niyg=niyg+1
                DZNM(ni,neax)=YGF(niyg,ng)
              ENDDO !ni
            ENDDO !neax
          ENDIF !NEAXT
        ELSE IF(FLUX) THEN !and not PETROV
          DO neax=1,NEAXT
C***        Evaluate residual multipliers (including Gauss weights)
            DO ni=1,NITE
              niyg=niyg+1
              DZNM(ni,neax)=YGF(niyg,ng)
            ENDDO !ni
          ENDDO !neax
        ELSE !higher derivative discontinuity
          niyg=niyg+1
          RWG(1)=YGF(niyg,ng)
        ENDIF !PETROV

C*** Elliptic eikonal equation for wavefront path determination.
C***   CG(1,ng) is the time constant source term,
C***   CG(2,ng) is the dimensionless coefficient of advection, co,
C***   CG(3..,ng) are coupling coefficients.
C***   FACTOR is the continuation parameter.

        IF(FLUX) THEN !Flux Term.
          nhx=1
C         First derivative of dependent var wrt out of face xi.
          DZXI(NIEF(0),1)=
     '      DDOT(NSTNAT,PG(1,1,ng,nbh_f),1,ZDF(1,2,1,nhx),1)
C         First derivatives wrt in face xi.
          DO ni_f=1,NITF
            DZXI(NIEF(ni_f),1)=
     '        DDOT(NSTNAT,PG(1,NU1(ni_f),ng,nbh_f),1,ZDF(1,1,1,nhx),1)
          ENDDO !ni
          IF(NEAXT.EQ.2) THEN
C           First derivative of dependent var wrt out of face xi.
            DZXI(NIEF(0),2)=
     '        DDOT(NSTNAT,PG(1,1,ng,nbh_f),1,ZDF(1,2,2,nhx),1)
C           First derivatives wrt in face xi.
            DO ni_f=1,NITF
              DZXI(NIEF(ni_f),2)=DZXI(NIEF(ni_f),1)
            ENDDO !ni
          ENDIF
        ENDIF !FLUX

        IF(PETROV) THEN !Petrov Flux Term.
C          IF(EXTERN.OR.(MATERIALPET.AND..NOT.PETLIMMAT)) THEN !variable wghts
          IF(EXTERN) THEN !derivative weights
C           Find temporary vector of coupling and grad(time), GUxDZ(ni).
            DO ni=1,NITE
              SUM=0.0d0
              DO mi=1,NITE
                SUM=SUM+DZXI(mi,1)*GU(mi,ni)
              ENDDO
              GUxDZ(ni)=SUM
            ENDDO
C           Coupling norm of first derivatives.
            DNORM=0.0d0
            DO ni=1,NITE
              DNORM=DNORM+GUxDZ(ni)*DZXI(ni,1)
            ENDDO
            DZN(1)=GUxDZ(NIEF(0))
C           Calculate multipliers for Galerkin and derivative weights.
            CALL PETROV_PREP(NITE,CG(1,ng),CONT0,DNORM,DZXI(1,1),
     '        GU,GUxDZ,%VAL(NULL),%VAL(NULL),%VAL(NULL),
     '        MGAL,MPET,MATERIALPET,.FALSE.,ERROR,*9999)
            ZGN=DZN(1)*RWG(1)
          ELSE !only Galerkin weights
            ZGN=0.0d0
            DO neax=1,NEAXT
              DZN(neax)=0.0d0
              DO ni=1,NITE
                DZN(neax)=DZN(neax)+DZNM(ni,neax)*DZXI(ni,neax)
              ENDDO
              ZGN=ZGN+DZN(neax)*RWG(neax)
            ENDDO
            MGAL=1.0d0
            MPET=0.0d0
          ENDIF

          IF(MPET.NE.0.0d0) THEN
            TEMP=MPET*NORMALSIGN1*ZGN
            IF(MATERIALPET) THEN !deriv based on material direction
              DO ni_f=0,NITF
                ni=NIEF(ni_f)
                PETM(ni_f)=TEMP*GUxDZ(ni)
              ENDDO
            ELSE !deriv based on xi direction
              DO ni_f=0,NITF
                ni=NIEF(ni_f)
                PETM(ni_f)=TEMP*DZXI(ni,1)
              ENDDO
            ENDIF !material deriv

            IF(EXTERN.AND.(NORMALSIGN1.GT.0.NEQV.DZN(1).GT.0.0d0)) THEN
C             Apply inflow boundary condition
              TEMP=CG(2,ng)*FACTOR*DSQRT(MNCPT) !Peclet number
              PPET=PETRATIO*TEMP
              PPET=TEMP*(1.0d0+PPET/(PETLIMIT*PPET+1.0d0))
C              MRESI=CONT0(1)*PPET*RWG(1)/2
              MRESI=PPET*RWG(1)/2
              TEMP=DNORM+PETSMOOTH*CG(1,ng)*CG(1,ng)/(CG(2,ng)*CG(2,ng))
              DIFF=TEMP-MNCPT*DZN(1)*DZN(1)
              DENOM1=DSQRT(TEMP)
              DENOM2=DSQRT(DIFF)
              RESI=MRESI*(DENOM1-DENOM2)
              TEMP=RESI/DENOM2
              TEMP1=RESI/DENOM1-TEMP
              TEMP2=MNCPT*DZN(1)*TEMP
              DO ni_f=0,NITF
                ni=NIEF(ni_f)
                PETM(ni_f)=
     '            PETM(ni_f)+TEMP1*GUxDZ(ni)+TEMP2*GU(ni,NIEF(0))
              ENDDO
C              TEMP2=MRESI*MNCPT*DZN(1)
C              DO ni_f=0,NITF
C                ni=NIEF(ni_f)
C                PETM(ni_f)=PETM(ni_f)+TEMP2*GU(ni,NIEF(0))
C              ENDDO
            ENDIF !inflow
          ENDIF !MPET!=0

C***      Make additions to the residual integrals.
!C         Galerkin weights.
!          CALL DAXPY(NSTB,MGAL*NORMALSIGN1*ZGN,
!     '      PG(1,1,ng,nbh_f),1,RDF(1,1,nhx),1)
C         Petrov derivatives.
          IF(MPET.NE.0.0d0) THEN
C           derivative weights wrt out of face xi.
            CALL DAXPY(NSTB,PETM(0),
     '        PG(1,1,ng,nbh_f),1,RDF(1,2,nhx),1)
C           derivative weights wrt in face xi.
            DO ni_f=1,NITF
              CALL DAXPY(NSTB,PETM(ni_f),
     '          PG(1,NU1(ni_f),ng,nbh_f),1,RDF(1,1,nhx),1)
            ENDDO !ELSE !can't choose a direction for derivative
          ENDIF

        ELSE IF(FLUX) THEN !first deriv discont
C         Evaluate integrand ZGN (with weights) from derivatives
          ZGN=0.0d0
          DO neax=1,NEAXT
            DO ni=1,NITE
              ZGN=ZGN+DZNM(ni,neax)*DZXI(ni,neax)
            ENDDO
          ENDDO
          ZGN=ZGN*FACTOR*NORMALSIGN1
          IF(NEAXT.EQ.2.OR.ZGN.LT.0d0) THEN!interior or inflow
C           Assemble into residual
            CALL DAXPY(NSTB,ZGN,
     '        PG(1,1,ng,nbh_f),1,RDF(1,1,nhx),1)
          ENDIF

CC!!! needs improving
CC         Evaluate integrand ZGN (with weights) from derivatives
C          ZGN=0.0d0
C          DO neax=1,NEAXT
C            ZGN=ZGN+DZNM(NIEF(0),neax)*DZXI(NIEF(0),neax)
C          ENDDO
C          ZGN=ZGN*FACTOR
CC         Only derivative weights wrt out of face xi.
C          DO idoxf=1,IDOXFT(1)
C            CALL DAXPY(NSTB,DZNM(NIEF(0),idoxf)*ZGN,
C     '        PG(1,1,ng,nbh_f),1,RDF(1,idoxf,nhx),1)
C          ENDDO !idoxf

        ELSE !higher derivative discontinuity term
          DO nhx=1,NHF
            nh=NH_LOC(nhx,nx)
            nb=NBHF(nh)
            NSTB=NST(nb)
            DO idoxf=1,IDOXFT(nhx)
C             Evaluate integrands ZGN (without weights) by interpolating
C             the necessary dependent variable derivatives.
              ZGN=DDOT(NSTNAT,PG(1,1,ng,nb),1,ZDF(1,idoxf,1,nhx),1)
C             multiply by Galerkin and other weights.
              CALL DAXPY(NSTB,ZGN*RWG(1),
     '          PG(1,1,ng,nb),1,RDF(1,idoxf,nhx),1)
            ENDDO !idoxf
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
          NSTB=NST(nb)
          DO idoxf=1,IDOXFT(nhx)
            WRITE(OP_STRING,
     '        '(/'' RDF(ns,'',I1,'','',I1,''): '',10E12.3,(/10E12.3))')
     '        idoxf,nhx,(RDF(ns,idoxf,nhx),ns=1,NSTB)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZFRF30')
      RETURN
 9999 CALL ERRORS('ZFRF30',ERROR)
      CALL EXITS('ZFRF30')
      RETURN 1
      END


