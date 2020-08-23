      SUBROUTINE GAUSS11(IBT,IDO,INP,nb1j,nb1jp,nbqh,
     '  nbuh,NGAP,D,DET,DET_ADAPT,PG_J,PG_Q,PG_U,XIG,XIG_J,XIG_Q,
     '  XIG_U,XIMIN,ERROR,*)

C#### Subroutine: GAUSS11
C###  Description:
C###    GAUSS11 defines the Gaussian quadrature coords XIG
C###    and evaluates the basis function Gauss point array PG for
C###    Lagrange / Hermite  tensor product type basis function when
C###    using the adaptive scheme of TELLES.

C**** The array D(ni) is used to determine the transformation to be
C**** used in each xi direction.  The xi coordinates of the nearest
C**** point to the singularity in ne are stored in XIMIN(ni).
C**** NBBEM indicates the basis function that needs to be defined.
C**** nb1jp indicates the basis function that is used to calculate the
C**** new basis function (family basis function for geometric variable).
C**** Use DET(nb1jp,0,ng,2) to store the determinant for this
C**** transformation

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  nb1j,nb1jp,nbqh,nbuh,NGAP(NIM,NBM)
      REAL*8 D(2),DET(NBFM,0:NNM,NGM,6),DET_ADAPT(NBFM,0:NNM,NGM),
     '  PG_J(NSM,NUM,NGM),PG_Q(NSM,NUM,NGM),PG_U(NSM,NUM,NGM),
     '  XIG(NIM,NGM,NBM),XIG_J(NIM,NGM),XIG_Q(NIM,NGM),XIG_U(NIM,NGM),
     '  XIMIN(2)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,nbuhp,nbqhp,ng,ng1,ng2,ng3,ni,nk,nn,ns,nu,NU1(0:3)
      REAL*8 A,B,BIGQ,C,DD,ETA,ETABAR,FACTOR,
     '  GAMMA,GAMMABAR,P,PSI1,PSI2_HERMITE,PSI5,Q,RBAR
      LOGICAL SECTOR

      DATA NU1/1,2,4,7/

C cpb 23/1/97 Adding multiprocssing adaptive integration. Store the
C transformed Gauss points in XIG_X and basis functions in PG_X where
C X=J for the geometry, Q for the normal derivative and U for the
C dependent variable. Store the transformed determinant in DET_ADAPT.

      CALL ENTERS('GAUSS11',*9999)

      nbuhp=NFBASE(1,nbuh)
      nbqhp=NFBASE(1,nbqh)

C***  Determine transformation in each direction
      DO ni=1,NIT(nb1jp)
        IF(D(ni).LT.0.05d0) THEN !Telles scheme not valid
          RBAR=1.0d0
C MLB July 2001 - Multiprocessing crashes sometimes when 2 strings
C  try to be written simultaneously so turning off warning with
C  multiprocessing for now.
C$        IF(.FALSE.) THEN
          WRITE(OP_STRING,'('' >>Warning: Distance too small for '
     '     //'adaptive integration scheme'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C$        ENDIF
        ELSE IF(D(ni).LE.1.3d0)THEN
          RBAR=0.85d0+0.24d0*DLOG(D(ni))
        ELSE IF(D(ni).LE.3.618d0)THEN
          RBAR=0.893d0+0.0832d0*DLOG(D(ni))
        ELSE
          RBAR=1.0d0
        ENDIF
        IF(RBAR.LT.1.0d0) THEN !Trans. is the identity if Rbar = 1
          ETABAR=2.0d0*XIMIN(ni)-1.0d0 !Put in (-1,1) coordinates
C*** Calculate the Telles parameters
          P=(4.0d0*RBAR*(1.0d0-RBAR)+3.0d0*(1.0d0-RBAR*RBAR))/
     '      (3.0d0*(1.0d0+RBAR+RBAR)*(1.0d0+RBAR+RBAR))
          Q=(((ETABAR*(3.0d0-RBAR-RBAR)-2.0d0*ETABAR**3/
     '       (1.0d0+RBAR+RBAR))/(1.0D0+RBAR+RBAR)-ETABAR))/
     '      (2.0D0+4.0D0*RBAR)
          FACTOR=DSQRT(Q*Q+P*P*P)
          GAMMABAR=(FACTOR-Q)**(1.0d0/3.0d0)-(FACTOR+Q)**(1.0d0/3.0d0)+
     '      ETABAR/(1.0d0+RBAR+RBAR)
          BIGQ=1.0d0+3.0d0*GAMMABAR*GAMMABAR
          A=(1.0d0-RBAR)/BIGQ
          B=-3.0d0*(1.0d0-RBAR)*GAMMABAR/BIGQ
          C=(RBAR+3.0d0*GAMMABAR*GAMMABAR)/BIGQ
          DD=-B
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' Xi='',I1,'' Telles parameters:'')') ni
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' D='',D12.4,'', Rbar='',D12.4,'
     '        //''', Etabar='',D12.4)') D(ni),RBAR,ETABAR
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' p='',D12.4,'', q='',D12.4)') P,Q
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Gammabar='',D12.4,'', Q='',D12.4)')
     '        GAMMABAR,BIGQ
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' a='',D12.4,'', b='',D12.4,'', c='','
     '        //'D12.4,'' d='',D12.4)') A,B,C,DD
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(ni.EQ.1)THEN
            DO ng=1,NGAP(1,nb1jp)
              GAMMA=2.0d0*XIG(1,ng,nb1jp)-1.0d0 !Change to (-1,1) limits
              ETA=((A*GAMMA+B)*GAMMA+C)*GAMMA+DD
              DET_ADAPT(nb1jp,0,ng)=DABS(((3.0d0*A*GAMMA+B+B)*GAMMA+C))*
     '          DET(nb1jp,0,ng,1)
              XIG_J(1,ng)=(ETA+1.0d0)/2.0d0 !Change to (0,1) limits

C LKC 23-AUG-2000 This is poked. Should never be looping over NIM
C  using nbuhp (should be the same as nbqhp)
C              IF(NIM.GT.1) THEN
              IF(NIT(nbuhp).GT.1) THEN
                DO ng1=1,NGAP(2,nb1jp)-1
                  ng2=ng+ng1*NGAP(1,nb1jp)
                  XIG_J(1,ng2)=XIG_J(1,ng)
                  DET_ADAPT(nb1jp,0,ng2)=DET_ADAPT(nb1jp,0,ng)
                ENDDO !ng1
              ENDIF
            ENDDO !ng
            IF(nbuhp.EQ.nb1jp) THEN
              DO ng=1,NGM
                XIG_U(1,ng)=XIG_J(1,ng)
              ENDDO !ng
            ELSE
              DO ng=1,NGAP(1,nbuhp)
                GAMMA=2.0d0*XIG(1,ng,nbuhp)-1.0d0 !Chng to (-1,1) limits
                ETA=((A*GAMMA+B)*GAMMA+C)*GAMMA+DD
                DET_ADAPT(nbuhp,0,ng)=DABS(((3.0d0*A*GAMMA+B+B)*GAMMA+
     '            C))*DET(nbuhp,0,ng,1)
                XIG_U(1,ng)=(ETA+1.0d0)/2.0d0 !Change to (0,1) limits
C LKC 23-AUG-2000 This is poked. Should never be looping over NIM
C  using nbuhp (should be the same as nbqhp)
C              IF(NIM.GT.1) THEN
                IF(NIT(nbuhp).GT.1) THEN
                  DO ng1=1,NGAP(2,nbuhp)-1
                    ng2=ng+ng1*NGAP(1,nbuhp)
                    XIG_U(1,ng2)=XIG_U(1,ng)
                    DET_ADAPT(nbuhp,0,ng2)=DET_ADAPT(nbuhp,0,ng)
                  ENDDO !ng1
                ENDIF
              ENDDO !ng
            ENDIF
            IF(nbqhp.EQ.nb1jp) THEN
              DO ng=1,NGM
                XIG_Q(1,ng)=XIG_J(1,ng)
              ENDDO !ng
            ELSE IF(nbqhp.EQ.nbuhp) THEN
              DO ng=1,NGM
                XIG_Q(1,ng)=XIG_U(1,ng)
              ENDDO !ng
            ELSE
              DO ng=1,NGAP(1,nbqhp)
                GAMMA=2.0d0*XIG(1,ng,nbqhp)-1.0d0 !Chng to (-1,1) limits
                ETA=((A*GAMMA+B)*GAMMA+C)*GAMMA+DD
                DET_ADAPT(nbqhp,0,ng)=DABS(((3.0d0*A*GAMMA+B+B)*GAMMA+
     '            C))*DET(nbqhp,0,ng,1)
                XIG_Q(1,ng)=(ETA+1.0d0)/2.0d0 !Change to (0,1) limits
C LKC 23-AUG-2000 This is poked. Should never be looping over NIM
C  using nbuhp (should be the same as nbqhp)
                IF(NIT(nbuhp).GT.1) THEN
                  DO ng1=1,NGAP(2,nbqhp)-1
                    ng2=ng+ng1*NGAP(1,nbqhp)
                    XIG_Q(1,ng2)=XIG_Q(1,ng)
                    DET_ADAPT(nbqhp,0,ng2)=DET_ADAPT(nbqhp,0,ng)
                  ENDDO !ng1
                ENDIF
              ENDDO !ng
            ENDIF
          ELSE !Need to get xi coordinates from the correct place.
            DO ng1=1,NGAP(ni,nb1jp)
              ng=1+(ng1-1)*NGAP(1,nb1jp)
              GAMMA=2.0d0*XIG(ni,ng,nb1jp)-1.0d0 !Chnge to (-1,1) limits
              ETA=((A*GAMMA+B)*GAMMA+C)*GAMMA+DD
              DET_ADAPT(nb1jp,0,ng)=DABS(((3.0d0*A*GAMMA+B+B)*GAMMA+C))*
     '          DET_ADAPT(nb1jp,0,ng)
              XIG_J(ni,ng)=(ETA+1.0d0)/2.0d0 !Chnge to (0,1) coord
              DO ng2=2,NGAP(1,nb1jp) !Fill in other entries.
                ng3=ng+ng2-1
                XIG_J(ni,ng3)=XIG_J(ni,ng)
                DET_ADAPT(nb1jp,0,ng3)=DABS(((3.0d0*A*GAMMA+B+B)*
     '            GAMMA+C))*DET_ADAPT(nb1jp,0,ng3)
              ENDDO !ng2
            ENDDO !ng1
            IF(nbuhp.EQ.nb1jp) THEN
              DO ng=1,NGM
                XIG_U(ni,ng)=XIG_J(ni,ng)
              ENDDO !ng
            ELSE
              DO ng1=1,NGAP(ni,nbuhp)
                ng=1+(ng1-1)*NGAP(1,nbuhp)
                GAMMA=2.0d0*XIG(ni,ng,nbuhp)-1.0d0 !Chg to (-1,1) limits
                ETA=((A*GAMMA+B)*GAMMA+C)*GAMMA+DD
                DET_ADAPT(nbuhp,0,ng)=DABS(((3.0d0*A*GAMMA+B+B)*GAMMA+
     '            C))*DET_ADAPT(nbuhp,0,ng)
                XIG_U(ni,ng)=(ETA+1.0d0)/2.0d0 !Chnge to (0,1) coord
                DO ng2=2,NGAP(1,nbuhp) !Fill in other entries.
                  ng3=ng+ng2-1
                  XIG_U(ni,ng3)=XIG_U(ni,ng)
                  DET_ADAPT(nbuhp,0,ng3)=DABS(((3.0d0*A*GAMMA+B+B)*
     '              GAMMA+C))*DET_ADAPT(nbuhp,0,ng3)
                ENDDO !ng2
              ENDDO !ng1
            ENDIF
            IF(nbqhp.EQ.nb1jp) THEN
              DO ng=1,NGM
                XIG_Q(ni,ng)=XIG_J(ni,ng)
              ENDDO !ng
            ELSE IF(nbqhp.EQ.nbuhp) THEN
              DO ng=1,NGM
                XIG_Q(ni,ng)=XIG_U(ni,ng)
              ENDDO !ng
            ELSE
              DO ng1=1,NGAP(ni,nbqhp)
                ng=1+(ng1-1)*NGAP(1,nbqhp)
                GAMMA=2.0d0*XIG(ni,ng,nbqhp)-1.0d0 !Chg to (-1,1) limits
                ETA=((A*GAMMA+B)*GAMMA+C)*GAMMA+DD
                DET_ADAPT(nbqhp,0,ng)=DABS(((3.0d0*A*GAMMA+B+B)*GAMMA+
     '            C))*DET_ADAPT(nbqhp,0,ng)
                XIG_Q(ni,ng)=(ETA+1.0d0)/2.0d0 !Chnge to (0,1) coord
                DO ng2=2,NGAP(1,nbqhp) !Fill in other entries.
                  ng3=ng+ng2-1
                  XIG_Q(ni,ng3)=XIG_Q(ni,ng)
                  DET_ADAPT(nbqhp,0,ng3)=DABS(((3.0d0*A*GAMMA+B+B)*
     '              GAMMA+C))*DET_ADAPT(nbqhp,0,ng3)
                ENDDO !ng2
              ENDDO !ng1
            ENDIF
          ENDIF
        ELSE !Rbar=1
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' Xi='',I1,'' Telles parameters:'')') ni
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' D='',D12.4,'', Rbar='',D12.4)') D(ni),
     '        RBAR
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(ni.EQ.1)THEN
            DO ng=1,NGAP(1,nb1jp)
              XIG_J(1,ng)=XIG(1,ng,nb1jp)
              DET_ADAPT(nb1jp,0,ng)=DET(nb1jp,0,ng,1)
C MLB July 2001 - Leo missed one! Needs the IF to stop 2d
C   problems crashing
              IF(NIT(nb1jp).GT.1) THEN
                DO ng1=1,NGAP(2,nb1jp)-1
                  ng2=ng+ng1*NGAP(1,nb1jp)
                  XIG_J(1,ng2)=XIG_J(1,ng)
                  DET_ADAPT(nb1jp,0,ng2)=DET_ADAPT(nb1jp,0,ng)
                ENDDO !ng1
              ENDIF
            ENDDO !ng
            IF(nbuhp.EQ.nb1jp) THEN
              DO ng=1,NGM
                XIG_U(1,ng)=XIG_J(1,ng)
              ENDDO !ng
            ELSE
              DO ng=1,NGAP(1,nbuhp)
                XIG_U(1,ng)=XIG(1,ng,nbuhp)
                DET_ADAPT(nbuhp,0,ng)=DET(nbuhp,0,ng,1)
                DO ng1=1,NGAP(2,nbuhp)-1
                  ng2=ng+ng1*NGAP(1,nbuhp)
                  XIG_U(1,ng2)=XIG_U(1,ng)
                  DET_ADAPT(nbuhp,0,ng2)=DET_ADAPT(nbuhp,0,ng)
                ENDDO !ng1
              ENDDO !ng
            ENDIF
            IF(nbqhp.EQ.nb1jp) THEN
              DO ng=1,NGM
                XIG_Q(1,ng)=XIG_J(1,ng)
              ENDDO !ng
            ELSE IF(nbqhp.EQ.nbuhp) THEN
              DO ng=1,NGM
                XIG_Q(1,ng)=XIG_U(1,ng)
              ENDDO !ng
            ELSE
              DO ng=1,NGAP(1,nbqhp)
                XIG_Q(1,ng)=XIG(1,ng,nbqhp)
                DET_ADAPT(nbqhp,0,ng)=DET(nbqhp,0,ng,1)
                DO ng1=1,NGAP(2,nbqhp)-1
                  ng2=ng+ng1*NGAP(1,nbqhp)
                  XIG_Q(1,ng2)=XIG_Q(1,ng)
                  DET_ADAPT(nbqhp,0,ng2)=DET_ADAPT(nbqhp,0,ng)
                ENDDO !ng1
              ENDDO !ng
            ENDIF
          ELSE !Need to get xi coordinates from the correct place
            DO ng1=1,NGAP(ni,nb1jp)
              ng=1+(ng1-1)*NGAP(1,nb1jp)
              XIG_J(ni,ng)=XIG(ni,ng,nb1jp)
              DO ng2=2,NGAP(1,nb1jp) !Fill in other entries.
                ng3=ng+ng2-1
                XIG_J(ni,ng3)=XIG_J(ni,ng)
              ENDDO !ng2
            ENDDO !ng1
            IF(nbuhp.EQ.nb1jp) THEN
              DO ng=1,NGM
                XIG_U(ni,ng)=XIG_J(ni,ng)
              ENDDO !ng
            ELSE
              DO ng1=1,NGAP(ni,nbuhp)
                ng=1+(ng1-1)*NGAP(1,nbuhp)
                XIG_U(ni,ng)=XIG(ni,ng,nbuhp)
                DO ng2=2,NGAP(1,nbuhp) !Fill in other entries.
                  ng3=ng+ng2-1
                  XIG_U(ni,ng3)=XIG_U(ni,ng)
                ENDDO !ng2
              ENDDO !ng1
            ENDIF
            IF(nbqhp.EQ.nb1jp) THEN
              DO ng=1,NGM
                XIG_Q(ni,ng)=XIG_J(ni,ng)
              ENDDO !ng
            ELSE IF(nbqhp.EQ.nbuhp) THEN
              DO ng=1,NGM
                XIG_Q(ni,ng)=XIG_U(ni,ng)
              ENDDO !ng
            ELSE
              DO ng1=1,NGAP(ni,nbqhp)
                ng=1+(ng1-1)*NGAP(1,nbqhp)
                XIG_Q(ni,ng)=XIG(ni,ng,nbqhp)
                DO ng2=2,NGAP(1,nbqhp) !Fill in other entries.
                  ng3=ng+ng2-1
                  XIG_Q(ni,ng3)=XIG_Q(ni,ng)
                ENDDO !ng2
              ENDDO !ng1
            ENDIF
          ENDIF
        ENDIF
      ENDDO !End of ni loop
      ng1=NGAP(1,nb1jp)
      IF(NIT(nb1jp).GT.1)THEN
        ng2=NGAP(2,nb1jp)
      ELSE
        ng2=1
      ENDIF
      SECTOR=.FALSE.
      DO ni=1,NIT(nb1jp)
        IF(IBT(1,ni,nb1jp).EQ.5.OR.IBT(1,ni,nb1jp).EQ.6) SECTOR=.TRUE.
      ENDDO !ni
      IF(IBT(1,1,nb1jp).LE.2.OR.SECTOR)THEN !Lagrange or Hermite
        DO J=1,ng2
          DO I=1,ng1
            ng=I+(J-1)*ng1
            ns=0
            DO nn=1,NNT(nb1jp)
              DO nk=1,NKT(nn,nb1jp)
                ns=ns+1
C cpb 26/5/95 Only need the nu=1,2,4 cases
C                DO nu=1,NUT(nb1jp)
                DO ni=0,NIT(nb1jp)
                  IF(SECTOR) THEN
                    PG_J(ns,NU1(ni),ng)=PSI5(IBT(1,1,nb1jp),
     '                IDO(1,1,0,nb1jp),INP(1,1,nb1jp),nb1jp,NU1(ni),
     '                nk,nn,XIG_J(1,ng))
                  ELSE
                    PG_J(ns,NU1(ni),ng)=PSI1(IBT(1,1,nb1jp),
     '                IDO(1,1,0,nb1jp),INP(1,1,nb1jp),nb1jp,NU1(ni),
     '                nk,nn,XIG_J(1,ng))
                  ENDIF
                ENDDO !ni
              ENDDO !nk
            ENDDO !nn
          ENDDO !I
        ENDDO !J
        IF(nbuh.NE.nb1j)THEN
          DO J=1,ng2
            DO I=1,ng1
              ng=I+(J-1)*ng1
              ns=0
              DO nn=1,NNT(nbuhp)
                DO nk=1,NKT(nn,nbuhp)
                  ns=ns+1
C cpb 26/5/95 Only need the nu=1,2,4 cases
C                  DO nu=1,NUT(nbuhp)
                  DO ni=0,NIT(nbuhp)
                    IF(SECTOR) THEN
                      PG_U(ns,NU1(ni),ng)=PSI5(IBT(1,1,nbuhp),
     '                  IDO(1,1,0,nbuhp),INP(1,1,nbuhp),nbuhp,NU1(ni),
     '                  nk,nn,XIG_U(1,ng))
                    ELSE
                      PG_U(ns,NU1(ni),ng)=PSI1(IBT(1,1,nbuhp),
     '                  IDO(1,1,0,nbuhp),INP(1,1,nbuhp),nbuhp,NU1(ni),
     '                  nk,nn,XIG_U(1,ng))
                    ENDIF
                  ENDDO !ni
                ENDDO !nk
              ENDDO !nn
            ENDDO !I
          ENDDO !J
        ELSE !copy the J case
          DO ng=1,ng1*ng2
            DO ni=0,NIT(nb1jp)
              nu=NU1(ni)
              DO ns=1,NST(nb1jp)
                PG_U(ns,nu,ng)=PG_J(ns,nu,ng)
              ENDDO !ns
            ENDDO !ni
          ENDDO !ng
        ENDIF
        IF(nbqh.NE.nb1j)THEN
          DO J=1,ng2
            DO I=1,ng1
              ng=I+(J-1)*ng1
              ns=0
              DO nn=1,NNT(nbqhp)
                DO nk=1,NKT(nn,nbqhp)
                  ns=ns+1
C cpb 26/5/95 Only need the nu=1,2,4 cases
C                  DO nu=1,NUT(nbqhp)
                  DO ni=0,NIT(nbqhp)
                    IF(SECTOR) THEN
                      PG_Q(ns,NU1(ni),ng)=PSI5(IBT(1,1,nbqhp),
     '                  IDO(1,1,0,nbqhp),INP(1,1,nbqhp),nbqhp,NU1(ni),
     '                  nk,nn,XIG_Q(1,ng))
                    ELSE
                      PG_Q(ns,NU1(ni),ng)=PSI1(IBT(1,1,nbqhp),
     '                  IDO(1,1,0,nbqhp),INP(1,1,nbqhp),nbqhp,NU1(ni),
     '                  nk,nn,XIG_Q(1,ng))
                    ENDIF
                  ENDDO !ni
                ENDDO !nk
              ENDDO !nn
            ENDDO !I
          ENDDO !J
        ELSE !copy the J case
          DO ng=1,ng1*ng2
            DO ni=0,NIT(nb1jp)
              nu=NU1(ni)
              DO ns=1,NST(nb1jp)
                PG_Q(ns,nu,ng)=PG_J(ns,nu,ng)
              ENDDO !ns
            ENDDO !ni
          ENDDO !ng
        ENDIF
      ELSE IF(IBT(1,1,nb1jp).EQ.3)THEN !Simplex
        IF(IBT(2,1,nb1jp).EQ.4)THEN !Special hermite simplex
          DO J=1,ng2
            DO I=1,ng1
              ng=I+(J-1)*ng1
              ns=0
              DO nn=1,NNT(nb1jp)
                DO nk=1,NKT(nn,nb1jp)
                  ns=ns+1
C cpb 26/5/95 Only need the nu=1,2,4 cases
C                  DO nu=1,NUT(nb1jp)
                  DO ni=0,NIT(nb1jp)
                    PG_J(ns,NU1(ni),ng)=PSI2_HERMITE(
     '                IDO(1,1,0,nb1jp),INP(1,1,nb1jp),nb1jp,NU1(ni),
     '                nk,nn,XIG_J(1,ng))
                  ENDDO !ni
                ENDDO !nk
              ENDDO !nn
            ENDDO !I
          ENDDO !J
          IF(nbuh.NE.nb1j)THEN
            DO J=1,ng2
              DO I=1,ng1
                ng=I+(J-1)*ng1
                ns=0
                DO nn=1,NNT(nbuhp)
                  DO nk=1,NKT(nn,nbuhp)
                    ns=ns+1
C cpb 26/5/95 Only need the nu=1,2,4 cases
C                    DO nu=1,NUT(nbuhp)
                    DO ni=0,NIT(nbuhp)
                      PG_U(ns,NU1(ni),ng)=PSI2_HERMITE(
     '                  IDO(1,1,0,nbuhp),INP(1,1,nbuhp),
     '                  nbuhp,NU1(ni),nk,nn,XIG_U(1,ng))
                    ENDDO !ni
                  ENDDO !nk
                ENDDO !nn
              ENDDO !I
            ENDDO !J
          ELSE !copy the J case
            DO ng=1,NGM
              DO nu=1,NUM
                DO ns=1,NSM
                  PG_U(ns,nu,ng)=PG_J(ns,nu,ng)
                ENDDO !ns
              ENDDO !nu
            ENDDO !ng
          ENDIF
          IF(nbqh.NE.nb1j)THEN
            DO J=1,ng2
              DO I=1,ng1
                ng=I+(J-1)*ng1
                ns=0
                DO nn=1,NNT(nbqhp)
                  DO nk=1,NKT(nn,nbqhp)
                    ns=ns+1
C cpb 26/5/95 Only need the nu=1,2,4 cases
C                    DO nu=1,NUT(nbqhp)
                    DO ni=1,NIT(nbqhp)
                      PG_Q(ns,NU1(ni),ng)=PSI2_HERMITE(
     '                  IDO(1,1,0,nbqhp),INP(1,1,nbqhp),
     '                  nbqhp,NU1(ni),nk,nn,XIG_Q(1,ng))
                    ENDDO !ni
                  ENDDO !nk
                ENDDO !nn
              ENDDO !I
            ENDDO !J
          ELSE !copy the J case
            DO ng=1,NGM
              DO nu=1,NUM
                DO ns=1,NSM
                  PG_Q(ns,nu,ng)=PG_J(ns,nu,ng)
                ENDDO !ns
              ENDDO !nu
            ENDDO !ng
          ENDIF
        ELSE
          WRITE(OP_STRING,'('' >>Simplex code not implemented'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ELSE
        WRITE(OP_STRING,'('' >>Add code to GAUSS11 for this element'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('GAUSS11')
      RETURN
9999  CALL ERRORS('GAUSS11',ERROR)
      CALL EXITS('GAUSS11')
      RETURN 1
      END


