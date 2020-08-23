      SUBROUTINE XGYG(NITB,nr,nx,CG,WG,XG,YG,MODDIFF,ERROR,*)

C#### Subroutine: XGYG
C###  Description:
C###    Given geometry XG and material parameters CG, XGYG calculates
C###    the appropriate expressions required for evaluation of a
C###    weighted residual of a differential equation, and stores the
C###    results in YG.  The purpose of the routine is to calculate
C###    expressions that need only be done once but are used several
C###    times in an iterative procedure.  The expressions do not depend
C###    on the current value of the solution but enable the required
C###    terms of the differential equation to be calculated from local
C###    element xi derivatives of the dependent variable.  The
C###    expressions calculated here do not depend on the interpolation
C###    functions (PG) as such terms would require a large amount of
C###    memory for storage.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NITB,nr,nx
      REAL*8 CG(NMM),WG,XG(NJM,NUM),YG(NIYGM)
      LOGICAL MODDIFF
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IEND,mi,mjj,ni,niyg,nj,njj,NJTR,NU1(3)
      REAL*8 COUP(3,3),CxXI_FIBRE(3,3),DIFFM1(3),DIFFMODM(3),
     '  DXIXRC(3,3),DXRCXI(3,3),D2XRCXI(3,3,3),D_XRC_FIBRE(3,3,3),
     '  RG,RWG,SUM,TEMP,TEMP2,U(3),XI_FIBRE(3,3),XRC_FIBRE(3,3)
      LOGICAL ISOTROPIC,PETROV

      DATA NU1/2,4,7/

      CALL ENTERS('XGYG',*9999)

      niyg=0
      NJTR=NJ_LOC(NJL_GEOM,0,nr)
      IF(ITYP3(nr,nx).EQ.1) THEN !isotropic
        ISOTROPIC=.TRUE.
      ELSE IF(ITYP3(nr,nx).EQ.2) THEN !anisotropic monodomain
        ISOTROPIC=.FALSE.
      ELSE
        ERROR='>>Only monodomain materials are implemented'
        GO TO 9999
      ENDIF !iso/aniso-tropic
      PETROV=ITYP15(nr,nx).EQ.1.OR.ITYP15(nr,nx).EQ.2

C     Calculate derivatives of r.c. coords wrt xi coords.
      IF(PETROV) THEN !need 2nd derivs
        CALL D2XRCDXI(NITB,nr,DXRCXI,D2XRCXI,XG,ERROR,*9999)
      ELSE
        CALL DXRCDXI(NITB,nr,DXRCXI,XG,ERROR,*9999)
      ENDIF

C     Invert dXrc/dxi.
C***  The determinant is the integration Jacobian.
      CALL INVERT(NJTR,DXRCXI,DXIXRC,RG)
      RWG=RG*WG
      niyg=niyg+1
      IF(niyg.LE.NIYGM) YG(niyg)=RWG

C*** Elliptic eikonal equation for wavefront path determination.
C***   CG(3..,ng) are space constants.
      IF(ISOTROPIC) THEN !isotropic
C***    Calculate coupling tensor in xi coords for advection term and
C***    for multiplying second xi derivs in diffusion term.
        TEMP=CG(3)**2
        DO ni=1,NITB
          DO mi=ni,NITB
            SUM=0.0d0
            DO njj=1,NJTR
              SUM=SUM+DXIXRC(mi,njj)*DXIXRC(ni,njj)
            ENDDO !njj
C           GU not used see below.
C            GU(mi,ni)=SUM
            COUP(mi,ni)=SUM*TEMP
          ENDDO !mi
          DO mi=1,ni-1
C           GU not used see below.
C            GU(mi,ni)=GU(ni,mi)
            COUP(mi,ni)=COUP(ni,mi)
          ENDDO !mi
        ENDDO !ni
      ELSE !anisotropic monodomain
C       Calculate material vectors in r.c. coords
        IF(PETROV) THEN
          CALL D_MAT_VEC(NITB,nr,DXRCXI,D2XRCXI,
     '      D_XRC_FIBRE,XRC_FIBRE,XG,ERROR,*9999)
        ELSE
          CALL MAT_VEC(NITB,nr,
     '      XRC_FIBRE(1,1),XRC_FIBRE(1,2),XRC_FIBRE(1,3),
     '      DXRCXI,XG,ERROR,*9999)
        ENDIF
C       Transform into xi coords (dxi/dnu).
        DO mjj=1,NITB
          DO ni=1,NITB
            XI_FIBRE(ni,mjj)=0.0d0
          ENDDO !ni
          DO njj=1,NJTR
            DO ni=1,NITB
              XI_FIBRE(ni,mjj)=
     '          XI_FIBRE(ni,mjj)+DXIXRC(ni,njj)*XRC_FIBRE(njj,mjj)
            ENDDO !ni
          ENDDO !mjj
        ENDDO !njj
C***    Calculate coupling tensor in xi coords for advection term and
C***    for multiplying second xi derivs in diffusion term.
        DO mjj=1,NITB
          TEMP=CG(2+mjj)**2
          DO ni=1,NITB
            CxXI_FIBRE(ni,mjj)=TEMP*XI_FIBRE(ni,mjj)
          ENDDO !ni
        ENDDO !mjj
        DO ni=1,NITB
          DO mi=1,ni-1
            COUP(mi,ni)=COUP(ni,mi)
          ENDDO !mi
          DO mi=ni,NITB
            SUM=0.0d0
            DO mjj=1,NITB
              SUM=SUM+XI_FIBRE(mi,mjj)*CxXI_FIBRE(ni,mjj)
            ENDDO !mjj
            COUP(mi,ni)=SUM
          ENDDO !mi
        ENDDO !ni
      ENDIF !iso/aniso-tropic
      DO mi=1,NITB
        DO ni=1,mi
          niyg=niyg+1
          IF(niyg.LE.NIYGM) YG(niyg)=COUP(mi,ni)
        ENDDO !ni
      ENDDO !mi

C***  Calculate multiplier of first derivs in diffusion

      IF(PETROV) THEN !Petrov-Galerkin
C       Initialize multiplier of first derivs for diffusion term.
        DO ni=1,NITB
          DIFFM1(ni)=0.0d0
        ENDDO !ni
C       Initialize temporary vector U(njj) to zero
        DO njj=1,NJTR
          U(njj)=0.0d0
        ENDDO !njj
        IF(ISOTROPIC) THEN !isotropic
C           This is to handle derivatives in coupling coefficients but
C           at present we don't know derivatves of material parameters.
C            DO mi=1,NITB
C              DO ni=1,NITB
C                DIFFM1(ni)=DIFFM1(ni)+D_CG(3,mi)*GU(ni,mi)
C              ENDDO !ni
C            ENDDO !mi
        ELSE !anisotropic
C           This is to handle derivatives in coupling coefficients but
C           at present we don't know derivatves of material parameters.
C            DO mjj=1,NITB
C              SUM=0.0d0
C              DO ni=1,NITB
C                SUM=SUM+D_CG(2+mjj,ni)*XI_FIBRE(ni,mjj)
C              ENDDO !ni
C              DO ni=1,NITB
C                DIFFM1(ni)=DIFFM1(ni)+XI_FIBRE(ni,mjj)*SUM
C              ENDDO !ni
C            ENDDO !mjj
C           This is to handle derivatives of fibre directions.
          DO mjj=1,NITB
            SUM=0.0d0
            DO ni=1,NITB
              DO njj=1,NJTR
                SUM=SUM+DXIXRC(ni,njj)*D_XRC_FIBRE(njj,mjj,ni)
              ENDDO !njj
            ENDDO !ni
            DO ni=1,NITB
              DIFFM1(ni)=DIFFM1(ni)+CxXI_FIBRE(ni,mjj)*SUM
            ENDDO !ni
          ENDDO !mjj
          DO ni=1,NITB
            DO mjj=1,NITB
              DO njj=1,NJTR
                U(njj)=
     '            U(njj)+CxXI_FIBRE(ni,mjj)*D_XRC_FIBRE(njj,mjj,ni)
              ENDDO !njj
            ENDDO !mjj
          ENDDO !ni
        ENDIF !iso/aniso-tropic
C       Include in temporary vector U(njj) the product of xi-coord
C       coupling tensor and second derivatives of r.c. coords wrt xi.
        DO ni=1,NITB
          DO mi=1,NITB
            DO njj=1,NJTR
              U(njj)=U(njj)-COUP(mi,ni)*D2XRCXI(njj,mi,ni)
            ENDDO !njj
          ENDDO !mi
        ENDDO !ni
C       Multiply U(njj) by dxi/dXrc
C       and include in the diffusion multiplier.
        DO njj=1,NJTR
          DO ni=1,NITB
            DIFFM1(ni)=DIFFM1(ni)+DXIXRC(ni,njj)*U(njj)
          ENDDO !ni
        ENDDO !njj
        DO ni=1,NITB
          niyg=niyg+1
          IF(niyg.LE.NIYGM) YG(niyg)=DIFFM1(ni)*RWG
        ENDDO !ni
      ENDIF !Petrov-Galerkin

      IF(MODDIFF) THEN
C       Initialize multiplier of first derivs for diffusion modifier
        DO ni=1,NITB
          DIFFMODM(ni)=0.0d0
        ENDDO !ni
        nj=NJ_LOC(NJL_FIEL,NJ_LOC(NJL_FIEL,0,nr),nr)
        TEMP=XG(nj,1)
        IF(TEMP.LT.1.0d0-ZERO_TOL) THEN
          CALL ASSERT(TEMP.NE.0.0d0,
     '      '>>Zero diffusion modifier',ERROR,*9999)
          DO ni=1,NITB
            TEMP2=XG(nj,NU1(ni))/TEMP
            DO mi=1,NITB
              DIFFMODM(mi)=DIFFMODM(mi)-COUP(mi,ni)*TEMP2
            ENDDO !mi
          ENDDO !ni
        ENDIF !modifying parameter not 1
        DO ni=1,NITB
          niyg=niyg+1
          IF(niyg.LE.NIYGM) YG(niyg)=DIFFMODM(ni)*RWG
        ENDDO !ni
      ENDIF !MODDIFF

      IF(niyg.GT.NIYGM) THEN
        IEND=0
        CALL APPENDC(IEND,' >>Increase NIYGM to ',ERROR)
        CALL APPENDI(IEND,niyg,ERROR)
        GOTO 9999
      ENDIF
      CALL EXITS('XGYG')
      RETURN
 9999 CALL ERRORS('XGYG',ERROR)
      CALL EXITS('XGYG')
      RETURN 1
      END


