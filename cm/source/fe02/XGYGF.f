      SUBROUTINE XGYGF(NEAXT,NIEF,NITF,nr,nx,CG,WG,XGF,YGF,ERROR,*)

C#### Subroutine: XGYGF
C###  Description:
C###    Given geometry XG and material parameters CG, XGYGF calculates
C###    the appropriate terms required for evaluation of a
C###    derivative expression on a face, and stores the
C###    results in YGF.  The purpose of the routine is to calculate
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
!     Parameter List
      INTEGER NEAXT,NITF,nr,nx
      REAL*8 CG(NMM),WG,XGF(NJM,NUM,2),YGF(NIYGFM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IEND,mi,mjj,neax,ni,NIEF(0:2),NITE,niyg,njj,NJTR,
     '  NORMALSIGN(2)
      REAL*8 CG2SQ,COUP(3,3),DET,DXIXRC(3,3),DXRCXI(3,3),FIBRE_XI(3,3),
     '  MNCPT(2),RG,RWG,SUM,TEMP,XI_FIBRE(3,3),XN_NORMAL(3),
     '  XRC_FIBRE(3,3),XRC_NORMAL(3)
      LOGICAL FLUX,ISOTROPIC,PETROV

      DATA NORMALSIGN/1,-1/ !for opposite normals for each adjacent element

      CALL ENTERS('XGYGF',*9999)

      niyg=0
      NITE=NITF+1
      NJTR=NJ_LOC(NJL_GEOM,0,nr)
      CALL ASSERT(NJTR.EQ.NITE,'>>NITF!=NJTR-1',ERROR,*9999)

      IF(ITYP3(nr,nx).EQ.1) THEN !isotropic
        ISOTROPIC=.TRUE.
      ELSE IF(ITYP3(nr,nx).EQ.2) THEN !anisotropic monodomain
        ISOTROPIC=.FALSE.
      ELSE
        ERROR='>>Only monodomain materials are implemented'
        GO TO 9999
      ENDIF !iso/aniso-tropic
      PETROV=ITYP15(nr,nx).EQ.1.OR.ITYP15(nr,nx).EQ.2
C!!!  need to check order of stabilizing derivs
      FLUX=PETROV.OR..TRUE.

C*** Elliptic eikonal equation for wavefront path determination.
C***   CG(2) is the dimensionless coefficient of advection, co,
C***   CG(3) are space constants.

      IF(FLUX) THEN !Flux term

        DO neax=1,NEAXT !loop over adjacent elements

C         Calculate DXRCXI, derivative of rect. Cart. coords. wrt xi
          CALL DXRCDXI(NITE,nr,DXRCXI,XGF(1,1,neax),ERROR,*9999)

          IF(.NOT.ISOTROPIC) THEN
C           Calculate material vectors in r.c. coords
            CALL MAT_VEC(NITE,nr,
     '        XRC_FIBRE(1,1),XRC_FIBRE(1,2),XRC_FIBRE(1,3),
     '        DXRCXI,XGF(1,1,neax),ERROR,*9999)
          ENDIF

C         Invert dXrc/dxi.
C         The determinant is the integration Jacobian.
          CALL INVERT(NJTR,DXRCXI,DXIXRC,RG)

          IF(ISOTROPIC) THEN !isotropic
C***        Calculate coupling tensor in xi coords for advection term and
C***        for multiplying second xi derivs in diffusion term.
            CG2SQ=CG(3)**2
            DO ni=1,NITE
              DO mi=ni,NITE
                SUM=0.0d0
                DO njj=1,NJTR
                  SUM=SUM+DXIXRC(mi,njj)*DXIXRC(ni,njj)
                ENDDO !njj
                COUP(mi,ni)=SUM*CG2SQ
              ENDDO !mi
              DO mi=1,ni-1
                COUP(mi,ni)=COUP(ni,mi)
              ENDDO !mi
            ENDDO !ni
          ELSE !anisotropic monodomain
C           Transform material vectors into xi coords (dxi/dnu).
            DO mjj=1,NITE
              DO ni=1,NITE
                XI_FIBRE(ni,mjj)=0.0d0
              ENDDO !ni
              DO njj=1,NJTR
                DO ni=1,NITE
                  XI_FIBRE(ni,mjj)=
     '              XI_FIBRE(ni,mjj)+DXIXRC(ni,njj)*XRC_FIBRE(njj,mjj)
                ENDDO !ni
              ENDDO !mjj
            ENDDO !njj
C           Scale by length constants
            DO mjj=1,NITE
              TEMP=CG(2+mjj)
              DO ni=1,NITE
                XI_FIBRE(ni,mjj)=TEMP*XI_FIBRE(ni,mjj)
              ENDDO !ni
            ENDDO !mjj
C***        Calculate coupling tensor in xi coords for advection term and
C***        for multiplying second xi derivs in diffusion term.
            DO ni=1,NITE
              DO mi=1,ni-1
                COUP(mi,ni)=COUP(ni,mi)
              ENDDO !mi
              DO mi=ni,NITE
                SUM=0.0d0
                DO mjj=1,NITE
                  SUM=SUM+XI_FIBRE(mi,mjj)*XI_FIBRE(ni,mjj)
                ENDDO !mjj
                COUP(mi,ni)=SUM
              ENDDO !mi
            ENDDO !ni
          ENDIF !iso/aniso-tropic

C         Choose the normal derivative in the dirn of increasing
C         cross face xi. i.e. XINORMAL(NIEF(0))=1
          RWG=NORMALSIGN(neax)*DABS(RG)*WG

          IF(NEAXT.EQ.1.OR..NOT.PETROV) THEN
C***        Calculate square of ratio of element length to space
C***        constant in the direction of the normal.
            IF(ISOTROPIC) THEN !isotropic
C             Calculate surface normal (including Jacobian) from
C             tangents to face.
              IF(NITF.EQ.2) THEN !cross product of tangents
                XRC_NORMAL(1)=DXRCXI(2,NIEF(1))*DXRCXI(3,NIEF(2))-
     '            DXRCXI(3,NIEF(1))*DXRCXI(2,NIEF(2))
                XRC_NORMAL(2)=DXRCXI(3,NIEF(1))*DXRCXI(1,NIEF(2))-
     '            DXRCXI(1,NIEF(1))*DXRCXI(3,NIEF(2))
                XRC_NORMAL(3)=DXRCXI(1,NIEF(1))*DXRCXI(2,NIEF(2))-
     '            DXRCXI(2,NIEF(1))*DXRCXI(1,NIEF(2))
              ELSE IF(NITF.EQ.1) THEN !rotate tangent through a right angle
                XRC_NORMAL(2)=DXRCXI(1,NIEF(1))
                XRC_NORMAL(1)=-DXRCXI(2,NIEF(1))
              ELSE !(NITF).EQ.0) unit magnitude
                XRC_NORMAL(1)=1.0d0
              ENDIF
C             Square of surface Jacobian
              SUM=0.0d0
              DO njj=1,NJTR
                SUM=SUM+XRC_NORMAL(njj)**2
              ENDDO !nj
C             Ratio
              MNCPT(neax)=RG*RG/(SUM*CG2SQ)
            ELSE !anisotropic monodomain
              CALL INVERT(NJTR,XI_FIBRE,FIBRE_XI,DET)
C             Calculate normal (including Jacobian) in fibre coords
              IF(NITF.EQ.2) THEN !cross product of tangents
                XN_NORMAL(1)=FIBRE_XI(2,NIEF(1))*FIBRE_XI(3,NIEF(2))-
     '            FIBRE_XI(3,NIEF(1))*FIBRE_XI(2,NIEF(2))
                XN_NORMAL(2)=FIBRE_XI(3,NIEF(1))*FIBRE_XI(1,NIEF(2))-
     '            FIBRE_XI(1,NIEF(1))*FIBRE_XI(3,NIEF(2))
                XN_NORMAL(3)=FIBRE_XI(1,NIEF(1))*FIBRE_XI(2,NIEF(2))-
     '            FIBRE_XI(2,NIEF(1))*FIBRE_XI(1,NIEF(2))
              ELSE IF(NITF.EQ.1) THEN !rotate tangent through a right angle
                XN_NORMAL(2)=FIBRE_XI(1,NIEF(1))
                XN_NORMAL(1)=-FIBRE_XI(2,NIEF(1))
              ELSE !(NITF).EQ.0) unit magnitude
                XN_NORMAL(1)=1.0d0
              ENDIF
C             Square of surface Jacobian
              SUM=0.0d0
              DO njj=1,NJTR
                SUM=SUM+XN_NORMAL(njj)**2
              ENDDO !nj
C             Ratio
              MNCPT(neax)=1.0d0/(SUM*DET*DET)
            ENDIF !iso/aniso-tropic
          ENDIF

          IF(PETROV) THEN
            IF(NEAXT.EQ.1) THEN !boundary face
              niyg=niyg+1
              IF(niyg.LE.NIYGFM) YGF(niyg)=RWG
C             Coupling tensor in xi coords.
              DO mi=1,NITE
                DO ni=1,mi
                  niyg=niyg+1
                  IF(niyg.LE.NIYGFM) YGF(niyg)=COUP(mi,ni)
                ENDDO !ni
              ENDDO !mi
              niyg=niyg+1
              IF(niyg.LE.NIYGFM) YGF(niyg)=MNCPT(neax)
            ELSE !interior face
              niyg=niyg+1
              IF(niyg.LE.NIYGFM) YGF(niyg)=RWG
C             Residual multipliers (not including Gauss weights)
              DO mi=1,NITE
                niyg=niyg+1
                IF(niyg.LE.NIYGFM) YGF(niyg)=COUP(mi,NIEF(0))
              ENDDO !mi
            ENDIF !NEAXT
          ELSE !not PETROV
C           Residual multipliers (including Gauss weights)
            DO mi=1,NITE
              niyg=niyg+1
              IF(niyg.LE.NIYGFM) YGF(niyg)=
     '          COUP(mi,NIEF(0))*RWG
C     '          COUP(mi,NIEF(0))*NORMALSIGN(neax)*DSQRT(MNCPT(neax))
            ENDDO !mi
          ENDIF !PETROV

        ENDDO !neax

        IF(.NOT.PETROV) THEN
C         Include the Peclet number in the residual.
          TEMP=0.0d0
          DO neax=1,NEAXT
            TEMP=TEMP+DSQRT(MNCPT(neax))
          ENDDO !neax
CC         Include a balancing factor of 1/12 for cubic Hermite elements.
C          TEMP=CG(2)*TEMP/(12*NEAXT)
C         Include a balancing factor of 1/6 for quadratic Lagrange elements.
          TEMP=CG(2)*TEMP/(6*NEAXT)
CC         Include a factor of Peclet*MNCPT/(120sqrt(15)) for cubic
CC         Hermite elements.  Take the square root because factor is used
CC         for weights and residual.
C          TEMP=0.0d0
C          DO neax=1,NEAXT
C            TEMP=TEMP+RWG(neax)*DSQRT(MNCPT(neax))
C          ENDDO !neax
C          TEMP=CG(2)*TEMP/(120*DSQRT(15.0d0/2))
C          TEMP=DSQRT(TEMP)
          niyg=0
          DO neax=1,NEAXT !loop over adjacent elements
            DO ni=1,NITE
              niyg=niyg+1
              IF(niyg.LE.NIYGFM) YGF(niyg)=YGF(niyg)*TEMP
            ENDDO !ni
          ENDDO !neax
        ENDIF !not PETROV

      ELSE !Derivative discontinuity term

C       Calculate DXRCXI, derivative of rect. Cart. coords. wrt xi
        CALL DXRCDXI(NITE,nr,DXRCXI,XGF(1,1,1),ERROR,*9999)

        IF(.NOT.ISOTROPIC) THEN
C         Calculate material vectors in r.c. coords
          CALL MAT_VEC(NITE,nr,
     '      XRC_FIBRE(1,1),XRC_FIBRE(1,2),XRC_FIBRE(1,3),
     '      DXRCXI,XGF(1,1,1),ERROR,*9999)
        ENDIF

C       Calculate normal (including Jacobian) from tangents to face.
        IF(NITF.EQ.2) THEN !cross product of tangents
          XRC_NORMAL(1)=DXRCXI(2,NIEF(1))*DXRCXI(3,NIEF(2))-
     '      DXRCXI(3,NIEF(1))*DXRCXI(2,NIEF(2))
          XRC_NORMAL(2)=DXRCXI(3,NIEF(1))*DXRCXI(1,NIEF(2))-
     '      DXRCXI(1,NIEF(1))*DXRCXI(3,NIEF(2))
          XRC_NORMAL(3)=DXRCXI(1,NIEF(1))*DXRCXI(2,NIEF(2))-
     '      DXRCXI(2,NIEF(1))*DXRCXI(1,NIEF(2))
        ELSE IF(NITF.EQ.1) THEN !rotate tangent through a right angle
          XRC_NORMAL(2)=DXRCXI(1,NIEF(1))
          XRC_NORMAL(1)=-DXRCXI(2,NIEF(1))
        ELSE !(NITF).EQ.0) unit magnitude
          XRC_NORMAL(1)=1.0d0
        ENDIF

C       Evaluate weighting const = Gauss weight * c0 * sqrt(n.Dn): RWG
C       (includes jacobian in the normal)
        IF(ISOTROPIC) THEN !Isotropic
          SUM=0.0d0
          DO njj=1,NJTR
            SUM=SUM+XRC_NORMAL(njj)*XRC_NORMAL(njj)
          ENDDO
C         Include coupling coefficient CG(3).
          SUM=SUM*CG(3)**2
        ELSE !Orthotropic monodomain
C         Calculate normal in fibre coords:
C         XN_NORMAL = XRC_NORMAL * DXRCXN
          DO mjj=1,NJTR
            SUM=0.0d0
            DO njj=1,NJTR
              SUM=SUM+XRC_NORMAL(njj)*XRC_FIBRE(njj,mjj)
            ENDDO
            XN_NORMAL(mjj)=SUM
          ENDDO
C         CG(3..) are space constants.
          SUM=0.0d0
          DO mjj=1,NJTR
            TEMP=CG(2+mjj)*XN_NORMAL(mjj)
            SUM=SUM+TEMP*TEMP
          ENDDO
C         Include c0 and Gauss weight.
C         (FACTOR, SE and upwind constants are included in ZPRP_FACE.)
          niyg=niyg+1
          IF(niyg.LE.NIYGFM) YGF(niyg)=CG(2)*DSQRT(SUM)*WG
        ENDIF !iso/aniso-tropic

      ENDIF !Flux / Derivative discontinuity

      IF(niyg.GT.NIYGFM) THEN
        IEND=0
        CALL APPENDC(IEND,' >>Increase NIYGFM to ',ERROR)
        CALL APPENDI(IEND,niyg,ERROR)
        GOTO 9999
      ENDIF
      CALL EXITS('XGYGF')
      RETURN
 9999 CALL ERRORS('XGYGF',ERROR)
      CALL EXITS('XGYGF')
      RETURN 1
      END


