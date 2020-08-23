      SUBROUTINE FACE_INT_LGE(IBT,IDO,IDOXFT,INP,JDOXFT,MYMS_F,
     '  NBH,NBHF,NEA,NEAXT,NF_EA,NHF,NIEF,NKB,NKHE,NORMALSIGN1,NNB,NNF,
     '  NPNE,nr,NSB,NSFE,NVHE,NYNP,NYNS_E,nx,SE,SM_F,SN_E,ERROR,*)

C#### Subroutine: FACE_INT_LGE
C###  Description:
C###    FACE_INT_LGE calculates solution and residual multipliers for
C###    FEM integrals on one face when changing between global and face
C###    values.
C***  Assumes nc=1.
C***  (Upwind terms will probably only be used for this matrix.)
C#### Variable: IDOXFT(nhx)
C###  Type: INTEGER
C###  Set_up: FACE_INT_LGE
C###  Description:
C###    IDOXFT(nhx) is the number of domain weighting functions for each
C###    dof in the face (corresponding to the number of derivatives
C###    across the face) for dependent variable nhx.
C#### Variable: idoxf
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    idoxf is the loop variable for out-of-face derivative order when
C###    using weighting functions in a surface of dimension
C###    NJ_LOC(NJL_GEOM,0,nr)-1.
C#### Variable: JDOXFT
C###  Type: INTEGER
C###  Set_up: FACE_INT_LGE
C###  Description:
C###    JDOXFT is the number of different stiffness matrices for
C###    different out-of-face derivative orders of the interpolation
C###    functions.
C#### Variable: jdoxf
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    jdoxf is the loop variable for out-of-face derivative order when
C###    evaluating interpolation functions in a surface of dimension
C###    NJ_LOC(NJL_GEOM,0,nr)-1.
C#### Variable: NYNS_E(0:ns_e,1:JDOXFT,1:NEAXT,nhx)
C###  Type: INTEGER
C###  Set_up: FACE_INT_LGE
C###  Description:
C###    NYNS_E(0,jdoxf,neax,nhx) is the total number of mesh dofs from
C###    the adjacent element neax contributing to the face integral
C###    terms for cross face derivative jdoxf of dependent variable
C###    nhx.  NYNS_E(ns_e,jdoxf,neax,nhx) are the global positions in
C###    YP(ny,1) of these dofs.
C###    (jdoxf = 1 for discontinuity terms.)
C#### Variable: NSFE(ns_e,1:JDOXFT,neax,nhx)
C###  Type: INTEGER
C###  Set_up: FACE_INT_LGE
C###  Description:
C###    NSFE(ns_e,jdoxf,nhx) are the corresponding face dofs for the
C###    elements of SN_E.  (jdoxf = 1 for discontinuity terms.)
C#### Variable: SN_E(ns_e,1:IDOXFT/JDOXFT,1:NEAXT,nhx)
C###  Type: REAL*8
C###  Set_up: FACE_INT_LGE
C###  Description:
C###    SN_E(ns_e,idoxf/jdoxf,neax,nhx) are the factors to apply to adjacent
C###    element dependent variable nhx dof ns_e to obtain the
C###    appropriate face nodal derivative discontinuities for cross face
C###    deriv idoxf/jdoxf of the weighting function (discontinuity
C###    terms) / dependent variable (flux terms) .
C#### Variable: MYMS_F(0:ns_f,1:IDOXFT,1:NEAXT,nhx)
C###  Type: INTEGER
C###  Set_up: FACE_INT_LGE
C###  Description:
C###    MYMS_F(0,idoxf,neax,nhx) is the total number of face residuals
C###    for the adjacent element neax and weighting function
C###    cross face derivative idoxf.  MYMS_E(ns_f,idoxf,neax,nhx)
C###    are the global positions in YP(ny,4) for these residuals.
C#### Variable: NORMALSIGN1
C###  Type: INTEGER
C###  Set_up: FACE_INT_LGE
C###  Description:
C###    For flux integration, the normal always crosses the face in the
C###    direction of increasing xi.  This variable contains a
C###    multiplier for the normal to change it to an outward facing
C###    normal from the first element adjacent to the face.
C#### Variable: SM_F(ns_f,1:IDOXFT,neax,nhx)
C###  Type: REAL*8
C###  Set_up: FACE_INT_LGE
C###  Description:
C###    SM_F(ns_f,idoxf,neax,nhx) are the factors to apply to face
C###    residual values ns_f of residual equation nhx when assembling
C###    into the global residual vector.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),IDOXFT(NHM),
     '  JDOXFT,INP(NNM,NIM,NBFM),MYMS_F(0:NSFM,2,2,NHM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM),NEA(2),NEAXT,NF_EA(2),
     '  NHF,NIEF(0:2),
     '  NKB(2,2,2,NNM,NBFM),NKHE(NKM,NNM,NHM,NEM),NORMALSIGN1,
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NPNE(NNM,NBFM,NEM),nr,
     '  NSB(NKM,NNM,NBFM),NSFE(NSM,2,2,NHM),
     '  NVHE(NNM,NBFM,NHM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM),NYNS_E(0:NSM,2,2,NHM),nx
      REAL*8 SE(NSM,NBFM,NEM),SM_F(NSFM,2,2,NHM),SN_E(NSM,2,2,NHM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER IDOI(3),IDOX,idoxf,INPX,INPXF,jdoxf,nb_e,nb_f,ne,neax,
     '  nf_e,nh,nhx,ni,nk,nk_e,nk_f,NORMALSIGN(3),nn_f,nn_e,np,
     '  ns_e,ns_f,ns_i(2),nv
      REAL*8 PH3D2(2,2,2),PH3D3(2,2,2),PL2D1(3,3),UPWHD(2)
      LOGICAL COLLAPSED,FLUX,HERMITE,PETROV,QUADLAG
C     For flux integration, the normal always crosses the face in the
C     direction of increasing xi.  This is compensated for by
C     NORMALSIGN1 to provide an integral with an outward facing normal.
C     NORMALSIGN depends on the cross face node position index INPXF.
C     This is applied to the weighting functions as the direction of
C     derivatives of the dependent variable may be needed in the
C     calculation of weighting functions.  It is assumed that if INPXF
C     is 1, xi=0; otherwise xi=1.
      DATA NORMALSIGN/-1,1,1/
C     Quadratic Lagrange basis function first derivatives PL2D1(J,K)
C     corresponding to node J, and at the position of node K.
C     Values at node 2 should not be used.  Signs are not modified.
      DATA PL2D1(1,3),PL2D1(2,3),PL2D1(3,3)/ 1d0,-4d0, 3d0/
      DATA PL2D1(1,1),PL2D1(2,1),PL2D1(3,1)/-3d0, 4d0,-1d0/
C     Hermite basis function derivatives PH3D*(I,J,K) corresponding to
C     derivative I of node J, and at the position of node K.
C     Signs are modified for direct insertion into upwind terms.
      DATA PH3D2/-6.0d0,-4.0d0,6.0d0,-2.0d0,-6.0d0,-2.0d0,6.0d0,-4.0d0/
      DATA PH3D3/12.0d0,6.0d0,-12.0d0,6.0d0,-12.d0,-6.0d0,12.0d0,-6.0d0/
C     Hermite upwind multipliers
      DATA UPWHD/13.1652d-3,-6.4508d-3/

      CALL ENTERS('FACE_INT_LGE',*9999)

      PETROV=ITYP15(nr,nx).NE.3

      IF(PETROV) THEN !Derivative Discontinuity Term
        FLUX=.TRUE.
        JDOXFT=2 !stiffness matrix depends on out-of-face deriv order
C       If this is an interior face flux term then fluxes on either side
C       of the face will cancel unless derivatives are discontinuous.
      ELSE !ITYP15(nr,nx)=1or2 !Flux Term
C!!!    should check on order of derivatives
        FLUX=.TRUE. !first derivative discontinuities
        IF(FLUX) THEN
          JDOXFT=2 !stiffness matrix depends on out-of-face deriv order
        ELSE
          JDOXFT=1 !stiffness matrix does not depend on deriv order
        ENDIF
      ENDIF

C***  Determine element nodal out-of-face derivative information.
C***  Derivatives in the face are interpolated from facial nodal values.

C     For upwinding terms, facial nodal values are calculated from
C     interpolation of adjacent element nodal values.

      DO nhx=1,NHF
        nh=NH_LOC(nhx,nx)
        nb_f=NBHF(nh)
        CALL ASSERT(NST(nb_f).LE.NSFM,'>>NSFM too small',ERROR,*9999)
C       If this is an interior face flux term then it will not be
C       evaluated unless derivatives are discontinuous.  The weights
C       must be continous, and so must not depend on derivatives.
        IF(PETROV.AND.NEAXT.EQ.1) THEN !petrov,exterior
          IDOXFT(nhx)=2 !values and derivs of weighting functions used
        ELSE IF(FLUX) THEN
          IDOXFT(nhx)=1 !Only Galerkin weights
        ELSE
          IDOXFT(nhx)=2 !values and derivs of weighting functions used
        ENDIF
        DO neax=1,NEAXT !loop over adjacent elements
          ne=NEA(neax)
          nb_e=NBH(nh,1,ne)
          nf_e=NF_EA(neax)
          ni=NIEF(0)
C         Defaults for interpolation
          QUADLAG=.FALSE. !quadratic Lagrange?
          HERMITE=.FALSE. !cubic or quadratic Hermite?
          COLLAPSED=.FALSE. !collapsed?
          nn_f=1 !any node
C         Determine type of interpolation in out-of-face direction.
          IF(IBT(1,ni,nb_e).EQ.1.AND.IBT(2,ni,nb_e).EQ.2) THEN
            QUADLAG=.TRUE.
          ELSE IF(IBT(1,ni,nb_e).EQ.2) THEN
            HERMITE=.TRUE.
          ELSE IF(IBT(1,ni,nb_e).EQ.5.AND.IBT(2,ni,nb_e).EQ.4) THEN
            HERMITE=.TRUE.
            COLLAPSED=.TRUE. !collapsed at xi=0
            nn_f=NNT(nb_f) !a non-collapsed node
          ELSE IF(IBT(1,ni,nb_e).EQ.6.AND.IBT(2,ni,nb_e).EQ.4) THEN
            HERMITE=.TRUE.
            COLLAPSED=.TRUE. !collapsed at xi=1
            !nn_f=1 !a non-collapsed node
          ENDIF
C         Find node index in cross face dirn for face nodes:
          INPXF=INP(NNF(1+nn_f,nf_e,nb_e),ni,nb_e)
          IF(neax.EQ.1) NORMALSIGN1=NORMALSIGN(INPXF)
          IF(FLUX.AND.HERMITE) THEN !Petrov Flux Term
            ns_f=0
            DO idoxf=1,IDOXFT(nhx) !cross face deriv order
              ns_i(idoxf)=0
            ENDDO
            DO nn_f=1,NNT(nb_f)
              nn_e=NNF(1+nn_f,nf_e,nb_e)
              np=NPNE(nn_e,nb_e,ne)
              nv=NVHE(nn_e,nb_e,nh,ne)
              DO nk_f=1,NKT(nn_f,nb_f)
                ns_f=ns_f+1
                IDOI(NIEF(1))=IDO(nk_f,nn_f,1,nb_f)
                IDOI(NIEF(2))=IDO(nk_f,nn_f,2,nb_f)
                DO jdoxf=1,JDOXFT !cross face deriv order
                  IDOI(NIEF(0))=jdoxf
                  nk_e=NKB(IDOI(1),IDOI(2),IDOI(3),nn_e,nb_e)
C                 Face element ns number.
                  NSFE(ns_f,jdoxf,neax,nhx)=ns_f
                  IF(nk_e.NE.0) THEN
                    nk=NKHE(nk_e,nn_e,nh,ne)
                    ns_e=NSB(nk_e,nn_e,nb_e)
                    IF(jdoxf.LE.IDOXFT(nhx)) THEN
                      idoxf=jdoxf
                      ns_i(idoxf)=ns_i(idoxf)+1
C                     Find position my in global array YP.
                      MYMS_F(ns_i(idoxf),idoxf,neax,nhx)=
     '                  NYNP(nk,nv,nh,np,1)
C                     Store scale factor corresponding to residual YP(4,my).
C                     (same sign in each adjacent element)
C                      IF(neax.EQ.1) THEN
                      SM_F(ns_i(idoxf),idoxf,neax,nhx)=
     '                  SE(ns_e,nb_e,ne)
C                      ELSE
C                        SM_F(ns_i(ijdoxf),ijdoxf,neax,nhx)=
C     '                    SE(ns_e,nb_e,ne)*NORMALSIGN1
C                      ENDIF
                    ENDIF
C                   Find position ny in global array YP.
                    NYNS_E(ns_f,jdoxf,neax,nhx)=
     '                NYNP(nk,nv,nh,np,2)
C                    IF(neax.EQ.1) THEN
C                   Calculate multipliers of YP(ny,1) for face integrals.
                    SN_E(ns_f,jdoxf,neax,nhx)=SE(ns_e,nb_e,ne)
C                    ELSE
CCC                     Calculate multipliers of YP(ny,1) for face
CCC                     integrals.  This means that ZPZF will provide a
CCC                     difference of values on either side of the face.
CCC                     Check if same dofs are used in each element
CC                      IF(NYNS_E(ns_f,ijdoxf,1,nhx)
CC     '                  .NE.NYNS_E(ns_f,ijdoxf,2,nhx).OR
CC     '                  .SN_E(ns_f,ijdoxf,1,nhx)
CC     '                  .NE.SE(ns_e,nb_e,ne)) THEN
CC                        INTEGRATE=.TRUE.
CC                      ENDIF
C                      SN_E(ns_f,ijdoxf,neax,nhx)=
C     '                  SE(ns_e,nb_e,ne)
CCC                      IF(SN_E(ns_f,ijdoxf,nhx).NE.SE(ns_e,nb_e,ne))
CCTHEN C                        ERROR='>>Differing scale factors on face'
CCC                        GO TO 9999
CCC                      ENDIF
C                    ENDIF
                  ELSE !nk_e!=0
                    NYNS_E(ns_f,jdoxf,neax,nhx)=0
                  ENDIF !nk_e
                ENDDO !jdoxf
              ENDDO !nk_f
            ENDDO !nn_f
            DO jdoxf=1,JDOXFT !cross face deriv order
              IF(jdoxf.LE.IDOXFT(nhx)) THEN
                idoxf=jdoxf
C               Number of residuals.
                MYMS_F(0,idoxf,neax,nhx)=ns_i(idoxf)
              ENDIF
C             Number of contributing adjacent element dofs
              NYNS_E(0,jdoxf,neax,nhx)=ns_f
            ENDDO !jdoxf

C          ELSE IF(FLUX.AND.HERMITE) THEN !Flux Term
C            ns_f=0
C            DO ijdoxf=1,IDOXFT(nhx) !cross face deriv order
C              ns_i(ijdoxf)=0
C            ENDDO
C            DO nn_f=1,NNT(nb_f)
C              nn_e=NNF(1+nn_f,nf_e,nb_e)
C              np=NPNE(nn_e,nb_e,ne)
C              nv=NVHE(nn_e,nb_e,nh,ne)
C              DO nk_f=1,NKT(nn_f,nb_f)
C                ns_f=ns_f+1
C                IDOI(NIEF(1))=IDO(nk_f,nn_f,1,nb_f)
C                IDOI(NIEF(2))=IDO(nk_f,nn_f,2,nb_f)
C                DO ijdoxf=1,JDOXFT !cross face deriv order
C                  IDOI(NIEF(0))=ijdoxf
C                  nk_e=NKB(IDOI(1),IDOI(2),IDOI(3),nn_e,nb_e)
CC                 Face element ns number.
C                  NSFE(ns_f,ijdoxf,neax,nhx)=ns_f
C                  IF(nk_e.NE.0) THEN
C                    nk=NKE(nk_e,nn_e,nb_e,ne)
C                    ns_e=NSB(nk_e,nn_e,nb_e)
C                    IF(ijdoxf.EQ.2) THEN !derivative
C                      ns_i(ijdoxf)=ns_i(ijdoxf)+1
CC                     Find position my in global array YP.
C                      MYMS_F(ns_i(ijdoxf),neax,1,nhx)=
C     '                  NYNP(nk,nv,nh,np,1)
C                      MYMS_F(ns_i(ijdoxf),neax,2,nhx)=
C     '                  NYNP(nk,nv,nh,np,1)
CC                     Store scale factor corresponding to residual YP(4,my).
CC                     (same sign in each adjacent element)
CC                      IF(neax.EQ.1) THEN
C                      SM_F(ns_i(ijdoxf),neax,1,nhx)=
C     '                  SE(ns_e,nb_e,ne)
C                      SM_F(ns_i(ijdoxf),neax,2,nhx)=
C     '                  SE(ns_e,nb_e,ne)
CC                      ELSE
CC                        SM_F(ns_i(ijdoxf),ijdoxf,neax,nhx)=
CC     '                    SE(ns_e,nb_e,ne)*NORMALSIGN1
CC                      ENDIF
C                    ENDIF
CC                   Find position ny in global array YP.
C                    NYNS_E(ns_f,ijdoxf,neax,nhx)=
C     '                NYNP(nk,nv,nh,np,2)
CC                    IF(neax.EQ.1) THEN
CC                   Calculate multipliers of YP(ny,1) for face integrals.
C                    SN_E(ns_f,ijdoxf,neax,nhx)=SE(ns_e,nb_e,ne)
CC                    ELSE
CCCC                     Calculate multipliers of YP(ny,1) for face
CCCC                     integrals.  This means that ZPZF will provide a
CCCC                     difference of values on either side of the face.
CCCC                     Check if same dofs are used in each element
CCC                      IF(NYNS_E(ns_f,ijdoxf,1,nhx)
CCC     '                  .NE.NYNS_E(ns_f,ijdoxf,2,nhx).OR
CCC     '                  .SN_E(ns_f,ijdoxf,1,nhx)
CCC     '                  .NE.SE(ns_e,nb_e,ne)) THEN
CCC                        INTEGRATE=.TRUE.
CCC                      ENDIF
CC                      SN_E(ns_f,ijdoxf,neax,nhx)=
CC     '                  SE(ns_e,nb_e,ne)
CCCC                      IF(SN_E(ns_f,ijdoxf,nhx).NE.SE(ns_e,nb_e,ne))
CCCTHEN C                        ERROR='>>Differing scale factors on face'
CCCC                        GO TO 9999
CCCC                      ENDIF
CC                    ENDIF
C                  ELSE !nk_e!=0
C                    NYNS_E(ns_f,ijdoxf,neax,nhx)=0
C                  ENDIF !nk_e
C                ENDDO !ijdoxf
C              ENDDO !nk_f
C            ENDDO !nn_f
C            DO ijdoxf=1,JDOXFT !cross face deriv order
C              IF(ijdoxf.EQ.2) THEN
CC               Number of residuals.
C                MYMS_F(0,neax,1,nhx)=ns_i(ijdoxf)
C                MYMS_F(0,neax,2,nhx)=ns_i(ijdoxf)
C              ENDIF
CC             Number of contributing adjacent element dofs
C              NYNS_E(0,ijdoxf,neax,nhx)=ns_f
C            ENDDO !ijdoxft
C
          ELSE IF(ITYP15(nr,nx).EQ.3.AND.QUADLAG)
     '        THEN !Derivative discontinuity term for quadratics
            ns_e=0
            DO nn_e=1,NNT(nb_e)
              nn_f=NNB(
     '          INP(nn_e,NIEF(1),nb_e),INP(nn_e,NIEF(2),nb_e),1,nb_f)
              INPX=INP(nn_e,NIEF(0),nb_e) !node index in cross face dirn
              np=NPNE(nn_e,nb_e,ne)
              nv=NVHE(nn_e,nb_e,nh,ne)
              DO nk_e=1,NKT(nn_e,nb_e)
                ns_e=ns_e+1
C               Find position ny in global array YP
                nk=NKHE(nk_e,nn_e,nh,ne)
                NYNS_E(ns_e,2,neax,nhx)=NYNP(nk,nv,nh,np,2)
C               Find corresponding face ns number
                nk_f=NKB(IDO(nk_e,nn_e,NIEF(1),nb_e),
     '            IDO(nk_e,nn_e,NIEF(2),nb_e),1,nn_f,nb_f)
                ns_f=NSB(nk_f,nn_f,nb_f)
                NSFE(ns_e,2,neax,nhx)=ns_f
C               Calculate multipliers of YP(ny,1)
                SN_E(ns_e,2,neax,nhx)=
     '            PL2D1(INPX,INPXF)*SE(ns_e,nb_e,ne)
                IF(INPX.EQ.INPXF) THEN !node in face
C                 Find position ny in global array YP
                  NYNS_E(ns_f,1,neax,nhx)=NYNP(nk,nv,nh,np,2)
C                 Find corresponding face ns number
                  NSFE(ns_f,1,neax,nhx)=ns_f
C                 Calculate multipliers of YP(ny,1)
                  SN_E(ns_f,1,neax,nhx)=SE(ns_e,nb_e,ne)
C                 Find position my in global array YP
                  MYMS_F(ns_f,1,neax,nhx)=NYNP(nk,nv,nh,np,1)
C                 Store scale factor corresponding to residual YP(4,my)
                  SM_F(ns_f,1,neax,nhx)=SE(ns_e,nb_e,ne)
                  IF(neax.EQ.2) THEN !check values are consistent
                    IF(MYMS_F(ns_f,1,1,nhx).NE.
     '                MYMS_F(ns_f,1,2,nhx)) THEN
                      ERROR='>>Differing degrees of freedom on face'
                      GO TO 9999
                    ENDIF
                    IF(SM_F(ns_f,1,1,nhx).NE.
     '                SM_F(ns_f,1,2,nhx)) THEN
                      ERROR='>>Differing scale factors on face'
                      GO TO 9999
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO !nk_e
            ENDDO !nn_e
C           Number of residuals.
            MYMS_F(0,1,neax,nhx)=NST(nb_f)
C           Number of contributing adjacent element dofs
            NYNS_E(0,1,neax,nhx)=NST(nb_f) !values depend on face only
            NYNS_E(0,2,neax,nhx)=ns_e !derivs depend on elements
C           Number of contributing adjacent element dofs

          ELSE IF(ITYP15(nr,nx).EQ.3.AND.HERMITE.AND..NOT.COLLAPSED)
     '        THEN !Derivative Discontinuity Term
            ns_e=0
            DO nn_e=1,NNT(nb_e)
              nn_f=NNB(
     '          INP(nn_e,NIEF(1),nb_e),INP(nn_e,NIEF(2),nb_e),1,nb_f)
              INPX=INP(nn_e,NIEF(0),nb_e) !node index in cross face dirn
              np=NPNE(nn_e,nb_e,ne)
              nv=NVHE(nn_e,nb_e,nh,ne)
              DO nk_e=1,NKT(nn_e,nb_e)
                ns_e=ns_e+1
C               Find position ny in global array YP
                nk=NKHE(nk_e,nn_e,nh,ne)
                NYNS_E(ns_e,1,neax,nhx)=NYNP(nk,nv,nh,np,2)
C               Find corresponding face ns number
                nk_f=NKB(IDO(nk_e,nn_e,NIEF(1),nb_e),
     '            IDO(nk_e,nn_e,NIEF(2),nb_e),1,nn_f,nb_f)
                ns_f=NSB(nk_f,nn_f,nb_f)
                NSFE(ns_e,1,neax,nhx)=ns_f
C               Calculate multipliers of YP(ny,1) for upwind terms
                IDOX=IDO(nk_e,nn_e,NIEF(0),nb_e) !der ind in cross face dirn
                SN_E(ns_e,1,neax,nhx)=
     '            PH3D3(IDOX,INPX,INPXF)*SE(ns_e,nb_e,ne)
                SN_E(ns_e,2,neax,nhx)=
     '            PH3D2(IDOX,INPX,INPXF)*SE(ns_e,nb_e,ne)
                IF(INPX.EQ.INPXF) THEN !node in face
C                 Find position my in global array YP
                  MYMS_F(ns_f,IDOX,neax,nhx)=NYNP(nk,nv,nh,np,1)
C                 Store scale factor corresponding to residual YP(4,my)
                  SM_F(ns_f,IDOX,neax,nhx)=
     '              SE(ns_e,nb_e,ne)*UPWHD(IDOX)*FACTOR
                  IF(neax.EQ.2) THEN !check values are consistent
                    IF(MYMS_F(ns_f,IDOX,1,nhx).NE.
     '                MYMS_F(ns_f,IDOX,2,nhx)) THEN
                      ERROR='>>Differing degrees of freedom on face'
                      GO TO 9999
                    ENDIF
                    IF(SM_F(ns_f,IDOX,1,nhx).NE.
     '                SM_F(ns_f,IDOX,2,nhx)) THEN
                      ERROR='>>Differing scale factors on face'
                      GO TO 9999
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO !nk_e
            ENDDO !nn_e
            DO idoxf=1,IDOXFT(nhx) !cross face deriv order
C             Number of residuals.
              MYMS_F(0,idoxf,neax,nhx)=NST(nb_f)
            ENDDO
C           Number of contributing adjacent element dofs
            NYNS_E(0,1,neax,nhx)=ns_e
          ELSE !Not Hermite
            ERROR='Basis function type not implemented'
            GO TO 9999
          ENDIF !Hermite
        ENDDO !neax
      ENDDO !nhx

      CALL EXITS('FACE_INT_LGE')
      RETURN
 9999 CALL ERRORS('FACE_INT_LGE',ERROR)
      CALL EXITS('FACE_INT_LGE')
      RETURN 1
      END


