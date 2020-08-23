      SUBROUTINE CALC_FV_GRID_COEF(NITB,nq,NQGP,NQGP_PIVOT,
     '  nr,nx_ext,NXQ,COEFFSEXT,CQ,GM,NQGW,XQ,PROPQ,BIDOMAIN,
     '  FIXQ,STATIONARY,ERROR,*)

C#### Subroutine: CALC_FV_GRID_COEF
C###  Description:
C###    CALC_FV_GRID_COEF is used for the grid-based Finite Volume
C###    scheme. It finds the coeffs of the non-zeros in the stiffenss
C###    matrix (NQGW) and the diagonal non-zero in the global mass 
C###    matrix (GM) for a single grid point, nq. The number of non-zeros
C###    on each row of the stiffness matrix is 3 in 1D, 5 in 2D and 7 in
C###    3D. 
C###    Important: The method is currently only implemented for grid points
C###    that are orthogonal in global xyz space. Code to test for this may
C###    be commented out to gain significant performance improvement. 
C###
C###    PROPQ is used to store geometric volume information and also to 
C###    store grid point derivatives and accumulate the secondary flux
C###    source values. The following arrangements are used in PROPQ(3,3,4,2)
C###    (these are expressed in terms of bidomain quantities).
C###    
C###    U   : U geometric volume face vector (i.e. vectors in the face)
C###    V   : V geometric volume face vector (i.e. vectors in the face)
C###    INU : product of intracellular conductivity tensor with face normal
C###          and U face vector and face area
C###    INV : product of intracellular conductivity tensor with face normal
C###          and V face vector and face area
C###    ENU : product of extracellular conductivity tensor with face normal
C###          and U face vector and face area
C###    ENV : product of extracellular conductivity tensor with face normal
C###          and V face vector and face area
C###    DV  : derivative approximations of transmembrane potential
C###    DE  : derivative approximations of extracellular potential
C###    LP  : finite volume dimensions in +dir
C###    LN  : finite volume dimensions in -dir
C###    SF*1: secondary flux contribution to first equation rhs 
C###          contributions are stored for the three possible positive directions
C###          (SFP*1) and the three possible negative directions (SFN*1).
C###    SF*2: secondary flux contribution to second equation rhs
C###          contributions are stored for the three possible positive directions
C###          (SFP*1) and the three possible negative directions (SFN*1).
C###
C###    PROPQ(i,j,k,l):
C###  
C###    l=1
C###    k=1           k=2          k=3    k=4
C###    +1            +2           +3             - face
C###    x   y   z     x   y   z    x y z
C###    U,  U,  U     U,  U,  U    U,U,U  -,-,-
C###    V,  V,  V     V,  V,  V    V,V,V  -,-,-
C###    INU,INU,INU   INV,INV,INV  -,-,-  -,-,-
C###    +1  +2  +3    +1  +2  +3                  - face
C###
C###    l=2
C###    k=1           k=2          k=3      k=4
C###    x   y   z     +1  +2  +3                  - direction/face
C###    DV, DV, DV    LP, LP, LP   SFP11,SFP21,SFP31  SFN11,SFN21,SFN31
C###    DE, DE, DE    LN, LN, LN   SFP12,SFP22,SFP32  SFN12,SFN22,SFN32
C###    ENU,ENU,ENU   ENV,ENV,ENV  -,  -,-  -,-,-
C###    +1  +2  +3    +1  +2  +3                  - face
C###    
C***  Created by Mark Trew 10 December 2002

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc' !NJL_FIBR
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER NITB,nq,NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),nr,
     '  nx_ext,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 COEFFSEXT(NQGM),CQ(NMM,NQM),GM(NZ_GM_M),NQGW(NQGM),
     '  XQ(NJM,NQM),PROPQ(3,3,4,2,NQM)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,FIXQ(NYQM,NIYFIXM,NXM),STATIONARY

!     Local variables
      INTEGER DOFF(2),D2NQGP(0:NQGM),lower,upper,PLoc,FDIR(NQGM),
     '  DIRF(-3:3),i,j,k,l,dir,nj1,nj,adj,nidx
      REAL*8 SUMSQ,
     '  FDIST(0:NQGM),LDIR(0:3),VOLUME,FAREA(-3:3),CMAM,
     '  SUMEXT,SUMINT,FNORM(3,NQGM),FANGLE(NJM,NQGM),
     '  DXDXI(3,3),DXDNU(3,3),DNUDX(3,3),DET,ETERM,ITERM,
     '  TSIGE(3,3),TSIGI(3,3),DOT_PROD,NV(3),UV(3),VV(3),LENGTH
      LOGICAL MONODOMAIN,TIMEVARY 
      CHARACTER ASSERTSTR*120

      CALL ENTERS('CALC_FV_GRID_COEF',*9999)

C logical problem characteristics
      TIMEVARY   = .FALSE.
      MONODOMAIN = .FALSE.
      IF (.NOT.STATIONARY) TIMEVARY   = .TRUE.
      IF (.NOT.BIDOMAIN)   MONODOMAIN = .TRUE.

C assertion:
C check that the number of Xi directions for the element is equal to
C     the number of global coordinates
      CALL ASSERT(NITB.EQ.NJ_LOC(NJL_GEOM,0,nr),
     '  '>>#Xi directions for element must equal #global'
     '  //' coordinates (for grid-based FE method)',ERROR,*9999)

C end assertion

C set material parameter (conductivities for generalise laplace, 
C membrane capacitance, cell membrane area to volume ratio and 
C conductivities for bidomain and monodomain)
      IF (MONODOMAIN) THEN
        IF (STATIONARY) THEN
          DOFF(1) = 0 ! no material parameters on time derivative
                      ! at most three conductivities
        ELSE
          DOFF(1) = 2 ! material parameters in mass matrix
                      ! at most three conductivities
        ENDIF
      ELSE  ! bidomain
        IF (STATIONARY) THEN
          DOFF(2) = 0 ! no material parameters on time derivative
          DOFF(1) = 3 ! three conductivities in extracellular dmn
                      ! three conductivities in intracellular dmn
        ELSE
          DOFF(2) = 2 ! material parameters in mass matrix
          DOFF(1) = 5 ! three conductivities in extracellular dmn
                      ! three conductivities in intracellular dm
        ENDIF
      ENDIF

C initialise NQGP (grid point support vector), NQGW (non-zero entries in row 
C of intracellular stiffness matrix), GM (non-zero entries in global mass 
C matrix), COEFFSEXT (non-zero entries in row of extracellular stiffness
C matrix), FNORM (face normal vector), NV, UV and VV (unit in face vectors),
C and PROPQ (grid point properties array)     
      NQGP(0,nq)=0
      D2NQGP(0) = 0 ! an NQGP for domain 2 - for bidomain, the 
                    ! intracellular domain 
      FDIST(0) = 0.0d0 ! distance to adjacent face 0 is 0 (boundary)
      DO i=1,NQGM ! initialise variables
        NQGP(i,nq)=0
        D2NQGP(i)=0
        COEFFSEXT(i)=0.0d0
        NQGW(i)=0.0d0
        GM(i+NQGM*(nq-1))=0.0d0    
        FNORM(1,i)=0.0d0
        FNORM(2,i)=0.0d0
        FNORM(3,i)=0.0d0
      ENDDO !i 

      DO i=1,3
        NV(i) = 0.0d0
        UV(i) = 0.0d0
        VV(i) = 0.0d0
      ENDDO
 
      DO i=1,3
        DO j=1,3
          DO k=1,4
            DO l=1,2
              PROPQ(i,j,k,l,nq) = 0.0d0
            ENDDO ! l
          ENDDO ! k
        ENDDO ! j
      ENDDO ! i

C grid points are stored in NQGP in order: -xi3,-xi2,-xi1,nq,xi1,xi2,xi3
      lower = -1
      upper = 1
      IF(NITB.GE.2) THEN
        lower = -2
        upper = 2
        IF(NITB.GE.3) THEN
          lower = -3
          upper = 3
        ENDIF
      ENDIF

C initialise the dimension in each coordinate direction array
      DO i=0,upper
        LDIR(i) = 0.0d0
      ENDDO

C pre-processing to look for face flux breaks and to determine 
C support vectors, face normal vectors, distances to adjacent grid
C points, the dimensions in a given direction and mappings from 
C face index to direction index and vice versa.
      DO dir = lower,upper
        DIRF(dir) = 0
        IF (dir.ne.0) THEN ! adjacent grid points
          IF (NXQ(dir,0,nq).GT.0) THEN ! adjacent grid pts exist
            DO adj = 1,NXQ(dir,0,nq) ! loop over all gps in specified direction
              NQGP(0,nq) = NQGP(0,nq)+1 ! increment the number of adjacent points
              NQGP(NQGP(0,nq),nq) = ABS(NXQ(dir,adj,nq)) ! store adjacent point
              IF (NXQ(dir,adj,nq).GT.0) THEN ! adjacent intracellular gps exist
                D2NQGP(0) = NQGP(0,nq) ! increment number of adj intracell. gps
                D2NQGP(NQGP(0,nq)) = NXQ(dir,adj,nq) ! store adj intracell. gp
              ENDIF

              SUMSQ = 0.0d0
              DO nj1=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,nj1,nr)
                SUMSQ = SUMSQ+
     '            (XQ(nj,NQGP(NQGP(0,nq),nq)) - XQ(nj,nq))**2
                FNORM(nj,NQGP(0,nq)) = XQ(nj,NQGP(NQGP(0,nq),nq)) - 
     '                               XQ(nj,nq)
              ENDDO
              FDIST(NQGP(0,nq)) = sqrt(SUMSQ)
              LDIR(ABS(dir)) = LDIR(ABS(dir))+0.5d0*FDIST(NQGP(0,nq))
              FDIR(NQGP(0,nq)) = dir ! map from face index to direction number
              DIRF(dir) = NQGP(0,nq) ! map from direction number to face index

C calculate face fibre/imbrication/sheet angle wrt xi1 using simple linear 
C interpolation
              DO nj1=1,NJ_LOC(NJL_FIBR,0,nr)
                nj=NJ_LOC(NJL_FIBR,nj1,nr)
                FANGLE(nj,NQGP(0,nq)) = 
     '            0.5d0*(XQ(nj,NQGP(NQGP(0,nq),nq))+XQ(nj,nq))
              ENDDO !nj1

            ENDDO !adj
          ENDIF
        ELSE ! grid point of interest
          NQGP(0,nq) = NQGP(0,nq)+1
          NQGP(NQGP(0,nq),nq) = nq
          PLoc = NQGP(0,nq)
          D2NQGP(0) = NQGP(0,nq)
          D2NQGP(NQGP(0,nq)) = nq

          SUMSQ = 0.0d0
          DO nj1=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,nj1,nr)
            SUMSQ = SUMSQ+
     '        (XQ(nj,NQGP(NQGP(0,nq),nq)) - XQ(nj,nq))**2
            FNORM(nj,NQGP(0,nq)) = XQ(nj,NQGP(NQGP(0,nq),nq)) - 
     '                             XQ(nj,nq)
          ENDDO
          FDIST(NQGP(0,nq)) = sqrt(SUMSQ)
          LDIR(ABS(dir)) = LDIR(ABS(dir))+0.5d0*FDIST(NQGP(0,nq))
          FDIR(NQGP(0,nq)) = dir ! map from face index to direction number
          DIRF(dir) = NQGP(0,nq) ! map from direction number to face index 

C calculate face fibre/imbrication/sheet angle wrt xi1 using simple linear 
C interpolation
          DO nj1=1,NJ_LOC(NJL_FIBR,0,nr)
            nj=NJ_LOC(NJL_FIBR,nj1,nr)
            FANGLE(nj,NQGP(0,nq)) = 
     '        0.5d0*(XQ(nj,NQGP(NQGP(0,nq),nq))+XQ(nj,nq))
          ENDDO !nj1

        ENDIF
      ENDDO ! direction loop

C assertion
C check that the grid points are orthogonal - note: this slows down the 
C system construction significantly and can be commented out when the 
C grid points are known to be orthogonally arranged
      DO i=1,NQGP(0,nq)
        DO j=1,NQGP(0,nq)
          IF (i.NE.PLoc.AND.j.NE.PLoc.AND.i.NE.j.AND.
     '        ABS(FDIR(i)).NE.ABS(FDIR(j))) THEN
            WRITE(ASSERTSTR,
     '        '(">>The grid points are not orthogonal in xi"'
     '        //',I2,"-xi",I2)') ABS(FDIR(i)),ABS(FDIR(j))
            CALL ASSERT(
     '        DABS(DOT_PROD(FNORM(1,i),FNORM(1,j))/
     '        (FDIST(i)*FDIST(j))).LT.LOOSE_TOL,
     '        ASSERTSTR,ERROR,*9999)
          ENDIF
        ENDDO ! supporting grid point loop
      ENDDO ! supporting grid point loop
C end assertion
      
C grid fv cell volume and face areas
      VOLUME = 1.0d0
      DO i=1,upper
        VOLUME = VOLUME*LDIR(i)
      ENDDO
      DO dir = lower,upper
        FAREA(dir) = 1.0d0
        DO i=1,upper
          IF (i.ne.ABS(dir)) THEN
            FAREA(dir) = FAREA(dir)*LDIR(i)
          ENDIF
        ENDDO
      ENDDO

C store the positive face distances in PROPQ
      DO dir = 1,upper
        PROPQ(1,dir,2,2,nq) = FDIST(DIRF(dir))
        PROPQ(2,dir,2,2,nq) = FDIST(DIRF(-dir))
      ENDDO

C grid point material parameters
      CMAM = 0.0d0
      IF (TIMEVARY) CMAM =  CQ(1,nq)*CQ(2,nq)

C initialise the intra and extracellular sums
      SUMINT = 0.0d0
      SUMEXT = 0.0d0

C loop over the supporting points
      DO i=1,NQGP(0,nq)

        IF (i.NE.PLoc) THEN ! surrounding points

C determine dxdxi and dxdnu and dnudx
          DO j=1,NITB
            DO k=1,NITB
              DXDXI(j,k) = 0.0d0
            ENDDO
            DXDXI(j,j) = 1.0d0
            DXDXI(j,ABS(FDIR(i))) = 
     '        DSIGN(FNORM(j,i),DBLE(FDIR(i)))
          ENDDO      
C note: MAT_VEC returns DNUDX not DXDNU.    
          CALL MAT_VEC(NITB,nr,DNUDX(1,1),DNUDX(1,2),
     '                DNUDX(1,3),DXDXI,FANGLE(1,i),ERROR,*9999)
          CALL INVERT(NITB,DNUDX,DXDNU,DET)
C assertion:
C check that the coordinate transformation array is full rank
C between the grid point of interest and the supporting grid point 
          WRITE(ASSERTSTR,
     '      '(">>DnuDx is rank deficient between grid points: "'
     '      //',I7," and ",I7)') nq,NQGP(i,nq)
          CALL ASSERT(DABS(DET).GT.LOOSE_TOL,ASSERTSTR,ERROR,*9999)
C end assertion

C calculate conductivity tensor in global coordinates midway between
C two adjacent grid points, i.e. conductivity tensor on volume face
C n.b. for monodomain problems the two conductivity tensors are not 
C required, but there is minimal overhead in calculating the second 
C tensor so no condition is used.
          DO j=1,NITB
            DO k=1,NITB
              TSIGE(j,k) = 0.0d0
              TSIGI(j,k) = 0.0d0
              DO l=1,NITB
                TSIGE(j,k) = TSIGE(j,k)+DNUDX(j,l)*DXDNU(l,k)*
     '            0.5d0*(CQ(DOFF(1)+l,nq)+CQ(DOFF(1)+l,NQGP(i,nq)))
                TSIGI(j,k) = TSIGI(j,k)+DNUDX(j,l)*DXDNU(l,k)*
     '            0.5d0*(CQ(DOFF(2)+l,nq)+CQ(DOFF(2)+l,NQGP(i,nq)))
              ENDDO
            ENDDO
          ENDDO

C Calculate face vector set - this should work for any oriented
C normal vector
          IF ((DABS(FNORM(1,i)).GT.DABS(FNORM(2,i))).AND.
     '        (DABS(FNORM(1,i)).GT.DABS(FNORM(3,i)))) THEN
            nidx = 1            
          ELSEIF ((DABS(FNORM(2,i)).GT.DABS(FNORM(1,i))).AND.
     '            (DABS(FNORM(2,i)).GT.DABS(FNORM(3,i)))) THEN
            nidx = 2
          ELSE
            nidx = 3
          ENDIF            
          LENGTH = 0.0d0
          DO j=1,NITB
            NV(j) = FNORM(j,i)/FDIST(i)
            IF (j.NE.nidx) THEN
              UV(j) = 1.0d0
              LENGTH = LENGTH + 1.0d0
            ENDIF
          ENDDO
          UV(nidx) = 0.0d0
          DO j=1,NITB
            IF (j.NE.nidx) THEN
              UV(nidx) = UV(nidx) - NV(j)/NV(nidx)
            ENDIF            
          ENDDO
          LENGTH = LENGTH + UV(nidx)**2

          LENGTH = SQRT(LENGTH)
          DO j=1,NITB
            UV(j) = UV(j)/LENGTH
          ENDDO
          CALL CROSS(NV,UV,VV)
          IF (FDIR(i).GT.0) THEN
            DO j=1,NITB
              PROPQ(1,j,FDIR(i),1,nq) = UV(j)
              PROPQ(2,j,FDIR(i),1,nq) = VV(j)
            ENDDO
          ENDIF

C Calculate conductivity-normal vector product and inner product this
C with face vectors and mulitply by face area - do not do this if 
C there is no flux between an intracellular face. 
          IF (FDIR(i).GT.0) THEN
            DO j=1,NITB
              DO k=1,NITB
                PROPQ(3,FDIR(i),1,2,nq) = PROPQ(3,FDIR(i),1,2,nq) + 
     '            TSIGE(k,j)*NV(k)*UV(j)*FAREA(FDIR(i))
                PROPQ(3,FDIR(i),2,2,nq) = PROPQ(3,FDIR(i),2,2,nq) + 
     '            TSIGE(k,j)*NV(k)*VV(j)*FAREA(FDIR(i))
c                IF (BIDOMAIN.AND.D2NQGP(i).NE.0) THEN
                IF (D2NQGP(i).NE.0) THEN
                  PROPQ(3,FDIR(i),1,1,nq) = PROPQ(3,FDIR(i),1,1,nq) + 
     '              TSIGI(k,j)*NV(k)*UV(j)*FAREA(FDIR(i))
                  PROPQ(3,FDIR(i),2,1,nq) = PROPQ(3,FDIR(i),2,1,nq) + 
     '              TSIGI(k,j)*NV(k)*VV(j)*FAREA(FDIR(i))
                ENDIF
              ENDDO     
            ENDDO
          ENDIF     

C accumulate primary flux terms
          ETERM = 0.0d0
          ITERM = 0.0d0
          DO j=1,NITB
            DO k=1,NITB
              ETERM = ETERM + TSIGE(j,k)*FNORM(k,i)*FNORM(j,i)
c              IF (BIDOMAIN) 
              ITERM = ITERM + TSIGI(j,k)*FNORM(k,i)*FNORM(j,i)
            ENDDO
          ENDDO
          ETERM = FAREA(FDIR(i))*ETERM/(FDIST(i)**3)
          ITERM = FAREA(FDIR(i))*ITERM/(FDIST(i)**3)

C add contribution to extracellular stiffness matrix 
C if no extracellular essential boundary condition specified
          IF (BIDOMAIN.AND.
     '        (.NOT.(FIXQ(nq,1,nx_ext).OR.
     '               FIXQ(nq,3,nx_ext)))) THEN
            COEFFSEXT(i) = ETERM

C add the intracellular contribution if the link to the ith
C intracellular point exists
            IF ((D2NQGP(i).NE.0)) THEN
               COEFFSEXT(i) = COEFFSEXT(i)+ITERM
            ENDIF
            SUMEXT = SUMEXT + COEFFSEXT(i)
          ENDIF

C add intracellular contribution to the stiffness matrix if 
C the link to the ith intracellular point exists 
          IF ((D2NQGP(i).NE.0)) THEN
            NQGW(i) = ITERM
            SUMINT = SUMINT + NQGW(i)
          ENDIF

        ENDIF
      ENDDO

C place diagonal terms in the stiffness matrices
      IF (BIDOMAIN) COEFFSEXT(PLoc) = -SUMEXT
      NQGW(PLoc) = -SUMINT

C mass matrix - note the same structure is used as grid FE with
C NQGM non-zeros per row. The NQGM-1 "non-zeros" in the GM row for 
C the grid FV method are filled with zeros.
      IF (TIMEVARY) GM(PLoc+NQGM*(nq-1)) = -1.0d0*CMAM*VOLUME

C the signs on the stiffness matrices are reversed to fit the grid FE
C and collocation conventions
      DO i=1,NQGP(0,nq)
        COEFFSEXT(i) = -1.0d0*COEFFSEXT(i)
        NQGW(i) = -1.0d0*NQGW(i)
      ENDDO

C essential or analytic extracellular potential boundary conditions
      IF (BIDOMAIN.AND.
     '    (FIXQ(nq,1,nx_ext).OR.FIXQ(nq,3,nx_ext))) THEN
        DO i=1,NQGP(0,nq)
          COEFFSEXT(i) = 0.0d0
        ENDDO
        COEFFSEXT(PLoc) = 1.0d0
      ENDIF

C unsupported grid point
      IF (NQGP(0,nq).EQ.1) THEN ! unsupported grid point
        NQGP(0,nq) = 1
        NQGP(1,nq) = nq
        NQGW(1) = 1.0d0 
        IF (BIDOMAIN) COEFFSEXT(1) = 1.0d0
      ELSEIF (TIMEVARY.AND.D2NQGP(0).EQ.1) THEN
C       The intracellular point is unsupported so its value remains
C       static (initial value), i.e. dV/dt = I. An alternative is to set
C       the weights in the NQGW array to be those for Laplace's equation
C       However, this would require the identification of the unsupported
C       intracellular point when the rhs is constructed. A logical 
C       alternative would be to set the extracellular stiffness matrix
C       COEFFSEXT to zero but this matrix is used to solve for the 
C       extracellular potentials. Conclusion: to set the transmembrane
C       potential of an intracellularly isolated grid point to be dependent
C       on the surrounding grid points will either require changes in the 
C       way the data structures are used or this subroutine will need to 
C       return a parameter indicating that nq is unsupported intracellularly.
        GM(PLoc+NQGM*(nq-1)) = 1.0d0  
      ENDIF 

C sort NQGP - this is not strictly necessary for the current suite of solvers
C but may give a small performance improvement.
      CALL ISORTP(NQGP(0,nq),NQGP(1,nq),NQGP_PIVOT(1,nq))

      CALL EXITS('CALC_FV_GRID_COEF')
      RETURN
 9999 CALL ERRORS('CALC_FV_GRID_COEF',ERROR)
      CALL EXITS('CALC_FV_GRID_COEF')
      RETURN 1
      END      


