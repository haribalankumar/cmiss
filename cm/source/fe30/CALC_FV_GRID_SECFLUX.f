      SUBROUTINE CALC_FV_GRID_SECFLUX(NITB,
     '  nr,NXQ,YQ1,YQ2,PROPQ,BIDOMAIN,ERROR,*)

C#### Subroutine: CALC_FV_GRID_SECFLUX
C###  Description:
C###    CALC_FV_GRID_SECFLUX is used for the grid-based Finite Volume
C###    scheme. It calculates dependent variable first derivative 
C###    approximations at each grid point using simple finite differences.
C###    The derivatives are then used with precalculated geometric information
C###    to determine the secondary flux (as souce terms) for the grid. 
C###    A regular orthogonal grid is assumed (this is tested for in 
C###    CALC_FV_GRID_COEF).
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
C***  Created by Mark Trew 5 March 2003

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc' !NJL_FIBR
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER NITB,nr,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 YQ1(NQM),YQ2(NQM),PROPQ(3,3,4,2,NQM)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN

!     Local variables
      INTEGER nq,nq2,i,j
      REAL*8 LENGTH

      CALL ENTERS('CALC_FV_GRID_SECFLUX',*9999)

C Parallel commands
C$OMP         PARALLEL DO
CC$OMP&          DEFAULT(none)
C$OMP&          PRIVATE(i,j,LENGTH,nq,nq2)
C$OMP&          SHARED(BIDOMAIN,NITB,NQR,nr,NXQ,PROPQ,YQ1,YQ2)

C Loop over all grid points
      DO nq = NQR(1,nr),NQR(2,nr)

C Loop over coordinate directions
        DO i = 1,NITB

C Calculate derivative approximations in direction dir - simple 1D differences
C are used since the grid is orthogonal. A one-sided difference with a second order
C truncation error is used on the boundaries when possible - otherwise a one-sided
C difference with a first order truncation error is used.
C MLT 20March03 May require more thought when intracellular domain is
C discontinuous - however, conceptual arguments can be formed either way.
          PROPQ(1,i,1,2,nq) = 0.0d0 ! first variable
          PROPQ(2,i,1,2,nq) = 0.0d0 ! second variable (IF BIDOMAIN?)
          LENGTH = MAX(PROPQ(1,i,2,2,nq),PROPQ(2,i,2,2,nq))
          IF (NXQ(-i,0,nq).NE.0.AND.NXQ(i,0,nq).NE.0) THEN ! internal gp derivative
            IF (NXQ(-i,0,nq).GT.0.AND.NXQ(i,0,nq).GT.0) THEN ! 2 point 2nd order approx
              PROPQ(1,i,1,2,nq) = (YQ1(ABS(NXQ(i,1,nq)))-
     '          YQ1(ABS(NXQ(-i,1,nq))))/(2.0d0*LENGTH)
              IF (BIDOMAIN) PROPQ(2,i,1,2,nq) = (YQ2(ABS(NXQ(i,1,nq)))-
     '          YQ2(ABS(NXQ(-i,1,nq))))/(2.0d0*LENGTH)
            ELSEIF (NXQ(-i,0,nq).GT.0) THEN ! 2 point 1st order approximation
              PROPQ(1,i,1,2,nq) = (YQ1(nq)-
     '          YQ1(ABS(NXQ(-i,1,nq))))/LENGTH
              IF (BIDOMAIN) PROPQ(2,i,1,2,nq) = (YQ2(nq)-
     '          YQ2(ABS(NXQ(-i,1,nq))))/LENGTH
            ELSEIF (NXQ(i,0,nq).GT.0) THEN ! 2 point 1st order approximation
              PROPQ(1,i,1,2,nq) = (YQ1(ABS(NXQ(i,1,nq)))-
     '          YQ1(nq))/LENGTH
              IF (BIDOMAIN) PROPQ(2,i,1,2,nq) = (YQ2(ABS(NXQ(i,1,nq)))-
     '          YQ2(nq))/LENGTH
            ENDIF ! if no 2 point derivative approx can be formed, defaults to zero
          ELSEIF ((NXQ(-i,0,nq).EQ.0)) THEN ! boundary point
            IF (NXQ(i,1,nq).GT.0.AND.NXQ(i,1,ABS(NXQ(i,1,nq))).GT.0) 
     '        THEN ! 3 point one-sided difference (if sufficient grid points and 
c     '             ! no intracellular breaks)
              PROPQ(1,i,1,2,nq) = (-3.0d0*YQ1(nq) + 
     '          4.0d0*YQ1(ABS(NXQ(i,1,nq)))-
     '          YQ1(ABS(NXQ(i,1,ABS(NXQ(i,1,nq))))))/
     '          (2.0d0*LENGTH)
              IF (BIDOMAIN) PROPQ(2,i,1,2,nq) = (-3.0d0*YQ2(nq) + 
     '          4.0d0*YQ2(ABS(NXQ(i,1,nq)))-
     '          YQ2(ABS(NXQ(i,1,ABS(NXQ(i,1,nq))))))/
     '          (2.0d0*LENGTH)
            ELSEIF (NXQ(i,1,nq).GT.0) THEN ! 2 point one-sided difference (if 
c     '                                     ! sufficient grid points and no intracellular
c     '                                     ! breaks)
              PROPQ(1,i,1,2,nq) = (-1.0d0*YQ1(nq) + 
     '          1.0d0*YQ1(ABS(NXQ(i,1,nq))))/LENGTH
              IF (BIDOMAIN) PROPQ(2,i,1,2,nq) = (-1.0d0*YQ2(nq) + 
     '          1.0d0*YQ2(ABS(NXQ(i,1,nq))))/LENGTH
            ENDIF ! if a one-sided boundary derivative cannot be formed then it
c     '            ! defaults to zero
          ELSE ! boundary point
            IF (NXQ(-i,1,nq).GT.0.AND.NXQ(-i,1,ABS(NXQ(-i,1,nq))).GT.0) 
     '        THEN ! 3 point one-sided difference (if sufficient grid points and 
c     '             ! no intracellular breaks)
              PROPQ(1,i,1,2,nq) = (3.0d0*YQ1(nq) -
     '          4.0d0*YQ1(ABS(NXQ(-i,1,nq)))+
     '          YQ1(ABS(NXQ(-i,1,ABS(NXQ(-i,1,nq))))))/
     '          (2.0d0*LENGTH)
              IF (BIDOMAIN) PROPQ(2,i,1,2,nq) = (3.0d0*YQ2(nq) - 
     '          4.0d0*YQ2(ABS(NXQ(-i,1,nq)))+
     '              YQ2(ABS(NXQ(-i,1,ABS(NXQ(-i,1,nq))))))/
     '              (2.0d0*LENGTH)
            ELSEIF (NXQ(-i,1,nq).GT.0) THEN ! 2 point one-sided difference (if 
c     '                                      ! sufficient grid points and no intracellular
c     '                                      ! breaks)
              PROPQ(1,i,1,2,nq) = (1.0d0*YQ1(nq) -
     '          1.0d0*YQ1(ABS(NXQ(-i,1,nq))))/LENGTH
              IF (BIDOMAIN) PROPQ(2,i,1,2,nq) = (1.0d0*YQ2(nq) - 
     '          1.0d0*YQ2(ABS(NXQ(-i,1,nq))))/LENGTH
            ENDIF ! if a one-sided boundary derivative cannot be formed then it
c     '            ! defaults to zero
          ENDIF
        ENDDO ! coordinate direction

C Sum up appropriate RHS source contributions for grid points with 
C gradients that have already been calculated.
C Note: secondary contributions are only added if there is a 
C fluxable face in the -dir direction from grid point nq. If there
C is no fluxable face (either external boundary or intracellular break)
C then no contribution to the appropriate terms are made.
C Note: only information in the -i direction will have been already
C calculated since lower numbered grid points only lie in - directions.
C
C Note: due to the fact that only the +ve direction information is stored
C       for a grid point, both nq and the adjacent grid point secondary flux
C       contributions must be stored at the same time - thus this loop will
C       not parallelise well.  

C Loop over coordinate directions
        IF (BIDOMAIN) THEN
          DO i = 1,NITB
            PROPQ(1,i,4,2,nq) = 0.0d0 ! -ve direction for grid point nq
            PROPQ(2,i,4,2,nq) = 0.0d0
            IF (ABS(NXQ(-i,0,nq)).NE.0.AND.ABS(NXQ(-i,1,nq)).LT.nq) THEN
              nq2 = ABS(NXQ(-i,1,nq))
              PROPQ(1,i,3,2,nq2) = 0.0d0 ! +ve direction for grid point nq2
              PROPQ(2,i,3,2,nq2) = 0.0d0            
              DO j=1,NITB

                IF (NXQ(-i,1,nq).GT.0) THEN ! Intracellular contribution
C First bidomain equation
                  PROPQ(1,i,4,2,nq) = PROPQ(1,i,4,2,nq) +
     '              0.5d0*(PROPQ(1,j,1,2,nq)+PROPQ(1,j,1,2,nq2)+
     '                     PROPQ(2,j,1,2,nq)+PROPQ(2,j,1,2,nq2))*
     '              (-PROPQ(3,i,1,1,nq2)*PROPQ(1,j,i,1,nq2)-
     '                PROPQ(3,i,2,1,nq2)*PROPQ(2,j,i,1,nq2))

                  PROPQ(1,i,3,2,nq2) = PROPQ(1,i,3,2,nq2) +
     '              0.5d0*(PROPQ(1,j,1,2,nq)+PROPQ(1,j,1,2,nq2)+
     '                     PROPQ(2,j,1,2,nq)+PROPQ(2,j,1,2,nq2))*
     '              (PROPQ(3,i,1,1,nq2)*PROPQ(1,j,i,1,nq2)+
     '               PROPQ(3,i,2,1,nq2)*PROPQ(2,j,i,1,nq2))

C Second bidomain equation
                  PROPQ(2,i,4,2,nq) = PROPQ(2,i,4,2,nq) +
     '              0.5d0*(PROPQ(1,j,1,2,nq)+PROPQ(1,j,1,2,nq2)+
     '                     PROPQ(2,j,1,2,nq)+PROPQ(2,j,1,2,nq2))*
     '              (-PROPQ(3,i,1,1,nq2)*PROPQ(1,j,i,1,nq2)-
     '                PROPQ(3,i,2,1,nq2)*PROPQ(2,j,i,1,nq2))

                  PROPQ(2,i,3,2,nq2) = PROPQ(2,i,3,2,nq2) +
     '              0.5d0*(PROPQ(1,j,1,2,nq)+PROPQ(1,j,1,2,nq2)+
     '                     PROPQ(2,j,1,2,nq)+PROPQ(2,j,1,2,nq2))*
     '              (PROPQ(3,i,1,1,nq2)*PROPQ(1,j,i,1,nq2)+
     '               PROPQ(3,i,2,1,nq2)*PROPQ(2,j,i,1,nq2))
                ENDIF ! Intracellular contribution

C Second bidomain equation - extracellular contribution
                PROPQ(2,i,4,2,nq) = PROPQ(2,i,4,2,nq) +
     '            0.5d0*(PROPQ(2,j,1,2,nq)+PROPQ(2,j,1,2,nq2))*
     '            (-PROPQ(3,i,1,2,nq2)*PROPQ(1,j,i,1,nq2)-
     '              PROPQ(3,i,2,2,nq2)*PROPQ(2,j,i,1,nq2))          

                PROPQ(2,i,3,2,nq2) = PROPQ(2,i,3,2,nq2) +
     '            0.5d0*(PROPQ(2,j,1,2,nq)+PROPQ(2,j,1,2,nq2))*
     '            (PROPQ(3,i,1,2,nq2)*PROPQ(1,j,i,1,nq2)+
     '             PROPQ(3,i,2,2,nq2)*PROPQ(2,j,i,1,nq2))

              ENDDO ! coordinate direction: j
            ENDIF ! not external grid point and lower number than current
          ENDDO ! coordiante direction: i

        ELSE ! Monodomain equation, Anisotropic diffusion or Laplace's equation
          DO i = 1,NITB
            PROPQ(1,i,4,2,nq) = 0.0d0 ! -ve direction for grid point nq
            IF (ABS(NXQ(-i,0,nq)).NE.0.AND.ABS(NXQ(-i,1,nq)).LT.nq) THEN
              nq2 = ABS(NXQ(-i,1,nq))
              PROPQ(1,i,3,2,nq2) = 0.0d0 ! +ve direction for grid point nq2
              DO j=1,NITB
                PROPQ(1,i,4,2,nq) = PROPQ(1,i,4,2,nq) +
     '            0.5d0*(PROPQ(1,j,1,2,nq)+PROPQ(1,j,1,2,nq2))*
     '            (-PROPQ(3,i,1,1,nq2)*PROPQ(1,j,i,1,nq2)-
     '              PROPQ(3,i,2,1,nq2)*PROPQ(2,j,i,1,nq2))

                PROPQ(1,i,3,2,nq2) = PROPQ(1,i,3,2,nq2) +
     '            0.5d0*(PROPQ(1,j,1,2,nq)+PROPQ(1,j,1,2,nq2))*
     '            (PROPQ(3,i,1,1,nq2)*PROPQ(1,j,i,1,nq2)+
     '             PROPQ(3,i,2,1,nq2)*PROPQ(2,j,i,1,nq2))
              ENDDO ! coordinate direction: j
            ENDIF ! not external grid point and lower number than current
          ENDDO ! coordiante direction: i

        ENDIF ! Bidomain or not

      ENDDO ! nq (grid point)
C$OMP END PARALLEL DO

      CALL EXITS('CALC_FV_GRID_SECFLUX')
      RETURN
 9999 CALL ERRORS('CALC_FV_GRID_SECFLUX',ERROR)
      CALL EXITS('CALC_FV_GRID_SECFLUX')
      RETURN 1
      END      

