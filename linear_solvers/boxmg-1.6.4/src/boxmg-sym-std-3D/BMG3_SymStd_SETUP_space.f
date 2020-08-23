      SUBROUTINE BMG3_SymStd_SETUP_space( 
     &                Nx, Ny, Nz, BMG_iPARMS,
     &                NOG, NF, NC, NSO, NSOR, NCI, NCBW, NCU,
     &                NCG_TEMP, NBMG_iWORK_PL, NBMG_rWORK_PL,
     &                BMG_IOFLAG
     &                )

C ========================================================================
C 
C   BMG3_SymStd_SETUP_space
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_SETUP_space is used to compute the storage requirements 
C   for the various arrays that are used by BMG3D.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:  2000/08/29 (JDM)
C
C =======================================================================
C   PARAMETERS:
C ========================
C
C   pNx   - index of x-direction dimensions
C   pNy   - index of y-direction dimensions
C   pNz   - index of z-direction dimensions
C   pN3   - 3D dimensions
C
C =======================================================================
C   INPUT:
C ========================
C
C   ----------------------
C   Fine Grid Dimensions:
C   ---------------------- 
C
C   Nx       Number of points in x-direction (excluding ghost points)
C   Ny       Number of points in y-direction (excluding ghost points)
C   Nz       Number of points in z-direction (excluding ghost points)
C
C  ------------------------
C   Discretization:
C  ------------------------
C
C   iSTNCL   Discrete operator index:
C            iSTNCL .EQ. 1 => a 7  point discretization
C            iSTNCL .NE. 1 => a 27 point discretization
C
C   iBC      Boundary condition index:
C            - unused
C
C   ------------------------
C   Smoothing:
C   ------------------------   
C
C   IRELAX   Relaxation index (refer to BMG3D).
C
C =======================================================================
C   OUTPUT:
C ===========================
C
C   NOG       Number of grids needed for the given (Nx,Ny,Nz)
C   NF        Storage for arrays Q and QF on all grids
C   NC        Storage for Q and QF on coarser grids
C             - Q the solution vector on all 3D grids
C             - QF the source vector on all 3D grids
C
C   NSO       Storage for the array SO
C             - SO holds the stencil on all 3D grids.
C   NSOR      Storage for the array SOR: 
C             - SOR holds the current residual, the reciprocal of 
C               the central stencil weight, and the LU decomposition
C               if iRELAX > 1.
C   NCI       Storage for the array CI
C             - CI holds the 3D interpolation weights on all 3D grids
C
C   NCBW      First dimension of ABD, which is set to the bandwidth
C             of the coarsest grid stencil (including the diagonal).
C   NCU       Second dimension of ABD, the dimension of the coarse grid.
C
C   --------------------------
C   Plane Relaxation:
C   --------------------------
C 
C   NF2       Storage for arrays Q2 and QF2 on all grids.
C
C   NSO3      Storage needed for SO2.
C   NSR3      Storage needed for SOR2.
C   NCI3      Storage needed for CI2.
C   NCBW2     First dimension of ABD2, which is the bandwidth
C             of the 2D coarsest grid stencils (including the diagonal).
C   NCU2      Second dimension of  ABD2, number of unknowns
C             on the largest planes coarsest grid.
C   NCPL2     Third dimension of ABD2, which is the number of 
C             coarse grid 2D planes
C
C =======================================================================
C   LOCAL:
C ========================
C
C   Nxc       Number of points in the x-direction on a coarser grid
C   Nyc       Number of points in the y-direction on a coarser grid
C   Nzc       Number of points in the z-direction on a coarser grid
C
C   NPLx      Number of (y,z) planes on a coarser grid
C   NPLy      Number of (x,z) planes on a coarser grid
C   NPLz      Number of (x,y) planes on a coarser grid
C
C   Nx2       -no idea
C   Ny2       -no idea
C
C   NC2       Storage for arrays Q2 and QF2 on coarser grids.
C
C   NSO2      Storage for SO2 on the largest fine-grid 2D plane
C   NSR2      Storage for SOR2 on the largest fine-grid 2D plane
C   NCI2      Storage for CI2 on the largest fine-grid 2D plane
C
C   kg        Current grid index 
C
C   iTEMP     Temp storage
C
C =======================================================================

      IMPLICIT NONE

C ---------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'
      
C ---------------------------
C     Parameter Declarations:
C
      INTEGER  pNx, pNy, pNz, pN3
      PARAMETER ( pNx=1, pNy=2, pNz=3, pN3=4 )

C ---------------------------
C     Argument Declarations:
C
      INTEGER  BMG_iPARMS(NBMG_iPARMS), 
     &         NBMG_iWORK_PL, NBMG_rWORK_PL, NCG_TEMP

      INTEGER  iBC, iCG, iSTNCL, iRELAX, NC, NCBW, NCU, NCI, NOG,
     &         NF, NSO, NSOR, Nx, Ny, Nz, NXYZc

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ---------------------------
C     Local Declarations:
C
      INTEGER  i, iSWITCH, kg, NBMG_iWORK_OP, NBMG_rWORK_OP,
     &         NCG_TEMP_ZO, NCG_TEMP_YO, NIGRD3,
     &         Nxc, Nyc, Nzc, NXYc

      INTEGER  BMG_pWORK_2D(NBMG_pWORK),
     &         NBMG_iWORK_2D, NBMG_rWORK_2D, 
     &         NFm_2D, NOGm_2D, NSOm_2D,
     &         pSI_2D, pSR_2D
      LOGICAL  BMG_InWORK_2D(NBMG_InWORK)

C =======================================================================

C =======================================================================
C --------------------------------------------
C     Unpack Parameters
C --------------------------------------------

      NXYZc   = BMG_iPARMS(id_BMG3_CG_MIN_DIM)
      IStncl  = BMG_iPARMS(id_BMG3_STENCIL)
      IRELAX  = BMG_iPARMS(id_BMG3_RELAX)
      IBC     = BMG_iPARMS(id_BMG3_BC)
      iCG     = BMG_iPARMS(id_BMG3_CG_TYPE)

      NXYc    = BMG_iPARMS(id_BMG2_CG_MIN_DIM)

C --------------------------------------------
C     Check coarsening criteria:
C --------------------------------------------

      IF(NXYZc.LT.3)THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Minimum coarse grid dimension is less than 3'
            WRITE(*,520) 'HAVE: NXYZc = ', NXYZc
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,11)
         RETURN
 
      ENDIF


      IF(NXYc.LT.3)THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Minimum coarse grid dimension is less than 3'
            WRITE(*,520) 'HAVE: NXYc = ', NXYc
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,12)
         RETURN

      ENDIF

C --------------------------------------------
C    Check relaxation index:
C --------------------------------------------

      IF ( iRELAX.NE.BMG_GS_RB_point 
     &    .AND. iRELAX.NE.BMG_GS_RB_planes_xy_yz_xz 
     &     .AND. iRELAX.NE.BMG_CG_relaxation 
     &     .AND. iRELAX.NE.BMG_GS_adaptive ) THEN 

         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The value of iRELAX is not supported ... '
            WRITE(*,520) '  HAVE: iRELAX = ', iRELAX
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,13)
         RETURN

      ENDIF
      
C --------------------------------------------
C =======================================================================

C =======================================================================
C -------------------------------------
C     Compute the number of 3D grids:
C -------------------------------------

      kg=1


 10   CONTINUE
      
         Nxc = (Nx-1)/2**kg+1
         Nyc = (Ny-1)/2**kg+1
         Nzc = (Nz-1)/2**kg+1

         IF (MIN(Nxc,Nyc,Nzc).GE.NXYZc) THEN   
            kg=kg+1
            GOTO 10
         ENDIF
        
 20   CONTINUE

      NOG=kg

C --------------------------------------------
C     Dimensions of the coarsest grid:
C --------------------------------------------
      
      IF ( NOG.LT.BMG_iPARMS(id_BMG3_MIN_NOG)) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The computed number of grids is too small!'
            WRITE(*,510) 'HAVE: BMG_iPARMS(id_BMG3_MIN_NOG)    = ', 
     &           BMG_iPARMS(id_BMG3_MIN_NOG)
            WRITE(*,510) 'HAVE: BMG_iPARMS(id_BMG3_CG_MIN_DIM) = ',
     &           BMG_iPARMS(id_BMG3_CG_MIN_DIM)
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,14)
         RETURN

      ELSEIF ( NOG.GT.1 ) THEN
         Nxc = (Nx-1)/2**(NOG-1)+1
         Nyc = (Ny-1)/2**(NOG-1)+1
         Nzc = (Nz-1)/2**(NOG-1)+1
      ELSEIF ( NOG.EQ.1 ) THEN
         Nxc = Nx
         Nyc = Ny
         Nzc = Nz
      ELSE
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The computed number of grids is less than 1!'
            WRITE(*,510) 'HAVE: Fine grid dimension Nx = ', Nx
            WRITE(*,510) 'HAVE: Fine grid dimension Ny = ', Ny
            WRITE(*,510) 'HAVE: Fine grid dimension Nz = ', Nz
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,15)
         RETURN

      ENDIF
         
C --------------------------------------------
C    Grid Indexing and pointers:
C --------------------------------------------

      NIGRD3 = NBMG_pIGRD*NOG

C --------------------------------------------
C     Storage space for all 3D grids:
C --------------------------------------------

      NF=(Nx+2)*(Ny+2)*(Nz+2)
      DO 25 kg=1, NOG-1
         NF = NF + (3+(Nx-1)/2**kg)*(3+(Ny-1)/2**kg)*(3+(Nz-1)/2**kg)
 25   CONTINUE

C ------------------------------------------
C     Storage for the coarse grids (3D):
C ------------------------------------------

      NC=NF-(Nx+2)*(Ny+2)*(Nz+2)

C ------------------------------------------
C     Storage for the 3D stencil:
C ------------------------------------------

      IF (iSTNCL.EQ.1) THEN
         NSO = 4*NF + 10*NC
      ELSE
         NSO = 14*NF
      ENDIF

C -------------------------------------------
C     Storage for the 3D interpolation:
C -------------------------------------------
      
      NCI = 26*NC

C -------------------------------------------
C     Solving the coarsest grid in 3D
C -------------------------------------------

      NCU  = Nxc*Nyc*Nzc
      NCBW = Nxc*(Nyc+1)+2

C -------------------------------------------
C     Zero workspace estimates:
C -------------------------------------------

      NBMG_iWORK_OP = 0 
      NBMG_rWORK_OP = 0

      NBMG_iWORK_PL = 0
      NBMG_rWORK_PL = 0 

C ------------------------------------------
C      Workspace for the operator setup:
C ------------------------------------------

      !
      ! Borrow space from plane relaxation work array
      !
      Nyc = (Ny-1)/2+1
      Nzc = (Nz-1)/2+1
      
      IF( iCG.EQ.BMG_CG_ITLI_IzIyIx) THEN
         NCG_TEMP_ZO = 14*(Nx+2)*(Ny+2)*(Nzc+3) 
         NCG_TEMP_YO = 14*(Nx+2)*(Nyc+3)*(Nzc+3)
      ELSE
         NCG_TEMP_ZO = 1
         NCG_TEMP_YO = 1
      ENDIF

      NCG_TEMP = NCG_TEMP_YO
      
      NBMG_iWORK_OP = NBMG_iWORK_PL_ptrs
      NBMG_rWORK_OP = NCG_TEMP_YO + NCG_TEMP_ZO

C ------------------------------------------
C     Work space for the Relaxations:
C ------------------------------------------

      IF ( iRELAX.EQ.BMG_GS_RB_point ) THEN
         
         NSOR = 2*NF

      ELSEIF ( iRELAX.EQ.BMG_GS_adaptive .OR. 
     &         iRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN

         NSOR = 2*NF            ! just to be safe for now, probably can
                       ! can get rid of this. 

         !
         ! Allocate space for pointers
         !
         NBMG_iWORK_PL = NBMG_iWORK_PL + NBMG_iWORK_PL_ptrs

         !
         ! Allocation for the pointer array BMG_pWORK_PL
         !  - one for every plane on every level 
         !    ( this is an overestimate, but should be easier )
         !
         NBMG_iWORK_PL = NBMG_iWORK_PL 
     &                 + NBMG_pWORK*NOG*((Nx+2)+(Ny+2)+(Nz+2))

         !
         ! Allocation for the parameter arrays BMG_iPARMS
         !  - one for each orientation
         !
         NBMG_iWORK_PL = NBMG_iWORK_PL + NBMG_iPARMS*NOG*3

         !
         ! We just need dimensional information at this stage
         !
         BMG_iPARMS(id_BMG2_POINTERS)=BMG_USE_pointers

         !
         ! Initialize irrelevant data
         !
         NOGm_2D = 1 
         NFm_2D = 1
         NSOm_2D = 1
         
         pSR_2D = 1
         pSI_2D = 1


         !
         ! NB: 
         ! 
         !
         ! For the first cut lets put it all in workspace and 
         ! not worry about the duplication of the solution, RHS 
         ! and residual space.
         !
         ! 
         ! At some point we'll consider having everything in worspace
         ! for the first plane on each level and in each orientation
         ! However, then we'll reuse the solution, RHS, and residual 
         ! vectors for subsequent planes on the same level and with 
         ! the same orientation.
         !

         !
         !  Initialize BMG_InWORK_PL
         !
         DO i=1, NBMG_InWORK
            BMG_InWORK_2D(i) = .FALSE.
         ENDDO

         !
         !  Override the BoxMG arguments
         !
         BMG_InWORK_2D(i_InWORK_SO)  = .TRUE.
         BMG_InWORK_2D(i_InWORK_U)   = .TRUE.
         BMG_InWORK_2D(i_InWORK_Q)   = .TRUE.
         BMG_InWORK_2D(i_InWORK_RES) = .TRUE.

         !
         ! (x,y) - planes
         !
         DO kg=2, NOG
            !
            !  Set stencil size
            !
            IF (kg.NE.NOG .OR. IStncl.NE.1 ) THEN
               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_9pt
            ELSE
               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_5pt
            ENDIF
            !
            Nxc = (Nx-1)/2**(NOG-kg)+1
            Nyc = (Ny-1)/2**(NOG-kg)+1
            Nzc = (Nz-1)/2**(NOG-kg)+1
            !
            ! Set error codes to zero
            !
            BMG_iPARMS(id_BMG2_Err_Code) = 0
            BMG_iPARMS(id_BMG2_Ext_Err_Code) = 0
            !
            CALL BMG2_SymStd_SETUP_PtrWork( 
     &                       Nxc, Nyc, BMG_iPARMS, 
     &                       NOGm_2D, NFm_2D, NSOm_2D,
     &                       NBMG_iWORK_2D, NBMG_rWORK_2D,
     &                       BMG_pWORK_2D, BMG_InWORK_2D,
     &                       pSR_2D, pSI_2D, BMG_IOFLAG 
     &                       )
            NBMG_iWORK_PL = NBMG_iWORK_PL + NBMG_iWORK_2D*Nzc
            NBMG_rWORK_PL = NBMG_rWORK_PL + NBMG_rWORK_2D*Nzc
            !
         ENDDO

         !
         ! (y,z) - planes
         !
         DO kg=2, NOG
            !
            !  Set stencil size
            !
            IF (kg.NE.NOG .OR. IStncl.NE.1 ) THEN
               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_9pt
            ELSE
               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_5pt
            ENDIF
            !
            Nxc = (Nx-1)/2**(NOG-kg)+1
            Nyc = (Ny-1)/2**(NOG-kg)+1
            Nzc = (Nz-1)/2**(NOG-kg)+1
            !
            ! Set error codes to zero
            !
            BMG_iPARMS(id_BMG2_Err_Code) = 0
            BMG_iPARMS(id_BMG2_Ext_Err_Code) = 0
            !
            CALL BMG2_SymStd_SETUP_PtrWork( 
     &                       Nyc, Nzc, BMG_iPARMS, 
     &                       NOGm_2D, NFm_2D, NSOm_2D,
     &                       NBMG_iWORK_2D, NBMG_rWORK_2D,
     &                       BMG_pWORK_2D, BMG_InWORK_2D,
     &                       pSR_2D, pSI_2D, BMG_IOFLAG 
     &                       )
            NBMG_iWORK_PL = NBMG_iWORK_PL + NBMG_iWORK_2D*Nxc
            NBMG_rWORK_PL = NBMG_rWORK_PL + NBMG_rWORK_2D*Nxc
            !
         ENDDO

         !
         ! (x,z) - planes
         !
         DO kg=2, NOG
            !
            !  Set stencil size
            !
            IF (kg.NE.NOG .OR. IStncl.NE.1 ) THEN
               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_9pt
            ELSE
               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_5pt
            ENDIF
            !
            Nxc = (Nx-1)/2**(NOG-kg)+1
            Nyc = (Ny-1)/2**(NOG-kg)+1
            Nzc = (Nz-1)/2**(NOG-kg)+1
            !
            ! Set error codes to zero
            !
            BMG_iPARMS(id_BMG2_Err_Code) = 0
            BMG_iPARMS(id_BMG2_Ext_Err_Code) = 0
            !
            CALL BMG2_SymStd_SETUP_PtrWork( 
     &                       Nxc, Nzc, BMG_iPARMS, 
     &                       NOGm_2D, NFm_2D, NSOm_2D,
     &                       NBMG_iWORK_2D, NBMG_rWORK_2D,
     &                       BMG_pWORK_2D, BMG_InWORK_2D,
     &                       pSR_2D, pSI_2D, BMG_IOFLAG 
     &                       )
            NBMG_iWORK_PL = NBMG_iWORK_PL + NBMG_iWORK_2D*Nyc
            NBMG_rWORK_PL = NBMG_rWORK_PL + NBMG_rWORK_2D*Nyc
            !
         ENDDO

      ELSEIF ( iRELAX.EQ.BMG_CG_relaxation ) THEN

         NSOR = 4*NF

      ENDIF

C ------------------------------------------
C     Dimension of shared workspace:
C ------------------------------------------

      NBMG_iWORK_PL = MAX( NBMG_iWORK_OP, NBMG_iWORK_PL )
      NBMG_rWORK_PL = MAX( NBMG_rWORK_OP, NBMG_rWORK_PL )

C ========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_space.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
