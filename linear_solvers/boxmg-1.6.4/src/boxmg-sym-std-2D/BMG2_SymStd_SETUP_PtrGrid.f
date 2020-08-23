      SUBROUTINE BMG2_SymStd_SETUP_PtrGrid( 
     &                Nx, Ny, BMG_iPARMS,
     &                NOGm, NFm, NCm, NCIm, NSOm, NSORm, NCBWm, NCUm,
     &                NOG, NF, NC, NCI, NSO, NSOR, NCBW, NCU, IGRD,
     &                BMG_IOFLAG
     &                )

C =======================================================================
C 
C   BMG2_SymStd_SETUP_PtrGrid
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_PtrGrid is used to compute the pointers into 
C   the various arrays for each grid within BOXMG.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/07/24 (JDM)
C
C =======================================================================
C   INPUT:
C ========================
C
C   ----------------------
C   Fine Grid Dimensions:
C   ---------------------- 
C
C   Nx         Number of points in x-direction (excluding ghost points)
C   Ny         Number of points in y-direction (excluding ghost points)
C
C   -----------------------
C   Coarsening Criteria:
C   -----------------------
C
C   NXYc       Minimum number of coarse grid points in either direction
C
C  ------------------------
C   Discretization:
C  ------------------------
C
C   iSTNCL     Discrete operator index:
C              iSTNCL .EQ. 1 => a 5  point discretization
C              iSTNCL .NE. 1 => a 9  point discretization
C
C   ------------------------
C   Smoothing:
C   ------------------------   
C
C   iRELAX     Relaxation index (refer to BOXMG).
C
C   ------------------------
C   Boundary Conditions:
C   ------------------------
C
C   iBC        Boundary condition index (refer to BOXMG,BOXMGP)
C
C   --------------------------------------------
C   (PRE)computed dimensions (which should be):
C   --------------------------------------------
C
C   NOG       Number of grids needed for the given (Nx,Ny)
C
C   NF        Storage for arrays Q and QF on all grids
C   NC        Storage for Q and QF on coarser grids
C             - Q the solution vector on all 2D grids
C             - QF the source vector on all 2D grids
C
C   NSO       Storage for the array SO
C             - SO holds the stencil on all 2D grids.
C   NSOR      Storage for the array SOR: 
C             - SOR holds the current residual, the reciprocal of 
C               the central stencil weight, and the LU decomposition
C               if iRELAX > 1.
C   NCI       Storage for the array CI
C             - CI holds the interpolation weights on all 2D grids
C
C   NCBW      First dimension of ABD, which is set to the bandwidth
C             of the coarsest grid stencil (including the diagonal).
C   NCU       Second dimension of ABD, the dimension of the coarse grid.
C
C =======================================================================
C   OUTPUT:
C ===========================
C
C   IGRD      Grid pointer array for the arrays internal to BOXMG
C
C =======================================================================
C   LOCAL:
C ========================
C
C   Nxc       Number of points in the x-direction on a coarser grid
C   Nyc       Number of points in the y-direction on a coarser grid
C
C   kg        Grid index
C
C =======================================================================

      IMPLICIT   NONE

      INCLUDE    'BMG_workspace.h'
      INCLUDE    'BMG_constants.h'
      INCLUDE    'BMG_parameters.h'

C ----------------------------
C     Parameter Declarations:
C
      INTEGER   PerGrid
      PARAMETER ( PerGrid = 9 )

C ---------------------------
C     Argument Declarations:
C
      INTEGER  BMG_iPARMS(NBMG_iPARMS) 

      INTEGER  NCBWm, NCIm, NCm, NCUm, NOGm, NFm, NSOm, NSORm

      INTEGER  iBC, IGRD(NOGm,PerGrid), iRELAX, iSTNCL, NC, NCBW, NCU,
     &         NCI, NOG, NF, NSO, NSOR, Nx, Ny, NXYc

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG) 

C ---------------------------
C     Local Declarations:
C
      INTEGER  CF, IERR, i, j, kg, Nxc, Nxg, NxNy, Nyc, Nyg

C ========================================================================

C ========================================================================
C --------------------------------------------
C     Unpack Parameters
C --------------------------------------------

      NXYc   = BMG_iPARMS(id_BMG2_CG_MIN_DIM)
      IStncl = BMG_iPARMS(id_BMG2_STENCIL)
      IRELAX = BMG_iPARMS(id_BMG2_RELAX)
      IBC    = BMG_iPARMS(id_BMG2_BC)

C --------------------------------------------
C =======================================================================

C --------------------------------------------
C    Compute the space requirements:
C --------------------------------------------

      CALL BMG2_SymStd_SETUP_space( 
     &          Nx, Ny, BMG_iPARMS,
     &          NOG, NF, NC, NCI, NSO, NSOR, NCBW, NCU
     &          )

C --------------------------------------------
C    Check the actual allocation:
C --------------------------------------------

      IERR = 0
      !
      !  The number of grids
      !
      IF ( NOG.GT.NOGm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'The number of grids is greater than NOGm'
         WRITE(*,510) 'HAVE: NOGm = ', NOGm
         WRITE(*,520) 'NEED: NOGm = ', NOG
      ENDIF

      !
      !  Total storage for a vector on all grids
      !
      IF ( NF.GT.NFm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Total grid space allocation is too small!'
         WRITE(*,510) 'HAVE: NF = ', NFm
         WRITE(*,520) 'NEED: NF = ', NF
      ENDIF

      !
      !  Total storage for a vector on all coarse grids
      !
      IF ( NC.GT.NCm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Total grid space allocation is too small!'
         WRITE(*,510) 'HAVE: NC = ', NCm
         WRITE(*,520) 'NEED: NC = ', NC
      ENDIF

      !
      !  The dimension of the stencil array
      !  -> NSO is IDISO 
      !
      IF ( NSO.GT.NSOm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Stencil allocation is too small!'
         WRITE(*,510) 'HAVE: NSO = ', NSOm
         WRITE(*,520) 'NEED: NSO = ', NSO
      ENDIF

      !
      !  Workspace for line relaxation's factorization
      !  -> NSOR is IDISOR
      !
      IF ( NSOR.GT.NSORm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Factorization allocation is too small!'
         WRITE(*,510) 'HAVE: NSOR = ', NSORm
         WRITE(*,520) 'NEED: NSOR = ', NSOR
      ENDIF

      !
      !  The dimension of the interpolation operator
      !  -> NCI is IDICI
      !
      IF ( NCI.GT.NCIm ) THEN
         IERR = IERR+1
         WRITE(*,500) 'Interpolation allocation is too small!'
         WRITE(*,510) 'HAVE: NCI = ', NCIm
         WRITE(*,520) 'NEED: NCI = ', NCI
      ENDIF

      !
      !  The bandwidth for the direct solve on the coarsest grid
      !  -> NCBW is NABD1 
      !
      IF ( NCBW.GT.NCBWm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Coarsest grid band width is too small!'
         WRITE(*,510) 'HAVE: NCBW = ', NCBWm
         WRITE(*,520) 'NEED: NCBW = ', NCBW
      ENDIF

      !
      !  The dimension of the coarsest grid needed for the direct solve
      !  -> NCU is NABD2
      !
      IF ( NCU.GT.NCUm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Coarsest grid space allocation is too small!'
         WRITE(*,510) 'HAVE: NCU = ', NCUm
         WRITE(*,520) 'NEED: NCU = ', NCU
      ENDIF
      
      !
      !  If there were errors STOP
      !
      IF ( IERR.GT.0 ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,530)IERR,'FATAL ERRORS in BMG2_SymStd_SETUP_PtrGrid,//
     &          // EXECUTION STOPPED'
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,7)
         RETURN

      ENDIF
      
C --------------------------------------------
C ========================================================================


C ========================================================================
C --------------------------------------------
C     Initialize the indexing array.
C --------------------------------------------

      DO j=1, PerGrid
         DO i=1, NOGm
            IGRD(i,j)= 0
         END DO
      END DO

C -------------------------------------------
C     Finest Grid:
C -------------------------------------------

      Nxg  = Nx+2
      Nyg  = Ny+2
      NxNy = Nxg*Nyg

      !
      ! Dimension; finest grid
      !
      IGRD(NOG,idL_BMG_Nx) = Nxg
      IGRD(NOG,idL_BMG_Ny) = Nyg
      !
      ! Pointer: vectors
      !
      IGRD(NOG,ipL_BMG_U)   = 1
      !
      ! Pointer: stencil
      !
      IGRD(NOG,ipL_BMG_SO)  = 1
      !
      ! Pointer: interpolation
      !
      IGRD(NOG,ipL_BMG_CI)   = 1
      !
      ! Pointer: relaxation workspace
      !
      IGRD(NOG,ipL_BMG_SOR)  = 1

C -------------------------------------------
C     Is there only one grid? => RETURN
C -------------------------------------------

      IF ( NOG.EQ.1 ) THEN
         RETURN
      ENDIF

C -------------------------------------------
C     First Coarse Grid:
C -------------------------------------------      

      !
      ! Pointer: vectors
      !
      IGRD(NOG-1,ipL_BMG_U) = IGRD(NOG,ipL_BMG_U) + NxNy
      !
      ! Pointer: stencil
      !
      IF (iSTNCL.EQ.1) THEN
         IGRD(NOG-1,ipL_BMG_SO) = IGRD(NOG,ipL_BMG_SO)+3*NxNy
      ELSE
         IGRD(NOG-1,ipL_BMG_SO) = IGRD(NOG,ipL_BMG_SO)+5*NxNy
      ENDIF
      !
      ! Pointer: interpolation
      !
      IGRD(NOG-1,ipL_BMG_CI) = 1
      !
      ! Pointer: relaxation workspace
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point
     &    .OR. IRELAX.EQ.BMG_GS_RB_x_lines ) THEN
         IGRD(NOG-1,ipL_BMG_SOR) = IGRD(NOG,ipL_BMG_SOR)+2*NxNy
      ELSE IF ( IRELAX.EQ.BMG_GS_RB_y_lines 
     &         .OR. IRELAX.EQ.BMG_GS_RB_x_y_lines ) THEN
         IGRD(NOG-1,ipL_BMG_SOR) = IGRD(NOG,ipL_BMG_SOR)+4*NxNy
      ELSE
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'A parameter is out of range!'
            WRITE(*,520) 'IRELAX is out of range [1,4] = ', IRELAX
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,8)
         RETURN

      ENDIF

C ----------------------------------
C     Remaining Coarse Grids:
C ----------------------------------

      CF = 1
      DO kg = NOG-1, 2, -1
         !
         ! Coarsening: previous grid
         !
         CF   = 2*CF
         Nxg  = (Nx-1)/CF + 3
         Nyg  = (Ny-1)/CF + 3
         NxNy = Nxg*Nyg
         !
         ! Dimension: previous grid
         !
         IGRD(kg,idL_BMG_Nx) = Nxg
         IGRD(kg,idL_BMG_Ny) = Nyg
         !
         ! Pointer: vectors
         !
         IGRD(kg-1,ipL_BMG_U) = IGRD(kg,ipL_BMG_U)+NxNy
         !
         ! Pointer: stencil
         !
         IGRD(kg-1,ipL_BMG_SO) = IGRD(kg,ipL_BMG_SO)+5*NxNy
         !
         ! Pointer: interpolation
         !
         IGRD(kg-1,ipL_BMG_CI) = IGRD(kg,ipL_BMG_CI)+8*NxNy
         !
         ! Pointer: relaxation workspace
         !
         IF ( IRELAX.EQ.BMG_GS_RB_point 
     &       .OR. IRELAX.EQ.BMG_GS_RB_x_lines ) THEN
            IGRD(kg-1,ipL_BMG_SOR) = IGRD(kg,ipL_BMG_SOR)+2*NxNy
         ELSE IF ( IRELAX.EQ.BMG_GS_RB_y_lines 
     &            .OR. IRELAX.EQ.BMG_GS_RB_x_y_lines ) THEN
            IGRD(kg-1,ipL_BMG_SOR) = IGRD(kg,ipL_BMG_SOR)+4*NxNy
         ELSE
            IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
               WRITE(*,500) 'A parameter is out of range!'
               WRITE(*,520) 'IRELAX is out of range [1,4] = ', IRELAX
            END IF

            CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,8)
            RETURN
         ENDIF
         !
      ENDDO

      !
      !  Coarsening: coarsest grid
      !
      CF=2*CF
      Nxc = (Nx-1)/CF + 3 
      Nyc = (Ny-1)/CF + 3
      !
      !  Dimension: coarsest grid
      !
      IGRD(1,idL_BMG_Nx) = Nxc
      IGRD(1,idL_BMG_Ny) = Nyc

C ---------------------------------------------
C     Check Dimensions
C ---------------------------------------------

C      Need to add code to double check the workspace pointers are 
C      in bounds.

C --------------------------------------------
C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_PtrGrid.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
