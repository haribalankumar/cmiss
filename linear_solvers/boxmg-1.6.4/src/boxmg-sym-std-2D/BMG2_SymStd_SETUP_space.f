      SUBROUTINE BMG2_SymStd_SETUP_space( 
     &                Nx, Ny, BMG_iPARMS, 
     &                NOG, NF, NC, NCI, NSO, NSOR, NCBW, NCU
     &                )

C =======================================================================
C 
C   BMG2_SymStd_SETUP_space
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_space is used to compute the space requirements
C   for the various arrays that are use by BOXMG.
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
C   Nx       Number of points in x-direction (excluding ghost points)
C   Ny       Number of points in y-direction (excluding ghost points)
C
C
C   -----------------------
C   BMG_iPARMS:
C   -----------------------
C
C   * Coarsening Criteria:
C
C   NXYc      Minimum number of coarse grid points in either direction
C
C   * Discretization:
C
C   iSTNCL   Discrete operator index:
C            iSTNCL .EQ. 1 => a 5  point discretization
C            iSTNCL .NE. 1 => a 9  point discretization
C
C   * Smoothing
C
C   iRELAX   Relaxation index (refer to BOXMG).
C
C   * Boundary Conditions:
C
C   iBC      Boundary condition index (refer to BOXMG,BOXMGP)
C
C =======================================================================
C   OUTPUT:
C ===========================
C
C   NOG       Number of grids needed for the given (Nx,Ny) based on NXYc
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

C ---------------------------
C     Argument Declarations:
C
C ---------------------------
C     Argument Declarations:
C
      INTEGER  BMG_iPARMS(NBMG_iPARMS) 

      INTEGER  iBC, iRELAX, iSTNCL, NC, NCBW, NCU, NCI, NOG, NF,
     &         NSO, NSOR, Nx, Ny, NXYc

C ---------------------------
C     Local Declarations:
C
      INTEGER  kg, Nxc, Nyc

C ========================================================================

C =======================================================================

C =======================================================================
C --------------------------------------------
C     Unpack Parameters
C --------------------------------------------

      NXYc    = BMG_iPARMS(id_BMG2_CG_MIN_DIM)
      IStncl  = BMG_iPARMS(id_BMG2_STENCIL)
      IRELAX  = BMG_iPARMS(id_BMG2_RELAX)
      IBC     = BMG_iPARMS(id_BMG2_BC)

C --------------------------------------------
C     Check coarsening criteria:
C --------------------------------------------

      IF(NXYc.LT.3)THEN
!         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Minimum coarse grid dimension is less than 3'
            WRITE(*,520) 'HAVE: NXYc = ', NXYc
!         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,4)
         RETURN

      ENDIF

C =======================================================================

C ========================================================================
C -----------------------------------
C     Compute the number of grids:
C -----------------------------------

      kg = 1
 10   CONTINUE
         Nxc = (Nx-1)/2**kg+1
         Nyc = (Ny-1)/2**kg+1
         IF ( MIN(Nxc,Nyc).GE.NXYc ) THEN
            kg = kg + 1
            GOTO 10
         ENDIF
 20   CONTINUE

      NOG = kg

C --------------------------------------------
C     Dimensions of the coarsest grid:
C --------------------------------------------
      
      IF ( NOG.LT.BMG_iPARMS(id_BMG2_MIN_NOG)) THEN
!         IF (BMG_IOFLAG(IBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The computed number of grids is too small!'
            WRITE(*,510) 'HAVE: BMG_iPARMS(id_BMG2_MIN_NOG)    = ', 
     &           BMG_iPARMS(id_BMG2_MIN_NOG)
            WRITE(*,510) 'HAVE: BMG_iPARMS(id_BMG2_CG_MIN_DIM) = ',
     &           BMG_iPARMS(id_BMG2_CG_MIN_DIM)
            WRITE(*,520) 'HAVE: NOG = ', NOG
!         End IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,5)
         RETURN

      ELSEIF ( NOG.GT.1 ) THEN
         Nxc = (Nx-1)/2**(NOG-1)+1
         Nyc = (Ny-1)/2**(NOG-1)+1
      ELSEIF ( NOG.EQ.1 ) THEN
         Nxc = Nx
         Nyc = Ny
      ELSE
!         IF (BMG_IOFLAG(IBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The computed number of grids is less than 1!'
            WRITE(*,510) 'HAVE: Fine grid dimension Nx = ', Nx
            WRITE(*,510) 'HAVE: Fine grid dimension Ny = ', Ny
            WRITE(*,520) 'HAVE: NOG = ', NOG
!         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,6)
         RETURN
         
      ENDIF
         
C -------------------------------------
C     Storage for all grids: NF
C -------------------------------------

      NF=(Nx+2)*(Ny+2)
      DO 30 kg=1, NOG-1
         NF=NF+(3+(Nx-1)/2**kg)*(3+(Ny-1)/2**kg)
 30   CONTINUE
      
C ------------------------------------------
C     Storage for the coarse grids: NC
C ------------------------------------------

      NC=NF-(Nx+2)*(Ny+2)

C ----------------------------------------
C     Interpolation operator:
C ----------------------------------------

      IF (NC.EQ.0) THEN
         NCI=2
      ELSE
         NCI=8*NC
      ENDIF
      
C -----------------------------------------
C     Stencil:
C -----------------------------------------

      IF ( ISTNCL.EQ.1 ) THEN
         NSO = 3*NF + 2*NC 
      ELSE
         NSO = 5*NF
      ENDIF

C ----------------------------------------
C     Storage for the workspace SOR:
C ----------------------------------------


      IF ( ABS(IBC).EQ.0 .OR. ABS(IBC).EQ.4 ) THEN   !!! Temporary 

         IF ( iRELAX.EQ.BMG_GS_RB_point 
     &       .OR. iRELAX.EQ.BMG_GS_RB_x_lines ) THEN
            NSOR = 2*NF
         ELSE 
            NSOR = 4*NF
         ENDIF

      ELSE

         IF ( IRELAX.EQ.BMG_GS_RB_point 
     &      .OR. IRELAX.EQ.BMG_GS_RB_x_lines ) THEN
            NSOR = 2*NF
         ELSE
            NSOR = 4*NF
         ENDIF
         
      ENDIF

C ----------------------------------------
C     Solving on the coarsest grid:
C ----------------------------------------

      IF (ABS(IBC).EQ.0 .OR. ABS(IBC).EQ.4) THEN 
         NCBW = Nxc+2
      ELSE
         NCBW = Nxc*Nyc
      ENDIF
      NCU = Nxc*Nyc

C ========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_space',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END

