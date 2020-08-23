      SUBROUTINE BMG3_SymStd_SETUP_PtrGrid( 
     &                Nx, Ny, Nz, BMG_iPARMS, BMG_iWORK_PL,
     &                NOGm, NFm, NCm, NSOm, NSORm, NCIm, NCBWm, NCUm,
     &                NBMG_iWORK_PLm, NBMG_rWORK_PLm,
     &                NOG, NF, NC, NSO, NSOR, NCI, NCBW, NCU, IGRD,
     &                NBMG_iWORK_PL, NBMG_rWORK_PL, BMG_IOFLAG
     &                ) 

C =======================================================================
C 
C   BMG3_SymStd_SETUP_PtrGrid
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_SETUP_PtrGrid is used to compute the pointers into 
C   the various arrays for each grid within BMG3D.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/08/29 (JDM)
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
C   Nz         Number of points in z-direction (excluding ghost points)
C
C   -----------------------
C   Coarsening Criteria:
C   -----------------------
C
C   NXYZc      Minimum number of coarse grid points in either direction
C
C  ------------------------
C   Discretization:
C  ------------------------
C
C   iSTNCL     Discrete operator index:
C              iSTNCL .EQ. 1 => a 7  point discretization
C              iSTNCL .NE. 1 => a 27 point discretization
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
      INCLUDE    'BMG_parameters.h'
      INCLUDE    'BMG_constants.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NBMG_iWORK_PL, NBMG_iWORK_PLm,
     &         NBMG_rWORK_PL, NBMG_rWORK_PLm

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_iWORK_PL(NBMG_iWORK_PLm)

      INTEGER  NCBWm, NCIm, NCm, NCUm, NOGm, NFm, NSOm, NSORm

      INTEGER  iBC, IGRD(NOGm,NBMG_pIGRD), iRELAX, iSTNCL,
     &         NC, NCBW, NCI, NCU, NOG, NF, NSO, NSOR,
     &         Nx, Ny, Nz, NXYZc

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ---------------------------
C     Local Declarations:
C
      INTEGER  CF, IERR, i, j, kg, NCG_TEMP, Nxc, Nxg, NxNyNz,
     &         Nyc, Nyg, Nzc, Nzg, pSI, pSR


C ========================================================================

C =======================================================================
C --------------------------------------------
C     Unpack Parameters
C --------------------------------------------

      NXYZc   = BMG_iPARMS(id_BMG3_CG_MIN_DIM)
      IStncl  = BMG_iPARMS(id_BMG3_STENCIL)
      IRELAX  = BMG_iPARMS(id_BMG3_RELAX)
      IBC     = BMG_iPARMS(id_BMG3_BC)

C --------------------------------------------
C =======================================================================

C ========================================================================
C --------------------------------------------
C    Initialize *_PL Free space pointers
C --------------------------------------------

      BMG_iWORK_PL(ip_BMG_iWORK_PL_SF) = NBMG_iWORK_PL_ptrs + 1
      BMG_iWORK_PL(ip_BMG_rWORK_PL_SF) = 1

C --------------------------------------------
C    Compute the space requirements:
C --------------------------------------------

      CALL BMG3_SymStd_SETUP_space( 
     &                 Nx, Ny, Nz, BMG_iPARMS,
     &                 NOG, NF, NC, NSO, NSOR, NCI, NCBW, NCU,
     &                 NCG_TEMP, NBMG_iWORK_PL, NBMG_rWORK_PL,
     &                 BMG_IOFLAG
     &                 )

      IF ( BMG_iPARMS(id_BMG3_Err_Code).NE.iZERO ) THEN
         RETURN
      END IF
      


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
      !  Workspace for relaxation (e.g., LU factorization)
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
      !  The dimension of the plane relaxation INTEGER workspace 
      !  -> NBMG_iWORK_PL
      !
      IF ( NBMG_iWORK_PL.GT.NBMG_iWORK_PLm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Plane relaxtion INTEGER workspace is too small!'
         WRITE(*,510) 'HAVE: NBMG_iWORK_PLm = ', NBMG_iWORK_PLm
         WRITE(*,520) 'NEED: NBMG_iWORK_PL  = ', NBMG_iWORK_PL
      ENDIF

      !
      !  The dimension of the plane relaxation REAL workspace 
      !  -> NBMG_rWORK_PL
      !
      IF ( NBMG_rWORK_PL.GT.NBMG_rWORK_PLm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Plane relaxtion REAL workspace is too small!'
         WRITE(*,510) 'HAVE: NBMG_rWORK_PLm = ', NBMG_rWORK_PLm
         WRITE(*,520) 'NEED: NBMG_rWORK_PL  = ', NBMG_rWORK_PL
      ENDIF

      !
      !  If there were errors STOP
      !
      IF ( IERR.GT.0 ) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,530) IERR, 'FATAL ERRORS, EXECUTATION STOPPED'
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,9)
         RETURN
         
      ENDIF
      
C --------------------------------------------
C ========================================================================
C --------------------------------------------
C     Initialize the indexing array.
C --------------------------------------------

      DO j=1, NBMG_pIGRD
         DO i=1, NOGm
            IGRD(i,j) = iZERO
         END DO
      END DO

C -------------------------------------------
C     Finest Grid and First Coarse Grid:
C -------------------------------------------      

      Nxg    = Nx+2
      Nyg    = Ny+2
      Nzg    = Nz+2
      NxNyNz = Nxg*Nyg*Nzg

      !
      ! Dimension; finest grid
      !
      IGRD(NOG,idL_BMG_Nx) = Nxg
      IGRD(NOG,idL_BMG_Ny) = Nyg      
      IGRD(NOG,idL_BMG_Nz) = Nzg
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
      IGRD(NOG-1,ipL_BMG_U) = IGRD(NOG,ipL_BMG_U) + NxNyNz
      !
      ! Pointer: stencil
      !
      IF (iSTNCL.EQ.1) THEN
         IGRD(NOG-1,ipL_BMG_SO) = IGRD(NOG,ipL_BMG_SO)+4*NxNyNz
      ELSE
         IGRD(NOG-1,ipL_BMG_SO) = IGRD(NOG,ipL_BMG_SO)+14*NxNyNz
      ENDIF
      !
      ! Pointer: interpolation
      !
      IGRD(NOG-1,ipL_BMG_CI) = 1
      !
      ! Pointer: relaxation workspace
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         IGRD(NOG-1,ipL_BMG_SOR)
     &        = IGRD(NOG,ipL_BMG_SOR)+2*NxNyNz
      ELSE IF ( IRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN
         IGRD(NOG-1,ipL_BMG_SOR) 
     &        = IGRD(NOG,ipL_BMG_SOR)+2*NxNyNz
      ELSE IF ( IRELAX.EQ.BMG_GS_adaptive ) THEN
         IGRD(NOG-1,ipL_BMG_SOR) 
     &        = IGRD(NOG,ipL_BMG_SOR)+2*NxNyNz
      ELSE IF ( IRELAX.EQ.BMG_CG_relaxation ) THEN
         IGRD(NOG-1,ipL_BMG_SOR) 
     &        = IGRD(NOG,ipL_BMG_SOR)+4*NxNyNz
      ELSE
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The value of iRELAX is not supported ... '
            WRITE(*,520) '  HAVE: iRELAX = ', iRELAX
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,10)
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
         Nzg  = (Nz-1)/CF + 3
         NxNyNz = Nxg*Nyg*Nzg
         !
         ! Dimension: previous grid
         !
         IGRD(kg,idL_BMG_Nx) = Nxg
         IGRD(kg,idL_BMG_Ny) = Nyg
         IGRD(kg,idL_BMG_Nz) = Nzg
         !
         ! Pointer: vectors
         !
         IGRD(kg-1,ipL_BMG_U) = IGRD(kg,ipL_BMG_U)+NxNyNz
         !
         ! Pointer: stencil
         !
         IGRD(kg-1,ipL_BMG_SO) = IGRD(kg,ipL_BMG_SO)+14*NxNyNz
         !
         ! Pointer: interpolation
         !
         IGRD(kg-1,ipL_BMG_CI) = IGRD(kg,ipL_BMG_CI)+26*NxNyNz
         !
         ! Pointer: relaxation workspace
         !
         IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
            IGRD(kg-1,ipL_BMG_SOR) 
     &           = IGRD(kg,ipL_BMG_SOR)+2*NxNyNz
         ELSE IF ( IRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN
            IGRD(kg-1,ipL_BMG_SOR) 
     &           = IGRD(kg,ipL_BMG_SOR)+2*NxNyNz
         ELSE IF ( IRELAX.EQ.BMG_GS_adaptive ) THEN
            IGRD(kg-1,ipL_BMG_SOR) 
     &           = IGRD(kg,ipL_BMG_SOR)+2*NxNyNz
         ELSE IF ( IRELAX.EQ.BMG_CG_relaxation ) THEN
            IGRD(kg-1,ipL_BMG_SOR) 
     &           = IGRD(kg,ipL_BMG_SOR)+4*NxNyNz
         ELSE
            IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
               WRITE(*,500) 'The value of iRELAX is not supported ... '
               WRITE(*,520) '  HAVE: iRELAX = ', iRELAX
            END IF
            
            CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,10)
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
      Nzc = (Nz-1)/CF + 3
      !
      !  Dimension: coarsest grid
      !
      IGRD(1,idL_BMG_Nx) = Nxc
      IGRD(1,idL_BMG_Ny) = Nyc
      IGRD(1,idL_BMG_Nz) = Nzc

C -------------------------------------------------------
C     Compute pointers for operator construction:
C -------------------------------------------------------

      BMG_iWORK_PL(ip_BMG_CGTEMP_yo) = 1
      BMG_iWORK_PL(ip_BMG_CGTEMP_zo) =
     &     BMG_iWORK_PL(ip_BMG_CGTEMP_yo) + NCG_TEMP

C -------------------------------------------------------
C     Compute pointers for plane relaxation
C -------------------------------------------------------

      IF ( iRELAX.EQ.BMG_GS_adaptive  .OR. 
     &     iRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN
         !
         ! Initialize free-space pointer shifts
         !
         pSI = BMG_iWORK_PL(ip_BMG_iWORK_PL_SF)
         pSR = BMG_iWORK_PL(ip_BMG_rWORK_PL_SF)
         !
         ! BMG_iPARMS_PL_**
         !
         BMG_iWORK_PL(ip_BMG_iPARMS_PL_xy) = pSI
         pSI = BMG_iWORK_PL(ip_BMG_iPARMS_PL_xy) + NBMG_iPARMS*NOG
         !
         BMG_iWORK_PL(ip_BMG_iPARMS_PL_yz) = pSI
         pSI = BMG_iWORK_PL(ip_BMG_iPARMS_PL_yz) + NBMG_iPARMS*NOG
         !
         BMG_iWORK_PL(ip_BMG_iPARMS_PL_xz) = pSI
         pSI = BMG_iWORK_PL(ip_BMG_iPARMS_PL_xz) + NBMG_iPARMS*NOG
         
         !
         ! BMG_pWORK_PL_**
         !
         BMG_iWORK_PL(ip_BMG_pWORK_PL_xy) = pSI
         pSI = BMG_iWORK_PL(ip_BMG_pWORK_PL_xy) 
     &       + NBMG_pWORK*NOG*(Nz+2)

         BMG_iWORK_PL(ip_BMG_pWORK_PL_yz) = pSI
         pSI = BMG_iWORK_PL(ip_BMG_pWORK_PL_yz) 
     &       + NBMG_pWORK*NOG*(Nx+2)

         BMG_iWORK_PL(ip_BMG_pWORK_PL_xz) = pSI
         pSI = BMG_iWORK_PL(ip_BMG_pWORK_PL_xz)
     &       + NBMG_pWORK*NOG*(Ny+2)

         !
         ! Update the free space counters
         !
         BMG_iWORK_PL(ip_BMG_iWORK_PL_SF) = pSI
         BMG_iWORK_PL(ip_BMG_rWORK_PL_SF) = pSR

      ENDIF

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_PtrGrid.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
