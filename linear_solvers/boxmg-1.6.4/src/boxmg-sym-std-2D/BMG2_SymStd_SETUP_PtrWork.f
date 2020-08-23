      SUBROUTINE BMG2_SymStd_SETUP_PtrWork( 
     &                Nx, Ny, BMG_iPARMS, 
     &                NOGm, NFm, NSOm, NIWORKm, NRWORKm, 
     &                BMG_pWORK, BMG_InWORK, pSR, pSI,
     &                BMG_IOFLAG
     &                )

C =======================================================================
C 
C   BMG2_SymStd_SETUP_PtrWork
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_PtrWork is used to compute the pointers into 
C   the user's work arrays for a direct call to BOXMG.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    1994/11/30 (JDM)
C   Modified:   1999/11/22 (JDM)
C   - defined a pointer array for the call to BOXMG
C   - added documentation
C   Modified:  2000/02/14 (JDM)
C   - added checks on NSOm, NFm, and NOGm
C   - added option to include SO, Q, and QF in the workspace.
C   Modified:  2000/07/24 (JDM)
C   - renamed boxspace.f -> BMG2_SymStd_SETUP_PtrWork in preparation
C     for rewriting the space allocation parts of boxmg.
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
C   -----------------------
C   Coarsening Criteria:
C   -----------------------
C
C   NXYc      Minimum number of coarse grid points in either direction
C
C  ------------------------
C   Discretization:
C  ------------------------
C
C   iSTNCL   Discrete operator index:
C            iSTNCL .EQ. 1 => a 5  point discretization
C            iSTNCL .NE. 1 => a 9  point discretization
C
C   ------------------------
C   Smoothing:
C   ------------------------   
C
C   iRELAX   Relaxation index (refer to BOXMG).
C
C   ------------------------
C   Boundary Conditions:
C   ------------------------
C
C   iBC      Boundary condition index (refer to BOXMG,BOXMGP)
C
C   ------------------------
C   Dimensions:
c   ------------------------
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
C   ------------------------
C   Workspace:
C   ------------------------
C
C   InWork    Logical array that controls which arrays are contained in rWORK
C
C   NInWork   Dimension of InWork
C
C   NBMG_pWORK    Dimension of pointer array BMG_pWORK
C
C =======================================================================
C   OUTPUT:
C ===========================
C
C   BMG_pWORK     Integer array of (fake) pointers to arrays in rWORK
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
C   NIGRD     INTEGER workspace for grid indexing
C
C =======================================================================

      IMPLICIT   NONE

      INCLUDE    'BMG_workspace.h'
      INCLUDE    'BMG_constants.h'
      INCLUDE    'BMG_parameters.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NFm, NIWORKm, NOGm, NRWORKm, NSOm

      INTEGER  BMG_iPARMS(NBMG_iPARMS), NC, NCBW, NCU, NCI, NOG, NF,
     &         NSO, NSOR, Nx, Ny, pSI, pSR, BMG_pWORK(NBMG_pWORK)
      LOGICAL  BMG_InWORK(NBMG_InWORK)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ---------------------------
C     Local Declarations:
C
      INTEGER  iBC, iRELAX, iSTNCL, NIGRD, NIWORK, NRWORK, NXYc,
     &         pSI_S, pSR_S
      LOGICAL  FAIL_NOG, FAIL_Q, FAIL_RES, FAIL_SO, FAIL_U, FAIL_WORK

C ========================================================================

C ========================================================================
C --------------------------------------------
C     Unpack Parameters:
C --------------------------------------------
      
      NXYc   = BMG_iPARMS(id_BMG2_CG_MIN_DIM)
      IStncl = BMG_iPARMS(id_BMG2_STENCIL)
      IBC    = BMG_iPARMS(id_BMG2_BC)
      IRELAX = BMG_iPARMS(id_BMG2_RELAX)

C ========================================================================

C ========================================================================
C --------------------------------------------
C    Compute the space requirements:
C --------------------------------------------
      
      CALL BMG2_SymStd_SETUP_space( 
     &          Nx, Ny, BMG_iPARMS,
     &          NOG, NF, NC, NCI, NSO, NSOR, NCBW, NCU
     &          )

C --------------------------------------------
C     Store the space requirements:
C --------------------------------------------

      BMG_iPARMS(id_BMG2_DIM_NOG)  = NOG
      BMG_iPARMS(id_BMG2_DIM_NF)   = NF
      BMG_iPARMS(id_BMG2_DIM_NC)   = NC
      BMG_iPARMS(id_BMG2_DIM_NCI)  = NCI
      BMG_iPARMS(id_BMG2_DIM_NSO)  = NSO
      BMG_iPARMS(id_BMG2_DIM_NSOR) = NSOR
      BMG_iPARMS(id_BMG2_DIM_NCBW) = NCBW
      BMG_iPARMS(id_BMG2_DIM_NCU)  = NCU

C --------------------------------------------
C    Grid Indexing and pointers:
C --------------------------------------------

      NIGRD = 9*NOG

C --------------------------------------------
C ========================================================================


C ========================================================================
C ------------------------------------------------------
C     Compute the pointers into the real work space:
C ------------------------------------------------------

      pSI_S = pSI
      pSR_S = pSR

      !
      ! Include the stencil, SO, in RWORK
      !
      IF ( BMG_InWORK(i_InWORK_SO) ) THEN 
         BMG_pWORK(ip_SO) = pSR
         pSR = BMG_pWORK(ip_SO) + NSO
      ELSE
         BMG_pWORK(ip_SO) = iZERO
      ENDIF

      !
      ! Include the solution, U, in RWORK
      !      
      IF ( BMG_InWORK(i_InWORK_U) ) THEN
         BMG_pWORK(ip_U) = pSR
         pSR = BMG_pWORK(ip_U) + NF
      ELSE
         BMG_pWORK(ip_U) = iZERO
      ENDIF

      !
      ! Include the right hand side, QF, in RWORK
      !
      IF ( BMG_InWORK(i_InWORK_Q) ) THEN
         BMG_pWORK(ip_Q) = pSR
         pSR = BMG_pWORK(ip_Q) + NF
      ELSE 
         BMG_pWORK(ip_Q) = iZERO
      ENDIF

      !
      ! Include the residual, RES, in RWORK
      !
      IF ( BMG_InWORK(i_InWORK_RES) ) THEN
         BMG_pWORK(ip_RES) = pSR
         pSR = BMG_pWORK(ip_RES) + NF
      ELSE 
         BMG_pWORK(ip_RES) = iZERO
      ENDIF

      !
      ! Include the the PCG workspace vector P in RWORK
      !
      IF (BMG_InWORK(i_InWORK_PCG_P)) THEN
         BMG_pWORK(ip_BMG_PCG_P) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_P) + NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_P) = iZERO
      ENDIF


      !
      ! Include the the PCG workspace vector R in RWORK
      !
      IF (BMG_InWORK(i_InWORK_PCG_R)) THEN
         BMG_pWORK(ip_BMG_PCG_R) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_R) +  NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_R) = iZERO
      ENDIF

      !
      ! Include the the PCG workspace vector Z in RWORK
      !
      IF (BMG_InWORK(i_InWORK_PCG_Z)) THEN
         BMG_pWORK(ip_BMG_PCG_Z) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_Z) + NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_Z) = iZERO
      ENDIF

      !
      BMG_pWORK(ip_SOR) = pSR
      BMG_pWORK(ip_CI)  = BMG_pWORK(ip_SOR) + NSOR
      BMG_pWORK(ip_CSO) = BMG_pWORK(ip_CI) + NCI
      BMG_pWORK(ip_CU)  = BMG_pWORK(ip_CSO) + NCBW*NCU

      psR = BMG_pWORK(ip_CU) + NCU
      NRWORK = pSR - pSR_S

C -------------------------------------------------------
C     Compute the pointers into the integer work space:
C -------------------------------------------------------

      BMG_pWORK(ip_iG) = pSI

      pSI = BMG_pWORK(ip_iG) + NIGRD
      NIWORK = pSI - pSI_S

C --------------------------------------------------------
C     Memory Allocation:
C --------------------------------------------------------   

      IF ( BMG_iPARMS(id_BMG2_POINTERS)
     &    .EQ.BMG_NO_pointers ) THEN   

         !
         !  Verify allocation for standard f77-style workspace
         !

         FAIL_NOG   = (NOG.GT.NOGm)
         FAIL_SO    = (.NOT.BMG_InWORK(i_InWORK_SO)).AND.(NSO.GT.NSOm)
         FAIL_U     = (.NOT.BMG_InWORK(i_InWORK_U)).AND.(NF.GT.NFm)
         FAIL_Q     = (.NOT.BMG_InWORK(i_InWORK_Q)).AND.(NF.GT.NFm)
         FAIL_RES   = (.NOT.BMG_InWORK(i_InWORK_RES)).AND.(NF.GT.NFm)
         FAIL_WORK  = (pSR-1.GT.NRWORKm).OR.(pSI-1.GT.NIWORKm)

         IF ( FAIL_NOG .OR. FAIL_SO .OR. FAIL_U .OR. FAIL_Q 
     &        .OR. FAIL_RES .OR. FAIL_WORK ) THEN 

           IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN

            WRITE(*,*) 
            WRITE(*,*) '*** ERROR: MEMORY ALLOCATION OUT OF BOUNDS ***'
            WRITE(*,*) 
               
            IF ( FAIL_NOG ) THEN 
               WRITE(*,*) 'There are not enough grids !!!! '
               WRITE(*,*) '  - required number of grids, NOG = ', NOG
               WRITE(*,*) '  - maximum number of grids, NOGm = ', NOGm
            ENDIF

            IF ( FAIL_SO ) THEN
               WRITE(*,*) 'The memory allocation for SO failed !!'
               WRITE(*,*) '  - required storage, NSO  = ', NSO
               WRITE(*,*) '  - maximum storage,  NSOm = ', NSOm
               WRITE(*,*) 
            ENDIF

            IF ( FAIL_U ) THEN
               WRITE(*,*) 'The memory allocation for U failed !!'
               WRITE(*,*) '  - required storage, NF  = ', NF
               WRITE(*,*) '  - maximum storage,  NFm = ', NFm
               WRITE(*,*) 
            ENDIF

            IF ( FAIL_Q ) THEN
               WRITE(*,*) 'The memory allocation for Q failed !!'
               WRITE(*,*) '  - required storage, NF  = ', NF
               WRITE(*,*) '  - maximum storage,  NFm = ', NFm
               WRITE(*,*) 
            ENDIF

            IF ( FAIL_RES ) THEN
               WRITE(*,*) 'The memory allocation for RES failed !!'
               WRITE(*,*) '  - required storage, NF  = ', NF
               WRITE(*,*) '  - maximum storage,  NFm = ', NFm
               WRITE(*,*) 
            ENDIF

            IF ( FAIL_WORK ) THEN
               WRITE(*,*) 'Workspace allocation failed !!'
               WRITE(*,*) 
               WRITE(*,*) ' - required REAL allocation, pSR-1  = ',
     &                    pSR - 1
               WRITE(*,*) ' - required REAL allocation, NRWORK = ',
     &                    NRWORK
               WRITE(*,*) ' - actual allocation, NRWORK = ', NRWORKm
     &                    
               WRITE(*,*) 
               WRITE(*,*) ' - required INTEGER allocation, pSI-1  = ',
     &                    pSI - 1
               WRITE(*,*) ' - required INTEGER allocation, NIWORK = ',
     &                    NIWORK
               WRITE(*,*) ' - actual allocation, NIWORK   = ', NIWORKm
               WRITE(*,*)
            ENDIF

            WRITE(*,*) ' *** CHECK WORKSPACE POINTERS *** '
            WRITE(*,*)
            WRITE(*,*) 'Initial pSI =   ', pSI_S
            WRITE(*,*) 'Initial pSR =   ', pSR_S
            WRITE(*,*)
            WRITE(*,*) 'BMG_pWORK(ip_SO)           ', BMG_pWORK(ip_SO)
            WRITE(*,*) 'BMG_pWORK(ip_U)            ', BMG_pWORK(ip_U)
            WRITE(*,*) 'BMG_pWORK(ip_Q)            ', BMG_pWORK(ip_Q)
            WRITE(*,*) 'BMG_pWORK(ip_RES)          ', BMG_pWORK(ip_RES)
            WRITE(*,*) 'BMG_pWORK(ip_SOR)          ', BMG_pWORK(ip_SOR)
            WRITE(*,*) 'BMG_pWORK(ip_CI)           ', BMG_pWORK(ip_CI)
            WRITE(*,*) 'BMG_pWORK(ip_CSO)          ', BMG_pWORK(ip_CSO)
            WRITE(*,*) 'BMG_pWORK(ip_CU)           ', BMG_pWORK(ip_CU)
            WRITE(*,*) 'BMG_pWORK(ip_BMG_PCG_P)    ',
     &                 BMG_pWORK(ip_BMG_PCG_P)
            WRITE(*,*) 'BMG_pWORK(ip_BMG_PCG_R)    ',
     &                 BMG_pWORK(ip_BMG_PCG_R)
            WRITE(*,*) 'BMG_pWORK(ip_BMG_PCG_Z)    ',
     &                 BMG_pWORK(ip_BMG_PCG_Z)
            WRITE(*,*) 
            WRITE(*,*) 'Required REAL allocation = ', NRWORK
            WRITE(*,*) 'Actual REAL Allocation   = ', NRWORKm
            WRITE(*,*) 
            WRITE(*,*) 'BMG_pWORK(ip_iG)    ', BMG_pWORK(ip_iG)
            WRITE(*,*)
            WRITE(*,*) 'Required INTEGER Allocation  = ', NIWORK
            WRITE(*,*) 'Actual INTEGER Allocation    = ', NIWORKm
            WRITE(*,*)
            WRITE(*,*) 

            WRITE(*,101) '****************************************'
            WRITE(*,102) 'CHECK PARAMETER VALUES'
            WRITE(*,103) '****************************************'
            
            WRITE(*,100) 'NF', NF 
            WRITE(*,100) 'NC', NC
            WRITE(*,100) 'NSO', NSO
            WRITE(*,100) 'NSOR', NSOR
            WRITE(*,100) 'NCI', NCI
            
            WRITE(*,100) 'NCU', NCU
            WRITE(*,100) 'NCBW', NCBW
           END IF

           CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,14)
           RETURN

         ENDIF

      ELSEIF ( BMG_iPARMS(id_BMG2_POINTERS)
     &        .EQ. BMG_USE_pointers      ) THEN

         !
         !  Return workspace dimensions for CRAY-style pointers/malloc
         !

         NIWORKm = NIWORK
         NRWORKm = NRWORK

      ELSE

        IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
         
         WRITE(*,*) 
         WRITE(*,*) '*** ERROR: MEMORY ALLOCATION MODE UNSPECIFIED ***'
         WRITE(*,*) 
         WRITE(*,*) 'BMG_iPARMS(id_BMG2_POINTERS) = ', 
     &              BMG_iPARMS(id_BMG2_POINTERS)
         WRITE(*,*)
         WRITE(*,*) 'Admissible values are: '
         WRITE(*,*) 
         WRITE(*,*) 'BMG_USE_pointers = ', BMG_USE_pointers
         WRITE(*,*) 'BMG_NO_pointers  = ', BMG_NO_pointers

        END IF
      
        CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,15)
        RETURN


      ENDIF

C -------------------------------------------------------
C =========================================================================

 100  FORMAT (8X,A,T24,I10)
      
 101  FORMAT (/,2X,A)
 102  FORMAT (10X,A)
 103  FORMAT (2X,A,/)
 104  FORMAT (/)

C ---------------------------

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_PtrWork.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END


