      SUBROUTINE BMG3_SymStd_SETUP_PtrWork( 
     &                Nx, Ny, Nz, BMG_iPARMS, 
     &                NOGm, NFm, NSOm, 
     &                NBMG_iWORKm, NBMG_rWORKm,
     &                NBMG_iWORK_PLm, NBMG_rWORK_PLm,
     &                BMG_pWORK, BMG_InWORK, pSR, pSI,
     &                BMG_IOFLAG
     &                )

C =======================================================================
C 
C   BMG3_SymStd_SETUP_PtrWork
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_SETUP_PtrWork is used to compute the pointers into 
C   the user's work arrays for a direct call to BMG3D.
C
C   Written:    1987/03/05  (JED)
C   - original name: STHMG3
C   Modified:   1988/01/21  (JED)
C   - details unknown
C   Modified:   1999/11/18  (JDM)
C   - new name: boxspace3
C   - renamed variables based on their use 
C   - shifted loop bounds in calculation of NCPL2 (formerly NABD32)
C   - renamed subroutine to BOXspace3
C   - defined a set of pointers for the call to BMG3
C   - added documentation
C   Modified:   2000/01/23  (JDM)
C   - added checks on NSOm, NFm, and NOGm
C   - added option to include SO, Q, and QF in the workspace
C   Modified:   2000/08/29  (JDM)
C   - new name: BMG3_SymStd_SETUP_PtrWork
C   - changelog maintained in CVS look there for details.
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
C   NX       Number of points in x-direction (excluding ghost points)
C   NY       Number of points in y-direction (excluding ghost points)
C   NZ       Number of points in z-direction (excluding ghost points)
C
C  ------------------------
C   Discretization:
C  ------------------------
C
C   IFD      Discrete operator index:
C            iFD .EQ. 1 => a 7  point discretization
C            iFD .NE. 1 => a 27 point discretization
C
C   ------------------------
C   Smoothing:
C   ------------------------   
C
C   IRELAX   Relaxation index (refer to BMG3D).
C
C   --------------------------
C   Work
C   --------------------------
C
C   IWORK     Integer work array (4*NOG)
C
C   -------------------------
C   Multigrid:
C   -------------------------
C
C   NOGm      Maximum number of grids.
C   NSOm      Maximum allowed space for 3D stencils (all grids).
C   NFm       Maximum allowed space for grid functions (all grids).
C
C   * note: NSOm and NFm may be zero if these arrays are to be stored
C           in the real work array RWORK.
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
C =======================================================================
C   LOCAL:
C ========================
C
C   Nxc       Number of points in the x-direction on a coarser grid
C   Nyc       Number of points in the y-direction on a coarser grid
C   Nzc       Number of points in the z-direction on a coarser grid
C
C =======================================================================

      IMPLICIT   NONE

      INCLUDE    'BMG_workspace.h'
      INCLUDE    'BMG_constants.h'
      INCLUDE    'BMG_parameters.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NFm, NBMG_iWORKm, NBMG_iWORK_PLm, NOGm,
     &         NBMG_rWORKm, NBMG_rWORK_PLm, NSOm

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_pWORK(NBMG_pWORK), 
     &         NC, NCBW, NCU, NCI, NOG, NF, NSO,
     &         NSOR, Nx, Ny, Nz, pSI, pSR

      LOGICAL  BMG_InWORK(NBMG_InWORK)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ---------------------------
C     Local Declarations:
C
      INTEGER  NCG_TEMP, NIGRD3, NBMG_iWORK, NBMG_iWORK_PL,
     &         NBMG_rWORK, NBMG_rWORK_PL, pSR_S, pSI_S
      LOGICAL  FAIL_NOG, FAIL_Q, FAIL_RES, FAIL_SO, FAIL_U,
     &         FAIL_WORK, FAIL_WORK_PL

C ========================================================================

C ========================================================================
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
C     Store the space requirements:
C --------------------------------------------

      BMG_iPARMS(id_BMG3_DIM_NOG)  = NOG
      BMG_iPARMS(id_BMG3_DIM_NF)   = NF
      BMG_iPARMS(id_BMG3_DIM_NC)   = NC
      BMG_iPARMS(id_BMG3_DIM_NCI)  = NCI
      BMG_iPARMS(id_BMG3_DIM_NSO)  = NSO
      BMG_iPARMS(id_BMG3_DIM_NSOR) = NSOR
      BMG_iPARMS(id_BMG3_DIM_NCBW) = NCBW
      BMG_iPARMS(id_BMG3_DIM_NCU)  = NCU

C --------------------------------------------
C    Grid Indexing and pointers:
C --------------------------------------------

      NIGRD3 = NBMG_pIGRD*NOG

C --------------------------------------------

C ========================================================================

C ------------------------------------------------------
C     Compute the pointers into the real work space:
C ------------------------------------------------------

      pSR_S=pSR
      pSI_S=pSI

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
      IF ( BMG_InWORK(i_InWORK_PCG_P) ) THEN
         BMG_pWORK(ip_BMG_PCG_P) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_P) + NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_P) = iZERO
      ENDIF


      !
      ! Include the the PCG workspace vector R in RWORK
      !
      IF ( BMG_InWORK(i_InWORK_PCG_R) ) THEN
         BMG_pWORK(ip_BMG_PCG_R) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_R) +  NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_R) = iZERO
      ENDIF

      !
      ! Include the the PCG workspace vector Z in RWORK
      !
      IF ( BMG_InWORK(i_InWORK_PCG_Z) ) THEN
         BMG_pWORK(ip_BMG_PCG_Z) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_Z) + NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_Z) = iZERO
      ENDIF

      BMG_pWORK(ip_SOR) = pSR
      BMG_pWORK(ip_CI)  = BMG_pWORK(ip_SOR) + NSOR
      BMG_pWORK(ip_CSO) = BMG_pWORK(ip_CI) + NCI
      BMG_pWORK(ip_CU)  = BMG_pWORK(ip_CSO) + NCBW*NCU

      pSR = BMG_pWORK(ip_CU) + NCU
      NBMG_rWORK = pSR - pSR_S

C -------------------------------------------------------
C     Compute the pointers into the integer work space:
C -------------------------------------------------------

      BMG_pWORK(ip_iG) = pSI

      pSI = pSI + NIGRD3
      NBMG_iWORK = pSI - pSI_S

C ========================================================================

C ========================================================================
C --------------------------------------------------------
C     Memory Allocation:
C --------------------------------------------------------   

      IF ( BMG_iPARMS(id_BMG3_POINTERS)
     &    .EQ.BMG_NO_pointers ) THEN

         !
         !  Verify allocation for standard f77-style workspace
         !

         FAIL_NOG     = (NOG.GT.NOGm)
         FAIL_SO      = (.NOT.BMG_InWORK(i_InWORK_SO)).AND.(NSO.GT.NSOm)
         FAIL_U       = (.NOT.BMG_InWORK(i_InWORK_U)).AND.(NF.GT.NFm)
         FAIL_Q       = (.NOT.BMG_InWORK(i_InWORK_Q)).AND.(NF.GT.NFm)
         FAIL_RES     = (.NOT.BMG_InWORK(i_InWORK_RES)).AND.(NF.GT.NFm)
         FAIL_WORK    = (NBMG_rWORK.GT.NBMG_rWORKm)
     &                  .OR. (NBMG_iWORK.GT.NBMG_iWORKm)
         FAIL_WORK_PL = (NBMG_rWORK_PL.GT.NBMG_rWORK_PLm)
     &                  .OR. (NBMG_iWORK_PL.GT.NBMG_iWORK_PLm)
         
         IF ( FAIL_NOG .OR. FAIL_SO .OR. FAIL_U .OR. FAIL_Q 
     &       .OR. FAIL_RES .OR. FAIL_WORK .OR. FAIL_WORK_PL ) THEN 
           IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN

            WRITE(*,*) 
            WRITE(*,*) '*** ERROR: MEMORY ALLOCATION OUT OF BOUNDS ***'
            WRITE(*,*) 
         
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
               WRITE(*,*) 'Insufficient Workspace allocated: 3D !!'
               WRITE(*,*) 
               WRITE(*,*) ' - required REAL allocation = ',
     &                    NBMG_rWORK
               WRITE(*,*) ' - actual REAL allocation   = ',
     &                    NBMG_rWORKm
               WRITE(*,*) 
               WRITE(*,*) ' - required INTEGER allocation = ',
     &                    NBMG_iWORK
               WRITE(*,*) ' - actual INTEGER allocation   = ',
     &                    NBMG_iWORKm
               WRITE(*,*)
            ENDIF

            IF ( FAIL_WORK_PL ) THEN
               WRITE(*,*) 'Insufficient workspace allocated: planes !!'
               WRITE(*,*) 
               WRITE(*,*) ' - required REAL allocation = ',
     &                    NBMG_rWORK_PL
               WRITE(*,*) ' - actual REAL allocation   = ', 
     &                    NBMG_rWORK_PLm
               WRITE(*,*) 
               WRITE(*,*) ' - required INTEGER allocation = ',
     &                    NBMG_iWORK_PL
               WRITE(*,*) ' - actual INTEGER allocation   = ',
     &                    NBMG_iWORK_PLm
               WRITE(*,*)
            ENDIF


            WRITE(*,*) ' *** CHECK WORKSPACE POINTERS *** '
            WRITE(*,*)
            WRITE(*,*) 'BMG_pWORK(ip_SO)        ', BMG_pWORK(ip_SO)
            WRITE(*,*) 'BMG_pWORK(ip_U)         ', BMG_pWORK(ip_U)
            WRITE(*,*) 'BMG_pWORK(ip_Q)         ', BMG_pWORK(ip_Q)
            WRITE(*,*) 'BMG_pWORK(ip_RES)       ', BMG_pWORK(ip_RES)
            WRITE(*,*) 'BMG_pWORK(ip_SOR)       ', BMG_pWORK(ip_SOR)
            WRITE(*,*) 'BMG_pWORK(ip_CI)        ', BMG_pWORK(ip_CI)
            WRITE(*,*) 'BMG_pWORK(ip_CSO)       ', BMG_pWORK(ip_CSO)
            WRITE(*,*) 'BMG_pWORK(ip_CU)        ', BMG_pWORK(ip_CU)
            WRITE(*,*) 'BMG_pWORK(ip_BMG_PCG_P) ',
     &                 BMG_pWORK(ip_BMG_PCG_P)
            WRITE(*,*) 'BMG_pWORK(ip_BMG_PCG_R) ',
     &                 BMG_pWORK(ip_BMG_PCG_R)
            WRITE(*,*) 'BMG_pWORK(ip_BMG_PCG_Z) ',
     &                 BMG_pWORK(ip_BMG_PCG_Z)
            WRITE(*,*) 
            WRITE(*,*) 'The required REAL allocation = ', NBMG_rWORK
            WRITE(*,*) 'Actual REAL Allocation       = ',
     &                 NBMG_rWORKm
            WRITE(*,*) 
            WRITE(*,*) 'BMG_pWORK(ip_iG)    ', BMG_pWORK(ip_iG)
            WRITE(*,*)
            WRITE(*,*) 'The required INTEGER Allocation  = ', NBMG_iWORK
            WRITE(*,*) 'Actual INTEGER Allocation        = ',
     &                 NBMG_iWORKm
            WRITE(*,*)
            WRITE(*,*) 

            WRITE(*,*) 'The required REAL allocation for planes = ',
     &                 NBMG_rWORK_PL
            WRITE(*,*) 'Actual REAL allocation for planes = ',
     &                 NBMG_rWORK_PLm
            WRITE(*,*) 
            WRITE(*,*) 'The required INTEGER Allocation for planes  = ',
     &                 NBMG_iWORK_PL
            WRITE(*,*) 'Actual INTEGER Allocation        = ',
     &                 NBMG_iWORK_PLm
            WRITE(*,*)
            WRITE(*,*) 

            WRITE(*,101) '****************************************'
            WRITE(*,102) 'CHECK PARAMETER VALUES'
            WRITE(*,103) '****************************************'
            
            WRITE(*,100) 'Vector on all grids, NF', NF 
            WRITE(*,100) 'Vector on coarse grids, NC', NC
            WRITE(*,100) 'Discrete operator on all grids, NSO', NSO
            WRITE(*,100) 'Workspace, at least, NSOR = 2*NF', NSOR
            WRITE(*,100) 'Interpolation operator coarse grids, NCI', NCI
            
            WRITE(*,100) 'Vector on the coarsest grid, NCU', NCU
            WRITE(*,100) 'Coarsest grid operator bandwidth, NCBW', NCBW
            
            WRITE(*,104) 
           END IF

           CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,16)
           RETURN

         ENDIF

      ELSEIF ( BMG_iPARMS(id_BMG3_POINTERS)
     &        .EQ. BMG_USE_pointers      ) THEN

         !
         !  Return workspace dimensions for CRAY-style pointers/malloc
         !

         NBMG_iWORKm = NBMG_iWORK
         NBMG_rWORKm = NBMG_rWORK
         
         NBMG_iWORK_PLm = NBMG_iWORK_PL
         NBMG_rWORK_PLm = NBMG_rWORK_PL

      ELSE
        IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN

         
         WRITE(*,*) 
         WRITE(*,*) '*** ERROR: MEMORY ALLOCATION MODE UNSPECIFIED ***'
         WRITE(*,*) 
         WRITE(*,*) 'BMG_iPARMS(id_BMG3_POINTERS) = ', 
     &              BMG_iPARMS(id_BMG3_POINTERS)
         WRITE(*,*)
         WRITE(*,*) 'Admissible values defined in BMG_parameters.h:'
         WRITE(*,*) 
         WRITE(*,*) 'BMG_USE_pointers = ', BMG_USE_pointers
         WRITE(*,*) 'BMG_NO_pointers  = ', BMG_NO_pointers
        END IF
        
        CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,17)
        RETURN

      ENDIF

C -------------------------------------------------------
C =========================================================================

 100  FORMAT (8X,A,I8)
      
 101  FORMAT (/,2X,A)
 102  FORMAT (10X,A)
 103  FORMAT (2X,A,/)
 104  FORMAT (/)

C ============================


      RETURN
      END


