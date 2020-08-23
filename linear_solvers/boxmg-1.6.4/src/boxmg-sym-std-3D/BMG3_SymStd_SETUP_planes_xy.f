      SUBROUTINE BMG3_SymStd_SETUP_planes_xy( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, IFD, IBC, IGRD, 
     &                SO, NSO, NOGm, NOG,
     &                BMG_iPARMS_PL_xy, BMG_pWORK_PL_xy,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

C ==========================================================================
C
C   BMG3_SymStd_SETUP_planes_xy.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_SETUP_planes_xy.f performs any necessary setup for 
C   (x,y) plane relaxation.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2003/02/03 (JDM)
C
C ==================================================================
C   INPUT:
C ========================
C
C
C
C ==================================================================
C   OUTPUT:
C ===========================
C
C
C
C ==================================================================
C   LOCAL:
C ========================
C
C
C
C ==========================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER  NSO, NBMG_iWORK_PL, NBMG_rWORK_PL,
     &         NOG, NOGm, Nxf, Nyf, Nzf

      INTEGER  BMG_iPARMS(NBMG_iPARMS),
     &         BMG_iPARMS_PL_xy(NBMG_iPARMS,NOG),
     &         BMG_iWORK_PL(NBMG_iWORK_PL),
     &         BMG_pWORK_PL_xy(NBMG_pWORK,Nzf,NOG),
     &         IBC, IFD, IGRD(NOGm,NBMG_pIGRD)
      REAL*8   BMG_rPARMS(NBMG_rPARMS),
     &         BMG_rWORK_PL(NBMG_rWORK_PL),
     &         SO(NSO)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ------------------------------------------------
C     Local Declarations
C
      INTEGER  i, id, iPL, kg, Nx, Ny, Nz, p_U, p_SO, p_SOR, p_CI
      LOGICAL  BMG_InWORK_PL(NBMG_InWORK)

      INTEGER  NFm_xy, NOGm_xy, NSOm_xy,
     &         NBMG_iWORK_PL_temp, NBMG_rWORK_PL_temp, 
     &         pSI, pSR

      INTEGER  KNOG, NC_xy, NCBW_xy, NCI_xy, NCU_xy, NF_xy, NOG_xy,
     &         NSO_xy, NSOR_xy, NStncl_3D, NStncl_2D 

      INTEGER  p_CI_xy, p_CSO_xy, p_CU_xy, p_IG_xy, p_SO_xy, p_SOR_xy,
     &         p_Q_xy, p_RES_xy, p_U_xy

C ======================================================================


C ======================================================================
    

      !
      !  Copy BMG_iPARMS for this level
      !  (necessary because of grid dimension info)
      !         
      DO id = 1, NBMG_iPARMS
         BMG_iPARMS_PL_xy(id,NOG) = BMG_iPARMS(id)
      ENDDO

      !
      !  Compute maximum dimensions for (x,y) - planes
      !
      BMG_iPARMS_PL_xy(id_BMG2_POINTERS,NOG)=BMG_USE_pointers
      !
      NOGm_xy = 1
      NFm_xy  = 1
      NSOm_xy = 1
      !
      pSR = 1
      pSI = 1
      !
      ! Set error codes to zero
      !
      BMG_iPARMS_PL_xy(id_BMG2_Err_Code, NOG) = 0
      BMG_iPARMS_PL_xy(id_BMG2_Ext_Err_Code, NOG) = 0
      !
      CALL BMG2_SymStd_SETUP_PtrWork( 
     &          Nxf-2, Nyf-2, BMG_iPARMS_PL_xy(1,NOG),
     &          NOGm_xy, NFm_xy, NSOm_xy,
     &          NBMG_iWORK_PL_temp,
     &          NBMG_rWORK_PL_temp,
     &          BMG_pWORK_PL_xy(1,1,NOG), 
     &          BMG_InWORK_PL, pSR, pSI,
     &          BMG_IOFLAG
     &          ) 

      !
      ! Store maximum dimensions for (x,y) - planes
      ! 
      NOGm_xy = BMG_iPARMS_PL_xy(id_BMG2_DIM_NOG,NOG)
      NFm_xy  = BMG_iPARMS_PL_xy(id_BMG2_DIM_NF,NOG)
      NSOm_xy = BMG_iPARMS_PL_xy(id_BMG2_DIM_NSO,NOG)
      BMG_iPARMS(id_BMG2_POINTERS)=BMG_NO_pointers

C         WRITE(*,*) 'Nxf, Nyf  = ', Nxf-2, Nyf-2
C         WRITE(*,*) 'NOGm_xy   = ', NOGm_xy
C         WRITE(*,*) 'NFm_xy    = ', NFm_xy
C         WRITE(*,*) 'NSOm_xy   = ', NSOm_xy

      !
      !   Set the starting point of free workspace
      !
      pSI = BMG_iWORK_PL(ip_BMG_iWORK_PL_SF)
      pSR = BMG_iWORK_PL(ip_BMG_rWORK_PL_SF)

C         WRITE(*,*) 'pSI = ', pSI
C         WRITE(*,*) 'pSR = ', pSR


      !
      !  Initialize BMG_InWORK_PL
      !
      DO i=1, NBMG_InWORK
         BMG_InWORK_PL(i) = .FALSE.
      ENDDO

      !
      !   Override the BoxMG arguments
      !
      BMG_InWORK_PL(i_InWORK_SO)  = .TRUE. 
      BMG_InWORK_PL(i_InWORK_U)   = .TRUE. 
      BMG_InWORK_PL(i_InWORK_Q)   = .TRUE.
      BMG_InWORK_PL(i_InWORK_RES) = .TRUE.

      !
      !   Loop of grid-levels
      !
      DO kg = NOG, 2, -1

         !
         !  Copy BMG_iPARMS for this level
         !  (necessary because of grid dimension info)
         !         
         DO id = 1, NBMG_iPARMS
            BMG_iPARMS_PL_xy(id,kg) = BMG_iPARMS(id)
         ENDDO

         !
         !  Modifications to relaxation parameters may
         !  be necessary if more exoctic relaxation options
         !  are considered.
         !

         !
         !  Call BoxMG to perform setup for planes on level kg
         !
         BMG_iPARMS_PL_xy(id_BMG2_SETUP,kg) = BMG_SETUP_only

         !
         !  Set stencil size
         !
         IF (kg.NE.NOG .OR. IFD.NE.1 ) THEN
            BMG_iPARMS_PL_xy(id_BMG2_STENCIL,kg) 
     &           = BMG_STENCIL_9pt
            NStncl_2D = 5
            NStncl_3D = 14
         ELSE
            BMG_iPARMS_PL_xy(id_BMG2_STENCIL,kg) 
     &           = BMG_STENCIL_5pt
            NStncl_2D = 3
            NStncl_3D = 4
         ENDIF


         !
         !  Get pointer data for 3D 
         !
         CALL BMG3_SymStd_GET_pointers( 
     &             kg, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &             )

         !
         !  Loop over interior planes (ignore ghost planes)
         !
         DO iPL = 2, Nz-1

C               WRITE(*,*) 'iPL, kg = ', iPL, kg
            !
            ! Set error codes to zero
            !
            BMG_iPARMS_PL_xy(id_BMG2_Err_Code, KG) = 0
            BMG_iPARMS_PL_xy(id_BMG2_Ext_Err_Code, KG) = 0
            !
            !  Create BMG_pWORK_PL_xy
            !
            CALL BMG2_SymStd_SETUP_PtrWork( 
     &                Nx-2, Ny-2, BMG_iPARMS_PL_xy(1,kg),
     &                NOGm_xy, NFm_xy, NSOm_xy,
     &                NBMG_iWORK_PL, NBMG_rWORK_PL,
     &                BMG_pWORK_PL_xy(1,iPL,kg),
     &                BMG_InWORK_PL, pSR, pSI,
     &                BMG_IOFLAG 
     &                ) 

C               WRITE(*,*) 'pSR = ', pSR
C               WRITE(*,*) 'pSI = ', pSI

            !
            !  Copy Stencil in place
            !
            CALL BMG3_SymStd_COPY_SO_xy(
     &                SO(p_SO),
     &                BMG_rWORK_PL(BMG_pWORK_PL_xy(ip_SO,iPL,kg)),
     &                iPL, Nx, Ny, Nz, NStncl_3D, NStncl_2D
     &                )

            
            !
            NOG_xy  = BMG_iPARMS_PL_xy(id_BMG2_DIM_NOG,kg)
            NF_xy   = BMG_iPARMS_PL_xy(id_BMG2_DIM_NF,kg)
            NC_xy   = BMG_iPARMS_PL_xy(id_BMG2_DIM_NC,kg)
            NSO_xy  = BMG_iPARMS_PL_xy(id_BMG2_DIM_NSO,kg)
            NSOR_xy = BMG_iPARMS_PL_xy(id_BMG2_DIM_NSOR,kg)
            NCI_xy  = BMG_iPARMS_PL_xy(id_BMG2_DIM_NCI,kg)
            NCBW_xy = BMG_iPARMS_PL_xy(id_BMG2_DIM_NCBW,kg)
            NCU_xy  = BMG_iPARMS_PL_xy(id_BMG2_DIM_NCU,kg)
            !
            p_SO_xy  = BMG_pWORK_PL_xy(ip_SO,iPL,kg)
            p_U_xy   = BMG_pWORK_PL_xy(ip_U,iPL,kg)
            p_Q_xy   = BMG_pWORK_PL_xy(ip_Q,iPL,kg)
            p_RES_xy = BMG_pWORK_PL_xy(ip_RES,iPL,kg)
            p_SOR_xy = BMG_pWORK_PL_xy(ip_SOR,iPL,kg)
            p_CI_xy  = BMG_pWORK_PL_xy(ip_CI,iPL,kg)
            p_CSO_xy = BMG_pWORK_PL_xy(ip_CSO,iPL,kg)
            p_CU_xy  = BMG_pWORK_PL_xy(ip_CU,iPL,kg)
            p_IG_xy  = BMG_pWORK_PL_xy(ip_IG,iPL,kg)

C   
C               WRITE(*,*) 
C               WRITE(*,*) '**** Calling Boxmg with ... '
C               WRITE(*,*) 
C               WRITE(*,*) ' NOG_xy  = ',
C        &                 BMG_iPARMS_PL_xy(id_BMG2_DIM_NOG,kg)
C               WRITE(*,*) ' NF_xy   = ',
C        &                 BMG_iPARMS_PL_xy(id_BMG2_DIM_NF,kg)
C               WRITE(*,*) ' NC_xy   = ', 
C        &                 BMG_iPARMS_PL_xy(id_BMG2_DIM_NC,kg)
C               WRITE(*,*) ' NSO_xy  = ', 
C        &                 BMG_iPARMS_PL_xy(id_BMG2_DIM_NSO,kg)
C               WRITE(*,*) ' NSOR_xy = ', 
C        &                 BMG_iPARMS_PL_xy(id_BMG2_DIM_NSOR,kg)
C               WRITE(*,*) ' NCI_xy  = ',
C        &                 BMG_iPARMS_PL_xy(id_BMG2_DIM_NCI,kg)
C               WRITE(*,*) ' NCBW_xy = ', 
C        &                 BMG_iPARMS_PL_xy(id_BMG2_DIM_NCBW,kg)
C               WRITE(*,*) ' NCU_xy  = ', 
C        &                 BMG_iPARMS_PL_xy(id_BMG2_DIM_NCU,kg)
C               WRITE(*,*) 
C               WRITE(*,*) ' p_SO_xy  = ', BMG_pWORK_PL_xy(ip_SO,iPL,kg)
C               WRITE(*,*) ' p_U_xy   = ', BMG_pWORK_PL_xy(ip_U,iPL,kg)
C               WRITE(*,*) ' p_Q_xy   = ', BMG_pWORK_PL_xy(ip_Q,iPL,kg)
C               WRITE(*,*) ' p_RES_xy = ', BMG_pWORK_PL_xy(ip_RES,iPL,kg)
C               WRITE(*,*) ' p_SOR_xy = ', BMG_pWORK_PL_xy(ip_SOR,iPL,kg)
C               WRITE(*,*) ' p_CI_xy  = ', BMG_pWORK_PL_xy(ip_CI,iPL,kg)
C               WRITE(*,*) ' p_CSO_xy = ', BMG_pWORK_PL_xy(ip_CSO,iPL,kg)
C               WRITE(*,*) ' p_CU_xy  = ', BMG_pWORK_PL_xy(ip_CU,iPL,kg)
C               WRITE(*,*) ' p_IG_xy  = ', BMG_pWORK_PL_xy(ip_IG,iPL,kg)

            CALL BMG2_SymStd_SOLVE_boxmg( 
     &                Nx-2, Ny-2, 
     &                BMG_iPARMS_PL_xy(1,kg), BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK_PL(p_U_xy), BMG_rWORK_PL(p_Q_xy),
     &                BMG_rWORK_PL(p_RES_xy), NF_xy, NC_xy,
     &                BMG_rWORK_PL(p_SO_xy), NSO_xy,
     &                BMG_rWORK_PL(p_SOR_xy), NSOR_xy,
     &                BMG_rWORK_PL(p_CI_xy), NCI_xy, 
     &                BMG_rWORK_PL(p_CSO_xy), BMG_rWORK_PL(p_CU_xy),
     &                NCBW_xy, NCU_xy, BMG_iWORK_PL(p_iG_xy),
     &                NOG_xy, KNOG
     &                )

C               WRITE(*,*) ' Effective maximum NOG = ', NOG_xy
C               WRITE(*,*) ' Computed NOG          = ', KNOG
C               WRITE(*,*) 

         ENDDO
         !
         ! We are finished with the setup for planes on level kg
         !
         BMG_iPARMS_PL_xy(id_BMG2_SETUP,kg) = BMG_SETUP_none
         !
      ENDDO
         
      !
      !   Store the starting point of free workspace
      !
      BMG_iWORK_PL(ip_BMG_iWORK_PL_SF) = pSI
      BMG_iWORK_PL(ip_BMG_rWORK_PL_SF) = pSR

      RETURN
      END
