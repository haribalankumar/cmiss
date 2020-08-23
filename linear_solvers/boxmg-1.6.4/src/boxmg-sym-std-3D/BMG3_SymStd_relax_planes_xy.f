      SUBROUTINE BMG3_SymStd_relax_planes_xy( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iPARMS_PL_xy, BMG_pWORK_PL_xy,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_IJK_MAP
     &                )

C ==========================================================================
C
C   BMG3_SymStd_relax_planes_xy
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_relax_planes_xy performs one relaxation step of 
C   red-black xy-plane relaxation (in the correct order)
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   2000/02/23  - written (M.Berndt)
C   2003/06/27  - major rewrite of memory management (JDM)
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
C ========================================================================

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
      INTEGER  kg, NCI, NBMG_iWORK_PL, NBMG_rWORK_PL,
     &         NF, NOG, NOGm, NSO, NSOR, Nxf, Nyf, Nzf

      INTEGER  BMG_iPARMS(NBMG_iPARMS), 
     &         BMG_iPARMS_PL_xy(NBMG_iPARMS,NOG),
     &         BMG_iWORK_PL(NBMG_iWORK_PL),
     &         BMG_pWORK_PL_xy(NBMG_pWORK,Nzf,NOG),
     &         IFD, IGRD(NOGm,NBMG_pIGRD), IRELAX_SYM, UPDOWN,
     &         BMG_IJK_MAP(*)

      REAL*8   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK_PL(NBMG_rWORK_PL), 
     &         CI(NCI), Q(NF), QF(NF), RES(NF), RES_L2,    
     &         SO(NSO), SOR(NSOR)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      
C -------------------------------------------------
C    Local Declarations
C
      INTEGER  iPL, iPL_beg, iPL_end, Nx, Ny, Nz,
     &         LSTART, LEND, LSTRIDE,
     &         p_CIC, p_SO, p_SOR, p_U

      INTEGER  NFm_xy, NOGm_xy, NSOm_xy

      INTEGER  NC_xy, NCBW_xy, NCI_xy, NCU_xy, NF_xy, NOG_xy,
     &         NSO_xy, NSOR_xy, NStncl_2D, NStncl_3D 

      INTEGER  p_CI_xy, p_CSO_xy, p_CU_xy, p_IG_xy, p_SO_xy, p_SOR_xy,
     &         p_Q_xy, p_RES_xy, p_U_xy

C ===========================================================================
      
      IRELAX_SYM = BMG_iPARMS_PL_xy(id_BMG2_RELAX_SYM,kg)

      IF ( ( IRELAX_SYM.EQ.BMG_RELAX_NONSYM ) 
     &   .OR. ( UPDOWN.EQ.BMG_DOWN )       ) THEN
         LSTART  =  2
         LEND    =  3
         LSTRIDE =  1 
      ELSE
         LSTART  =  3
         LEND    =  2
         LSTRIDE = -1
      ENDIF

      !
      !  Set stencil size
      !
      IF (kg.NE.NOG .OR. IFD.NE.1 ) THEN
         NStncl_2D = 5
         NStncl_3D = 14
      ELSE
         NStncl_2D = 3
         NStncl_3D = 4
      ENDIF

      !
      ! Store maximum dimensions for (x,y) - planes
      ! 
      NOGm_xy = BMG_iPARMS_PL_xy(id_BMG2_DIM_NOG,NOG)
      NFm_xy  = BMG_iPARMS_PL_xy(id_BMG2_DIM_NF,NOG)
      NSOm_xy = BMG_iPARMS_PL_xy(id_BMG2_DIM_NSO,NOG)

C -------------------------------------------------------------------------
C    >>>>>>>>>>>>>>>>>>>>>> BEGIN: xy sweep <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C -------------------------------------------------------------------------

      CALL BMG3_SymStd_GET_pointers( 
     &          kg, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CIC, Nx, Ny, Nz 
     &          )

      NOG_xy  = BMG_iPARMS_PL_xy(id_BMG2_DIM_NOG,kg)
      NF_xy   = BMG_iPARMS_PL_xy(id_BMG2_DIM_NF,kg)
      NC_xy   = BMG_iPARMS_PL_xy(id_BMG2_DIM_NC,kg)
      NSO_xy  = BMG_iPARMS_PL_xy(id_BMG2_DIM_NSO,kg)
      NSOR_xy = BMG_iPARMS_PL_xy(id_BMG2_DIM_NSOR,kg)
      NCI_xy  = BMG_iPARMS_PL_xy(id_BMG2_DIM_NCI,kg)
      NCBW_xy = BMG_iPARMS_PL_xy(id_BMG2_DIM_NCBW,kg)
      NCU_xy  = BMG_iPARMS_PL_xy(id_BMG2_DIM_NCU,kg)
         !
      DO iPL_beg=LSTART, LEND, LSTRIDE
         !
         iPL_end=2*((Nz-1-iPL_beg)/2)+iPL_beg
         !
         BMG_rPARMS(id_BMG2_STOP_TOL) = rZERO
         !
CC$OMP    PARALLEL DO 
CC$OMP&    DEFAULT(NONE)
CC$OMP&    PRIVATE(iPL,p_SO_xy,p_U_xy,p_Q_xy,p_RES_xy,p_SOR_xy,
CC$OMP&            p_CI_xy,p_CSO_xy,p_CU_xy,p_IG_xy)
CC$OMP&    SHARED(BMG_IOFLAG,BMG_rPARMS,BMG_iWORK_PL,BMG_rWORK_PL,
CC$OMP&           BMG_pWORK_PL_xy,BMG_iPARMS_PL_xy,kg,NStncl_3D,
CC$OMP&           NF_xy,NC_xy,NSO_xy,NSOR_xy, NCI_xy,NCBW_xy,
CC$OMP&           NCU_xy,NOG_xy,Nx,Ny,Nz,p_U,p_SO,Q,QF,SO,
CC$OMP&           iPL_beg,iPL_end)
         DO iPL=iPL_beg, iPL_end, 2

            p_SO_xy  = BMG_pWORK_PL_xy(ip_SO,iPL,kg)
            p_U_xy   = BMG_pWORK_PL_xy(ip_U,iPL,kg)
            p_Q_xy   = BMG_pWORK_PL_xy(ip_Q,iPL,kg)
            p_RES_xy = BMG_pWORK_PL_xy(ip_RES,iPL,kg)
            p_SOR_xy = BMG_pWORK_PL_xy(ip_SOR,iPL,kg)
            p_CI_xy  = BMG_pWORK_PL_xy(ip_CI,iPL,kg)
            p_CSO_xy = BMG_pWORK_PL_xy(ip_CSO,iPL,kg)
            p_CU_xy  = BMG_pWORK_PL_xy(ip_CU,iPL,kg)
            p_IG_xy  = BMG_pWORK_PL_xy(ip_IG,iPL,kg)

            ! solve the plane

            !
            ! Copy the iPL{th} plane of the 3D solution into 2D
            !
            CALL BMG3_SymStd_COPY_rV_32_xy( 
     &                Q(p_U), BMG_rWORK_PL(p_U_xy), iPL, Nx, Ny, Nz 
     &                )
            
            CALL BMG3_SymStd_COPY_RHS_xy( 
     &                SO(p_SO), Q(p_U), QF(p_U), BMG_rWORK_PL(p_Q_xy), 
     &                iPL, Nx, Ny, Nz, NStncl_3D
     &                )

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
     &                NOG_xy, NOG_xy
     &                )

            !
            ! Copy the 2D solution back in the iPL{th} plane of 3D
            !
            CALL BMG3_SymStd_COPY_rV_23_xy( 
     &                Q(p_U), BMG_rWORK_PL(p_U_xy), iPL, Nx, Ny, Nz 
     &                )
            !
         ENDDO
CC$OMP    END PARALLEL DO
         !
      ENDDO

      IF ( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) ) THEN
         CALL BMG3_SymStd_residual( 
     &             kg, NOG, IFD, Q(p_U), QF(p_U), SO(p_SO),
     &             RES(p_U), Nx, Ny, Nz, NStncl_3D,
     &             BMG_IJK_MAP
     &             )
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES(p_U), Nx, Ny, Nz, RES_L2
     &             )
         WRITE(*,100) 'Residual (l2-norm) after xy sweep = ', RES_L2
      ENDIF

C -------------------------------------------------------------------------
C    >>>>>>>>>>>>>>>>>>>>>> END: xy sweep <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C -------------------------------------------------------------------------

C =========================================================================

 100  FORMAT(1X,A,1X,1P,E15.7)

C ==============================

      RETURN
      END









