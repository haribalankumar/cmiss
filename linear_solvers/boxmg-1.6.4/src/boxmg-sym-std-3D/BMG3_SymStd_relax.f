      SUBROUTINE BMG3_SymStd_relax( 
     &                kg, kf, IGRD,
     &                Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IFD, NStncl, IRELAX, IRELAX_SYM,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_IJK_MAP
     &                )


C ==========================================================================
C
C   BMG3_SymStd_relax
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_relax performs one relaxation step (either colored
C   point Gauss Seidel, or alternating plane red-black Gauss-Seidel)
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   2000/02/23  - written (M.Berndt)
C
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
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ------------------------------------------------
C     Argument Declarations
C 
      INTEGER  NCI, NSO, NSOR, kg, kf, NBMG_iWORK_PL,
     &         NBMG_rWORK_PL, NOG, NF, NOGm, NStncl

      INTEGER  BMG_iPARMS(NBMG_iPARMS),
     &         BMG_iWORK_PL(NBMG_iWORK_PL),
     &         IFD, IGRD(NOGm,NBMG_pIGRD),
     &         IRELAX_SYM, UPDOWN, IRELAX,
     &         BMG_IJK_MAP(*)
      REAL*8   BMG_rPARMS(NBMG_rPARMS),
     &         BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), q(NF), qf(NF), RES(NF),
     &         RES_L2, SO(NSO), SOR(NSOR)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  NSORv, Nx, Nxf, Ny, Nyf, Nz, Nzf, p_CI, p_SO, p_SOR, p_U,
     &         p_IJK_MAP
      REAL*8  starttime, endtime, OMP_GET_WTIME
      EXTERNAL OMP_GET_WTIME

C ========================================================================

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      !
      ! Collect fine-grid dimensions
      !
      CALL BMG3_SymStd_GET_pointers( 
     &          NOG, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, Nxf, Nyf, Nzf 
     &          )

      CALL BMG3_SymStd_GET_pointers( 
     &          kg, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &          )
      

C      starttime = OMP_GET_WTIME()

      IF ( IRELAX.EQ.BMG_GS_RB_point .OR.
     &     IRELAX.EQ.BMG_GS_adaptive ) THEN
         !
         ! Gauss-Seidel relaxation
         !
         p_IJK_MAP = BMG_IJK_MAP(kg)
         CALL BMG3_SymStd_relax_ptwise_GS( 
     &             kg, SO(p_SO), qf(p_U), q(p_U), RES(p_U), SOR(p_SOR), 
     &             Nx, Ny, Nz, RES_L2, BMG_IOFLAG,
     &             KF, IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,
     &             BMG_IJK_MAP(p_IJK_MAP)
     &             )

C         endtime = OMP_GET_WTIME()
C         PRINT *,kg,NOG,' Computing Point Relax Time = ', 
C     &        endtime-starttime
         !
      ELSEIF ( IRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN

         !
         ! Alternating plane relaxation
         !
         IF ( ( IRELAX_SYM .EQ.BMG_RELAX_NONSYM )
     &      .OR. ( UPDOWN.EQ.BMG_DOWN )         )THEN
            !
            ! On the way down, or in the non-symmetric case, 
            ! relax in the order xy, yz, xz
            !
            CALL BMG3_SymStd_relax_planes_xy(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xy)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xy)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_IJK_MAP
     &                )

            CALL BMG3_SymStd_relax_planes_yz(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_yz)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_yz)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_IJK_MAP
     &                )
            !
            CALL BMG3_SymStd_relax_planes_xz(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xz)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xz)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_IJK_MAP
     &                )

C            endtime = OMP_GET_WTIME()
C            PRINT *,kg,NOG,' Computing Plane Relax Time = ', 
C     &           endtime-starttime
            !
         ELSE

            !
            ! on the way up, in the symmetric case, use
            ! the opposite order xz, yz, xy
            !
            CALL BMG3_SymStd_relax_planes_xz(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xz)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xz)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_IJK_MAP
     &                )

              !
            CALL BMG3_SymStd_relax_planes_yz(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_yz)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_yz)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_IJK_MAP
     &                )

            !
            CALL BMG3_SymStd_relax_planes_xy(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xy)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xy)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_IJK_MAP
     &                )

C            endtime = OMP_GET_WTIME()
C            PRINT *,kg,NOG,' Computing Plane Relax Time = ', 
C     &           endtime-starttime

            !
         ENDIF

      ELSEIF ( IRELAX.EQ.BMG_CG_relaxation ) THEN

         !
         !  Use a few iterations of Conjugate Gradients as smoother
         !
         CALL BMG3_SymStd_relax_CG( 
     &             kg, SO(p_SO), qf(p_U), q(p_U), RES(p_U), SOR(p_SOR), 
     &             Nx, Ny, Nz, RES_L2, BMG_IOFLAG, KF, IFD, NStncl, 
     &             BMG_IJK_MAP
     &             )         

      ENDIF
    

C ==========================================================================

      RETURN
      END









