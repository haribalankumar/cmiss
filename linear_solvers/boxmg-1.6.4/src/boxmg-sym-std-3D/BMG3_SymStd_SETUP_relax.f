      SUBROUTINE BMG3_SymStd_SETUP_relax( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                IFD, IBC, IRELAX, NOGm, NOG, IGRD,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

C ==========================================================================
C
C   BMG3_SymStd_SETUP_relax.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_SETUP_relax.f performs any necessary setup for 
C   the chosen relaxation scheme.  In particular, it calls routines
C   that setup the operators for plane relaxation.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/11/15 (JDM)
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

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'


C ---------------------------
C    Argument Declarations:
C
      INTEGER  NBMG_iWORK_PL, NBMG_rWORK_PL, NCI, NOGm, NSO, NSOR

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_iWORK_PL(NBMG_iWORK_PL),
     &         IBC, IFD, IGRD(24*NOGm), IRELAX, NOG
      REAL*8   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), SO(NSO), SOR(NSOR)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ---------------------------
C    Local Declarations:
C
      INTEGER  K, NSORv, NStncl, Nx, Nxf, Ny, Nyf, Nz, Nzf, 
     &         p_CI, p_SO, p_SOR, p_U

C ==========================================================================
      
      !
      ! Sanity check
      !
      IF ( NOG.EQ.1) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN 
            WRITE(*,500) 'There is only 1 grid!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,6)
         RETURN

      ELSE IF ( NOG.EQ.0 ) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN 
            WRITE(*,500) 'There are no grids?'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,7)
         RETURN
         
      ELSE IF ( NOG.LT.0 ) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN 
            WRITE(*,500) 'There number of grids is negative!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,8)
         RETURN
         
      ENDIF

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      IF ( IRELAX.EQ.BMG_GS_RB_point .OR. 
     &     IRELAX.EQ.BMG_GS_adaptive ) THEN

         !
         ! Loop over grids
         !
         DO K = NOG, 2, -1
            
            IF ( K.NE.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
               NStncl=14
            ELSE
               NStncl=4
            ENDIF

            CALL BMG3_SymStd_GET_pointers( 
     &                K, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &                )

            CALL BMG3_SymStd_SETUP_recip(
     &                SO(p_SO), SOR(p_SOR), Nx, Ny, Nz, Nstncl, NSORv
     &                )
            
         ENDDO

      ENDIF

      IF ( IRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz .OR. 
     &     IRELAX.EQ.BMG_GS_adaptive ) THEN

         !
         ! Collect fine-grid dimensions (Nxf,Nyf,Nzf)
         !
         CALL BMG3_SymStd_GET_pointers( 
     &             NOG, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, Nxf, Nyf, Nzf 
     &             )

         CALL BMG3_SymStd_SETUP_planes_xy( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             Nxf, Nyf, Nzf, IFD, IBC, IGRD, 
     &             SO, NSO, NOGm, NOG,
     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xy)),
     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xy)),
     &             BMG_iWORK_PL, NBMG_iWORK_PL,
     &             BMG_rWORK_PL, NBMG_rWORK_PL
     &             )

         CALL BMG3_SymStd_SETUP_planes_yz( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             Nxf, Nyf, Nzf, IFD, IBC, IGRD, 
     &             SO, NSO, NOGm, NOG,
     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_yz)),
     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_yz)),
     &             BMG_iWORK_PL, NBMG_iWORK_PL,
     &             BMG_rWORK_PL, NBMG_rWORK_PL
     &             )

        CALL BMG3_SymStd_SETUP_planes_xz( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             Nxf, Nyf, Nzf, IFD, IBC, IGRD, 
     &             SO, NSO, NOGm, NOG,
     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xz)),
     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xz)),
     &             BMG_iWORK_PL, NBMG_iWORK_PL,
     &             BMG_rWORK_PL, NBMG_rWORK_PL
     &             )

      ENDIF

      IF ( iRELAX.EQ.BMG_CG_relaxation ) THEN

         !
         ! Loop over grids
         !
         DO K = NOG, 2, -1
            
            IF ( K.NE.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
               NStncl=14
            ELSE
               NStncl=4
            ENDIF

            CALL BMG3_SymStd_GET_pointers( 
     &                K, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &                )

            CALL BMG3_SymStd_SETUP_cg_recip(
     &                SO(p_SO), SOR(p_SOR), Nx, Ny, Nz, Nstncl
     &                )
            
         ENDDO

      ENDIF


      
C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_relax.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C =====================

      RETURN
      END


