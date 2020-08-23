      SUBROUTINE BMG2_SymStd_updown(
     &                K, KF, UPDOWN, IFD, IBC, IU, ID,
     &                IRELAX, IRELAX_SYM, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm, RES_L2, BMG_iPARMS
     &                )

C ==========================================================================
C
C   BMG2_SymStd_updown.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_updown performs all the necessary tasks associated
C   with moving "down" to a coarser grid or moving "up" to 
C   a finer grid.  It is necessary because f77 does not support
C   recursion.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/02/17 (JDM)
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
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  NC, NCI, NF, NOGm, NSO, NSOR

      INTEGER  IBC, ID, IFD, iGRD(NOGm,9), IRELAX, IRELAX_SYM, IU, 
     &         K, KF, UPDOWN
      REAL*8   CI(NCI), Q(NF), QF(NF), RES(NF), RES_L2, 
     &         SO(NSO), SOR(NSOR)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_iPARMS)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, NSORv, NStncl, Nx, Nxc, Ny, Nyc, p_CI, p_CIC,
     &         p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC, myid
      REAL*8   TAUERR     !! Temporary 

C =========================================================================


      !
      ! How many points in the stencil?
      !
      IF ( K.EQ.KF .AND. IFD.EQ.1 ) THEN
         NStncl=3
      ELSE
         NStncl=5
      ENDIF

      !
      !  Number of temporary vectors in SOR
      !
      IF ( iRELAX.EQ.BMG_GS_RB_point 
     &    .OR. iRELAX.EQ.BMG_GS_RB_x_lines ) THEN
         NSORv = 2
      ELSE 
         NSORv = 4
      ENDIF

      !  
      !  Get IGRD pointers
      !
      CALL BMG2_SymStd_GET_pointers(
     &          K, iGRD, NOGm, Nx, Ny, 
     &          p_U, p_SO, p_SOR, p_CI
     &          )
      CALL BMG2_SymStd_GET_pointers(
     &          K-1, iGRD, NOGm, Nxc, Nyc,
     &          p_UC, p_SOC, p_SORC, p_CIC
     &          )

      IF (UPDOWN.EQ.BMG_DOWN) THEN
         !
         ! Relaxation
         !
         IF ( IBC.EQ.BMG_BCs_definite .OR.
     &        IBC.EQ.BMG_BCs_indef_nonper ) THEN
            !
            DO i=1, ID
               CALL BMG2_SymStd_relax ( 
     &              K, SO(p_SO), QF(p_U), Q(p_U), 
     &              RES(p_U), SOR(p_SOR), NSOR, 
     &              Nx, Ny, KF, IFD, NStncl, IBC, 
     &              IRELAX, IRELAX_SYM, BMG_DOWN,
     &              BMG_IOFLAG, BMG_iPARMS 
     &              )
               
               IF( BMG_IOFLAG(iBMG2_BUG_RES_RELAX) ) THEN
                  CALL BMG2_SymStd_residual( 
     &                 K, SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &                 Nx, Ny, KF, IFD, NStncl
     &                 )
                  CALL BMG2_SymStd_UTILS_norm_l2(
     &                 RES(p_U), Nx, Ny, RES_L2
     &                 )
                  WRITE (*,100) K, RES_L2
               ENDIF
            END DO
            !
         ELSE
            ! Relaxation + residual
            DO i=1, ID
               !
               CALL BMG2_PerSymStd_relax(
     &                   K, SO(p_SO), QF(p_U), Q(p_U), SOR(p_SOR),
     &                   Nx, Ny, RES_L2, KF, IFD, NStncl, NSORv,
     &                   IRELAX, IBC
     &                   )
               !
               IF( BMG_IOFLAG(iBMG2_BUG_RES_RELAX) ) THEN
                  WRITE (*,100) K, RES_L2
               ENDIF
            ENDDO
         ENDIF
         !
         ! Restriction
         !
         IF ( IBC.EQ.BMG_BCs_definite .OR. 
     &        IBC.EQ.BMG_BCs_indef_nonper ) THEN
            ! Residual calculation
            CALL BMG2_SymStd_residual(
     &                K, SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &                Nx, Ny, KF, IFD, NStncl
     &                )
            ! Restriction
            CALL BMG2_SymStd_restrict( 
     &                K, K-1, RES(p_U), QF(p_UC), CI(p_CIC),
     &                Nx, Ny, Nxc, Nyc, IBC
     &                )
            
         ELSE
            ! Restriction
            CALL BMG2_PerSymStd_residual(
     &                K, K-1, iONE, 
     &                SO(p_SO), SO(p_SOC), QF(p_U), QF(p_UC), 
     &                Q(p_U), Q(p_UC), CI(p_CIC), SOR(p_SOR), 
     &                Nxc, Nyc, Nx, Ny, TAUERR, NStncl, NSORv, IBC
     &                )
         ENDIF
         ! Zero the initial guess
         CALL BMG2_SymStd_UTILS_rV_zero(
     &             Q(p_UC), Nxc, Nyc 
     &             )

      ELSE IF (UPDOWN.EQ.BMG_UP) THEN
         !
         !  Interpolation
         ! 
         IF ( IBC.EQ.BMG_BCs_definite .OR. 
     &        IBC.EQ.BMG_BCs_indef_nonper ) THEN
            !
            CALL BMG2_SymStd_interp_add( 
     &                K-1, K, Q(p_U), Q(p_UC), RES(p_U), SO(p_SO),
     &                CI(p_CIC), Nxc, Nyc, Nx, Ny, NStncl, IBC                
     &                )
            !
         ELSE
            !
            CALL BMG2_PerSymStd_interp_add(
     &                K-1, K, Q(p_U), Q(p_UC), SOR(p_SOR), 
     &                CI(p_CIC), Nxc, Nyc, Nx, Ny, NSORv, IBC 
     &                )
            !
         ENDIF

         !
         ! Relaxation
         !
         IF ( IBC.EQ.BMG_BCs_definite .OR.
     &        IBC.EQ.BMG_BCs_indef_nonper ) THEN
            !
            !
            DO i=1, IU
               CALL BMG2_SymStd_relax( 
     &              K, SO(p_SO), QF(p_U), Q(p_U),
     &              RES(p_U), SOR(p_SOR), NSOR,
     &              Nx, Ny, KF, IFD, NStncl, IBC, 
     &              IRELAX, IRELAX_SYM, BMG_UP,
     &              BMG_IOFLAG, BMG_iPARMS
     &              )
                  
               IF( BMG_IOFLAG(iBMG2_BUG_RES_RELAX) ) THEN
                  CALL BMG2_SymStd_residual( 
     &                 K, SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &                 Nx, Ny, KF, IFD, NStncl
     &                 )
                  CALL BMG2_SymStd_UTILS_norm_l2( 
     &                 RES(p_U), Nx, Ny, RES_L2
     &                 )
                  WRITE (*,100) K, RES_L2
               ENDIF
            ENDDO
            !
         ELSE
            ! Relaxation + residual
            DO i=1, IU
               CALL BMG2_PerSymStd_relax(
     &                   K, SO(p_SO), QF(p_U), Q(p_U), SOR(p_SOR),
     &                   Nx, Ny, RES_L2, KF, IFD, NStncl, NSORv, 
     &                   IRELAX, IBC
     &                   )
               IF( BMG_IOFLAG(iBMG2_BUG_RES_RELAX) ) THEN
                  WRITE (*,100) K, RES_L2
               ENDIF
            ENDDO
         ENDIF
         
      ELSE
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) '*** UPDOWN out of range: '
            WRITE(*,510) 'HAVE: UPDOWN = ', UPDOWN
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,10)
         RETURN
         
      ENDIF

C ==========================================================================

 100  FORMAT (' LEVEL',I2,' RESIDUAL NORM= ',1P,E10.3)

C -----------------------------                                                 
 500    FORMAT (/,'FATAL ERROR: BMG2_SymStd_updown.f',/,5X,A)
 510    FORMAT (5X,A,1X,I2)

C =============================

      RETURN
      END

