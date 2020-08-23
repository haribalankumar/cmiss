      SUBROUTINE BMG2_SymStd_SETUP_relax( 
     &                IFD, IBC, IRELAX, SO, NSO, SOR, NSOR, 
     &                IGRD, NOGm, NOG, BMG_IOFLAG, BMG_iPARMS
     &                )

C ========================================================================
C
C   BMG2_SymStd_SETUP_relax.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_relax.f performs any necessary setup for 
C   the chosen relaxation scheme.  In particular, it performs
C   calls routines to perform the LU factorization of the Tridiagonal 
C   systems that arise in line relaxation.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/02/28 (JDM)
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

      INTEGER  NSO, NSOR, NOGm

      INTEGER  IBC, IFD, IGRD(NOGm,9), IRELAX, NOG
      REAL*8   SO(NSO), SOR(NSOR)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_iPARMS)

C ---------------------------
C    Local Declarations:
C
      INTEGER  K, NSORv, NStncl, Nx, Ny, p_CI, p_SO, p_SOR, p_U

C ========================================================================
      
      !
      ! Sanity check
      !
      IF ( NOG.EQ.1) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'There is only 1 grid!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,11)
         RETURN

      ELSE IF ( NOG.EQ.0 ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'There are no grids?'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,12)
         RETURN
          
      ELSE IF ( NOG.LT.0 ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'There number of grids is negative!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,13)
         RETURN

      ENDIF


      !
      !  Number of temporary vectors in SOR
      !
      IF ( iRELAX.EQ.BMG_GS_RB_point 
     &   .OR. iRELAX.EQ.BMG_GS_RB_x_lines ) THEN
         NSORv = 2
      ELSE 
         NSORv = 4
      ENDIF

      !
      ! Loop over grids
      !
      DO  K = NOG, 2, -1
         !
         ! Determine the number of points in the stencil
         !
         IF (K.NE.NOG .OR. IFD.NE.1) THEN
            NStncl=5
         ELSE
            NStncl=3
         ENDIF

         !
         ! (fake) memory pointers
         !
         CALL BMG2_SymStd_GET_pointers(
     &             K, IGRD, NOGm, Nx, Ny, 
     &             p_U, p_SO, p_SOR, p_CI
     &             )

         IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
            !
            CALL BMG2_SymStd_SETUP_recip(
     &                SO(p_SO), SOR(p_SOR), Nx, Ny, NStncl, NSORv
     &                )
         ELSE IF ( IRELAX.EQ.BMG_GS_RB_x_lines ) THEN
            !
            p_SOR = p_SOR + (ipL_BMG_LUL1-1)*Nx*Ny 
            CALL BMG2_SymStd_SETUP_lines_x( 
     &                SO(p_SO), SOR(p_SOR), Nx, Ny, NStncl 
     &                )
            !
         ELSE IF ( IRELAX.EQ.BMG_GS_RB_y_lines ) THEN
            !
            p_SOR = p_SOR + (ipL_BMG_LUL2-1)*Nx*Ny
            CALL BMG2_SymStd_SETUP_lines_y(
     &                SO(p_SO), SOR(p_SOR), Nx, Ny, NStncl
     &                )
            !
         ELSE IF ( IRELAX.EQ.BMG_GS_RB_x_y_lines ) THEN
            !
            p_SOR = p_SOR + (ipL_BMG_LUL1-1)*Nx*Ny
            CALL BMG2_SymStd_SETUP_lines_x( 
     &                SO(p_SO), SOR(p_SOR), Nx, Ny, NStncl 
     &                )
            p_SOR = p_SOR + (ipL_BMG_LUL2-1)*Nx*Ny 
            CALL BMG2_SymStd_SETUP_lines_y( 
     &                SO(p_SO), SOR(p_SOR), Nx, Ny, NStncl
     &                )
            !
         ENDIF

      ENDDO

C ========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_relax.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)


C =====================

      RETURN
      END
