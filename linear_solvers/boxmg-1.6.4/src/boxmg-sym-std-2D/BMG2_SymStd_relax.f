      SUBROUTINE BMG2_SymStd_relax(
     &                K, SO, QF, Q, RES, SOR, NSOR, Nx, Ny, 
     &                KF, IFD, NStncl, IBC, IRELAX, IRELAX_SYM, UPDOWN,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )

C ==========================================================================
C
C
C***BEGIN PROLOGUE  BMG2_SymStd_relax
C***SUBSIDIARY
C***PURPOSE  MGRLX performs point or line relaxation
C            on the fine grid depending on the value of irelax.
C
C***LIBRARY   SLATEC
C***AUTHOR  DENDY, J. E. JR.
C             LOS ALAMOS NATIONAL LABORATORY
C           VOYTKO, M. H.
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C
C          Copyright, 1988. The Regents of the University of California.
C          This software was produced under a U.S. Government contract
C          (W-7405-ENG-36) by Los Alamso National Laboratory, which is
C          operated by the University of California for the U.S. Depart-
C          ment of Energy. The U.S. Government is licensed to use,
C          reproduce, and distribute this software. Permission is
C          granted to the public to copy and use this software without
C          charge, provided that this Notice and any statement of
C          authorship are reproduced on all copies. Neither the
C          Government nor the University makes any warranty, expressed
C          or implied, or assumes any liability or responsibility for
C          the use of this software.
C
C***PARAMETERS
C***INPUT
C   K         K is the grid number.
C   SO        Refer to BOXMG.
C   QF        Refer to BOXMG.
C   SOR       Refer to BOXMG.
C   Nx        Number of grid points in x direction, including
C             two fictitious points.
C   Ny        Number of grid points in y direction, including
C             two fictitious points.
C   KF        Index of the finest grid
C   IFD       Stencil parameter
C   IRELAX    Refer to BOXMG.
C***INPUT/OUTPUT
C   Q         Refer to BOXMG.
C***OUTPUT
C   ERR       ERR is the l2 error of the resiuals.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830925  DATE WRITTEN
C   900627  modified to conform to the 4/10/90 "Guide to the SLATEC
C           Common Mathematical Library" by Victor A. Bandy
C
C   991206  Changed all gotos to if statements, added indentation, and 
C           changed do-loops to end on enddo (M. Berndt)
C
C   000218  Removed the residual calculation form this subroutine. 
C           The residual caclulation is now in BMG2_SymStd_resl2 (M.Berndt)
C
C***END PROLOGUE  MGRLX

C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER Nx, Ny, NSOR, NStncl

      INTEGER IBC, IFD, IRELAX, IRELAX_SYM, K, KF, UPDOWN
      REAL*8  Q(Nx*Ny), QF(Nx*Ny), RES(Nx*Ny),
     &        SO(Nx*Ny*NStncl), SOR(NSOR)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER BMG_iPARMS(NBMG_iPARMS)

C ----------------------------
C     Local Declarations
C 
      INTEGER p_SOR, NSORv

C =========================================================================

C     NB: Workspace double time
C
C     RES is used for temporary workspace in the y-lines case.
C     It is included in the call to x-lines for consistency of
C     the calling sequence only.
C

C     -------------------------------------------------------------------
      
      !
      !  Number of temporary vectors in SOR
      !
      IF ( iRELAX.EQ.BMG_GS_RB_point 
     &    .OR. iRELAX.EQ.BMG_GS_RB_x_lines ) THEN
         NSORv = 2
      ELSE 
         NSORv = 4
      ENDIF


      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         !
         CALL BMG2_SymStd_relax_GS( 
     &             K, SO, QF, Q, SOR, Nx, Ny, 
     &             KF, IFD, NStncl, NSORv, IRELAX_SYM, UPDOWN 
     &             )
         !
      ELSE IF ( IRELAX.EQ.BMG_GS_RB_x_lines ) THEN
         !
         p_SOR = (ipL_BMG_LUL1-1)*Nx*Ny + 1
         CALL BMG2_SymStd_relax_lines_x ( 
     &             K, SO, QF, Q, SOR(p_SOR), RES, 
     &             Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &             )
         !
      ELSE IF ( IRELAX.EQ.BMG_GS_RB_y_lines ) THEN
         !
         p_SOR = (ipL_BMG_LUL2-1)*Nx*Ny + 1
         CALL BMG2_SymStd_relax_lines_y ( 
     &             K, SO, QF, Q, SOR(p_SOR), RES, 
     &             Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &             )
         !
      ELSE IF ( IRELAX.EQ.BMG_GS_RB_x_y_lines ) THEN
         ! 
         IF ( IRELAX_SYM.EQ.BMG_RELAX_NONSYM
     &       .OR.( UPDOWN.EQ.BMG_DOWN
     &            .AND.IRELAX_SYM.EQ.BMG_RELAX_SYM ) ) THEN
            !
            p_SOR = (ipL_BMG_LUL1-1)*Nx*Ny + 1
            CALL BMG2_SymStd_relax_lines_x ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &                )
            p_SOR = (ipL_BMG_LUL2-1)*Nx*Ny + 1
            CALL BMG2_SymStd_relax_lines_y ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &                )
            !
         ELSE ! on the way up relax in the opposite order
            !
            p_SOR = (ipL_BMG_LUL2-1)*Nx*Ny + 1
            CALL BMG2_SymStd_relax_lines_y ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &                )
            p_SOR = (ipL_BMG_LUL1-1)*Nx*Ny + 1
            CALL BMG2_SymStd_relax_lines_x ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &                )
            !
         ENDIF
         !
      ENDIF

      RETURN
      END
