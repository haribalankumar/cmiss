      SUBROUTINE BMG2_SymStd_ErrTrap( 
     &                BMG_iPARMS, IERROR
     &                )

C ==========================================================================
C
C   BMG2_SymStd_ErrTrap.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_ErrTrap.f 
C   implements the a user defined code behavior in case of an error.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2003/09/10 (MB)
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
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_constants.h'

C ----------------------------
C     Argument Declarations
C 

      INTEGER  IERROR
      INTEGER  BMG_iPARMS(NBMG_iPARMS) 

C ----------------------------
C     Local Declarations
C

C ==========================================================================
      
      !
      ! first test if an error has occured
      !
      IF (IERROR .ne. iZERO) THEN
         !
         ! store the error code in BMG_iPARMS
         !
         BMG_iPARMS(id_BMG2_Err_Code) = IERROR
      ENDIF

C ===========================================

      RETURN
      END
