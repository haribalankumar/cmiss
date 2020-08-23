
      SUBROUTINE BMG_ALLOC(
     &               IWORK_PTR, NIWORK, IWORK_PL_PTR, NIWORK_PL,
     &               RWORK_PTR, NRWORK, RWORK_PL_PTR, NRWORK_PL,
     &               ERROR, *)

C#### Subroutine: BMG_ALLOC
C###  Description:
C###    BMG_ALLOC allocates workspace memory needed to call Boxmg.
C###  Written by Travis Austin 18/01/05

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
!     Parameter List
      CHARACTER ERROR*(*)
      INTEGER NIWORK,NRWORK,NIWORK_PL,NRWORK_PL
      POINTER IWORK_PTR,RWORK_PTR,IWORK_PL_PTR,RWORK_PL_PTR
!     Local Variables
      LOGICAL INIT

      INIT = .TRUE.

!     Allocate INTEGER arrays
      CALL ALLOCATE_MEMORY(
     &              NIWORK,1,INTTYPE,IWORK_PTR,INIT,ERROR,*9999) 
      CALL ALLOCATE_MEMORY(
     &              NIWORK_PL,1,INTTYPE,IWORK_PL_PTR,INIT,ERROR,*9999) 

!     Allocate REAL arrays
      CALL ALLOCATE_MEMORY(
     &              NRWORK,1,DPTYPE,RWORK_PTR,INIT,ERROR,*9999) 
      CALL ALLOCATE_MEMORY(
     &              NRWORK_PL,1,DPTYPE,RWORK_PL_PTR,INIT,ERROR,*9999) 
      RETURN

 9999 CALL ERRORS('BMG_ALLOC',ERROR)
      RETURN 1
      END

      SUBROUTINE BMG_FREE( 
     &              IWORK_PTR,IWORK_PL_PTR,
     &              RWORK_PTR,RWORK_PL_PTR,
     &              ERROR,*)

C#### Subroutine: BMG_FREE
C###  Description:
C###    BMG_FREE frees memory allocated with BMG_ALLOC.
C###  Written by Travis Austin 18/01/05
 
      IMPLICIT NONE 
!     Parameter List
      CHARACTER ERROR*(*)
      POINTER IWORK_PTR,RWORK_PTR,
     &        IWORK_PL_PTR,RWORK_PL_PTR

      IF(IWORK_PTR.NE.0) CALL FREE_MEMORY(IWORK_PTR,ERROR,*9999)
      IF(RWORK_PTR.NE.0) CALL FREE_MEMORY(RWORK_PTR,ERROR,*9999)

      IF(IWORK_PL_PTR.NE.0) CALL FREE_MEMORY(IWORK_PL_PTR,ERROR,*9999)
      IF(RWORK_PL_PTR.NE.0) CALL FREE_MEMORY(RWORK_PL_PTR,ERROR,*9999)
      
      RETURN

 9999 CALL ERRORS('BMG_FREE',ERROR)
      RETURN 1
      END

