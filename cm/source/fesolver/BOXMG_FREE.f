      SUBROUTINE BOXMG_FREE(NX,FIRST_A,ERROR,*)

C#### Subroutine: BOXMG_FREE
C###  Description:
C###    BOXMG_FREE deallocates memory for the BOXMG arrays.
C###
      IMPLICIT NONE 

      INCLUDE 'bmg00.cmn' 
      INCLUDE 'ptr00.cmn'
!      INCLUDE 'solv00.cmn' 

!     Parameter List
      INTEGER NX
      LOGICAL FIRST_A
      CHARACTER ERROR*(*)

      CALL ENTERS('BOXMG_FREE',*9999)

      CALL ASSERT(nx.LE.9,'>>Can only have nx.le.9  at present',
     &  ERROR,*9999)

!     Only deallocate if currently allocated
      IF(SOLV_ALLOC(NX)) THEN

         CALL BMG_FREE(
     &           BMG_IWORK_PTR(NX),BMG_RWORK_PTR(NX),
     &           BMG_IWORK_PL_PTR(NX),BMG_RWORK_PL_PTR(NX),
     &           ERROR,*9999)

      ENDIF

C     Flag memory as free'ed
      IF(FIRST_A) SOLV_ALLOC(NX)=.FALSE.
      
      CALL EXITS('BOXMG_FREE') 
      RETURN 

 9999 CALL ERRORS('BOXMG_FREE',ERROR) 
      CALL EXITS('BOXMG_FREE') 
      RETURN 1 
      END 

