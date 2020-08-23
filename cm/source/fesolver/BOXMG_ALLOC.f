      SUBROUTINE BOXMG_ALLOC(NX,ERROR,*)

C#### Subroutine: BOXMG_ALLOC
C###  Description:
C###    BOXMG_ALLOC allocates memory for the BOXMG arrays needed to 
C###    construct coarse grid operators.  On completion the SOLV_ALLOC 
C###    variable is set to true to flag that memory is allocated.

      IMPLICIT NONE 

      INCLUDE 'bmg00.cmn' 
      INCLUDE 'ptr00.cmn'

!     Parameter List
      INTEGER NX
      CHARACTER ERROR*(*)

      CALL ENTERS('BOXMG_ALLOC',*9999)

      CALL ASSERT(nx.LE.9,'>>Can only have nx.le.9  at present',
     &  ERROR,*9999)

      CALL BMG_ALLOC(
     &        BMG_iWORK_PTR(NX),NBMG_iWORK(NX),
     &        BMG_iWORK_PL_PTR(NX),NBMG_IWORK_PL(NX),
     &        BMG_rWORK_PTR(NX),NBMG_rWORK(NX),
     &        BMG_rWORK_PL_PTR(NX),NBMG_RWORK_PL(NX),
     &        ERROR,*9999) 

C     Flag memory as allocated
      SOLV_ALLOC(NX)=.TRUE. 

      CALL EXITS('BOXMG_ALLOC') 
      RETURN 

 9999 CALL ERRORS('BOXMG_ALLOC',ERROR) 
      CALL EXITS('BOXMG_ALLOC') 
      RETURN 1 
      END 


