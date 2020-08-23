      SUBROUTINE NPLINK(NPN,NPL,ERROR,*)

C#### Subroutine: NPLINK
C###  Description:
C###    Creates an array listing the neighbours of given node
C###    in each Xi-direction.

C#### Variable: NPN(np,ni)
C###  Type: INTEGER
C###  Set_up: NPLINK
C###  Description:
C###    NPN(np,ni) lists the neigbouring node to np in the 
C###    ni xi-direction, where ni=-3,-2,-1,1,2,3

C Author: Duane Malcolm
C Created: 17 September 2003

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'mxch.inc'
      
!     Parameter List
      INTEGER NPN(NPM,-3:3),NPL(5,0:3,NLM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni,nl,np

      CALL ENTERS('NPLINK',*9999)
      
      DO np=1,NPM
        DO ni=-3,3
          NPN(np,ni)=0
        ENDDO
      ENDDO
      
      DO nl=1,NLT
        IF ((NPL(1,1,nl).EQ.1).OR.(NPL(1,1,nl).EQ.4)) THEN
        ! linear or cubic-Hermite i.e., four local nodes
          NPN(NPL(2,1,nl),NPL(1,0,nl))=NPL(3,1,nl)
          NPN(NPL(3,1,nl),-NPL(1,0,nl))=NPL(2,1,nl)
        ELSE
          CALL ASSERT(.FALSE.,
     &      '>> Zero implementation for this basis function',
     &      ERROR,*9999)
        ENDIF
      ENDDO

      CALL EXITS('NPLINK')
      RETURN
 9999 CALL ERRORS('NPLINK',ERROR)
      CALL EXITS('NPLINK')
      RETURN 1
      END


