      LOGICAL FUNCTION SPACE(ELEM,CIRC,FSPACE,GAPSQ,NODE,RSQ)

C#### Function: SPACE
C###  Description:
C###    Determines whether there is space available to insert
C###    an internal point. Genmesh function.
CC JMB 20-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4)
      REAL*8 CIRC(LDCIRC,3),FSPACE(0:N_GM),GAPSQ,NODE(LDNODE,3),
     '  RSQ

      CALL GET_SPACE(ELEM,CIRC,FSPACE,NODE)
      IF(RSQ.GE.GAPSQ)THEN
        SPACE=.TRUE.
      ELSE
        SPACE=.FALSE.
      ENDIF !rsq

      RETURN
      END


