      LOGICAL FUNCTION CELEM(CHAR,SET)

C#### Function: CELEM
C###  Type: LOGICAL
C###  Description:
C###    CELEM returns .TRUE. if CHAR is an element of the set
C###    of characters SET.

      IMPLICIT NONE
!     Parameter List
      CHARACTER CHAR*1,SET*(*)
!     Local Variables
      INTEGER n1
      CHARACTER C1*1,C2*1

      CELEM=.FALSE.
      DO n1=1,LEN(SET)
        CALL CUPPER(CHAR,C1)
        CALL CUPPER(SET(n1:n1),C2)
C        IF(CUPPER(CHAR).EQ.CUPPER(SET(n1:n1))) THEN
        IF(C1.EQ.C2) THEN
          CELEM=.TRUE.
          RETURN
        ENDIF
      ENDDO
      RETURN
C GMH 12/12/96 This does not make sense
C 9999 WRITE(OP_STRING,'('' >>ERROR: RETURN FROM CUPPER'',A)') ERROR
C      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
C 9998 RETURN
      END


