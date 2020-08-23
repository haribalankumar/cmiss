      LOGICAL FUNCTION ENTRIES_MATCH(ARRAY,ARRAY_SIZE,ENTRIES_TO_FIND,
     '  NUM_ENTRIES_TO_FIND,ENTRY_FOUND)

C#### Function: ENTRIES_MATCH
C###  Type: LOGICAL
C###  Description:
C###    ENTRIES_MATCH finds out whether all of ENTRIES_TO_FIND are
C###    present in ARRAY

      IMPLICIT NONE
!     Parameter list
      INTEGER ARRAY(*),ARRAY_SIZE,ENTRIES_TO_FIND(*),
     '  NUM_ENTRIES_TO_FIND
      LOGICAL ENTRY_FOUND(*)
!     Local variables
      INTEGER i,j

      DO i=1,NUM_ENTRIES_TO_FIND
        ENTRY_FOUND(i)=.FALSE.
      ENDDO
      DO i=1,NUM_ENTRIES_TO_FIND
        DO j=1,ARRAY_SIZE
          IF(ARRAY(j).EQ.ENTRIES_TO_FIND(i)) THEN
            ENTRY_FOUND(i)=.TRUE.
            GOTO 10
          ENDIF
        ENDDO
 10     CONTINUE
      ENDDO
      ENTRIES_MATCH=.TRUE.
      DO i=1,NUM_ENTRIES_TO_FIND
        IF(.NOT.ENTRY_FOUND(i)) THEN
          ENTRIES_MATCH=.FALSE.
          GOTO 20
        ENDIF
      ENDDO
 20   CONTINUE
      RETURN
      END



