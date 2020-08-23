      SUBROUTINE COPY_AND_INIT_INT(NEW_LENGTH,OLD_LENGTH,
     '  NEW_ARRAY,OLD_ARRAY,INIT_VALUE,ERROR,*)

C#### Subroutine: COPY_AND_INIT_INT
C###  Description:
C###    Copies data from one integer array to another.

      IMPLICIT NONE
!     Parameter List
      INTEGER NEW_LENGTH,OLD_LENGTH,
     '  NEW_ARRAY(NEW_LENGTH),OLD_ARRAY(OLD_LENGTH),INIT_VALUE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER count

      CALL ENTERS('COPY_AND_INIT_INT',*9999)
      IF(NEW_LENGTH.GT.OLD_LENGTH) THEN
        DO count=1,OLD_LENGTH
          NEW_ARRAY(count)=OLD_ARRAY(count)
        ENDDO !count
        DO count=OLD_LENGTH+1,NEW_LENGTH
          NEW_ARRAY(count)=INIT_VALUE
        ENDDO !count
      ELSE
        DO count=1,NEW_LENGTH
          NEW_ARRAY(count)=OLD_ARRAY(count)
        ENDDO !count
      ENDIF

      CALL EXITS('COPY_AND_INIT_INT')
      RETURN
 9999 CALL ERRORS('COPY_AND_INIT_INT',ERROR)
      CALL EXITS('COPY_AND_INIT_INT')
      RETURN 1
      END


