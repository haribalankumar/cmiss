      SUBROUTINE APPENDC(IEND,STRING1,STRING2)

C#### Subroutine: APPENDC
C###  Description:
C###    APPENDC writes a fortran character STRING1 into fortran
C###    character STRING2 after position IEND and sets IEND to the last
C###    character of the copied string in STRING2.

      IMPLICIT NONE
!     Parameter List
      INTEGER IEND
      CHARACTER STRING1*(*),STRING2*(*)
!     Local Variables
      INTEGER IBEG,LENGTH1,LENGTH2

C      CALL ENTERS('APPENDC',*9999)

      LENGTH2=LEN(STRING2)
      IF(IEND.LT.LENGTH2) THEN
        LENGTH1=LEN(STRING1)
        IBEG=IEND+1
        IEND=MIN(IEND+LENGTH1,LENGTH2)
        STRING2(IBEG:IEND)=STRING1
      ENDIF

C      CALL EXITS('APPENDC')
      RETURN
      END


