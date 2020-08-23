      LOGICAL FUNCTION CBBREV(CO,SUBSTR,MNCH,N1CO,N2CO,N3CO)

C#### Function: CBBREV
C###  Type: LOGICAL
C###  Description:
C###    CBBREV returns .TRUE. if CO(N1CO..N2CO) contains an
C###    abbreviation of SUBSTR which is at least MNCH characters long.
C###    As a side effect, if CBBREV is .TRUE., N3CO is position in
C###    CO of SUBSTR.

      IMPLICIT NONE
!     Parameter List
      INTEGER MNCH,N1CO,N2CO,N3CO
      CHARACTER CO(*)*(*),SUBSTR*(*)
!     Local Variables
      INTEGER noco_local
      LOGICAL ABBREV

      CBBREV=.FALSE.
      DO noco_local=N1CO,N2CO
        IF(ABBREV(CO(noco_local),SUBSTR,MNCH)) THEN
          N3CO=noco_local
          CBBREV=.TRUE.
          RETURN
        ENDIF
      ENDDO

      RETURN
      END


