      SUBROUTINE PARSE_DATA(nd0,nd1,noco,NTCO,CO,ERROR,*)

C#### Subroutine: PARSE_DATA
C###  Description:
C###    PARSE_DATA parses command CO for keyword 'data'
C###    and returns first and last data point numbers as nd0,nd1

C!!!  PARSILG should be extended to handle DATA pts; this routine could
C!!!  then call PARSILG.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'

!     Parameter List
      INTEGER nd0,nd1,noco,NTCO
      CHARACTER CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,nogrda,n1grda
      LOGICAL CBBREV,CONTINUE
      CHARACTER LABEL*30

      CALL ENTERS('PARSE_DATA',*9999)
      n1grda=0
      IF(CBBREV(CO,'DATA',2,NOCO+1,NTCO,N3CO)) THEN
        LABEL=CO(N3CO+1)
        CALL STRING_TRIM(LABEL,IBEG,IEND)
        CONTINUE=.TRUE.
        nogrda=0
        DO WHILE(CONTINUE.AND.nogrda.LT.NTGRDA)
          nogrda=nogrda+1
          IF(LAGRDA(nogrda)(IBEG:IEND).EQ.LABEL(IBEG:IEND)) THEN
            n1grda=nogrda
            CONTINUE=.FALSE.
          ENDIF
        ENDDO !while continue
        IF(n1grda.EQ.0) THEN
          ERROR=' >>Group name '//LABEL(IBEG:IEND)//' not defined'
          GO TO 9999
        ENDIF
      ELSE IF(NTGRDA.GT.0) THEN
        n1grda=1
      ENDIF

      IF(n1grda.GT.0) THEN
        nd0=LIGRDA(1,n1grda)
        nd1=LIGRDA(2,n1grda)
      ELSE
        nd0=1
        nd1=NDT
      ENDIF

      CALL EXITS('PARSE_DATA')
      RETURN
 9999 CALL ERRORS('PARSE_DATA',ERROR)
      CALL EXITS('PARSE_DATA')
      RETURN 1
      END


