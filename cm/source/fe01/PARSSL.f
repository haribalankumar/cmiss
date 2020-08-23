      SUBROUTINE PARSSL(STRING,NXSL,NTSL,SL,ERROR,*)

C#### Subroutine: PARSSL
C###  Description:
C###    PARSSL converts the character string STRING into a list of
C###    NTSL STRING variables SL. The list may be composed of items
C###    separated by commas.

      IMPLICIT NONE
!     Parameter List
      INTEGER NTSL,NXSL
      CHARACTER ERROR*(*),SL(NXSL)*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IEND,n1ch,n2ch
      LOGICAL FOUND_END

      CALL ENTERS('PARSSL',*9999)
      CALL STRING_TRIM(STRING,IBEG,IEND)
      n1ch=IBEG-1
      NTSL=0
      IF(STRING.NE.' ') THEN
        DO 2 n2ch=IBEG,IEND+1
          FOUND_END=n2ch.GT.IEND
          IF(.NOT.FOUND_END) FOUND_END=STRING(n2ch:n2ch).EQ.','
          IF(FOUND_END) THEN
            IF(NTSL.GE.NXSL) THEN
              ERROR='>>Number of items exceeds storage'
              GOTO 9999
            ENDIF
            NTSL=NTSL+1
            SL(NTSL)=STRING(n1ch+1:n2ch-1)
            n1ch=n2ch
          ENDIF
 2      CONTINUE
      ENDIF

      CALL EXITS('PARSSL')
      RETURN
 9999 CALL ERRORS('PARSSL',ERROR)
      CALL EXITS('PARSSL')
      RETURN 1
      END


