      SUBROUTINE USER(STRING,DE,ERROR,*)

C#### Subroutine: USER
C###  Description:
C###    USER substitutes user defined names within STRING, separated
C###    by delimiters found in DE, by the strings they represent.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'user00.cmn'
!     Parameter List
      CHARACTER DE*(*),ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER CLOCAT,IBEG,IEND,l0,l1,l2,nous
      LOGICAL CELEM

      CALL ENTERS('USER',*9999)
      DO 2 nous=1,NTUS
        IF(ACUS(nous)) THEN
          l0=1
C new MPN 28Jul97: bug fix
          CALL STRING_TRIM(STRING,IBEG,IEND)
          DO 1 l1=IBEG,IEND
c old          DO 1 l1=1,LEN(STRING)
            IF(STRING(l1:l1+2).EQ.'   ') GO TO 2
            l2=CLOCAT(USID(nous)(:LNUSID(nous)),STRING(l0:))
            l0=l0+l2-1
            IF(l2.LE.0) THEN
              GOTO 2
            ELSE IF(l0.EQ.1) THEN
              IF(CELEM(STRING(l0+LNUSID(nous):l0+LNUSID(nous)),DE)) THEN
                STRING=US(nous)(:LNUS(nous))
     '            //STRING(l0+LNUSID(nous):)
                l0=l0+LNUS(nous)-1
              ENDIF
            ELSE IF(l0.GT.1) THEN
              IF(CELEM(STRING(l0-1:l0-1),DE).AND.
     '          CELEM(STRING(l0+LNUSID(nous):l0+LNUSID(nous)),DE)) THEN
                STRING=STRING(:l0-1)
     '            //US(nous)(:LNUS(nous))
     '            //STRING(l0+LNUSID(nous):)
                l0=l0+LNUS(nous)-1
              ENDIF
              l0=l0+1
            ENDIF
 1        CONTINUE
        ENDIF
 2    CONTINUE

      CALL EXITS('USER')
      RETURN
 9999 CALL ERRORS('USER',ERROR)
      CALL EXITS('USER')
      RETURN 1
      END


