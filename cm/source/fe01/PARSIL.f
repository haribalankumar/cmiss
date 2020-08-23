      SUBROUTINE PARSIL(STRING,NXIL,NTIL,IL,ERROR,*)

C#### Subroutine: PARSIL
C###  Description:
C###    PARSIL converts the character string STRING
C###    into a list of NTIL integer variables IL.
C**** The list may be composed of items separated by commas/spaces
C**** and/or items in an implied do list.
C**** An implied do list is structured as:
C**** ISTART..IFINISH:IINCR
C**** where ISTART and IFINISH represent the lower and upper bounds
C**** of the list while IINCR is the increment.
C**** If IINCR is omitted a value of 1 is assumed.
C**** Rewritten by MPN 24-May-95

      IMPLICIT NONE
!     Parameter List
      INTEGER NTIL,NXIL,IL(NXIL)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER CLOCAT,ERR,IBEG,IEND,IFINISH,IINCR,
     '  IPOS1,IPOS2,IPOS3,IPOS4,IPOS_comma,IPOS_space,
     '  ISTART,ival
      LOGICAL CONTINUE

      CALL ENTERS('PARSIL',*9999)
      NTIL=0
C KAT 10Oct98: not necessary
C      DO noil=1,NXIL !initialise list
C        IL(noil)=0
C      ENDDO
      IF(STRING.NE.' ') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        IPOS1=IBEG
        DO WHILE(IPOS1.LE.IEND)
          IPOS_comma=CLOCAT(',',STRING(IPOS1:IEND))
          IPOS_space=CLOCAT(' ',STRING(IPOS1:IEND))
          IF(IPOS_comma.EQ.0.AND.IPOS_space.EQ.0) THEN
            IPOS2=IEND !no comma/space so last entry in list
          ELSE IF(IPOS_comma.EQ.0) THEN
            IPOS2=IPOS1+IPOS_space-2
          ELSE IF(IPOS_space.EQ.0) THEN
            IPOS2=IPOS1+IPOS_comma-2
          ELSE
            IPOS2=IPOS1+MIN(IPOS_comma,IPOS_space)-2
          ENDIF
          IPOS3=IPOS1+CLOCAT('..',STRING(IPOS1:IPOS2))
          IF(IPOS3.EQ.IPOS1) THEN !single entry
C LKC 16-NOV-1998 Change to allow exact storage
C            IF(NTIL.GE.NXIL) THEN

C LKC 14-MAY-1999 Move below to grab errors
C            IF(NTIL.GT.NXIL) THEN
C              ERROR=' Number of items exceeds storage'
C              GOTO 9999
C            ENDIF
            NTIL=NTIL+1
            IF(NTIL.GT.NXIL) THEN
              ERROR=' Number of items exceeds storage'
              GOTO 9999
            ENDIF

            CALL INTFROMCHAR(IL(NTIL),STRING(IPOS1:IPOS2),ERR)
            IF(ERR.NE.0) GOTO 9998
          ELSE !list entry
            IPOS4=IPOS1+CLOCAT(':',STRING(IPOS1:IPOS2))
            IF(IPOS4.EQ.IPOS1) THEN !no increment entered
              IINCR=1
              IPOS4=IPOS2
            ELSE !set the increment
              CALL INTFROMCHAR(IINCR,STRING(IPOS4:IPOS2),ERR)
              IF(ERR.NE.0) GOTO 9998
              IPOS4=IPOS4-2
            ENDIF
            CALL INTFROMCHAR(ISTART,STRING(IPOS1:IPOS3-2),ERR)
            IF(ERR.NE.0) GOTO 9998
            CALL INTFROMCHAR(IFINISH,STRING(IPOS3+1:IPOS4),ERR)
            IF(ERR.NE.0) GOTO 9998
            DO ival=ISTART,IFINISH,IINCR
C LKC 16-NOV-1998 Change to allow exact storage
C              IF(NTIL.GE.NXIL) THEN
C LKC 14-MAY-1999 Move below to grab errors
C              IF(NTIL.GT.NXIL) THEN
C                ERROR=' Number of items exceeds storage'
C                GOTO 9999
C              ENDIF
              NTIL=NTIL+1
              IF(NTIL.GT.NXIL) THEN
                ERROR=' Number of items exceeds storage'
                GOTO 9999
              ENDIF
              IL(NTIL)=ival
            ENDDO !ival
          ENDIF
          IPOS1=IPOS2+2 !find position of next entry
          CONTINUE=.TRUE.
          DO WHILE(IPOS1.LE.IEND.AND.CONTINUE)
            IF(STRING(IPOS1:IPOS1).EQ.' '.OR.
     '        STRING(IPOS1:IPOS1).EQ.',') THEN !trim leading blanks/commas
              IPOS1=IPOS1+1
            ELSE
              CONTINUE=.FALSE.
            ENDIF
          ENDDO
        ENDDO !while(IPOS.LE.IEND)
      ENDIF

      CALL EXITS('PARSIL')
      RETURN
 9998 ERROR=' '
 9999 CALL ERRORS('PARSIL',ERROR)
      CALL EXITS('PARSIL')
      RETURN 1
      END


