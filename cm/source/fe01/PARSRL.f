      SUBROUTINE PARSRL(STRING,NXRL,NTRL,RL,ERROR,*)

C#### Subroutine: PARSRL
C###  Description:
C###    PARSRL converts the character string STRING into a list of NTRL
C###    REAL*8 variables RL. The list may be composed of single items
C###    and/or items in an implied do list separated by commas and/or
C###    spaces. An implied do list is structured as:
C###    RSTART..RFINISH:RINCR where RSTART and RFINISH represent the
C###    lower and upper bounds of the list while RINCR is the increment.
C###    If RINCR is omitted a value of 1.0d0 is assumed.

C**** Rewritten by MPN 2-Aug-95

      IMPLICIT NONE
!     Parameter List
      INTEGER NTRL,NXRL
      REAL*8 RL(NXRL)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER CLOCAT,IBEG,IEND,
     '  IPOS1,IPOS2,IPOS3,IPOS4,IPOS_comma,IPOS_space
      REAL*8 RFINISH,RINCR,RSTART,rval
      LOGICAL CONTINUE

      CALL ENTERS('PARSRL',*9999)
      NTRL=0
C KAT 10Oct98: not necessary
C      DO norl=1,NXRL !initialise list
C        RL(norl)=0.0d0
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
            IF(NTRL.GE.NXRL) THEN
              ERROR=' Number of items exceeds storage'
              GOTO 9999
            ENDIF
            NTRL=NTRL+1
            CALL PARSRE(STRING(IPOS1:IPOS2),RL(NTRL),ERROR,*9999)
          ELSE !list entry
            IPOS4=IPOS1+CLOCAT(':',STRING(IPOS1:IPOS2))
            IF(IPOS4.EQ.IPOS1) THEN !no increment entered
              RINCR=1.0d0
              IPOS4=IPOS2
            ELSE !set the increment
              CALL PARSRE(STRING(IPOS4:IPOS2),RINCR,ERROR,*9999)
              IPOS4=IPOS4-2
            ENDIF
            CALL PARSRE(STRING(IPOS1:IPOS3-2),RSTART,ERROR,*9999)
            CALL PARSRE(STRING(IPOS3+1:IPOS4),RFINISH,ERROR,*9999)
            DO rval=RSTART,RFINISH,RINCR
              IF(NTRL.GE.NXRL) THEN
                ERROR='>>Number of items exceeds storage'
                GOTO 9999
              ENDIF
              NTRL=NTRL+1
              RL(NTRL)=rval
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

      CALL EXITS('PARSRL')
      RETURN
 9999 CALL ERRORS('PARSRL',ERROR)
      CALL EXITS('PARSRL')
      RETURN 1
      END


