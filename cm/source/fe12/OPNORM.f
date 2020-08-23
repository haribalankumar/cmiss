      SUBROUTINE OPNORM(NEELEM,NRLIST,NW,ERROR,*)

C#### Subroutine: OPNORM
C###  Description:
C###    OPNORM outputs normal reversals for BEM problems

      IMPLICIT NONE

      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NRLIST(0:NRM),NW(NEM,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,ne,nee,nr,nrr
      CHARACTER CHAR1*6

      CALL ENTERS('OPNORM',*9999)

      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr)
        WRITE(CHAR1,'(I2)') nr
        CALL STRING_TRIM(CHAR1,IBEG,IEND)
        WRITE(OP_STRING,'(''  Region '//CHAR1(IBEG:IEND)//': '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        DO nee=1,NEELEM(0,nr)
          ne=NEELEM(nee,nr)
          WRITE(CHAR1,'(I6)') ne
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
          IF(NW(ne,3).EQ.0) THEN
            WRITE(OP_STRING,'(''  Element '//CHAR1(IBEG:IEND)
     '        //' normal is Standard '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSEIF(NW(ne,3).EQ.1) THEN
            WRITE(OP_STRING,'(''  Element '//CHAR1(IBEG:IEND)
     '        //' normal is Reversed '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'(''  Element '//CHAR1(IBEG:IEND)
     '        //' normal is incorrectly set up '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
      ENDDO

      CALL EXITS('OPNORM')
      RETURN
 9999 CALL ERRORS('OPNORM',ERROR)
      CALL EXITS('OPNORM')
      RETURN 1
      END


