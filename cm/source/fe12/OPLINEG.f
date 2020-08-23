      SUBROUTINE OPLINEG(ERROR,*)

C#### Subroutine: OPLINEG
C###  Description:
C###    OPELING outputs line groups.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'grou00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n,nogrli,nog,NOGT
!     Functions
      INTEGER ILISTMBR

      CALL ENTERS('OPLINEG',*9999)
C TVK 28July1999: Dynamic Groups (whole routine is new)
      WRITE(OP_STRING,'('' Number of line groups = '',I2)') NTGRLI
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nogrli=1,NTGRLI
        NOGT=NLIGRLI(nogrli)
        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'' Total= '','
     '    //'I6,'' :'')') nogrli,LAGRLI(nogrli),NOGT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nog=1,NOGT/10+1
        WRITE(OP_STRING,'((10(1X,I8)))')
     '      (ILISTMBR(%VAL(LIGRLI_PTR(nogrli)),n),
     '      n=(nog-1)*10+1,MIN(nog*10,NOGT))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO

      CALL EXITS('OPLINEG')
      RETURN
 9999 CALL ERRORS('OPLINEG',ERROR)
      CALL EXITS('OPLINEG')
      RETURN 1
      END


