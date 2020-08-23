      SUBROUTINE OPGRIDG(ERROR,*)

C#### Subroutine: OPGRIDG
C###  Description:
C###    OPGRIDG outputs grid point groups.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'grou00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n,nogrgr,nog,NOGT
!     Functions
      INTEGER ILISTMBR

      CALL ENTERS('OPGRIDG',*9999)

      WRITE(OP_STRING,'('' Number of grid groups = '',I2)') NTGRGR
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nogrgr=1,NTGRGR
C KAT 11May99: dynamic groups
C        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'' Total= '','
C     '    //'I6,'' :'')') nogrgr,LAGRGR(nogrgr),LIGRGR(0,nogrgr)
        NOGT=NLIGRGR(nogrgr)
        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'' Total= '','
     '    //'I6,'' :'')') nogrgr,LAGRGR(nogrgr),NOGT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C KAT 11May99: dynamic groups
C        NOGT=LIGRGR(0,nogrgr)
C        DO nog=1,INT(NOGT/10)+1
C          WRITE(OP_STRING,'((10(1X,I6)))')
C     '      (LIGRGR(n,nogrgr),n=(nog-1)*10+1,MIN(nog*10,NOGT))
C          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        ENDDO
        DO nog=1,NOGT/10+1
          WRITE(OP_STRING,'((10(1X,I8)))')
     '      (ILISTMBR(%VAL(LIGRGR_PTR(nogrgr)),n),
     '      n=(nog-1)*10+1,MIN(nog*10,NOGT))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO

      CALL EXITS('OPGRIDG')
      RETURN
 9999 CALL ERRORS('OPGRIDG',ERROR)
      CALL EXITS('OPGRIDG')
      RETURN 1
      END


