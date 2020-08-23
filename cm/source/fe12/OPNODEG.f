      SUBROUTINE OPNODEG(ERROR,*)

C#### Subroutine: OPNODEG
C###  Description:
C###    OPNODEG outputs node groups.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'grou00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER N,nogrno,nog,NOGT
!     Functions
      INTEGER ILISTMBR

      CALL ENTERS('OPNODEG',*9999)

      WRITE(OP_STRING,'('' Number of node groups = '',I2)') NTGRNO
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nogrno=1,NTGRNO
C TVK 05July1999: Dynamic Groups
C        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'' Total= '','
C     '    //'I5,'' :'')') nogrno,LAGRNO(nogrno),LIGRNO(0,nogrno)
        NOGT=NLIGRNO(nogrno)
        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'' Total= '','
     '    //'I6,'' :'')') nogrno,LAGRNO(nogrno),NOGT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C TVK 05July1999: Dynamic Groups
C        WRITE(OP_STRING,'((20I5))')
C     '    (LIGRNO(N,nogrno),N=1,LIGRNO(0,nogrno))
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C      ENDDO
      DO nog=1,NOGT/10+1
        WRITE(OP_STRING,'((10(1X,I8)))')
     '      (ILISTMBR(%VAL(LIGRNO_PTR(nogrno)),n),
     '      n=(nog-1)*10+1,MIN(nog*10,NOGT))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO

      CALL EXITS('OPNODEG')
      RETURN
 9999 CALL ERRORS('OPNODEG',ERROR)
      CALL EXITS('OPNODEG')
      RETURN 1
      END


