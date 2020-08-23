      SUBROUTINE OPELEMG(ERROR,*)

C#### Subroutine: OPELEMG
C###  Description:
C###    OPELEMG outputs element groups.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'grou00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n,nogrel,nog,NOGT
!     Functions
      INTEGER ILISTMBR

      CALL ENTERS('OPELEMG',*9999)

      WRITE(OP_STRING,'('' Number of element groups = '',I2)') NTGREL
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nogrel=1,NTGREL
C TVK 14July1999: Dynamic Groups
C        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'' Total= '','
C     '    //'I5,'' :'')') nogrno,LAGRNO(nogrno),LIGRNO(0,nogrno)
        NOGT=NLIGREL(nogrel)
        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'' Total= '','
     '    //'I6,'' :'')') nogrel,LAGREl(nogrel),NOGT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C TVK 14July1999: Dynamic Groups
C        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'
C     '    //''' Total= '',I5,'' :'')')
C     '    nogrel,LAGREL(nogrel),LIGREL(0,nogrel)
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'((20I5))') (LIGREL(n,nogrel),
C     '    n=1,LIGREL(0,nogrel))
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C      ENDDO
      DO nog=1,NOGT/10+1
        WRITE(OP_STRING,'((10(1X,I8)))')
     '      (ILISTMBR(%VAL(LIGREL_PTR(nogrel)),n),
     '      n=(nog-1)*10+1,MIN(nog*10,NOGT))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO

      CALL EXITS('OPELEMG')
      RETURN
 9999 CALL ERRORS('OPELEMG',ERROR)
      CALL EXITS('OPELEMG')
      RETURN 1
      END


