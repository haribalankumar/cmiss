      SUBROUTINE OPFACEG(ERROR,*)

C#### Subroutine: OPFACEG
C###  Description:
C###    OPFACEG outputs face groups.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'grou00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n,nogrfa,nog,NOGT
!     Functions
      INTEGER ILISTMBR

      CALL ENTERS('OPFACEG',*9999)
C TVK 09July1999: Dynamic Groups (whole routine is new)
      WRITE(OP_STRING,'('' Number of face groups = '',I2)') NTGRFA
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nogrfa=1,NTGRFA
        NOGT=NLIGRFA(nogrfa)
        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'' Total= '','
     '    //'I6,'' :'')') nogrfa,LAGRFA(nogrfa),NOGT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nog=1,NOGT/10+1
        WRITE(OP_STRING,'((10(1X,I8)))')
     '      (ILISTMBR(%VAL(LIGRFA_PTR(nogrfa)),n),
     '      n=(nog-1)*10+1,MIN(nog*10,NOGT))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO

      CALL EXITS('OPFACEG')
      RETURN
 9999 CALL ERRORS('OPFACEG',ERROR)
      CALL EXITS('OPFACEG')
      RETURN 1
      END


