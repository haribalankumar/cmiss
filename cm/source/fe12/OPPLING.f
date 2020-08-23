      SUBROUTINE OPPLING(ERROR,*)

C#### Subroutine: OPPLING
C###  Description:
C###    OPPLING outputs polyline groups.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'grou00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n,nogrpl

      CALL ENTERS('OPPLING',*9999)

      WRITE(OP_STRING,'('' Number of polyline groups = '',I2)')NTGRPL
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nogrpl=1,NTGRPL
        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'' Total= '','
     '    //'I5,'' :'')') nogrpl,LAGRPL(nogrpl),LIGRPL(0,nogrpl)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'((20I5))')
     '    (LIGRPL(n,nogrpl),n=1,LIGRPL(0,nogrpl))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO

      CALL EXITS('OPPLING')
      RETURN
 9999 CALL ERRORS('OPPLING',ERROR)
      CALL EXITS('OPPLING')
      RETURN 1
      END


