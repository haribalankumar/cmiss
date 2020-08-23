      SUBROUTINE WRITES(id,OP_STRING,ERROR,*)

C#### Subroutine: WRITES
C###  Description:
C###    <HTML>
C###    WRITES writes strings to output identified by id as follows:
C###    <PRE>
C###    id=2 (=IOOP) listing
C###      3  (=IODI) diagnostics
C###      4  (=IOTR) trace
C###      5  (=IOER) errors
C###      6  (=IOH1) ? help
C###      7  (=IOH2) ?? help
C###      8  (=IOH3) ??? help
C###      9
C###    id>9 is used for file output
C###    </PRE></HTML>
C**** Note: WRITES cannot call Enters or Exits (since would get into
C**** an infinite loop).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'disp00.cmn'
      INCLUDE 'fsklib.inc'
C$    INCLUDE 'openmp.inc'
C$    INCLUDE 'mp00.cmn'
!     Parameter List
      INTEGER id
      CHARACTER ERROR*(*),OP_STRING(100)*(255)
!     Local Variables
      INTEGER i,IEND,j,NUM_BLANKS,NUM_RECORDS
!     Functions
      INTEGER LEN_TRIM1

C$    IF((OMP_GET_THREAD_NUM().EQ.THREAD_NUM).OR.
C$   &   (THREAD_NUM.EQ.-1)) THEN
C ***   Calculate number of records in OP_STRING
        IF(OS_TYPE(1:3).EQ.'VMS') THEN
          i=1
          DO WHILE (OP_STRING(i)(1:1).NE.CHAR(0).AND.I.LT.100)
          i=i+1
          ENDDO
          NUM_RECORDS=i-1
          IF(i.EQ.100) NUM_RECORDS=100
        ELSE IF(OS_TYPE(1:4).EQ.'UNIX') THEN
          i=1
          NUM_BLANKS=0
          DO WHILE (i.LT.100.AND.NUM_BLANKS.LT.2)
            i=i+1
            IF(LEN_TRIM1(OP_STRING(i)).EQ.1) THEN
              NUM_BLANKS=NUM_BLANKS+1
            ELSE
              NUM_BLANKS=0
            ENDIF
          ENDDO
          NUM_RECORDS=i-NUM_BLANKS
          IF(i.EQ.100) NUM_RECORDS=100
        ENDIF

        DO i=1,NUM_RECORDS
C         Only trimming to min length of 1 for FORTRAN 77
          IEND=LEN_TRIM1(OP_STRING(i))
          CALL WRITE_LINE(id,OP_STRING(i)(:IEND),ERROR,*9999)
        ENDDO

C *** Reset OP_STRING
        IF(OS_TYPE(1:3).EQ.'VMS') THEN
          DO i=1,NUM_RECORDS
            OP_STRING(i)(1:1)=CHAR(0)
          ENDDO
        ELSE IF(OS_TYPE(1:4).EQ.'UNIX') THEN
          DO i=1,NUM_RECORDS
            DO j=1,MXCH
              OP_STRING(i)(j:j)=' '
            ENDDO
          ENDDO
        ENDIF

C$    ENDIF
      RETURN
 9999 CALL ERRORS('WRITES',ERROR)
      RETURN 1
      END


