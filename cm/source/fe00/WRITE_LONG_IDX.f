      SUBROUTINE WRITE_LONG_IDX(DATA_TYPE,DELTA,INDEX,IUNIT,NUM_FIRST,
     '  NUM_REPEAT,VECTOR_INT,VECTOR_DP,FMT_FIRST,FMT_REPEAT,ERROR,*)

C#### Subroutine: WRITE_LONG_IDX
C###  Description:
C###    WRITE_LONG_IDX writes a VECTOR indexed by INDEX to the given
C###    IUNIT.  The first line format is initially used, followed by
C###    the second line format - repeated as many times as necessary.
C###    NUM_FIRST is number of data items in first FMT
C###    NUM_REPEAT is number of data items in repeat FMT
C###    INDEX(0) is the number of items, INDEX(i) are the indicies of
C###    the vector.
C###  See-Also: WRITE_LONG

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER DATA_TYPE,DELTA,INDEX(0:*),IUNIT,NUM_FIRST,NUM_REPEAT,
     '  VECTOR_INT(*)
      REAL*8 VECTOR_DP(*)
      CHARACTER ERROR*(*),FMT_FIRST*(*),FMT_REPEAT*(*)
!     Local Variables
      INTEGER current,count,NUM_TO_DO

      CALL ENTERS('WRITE_LONG_IDX',*9999)

      CALL ASSERT((DATA_TYPE.EQ.INTTYPE).OR.(DATA_TYPE.EQ.DPTYPE),
     '  '>>Data type not implemented',ERROR,*9999)
      NUM_TO_DO=INDEX(0)
      IF(DATA_TYPE.EQ.INTTYPE) THEN
        WRITE(OP_STRING,FMT_FIRST)
     '    (VECTOR_INT((INDEX(count)-1)*DELTA+1),count=1,
     '    MIN(NUM_FIRST,INDEX(0)))
      ELSEIF(DATA_TYPE.EQ.DPTYPE) THEN
        WRITE(OP_STRING,FMT_FIRST)
     '    (VECTOR_DP((INDEX(count)-1)*DELTA+1),count=1,
     '    MIN(NUM_FIRST,INDEX(0)))
      ENDIF !DATA_TYPE
      CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
      NUM_TO_DO=INDEX(0)-NUM_FIRST
      current=NUM_FIRST+1
      DO WHILE(NUM_TO_DO.GT.0) !more stuff to do
        IF(DATA_TYPE.EQ.INTTYPE) THEN
          WRITE(OP_STRING,FMT_REPEAT)
     '      (VECTOR_INT((INDEX(count)-1)*DELTA+1),count=current,
     '      MIN(current+NUM_REPEAT-1,INDEX(0)))
        ELSEIF(DATA_TYPE.EQ.DPTYPE) THEN
          WRITE(OP_STRING,FMT_REPEAT)
     '      (VECTOR_DP((INDEX(count)-1)*DELTA+1),count=current,
     '      MIN(current+NUM_REPEAT-1,INDEX(0)))
        ENDIF !DATA_TYPE
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        current=current+NUM_REPEAT
        NUM_TO_DO=NUM_TO_DO-NUM_REPEAT
      ENDDO !NUM_TO_DO > 0

      CALL EXITS('WRITE_LONG_IDX')
      RETURN
 9999 CALL ERRORS('WRITE_LONG_IDX',ERROR)
      CALL EXITS('WRITE_LONG_IDX')
      RETURN 1
      END


