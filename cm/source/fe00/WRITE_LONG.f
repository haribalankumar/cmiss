      SUBROUTINE WRITE_LONG(DATA_TYPE,FIRST,DELTA,IUNIT,LAST,NUM_FIRST,
     '  NUM_REPEAT,VECTOR_INT,VECTOR_DP,FMT_FIRST,FMT_REPEAT,ERROR,*)

C#### Subroutine: WRITE_LONG
C###  Description:
C###    WRITE_LONG writes VECTOR to the given IUNIT.  The first
C###    line format is initially used, followed by the second line
C###    format - repeated as many times as necessary.
C###    NUM_FIRST is number of data items in first FMT
C###    NUM_REPEAT is number of data items in repeat FMT
C###    FIRST and LAST are the extents of the data
C###  See-Also: WRITE_LONG_IDX

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER DATA_TYPE,FIRST,DELTA,IUNIT,LAST,
     '  NUM_FIRST,NUM_REPEAT,VECTOR_INT(*)
      REAL*8 VECTOR_DP(*)
      CHARACTER ERROR*(*),FMT_FIRST*(*),FMT_REPEAT*(*)
!     Local Variables
      INTEGER current,final,count

      CALL ENTERS('WRITE_LONG',*9999)
      CALL ASSERT((DATA_TYPE.EQ.INTTYPE).OR.(DATA_TYPE.EQ.DPTYPE),
     '  '>>Data type not implemented',ERROR,*9999)
      current=FIRST
      final=current+(NUM_FIRST-1)*DELTA
      IF(final.GT.LAST) final=LAST
      IF(DATA_TYPE.EQ.INTTYPE) THEN
        WRITE(OP_STRING,FMT_FIRST)
     '    (VECTOR_INT(count),count=current,final,DELTA)
      ELSEIF(DATA_TYPE.EQ.DPTYPE) THEN
        WRITE(OP_STRING,FMT_FIRST)
     '    (VECTOR_DP(count),count=current,final,DELTA)
      ENDIF !DATA_TYPE
      CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
      DO WHILE(final.LT.LAST) !more stuff to do
        current=final+DELTA
        final=final+NUM_REPEAT*DELTA
        IF(final.GT.LAST) final=LAST
        IF(DATA_TYPE.EQ.INTTYPE) THEN
          WRITE(OP_STRING,FMT_REPEAT)
     '      (VECTOR_INT(count),count=current,final,DELTA)
        ELSEIF(DATA_TYPE.EQ.DPTYPE) THEN
          WRITE(OP_STRING,FMT_REPEAT)
     '      (VECTOR_DP(count),count=current,final,DELTA)
        ENDIF !DATA_TYPE
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
      ENDDO !final.LT.LAST

      CALL EXITS('WRITE_LONG')
      RETURN
 9999 CALL ERRORS('WRITE_LONG',ERROR)
      CALL EXITS('WRITE_LONG')
      RETURN 1
      END


