
      SUBROUTINE ALLOCATE_MEMORY(N,NMIN,TYPE,PTR,MEM_INIT,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER N,NMIN,TYPE,MEM_INIT
      POINTER PTR
      CHARACTER*(*) ERROR
!     Local variables
      POINTER NMEM,NSIZE
!     Functions
      POINTER FMALLOC
      EXTERNAL FMALLOC


      CALL ASSERT(N.GT.0,'>>Trying to allocate a negative amount of '
     '  //'memory',ERROR,*9999)

      IF(TYPE.EQ.DPTYPE) THEN
        NSIZE=8
      ELSE IF(TYPE.EQ.CHARTYPE) THEN
        NSIZE=1
      ELSE
        NSIZE=4
      ENDIF
      NMEM=N
      
      PTR=FMALLOC(NMEM,NSIZE,MEM_INIT)
      IF(PTR.LE.0) THEN
        ERROR='>>Error allocating memory'
        GOTO 9999
      ENDIF

      RETURN

 9999 CONTINUE
      RETURN 1
      END


      SUBROUTINE FREE_MEMORY(PTR,ERROR,*)

      IMPLICIT NONE
      POINTER PTR
      CHARACTER*(*) ERROR

      IF(PTR.NE.0) THEN
        CALL FFREE(PTR)
        PTR=0
      ELSE
        ERROR='>>ERROR: pointer to free is NULL'
        RETURN 1
      ENDIF

      RETURN
      END

      SUBROUTINE ASSERT(LOG,BLURB,ERROR,*)

      IMPLICIT NONE
!     Routine arguments
      LOGICAL LOG
      CHARACTER*(*) BLURB,ERROR

      IF(.NOT.LOG) THEN
        ERROR=BLURB
        RETURN 1
      ENDIF

      RETURN
      END

      SUBROUTINE WRITES(IOOP,OP_STRING,ERROR,*)

      IMPLICIT NONE
!     Routine arguments
      INTEGER IOOP
      CHARACTER*(*) OP_STRING,ERROR
!     Local variables
      INTEGER IERR
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      WRITE(*,'(A)',IOSTAT=IERR)OP_STRING(:LEN_TRIM(OP_STRING))
      IF(IERR.NE.0) THEN
        ERROR='>>Error writing to stdout'
        RETURN 1
      ENDIF

      RETURN
      END

      SUBROUTINE ERRORS(ROUTINE,ERROR)

      IMPLICIT NONE
      CHARACTER*(*) ROUTINE,ERROR

      INTEGER LENGTH
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM

      LENGTH=LEN_TRIM(ERROR)
      IF(LENGTH.GT.0) THEN
        WRITE(0,'(3A)') ROUTINE,': ',ERROR(:LENGTH)
        ERROR=' '
      ELSE
        WRITE(0,'(2A)') '>>',ROUTINE
      ENDIF

      RETURN
      END

      INTEGER FUNCTION LEN_TRIM(STRING)

      IMPLICIT NONE
      CHARACTER*(*) STRING

      INTEGER LNBLNK
      EXTERNAL LNBLNK

      LEN_TRIM=MAX(1,LNBLNK(STRING))

      RETURN
      END
