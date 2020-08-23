      SUBROUTINE SET_NUM_THREADS(STRING,ERROR,*)

C#### Subroutine: SET_NUM_THREADS
C###  Description:
C###    If allowed, sets the maximum number of threads to use. Can be
C###    used more than once to dynamically change the number of threads
C###    used during a problem.

C *** Created : DPN 03 July 1999

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,NUM_THREADS
      LOGICAL CBBREV
C$    INTEGER OMP_GET_MAX_THREADS
C$    LOGICAL OMP_IN_PARALLEL
C$    EXTERNAL OMP_IN_PARALLEL,OMP_GET_MAX_THREADS

      CALL ENTERS('SET_NUM_THREADS',*9999)
      IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C#### Command: SET num_threads <#[2]>
C###  Parameter:      <#[2]>
C###    Specify the maximum number of threads to use.
C###  Description:
C###    If allowed, sets the maximum number of threads to use. Can be
C###    used more than once to dynamically change the number of threads
C###    used during a problem. If the environment variable
C###    OMP_DYNAMIC is set to false, the the maximum number of threads
C###    will always be used.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<#>[2]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('fe01','doc','SET_NUM_THREADS',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'NUM_THREADS',4,noco,NTCO,N3CO)) THEN
          IF(NTCO.GE.N3CO+1) THEN
            NUM_THREADS=IFROMC(CO(N3CO+1))
            CALL ASSERT(NUM_THREADS.GT.0,
     '        'Number of threads must be > 0',ERROR,*9999)
          ELSE
            NUM_THREADS=2
          ENDIF
C ***     Should only set the number of threads when not in a region
C         executing in parallel!!
C$        IF(OMP_IN_PARALLEL()) THEN
C$          ERROR='Cannot set num_threads in parallel'
C$          GOTO 9999
C$        ELSE
C$          CALL OMP_SET_NUM_THREADS(NUM_THREADS)
C$        ENDIF
        ENDIF
      ENDIF

      CALL EXITS('SET_NUM_THREADS')
      RETURN
 9999 CALL ERRORS('SET_NUM_THREADS',ERROR)
      CALL EXITS('SET_NUM_THREADS')
      RETURN 1
      END


