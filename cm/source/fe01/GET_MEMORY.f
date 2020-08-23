      SUBROUTINE GET_MEMORY(N,ITEM_TYPE,ADDR_PTR,INITIALISE,ERROR,*)

C#### Subroutine: GET_MEMORY
C###  Description:
C###    GET_MEMORY dynamically allocates N items of type ITEM_TYPE
C###    in memory. ITEM_TYPES are 1/2/3/4/5/6/7/8 for INTEGER,REAL,
C###    REAL*8,CHARACTER,LOGICAL,INTEGER*2,COMPLEX and COMPLEX*16.
C###    Constants are setup in mach00.cmn for these types. If
C###    INITIALISE is true the variables are initialised to a
C###    special value (currently only for INTEGER's, REAL's and
C###    REAL*8's).
C###  See-Also: ALLOCATE_MEMORY,FREE_MEMORY

      IMPLICIT NONE
!     Parameter List
      INTEGER N,ITEM_TYPE
      INTEGER*4 ADDR_PTR
      LOGICAL INITIALISE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR_SIZE
      PARAMETER(CERROR_SIZE=50)
      INTEGER CERROR(CERROR_SIZE),CERRLEN,ERR,INIT


      CALL ENTERS('GET_MEMORY',*9999)

      CALL ASSERT(N.GE.0,'>>N < 0, Integer overflow',ERROR,*9999)
      ERR=0
      IF(INITIALISE) THEN
        INIT=1
      ELSE
        INIT=0
      ENDIF
      CALL MALLOCMEMORY(N,ITEM_TYPE,INIT,ADDR_PTR,
     '  ERR,CERROR,CERROR_SIZE)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        GOTO 9999
      ENDIF

      CALL EXITS('GET_MEMORY')
      RETURN
 9999 CALL ERRORS('GET_MEMORY',ERROR)
      CALL EXITS('GET_MEMORY')
      RETURN 1
      END


      SUBROUTINE BIG_GET_MEMORY(N1,N2,ITEM_TYPE,ADDR_PTR,INITIALISE,
     '  ERROR,*)

C#### Subroutine: BIG_GET_MEMORY
C###  Description:
C###    BIG_GET_MEMORY dynamically allocates N1*N2 items of type ITEM_TYPE
C###    in memory. ITEM_TYPES are 1/2/3/4/5/6/7/8 for INTEGER,REAL,
C###    REAL*8,CHARACTER,LOGICAL,INTEGER*2,COMPLEX and COMPLEX*16.
C###    Constants are setup in mach00.cmn for these types. If
C###    INITIALISE is true the variables are initialised to a
C###    special value (currently only for INTEGER's, REAL's and
C###    REAL*8's). The only difference between BIG_GET_MEMORY and 
C###    GET_MEMORY is that the amount of requested memory is passed as two 
C###    integer variables to overcome possible 32 bit overflow problems when
C###    large amounts of memory are requested in 64 bit executables.
C###  See-Also: GET_MEMORY,ALLOCATE_MEMORY,FREE_MEMORY

      IMPLICIT NONE
!     Parameter List
      INTEGER N1,N2,ITEM_TYPE
      INTEGER*4 ADDR_PTR
      LOGICAL INITIALISE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR_SIZE
      PARAMETER(CERROR_SIZE=50)
      INTEGER CERROR(CERROR_SIZE),CERRLEN,ERR,INIT


      CALL ENTERS('BIG_GET_MEMORY',*9999)

      CALL ASSERT(N1.GE.0,'>>N1 < 0, Integer overflow',ERROR,*9999)
      CALL ASSERT(N2.GE.0,'>>N2 < 0, Integer overflow',ERROR,*9999)
      ERR=0
      IF(INITIALISE) THEN
        INIT=1
      ELSE
        INIT=0
      ENDIF
      CALL BIGMALLOCMEMORY(N1,N2,ITEM_TYPE,INIT,ADDR_PTR,
     '  ERR,CERROR,CERROR_SIZE)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        GOTO 9999
      ENDIF

      CALL EXITS('BIG_GET_MEMORY')
      RETURN
 9999 CALL ERRORS('BIG_GET_MEMORY',ERROR)
      CALL EXITS('BIG_GET_MEMORY')
      RETURN 1
      END


