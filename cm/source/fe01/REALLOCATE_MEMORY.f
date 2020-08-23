      SUBROUTINE REALLOCATE_MEMORY(N,ITEM_TYPE,ADDR_PTR,ERROR,*)

C#### Subroutine: REALLOCATE_MEMORY
C###  Description:
C###    REALLOCATE_MEMORY dynamically allocates the maximum of (N,MINIMUM)
C###    items of type ITEM_TYPE in memory. If the address pointer is
C###    not zero REALLOCATE_MEMORY will free the memory first before
C###    allocating. ITEM_TYPES are 1/2/3/4/5/6/7/8/9 for INTEGER, REAL,
C###    REAL*8, CHARACTER, LOGICAL, INTEGER*2, COMPLEX, COMPLEX*16, and
C###    POINTER.
C###    Constants are setup in mach00.cmn for these types. If
C###    INITIALISE is true the variables are initialised to a special
C###    value (currently only for INTEGER's, REAL's and REAL*8's).
C###  See-Also: GET_MEMORY,FREE_MEMORY

      IMPLICIT NONE
!     Parameter List
      INTEGER N,ITEM_TYPE
      INTEGER*4 ADDR_PTR
!      LOGICAL
      CHARACTER ERROR*(*)
      INTEGER CERROR_SIZE
      PARAMETER(CERROR_SIZE=50)
      INTEGER CERROR(CERROR_SIZE),CERRLEN,ERR
!     Local Variables

      CALL ENTERS('REALLOCATE_MEMORY',*9999)

      CALL ASSERT(N.GE.0,'>>N < 0, Integer overflow',ERROR,*9999)
      ERR=0
      CALL REALLOCMEMORY(N,ITEM_TYPE,ADDR_PTR,
     '  ERR,CERROR,CERROR_SIZE)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        GOTO 9999
      ENDIF

      CALL EXITS('REALLOCATE_MEMORY')
      RETURN
 9999 CALL ERRORS('REALLOCATE_MEMORY',ERROR)
      CALL EXITS('REALLOCATE_MEMORY')
      RETURN 1
      END


