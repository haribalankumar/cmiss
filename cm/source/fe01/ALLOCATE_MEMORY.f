      SUBROUTINE ALLOCATE_MEMORY(N,MINIMUM,ITEM_TYPE,ADDR_PTR,
     '  INITIALISE,ERROR,*)

C#### Subroutine: ALLOCATE_MEMORY
C###  Description:
C###    ALLOCATE_MEMORY dynamically allocates the maximum of (N,MINIMUM)
C###    items of type ITEM_TYPE in memory. If the address pointer is
C###    not zero ALLOCATE_MEMORY will free the memory first before
C###    allocating. ITEM_TYPES are 1/2/3/4/5/6/7/8 for INTEGER, REAL,
C###    REAL*8, CHARACTER, LOGICAL, INTEGER*2, COMPLEX and COMPLEX*16.
C###    Constants are setup in mach00.cmn for these types. If
C###    INITIALISE is true the variables are initialised to a special
C###    value (currently only for INTEGER's, REAL's and REAL*8's).
C###  See-Also: GET_MEMORY,FREE_MEMORY

      IMPLICIT NONE
!     Parameter List
      INTEGER N,MINIMUM,ITEM_TYPE
      INTEGER*4 ADDR_PTR
      LOGICAL INITIALISE
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('ALLOCATE_MEMORY',*9999)

      IF(ADDR_PTR.NE.0) CALL FREE_MEMORY(ADDR_PTR,ERROR,*9999)

      CALL GET_MEMORY(MAX(N,MINIMUM),ITEM_TYPE,ADDR_PTR,INITIALISE,
     '  ERROR,*9999)

      CALL EXITS('ALLOCATE_MEMORY')
      RETURN
 9999 CALL ERRORS('ALLOCATE_MEMORY',ERROR)
      CALL EXITS('ALLOCATE_MEMORY')
      RETURN 1
      END

      SUBROUTINE BIG_ALLOCATE_MEMORY(N1,N2,MINIMUM,ITEM_TYPE,ADDR_PTR,
     '  INITIALISE,ERROR,*)

C#### Subroutine: BIG_ALLOCATE_MEMORY
C###  Description:
C###    BIG_ALLOCATE_MEMORY dynamically allocates the maximum of (N1*N2,MINIMUM)
C###    items of type ITEM_TYPE in memory. If the address pointer is
C###    not zero BIG_ALLOCATE_MEMORY will free the memory first before
C###    allocating. ITEM_TYPES are 1/2/3/4/5/6/7/8 for INTEGER, REAL,
C###    REAL*8, CHARACTER, LOGICAL, INTEGER*2, COMPLEX and COMPLEX*16.
C###    Constants are setup in mach00.cmn for these types. If
C###    INITIALISE is true the variables are initialised to a special
C###    value (currently only for INTEGER's, REAL's and REAL*8's).
C###    The only difference between BIG_ALLOCATE_MEMORY and ALLOCATE_MEMORY
C###    is that the amount of requested memory is passed as two integer 
C###    variables to overcome possible 32 bit overflow problems when
C###    large amounts of memory are reqested in 64 bit executables.
C###  See-Also: ALLOCATE_MEMORY,GET_MEMORY,FREE_MEMORY

      IMPLICIT NONE
!     Parameter List
      INTEGER N1,N2,MINIMUM,ITEM_TYPE
      INTEGER*4 ADDR_PTR
      LOGICAL INITIALISE
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('BIG_ALLOCATE_MEMORY',*9999)

      IF(ADDR_PTR.NE.0) CALL FREE_MEMORY(ADDR_PTR,ERROR,*9999)

      CALL BIG_GET_MEMORY(MAX(N1,MINIMUM),MAX(N2,MINIMUM),ITEM_TYPE,
     '  ADDR_PTR,INITIALISE,
     '  ERROR,*9999)

      CALL EXITS('BIG_ALLOCATE_MEMORY')
      RETURN
 9999 CALL ERRORS('BIG_ALLOCATE_MEMORY',ERROR)
      CALL EXITS('BIG_ALLOCATE_MEMORY')
      RETURN 1
      END


