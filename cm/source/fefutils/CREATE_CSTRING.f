      SUBROUTINE CREATE_CSTRING(C_STRING_PTR,F_STRING)

C#### Subroutine: CREATE_CSTRING
C###  Description:
C###    CREATE_CSTRING allocates memory to create a C
C###    string from a fortran string and sets a pointer to the location.
C###    If the memory allocation fails the pointer is set to NULL.
C###  See-Also: DESTROY_CSTRING
      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='CREATE_CSTRING')
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER*4 C_STRING_PTR
      CHARACTER F_STRING*(*)
!     Local Variables
      CHARACTER ERROR*100

C      CALL ENTER(ROUTINENAME)
      CALL ENTERS(ROUTINENAME,*1)
 1    CONTINUE

      CALL GET_MEMORY(LEN(F_STRING)+1,CHARTYPE,C_STRING_PTR,.FALSE.,
     '  ERROR,*9999)
      CALL F2CSTRING(%VAL(C_STRING_PTR),F_STRING)

C      CALL EXIT(ROUTINENAME)
      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORIN(ROUTINENAME)
      CALL EXITS(ROUTINENAME)
      RETURN
      END


