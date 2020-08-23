      SUBROUTINE SETEXAMPLEDIR(EXAMPLENAME,COMFILENAME,ERROR,*)

C#### Subroutine: SETEXAMPLEDIR
C###  Description:
C###    SETEXAMPLEDIR sets the example directory and subdirectory based
C###    on the given example number.  COMFILENAME is set to the name of
C###    the default comfile.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='SETEXAMPLEDIR')
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'disp00.cmn'
!     Parameter List
      CHARACTER EXAMPLENAME*(*),COMFILENAME*(*),ERROR*(*)
!     Local Variables
      INTEGER ERR,LENGTH
      INTEGER*4 COMNAME_PTR,EXAMPLE_DIR_PTR,EXAMPLENAME_PTR
!     Functions
      INTEGER C_STRLEN,LEN_TRIM1
      INTEGER*4 RESOLVE_EXAMPLE_PATH

      COMNAME_PTR=0
      EXAMPLE_DIR_PTR=0
      EXAMPLENAME_PTR=0
      CALL ENTERS(ROUTINENAME,*9999)

C Find the example subdirectory from the example number
C KAT 20/4/00: Using script in example tree to get the path if
C     EXAMPLENAME is not blank.
      IF(EXAMPLENAME.NE.' ') THEN

        CALL CREATE_CSTRING(EXAMPLENAME_PTR,EXAMPLENAME)
        IF(EXAMPLENAME_PTR.EQ.0) GOTO 9998

        EXAMPLE_DIR_PTR=RESOLVE_EXAMPLE_PATH(%VAL(EXAMPLE_ROOT_PTR),
     '    %VAL(EXAMPLENAME_PTR),COMNAME_PTR,%VAL(0))
        IF(EXAMPLE_DIR_PTR.EQ.0) GOTO 9998

        CALL DESTROY_CSTRING(EXAMPLENAME_PTR)

        LENGTH = C_STRLEN(%VAL(EXAMPLE_DIR_PTR))
        CALL ASSERT(LENGTH.LE.LEN(EXAMPLES_DIR),
     '    'Example directory is too long',ERROR,*9999)

        CALL C2FSTRING(%VAL(EXAMPLE_DIR_PTR),LENGTH,EXAMPLES_DIR)

C       Set COMFILENAME to the default comfile name.
        IF(COMNAME_PTR.EQ.0) THEN
C         No default name available
          COMFILENAME=' '
        ELSE
          LENGTH = C_STRLEN(%VAL(COMNAME_PTR))
          IF(LENGTH.GT.LEN(COMFILENAME)) THEN
            CALL FLAG_ERROR(0,'Default comfile name is too long')
          ENDIF
          CALL C2FSTRING(%VAL(COMNAME_PTR),LENGTH,COMFILENAME)
        ENDIF

        CALL DESTROY_EXAMPLE_PATH(EXAMPLE_DIR_PTR,COMNAME_PTR)

      ELSE !Just use the example path

        LENGTH = C_STRLEN(%VAL(EXAMPLE_ROOT_PTR))
        CALL ASSERT(LENGTH.LT.LEN(EXAMPLES_DIR),
     '    'Example directory is too long',ERROR,*9999)
        CALL C2FSTRING(%VAL(EXAMPLE_ROOT_PTR),LENGTH,EXAMPLES_DIR)
        LENGTH=LENGTH+1
        EXAMPLES_DIR(LENGTH:LENGTH)='/'
        COMFILENAME=' '
      ENDIF

C     CALL STRING_TRIM(EXAMPLENUM,IBEG1,IEND1)
C      SUBDIR=EXAMPLENUM(IBEG1:IBEG1)
C      DO i=IBEG1+1,IEND1
C        CALL STRING_TRIM(SUBDIR,IBEG2,IEND2)
C        SUBDIR=SUBDIR(1:IEND2)//'/'
C     '    //EXAMPLENUM(IBEG1:IBEG1+i-1)
C      ENDDO !subdir='3/31'
C      CALL STRING_TRIM(SUBDIR,IBEG2,IEND2)

CC Locate the root location for the CMISS examples
C      CALL STRING_TRIM(EXAMPLE_PATH,IBEG3,IEND3)

CC Append the example subdirectory to the cmiss examples root directory
CC KAT 31/3/00: if subdirectory is not blank
C      EXAMPLES_DIR=EXAMPLE_PATH(IBEG3:IEND3)//'/'
C      IF(SUBDIR(IBEG2:IEND2).NE.' ') THEN
C        CALL STRING_TRIM(EXAMPLES_DIR,IBEG4,IEND4)
C        EXAMPLES_DIR(IEND4+1:)=SUBDIR(IBEG2:IEND2)//'/'
C      ENDIF
CC      EXAMPLES_DIR=EXAMPLE_PATH(IBEG3:IEND3)//'/'
CC     '  //SUBDIR(IBEG2:IEND2)//'/'
C      CALL STRING_TRIM(EXAMPLES_DIR,IBEG4,IEND4)

C Set the interpreter's user variable.
      IF(.NOT.CMGUI_LINK) THEN
        LENGTH=LEN_TRIM1(EXAMPLES_DIR)
        CALL SET_USER_CHARACTER('cmiss::example',
     '    EXAMPLES_DIR(:LENGTH),ERR)
        IF(ERR.NE.0) THEN
          CALL FLAG_ERROR(-1,'Unable to set user var "example"')
          GOTO 9998
        ENDIF
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' EXAMPLENAME='',A,'
     '    //'/,'' EXAMPLES_DIRECTORY='',A)')
     '    EXAMPLENAME, EXAMPLES_DIR(:LENGTH)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9998 ERROR=' '
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      IF(EXAMPLENAME_PTR.NE.0) THEN
        CALL DESTROY_CSTRING(EXAMPLENAME_PTR)
      ENDIF
      IF(EXAMPLE_DIR_PTR.NE.0) THEN
        CALL DESTROY_EXAMPLE_PATH(EXAMPLE_DIR_PTR,COMNAME_PTR)
      ENDIF
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


