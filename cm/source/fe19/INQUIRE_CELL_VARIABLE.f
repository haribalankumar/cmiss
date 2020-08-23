      SUBROUTINE INQUIRE_CELL_VARIABLE(CELL_ICQS_NAMES,
     '  CELL_RCQS_NAMES,CELL_YQS_NAMES,STRING,ERROR,*)

C#### Subroutine: INQUIRE_CELL_VARIABLE
C###  Description:
C###    <HTML>
C###    INQUIRE_CELL_VARIABLE searches for a specified CELLVARNAME
C###    in the character arrays that contain cell variable names.
C###    If desired the cell array name (CELLARRAYNAME = YQS, ICQS or
C###    RCQS) that contains the named variable and the index
C###    of the cell variable (CELLVARINDEXNAME) in the array can be
C###    returned into user defined variables (which are should not
C###    be defined before this routine is called).
C###    </HTML>
C*** Created by Martyn Nash July 2000

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      CHARACTER CELL_ICQS_NAMES(NQIM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH),
     '  ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER CELLVARINDEX,ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  N3CO,niq,nqv,NUMRETVARS
      LOGICAL CBBREV,RETVAR,FOUNDVAR
      CHARACTER CELLARRAYNAME*4,CELL_NAME_UC*(CELL_NAME_LENGTH),
     '  CELLVARNAME*(CELL_NAME_LENGTH),
     '  CELLVARNAME_UC*(CELL_NAME_LENGTH),
     '  RETVARNAMES(2)*(CELL_NAME_LENGTH)

      CALL ENTERS('INQUIRE_CELL_VARIABLE',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------
C#### Command: FEM inquire cell_variable CELLVARNAME
C###    Specifies the list of grid point variable indices.
C###  Parameter:      <return_variables CELLARRAYNAME,CELLVARINDEXNAME>
C###    Specify the value to return and the user variable name to use.
C###  Description:
C###    Searches for a specified CELLVARNAME
C###    in the character arrays that contain cell variable names.
C###    If desired the cell array name (CELLARRAYNAME = YQS, ICQS or
C###    RCQS) that contains the named variable and the index
C###    of the cell variable (CELLVARINDEXNAME) in the array can be
C###    returned into user defined variables.  It is assumed that
C###    these two user variables are not defined before this command
C###    is issued.

        OP_STRING(1)=STRING(IBEG:IEND)//' CELLVARNAME'
        OP_STRING(2)=BLANK(1:15)
     '    //'<return_variables CELLARRAYNAME,CELLVARINDEXNAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe19','doc','INQUIRE_CELL_VARIABLE',ERROR,*9999)
      ELSE
        CALL ASSERT(NTCO.GT.noco,'>>Must supply a cell '
     '    //'variable name',ERROR,*9999)
        CALL STRING_TRIM(CO(noco+1),IBEG1,IEND1)
        CELLVARNAME=CO(noco+1)(IBEG1:IEND1)

C       Check for return variables
        IF(CBBREV(CO,'RETURN_VARIABLES',3,noco+2,NTCO,N3CO)) THEN
          RETVAR=.TRUE.
C         Set the name of the variable to which the value is returned
          CALL PARSSL(CO(N3CO+1),2,NUMRETVARS,RETVARNAMES,ERROR,*9999)
          CALL ASSERT(NUMRETVARS.EQ.2,'>>Must supply return user '
     '      //'variable names for cell array name AND variable index',
     '      ERROR,*9999)
        ELSE
          RETVAR=.FALSE.
        ENDIF

C       Search for the array name and index of the specified
C       cell variable
        CALL CUPPER(CELLVARNAME,CELLVARNAME_UC)
        CALL STRING_TRIM(CELLVARNAME_UC,IBEG1,IEND1)
        FOUNDVAR=.FALSE.
        nqv=0
        DO WHILE(nqv.LT.NQVT.AND..NOT.FOUNDVAR)
          nqv=nqv+1
C         Search CELL_YQS_NAMES
          niq=0
          DO WHILE(niq.LT.NIQST.AND..NOT.FOUNDVAR)
            niq=niq+1
            CALL CUPPER(CELL_YQS_NAMES(niq,nqv),CELL_NAME_UC)
            CALL STRING_TRIM(CELL_NAME_UC,IBEG2,IEND2)
            IF(CELLVARNAME_UC(IBEG1:IEND1).EQ.
     '        CELL_NAME_UC(IBEG2:IEND2)) THEN
              CELLARRAYNAME='YQS'
              CELLVARINDEX=niq
              FOUNDVAR=.TRUE.
            ENDIF
          ENDDO !while niq.LE.NIQST
C         Search CELL_ICQS_NAMES
          niq=0
          DO WHILE(niq.LT.NQIT.AND..NOT.FOUNDVAR)
            niq=niq+1
            CALL CUPPER(CELL_ICQS_NAMES(niq,nqv),CELL_NAME_UC)
            CALL STRING_TRIM(CELL_NAME_UC,IBEG2,IEND2)
            IF(CELLVARNAME_UC(IBEG1:IEND1).EQ.
     '        CELL_NAME_UC(IBEG2:IEND2)) THEN
              CELLARRAYNAME='ICQS'
              CELLVARINDEX=niq
              FOUNDVAR=.TRUE.
            ENDIF
          ENDDO !while niq.LE.NQIT
C         Search CELL_ICQS_NAMES
          niq=0
          DO WHILE(niq.LT.NQRT.AND..NOT.FOUNDVAR)
            niq=niq+1
            CALL CUPPER(CELL_RCQS_NAMES(niq,nqv),CELL_NAME_UC)
            CALL STRING_TRIM(CELL_NAME_UC,IBEG2,IEND2)
            IF(CELLVARNAME_UC(IBEG1:IEND1).EQ.
     '        CELL_NAME_UC(IBEG2:IEND2)) THEN
              CELLARRAYNAME='RCQS'
              CELLVARINDEX=niq
              FOUNDVAR=.TRUE.
            ENDIF
          ENDDO !while niq.LE.NQRT
        ENDDO !while nqv.LE.NQVT
        CALL ASSERT(FOUNDVAR,'Cell variable "'//CELLVARNAME(IBEG1:IEND1)
     '    //'" was not found in YQS, ICQS or RCQS',ERROR,*9999)

        CALL STRING_TRIM(CELLARRAYNAME,IBEG1,IEND1)
        IF(RETVAR) THEN
C         Assign the cell array na,e to the given user variable
          CALL STRING_TRIM(RETVARNAMES(1),IBEG2,IEND2)
          CALL SET_USER_CHARACTER(RETVARNAMES(1)(IBEG2:IEND2),
     '      CELLARRAYNAME(IBEG1:IEND1),ERR)
          IF(ERR.NE.0) THEN
            ERROR='Unable to set user var "'//
     '        RETVARNAMES(1)(IBEG2:IEND2)//'"'
            GOTO 9999
          ENDIF
C         Assign the cell variable index to the given user variable
          CALL STRING_TRIM(RETVARNAMES(2),IBEG2,IEND2)
          CALL SET_USER_INTEGER(RETVARNAMES(2)(IBEG2:IEND2),
     '      CELLVARINDEX,ERR)
          IF(ERR.NE.0) THEN
            ERROR='Unable to set user var "'//
     '        RETVARNAMES(2)(IBEG2:IEND2)//'"'
            GOTO 9999
          ENDIF
        ELSE !not RETVAR
C         Output the cell variable index to screen
          CALL STRING_TRIM(CELLVARNAME,IBEG2,IEND2)
          WRITE(OP_STRING,
     '      '('' Cell variable "'//CELLVARNAME(IBEG2:IEND2)
     '      //'" is in position '',I3,'' of the cell array '
     '      //CELLARRAYNAME(IBEG1:IEND1)//''')') CELLVARINDEX
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF !RETVAR

      ENDIF

      CALL EXITS('INQUIRE_CELL_VARIABLE')
      RETURN
 9999 CALL ERRORS('INQUIRE_CELL_VARIABLE',ERROR)
      CALL EXITS('INQUIRE_CELL_VARIABLE')
      RETURN 1
      END


