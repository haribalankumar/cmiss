      SUBROUTINE IPCELL_CELLML(CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '  CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,nr,CELL_RCQS_VALUE,
     '  CELL_YQS_VALUE,CELL_ICQS_NAMES,
     '  CELL_RCQS_NAMES,CELL_YQS_NAMES,CQ,ERROR,*)

C#### Subroutine: IPCELL_CELLML
C###  Description:
C###    IPCELL_CELLML inputs all cell model parameters for
C###    cellular models defined in cellML.

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'cell02.cmn'
      INCLUDE 'file01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'

      INCLUDE 'cellml.cmn'

!     Parameter List
      INTEGER CELL_ICQS_VALUE(NQIM,NQVM),CELL_ICQS_SPATIAL(NQIM,NQVM),
     '  CELL_RCQS_SPATIAL(NQRM,NQVM),CELL_YQS_SPATIAL(NIQSM,NQVM),nr
      REAL*8 CELL_RCQS_VALUE(NQRM,NQVM),CELL_YQS_VALUE(NIQSM,NQVM),
     '  CQ(NMM,NQM)
      CHARACTER CELL_ICQS_NAMES(NQIM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH),
     '  ERROR*(*)
!     Local Variables
      INTEGER ERR,IBEG,ICHAR,IEND,INFO,NOQUES,variant,VIBEG,VIEND,
     '  ERROR_CODE,NSTATE,NODE,NDERIVED,NPARAMETERS,INDEX1,INDEX2,
     '  variable,NIBEG,NIEND,NFCN,NFCNS,NCONSTS,i,j,nq
      CHARACTER CHAR*100,VCHAR*15,CHAR11*11,NAME*(CELL_NAME_LENGTH)
      LOGICAL TERMINATE_CALLED,INITIALISED,SAVE_FILES,DEBUG
      REAL*8 INITIAL_VALUE

      CALL ENTERS('IPCELL_CELLML',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

C *** Initialise parameters
      NIQST = 1
      NQIT = 1
      NQRT = 1
      NQVT = 1

C *** Need to initialise the cell arrays
      DO i=1,NQVM
        DO j=1,NQIM
          CELL_ICQS_VALUE(j,i) = 0
          CELL_ICQS_SPATIAL(j,i) = 0
          CELL_ICQS_NAMES(j,i) = 'unused'
        ENDDO
        DO j=1,NQRM
          CELL_RCQS_VALUE(j,i) = 0.0d0
          CELL_RCQS_SPATIAL(j,i) = 0
          CELL_RCQS_NAMES(j,i) = 'unused'
        ENDDO
        DO j=1,NIQSM
          CELL_YQS_VALUE(j,i) = 0.0d0
          CELL_YQS_SPATIAL(j,i) = 0
          CELL_YQS_NAMES(j,i) = 'unused'
        ENDDO
      ENDDO

C *** Initialise the CellML processor
      TERMINATE_CALLED = .FALSE.
      INITIALISED = .FALSE.
      CALL CELLML_INITIALISE(ERROR_CODE)
      CALL ASSERT(ERROR_CODE.EQ.0,
     '  'Failed to initialise CellML processor',ERROR,*9999)
      INITIALISED = .TRUE.

C *** First prompt for the number of models to be used - i.e. the number
C     of variants in the distributed simulation, so each variant can
C     have its own RHS routine and we can hack up the integrator to call
C     the appropriate one ??
      IDEFLT(1) = 1
      WRITE(CHAR,'(I1)') IDEFLT(1)
      CALL STRING_TRIM(CHAR,IBEG,IEND)
      FORMAT = '($,'' Enter the number of model variants ['
     '  //CHAR(IBEG:IEND)//']: '',I2)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
      CELL_NUM_VARIANTS = IDATA(1)

      CALL ASSERT(CELL_NUM_VARIANTS.LE.CELLML_MAX_MODELS,
     '  'Increase CELLML_MAX_MODELS in cellml.cmn',ERROR,*9999)

      CALL ASSERT(CELL_NUM_VARIANTS.LE.NQVM,'>>Increase NQVM',ERROR,
     '  *9999)

C *** Variant loop....
      DO variant = 1,CELL_NUM_VARIANTS

        WRITE(VCHAR,'(I12)') variant
        CALL STRING_TRIM(VCHAR,VIBEG,VIEND)

C       Default array sizes
        CELL_NUM_STATE(variant) = 1
        CELL_NUM_ODE(variant) = 1
        CELL_NUM_DERIVED(variant) = 0
        CELL_NUM_MODEL(variant) = 0
        CELL_NUM_CONTROL(variant) = 1 ! ODE is reserved
        CELL_NUM_PARAMETERS(variant) = 2
        CELL_NUM_PROTOCOL(variant) = 10
        CELL_NUM_AII(variant) = 0
        CELL_NUM_AIO(variant) = 0
        CELL_NUM_ARI(variant) = 0
        CELL_NUM_ARO(variant) = 0

        CELLML_CONTAINS_VM(variant) = .FALSE.
        CELLML_CONTAINS_CM(variant) = .FALSE.
        CELLML_CONTAINS_AM(variant) = .FALSE.
        CELLML_CONTAINS_ISTIM(variant) = .FALSE.
        CELLML_ISTIM_ARRAY_NAME(variant) = ' '
        CELLML_ISTIM_ARRAY_INDEX(variant) = 0

C *** *** Get the URI of the file containing the cellML model
C       description
        FORMAT =
     '    '($,'' Enter the URI of the cellML file for model variant '
     '    //VCHAR(VIBEG:VIEND)//': '',A)'
        CDEFLT(1) = ' '
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        CALL STRING_TRIM(CDATA(1),IBEG,IEND)
        CELLML_URIS(variant) = CDATA(1)(IBEG:IEND)

C ***   Parse the cellML file - if we're gonna have multiple models in
C       memory at the same time, we're most likely gonna need some kind
C       of handle on each one...which is the variant number ??
        CALL CELLML_PARSE(variant,ERROR_CODE)
        CALL ASSERT(ERROR_CODE.EQ.0,'The file failed to parse',ERROR,
     '    *9999)

C ***   And create the math object for the model
        CALL CELLML_CREATE_MATH(variant,ERROR_CODE)
        CALL ASSERT(ERROR_CODE.EQ.0,'Error creating math object',ERROR,
     '    *9999)

C *** *** Set up the ODE variables - these need to be put
C       first into the YQS array in order for them to be safely
C       integrated.
        CALL CELLML_GET_NUM_ODE(variant,NODE)
        NSTATE = NODE
c        CALL ASSERT(NODE.GE.1,
c     '    '>>Must have at least one differential variable',ERROR,*9999)
c        CALL ASSERT(NODE.LE.NIQSM,'>>Increase NIQSM',ERROR,*9999)

        INDEX1 = 0
        INDEX2 = 2
        DO variable=1,NODE
          ! Get the variable name
          CALL CELLML_GET_ODE_VARIABLE_NAME(variant,variable,
     '      NAME)
          CALL STRING_TRIM(NAME,NIBEG,NIEND)
          ADEFLT(1) = 'N'
          FORMAT='(/'' ODE Variable '//NAME(NIBEG:NIEND)//':'''//
     '      '/$,''   Membrane potential [N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELLML_CONTAINS_VM(variant) = .TRUE.
            INDEX1 = 1
          ELSE
            INDEX1 = INDEX2
            INDEX2 = INDEX2 + 1
          ENDIF
          IF(INDEX1.LE.NIQSM) THEN
            CELL_YQS_NAMES(INDEX1,variant) = NAME(NIBEG:NIEND)
            CALL CELLML_GET_ODE_VAR_INIT_VALUE(variant,variable,
     '        INITIAL_VALUE)
            RDEFLT(1) = INITIAL_VALUE
            WRITE(CHAR11,'(E11.4)') RDEFLT(1)
            CALL STRING_TRIM(CHAR11,IBEG,IEND)
            FORMAT='($,''   Initial value ['//CHAR11(IBEG:IEND)//
     '        ']: '',E11.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            CELL_YQS_VALUE(INDEX1,variant) = RDATA(1)
            ADEFLT(1) = 'N'
            FORMAT = '($,''   Spatially varying [N]? '',A)'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF (ADATA(1).EQ.'Y') THEN
              CELL_YQS_SPATIAL(INDEX1,variant) = -1
            ELSE
              CELL_YQS_SPATIAL(INDEX1,variant) = 0
            ENDIF
          ENDIF ! < NIQSM
          CALL CELLML_SET_VARIABLE_ARRAY(variant,NAME(NIBEG:NIEND),'Y')
          CALL CELLML_SET_VARIABLE_ARRAY_INDEX(variant,
     '      NAME(NIBEG:NIEND),INDEX1)
        ENDDO !NUM_ODE

        IF(CELLML_CONTAINS_VM(variant)) THEN
          CELL_NUM_ODE(variant) = NODE
          CELL_NUM_STATE(variant) = NODE
        ELSE
          CELL_NUM_ODE(variant) = NODE+1
          CELL_NUM_STATE(variant) = NODE+1
          IF(Vm.LE.NIQSM) THEN
            CELL_YQS_NAMES(Vm,variant) = 'Vm'
            CELL_YQS_SPATIAL(Vm,variant) = 0
            CELL_YQS_VALUE(Vm,variant) = 200.0d0
          ENDIF ! < NIQSM
        ENDIF

        !Need to get this right to get the derived variables in
        !the right place
        NODE = CELL_NUM_STATE(variant)

C ***   Set up the function variables - if the user wants to
C       save them they need to go into the derived array, otherwise they
C       are local variables in the cellular RHS routine
        CALL CELLML_GET_NUM_FUNCTIONS(variant,NFCNS)

        NDERIVED = 0
        NFCN = 0
        DO variable=1,NFCNS
          ! Get the variable name
          CALL CELLML_GET_FUNCTION_VAR_NAME(variant,variable,
     '      NAME)
          CALL STRING_TRIM(NAME,NIBEG,NIEND)
          ADEFLT(1) = 'N'
          FORMAT = '(/$,'' Save function variable '//NAME(NIBEG:NIEND)//
     '      ' [N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            NDERIVED = NDERIVED+1
            IF(NODE+NDERIVED.LE.NIQSM) THEN
              CELL_YQS_NAMES(NODE+NDERIVED,variant) = NAME(NIBEG:NIEND)
            ENDIF ! < NIQSM
            CALL CELLML_SET_VARIABLE_ARRAY(variant,NAME(NIBEG:NIEND),
     '        'DERIVED')
            CALL CELLML_SET_VARIABLE_ARRAY_INDEX(variant,
     '        NAME(NIBEG:NIEND),NDERIVED)
          ELSE
            NFCN = NFCN + 1
            CALL CELLML_SET_VARIABLE_ARRAY(variant,NAME(NIBEG:NIEND),
     '        'FCN')
            CALL CELLML_SET_VARIABLE_ARRAY_INDEX(variant,
     '        NAME(NIBEG:NIEND),NFCN)
          ENDIF
        ENDDO !NFCNS

        IF(NODE+NDERIVED.GT.NIQSM) THEN
          CALL FLAG_ERROR(0,' ')
          CALL WRITE_CHAR(IOER,'Increase NIQSM to at least ',ERR)
          CALL WRITE_INT(IOER,NODE+NDERIVED,ERR)
          CALL WRITE_CHAR(IOER,NEWLINE,ERR)
          GOTO 9998
        ENDIF

        CELL_NUM_DERIVED(variant) = NDERIVED

C ***   Set up the constant variables - if the user wants to save them,
C       they need to be in the DERIVED array; otherwise they go in the
C       PARAM array. - Oops, if you put a parameter in the DERIVED array
C       you cannot specify spatial variation...
        CALL CELLML_GET_NUM_CONSTANTS(variant,NCONSTS)

        INDEX1 = 0
        INDEX2 = 3
        NPARAMETERS = NCONSTS+2
        IF(NPARAMETERS.GT.NQRM) THEN
          CALL FLAG_ERROR(0,' ')
          CALL WRITE_CHAR(IOER,'Increase NQRM to at least ',ERR)
          CALL WRITE_INT(IOER,NPARAMETERS,ERR)
          CALL WRITE_CHAR(IOER,NEWLINE,ERR)
          GOTO 9998
        ENDIF
        DO variable=1,NCONSTS
          ! Get the variable name
          CALL CELLML_GET_CONST_VAR_NAME(variant,variable,
     '      NAME)
          CALL STRING_TRIM(NAME,NIBEG,NIEND)
          ADEFLT(1) = 'N'
          FORMAT='(/'' Constant Variable '//NAME(NIBEG:NIEND)//':'''//
     '      '/$,''   Specific membrane capacitance (Cm) [N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELLML_CONTAINS_CM(variant) = .TRUE.
            INDEX1 = 1
            NPARAMETERS = NPARAMETERS-1
          ELSE
            ADEFLT(1) = 'N'
            FORMAT=
     '        '($,''   Surface to volume ratio (Am) [N]? '',A)'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF (ADATA(1).EQ.'Y') THEN
              CELLML_CONTAINS_AM(variant) = .TRUE.
              INDEX1 = 2
              NPARAMETERS = NPARAMETERS-1
            ELSE
              INDEX1 = INDEX2
              INDEX2 = INDEX2 + 1
            ENDIF
          ENDIF
          CELL_RCQS_NAMES(INDEX1,variant) = NAME(NIBEG:NIEND)
          ADEFLT(1) = 'N'
          FORMAT='($,''   Total electrical stimulus current [N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELLML_CONTAINS_ISTIM(variant) = .TRUE.
            CELLML_ISTIM_ARRAY_NAME(variant) = 'PARAM'
            CELLML_ISTIM_ARRAY_INDEX(variant) = INDEX1
            CELL_RCQS_VALUE(INDEX1,variant) = 0.0d0
            CELL_RCQS_SPATIAL(INDEX1,variant) = 0
          ELSE
            CALL CELLML_GET_CONST_VAR_INIT_VALUE(variant,
     '        variable,INITIAL_VALUE)
            RDEFLT(1) = INITIAL_VALUE
            WRITE(CHAR11,'(E11.4)') RDEFLT(1)
            CALL STRING_TRIM(CHAR11,IBEG,IEND)
            FORMAT='($,''   Initial value ['//CHAR11(IBEG:IEND)//
     '        ']: '',E11.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            CELL_RCQS_VALUE(INDEX1,variant) = RDATA(1)
            ADEFLT(1) = 'N'
            FORMAT = '($,''   Spatially varying [N]? '',A)'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF (ADATA(1).EQ.'Y') THEN
              CELL_RCQS_SPATIAL(INDEX1,variant) = -1
            ELSE
              CELL_RCQS_SPATIAL(INDEX1,variant) = 0
            ENDIF
          ENDIF
          CALL CELLML_SET_VARIABLE_ARRAY(variant,NAME(NIBEG:NIEND),
     '      'PARAM')
          CALL CELLML_SET_VARIABLE_ARRAY_INDEX(variant,
     '      NAME(NIBEG:NIEND),INDEX1)
        ENDDO !NCONSTS

        CELL_NUM_PARAMETERS(variant) = NPARAMETERS

        CALL IPCELL_COMMON(variant,ERROR,*9999)

C ***   Set the default/standard parameter names
        CELL_ICQS_NAMES(CELL_CONTROL_OFFSET(variant)+ODE-1,variant) =
     '    'ODEswitch'
        CELL_ICQS_SPATIAL(CELL_CONTROL_OFFSET(variant)+ODE-1,variant) =
     '    0
        CELL_ICQS_VALUE(CELL_CONTROL_OFFSET(variant)+ODE-1,variant) = 0
        IF (.NOT.CELLML_CONTAINS_CM(variant)) THEN
          CELL_RCQS_NAMES(CELL_PARAMETERS_OFFSET(variant)+Cm-1,variant)
     '      = 'Cm'
          RDEFLT(1) = 0.01d0
          FORMAT=
     '      '(/$,'' Membrane capcitance value, Cm (uF/mm^2) '
     '      //'[0.01]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(variant)+Cm-1,variant)
     '      = RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT='($,'' Membrane capcitance spatially varying '
     '      //'[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PARAMETERS_OFFSET(variant)+Cm-1,
     '        variant) = -1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PARAMETERS_OFFSET(variant)+Cm-1,
     '        variant) = 0
          ENDIF
        ENDIF
        IF (.NOT.CELLML_CONTAINS_AM(variant)) THEN
          CELL_RCQS_NAMES(CELL_PARAMETERS_OFFSET(variant)+Am-1,variant)
     '      = 'Am'
          RDEFLT(1) = 200.0d0
          FORMAT=
     '      '(/$,'' Surface-to-volume ratio, Am (1/mm) '
     '      //'[200]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(variant)+Am-1,variant)
     '      = RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT='($,'' Surface-to-volume ratio spatially varying '
     '      //'[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PARAMETERS_OFFSET(variant)+Am-1,
     '        variant) = -1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PARAMETERS_OFFSET(variant)+Am-1,
     '        variant) = 0
          ENDIF
        ENDIF
        IF (CELLML_CONTAINS_ISTIM(variant)) THEN
          ! Need to define the stimulus protocol
          !Pseudo stimulus
          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+PseudoIs-1,
     '      variant) ='PseudoIs'
          CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+PseudoIs-1,
     '      variant) = 0
          CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+PseudoIs-1,
     '      variant) =0.0d0
          !Stimulus 1 - start
          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+Is1start-1,
     '      variant) ='Is1start'
          RDEFLT(1) = 0.0d0
          FORMAT='(/'' Define stimulus protocol:'''//
     '      '/$,''   Stimulus 1 start time (ms) [0.0]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+Is1start-1,
     '      variant) = RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT='($,''   Stimulus 1 start time spatially varying '
     '      //'[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is1start-1,
     '        variant) =-1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is1start-1,
     '        variant) =0
          ENDIF
          !Stimulus 1 - stop
          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+Is1stop-1,
     '      variant) ='Is1stop'
          RDEFLT(1) = 0.0d0
          FORMAT='($,''   Stimulus 1 stop time (ms) [0.0]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+Is1stop-1,
     '      variant) =RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT='($,''   Stimulus 1 stop time spatially varying '
     '      //'[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is1stop-1,
     '        variant) =-1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is1stop-1,
     '        variant) =0
          ENDIF
          !Stimulus 1 - magnitude
          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+Is1current-1,
     '      variant) ='Is1current'
          RDEFLT(1) = 0.0d0
          FORMAT=
     '      '($,''   Stimulus 1 current (uA/mm^3) [0.0]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+Is1current-1,
     '      variant) =RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT='($,''   Stimulus 1 current spatially varying '
     '      //'[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is1current
     '        -1,variant)=-1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is1current
     '        -1,variant)=0
          ENDIF
          !Stimulus 2 - start
          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+Is2start-1,
     '      variant) ='Is2start'
          RDEFLT(1) = 0.0d0
          FORMAT='($,''   Stimulus 2 start time (ms) [0.0]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+Is2start-1,
     '      variant) =RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT='($,''   Stimulus 2 start time spatially varying '
     '      //'[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is2start-1,
     '        variant) =-1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is2start-1,
     '        variant) =0
          ENDIF
          !Stimulus 2 - stop
          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+Is2stop-1,
     '      variant) ='Is2stop'
          RDEFLT(1) = 0.0d0
          FORMAT='($,''   Stimulus 2 stop time (ms) [0.0]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+Is2stop-1,
     '      variant) =RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT='($,''   Stimulus 2 stop time spatially varying '
     '      //'[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is2stop-1,
     '        variant) =-1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is2stop-1,
     '        variant) =0
          ENDIF
          !Stimulus 2 - magnitude
          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+Is2current-1,
     '      variant) ='Is2current'
          RDEFLT(1) = 0.0d0
          FORMAT=
     '      '($,''   Stimulus 2 current (uA/mm^3) [0.0]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+Is2current-1,
     '      variant) =RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT='($,''   Stimulus 2 current spatially varying '
     '      //'[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is2current
     '        -1,variant)=-1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+Is2current
     '        -1,variant)=0
          ENDIF
          !Frequency stimulus - period
          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+IsFreqPeriod-1,
     '      variant) ='IsFreqPeriod'
          RDEFLT(1) = 0.0d0
          FORMAT=
     '      '($,''   Frequency stimulus period (ms) [0.0]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+IsFreqPeriod-1,
     '      variant) =RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT='($,''   Frequency stimulus period spatially varying '
     '      //'[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+IsFreqPeriod
     '        -1,variant) =-1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+IsFreqPeriod
     '        -1,variant) =0
          ENDIF
          !Frequency stimulus - duration
          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+IsFreqDuration
     '      -1,variant)='IsFreqDuration'
          RDEFLT(1) = 0.0d0
          FORMAT=
     '      '($,''   Frequency stimulus duration (ms) [0.0]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+IsFreqDuration
     '      -1,variant)=RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT=
     '      '($,''   Frequency stimulus duration spatially varying '//
     '      '[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)
     '        +IsFreqDuration-1,variant) =-1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)
     '        +IsFreqDuration-1,variant) =0
          ENDIF
          !Frequency stimulus - magnitude
          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+IsFreqMag-1,
     '      variant) ='IsFreqMag'
          RDEFLT(1) = 0.0d0
          FORMAT=
     '      '($,''   Frequency stimulus magnitude (uA/mm^3) '//
     '      '[0.0]: '',E11.4)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+IsFreqMag-1,
     '      variant) =RDATA(1)
          ADEFLT(1) = 'N'
          FORMAT=
     '      '($,''   Frequency stimulus magnitude spatially varying '//
     '      '[N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF (ADATA(1).EQ.'Y') THEN
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+IsFreqMag-1,
     '        variant) =-1
          ELSE
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+IsFreqMag-1,
     '        variant) =0
          ENDIF
        ELSE
          DO i=1,CELL_NUM_PROTOCOL(variant)
            CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(variant)+i-1,variant) =
     '        'Unused'
            CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(variant)+i-1,variant)
     '        = 0
            CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(variant)+i-1,variant) =
     '        0.0d0
          ENDDO
        ENDIF

        ADEFLT(1) = 'N'
        FORMAT='(/$,'' Debug code [N]? '',A)'
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF (ADATA(1).EQ.'Y') THEN
          DEBUG = .TRUE.
        ELSE
          DEBUG = .FALSE.
        ENDIF

        ADEFLT(1) = 'N'
        FORMAT='(/$,'' Save temporary files [N]? '',A)'
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF (ADATA(1).EQ.'Y') THEN
          SAVE_FILES = .TRUE.
        ELSE
          SAVE_FILES = .FALSE.
        ENDIF

C ***   Now have enough information to create the RHS routine
        CALL CELLML_CREATE_RHS_ROUTINE(variant,DEBUG,SAVE_FILES,
     '    ERROR_CODE)
        CALL ASSERT(ERROR_CODE.EQ.0,
     '    'Failed to create RHS routine',ERROR,*9999)

      ENDDO
C *** ... end variant loop

C *** Distributed modelling needs Cm and Am in the CQ array
C ??? Just use variant 1 ???
      DO nq=NQR(1,nr),NQR(2,nr)
        CQ(1,nq)=CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(1)+Cm-1,1)
        CQ(2,nq)=CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(1)+Am-1,1)
      ENDDO

C *** Terminate the CellML processor
      CALL CELLML_TERMINATE(ERROR_CODE)
      TERMINATE_CALLED=.TRUE.
      CALL ASSERT(ERROR_CODE.EQ.0,
     '  'Failed to terminate CellML processor',ERROR,*9999)

      CALL EXITS('IPCELL_CELLML')
      RETURN

 9998 ERROR=' '
 9999 CALL ERRORS('IPCELL_CELLML',ERROR)

C *** Terminate the CellML processor if it was initialised and not
C     terminated.
      IF (INITIALISED.AND..NOT.TERMINATE_CALLED) THEN
        CALL CELLML_TERMINATE(ERROR_CODE)
        TERMINATE_CALLED=.TRUE.
        CALL ASSERT(ERROR_CODE.EQ.0,
     '    'Failed to terminate CellML processor',ERROR,*9999)
      ENDIF

      CALL EXITS('IPCELL_CELLML')
      RETURN 1
      END


