      SUBROUTINE OPCELL(CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '  CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,CELL_RCQS_VALUE,
     '  CELL_YQS_VALUE,CELL_ICQS_NAMES,CELL_RCQS_NAMES,CELL_YQS_NAMES,
     '  ERROR,*)

C#### Subroutine: OPCELL
C###  Description:
C###    OPCELL outputs a summary of cellular parameter numbers.

      IMPLICIT NONE

      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER CELL_ICQS_VALUE(NQIM,NQVM),CELL_ICQS_SPATIAL(NQIM,NQVM),
     '  CELL_RCQS_SPATIAL(NQRM,NQVM),CELL_YQS_SPATIAL(NIQSM,NQVM)
      REAL*8 CELL_RCQS_VALUE(NQRM,NQVM),CELL_YQS_VALUE(NIQSM,NQVM)
      CHARACTER CELL_ICQS_NAMES(NQIM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH),ERROR*(*)
!     Local Variables
      INTEGER i,j
      CHARACTER SPATIAL*1

      CALL ENTERS('OPCELL',*9999)

C *** write out the model size information
      WRITE(OP_STRING,
     '  '('' The number of cell model variants is: '',I12)')
     '  CELL_NUM_VARIANTS
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' The number of state variables is: '',I12)')
     '  CELL_NUM_STATE(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' The number of ODE variables is: '',I12)')
     '  CELL_NUM_ODE(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' The number of derived variables is: '',I12)')
     '  CELL_NUM_DERIVED(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '('' The number of cell model parameters is: '',I12)')
     '  CELL_NUM_MODEL(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' The number of cell control parameters '
     '  //'is: '',I12)') CELL_NUM_CONTROL(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' The number of cell material parameters '
     '  //'is: '',I12)') CELL_NUM_PARAMETERS(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' The number of cell protocol parameters '
     '  //'is: '',I12)') CELL_NUM_PROTOCOL(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' The number of additional integer input '
     '  //'parameters is: '',I12)') CELL_NUM_AII(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' The number of additional integer output '
     '  //'parameters is: '',I12)') CELL_NUM_AIO(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '('' The number of additional real input parameters '
     '  //'is: '',I12)') CELL_NUM_ARI(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '('' The number of additional real output parameters '
     '  //'is: '',I12)') CELL_NUM_ARO(1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C *** Write out the parameter values and names
      DO i=1,CELL_NUM_VARIANTS
        WRITE(OP_STRING,'('' Variant number '',I3)') i
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''   State Variables:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO j=1,CELL_NUM_STATE(i)
          IF(CELL_YQS_SPATIAL(CELL_STATE_OFFSET(i)-1+j,i).NE.0) THEN
            WRITE(SPATIAL,'(''*'')')
          ELSE
            WRITE(SPATIAL,'('' '')')
          ENDIF
          WRITE(OP_STRING,'(''     '',A,''= '',D15.8,A)')
     '      CELL_YQS_NAMES(CELL_STATE_OFFSET(i)-1+j,i),
     '      CELL_YQS_VALUE(CELL_STATE_OFFSET(i)-1+j,i),SPATIAL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(''   Model Parameters:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO j=1,CELL_NUM_MODEL(i)
          IF(CELL_ICQS_SPATIAL(CELL_MODEL_OFFSET(i)-1+j,i).NE.0) THEN
            WRITE(SPATIAL,'(''*'')')
          ELSE
            WRITE(SPATIAL,'('' '')')
          ENDIF
          WRITE(OP_STRING,'(''     '',A,'' = '',I12,A)')
     '      CELL_ICQS_NAMES(CELL_MODEL_OFFSET(i)-1+j,i),
     '      CELL_ICQS_VALUE(CELL_MODEL_OFFSET(i)-1+j,i),SPATIAL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(''   Control Parameters:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO j=1,CELL_NUM_CONTROL(i)
          IF(CELL_ICQS_SPATIAL(CELL_CONTROL_OFFSET(i)-1+j,i).NE.0) THEN
            WRITE(SPATIAL,'(''*'')')
          ELSE
            WRITE(SPATIAL,'('' '')')
          ENDIF
          WRITE(OP_STRING,'(''     '',A,'' = '',I12,A)')
     '      CELL_ICQS_NAMES(CELL_CONTROL_OFFSET(i)-1+j,i),
     '      CELL_ICQS_VALUE(CELL_CONTROL_OFFSET(i)-1+j,i),SPATIAL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(''   Material Parameters:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO j=1,CELL_NUM_PARAMETERS(i)
          IF(CELL_RCQS_SPATIAL(CELL_PARAMETERS_OFFSET(i)-1+j,i).NE.0)
     '      THEN
            WRITE(SPATIAL,'(''*'')')
          ELSE
            WRITE(SPATIAL,'('' '')')
          ENDIF
          WRITE(OP_STRING,'(''     '',A,'' = '',D15.8,A)')
     '      CELL_RCQS_NAMES(CELL_PARAMETERS_OFFSET(i)-1+j,i),
     '      CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(i)-1+j,i),SPATIAL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(''   Protocol Parameters:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO j=1,CELL_NUM_PROTOCOL(i)
          IF(CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(i)-1+j,i).NE.0) THEN
            WRITE(SPATIAL,'(''*'')')
          ELSE
            WRITE(SPATIAL,'('' '')')
          ENDIF
          WRITE(OP_STRING,'(''     '',A,'' = '',D15.8,A)')
     '      CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(i)-1+j,i),
     '      CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(i)-1+j,i),SPATIAL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(''   Additional Integer Input Parameters:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO j=1,CELL_NUM_AII(i)
          IF(CELL_ICQS_SPATIAL(CELL_AII_OFFSET(i)-1+j,i).NE.0) THEN
            WRITE(SPATIAL,'(''*'')')
          ELSE
            WRITE(SPATIAL,'('' '')')
          ENDIF
          WRITE(OP_STRING,'(''     '',A,'' = '',I12,A)')
     '      CELL_ICQS_NAMES(CELL_AII_OFFSET(i)-1+j,i),
     '      CELL_ICQS_VALUE(CELL_AII_OFFSET(i)-1+j,i),SPATIAL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(''   Additional Real Input Parameters:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO j=1,CELL_NUM_ARI(i)
          IF(CELL_RCQS_SPATIAL(CELL_ARI_OFFSET(i)-1+j,i).NE.0) THEN
            WRITE(SPATIAL,'(''*'')')
          ELSE
            WRITE(SPATIAL,'('' '')')
          ENDIF
          WRITE(OP_STRING,'(''     '',A,'' = '',D15.8,A)')
     '      CELL_RCQS_NAMES(CELL_ARI_OFFSET(i)-1+j,i),
     '      CELL_RCQS_VALUE(CELL_ARI_OFFSET(i)-1+j,i),SPATIAL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO !i=1,CELL_NUM_VARIANTS


      CALL EXITS('OPCELL')
      RETURN
 9999 CALL ERRORS('OPCELL',ERROR)
      CALL EXITS('OPCELL')
      RETURN 1
      END


