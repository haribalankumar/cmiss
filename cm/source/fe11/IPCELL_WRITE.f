      SUBROUTINE IPCELL_WRITE(CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '  CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,nr,CELL_RCQS_VALUE,
     '  CELL_YQS_VALUE,CELL_ICQS_NAMES,
     '  CELL_RCQS_NAMES,CELL_YQS_NAMES,CQ,ERROR,*)

C#### Subroutine: IPCELL_WRITE
C###  Description:
C###    IPCELL_WRITE outputs all material and cell model parameters for
C###    user defined cellular models.

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'cell02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
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
      INTEGER i,IBEG,IEND,j
      CHARACTER VALUE*30

      CALL ENTERS('IPCELL_WRITE',*9999)

C     Write out the problem sizes
      WRITE(IFILE,'('' The number of cell model variants is: '',I5)')
     '  CELL_NUM_VARIANTS
      WRITE(IFILE,'('' The number of state variables is: '',I5)')
     '  CELL_NUM_STATE(1)
      WRITE(IFILE,'('' The number of ODE variables is: '',I5)')
     '  CELL_NUM_ODE(1)
      WRITE(IFILE,'('' The number of derived variables is: '',I5)')
     '  CELL_NUM_DERIVED(1)
      WRITE(IFILE,'('' The number of cell model parameters is: '',I5)')
     '  CELL_NUM_MODEL(1)
      WRITE(IFILE,
     '  '('' The number of cell control parameters is: '',I5)')
     '  CELL_NUM_CONTROL(1)
      WRITE(IFILE,
     '  '('' The number of cell material parameters is: '',I5)')
     '  CELL_NUM_PARAMETERS(1)
      WRITE(IFILE,
     '  '('' The number of cell protocol parameters is: '',I5)')
     '  CELL_NUM_PROTOCOL(1)
      WRITE(IFILE,
     '  '('' The number of additional integer input '//
     '  'parameters is: '',I5)') CELL_NUM_AII(1)
      WRITE(IFILE,
     '  '('' The number of additional integer output '//
     '  'parameters is: '',I5)') CELL_NUM_AIO(1)
      WRITE(IFILE,
     '  '('' The number of additional real input '//
     '  'parameters is: '',I5)') CELL_NUM_ARI(1)
      WRITE(IFILE,
     '  '('' The number of additional real output '//
     '  'parameters is: '',I5)') CELL_NUM_ARO(1)

      IF (CELL_NUM_STATE(1).GT.0) THEN
        WRITE(IFILE,*)
        WRITE(IFILE,'('' State variables:'')')
        DO i=1,CELL_NUM_VARIANTS
          WRITE(IFILE,'('' Cell variant '',I5,'':'')') i
          DO j=1,CELL_NUM_STATE(1)
            WRITE(VALUE,'(D12.4)')
     '        CELL_YQS_VALUE(CELL_STATE_OFFSET(1)+j-1,i)
            CALL STRING_TRIM(VALUE,IBEG,IEND)
            IF(CELL_YQS_SPATIAL(CELL_STATE_OFFSET(1)+j-1,i).EQ.0) THEN
              WRITE(IFILE,'(I5,'' '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_YQS_NAMES(CELL_STATE_OFFSET(1)+j-1,i)
            ELSE
              WRITE(IFILE,'(I5,''* '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_YQS_NAMES(CELL_STATE_OFFSET(1)+j-1,i)
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_STATE
      IF (CELL_NUM_MODEL(1).GT.0) THEN
        WRITE(IFILE,*)
        WRITE(IFILE,'('' Model variables:'')')
        DO i=1,CELL_NUM_VARIANTS
          WRITE(IFILE,'('' Cell variant '',I5,'':'')') i
          DO j=1,CELL_NUM_MODEL(1)
            WRITE(VALUE,'(I5)')
     '        CELL_ICQS_VALUE(CELL_MODEL_OFFSET(1)+j-1,i)
            CALL STRING_TRIM(VALUE,IBEG,IEND)
            IF(CELL_ICQS_SPATIAL(CELL_MODEL_OFFSET(1)+j-1,i).EQ.0) THEN
              WRITE(IFILE,'(I5,'' '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_ICQS_NAMES(CELL_MODEL_OFFSET(1)+j-1,i)
            ELSE
              WRITE(IFILE,'(I5,''* '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_ICQS_NAMES(CELL_MODEL_OFFSET(1)+j-1,i)
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_MODEL
      IF (CELL_NUM_CONTROL(1).GT.0) THEN
        WRITE(IFILE,*)
        WRITE(IFILE,'('' Control variables:'')')
        DO i=1,CELL_NUM_VARIANTS
          WRITE(IFILE,'('' Cell variant '',I5,'':'')') i
          DO j=1,CELL_NUM_CONTROL(1)
            WRITE(VALUE,'(I5)')
     '        CELL_ICQS_VALUE(CELL_CONTROL_OFFSET(1)+j-1,i)
            CALL STRING_TRIM(VALUE,IBEG,IEND)
            IF(CELL_ICQS_SPATIAL(CELL_CONTROL_OFFSET(1)+j-1,i).EQ.0)
     '        THEN
              WRITE(IFILE,'(I5,'' '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_ICQS_NAMES(CELL_CONTROL_OFFSET(1)+j-1,i)
            ELSE
              WRITE(IFILE,'(I5,''* '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_ICQS_NAMES(CELL_CONTROL_OFFSET(1)+j-1,i)
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_CONTROL
      IF (CELL_NUM_PARAMETERS(1).GT.0) THEN
        WRITE(IFILE,*)
        WRITE(IFILE,'('' Parameter variables:'')')
        DO i=1,CELL_NUM_VARIANTS
          WRITE(IFILE,'('' Cell variant '',I5,'':'')') i
          DO j=1,CELL_NUM_PARAMETERS(1)
            WRITE(VALUE,'(D12.4)')
     '        CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(1)+j-1,i)
            CALL STRING_TRIM(VALUE,IBEG,IEND)
            IF(CELL_RCQS_SPATIAL(CELL_PARAMETERS_OFFSET(1)+j-1,i).EQ.0)
     '        THEN
              WRITE(IFILE,'(I5,'' '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_RCQS_NAMES(CELL_PARAMETERS_OFFSET(1)+j-1,i)
            ELSE
              WRITE(IFILE,'(I5,''* '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_RCQS_NAMES(CELL_PARAMETERS_OFFSET(1)+j-1,i)
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_PARAMETERS
      IF (CELL_NUM_PROTOCOL(1).GT.0) THEN
        WRITE(IFILE,*)
        WRITE(IFILE,'('' Protocol variables:'')')
        DO i=1,CELL_NUM_VARIANTS
          WRITE(IFILE,'('' Cell variant '',I5,'':'')') i
          DO j=1,CELL_NUM_PROTOCOL(1)
            WRITE(VALUE,'(D12.4)')
     '        CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(1)+j-1,i)
            CALL STRING_TRIM(VALUE,IBEG,IEND)
            IF(CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(1)+j-1,i).EQ.0)
     '        THEN
              WRITE(IFILE,'(I5,'' '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(1)+j-1,i)
            ELSE
              WRITE(IFILE,'(I5,''* '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(1)+j-1,i)
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_PROTOCOL
      IF (CELL_NUM_AII(1).GT.0) THEN
        WRITE(IFILE,*)
        WRITE(IFILE,'('' Additional integer input variables:'')')
        DO i=1,CELL_NUM_VARIANTS
          WRITE(IFILE,'('' Cell variant '',I5,'':'')') i
          DO j=1,CELL_NUM_AII(1)
            WRITE(VALUE,'(I5)')
     '        CELL_ICQS_VALUE(CELL_AII_OFFSET(1)+j-1,i)
            CALL STRING_TRIM(VALUE,IBEG,IEND)
            IF(CELL_ICQS_SPATIAL(CELL_AII_OFFSET(1)+j-1,i).EQ.0)
     '        THEN
              WRITE(IFILE,'(I5,'' '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_ICQS_NAMES(CELL_AII_OFFSET(1)+j-1,i)
            ELSE
              WRITE(IFILE,'(I5,''* '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_ICQS_NAMES(CELL_AII_OFFSET(1)+j-1,i)
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_AII
      IF (CELL_NUM_ARI(1).GT.0) THEN
        WRITE(IFILE,*)
        WRITE(IFILE,'('' Additional real input variables:'')')
        DO i=1,CELL_NUM_VARIANTS
          WRITE(IFILE,'('' Cell variant '',I5,'':'')') i
          DO j=1,CELL_NUM_ARI(1)
            WRITE(VALUE,'(D12.4)')
     '        CELL_RCQS_VALUE(CELL_ARI_OFFSET(1)+j-1,i)
            CALL STRING_TRIM(VALUE,IBEG,IEND)
            IF(CELL_RCQS_SPATIAL(CELL_ARI_OFFSET(1)+j-1,i).EQ.0)
     '        THEN
              WRITE(IFILE,'(I5,'' '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_RCQS_NAMES(CELL_ARI_OFFSET(1)+j-1,i)
            ELSE
              WRITE(IFILE,'(I5,''* '//VALUE(IBEG:IEND)//' '',A)') j,
     '          CELL_RCQS_NAMES(CELL_ARI_OFFSET(1)+j-1,i)
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_ARI

      CALL EXITS('IPCELL_WRITE')
      RETURN
 9999 CALL ERRORS('IPCELL_WRITE',ERROR)
      CALL EXITS('IPCELL_WRITE')
      RETURN 1
      END


