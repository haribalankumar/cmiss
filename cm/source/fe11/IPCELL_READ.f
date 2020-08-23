      SUBROUTINE IPCELL_READ(CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '  CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,nr,CELL_RCQS_VALUE,
     '  CELL_YQS_VALUE,CELL_ICQS_NAMES,
     '  CELL_RCQS_NAMES,CELL_YQS_NAMES,CQ,ERROR,*)

C#### Subroutine: IPCELL_READ
C###  Description:
C###    IPCELL inputs all material and cell model parameters for
C###    user defined cellular models. It also reads in each parameter's
C###    name from the ipcell file. Define materials must also be
C###    called for continuum modelling problems.

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
      INTEGER CLOCAT,CURRENT_OFFSET,i,IBEG,IBEG1,IEND,IEND1,IFROMC,j,
     '  IVALUE,nq,POS1,POS2,POS3,VARNUM,variant
      REAL*8 RVALUE
      CHARACTER VALUE*30,VARIABLE*30,LINE*132
      LOGICAL IVALID,RVALID

      CALL ENTERS('IPCELL_READ',*9999)

C     Read in the problem sizes
C *** DPN 24 April 2001 - this is bad, Fortran spec says you cannot read
C     in like this...
C      READ(IFILE,'('' The number of cell model variants is: '',I12)')
C     '  CELL_NUM_VARIANTS
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_VARIANTS = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of state variables is: '',I12)')
C     '  CELL_NUM_STATE
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_STATE(1) = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of ODE variables is: '',I12)')
C     '  CELL_NUM_ODE
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_ODE(1) = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of derived variables is: '',I12)')
C     '  CELL_NUM_DERIVED
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_DERIVED(1) = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of cell model parameters is: '',I12)')
C     '  CELL_NUM_MODEL
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_MODEL(1) = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of cell control parameters '
C     '  //'is: '',I12)') CELL_NUM_CONTROL
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_CONTROL(1) = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of cell material parameters '
C     '  //'is: '',I12)') CELL_NUM_PARAMETERS
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_PARAMETERS(1) = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of cell protocol parameters '
C     '  //'is: '',I12)') CELL_NUM_PROTOCOL
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_PROTOCOL(1) = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of additional integer input '
C     '  //'parameters is: '',I12)') CELL_NUM_AII
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_AII(1) = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of additional integer output '
C     '  //'parameters is: '',I12)') CELL_NUM_AIO
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_AIO(1) = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of additional real input parameters '
C     '  //'is: '',I12)') CELL_NUM_ARI
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
      CELL_NUM_ARI(1) = IFROMC(VARIABLE(IBEG1:IEND1))
C      READ(IFILE,'('' The number of additional real output parameters '
C     '  //'is: '',I12)') CELL_NUM_ARO
      READ(IFILE,'(A)') LINE
      CALL STRING_TRIM(LINE,IBEG,IEND)
      POS1=CLOCAT(':',LINE)
      VARIABLE=LINE(POS1+1:IEND)
      CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)

      CELL_NUM_ARO(1) = IFROMC(VARIABLE(IBEG1:IEND1))

C *** Initialise parameters
      NIQST = 1
      NQIT = 1
      NQRT = 1
      NQVT = 1

      CALL IPCELL_COMMON(1,ERROR,*9999)
      ! Set all array sizes equal
      DO variant=2,CELL_NUM_VARIANTS
        CELL_NUM_STATE(variant) = CELL_NUM_STATE(1)
        CELL_NUM_ODE(variant) = CELL_NUM_ODE(1)
        CELL_NUM_DERIVED(variant) = CELL_NUM_DERIVED(1)
        CELL_NUM_MODEL(variant) = CELL_NUM_MODEL(1)
        CELL_NUM_CONTROL(variant) = CELL_NUM_CONTROL(1)
        CELL_NUM_PARAMETERS(variant) = CELL_NUM_PARAMETERS(1)
        CELL_NUM_PROTOCOL(variant) = CELL_NUM_PROTOCOL(1)
        CELL_NUM_AII(variant) = CELL_NUM_AII(1)
        CELL_NUM_AIO(variant) = CELL_NUM_AIO(1)
        CELL_NUM_ARI(variant) = CELL_NUM_ARI(1)
        CELL_NUM_ARO(variant) = CELL_NUM_ARO(1)
        CALL IPCELL_COMMON(variant,ERROR,*9999)
      ENDDO

      IF (CELL_NUM_STATE(1).GT.0) THEN
C       Skip blank line
        READ(IFILE,'(A)') LINE
C       Skip state variable data line
        READ(IFILE,'(A)') LINE
        DO i=1,CELL_NUM_VARIANTS
C         Skip cell variant number line
          READ(IFILE,'(A)') LINE
          DO j=1,CELL_NUM_STATE(1)
            READ(IFILE,'(A)') LINE
            CALL STRING_TRIM(LINE,IBEG,IEND)
            POS1=CLOCAT(' ',LINE(IBEG:IEND))
C ***       CHAR(9) = TAB
            IF(POS1.EQ.0) POS1=CLOCAT(CHAR(9),LINE(IBEG:IEND))
            IF(POS1.NE.0) THEN
              VARIABLE=LINE(IBEG:IBEG+POS1-1)
              CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
              IF(VARIABLE(IEND1:IEND1).EQ.'*') THEN
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1-1))
                CELL_YQS_SPATIAL(CELL_STATE_OFFSET(1)+VARNUM-1,i)=-1
              ELSE
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1))
                CELL_YQS_SPATIAL(CELL_STATE_OFFSET(1)+VARNUM-1,i)=0
              ENDIF
              POS2=CLOCAT(' ',LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) POS2=CLOCAT(CHAR(9),LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) THEN
                VALUE=LINE(IBEG+POS1:IEND)
              ELSE
                VALUE=LINE(IBEG+POS1:IBEG+POS1+POS2-1)
              ENDIF
              CALL STRING_TRIM(VALUE,IBEG1,IEND1)
              IF(RVALID(VALUE(IBEG1:IEND1))) THEN
                READ(VALUE(IBEG1:IEND1),*) RVALUE
                CELL_YQS_VALUE(CELL_STATE_OFFSET(1)+VARNUM-1,i)=RVALUE
C               Read in the parameter's name
                POS3=CLOCAT(' ',LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) POS3=CLOCAT(CHAR(9),
     '            LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) THEN
                  VALUE=LINE(IBEG+POS1+POS2:IEND)
                ELSE
                  VALUE=LINE(IBEG+POS1+POS2:IBEG+POS1+POS2+POS3-1)
                ENDIF
                CALL STRING_TRIM(VALUE,IBEG1,IEND1)
                IF(IBEG1-IEND1.GT.CELL_NAME_LENGTH) IEND1=IBEG1+
     '            CELL_NAME_LENGTH
                CELL_YQS_NAMES(CELL_STATE_OFFSET(1)+VARNUM-1,i) =
     '            VALUE(IBEG1:IEND1)
              ELSE
                ERROR='>>'//VALUE(IBEG1:IEND1)
     '            //' is an invalid real number'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>Variable number not found'
              GOTO 9999
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_STATE
      IF (CELL_NUM_MODEL(1).GT.0) THEN
C       Skip blank line
        READ(IFILE,'(A)') LINE
C       Skip model data line
        READ(IFILE,'(A)') LINE
        DO i=1,CELL_NUM_VARIANTS
C         Skip cell variant number line
          READ(IFILE,'(A)') LINE
          DO j=1,CELL_NUM_MODEL(1)
            READ(IFILE,'(A)') LINE
            CALL STRING_TRIM(LINE,IBEG,IEND)
            POS1=CLOCAT(' ',LINE(IBEG:IEND))
            IF(POS1.EQ.0) POS1=CLOCAT(CHAR(9),LINE(IBEG:IEND))
            IF(POS1.NE.0) THEN
              VARIABLE=LINE(IBEG:IBEG+POS1-1)
              CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
              IF(VARIABLE(IEND1:IEND1).EQ.'*') THEN
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1-1))
                CELL_ICQS_SPATIAL(CELL_MODEL_OFFSET(1)+VARNUM-1,i)=1
              ELSE
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1))
                CELL_ICQS_SPATIAL(CELL_MODEL_OFFSET(1)+VARNUM-1,i)=0
              ENDIF
              POS2=CLOCAT(' ',LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) POS2=CLOCAT(CHAR(9),LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) THEN
                VALUE=LINE(IBEG+POS1:IEND)
              ELSE
                VALUE=LINE(IBEG+POS1:IBEG+POS1+POS2-1)
              ENDIF
              CALL STRING_TRIM(VALUE,IBEG1,IEND1)
              IF(IVALID(VALUE(IBEG1:IEND1))) THEN
                READ(VALUE(IBEG1:IEND1),*) IVALUE
                CELL_ICQS_VALUE(CELL_MODEL_OFFSET(1)+VARNUM-1,i)=IVALUE
C               Read in the parameter's name
                POS3=CLOCAT(' ',LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) POS3=CLOCAT(CHAR(9),
     '            LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) THEN
                  VALUE=LINE(IBEG+POS1+POS2:IEND)
                ELSE
                  VALUE=LINE(IBEG+POS1+POS2:IBEG+POS1+POS2+POS3-1)
                ENDIF
                CALL STRING_TRIM(VALUE,IBEG1,IEND1)
                IF(IBEG1-IEND1.GT.CELL_NAME_LENGTH) IEND1=IBEG1+
     '            CELL_NAME_LENGTH
                CELL_ICQS_NAMES(CELL_MODEL_OFFSET(1)+VARNUM-1,i) =
     '            VALUE(IBEG1:IEND1)
              ELSE
                ERROR='>>'//VALUE(IBEG1:IEND1)
     '            //' is an invalid integer number'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>Variable number not found'
              GOTO 9999
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_MODEL
      IF (CELL_NUM_CONTROL(1).GT.0) THEN
C       Skip blank line
        READ(IFILE,'(A)') LINE
C       Skip control data line
        READ(IFILE,'(A)') LINE
        DO i=1,CELL_NUM_VARIANTS
C         Skip cell variant number line
          READ(IFILE,'(A)') LINE
          DO j=1,CELL_NUM_CONTROL(1)
            READ(IFILE,'(A)') LINE
            CALL STRING_TRIM(LINE,IBEG,IEND)
            POS1=CLOCAT(' ',LINE(IBEG:IEND))
            IF(POS1.EQ.0) POS1=CLOCAT(CHAR(9),LINE(IBEG:IEND))
            IF(POS1.NE.0) THEN
              VARIABLE=LINE(IBEG:IBEG+POS1-1)
              CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
              IF(VARIABLE(IEND1:IEND1).EQ.'*') THEN
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1-1))
                CELL_ICQS_SPATIAL(CELL_CONTROL_OFFSET(1)+VARNUM-1,i)=1
              ELSE
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1))
                CELL_ICQS_SPATIAL(CELL_CONTROL_OFFSET(1)+VARNUM-1,i)=0
              ENDIF
              POS2=CLOCAT(' ',LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) POS2=CLOCAT(CHAR(9),LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) THEN
                VALUE=LINE(IBEG+POS1:IEND)
              ELSE
                VALUE=LINE(IBEG+POS1:IBEG+POS1+POS2-1)
              ENDIF
              CALL STRING_TRIM(VALUE,IBEG1,IEND1)
              IF(IVALID(VALUE(IBEG1:IEND1))) THEN
                READ(VALUE(IBEG1:IEND1),*) IVALUE
                CELL_ICQS_VALUE(CELL_CONTROL_OFFSET(1)+VARNUM-1,i)=
     '            IVALUE
C               Read in the parameter's name
                POS3=CLOCAT(' ',LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) POS3=CLOCAT(CHAR(9),
     '            LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) THEN
                  VALUE=LINE(IBEG+POS1+POS2:IEND)
                ELSE
                  VALUE=LINE(IBEG+POS1+POS2:IBEG+POS1+POS2+POS3-1)
                ENDIF
                CALL STRING_TRIM(VALUE,IBEG1,IEND1)
                IF(IBEG1-IEND1.GT.CELL_NAME_LENGTH) IEND1=IBEG1+
     '            CELL_NAME_LENGTH
                CELL_ICQS_NAMES(CELL_CONTROL_OFFSET(1)+VARNUM-1,i) =
     '            VALUE(IBEG1:IEND1)
              ELSE
                ERROR='>>'//VALUE(IBEG1:IEND1)
     '            //' is an invalid integer number'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>Variable number not found'
              GOTO 9999
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_CONTROL
      IF (CELL_NUM_PARAMETERS(1).GT.0) THEN
C       Skip blank line
        READ(IFILE,'(A)') LINE
C       Skip material data line
        READ(IFILE,'(A)') LINE
        DO i=1,CELL_NUM_VARIANTS
C         Skip cell variant number line
          READ(IFILE,'(A)') LINE
          DO j=1,CELL_NUM_PARAMETERS(1)
            READ(IFILE,'(A)') LINE
            CALL STRING_TRIM(LINE,IBEG,IEND)
            POS1=CLOCAT(' ',LINE(IBEG:IEND))
            IF(POS1.EQ.0) POS1=CLOCAT(CHAR(9),LINE(IBEG:IEND))
            IF(POS1.NE.0) THEN
              VARIABLE=LINE(IBEG:IBEG+POS1-1)
              CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
              IF(VARIABLE(IEND1:IEND1).EQ.'*') THEN
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1-1))
                CELL_RCQS_SPATIAL(CELL_PARAMETERS_OFFSET(1)+VARNUM-1,i)=
     '            1
              ELSE
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1))
                CELL_RCQS_SPATIAL(CELL_PARAMETERS_OFFSET(1)+VARNUM-1,i)=
     '            0
              ENDIF
              POS2=CLOCAT(' ',LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) POS2=CLOCAT(CHAR(9),LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) THEN
                VALUE=LINE(IBEG+POS1:IEND)
              ELSE
                VALUE=LINE(IBEG+POS1:IBEG+POS1+POS2-1)
              ENDIF
              CALL STRING_TRIM(VALUE,IBEG1,IEND1)
              IF(RVALID(VALUE(IBEG1:IEND1))) THEN
                READ(VALUE(IBEG1:IEND1),*) RVALUE
                CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(1)+VARNUM-1,i)=
     '            RVALUE
C               Read in the parameter's name
                POS3=CLOCAT(' ',LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) POS3=CLOCAT(CHAR(9),
     '            LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) THEN
                  VALUE=LINE(IBEG+POS1+POS2:IEND)
                ELSE
                  VALUE=LINE(IBEG+POS1+POS2:IBEG+POS1+POS2+POS3-1)
                ENDIF
                CALL STRING_TRIM(VALUE,IBEG1,IEND1)
                IF(IBEG1-IEND1.GT.CELL_NAME_LENGTH) IEND1=IBEG1+
     '            CELL_NAME_LENGTH
                CELL_RCQS_NAMES(CELL_PARAMETERS_OFFSET(1)+VARNUM-1,i) =
     '            VALUE(IBEG1:IEND1)
              ELSE
                ERROR='>>'//VALUE(IBEG1:IEND1)
     '            //' is an invalid real number'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>Variable number not found'
              GOTO 9999
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_PARAMETERS
      IF (CELL_NUM_PROTOCOL(1).GT.0) THEN
C       Skip blank line
        READ(IFILE,'(A)') LINE
C       Skip protocol data line
        READ(IFILE,'(A)') LINE
        DO i=1,CELL_NUM_VARIANTS
C         Skip cell variant number line
          READ(IFILE,'(A)') LINE
          DO j=1,CELL_NUM_PROTOCOL(1)
            READ(IFILE,'(A)') LINE
            CALL STRING_TRIM(LINE,IBEG,IEND)
            POS1=CLOCAT(' ',LINE(IBEG:IEND))
            IF(POS1.EQ.0) POS1=CLOCAT(CHAR(9),LINE(IBEG:IEND))
            IF(POS1.NE.0) THEN
              VARIABLE=LINE(IBEG:IBEG+POS1-1)
              CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
              IF(VARIABLE(IEND1:IEND1).EQ.'*') THEN
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1-1))
                CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(1)+VARNUM-1,i)=1
              ELSE
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1))
                CELL_RCQS_SPATIAL(CELL_PROTOCOL_OFFSET(1)+VARNUM-1,i)=0
              ENDIF
              POS2=CLOCAT(' ',LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) POS2=CLOCAT(CHAR(9),LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) THEN
                VALUE=LINE(IBEG+POS1:IEND)
              ELSE
                VALUE=LINE(IBEG+POS1:IBEG+POS1+POS2-1)
              ENDIF
              CALL STRING_TRIM(VALUE,IBEG1,IEND1)
              IF(RVALID(VALUE(IBEG1:IEND1))) THEN
                READ(VALUE(IBEG1:IEND1),*) RVALUE
                CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(1)+VARNUM-1,i)=
     '            RVALUE
C               Read in the parameter's name
                POS3=CLOCAT(' ',LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) POS3=CLOCAT(CHAR(9),
     '            LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) THEN
                  VALUE=LINE(IBEG+POS1+POS2:IEND)
                ELSE
                  VALUE=LINE(IBEG+POS1+POS2:IBEG+POS1+POS2+POS3-1)
                ENDIF
                CALL STRING_TRIM(VALUE,IBEG1,IEND1)
                IF(IBEG1-IEND1.GT.CELL_NAME_LENGTH) IEND1=IBEG1+
     '            CELL_NAME_LENGTH
                CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(1)+VARNUM-1,i) =
     '            VALUE(IBEG1:IEND1)
              ELSE
                ERROR='>>'//VALUE(IBEG1:IEND1)
     '            //' is an invalid real number'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>Variable number not found'
              GOTO 9999
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_PROTOCOL
      IF (CELL_NUM_AII(1).GT.0) THEN
C       Skip blank line
        READ(IFILE,'(A)') LINE
C       Skip aii data line
        READ(IFILE,'(A)') LINE
        DO i=1,CELL_NUM_VARIANTS
C         Skip cell variant number line
          READ(IFILE,'(A)') LINE
          DO j=1,CELL_NUM_AII(1)
            READ(IFILE,'(A)') LINE
            CALL STRING_TRIM(LINE,IBEG,IEND)
            POS1=CLOCAT(' ',LINE(IBEG:IEND))
            IF(POS1.EQ.0) POS1=CLOCAT(CHAR(9),LINE(IBEG:IEND))
            IF(POS1.NE.0) THEN
              VARIABLE=LINE(IBEG:IBEG+POS1-1)
              CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
              IF(VARIABLE(IEND1:IEND1).EQ.'*') THEN
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1-1))
                CELL_ICQS_SPATIAL(CELL_AII_OFFSET(1)+VARNUM-1,i)=1
              ELSE
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1))
                CELL_ICQS_SPATIAL(CELL_AII_OFFSET(1)+VARNUM-1,i)=0
              ENDIF
              POS2=CLOCAT(' ',LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) POS2=CLOCAT(CHAR(9),LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) THEN
                VALUE=LINE(IBEG+POS1:IEND)
              ELSE
                VALUE=LINE(IBEG+POS1:IBEG+POS1+POS2-1)
              ENDIF
              CALL STRING_TRIM(VALUE,IBEG1,IEND1)
              IF(IVALID(VALUE(IBEG1:IEND1))) THEN
                READ(VALUE(IBEG1:IEND1),*) IVALUE
                CELL_ICQS_VALUE(CELL_AII_OFFSET(1)+VARNUM-1,i)=IVALUE
C               Read in the parameter's name
                POS3=CLOCAT(' ',LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) POS3=CLOCAT(CHAR(9),
     '            LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) THEN
                  VALUE=LINE(IBEG+POS1+POS2:IEND)
                ELSE
                  VALUE=LINE(IBEG+POS1+POS2:IBEG+POS1+POS2+POS3-1)
                ENDIF
                CALL STRING_TRIM(VALUE,IBEG1,IEND1)
                IF(IBEG1-IEND1.GT.CELL_NAME_LENGTH) IEND1=IBEG1+
     '            CELL_NAME_LENGTH
                CELL_ICQS_NAMES(CELL_AII_OFFSET(1)+VARNUM-1,i) =
     '            VALUE(IBEG1:IEND1)
              ELSE
                ERROR='>>'//VALUE(IBEG1:IEND1)
     '            //' is an invalid integer number'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>Variable number not found'
              GOTO 9999
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_AII
      IF (CELL_NUM_ARI(1).GT.0) THEN
C       Skip blank line
        READ(IFILE,'(A)') LINE
C       Skip ari data line
        READ(IFILE,'(A)') LINE
        DO i=1,CELL_NUM_VARIANTS
C       Skip cell variant number line
          READ(IFILE,'(A)') LINE
          DO j=1,CELL_NUM_ARI(1)
            READ(IFILE,'(A)') LINE
            CALL STRING_TRIM(LINE,IBEG,IEND)
            POS1=CLOCAT(' ',LINE(IBEG:IEND))
            IF(POS1.EQ.0) POS1=CLOCAT(CHAR(9),LINE(IBEG:IEND))
            IF(POS1.NE.0) THEN
              VARIABLE=LINE(IBEG:IBEG+POS1-1)
              CALL STRING_TRIM(VARIABLE,IBEG1,IEND1)
              IF(VARIABLE(IEND1:IEND1).EQ.'*') THEN
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1-1))
                CELL_RCQS_SPATIAL(CELL_ARI_OFFSET(1)+VARNUM-1,i)=1
              ELSE
                VARNUM=IFROMC(VARIABLE(IBEG1:IEND1))
                CELL_RCQS_SPATIAL(CELL_ARI_OFFSET(1)+VARNUM-1,i)=0
              ENDIF
              POS2=CLOCAT(' ',LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) POS2=CLOCAT(CHAR(9),LINE(IBEG+POS1:IEND))
              IF(POS2.EQ.0) THEN
                VALUE=LINE(IBEG+POS1:IEND)
              ELSE
                VALUE=LINE(IBEG+POS1:IBEG+POS1+POS2-1)
              ENDIF
              CALL STRING_TRIM(VALUE,IBEG1,IEND1)
              IF(RVALID(VALUE(IBEG1:IEND1))) THEN
                READ(VALUE(IBEG1:IEND1),*) RVALUE
                CELL_RCQS_VALUE(CELL_ARI_OFFSET(1)+VARNUM-1,i)=RVALUE
C               Read in the parameter's name
                POS3=CLOCAT(' ',LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) POS3=CLOCAT(CHAR(9),
     '            LINE(IBEG+POS1+POS2+1:IEND))
                IF(POS3.EQ.0) THEN
                  VALUE=LINE(IBEG+POS1+POS2:IEND)
                ELSE
                  VALUE=LINE(IBEG+POS1+POS2:IBEG+POS1+POS2+POS3-1)
                ENDIF
                CALL STRING_TRIM(VALUE,IBEG1,IEND1)
                IF(IBEG1-IEND1.GT.CELL_NAME_LENGTH) IEND1=IBEG1+
     '            CELL_NAME_LENGTH
                CELL_RCQS_NAMES(CELL_ARI_OFFSET(1)+VARNUM-1,i) =
     '            VALUE(IBEG1:IEND1)
              ELSE
                ERROR='>>'//VALUE(IBEG1:IEND1)
     '            //' is an invalid real number'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>Variable number not found'
              GOTO 9999
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDIF !CELL_NUM_ARI

      IF(DOP) THEN
        DO i=1,CELL_NUM_VARIANTS
          WRITE(OP_STRING,'('' Parameter names for variant '',I3)') i
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' CELL_YQS_NAMES: '',4A)')
     '      (CELL_YQS_NAMES(j,i),j=1,CELL_NUM_STATE(1))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' CELL_ICQS_NAMES: '',4A)')
     '      (CELL_ICQS_NAMES(j,i),j=1,CELL_NUM_MODEL(1)
     '      +CELL_NUM_CONTROL(1)+CELL_NUM_AII(1)+1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' CELL_RCQS_NAMES: '',4A)')
     '      (CELL_RCQS_NAMES(j,i),j=1,CELL_NUM_PARAMETERS(1)+
     '      CELL_NUM_PROTOCOL(1)+CELL_NUM_ARI(1))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      DO nq=NQR(1,nr),NQR(2,nr)
C *** DPN 14 March 2000 - fixing up units
c       CQ(1,nq)=CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(1)+Cm-1,1)*1.0d3
        CQ(1,nq)=CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(1)+Cm-1,1)
        CQ(2,nq)=CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(1)+Am-1,1)
      ENDDO

      CALL EXITS('IPCELL_READ')
      RETURN
 9999 CALL ERRORS('IPCELL_READ',ERROR)
      CALL EXITS('IPCELL_READ')
      RETURN 1
      END


