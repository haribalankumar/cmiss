      SUBROUTINE EVEVENT(NBJ,NDDATA,NELEC,NTIME,LD,PHI_H,
     '  WD,XID,ZD,STRING,ERROR,*)

C#### Subroutine: EVEVENT
C###  Description:
C###    <HTML>
C###    EVELEC evaluates event times from signals stored in the PHI_H
C###    array using the unemap event detection algorithms.
C###    <BR><BR>
C###    Requires libunemap.a to be linked in.
C###    </HTML>

C***     Created by Leo Cheng on June 2002

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'sign00.cmn'
      INCLUDE 'tol00.cmn'

 !    Parameter List
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),NELEC,NTIME
      REAL*8 PHI_H(NY_TRANSFER_M,NTSM),WD(NJM,NDM),
     '  XID(NIM,NDM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER CERROR(50),CERRLEN,EDA,ERR,IBEG,IBEG1,IEND,IEND1,
     &  OBJECTIVE_TYPE,THRESHOLD,N3CO,noelec,noevent,nr,NUMTIMEDATA,nts,
     '  SEPARATION,WIDTH
      REAL*4 LEVEL,OBJECTIVE(NTSM)
      REAL*8 FREQUENCY,SIGNALMAX(9),SIGNALMIN(9),TIME
      CHARACTER FILEFORMAT*6,SIGNALFILE*(MXCH)
      LOGICAL CBBREV,DEBUG,ENDFILE
      ! Save event times for the given electrode
      LOGICAL SAVE_EVENT_TIMES
      CHARACTER VARIABLE_NAME*(MXCH),SAVE_STRING*(MXCH)
      INTEGER LEN_TRIM,SAVE_ELECTRODE

!     Functions
      INTEGER IFROMC
      REAL*8 RFROMC

C***  Declaration of the event times and the number of event times
      INTEGER NEVENTSM
      PARAMETER(NEVENTSM=10)
      INTEGER EVENT_TIMES(NEVENTSM,NELEC)
      INTEGER EVENTS_FOUND(NELEC)

C*** Declaration of strings and constants for the detection algorithms
      INTEGER EDA_THRESHOLD,EDA_INTERVAL,EDA_LEVEL
      INTEGER ABSOLUTE_SLOPE,NEGATIVE_SLOPE,POSITIVE_SLOPE,
     '  ABSOLUTE_VALUE,NEGATIVE_VALUE,POSITIVE_VALUE
      CHARACTER EDA_LABELS(0:2)*9,OBJECTIVE_LABELS(0:5)*14

      PARAMETER(EDA_INTERVAL=0,EDA_THRESHOLD=1,EDA_LEVEL=2)
      PARAMETER(
     '  ABSOLUTE_SLOPE=0,NEGATIVE_SLOPE=1,POSITIVE_SLOPE=2,
     '  ABSOLUTE_VALUE=3,NEGATIVE_VALUE=4,POSITIVE_VALUE=5)

      DATA EDA_LABELS /
     '  'Interval', 'Threshold', 'Level'/
      DATA OBJECTIVE_LABELS /
     '  'Absolute Slope','Negative Slope','Positive Slope',
     '  'Absolute Value','Negative Value','Positive Value'/


      CALL ENTERS('EVEVENT',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C--------------------------------------------------------------------

C#### Command: FEM evaluate events<;FILENAME>
C###  Parameter:        <method (interval/level/threshold)[interval]>
C###    The algorithm to search for event times.
C###  Parameter:        <objective (slope/magnitude)[slope]>
C###    The objective function.
C###  Parameter:        <slope (absolute/positive/negative)[negative]>
C###    The value of the objective function.
C###  Parameter:        <width (#)[3]>
C###    The average width of the finite difference interval. Used with
C###    all algorithms
C###  Parameter:        <level (#)[0]>
C###    Level to select event. Used the LEVEL algorithm
C###  Parameter:        <percentage (#)[90]>
C###    Threshold percentage. Used the THRESHOLD algorithm
C###  Parameter:        <separation (#)[100]>
C###    Minimum separation between events. Used with the THRESHOLD
C###    algorithm
C###  Parameter:        <frequency (#)[1.0]>
C###    The frequency of the signal information is kHz.
C###  Parameter:        <save (VARIABLE_NAME)[]>
C###    Give the interpreter variable name to store the event
C###    times in a space separated format
C###  Parameter:          <electrode (#)[1]>
C###    When saving event times, give the electrode to save.
C###  Parameter:        <(ascii/binary)[ascii]>
C###    Output the events files as ascii or binary format.
C###  Description:
C###    Evaluates event times from the PHI array
C###    using a moving finite difference
C###    window to approximate the derivatives.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//
     '    '<method (interval/level/threshold)[interval]>'
        OP_STRING(3)=BLANK(1:15)//
     '    '<objective (value/slope) [slope]>'
        OP_STRING(4)=BLANK(1:15)//
     '    '<magnitude (negative/positive/value) [negative]>'
        OP_STRING(5)=BLANK(1:15)//'<width #[3]>'
        OP_STRING(6)=BLANK(1:15)//'<frequency #[1.0]>'
        OP_STRING(7)=BLANK(1:15)//'<level #[0]>'
        OP_STRING(8)=BLANK(1:15)//'<percentage #[90]>'
        OP_STRING(9)=BLANK(1:15)//'<separation #[100]>'
        OP_STRING(10)=BLANK(1:15)//'<save VARIABLE_NAME[]>'
        OP_STRING(11)=BLANK(1:15)//'  <electrode #[1]>'
        OP_STRING(12)=BLANK(1:15)//'<(ascii/binary) [ascii]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVELEC',ERROR,*9999)
      ELSE

        IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'INTERVAL',3,noco+2,NTCO,N3CO)) THEN
            EDA=EDA_INTERVAL
          ELSEIF(CBBREV(CO,'THRESHOLD',3,noco+2,NTCO,N3CO)) THEN
            EDA=EDA_THRESHOLD
          ELSEIF(CBBREV(CO,'LEVEL',3,noco+2,NTCO,N3CO)) THEN
            EDA=EDA_LEVEL
          ELSE
            ERROR="Unknown detection algorithm"
            GOTO 9999
          ENDIF
        ENDIF

        CALL ASSERT(NELEC.GT.0,'No electrodes exist'
     '    ,ERROR,*9999)
        CALL ASSERT(NTIME.GT.3,
     '    'There must be more than 3 timesteps',ERROR,*9999)

C*** Fileformat
        IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF

        IF(CBBREV(CO,'FREQENCY',3,noco+1,NTCO,N3CO)) THEN
          FREQUENCY=RFROMC(CO(N3CO+1))
        ELSE
          FREQUENCY=1.D0
        ENDIF

        IF(NTCOQU(noco).GT.0) THEN
          SIGNALFILE=COQU(noco,1)
        ELSE
          SIGNALFILE='events'
        ENDIF

        ERR=0
        DEBUG=.FALSE.
C        DEBUG=.TRUE.
C        ISIZE_PHI(1)=1

        IF(CBBREV(CO,'METHOD_OF_DETECTION',3,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'INTERVAL',3,noco+2,NTCO,N3CO)) THEN
            EDA=EDA_INTERVAL
          ELSEIF(CBBREV(CO,'THRESHOLD',3,noco+2,NTCO,N3CO)) THEN
            EDA=EDA_THRESHOLD
          ELSEIF(CBBREV(CO,'LEVEL',3,noco+2,NTCO,N3CO)) THEN
            EDA=EDA_LEVEL
          ELSE
            ERROR="Unknown detection algorithm"
            GOTO 9999
          ENDIF
        ELSE
          EDA=EDA_INTERVAL
        ENDIF

        IF(CBBREV(CO,'LEVEL',3,noco+1,NTCO,N3CO)) THEN
          LEVEL=SNGL(RFROMC(CO(N3CO+1)))
        ELSE
          LEVEL=0.0
        ENDIF

        IF(CBBREV(CO,'PERCENTAGE',3,noco+1,NTCO,N3CO)) THEN
          THRESHOLD=IFROMC(CO(N3CO+1))
        ELSE
          THRESHOLD=90
        ENDIF

        IF(CBBREV(CO,'SEPARATION',3,noco+1,NTCO,N3CO)) THEN
          SEPARATION=IFROMC(CO(N3CO+1))
        ELSE
          SEPARATION=100
        ENDIF

        IF(CBBREV(CO,'WIDTH',3,noco+1,NTCO,N3CO)) THEN
          WIDTH=IFROMC(CO(N3CO+1))
        ELSE
          WIDTH=4
        ENDIF

        IF(CBBREV(CO,'SAVE',4,noco+1,NTCO,N3CO)) THEN
          SAVE_EVENT_TIMES=.TRUE.
          VARIABLE_NAME=CO(N3CO+1)
          IF(LEN_TRIM(VARIABLE_NAME).EQ.0) THEN
            ERROR='Need to give a variable name to save to'
            GOTO 9999
          ENDIF
          IF(CBBREV(CO,'ELECTRODE',4,noco+1,NTCO,N3CO)) THEN
            SAVE_ELECTRODE=IFROMC(CO(N3CO+1))
          ELSE
            SAVE_ELECTRODE=1
          ENDIF
        ELSE
          SAVE_EVENT_TIMES=.FALSE.
        ENDIF

        IF(CBBREV(CO,'OBJECTIVE',3,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'SLOPE',3,noco+2,NTCO,N3CO)) THEN
            IF(CBBREV(CO,'MAGNITUDE',3,noco+1,NTCO,N3CO)) THEN

              IF(CBBREV(CO,'ABSOLUTE',3,noco+2,NTCO,N3CO)) THEN
                OBJECTIVE_TYPE=ABSOLUTE_SLOPE
              ELSEIF(CBBREV(CO,'NEGATIVE',3,noco+2,NTCO,N3CO)) THEN
                OBJECTIVE_TYPE=NEGATIVE_SLOPE
              ELSEIF(CBBREV(CO,'POSITIVE',3,noco+2,NTCO,N3CO)) THEN
                OBJECTIVE_TYPE=POSITIVE_SLOPE
              ELSE
                ERROR='Invalide objective type'
                GOTO 9999
              ENDIF
            ENDIF
          ELSEIF(CBBREV(CO,'VALUE',3,noco+2,NTCO,N3CO)) THEN
            IF(CBBREV(CO,'MAGNITUDE',3,noco+2,NTCO,N3CO)) THEN
              IF(CBBREV(CO,'ABSOLUTE',3,noco+2,NTCO,N3CO)) THEN
                OBJECTIVE_TYPE=ABSOLUTE_VALUE
              ELSEIF(CBBREV(CO,'NEGATIVE',3,noco+2,NTCO,N3CO)) THEN
                OBJECTIVE_TYPE=NEGATIVE_VALUE
              ELSEIF(CBBREV(CO,'POSITIVE',3,noco+2,NTCO,N3CO)) THEN
                OBJECTIVE_TYPE=POSITIVE_VALUE
              ELSE
                ERROR='Invalide objective type'
                GOTO 9999
              ENDIF
            ENDIF
          ENDIF ! objective/value
        ELSE
          OBJECTIVE_TYPE=NEGATIVE_SLOPE
        ENDIF


        WRITE(OP_STRING(1),'('' Algorithm   : '',A)') EDA_LABELS(EDA)
        WRITE(OP_STRING(2),'('' Width       : '',I6)') WIDTH
        IF(EDA.EQ.EDA_INTERVAL) THEN
          WRITE(OP_STRING(3),'('' Objective   : '',A)')
     '      OBJECTIVE_LABELS(OBJECTIVE_TYPE)
        ELSEIF(EDA.EQ.EDA_THRESHOLD) THEN
          WRITE(OP_STRING(3),'('' Objective   : '',A)')
     '      OBJECTIVE_LABELS(OBJECTIVE_TYPE)
          WRITE(OP_STRING(4),'('' Separation  : '',I6)') SEPARATION
          WRITE(OP_STRING(5),'('' Threshold % : '',I6)') THRESHOLD
        ELSEIF(EDA.EQ.EDA_LEVEL) THEN
          WRITE(OP_STRING(3),'('' Level       : '',F6.2)') LEVEL
        ENDIF
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        WRITE(OP_STRING(1),'('' '')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        DO noelec=1,NELEC
          DO nts=1,NTIME
            OBJECTIVE(nts)=SNGL(PHI_H(noelec,nts))
          ENDDO
C!!! This is here because the c portion is wrong!!!
C          OBJECTIVE(NTIME+1)=-0.00012345

          IF(DEBUG) THEN
            WRITE(*,*) 'Phi vector'
            DO nts=1,NTIME
              WRITE(*,*) nts, OBJECTIVE(nts)
            ENDDO
          ENDIF


          CALL CMUNEMAPWRAPPER(
     '      EDA,
     '      EVENT_TIMES(1,noelec),
     '      EVENTS_FOUND(noelec),
     '      NEVENTSM,
     '      SEPARATION,
     '      NTIME,
     '      NELEC,
     '      OBJECTIVE_TYPE,
     '      WIDTH,
     '      THRESHOLD,
     '      OBJECTIVE,
     '      LEVEL,
     '      SNGL(FREQUENCY)*1000.0,
     '      ERR,
     '      CERROR)

          IF(ERR.NE.0) THEN
            CALL CSTRINGLEN(CERRLEN,CERROR)
            CALL C2FSTRING(CERROR,CERRLEN,ERROR)
            GOTO 9999
          ENDIF

        ENDDO !noelec

        WRITE(OP_STRING,'('' Electrode  #Events  Time'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO noelec=1,NELEC

          WRITE(OP_STRING,'(I8,I8,20F10.1)')
     '      noelec,EVENTS_FOUND(noelec),
     '      (EVENT_TIMES(noevent,noelec)/FREQUENCY,
     '      noevent=1,EVENTS_FOUND(noelec))

          IF(SAVE_EVENT_TIMES.AND.noelec.EQ.SAVE_ELECTRODE) THEN
            !Save the event time string
            WRITE(SAVE_STRING,'(20E15.8)') (EVENT_TIMES(noevent,noelec)
     &        /FREQUENCY,noevent=1,EVENTS_FOUND(noelec))
            CALL STRING_TRIM(VARIABLE_NAME,IBEG,IEND)
            CALL STRING_TRIM(SAVE_STRING,IBEG1,IEND1)
            CALL SET_USER_CHARACTER(VARIABLE_NAME(IBEG:IEND),
     &        SAVE_STRING(IBEG1:IEND1),ERR)
          ENDIF

C*** Fill ZD with the first event time for the signal file
          IF(EVENTS_FOUND(noelec).GT.0) THEN
            ZD(NJT+1,noelec)=EVENT_TIMES(1,noelec)/FREQUENCY
          ELSE
            ZD(NJT+1,noelec)=0.d0
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C LKC 23-SEP-2002 Need to set up the fields
          nr=1
          NJ_LOC(NJL_FIEL,0,nr)=1
          NJ_LOC(NJL_FIEL,1,nr)=4

        ENDDO

C*** Need to write the solutions to a file - use a binsig
C*** file for now, where event time is the potential value.
        SIGNAL_HEADER(IOFILE1)='Computed event times '
        SIGNAL_HOWGENERATED(IOFILE1)=0

        TIME=0.0d0
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,SIGNALFILE,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

        SIGNAL_NUMREGIONS(IOFILE1)=1
        SIGNAL_REGNAME(1,IOFILE1)=' '
        SIGNAL_NUMELEC(1,IOFILE1)=NDT
        SIGNAL_REGTYPE(1,IOFILE1)=0 !irregular rectangular cartesian
        SIGNAL_ELEMLOC(IOFILE1)=0
        SIGNAL_NUMXI(1,IOFILE1)=0
        CALL ASSERT(SIGNAL_NUMELEC(1,IOFILE1).GE.1,
     '    'No electodes to export',ERROR,*9999)

C*** Write electrode data
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,SIGNALFILE,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Write signal data
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,
     '    SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'WRITE',
     '    FILEFORMAT,SIGNALFILE,'SIGNAL_DATA',ENDFILE,
     '    .TRUE.,ERROR,*9999)

C*** Close signal
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,SIGNALFILE,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)

      ENDIF

      CALL EXITS('EVEVENT')
      RETURN
 9999 CALL ERRORS('EVEVENT',ERROR)
      CALL EXITS('EVEVENT')
      RETURN 1
      END


