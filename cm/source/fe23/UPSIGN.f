      SUBROUTINE UPSIGN(LD,NBJ,NDDATA,NEELEM,WD,XID,ZD,
     '  STRING,ERROR,*)

C#### Subroutine: UPSIGN
C###  Description:
C###    UPSIGN updates a signal file ie. rewrites the signal file
C###    with new parameters etc.

C***  LKC 12-JUN-1998 big overhaul
C***  LKC 18-OCT-1998 select a time period

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'sign00.cmn'

!     Parameter List
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     '  NEELEM(0:NE_R_M,0:NRM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),XID_TMP(NIM,NDM),
     '  ZD(NJM,NDM),ZD_TMP(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,N3CO,nb,nr,notime,
     '  NUMTIMEDATA,NUMTIMEDATA1
      REAL*8 SIGNALMAX(9),SIGNALMIN(9),TIME,TEND,TSTART
      CHARACTER*6 FILEFORMAT,INFILE*100,OUTFILE*100,ERROR_DUMMY*255
      LOGICAL CBBREV,ENDFILE

!     Functions
      REAL*8 RFROMC

      CALL ENTERS('UPSIGN',*9999)

      CALL STRING_TRIM(FILE00,IBEG1,IEND1)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update signal
C###  Parameter:      <infile  FILENAME[$current]>
C###  Specify the infile name
C###  Parameter:      <outfile FILENAME[$current_new]>
C###  Specify the outfile name
C###  Parameter:      <tstart  #[beginning]>
C###  Specify the start time
C###  Parameter:      <tend    #[end]>
C###  Specify the finish time
C###  Parameter:      <(ascii/binary)[ascii]>
C###  Specify the file format
C###  Description: updates a signal file ie. rewrites the signal file
C###    with new parameters.
C###

        OP_STRING(1)=STRING(1:IEND)//' <infile FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//' <outfile   FILENAME['
     '    //FILE00(IBEG1:IEND1)//'_new]>'
        OP_STRING(3)=BLANK(1:15)//' <tstart    #[beginning]>'
        OP_STRING(4)=BLANK(1:15)//' <tend      #[end]>'
        OP_STRING(5)=BLANK(1:15)//' <(ascii/binary)[ascii]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPSIGN',ERROR,*9999)
      ELSE
        nr=1 !temp until data is region dependant
        CALL ASSERT(NDDATA(0,nr).GT.0,
     '    '>> Define data first',ERROR,*9999)
        CALL ASSERT(CALC_XI,'>> Calculate Xi first',ERROR,*9999)

        IF(CBBREV(CO,'INFILE',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          INFILE=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          INFILE=FILE00(IBEG1:IEND1)
        ENDIF

        IF(CBBREV(CO,'OUTFILE',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          OUTFILE=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          OUTFILE=FILE00(IBEG1:IEND1)//'_new'
        ENDIF

        CALL STRING_TRIM(INFILE,IBEG,IEND)
        CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)
        IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
     '    INFILE(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
          ERROR='>>Input and output file names must be different'
          GOTO 9999
        ENDIF

C LKC 12-JUN-1998 These are just repeats ?
C
C        CALL STRING_TRIM(INFILE,IBEG,IEND)
C        CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)
C        IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
C     '    INFILE(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
C          ERROR='>>Input and output file names must be different'
C          GOTO 9999
C        ENDIF
C        CALL STRING_TRIM(INFILE,IBEG,IEND)
C        CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)
C        IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
C     '    INFILE(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
C          ERROR='>>Input and output file names must be different'
C          GOTO 9999
C        ENDIF

        IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF

        IF(CBBREV(CO,'TSTART',2,noco+1,NTCO,N3CO)) THEN
          TSTART=RFROMC(CO(N3CO+1))
        ELSE
          TSTART=-RMAX
        ENDIF

        IF(CBBREV(CO,'TEND',2,noco+1,NTCO,N3CO)) THEN
          TEND=RFROMC(CO(N3CO+1))
        ELSE
          TEND=RMAX
        ENDIF

C cpb 7/11/95 This is only implemented quickly at the moment. As a
C result its only function is to allow the updating of xi positions
C in a binary signal file. Other things may follow.

C*** Open up the input signal file and reset to beginning of signal
C*** data

C MLB time was not set
        TIME=0.0d0

        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C LKC 9-FEB-1998 read ZD and XID into temp arrays to preserve
C   existing information
C
C  LKC 12-JUN-1998 new read - previously no electrode data was read !?!
C        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
C     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
C     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID_TMP,ZD_TMP,'READ',FILEFORMAT,INFILE,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)



C***    IF ASCII need to calculate NUMTIMEDATA
        IF(FILEFORMAT.EQ.'ASCII') THEN
          NUMTIMEDATA=0
          ENDFILE=.FALSE.
          DO WHILE(.NOT.ENDFILE)
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
          ENDDO
        ENDIF

        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C*** Open up the ouput signal file
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C*** Write out the electrode data
        SIGNAL_HOWGENERATED(IOFILE2)=SIGNAL_HOWGENERATED(IOFILE1)
        SIGNAL_HEADER(IOFILE2)=SIGNAL_HEADER(IOFILE1)
        SIGNAL_NUMREGIONS(IOFILE2)=SIGNAL_NUMREGIONS(IOFILE1)



C LKC 31-AUG-1999 Trying to accomadate multiple regions
C   If multiple regions assumes that the electrodes are
C   sequential from region 1 up to region N.

        IF(CALC_XI) THEN

C LKC 31-AUG-1999 Trying to accomadate multiple regions
C          CALL ASSERT(SIGNAL_NUMREGIONS(IOFILE2).EQ.1,
C     '      '>> More the 1 reg, data has no reg assoc.',ERROR,*9999)
C            nb=NBJ(1,NEELEM(1,1)) !nr=1

          SIGNAL_ELEMLOC(IOFILE2)=1
          DO nr=1,SIGNAL_NUMREGIONS(IOFILE2)

C LKC 31-AUG-1999 Trying to accomadate multiple regions
C            CALL ASSERT(NDT.EQ.SIGNAL_NUMELEC(nr,IOFILE1),
C     '        '>> Num Data <> Num Electrodes',ERROR,*9999)

            CALL ASSERT(NDDATA(0,nr).EQ.SIGNAL_NUMELEC(nr,IOFILE1),
     '        '>> Num Data <> Num Electrodes',ERROR,*9999)

            nb=NBJ(1,NEELEM(1,nr))
            SIGNAL_NUMXI(nr,IOFILE2)=NIT(nb)
          ENDDO !nr
        ENDIF
        DO nr=1,SIGNAL_NUMREGIONS(IOFILE2)
          SIGNAL_NUMELEC(nr,IOFILE2)=SIGNAL_NUMELEC(nr,IOFILE1)
          SIGNAL_REGNAME(nr,IOFILE2)=SIGNAL_REGNAME(nr,IOFILE1)
          SIGNAL_REGTYPE(nr,IOFILE2)=SIGNAL_REGTYPE(nr,IOFILE1)
        ENDDO !nr

        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Loop over the signals

        DO notime=1,NUMTIMEDATA

          IF(TIME.GE.TSTART.AND.TIME.LE.TEND) THEN
C***        Read in the signal values
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C***        Write out the signal data
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
          ENDIF ! tstart <= time <= tend

        ENDDO !notime

C*** Close input and output files
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,OUTFILE,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)

      ENDIF

      CALL EXITS('UPSIGN')
      RETURN
 9999 CALL ERRORS('UPSIGN',ERROR)
C LKC 12-JUN-1998
C      IF(ISBINFILEOPEN(IOFILE1))
C     '  CALL BINARYCLOSEFILE(IOFILE1,ERR,CERROR)
C      IF(ISBINFILEOPEN(IOFILE2))
C     '  CALL BINARYCLOSEFILE(IOFILE2,ERR,CERROR)
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE,
     '  ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9998)
 9998 CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,OUTFILE,
     '  ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9997)
 9997 CALL EXITS('UPSIGN')
      RETURN 1
      END


