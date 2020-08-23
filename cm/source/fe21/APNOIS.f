      SUBROUTINE APNOIS(LD,NBJ,NDDATA,WD,XID,ZD,STRING,ERROR,*)
C#### SUBROUTINE: APNOIS
C###  Description:
C###    <HTML>
C###    APNOIS applys random noise to a signal file. The noise/error
C###    can either be applied to the signal values, electode
C###    locations or both.
C###    <BR><BR>
C###    NOISE_LEV is the RMS noise power in millivolts.
C###    NOIS_DIST is average distance the electrodes are displaced.
C###
C###    <BR><BR>
C###    The IO files are stored such that the
C###    <UL>
C###      <LI> output 'noisy' signal file is IOFILE1
C###      <LI> input 'clean' signal file is IOFILE2
C###      <LI> summary (if applicable) is IOFILE3
C###    </UL>
C###    The noise is applied to the clean signal and is stored
C###       in ZD(NJT+1,nd) and
C###    The noisy itself placed in  ZD(NJT+2,nd)
C###    </HTML>
C***  Create By   : LKC 12 September 1997
C***  Last Change : 5-OCT-97 fixing sd calculations
C***              : 8-OCT-97 adding new noise generators
C***              : 12-JAN-98 errors for elec locats.
C***              : 21-AUG-00 error for sd value.
C***              : 25-AUG-00 adding random start option

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'nois00.cmn'
      INCLUDE 'sign00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Paramter List
      INTEGER
     '  LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM)
      REAL*8 WD(NJM,NDM),
     '  XID(NIM,NDM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,INIT_SEED1,INIT_SEED2,INIT_SEED3,
     '  N3CO,nd,nj,noelec,notime,NUMTIMEDATA,NUMTIMEDATA1,nodata,
     '  nr,NRAND,SEED1,SEED2,SEED3
      CHARACTER FILEFORMAT*6,INFILE*(MXCH),OUTFILE*(MXCH),LEV*12,
     '  ERROR1*255,SUMMARYFILE*(MXCH)
      REAL*8 ELECT_ERR(NJM,6),MEAN(3),SD,SDTEMP,RMS(2),SIGNALMAX(9),
     '  SIGNALMIN(9),TIME,TMP
      LOGICAL APPLYTODATA,CBBREV,ENDFILE,FIRST,OPFILE,RANDOM_START
      SAVE FIRST,SEED1,SEED2,SEED3
      DATA FIRST / .TRUE. /

!     Functions
      REAL*8 RAN0,RAN1,RAN10

      CALL ENTERS('APNOIS',*9999)
      CALL STRING_TRIM(FILE00,IBEG,IEND)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM apply noise<;FILENAME>
C###  Description:
C###    Applies random noise to a signal file. The noise/error
C###    can be applied to either the signal values, electrode
C###    locations or both.
C###
C###    Summary information
C###    about the noise is displayed either the screen or file
C###    as specified by <;FILENAME> tag. The outfile
C###       <outfile FILENAME[$current_nsLEV]>
C###    contains the 'dirty' signal. This is written either to a
C###    .binsig or .ipsign file according to the
C###    file format specified
C###
C###  Parameter:      <signal FILENAME[$current]>
C###    Input signal filename
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the element file is stored as a binary or
C###    ascii file.
C###  Parameter:      <outfile FILENAME[$current_nsLEV]>
C###    Specifies the output filename with noise added
C###  Parameter:      <set_random_start>
C###    Specifies whether the starting seeds are not set to 1 by default

        OP_STRING(1)=STRING(1:IEND1)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<signal FILENAME['//
     '    FILE00(IBEG:IEND)//']>'
        OP_STRING(3)=BLANK(1:15)//'<outfile FILENAME['//
     '    FILE00(IBEG:IEND)//'_nsLVL]>'
        OP_STRING(4)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(5)=BLANK(1:15)//'<set_random_start>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM apply noise<;FILENAME> data
C###  Description:
C###    Applies random noise to a data set.

        OP_STRING(1)=STRING(1:IEND1)//'<;FILENAME> data'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','ADDSIG',ERROR,*9999)
      ELSE
        CALL ASSERT(CALL_NOIS,'>>Define noise first',ERROR,*9999)
        CALL ASSERT(DEFNOIS_LOC.LE.3,'>>Update code ',ERROR,*9999)
        CALL ASSERT(DEFNOIS_LOC.GE.1,'>>Update code ',ERROR,*9999)

        IF(CBBREV(CO,'DATA',3,noco+1,NTCO,N3CO)) THEN
          APPLYTODATA=.TRUE.
        ELSE
          APPLYTODATA=.FALSE.
          CALL ASSERT(NJM.GE.NJT+2,'>>Increase NJM',ERROR,*9999)
        ENDIF


C*** Warning LKC JAN-98
        IF(DEFNOIS_TYPE.EQ.2.OR.DEFNOIS_TYPE.EQ.3)THEN
          OP_STRING(1)='Warning : Coefficients not set'
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        ENDIF

C***    Summary output location
        IF(NTCOQU(noco).GT.0) THEN
          SUMMARYFILE=COQU(noco,1)
          OPFILE=.TRUE.
          IOFI=IOFILE3
          CALL STRING_TRIM(SUMMARYFILE,IBEG1,IEND1)
          CALL OPENF(IOFI,'DISK',SUMMARYFILE(IBEG1:IEND1)//'.opnois',
     '      'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

C***    Set default output extension  (remove this ??)
        IF(DEFNOIS_TYPE.EQ.1) THEN !absolute noise
CC JMBs 30-MAY-2000
C          IF(NOIS_LEV.LT.10) THEN !one s.f.
C            WRITE(LEV,'(I1)') NOIS_LEV
C          ELSE
C            WRITE(LEV,'(I2)') NOIS_LEV
C          ENDIF
          WRITE(LEV,'(D12.4)') NOIS_LEV
CC JMBe
        ENDIF


        IF(APPLYTODATA) THEN !apply noise to data

          nr=1
          SEED1=1 !initialise seeds
          SEED2=1
          SEED3=1

          RMS(1)=0.D0
          DO nd=1,NDDATA(0,nr)
            DO nj=1,NJT
              IF(NOIS_DIST.EQ.1)THEN
                TMP=RAN10(seed1,seed2,seed3)
                RMS(1)=RMS(1)+TMP*TMP
                ZD(nj,nd)=ZD(nj,nd)+TMP*NOIS_DISP
              ELSE IF(NOIS_DIST.EQ.2) THEN
                TMP=RAN0(seed1)
                RMS(1)=RMS(1)+TMP*TMP
                ZD(nj,nd)=ZD(nj,nd)+TMP*NOIS_DISP*2.D0
              ELSE IF(NOIS_DIST.EQ.3) THEN
                TMP=RAN1(seed1)
                RMS(1)=RMS(1)+TMP*TMP
                ZD(nj,nd)=ZD(nj,nd)+TMP*NOIS_DISP*2.D0
              ELSE
                CALL ASSERT(NOIS_DIST.LE.3,'>>Update Code',ERROR,*9999)
              ENDIF
            ENDDO
            RMS(1)=SQRT(RMS(1)/NDDATA(0,nr)/NJT)
          ENDDO
          WRITE(OP_STRING(1),'('' Summary of Error Applied to Data'')')
          WRITE(OP_STRING(2),
     '      '(''   Data points         = '',I12)')   NDDATA(0,nr)
          WRITE(OP_STRING(3),
     '      '(''   Total Components    = '',I12)')   NDDATA(0,nr)*NJT
          WRITE(OP_STRING(4),
     '      '(''   Expected RMS        = '',D12.4)') NOIS_DISP
          WRITE(OP_STRING(5),
     '      '(''   Actual  RMS         = '',D12.4)') RMS(1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE !apply noise to signal

C*** Set infile
          IF(CBBREV(CO,'SIGNAL',3,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            INFILE=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            INFILE=FILE00(IBEG:IEND)
          ENDIF

C*** Set outfile  (remove this ??)
          IF(CBBREV(CO,'OUTFILE',3,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            OUTFILE=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            OUTFILE=FILE00(IBEG:IEND)//'_ns'//LEV
          ENDIF

C*** Set Filetype
          IF(CBBREV(CO,'BINARY',3,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF

C*** Set starting seed format
          IF(CBBREV(CO,'SET_RANDOM_START',3,noco+1,NTCO,N3CO)) THEN
            RANDOM_START=.TRUE.
          ELSE
            RANDOM_START=.FALSE.
          ENDIF

C LKC 8-JAN-98 Initialise time to some value
          TIME=0.D0

C***    Open input signal file
          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '      'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C***    Read the electrode data
          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '      'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C***    IF ASCII need to calculate NUMTIMEDATA
          IF(FILEFORMAT.EQ.'ASCII') THEN
            NUMTIMEDATA=0
            ENDFILE=.FALSE.
            DO WHILE(.NOT.ENDFILE)
              CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
            ENDDO
          ENDIF

          SIGNAL_HOWGENERATED(IOFILE1)=SIGNAL_HOWGENERATED(IOFILE2)
          SIGNAL_HEADER(IOFILE1)=SIGNAL_HEADER(IOFILE2)
          SIGNAL_NUMREGIONS(IOFILE1)=SIGNAL_NUMREGIONS(IOFILE2)
          SIGNAL_ELEMLOC(IOFILE1)=SIGNAL_ELEMLOC(IOFILE2)
          DO nr=1,SIGNAL_NUMREGIONS(IOFILE2)
            SIGNAL_NUMELEC(nr,IOFILE1)=SIGNAL_NUMELEC(nr,IOFILE2)
            SIGNAL_REGNAME(nr,IOFILE1)=SIGNAL_REGNAME(nr,IOFILE2)
            SIGNAL_REGTYPE(nr,IOFILE1)=SIGNAL_REGTYPE(nr,IOFILE2)
            IF(SIGNAL_ELEMLOC(IOFILE2).GT.0) THEN
              SIGNAL_NUMXI(nr,IOFILE1)=SIGNAL_NUMXI(nr,IOFILE2)
            ENDIF
          ENDDO !nr

C*** The starting seeds by default are initialised to 1.
C*** Subsequent calls to the routine can enable the starting
C*** seeds as they were last set.
          IF(.NOT.RANDOM_START.OR.FIRST) THEN
            SEED1=1 !initialise seeds
            SEED2=1
            SEED3=1
            FIRST=.FALSE.
          ENDIF
          INIT_SEED1=SEED1
          INIT_SEED2=SEED2
          INIT_SEED3=SEED3

C*** Noise for the electrode locations
          IF(DEFNOIS_LOC.EQ.2.OR.DEFNOIS_LOC.EQ.3) THEN

            OP_STRING(1)='Warning : Need to update projections'
            OP_STRING(2)='Warning : Need to update xi positions'
            CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C*** Reset the input file
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '        'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C*** Open the output file
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '        'OPEN',ENDFILE,.TRUE.,ERROR,*9999)


            OP_STRING(1)=' '
            OP_STRING(2)=' Summary of Error Applied to Electodes'
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

            DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
              NRAND=0 !initialise
              MEAN(1)=0.D0
              MEAN(2)=0.D0
              MEAN(3)=0.D0
              RMS(1)=0.D0

              DO nodata=1,SIGNAL_NUMELEC(nr,IOFILE2)
                nd=NDDATA(nodata,nr)
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  IF(NOIS_DIST.EQ.1)THEN
                    ELECT_ERR(nj,1)=RAN10(seed1,seed2,seed3)
     '              *DBLE(NOIS_DISP)
                  ELSE IF(NOIS_DIST.EQ.2) THEN
                    ELECT_ERR(nj,1)=RAN0(seed1)*DBLE(NOIS_DISP)
                  ELSE IF(NOIS_DIST.EQ.3) THEN
                    ELECT_ERR(nj,1)=RAN1(seed1)*DBLE(NOIS_DISP)
                  ELSE
                    CALL ASSERT(NOIS_DIST.LE.3,'>>Update Code',
     '                ERROR,*9999)
C                   ! should not get here unless ipnois changed
                  ENDIF

                  ZD(nj,nd)=ZD(nj,nd)+ELECT_ERR(nj,1)
                  NRAND=NRAND+1
                  MEAN(nj)=MEAN(nj)+DABS(ELECT_ERR(nj,1))
                  RMS(1)=RMS(1)+ELECT_ERR(nj,1)**2
                ENDDO !nj
              ENDDO !nd
              RMS(1)=DSQRT(RMS(1)/DBLE(NRAND))


C*** Output summary
              WRITE(OP_STRING(1),
     '          '(''   Region              = '',I12)') nr

              WRITE(OP_STRING(2),
     '          '(''   Number of samples   = '',I12)') NRAND
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                WRITE(OP_STRING(1),
     '            '(''    In direction xj('',I1,'')'')') nj
                WRITE(OP_STRING(2),
     '            '(''      Mean Absolute    = '',D12.4)')
     '            MEAN(nj)/NRAND
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !nj

              WRITE(OP_STRING(1),
     '        '(''   Expected RMS        = '',D12.4)') DBLE(NOIS_DISP)
              WRITE(OP_STRING(2),
     '          '(''   Actual  RMS         = '',D12.4)') RMS(1)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !nr

C*** Write the electrode data
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '        'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

            DO notime=1,NUMTIMEDATA !Loop over the times
C*** Read and write the signal data
              CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA1,0,
     '          SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'WRITE',
     '          FILEFORMAT,OUTFILE,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            ENDDO !notime

C*** Close output signals
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,OUTFILE,
     '        ' ',ENDFILE,.TRUE.,ERROR1,*9998)

          ENDIF !applying electrode noise

          IF(DEFNOIS_LOC.EQ.1.OR.DEFNOIS_LOC.EQ.3) THEN

C*** Reset the input file
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '        'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C*** Open the output file
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '        'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C*** Write the electrode data
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '        'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Noise for the signal data
            SEED1=INIT_SEED1 !initialise seeds
            SEED2=INIT_SEED2
            SEED3=INIT_SEED3
            NRAND=0 !initialise
            MEAN(1)=0.D0
            RMS(1)=0.D0
            DO notime=1,NUMTIMEDATA !Loop over the times
C***        Read the signal data in
              CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Add the noise to the clean signal
C*** Additional factor of 2 for the normal distributions
              DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
                DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE1)
                  nd=NDDATA(noelec,nr)
                  IF(NOIS_DIST.EQ.1)THEN
                    ZD(NJT+2,nd)=RAN10(seed1,seed2,seed3)*NOIS_LEV
                  ELSE IF(NOIS_DIST.EQ.2) THEN
                    ZD(NJT+2,nd)=RAN0(seed1)*NOIS_LEV*2.D0
                  ELSE IF(NOIS_DIST.EQ.3) THEN
                    ZD(NJT+2,nd)=RAN1(seed1)*NOIS_LEV*2.D0
                  ELSE
                    CALL ASSERT(NOIS_DIST.LE.3,'>>Update Code',
     '                ERROR,*9999)
C                 ! should not get here unless ipnois changed
                  ENDIF
                  MEAN(1)=MEAN(1)+ZD(NJT+2,nd) !sum noise
                  ZD(NJT+1,nd)=ZD(NJT+1,nd)+ZD(NJT+2,nd) !dirty orig sig
                  RMS(1)=RMS(1)+ZD(NJT+2,nd)*ZD(NJT+2,nd) !sqr noise
                  NRAND=NRAND+1
                ENDDO !noelec
              ENDDO !nr
C*** Write the new noisy signal
              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA1,0,
     '          SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'WRITE',
     '          FILEFORMAT,OUTFILE,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            ENDDO !notime
            MEAN(1)=MEAN(1)/DBLE(NRAND)
            RMS(1)=SQRT(RMS(1)/DBLE(NRAND))

C*** Recalculating the RN's because they are not stored
C*** for each time period. This is required for the calculation
C*** of the SD as the mean of the noise is not known a priori
            SEED1=INIT_SEED1 !reset the seeds
            SEED2=INIT_SEED2
            SEED3=INIT_SEED3

C***      Standard Deviation Calculation
            SD=0.D0
C         !Loops separately for each noise type to prevent multiple ifs
            IF(NOIS_DIST.EQ.1)THEN
              DO notime=1,NUMTIMEDATA
                DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
                  DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE1)
                    nd=NDDATA(noelec,nr)
                    SDTEMP=RAN10(seed1,seed2,seed3)*NOIS_LEV-MEAN(1)
                    SD=SD+SDTEMP*SDTEMP
                  ENDDO
                ENDDO
              ENDDO !notime
            ELSE IF(NOIS_DIST.EQ.2) THEN
              DO notime=1,NUMTIMEDATA
                DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
                  DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE1)
                    nd=NDDATA(noelec,nr)
                    SDTEMP=RAN0(seed1)*NOIS_LEV-MEAN(1)
                    SD=SD+SDTEMP*SDTEMP
                  ENDDO
                ENDDO
              ENDDO !notime
            ELSE IF(NOIS_DIST.EQ.3) THEN
              DO notime=1,NUMTIMEDATA
                DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
                  DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE1)
                    nd=NDDATA(noelec,nr)
                    SDTEMP=RAN1(seed1)*NOIS_LEV-MEAN(1)
                    SD=SD+SDTEMP*SDTEMP
                  ENDDO
                ENDDO
              ENDDO !notime
            ELSE
              CALL ASSERT(NOIS_DIST.LE.3,'>>Unknown noise type',
     '          ERROR,*9999)
            ENDIF
            SD=DSQRT(1.D0/DBLE(NRAND-1)*SD)

C***      Write summary information for signal noise
            OP_STRING(1)=' '
            OP_STRING(2)='Summary of Noise Applied to Signal'
            WRITE(OP_STRING(3),'(''  Noise Type     ='',I12)')
     '        DEFNOIS_TYPE
            WRITE(OP_STRING(4),'(''  Noise Samples  ='',I12)') NRAND
            WRITE(OP_STRING(5),'(''  Mean Noise     ='',D12.4)') MEAN(1)
            WRITE(OP_STRING(6),'(''  Applied RMS    ='',D12.4)')
     '        NOIS_LEV
            WRITE(OP_STRING(7),'(''  Actual RMS     ='',D12.4)') RMS(1)
            WRITE(OP_STRING(8),'(''  Std. Deviation ='',D12.4)') SD
            OP_STRING(9)='  '
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C*** Close output signal
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,OUTFILE,
     '        ' ',ENDFILE,.TRUE.,ERROR1,*9998)
          ENDIF !DEFNOIS_LOC

C*** Close input signal
          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE,
     '      ' ',ENDFILE,.TRUE.,ERROR1,*9998)

        ENDIF !apply noise to signal/data

C***    Close the files
        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF

      ENDIF

      CALL EXITS('APNOIS')
      RETURN
 9999 CALL ERRORS('APNOIS',ERROR)
C***   Close the Files if and error has occurs
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,OUTFILE,
     '  ' ',ENDFILE,.TRUE.,ERROR1,*9998)
 9998 CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE,
     '  ' ',ENDFILE,.TRUE.,ERROR1,*9997)
 9997 CALL EXITS('APNOIS')
      RETURN 1
      END


