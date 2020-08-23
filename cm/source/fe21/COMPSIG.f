      SUBROUTINE COMPSIG(LD,LIST,NBJ,NDDATA,NRLIST,WD,XID,ZD,ZD2,
     '  STRING,ERROR,*)

C#### Subroutine: COMPSIG
C###  Description:
C###    <HTML>
C###    COMPSIG compares signal files with respect to some measures.
C###    The signals can be compared at a given time instance, a
C###    specified period or the entire signal. Data samples can be
C###    skipped on the MASTER file with the use of SKIPSAMPLE.
C###
C###    The following types of measures are used to compare the
C###    signals.
C###    <OL>
C###      <LI> RMS
C###      <LI> relative RMS
C###      <LI> similiarity index (SI)
C###      <LI> correlation coefficient (CC)
C###      <LI> max/min signal value
C###      <LI> max/min differences between signals (ABS,%)
C###      <LI> max/min signal location
C###      <LI> amount of movement of max/min signal location
C###    </OL>
C###    </HTML>

C***  LKC 13-NOV-97 Added Max/Min comparisons
C***  LKC 23-NOV-97 Added Integral over a region
C***  LKC 19-JUN-99 Added sample skipping

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'sign00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER LD(NDM),LIST(0:NLISTM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     '  NRLIST(0:NRM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM),ZD2(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Parameters
      INTEGER NEXCEPTM
      PARAMETER(NEXCEPTM=300)

!     Local Variables
      INTEGER COUNT,IBEG,IBEG1,IEND,IEND1,N3CO,nd,nexcept,nj,nodata,
     '  nonr,nr,NUMTIMEDATA,NUMTIMEDATA1,notime,LIST2(0:NEXCEPTM),
     '  NSAMPLE,nsig,noskip,
     '  maxsig_nd(2),minsig_nd(2),SKIPSAMPLE
      REAL*8 COMPARE,COMPARESUM,COMPARESQUAREDSUM,DIFF,DIFFSQUAREDSUM,
     '  END_TIME,MASTER,MASTERSUM,MASTERSQUAREDSUM,
     '  maxposn_diff,maxsig_diff,maxsig_per,maxsig_total(2),
     '  minposn_diff,minsig_diff,minsig_per,minsig_total(2),
     '  PRODUCTSUM,RELRMS,RFROMC,
     '  RMS,SI,SIGNALMAX(9),SIGNALMIN(9),START_TIME,TIME,TMP
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,ENDFILE,ENTIRE_SIGNAL,
     '  GEOMETRY,INTEGRAL,OPFILE
      CHARACTER BASIS*9,COMPAREFILE*(MXCH),ERROR_DUMMY*255,FILEFORMAT*6,
     '  MASTERFILE*(MXCH),OUTPUTFILE*(MXCH)
!     Function
      INTEGER IFROMC


      CALL ENTERS('COMPSIG',*9999)

      CALL STRING_TRIM(FILE00,IBEG1,IEND1)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

c cpb 29/8/00 Adding electrode and time bases.

C---------------------------------------------------------------------

C#### Command: FEM compare signal;FILENAME
C###  Description:
C###    Compares the master signal file with the compare signal file
C###    according to some measure. The comparison can occur over
C###    electrodes (electrode basis) at a time instance (TIME), the
C###    entire signal (INTEGRAL), or a specified time period
C###    (START_TIME # END_TIME #), or the comparison can occur over
C###    time (time basis) for a given list of electrodes, or a single
C###    comparision can over over both time and electordes (both basis).
C###  Parameter:      <basis (electrode/time/both)[electrodes]>
C###    Specify wether the comparison basis is over electrode(s) for a
C###    given time (or times), or whether the comparison basis is over
C###    time for a given electrode(s), or whether the comparision is
C###    over both time and electrodes.
C###  Parameter:      <masterfile FILENAME[$current]>
C###    The main signal file to compare to.
C###  Parameter:      <comparefile FILENAME[$current]>
C###    The comparison signal file
C###  Parameter:      <electrodes (#s/all)[all]>
C###    List of electrode to compare.
C###  Parameter:      <except (#s/none)[none]>
C###    If all electrodes are used then some can be excluded using
C###    this tag.
C###  Parameter:      <(time #/integral/start_time # end_time #)[time 0.0]>
C###    The type of comparison to use for time based comparisons. For
C###    electrode based comparisons the start_time and end_time can be
C###    used to specify the time for the comparison. If omitted
C###    they default to the entire signal.
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Whether the files are ascii (.ipsign) or binary (.binsig)
C###  Parameter:      <skipsample #[1]>
C###    Whether to skip time samples in the master file for an integral
C###    comparison or a time based comparison.
C###  Parameter:      <region (#s/all)[1]>
C###    Which region to compare. Not valid at the moment.


        OP_STRING(1)=STRING(1:IEND)//';FILENAME'
        OP_STRING(2)=BLANK(1:15)//
     '    '<basis (electrode/time/both)[electrode]>'
        OP_STRING(3)=BLANK(1:15)//
     '    '<masterfile FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(4)=BLANK(1:15)//
     '    '<comparefile FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(5)=BLANK(1:15)//'<electrodes (#s/all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<except (#s/none)[none]>'
        OP_STRING(7)=BLANK(1:15)//
     '    '<(time #/integral/start_time # end_time #)[time 0.0]>'
        OP_STRING(8)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(9)=BLANK(1:15)//'<skipsample #[1]>'
        OP_STRING(10)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM compare signal;FILENAME geometry
C###  Description:
C###    Compares the electrode posistions of the  master signal
C###    file with the compare signal file
C###  Parameter:      <masterfile FILENAME[$current]>
C###    The main signal file to compare to.
C###  Parameter:      <comparefile FILENAME[$current]>
C###    The comparison signal file
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Whether the files are ascii (.ipsign) or binary (.binsig)
C###  Parameter:      <region (#s/all)[1]>
C###    Which region to compare. Not valid at the moment.


        OP_STRING(1)=STRING(1:IEND)//';FILENAME geometry'
        OP_STRING(2)=BLANK(1:15)//
     '    '<masterfile FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//
     '    '<comparefile FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(4)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(5)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------



      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','COMPSIG',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          OUTPUTFILE=COQU(noco,1)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL STRING_TRIM(OUTPUTFILE,IBEG1,IEND1)
          CALL OPENF(IOFI,'DISK',OUTPUTFILE(IBEG1:IEND1)//'.opcomp',
     '      'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)


        IF(CBBREV(CO,'BASIS',1,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'ELECTRODE',1)) THEN
            BASIS='ELECTRODE'
          ELSE IF(ABBREV(CO(N3CO+1),'TIME',1)) THEN
            BASIS='TIME'
CC JMB 13-SEP-2000 Adding combined electrode + time basis
          ELSE IF(ABBREV(CO(N3CO+1),'BOTH',1)) THEN
            BASIS='BOTH'
          ELSE
            ERROR='>>Invalid basis for comparison'
            GOTO 9999
          ENDIF
        ELSE
          BASIS='ELECTRODE'
        ENDIF

        GEOMETRY=.FALSE.
        IF(CBBREV(CO,'GEOMETRY',3,noco+1,NTCO,N3CO)) THEN
          GEOMETRY=.TRUE.
          BASIS='NEITHER'
        ENDIF


        IF(CBBREV(CO,'MASTERFILE',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          MASTERFILE=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          MASTERFILE=FILE00(IBEG1:IEND1)
        ENDIF
        IF(CBBREV(CO,'COMPAREFILE',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          COMPAREFILE=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          COMPAREFILE=FILE00(IBEG1:IEND1)
        ENDIF


C cpb 19/11/97 Adding integral
C LKC 23-NOV-97 Adding integral for a time period

        IF(BASIS(1:9).EQ.'ELECTRODE') THEN
          IF(CBBREV(CO,'TIME',2,noco+1,NTCO,N3CO)) THEN
            INTEGRAL=.FALSE.
            ENTIRE_SIGNAL=.FALSE.
            TIME=RFROMC(CO(N3CO+1))
          ELSE IF(CBBREV(CO,'INTEGRAL',2,noco+1,NTCO,N3CO)) THEN
            INTEGRAL=.TRUE.
            ENTIRE_SIGNAL=.TRUE.
            TIME=0.0d0
          ELSE IF(CBBREV(CO,'START_TIME',2,noco+1,NTCO,N3CO)) THEN
C**** Integrates a given time inteval
            INTEGRAL=.TRUE.
            ENTIRE_SIGNAL=.FALSE.
            IF(N3CO.LT.NTCO) THEN
              START_TIME=RFROMC(CO(N3CO+1))
              IF(CBBREV(CO,'END_TIME',2,noco+1,NTCO,N3CO)) THEN
                IF(N3CO.LT.NTCO) THEN
                  END_TIME=RFROMC(CO(N3CO+1))
                ELSE
                  ERROR='>>No END_TIME supplied'
                  GOTO 9999
                ENDIF
              ELSE
                ERROR='>>No END_TIME supplied'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>No START_TIME supplied'
              GOTO 9999
            ENDIF
          ELSE
            INTEGRAL=.FALSE. !default to TIME=0.0
            ENTIRE_SIGNAL=.FALSE.
            TIME=0.0d0
          ENDIF

          IF(CBBREV(CO,'SKIPSAMPLE',2,noco+1,NTCO,N3CO)) THEN
            IF(INTEGRAL) THEN
              SKIPSAMPLE=IFROMC(CO(N3CO+1))
            ELSE
              OP_STRING(1)=
     '          'WARNING: SKIPSAMPLE only valid to integral comparison'
              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
              SKIPSAMPLE=1
            ENDIF
          ELSE
            SKIPSAMPLE=1
          ENDIF
        ELSE IF(BASIS(1:4).EQ.'TIME') THEN
          IF(CBBREV(CO,'START_TIME',2,noco+1,NTCO,N3CO)) THEN
C**** Integrates a given time inteval
            ENTIRE_SIGNAL=.FALSE.
            IF(N3CO.LT.NTCO) THEN
              START_TIME=RFROMC(CO(N3CO+1))
              IF(CBBREV(CO,'END_TIME',2,noco+1,NTCO,N3CO)) THEN
                IF(N3CO.LT.NTCO) THEN
                  END_TIME=RFROMC(CO(N3CO+1))
                ELSE
                  ERROR='>>No END_TIME supplied'
                  GOTO 9999
                ENDIF
              ELSE
                ERROR='>>No END_TIME supplied'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>No START_TIME supplied'
              GOTO 9999
            ENDIF
          ELSE
            ENTIRE_SIGNAL=.TRUE.
          ENDIF
          IF(CBBREV(CO,'SKIPSAMPLE',2,noco+1,NTCO,N3CO)) THEN
            SKIPSAMPLE=IFROMC(CO(N3CO+1))
          ELSE
            SKIPSAMPLE=1
          ENDIF
        ENDIF

        IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF

C*** Open up the two signal files
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFILE,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFILE,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C*** Read in the electrode data
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFILE,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C LKC 25-SEP-2000 Change ZD to ZD2 to electrode positions can
C       be compared
C        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
C     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFILE,
C     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD2,'READ',FILEFORMAT,COMPAREFILE,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)



C*** Electrode position comparison ...

        IF(GEOMETRY) THEN

          WRITE(*,*)
          WRITE(*,*) 'comparing electrode geometry'
          WRITE(*,*)

          MAXSIG_ND(1)=0
          MINSIG_ND(1)=0

          MAXSIG_DIFF=-RMAX !initialise max and min displacements
          MINSIG_DIFF=RMAX
          COUNT=0
          RMS=0.D0
          DO nd=1,NDT
            COUNT=COUNT+1 !currently COUNT should == NDT
            DIFFSQUAREDSUM=0.D0
            DO nj=1,NJT
              DIFFSQUAREDSUM=DIFFSQUAREDSUM+(ZD(nj,nd)-ZD2(nj,nd))**2
            ENDDO !nj
            TMP=SQRT(DIFFSQUAREDSUM)

            IF(TMP.LT.MINSIG_DIFF) THEN ! set max and min signals
              MINSIG_DIFF=TMP
              MINSIG_ND(1)=COUNT
            ENDIF
            IF(TMP.GT.MAXSIG_DIFF) THEN
              MAXSIG_DIFF=TMP
              MAXSIG_ND(1)=COUNT
            ENDIF

            RMS=RMS+DIFFSQUAREDSUM
          ENDDO !nd
          RMS=SQRT(RMS/DBLE(COUNT))




          WRITE(OP_STRING(1),'('' '')')
          WRITE(OP_STRING(2),'('' Geometric data comparison'')')
          WRITE(OP_STRING(3),'('' '')')
          WRITE(OP_STRING(4),'('' RMS displacement     :'',F12.3)')
     '      RMS
          WRITE(OP_STRING(5),'('' Maximum displacement :'',F12.3)')
     '      MAXSIG_DIFF
          WRITE(OP_STRING(6),'(''   at Electrode       :'',I12)')
     '      MAXSIG_ND(1)
          WRITE(OP_STRING(7),'('' Minimum displacement :'',F12.3)')
     '      MINSIG_DIFF
          WRITE(OP_STRING(8),'(''   at Electrode       :'',I12)')
     '      MINSIG_ND(1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)


        ENDIF



        IF(BASIS(1:9).EQ.'ELECTRODE') THEN
C*** Check the number of electodes are the same
          DO nonr=1,NRLIST(0)
            nr=NRLIST(nonr)
            CALL ASSERT(
     '        SIGNAL_NUMELEC(nr,IOFILE2).EQ.SIGNAL_NUMELEC(nr,IOFILE1),
     '        '>>Different number of electrodes',ERROR,*9999)
          ENDDO !nonr
        ENDIF

C***    IF ASCII need to calculate NUMTIMEDATA
        IF(FILEFORMAT(1:5).EQ.'ASCII'.AND.(INTEGRAL.OR.
     '    BASIS.EQ.'TIME'.OR.BASIS.EQ.'BOTH')) THEN
          NUMTIMEDATA=0
          NUMTIMEDATA1=0
          ENDFILE=.FALSE.
          DO WHILE(.NOT.ENDFILE)
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
          ENDDO
          ENDFILE=.FALSE.
          DO WHILE(.NOT.ENDFILE)
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            IF(.NOT.ENDFILE) NUMTIMEDATA1=NUMTIMEDATA1+1
          ENDDO
        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN
          CALL ASSERT(NUMTIMEDATA.GE.1,
     '      '>> No signal in master file',ERROR,*9999)
          CALL ASSERT(NUMTIMEDATA1.GE.1,
     '      '>> No signal in comparison file',ERROR,*9999)
        ENDIF

        nr=NRLIST(1)
        nexcept=1
        nd=1
        IF(CBBREV(CO,'ELECTRODES',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'ALL',1)) THEN

            IF (CBBREV(CO,'EXCEPT',3,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),NEXCEPTM,LIST2(0),LIST2(1),
     '          ERROR,*9999)

              LIST(0)=NDDATA(0,nr)
              DO nodata=1,LIST(0)
                IF(LIST2(nexcept).NE.NDDATA(nodata,nr)) THEN
                  LIST(nd)=NDDATA(nodata,nr)
                  nd=nd+1
                ELSE
                  nexcept=nexcept+1
                ENDIF
              ENDDO
              LIST(0)=NDDATA(0,nr)-LIST2(0)
            ELSE
              LIST(0)=NDDATA(0,nr)
              DO nodata=1,LIST(0)
                LIST(nodata)=NDDATA(nodata,nr)
              ENDDO
            ENDIF
          ELSE !specified electrodes
            CALL PARSIL(CO(N3CO+1),NLISTM,LIST(0),LIST(1),ERROR,*9999)

            DO nodata=1,LIST(0)
              nd=LIST(nodata)
              IF(nd.GT.0.AND.nd.LE.NDDATA(0,nr)) THEN
                LIST(nodata)=NDDATA(nd,nr)
              ELSE
                ERROR='>>Invalid electrode number '
                GOTO 9999
              ENDIF
            ENDDO
          ENDIF
        ELSE !default ALL Electrodes
          IF (CBBREV(CO,'EXCEPT',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NEXCEPTM,LIST2(0),LIST2(1),
     '        ERROR,*9999)

            LIST(0)=NDDATA(0,nr)
            DO nodata=1,LIST(0)
              IF(LIST2(nexcept).NE.NDDATA(nodata,nr)) THEN
                LIST(nd)=NDDATA(nodata,nr)
                nd=nd+1
              ELSE
                nexcept=nexcept+1
              ENDIF
            ENDDO
            LIST(0)=NDDATA(0,nr)-LIST2(0)
          ELSE
            LIST(0)=NDDATA(0,nr)
            DO nodata=1,LIST(0)
              LIST(nodata)=NDDATA(nodata,nr)
            ENDDO
          ENDIF !EXCEPT
        ENDIF !ELECTRODES

        NSAMPLE=0

        IF(BASIS(1:9).EQ.'ELECTRODE') THEN

C*** Reset both files
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFILE,
     '      'RESET',ENDFILE,.TRUE.,ERROR,*9999)
          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFILE,
     '      'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C*** Get the two sets of signal data at the correct time or add up
C*** all the times (integral).

C*** Compare the entire signal range
          IF(INTEGRAL.AND.ENTIRE_SIGNAL) THEN
            CALL ASSERT(NJT+3.LE.NJM,'>>Increase NJM',ERROR,*9999)

            IF(NUMTIMEDATA.NE.NUMTIMEDATA1) THEN
              WRITE(OP_STRING(1),'('' Master length      = '',I8)')
     '          NUMTIMEDATA
              WRITE(OP_STRING(2),'('' Comparision length = '',I8)')
     '          NUMTIMEDATA1
              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
              ERROR='>> # Time samples not equal for integral'
              GOTO 9999
            ENDIF


            NSAMPLE=NUMTIMEDATA
            DO nodata=1,LIST(0)
              nd=LIST(nodata)
              ZD(NJT+1,nd)=0.0d0
              ZD(NJT+2,nd)=0.0d0
            ENDDO !nodata

            DO notime=1,NUMTIMEDATA

              DO noskip=1,SKIPSAMPLE
                CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,2,
     '            SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '            FILEFORMAT,MASTERFILE,
     '            'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              ENDDO
              DO nodata=1,LIST(0)
                nd=LIST(nodata)
                ZD(NJT+1,nd)=ZD(NJT+1,nd)+ZD(NJT+3,nd)
              ENDDO !nodata

              CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,2,
     '          SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '          FILEFORMAT,COMPAREFILE,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              DO nodata=1,LIST(0)
                nd=LIST(nodata)
                ZD(NJT+2,nd)=ZD(NJT+2,nd)+ZD(NJT+3,nd)
              ENDDO !noskip
            ENDDO !notime

C*** Compare a specified time period
          ELSE IF(INTEGRAL.AND..NOT.ENTIRE_SIGNAL) THEN
            CALL ASSERT(NJT+3.LE.NJM,'>>Increase NJM',ERROR,*9999)
            DO nodata=1,LIST(0)
              nd=LIST(nodata)
              ZD(NJT+1,nd)=0.0d0
              ZD(NJT+2,nd)=0.0d0
            ENDDO !nodata

            NUMTIMEDATA=0
            ENDFILE=.FALSE.
            TIME=START_TIME !  Get to the correct time
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,2,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '        FILEFORMAT,MASTERFILE,
     '        'SIGNAL_DATA',ENDFILE,.FALSE.,ERROR,*9999)

            IF(SKIPSAMPLE.EQ.1) THEN !else don't do 1st one
              DO nodata=1,LIST(0)
                nd=LIST(nodata)
                ZD(NJT+1,nd)=ZD(NJT+1,nd)+ZD(NJT+3,nd)
              ENDDO !nodata
            ENDIF
            CALL ASSERT(.NOT.ENDFILE,
     '        '>>Could not find specified time in MASTER file '
     '        ,ERROR,*9999)

            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,2,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '        FILEFORMAT,COMPAREFILE,
     '        'SIGNAL_DATA',ENDFILE,.FALSE.,ERROR,*9999)
            IF(SKIPSAMPLE.EQ.1) THEN !else don't do 1st one
              DO nodata=1,LIST(0)
                nd=LIST(nodata)
                ZD(NJT+2,nd)=ZD(NJT+2,nd)+ZD(NJT+3,nd)
              ENDDO !nodata
            ENDIF

            CALL ASSERT(.NOT.ENDFILE,
     '        '>>Could not find specified time in COMPARE file '
     '        ,ERROR,*9999)
            NUMTIMEDATA=NUMTIMEDATA+1

C           Continue until END_TIME
            DO WHILE((.NOT.ENDFILE).AND.(TIME.LT.END_TIME))
              DO noskip=1,SKIPSAMPLE
                CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,2,
     '            SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '            FILEFORMAT,MASTERFILE,'SIGNAL_DATA',ENDFILE,
     '            .TRUE.,ERROR,*9999)
              ENDDO
              NSAMPLE=NSAMPLE+1
              DO nodata=1,LIST(0)
                nd=LIST(nodata)
                ZD(NJT+1,nd)=ZD(NJT+1,nd)+ZD(NJT+3,nd)
              ENDDO !nodata
              CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,2,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFILE,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              DO nodata=1,LIST(0)
                nd=LIST(nodata)
                ZD(NJT+2,nd)=ZD(NJT+2,nd)+ZD(NJT+3,nd)
              ENDDO !nodata
            ENDDO ! WHILE not ENDFILE && time < ENDTIME

C*** set the actual end_time
            END_TIME=TIME
            IF(ENDFILE) THEN
              WRITE(OP_STRING,
     '          '('' >>WARNING: EOF occurred before END_TIME'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

C*** Compare at a given time instance
          ELSE
            CALL ASSERT(NJT+2.LE.NJM,'>>Increase NJM',ERROR,*9999)
            NSAMPLE=1
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFILE,
     '        'SIGNAL_DATA',ENDFILE,.FALSE.,ERROR,*9999)
            CALL ASSERT(.NOT.ENDFILE,
     '        '>>Could not findtime in MASTER file ',ERROR,*9999)
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,1,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFILE,
     '        'SIGNAL_DATA',ENDFILE,.FALSE.,ERROR,*9999)
            CALL ASSERT(.NOT.ENDFILE,
     '        '>>Could not find time in COMPARE file ',ERROR,*9999)
          ENDIF

C*** Calculate RMS, rel RMS, similarity index
          MASTERSUM=0.0d0
          MASTERSQUAREDSUM=0.0d0
          COMPARESUM=0.0d0
          COMPARESQUAREDSUM=0.0d0
          PRODUCTSUM=0.0d0
          DIFFSQUAREDSUM=0.0d0
          COUNT=0
          DO nodata=1,LIST(0)
            nd=LIST(nodata)
            MASTER=ZD(NJT+1,nd)
            COMPARE=ZD(NJT+2,nd)
            MASTERSUM=MASTERSUM+MASTER
            COMPARESUM=COMPARESUM+COMPARE
            MASTERSQUAREDSUM=MASTERSQUAREDSUM+MASTER**2
            COMPARESQUAREDSUM=COMPARESQUAREDSUM+COMPARE**2
            PRODUCTSUM=PRODUCTSUM+MASTER*COMPARE
            DIFFSQUAREDSUM=DIFFSQUAREDSUM+(MASTER-COMPARE)**2
            COUNT=COUNT+1
          ENDDO !nd


C LKC 18-OCT-2000 adding check for no data to prevent divide by 0
C          CALL ASSERT(DABS(MASTERSQUAREDSUM).GT.ZERO_TOL,
C     '      '>> No data in Master signal file',ERROR,*9999)
C GBS 25-FEB-2002 changed not to exit if sum is zero.
          RMS=DSQRT(DIFFSQUAREDSUM/DBLE(COUNT))
          IF (DABS(MASTERSQUAREDSUM).GT.ZERO_TOL) THEN
            RELRMS=DSQRT(DIFFSQUAREDSUM/MASTERSQUAREDSUM)
          ELSE
            RELRMS=0.0d0
          ENDIF

C LKC 18-DEC-2000 adding same check for no data in comp signal
C          CALL ASSERT(DABS(COMPARESQUAREDSUM).GT.ZERO_TOL,
C     '      '>> No data in Master signal file',ERROR,*9999)

C!!! LKC 21-NOV-97
C    if the dsqrt() part.EQ.0.0 are the signals exactly the same ?
C    results in divide by 0.0

C LKC 11-MAY-1999 Can't do SI for 1 electrode

          IF(LIST(0).GT.1) THEN
            SI=(PRODUCTSUM-(MASTERSUM*COMPARESUM/DBLE(COUNT)))/
     '        DSQRT((MASTERSQUAREDSUM-MASTERSUM**2/DBLE(COUNT))*
     '        (COMPARESQUAREDSUM-COMPARESUM**2/DBLE(COUNT)))
          ELSE
            SI=0.0d0
          ENDIF

C LKC Note COUNT == LIST(0) can remove count?


C LKC 10-NOV-97 tidy up this so less WRITES statements
C        WRITE(OP_STRING,'(/'' Comparison metrics :'')')
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' At time          = '',D12.5)') TIME
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' RMS              = '',D12.5)') RMS
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Relative RMS     = '',D12.5)') RELRMS
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Similarity Index = '',D12.5/)') SI
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C*** Output comparison metrices
          OP_STRING(1)=' '
          WRITE(OP_STRING(2),'('' Comparison metrics :'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          IF(INTEGRAL) THEN
            IF(ENTIRE_SIGNAL) THEN
              WRITE(OP_STRING(1),'('' Integral          '')')
            ELSE
              WRITE(OP_STRING(1),'('' From time          = '',D12.5,'//
     '          ''' to '',D12.5)') START_TIME,END_TIME
            ENDIF
          ELSE
            WRITE(OP_STRING(1),'('' At time            = '',D12.5)')
     '        TIME
          ENDIF
          WRITE(OP_STRING(2),'('' Num. Time Samples  = '',I12)')
     '      NSAMPLE
          WRITE(OP_STRING(3),'('' Num. Electrodes    = '',I12)')
     '      LIST(0)
          WRITE(OP_STRING(4),'('' RMS                = '',D12.5)')
     '      RMS
          WRITE(OP_STRING(5),'('' Relative RMS       = '',D12.5)')
     '      RELRMS
          WRITE(OP_STRING(6),'('' Similarity Index   = '',D12.5)')
     '      SI
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C LKC 13-NOV-97 New section of comparisons
C*** Finding maximum and minimum for master THEN compare signal
          DO nsig=1,2 !master sig then compare sig
            maxsig_total(nsig)=-RMAX !initialise
            minsig_total(nsig)=RMAX
            maxsig_nd(nsig)=0
            minsig_nd(nsig)=0
            DO nodata=1,LIST(0)
              nd=LIST(nodata)
              IF(ZD(NJT+nsig,nd).GT.maxsig_total(nsig)) THEN
                maxsig_total(nsig)=ZD(NJT+nsig,nd)
                maxsig_nd(nsig)=nd
              ENDIF
              IF(ZD(NJT+nsig,nd).LT.minsig_total(nsig)) THEN
                minsig_total(nsig)=ZD(NJT+nsig,nd)
                minsig_nd(nsig)=nd
              ENDIF
            ENDDO !nd
          ENDDO !nsig

C*** Calculate differences in ranges
          maxsig_diff=DABS(maxsig_total(1)-maxsig_total(2))
          minsig_diff=DABS(minsig_total(1)-minsig_total(2))

C*** Calculate the percetange differences

C LKC 27-FEB-1999 Changes to check for divide by 0

          IF(maxsig_total(2) .GE.maxsig_total(1)) THEN
            IF(DABS(maxsig_total(1)).GT.ZERO_TOL) THEN
              maxsig_per=DABS(maxsig_diff/maxsig_total(1)*100.0d0)
            ELSE
              maxsig_per=0.0d0
            ENDIF
          ELSE
            IF(DABS(maxsig_total(1)).GT.ZERO_TOL) THEN
              maxsig_per=-DABS(maxsig_diff/maxsig_total(1)*100.0d0)
            ELSE
              maxsig_per=0.0d0
            ENDIF
          ENDIF
          IF(minsig_total(2) .GE. minsig_total(1)) THEN
            IF(DABS(minsig_total(1)).GT.ZERO_TOL) THEN
              minsig_per=DABS(minsig_diff/minsig_total(1)*100.0d0)
            ELSE
              minsig_per=0.0d0
            ENDIF
          ELSE
            IF(DABS(minsig_total(1)).GT.ZERO_TOL) THEN
              minsig_per=-DABS(minsig_diff/minsig_total(1)*100.0d0)
            ELSE
              minsig_per=0.0d0
            ENDIF
          ENDIF

C*** Calculate the movement of the max/min positions
          maxposn_diff=0.0d0
          minposn_diff=0.0d0
          DO nj=1,NJT
            maxposn_diff=maxposn_diff+(ZD(nj,maxsig_nd(1))-
     '        ZD(nj,maxsig_nd(2)))**2
            minposn_diff=minposn_diff+(ZD(nj,minsig_nd(1))-
     '        ZD(nj,minsig_nd(2)))**2
          ENDDO
          maxposn_diff=DSQRT(maxposn_diff)
          minposn_diff=DSQRT(minposn_diff)

C*** Output max/min location/values ...
          OP_STRING(1)=' '
          OP_STRING(2)=BLANK(1:19)
     '      //' Master    Comparison   ABS Diff     % Diff'
          OP_STRING(1)=' '
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C*** Output 'Maximum' information
          WRITE(OP_STRING(1),'('' MAXIMUM Signal '','
     '      //'3D12.4,2X,F8.2)')
     '      maxsig_total(1),maxsig_total(2),maxsig_diff,maxsig_per
          WRITE(OP_STRING(2),'('' Elect. Number  '',2I12)')
     '      maxsig_nd(1),maxsig_nd(2)
          WRITE(OP_STRING(3),'('' Elect. Position'','
     '      //'''                        '' ,D12.4)') maxposn_diff

          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          DO nj=1,NJT
            WRITE(OP_STRING(nj),'(''       Xj('',I1,'')    '''
     '        //',D12.4,D12.4)')
     '        nj,ZD(nj,maxsig_nd(1)),ZD(nj,maxsig_nd(2))
          ENDDO
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C Output 'Minimum' information
          OP_STRING(1)=' '
          WRITE(OP_STRING(2),
     '      '('' MINIMUM Sig.   '',3D12.4,2X,F8.2)')
     '      minsig_total(1),minsig_total(2),minsig_diff,minsig_per
          WRITE(OP_STRING(3),'('' Elect. Number  '',2I12)')
     '      minsig_nd(1),minsig_nd(2)
          WRITE(OP_STRING(4),'('' Elect. Position'''
     '      //',''                        '',D12.4)') minposn_diff
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          DO nj=1,NJT
            WRITE(OP_STRING(nj),'(''       Xj('',I1,'')    '''
     '        //',2D12.4)')
     '        nj,ZD(nj,minsig_nd(1)),ZD(nj,minsig_nd(2))
          ENDDO
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          OP_STRING(1)=' '
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C*** TIME BASIS
        ELSE IF(BASIS(1:4).EQ.'TIME') THEN

          CALL ASSERT(NJT+2.LE.NJM,'>>Increase NJM',ERROR,*9999)
          IF(ENTIRE_SIGNAL) THEN
            IF(NUMTIMEDATA.NE.NUMTIMEDATA1) THEN
              WRITE(OP_STRING(1),'('' Master length      = '',I8)')
     '          NUMTIMEDATA
              WRITE(OP_STRING(2),'('' Comparision length = '',I8)')
     '          NUMTIMEDATA1
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ERROR='>>Number time samples not equal for comparison'
              GOTO 9999
            ENDIF
            NSAMPLE=NUMTIMEDATA
          ENDIF

          DO nodata=1,LIST(0)
            nd=LIST(nodata)

C*** Reset both files
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFILE,
     '        'RESET',ENDFILE,.TRUE.,ERROR,*9999)
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFILE,
     '        'RESET',ENDFILE,.TRUE.,ERROR,*9999)

            MASTERSUM=0.0d0
            MASTERSQUAREDSUM=0.0d0
            COMPARESUM=0.0d0
            COMPARESQUAREDSUM=0.0d0
            PRODUCTSUM=0.0d0
            DIFFSQUAREDSUM=0.0d0
            COUNT=0

            IF(ENTIRE_SIGNAL) THEN
              DO notime=1,NUMTIMEDATA

                DO noskip=1,SKIPSAMPLE
                  CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,
     '              SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '              FILEFORMAT,MASTERFILE,
     '              'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
                ENDDO

                CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,1,
     '            SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '            FILEFORMAT,COMPAREFILE,
     '            'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

                MASTER=ZD(NJT+1,nd)
                COMPARE=ZD(NJT+2,nd)
                MASTERSUM=MASTERSUM+MASTER
                COMPARESUM=COMPARESUM+COMPARE
                MASTERSQUAREDSUM=MASTERSQUAREDSUM+MASTER**2
                COMPARESQUAREDSUM=COMPARESQUAREDSUM+COMPARE**2
                PRODUCTSUM=PRODUCTSUM+MASTER*COMPARE
                DIFFSQUAREDSUM=DIFFSQUAREDSUM+(MASTER-COMPARE)**2
                COUNT=COUNT+1

              ENDDO !notime

            ELSE

              NUMTIMEDATA=0
              ENDFILE=.FALSE.
              TIME=START_TIME !Get to the correct time
              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,
     '          SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '          FILEFORMAT,MASTERFILE,'SIGNAL_DATA',ENDFILE,
     '          .FALSE.,ERROR,*9999)
              CALL ASSERT(.NOT.ENDFILE,
     '          '>>Could not find specified time in master file',
     '          ERROR,*9999)
              CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,1,
     '          SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '          FILEFORMAT,COMPAREFILE,'SIGNAL_DATA',ENDFILE,
     '          .FALSE.,ERROR,*9999)
              CALL ASSERT(.NOT.ENDFILE,
     '          '>>Could not find specified time in compare file',
     '          ERROR,*9999)
              IF(SKIPSAMPLE.EQ.1) THEN !else don't do 1st one
                MASTER=ZD(NJT+1,nd)
                COMPARE=ZD(NJT+2,nd)
                MASTERSUM=MASTERSUM+MASTER
                COMPARESUM=COMPARESUM+COMPARE
                MASTERSQUAREDSUM=MASTERSQUAREDSUM+MASTER**2
                COMPARESQUAREDSUM=COMPARESQUAREDSUM+COMPARE**2
                PRODUCTSUM=PRODUCTSUM+MASTER*COMPARE
                DIFFSQUAREDSUM=DIFFSQUAREDSUM+(MASTER-COMPARE)**2
                COUNT=COUNT+1
              ENDIF

              NUMTIMEDATA=NUMTIMEDATA+1

C             Continue until END_TIME
              DO WHILE((.NOT.ENDFILE).AND.(TIME.LT.END_TIME))
                DO noskip=1,SKIPSAMPLE
                  CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,
     '              SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '              FILEFORMAT,MASTERFILE,'SIGNAL_DATA',ENDFILE,
     '              .TRUE.,ERROR,*9999)
                ENDDO
                NSAMPLE=NSAMPLE+1
                CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,1,
     '            SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,
     '            COMPAREFILE,'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
                MASTER=ZD(NJT+1,nd)
                COMPARE=ZD(NJT+2,nd)
                MASTERSUM=MASTERSUM+MASTER
                COMPARESUM=COMPARESUM+COMPARE
                MASTERSQUAREDSUM=MASTERSQUAREDSUM+MASTER**2
                COMPARESQUAREDSUM=COMPARESQUAREDSUM+COMPARE**2
                PRODUCTSUM=PRODUCTSUM+MASTER*COMPARE
                DIFFSQUAREDSUM=DIFFSQUAREDSUM+(MASTER-COMPARE)**2
                COUNT=COUNT+1
              ENDDO ! WHILE not ENDFILE && time < ENDTIME

C*** set the actual end_time
              END_TIME=TIME
              IF(ENDFILE) THEN
                WRITE(OP_STRING,
     '            '('' >>WARNING: EOF occurred before END_TIME'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF

            ENDIF

            RMS=DSQRT(DIFFSQUAREDSUM/DBLE(COUNT))
            RELRMS=DSQRT(DIFFSQUAREDSUM/MASTERSQUAREDSUM)
            IF(COUNT.GT.1) THEN
              SI=(PRODUCTSUM-(MASTERSUM*COMPARESUM/DBLE(COUNT)))/
     '          DSQRT((MASTERSQUAREDSUM-MASTERSUM**2/DBLE(COUNT))*
     '          (COMPARESQUAREDSUM-COMPARESUM**2/DBLE(COUNT)))
            ELSE
              SI=0.0d0
            ENDIF

C*** Output comparison metrices
            OP_STRING(1)=' '
            WRITE(OP_STRING(2),'('' Comparison metrics : Electrode '','
     '        //'I5)') nd
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

            IF(ENTIRE_SIGNAL) THEN
              WRITE(OP_STRING(1),'('' Entire signal'')')
            ELSE
              WRITE(OP_STRING(1),'('' From time          = '',D12.5,'//
     '          ''' to '',D12.5)') START_TIME,END_TIME
            ENDIF
            WRITE(OP_STRING(2),'('' Num. Time Samples  = '',I12)')
     '        NSAMPLE
            WRITE(OP_STRING(4),'('' RMS                = '',D12.5)')
     '        RMS
            WRITE(OP_STRING(5),'('' Relative RMS       = '',D12.5)')
     '        RELRMS
            WRITE(OP_STRING(6),'('' Similarity Index   = '',D12.5)')
     '        SI
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C KAT 25Sep00: Reset done at start of loop
CCC JMBs 13-SEP-2000 Signal files should be reset for each data set
C            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
C     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFILE,
C     '        'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,1,SIGNALMAX,
C     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFILE,
C     '        'RESET',ENDFILE,.TRUE.,ERROR,*9999)
CCC JMBe

          ENDDO !nodata

        ELSEIF(BASIS(1:4).EQ.'BOTH') THEN
C*** Calculate the Frobenius norm of the entire signal over the electrodes
          CALL ASSERT(NJT+2.LE.NJM,'>>Increase NJM',ERROR,*9999)

          IF(NUMTIMEDATA.NE.NUMTIMEDATA1) THEN
            WRITE(OP_STRING(1),'('' Master length      = '',I8)')
     '        NUMTIMEDATA
            WRITE(OP_STRING(2),'('' Comparision length = '',I8)')
     '        NUMTIMEDATA1
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ERROR='>> # Time samples not equal for integral'
            GOTO 9999
          ENDIF

C*** Reset both files
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFILE,
     '      'RESET',ENDFILE,.TRUE.,ERROR,*9999)
          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFILE,
     '      'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C*** More efficient to first loop over time then data
          MASTERSQUAREDSUM=0.0d0
          DIFFSQUAREDSUM=0.0d0

C*** Compute the max/min errors at all times steps for all electrodes
            maxsig_diff=-RMAX
            minsig_diff=RMAX

          DO notime=1,NUMTIMEDATA
C*** Masterfile
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Comparefile
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,1,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

            DO nodata=1,LIST(0)
              nd=LIST(nodata)

              MASTER=ZD(NJT+1,nd)
              COMPARE=ZD(NJT+2,nd)

              DIFF=DABS(MASTER-COMPARE)
              IF(DIFF.GT.maxsig_diff) THEN
                maxsig_diff=DIFF
C                WRITE(*,*) 'max',nd,MASTER,COMPARE
              ENDIF

              IF(DIFF.LT.minsig_diff) THEN
                minsig_diff=DIFF
C                WRITE(*,*) 'min',nd,MASTER,COMPARE
              ENDIF

              MASTERSQUAREDSUM=MASTERSQUAREDSUM+MASTER**2
              DIFFSQUAREDSUM=DIFFSQUAREDSUM+(MASTER-COMPARE)**2
            ENDDO !nodata
          ENDDO !notime
          RMS=DSQRT(DIFFSQUAREDSUM/DBLE(NUMTIMEDATA*LIST(0)))
          RELRMS=DSQRT(DIFFSQUAREDSUM/MASTERSQUAREDSUM)

          OP_STRING(1)=' '
          WRITE(OP_STRING(2),
     '      '('' Comparison metrics : Electrode + Time'')')
          WRITE(OP_STRING(3),'('' Entire signal'')')
          WRITE(OP_STRING(4),'('' Num. Time Samples  = '',I12)')
     '      NUMTIMEDATA
          WRITE(OP_STRING(5),'('' Num. Electrodes    = '',I12)') LIST(0)
          WRITE(OP_STRING(6),'('' RMS                = '',D12.5)') RMS
          WRITE(OP_STRING(7),'('' Relative RMS       = '',D12.5)')
     '      RELRMS

C LKC 22-NOV-2001 Adding max/min differences
          WRITE(OP_STRING(8),'('' Max. Differences   = '',D12.5)')
     '      maxsig_diff
          WRITE(OP_STRING(9),'('' Min. Differences   = '',D12.5)')
     '      minsig_diff

          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF

        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,MASTERFILE,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,COMPAREFILE,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)

      ENDIF

      CALL EXITS('COMPSIG')
      RETURN
 9999 CALL ERRORS('COMPSIG',ERROR)
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,MASTERFILE,
     '  ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9998)
 9998 CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,COMPAREFILE,
     '  ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9997)
 9997 CALL EXITS('COMPSIG')
      RETURN 1
      END


