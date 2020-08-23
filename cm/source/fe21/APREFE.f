      SUBROUTINE APREFE(LD,NBJ,NDDATA,NHQ,NP_INTERFACE,NPNY,NQNY,
     '  NRLIST,NRLIST2,NXLIST,NYNP,NYNR,NYQNR,WD,XID,
     '  YP,YQ,YQS,ZD,STRING,ERROR,*)

C#### Subroutine: APREFE
C###  Description:
C###    APREFE uses the reference location given in ref000.cmn to
C###    re-reference a set of signals to.  The average of the signals
C###    at each reference location is subtracted from all the current
C###    signals in the history/signal file.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ref000.cmn'
      INCLUDE 'sign00.cmn'
!     Parameter List
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     '  NHQ(NRM,NXM),NP_INTERFACE(0:NPM,0:3),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NQNY(2,NYQM,0:NRCM,NXM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),ZD(NJM,NDM)
      CHARACTER STRING*(MXCH),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,N3CO,na,nc,nd,NIQLIST(0:1),
     '  NIQSLIST(0:1),NIYLIST(0:16),noelec,no_nr,no_ref_loc,nony,
     '  notime,nr,NUMTIMEDATA,NUMTIMEDATA1,nx,nxc,ny
      REAL*8 REF_VAL,SIGNALMAX(9),SIGNALMIN(9),TIME,YPMAX(16),YPMIN(16)
      CHARACTER ERROR1*255,FILEFORMAT*6,INFILE1*(MXCH),
     '  OUTFILE*(MXCH)
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,ENDFILE,HISTORY,
     '  YPDATA,YQDATA,YQSDATA

      CALL ENTERS('APREFE',*9999)

      CALL STRING_TRIM(FILE00,IBEG1,IEND1)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM apply reference history
C###  Parameter:      <infile FILENAME[$current]>
C###    The name of the history input file
C###  Parameter:      <outfile FILENAME[$current_new]>
C###    The output history file
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:      <using (fit/solve)[solve]>
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description:
C###    Creates a history file with the specified outfile filename from
C###    the input history file after subtracting the signal
C###    corresponding to the location(s) given in NP_REF_LOC.

        OP_STRING(1)=STRING(1:IEND)//' history <infile FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<outfile FILENAME['
     '    //FILE00(IBEG1:IEND1)//'_new]>'
        OP_STRING(3)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(4)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM apply reference signal
C###  Parameter:      <infile FILENAME[$current]>
C###    The name of the input signal file
C###  Parameter:      <outfile FILENAME[$current_new]>
C###    The output signal file
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Description:
C###    Creates a signal file with the specified outfile filename from
C###    the input signal file after subtracting the signal
C###    corresponding to the location(s) given in NP_REF_LOC.

        OP_STRING(1)=STRING(1:IEND)//' signal <infile FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<outfile FILENAME['
     '    //FILE00(IBEG1:IEND1)//'_new]>'
        OP_STRING(3)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','ADDSIG',ERROR,*9999)
      ELSE
        CALL ASSERT(NJM.GE.NJT+1,'>>Increase NJM to NJT+2',
     '    ERROR,*9997)

C LKC 11-JUN-2002
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

C GBS 7-JUN-2000 Check that reference has been defined
        CALL ASSERT(NP_REF_LOC(0).GT.0,'>>Define reference first',
     '    ERROR,*9997)
        IF(CBBREV(CO,'HISTORY',2,noco+1,NTCO,N3CO)) THEN
          HISTORY=.TRUE.
        ELSE
          HISTORY=.FALSE.
        ENDIF
        IF(CBBREV(CO,'INFILE',4,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          INFILE1=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          INFILE1=FILE00(IBEG1:IEND1)
        ENDIF

        IF(CBBREV(CO,'OUTFILE',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          OUTFILE=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          OUTFILE=FILE00(IBEG1:IEND1)//'_new'
        ENDIF

        CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)
        CALL STRING_TRIM(INFILE1,IBEG,IEND)
        IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
     '    INFILE1(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
          ERROR='>>Input and output file names must be different'
          GOTO 9997
        ENDIF

        IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF



        IF(HISTORY) THEN

C GBS 7-JUN-2000 Added <using fit/solve>
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
              CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '          ERROR,*9997)
            ELSE
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this solve class',ERROR,*9997)
            ENDIF
          ELSE
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '        ERROR,*9997)
          ENDIF

          NIYLIST(0)=0
          NIQLIST(0)=0
          NIQSLIST(0)=0
          na=1

C***    Open the INPUT history file
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '      NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,
     '      INFILE1,'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

          IF(FILEFORMAT.EQ.'ASCII') THEN !need to find NUMTIMEDATA
            NUMTIMEDATA=0
            ENDFILE=.FALSE.
            DO WHILE(.NOT.ENDFILE)
              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,
     '          NRLIST2,NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),
     '          NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,YPMIN,
     '          YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,INFILE1,'TIME_DATA',
     '          ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
              IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
            ENDDO
C*** RESET the input file
            CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '        NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '        TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '        FILEFORMAT,INFILE1,'RESET',ENDFILE,.TRUE.,YPDATA,YQDATA,
     '        YQSDATA,ERROR,*9999)
          ENDIF
          CALL ASSERT(NUMTIMEDATA.GT.0,
     '        '>>No signal, NUMTIMEDATA.LE.0',ERROR,*9999)

C***  Set up arrays for writing out.
          NRLIST(0)=NRLIST_HIST(0,IOFILE1)
          NRLIST2(0)=NRLIST_HIST(0,IOFILE1)
          DO no_nr=1,NRLIST_HIST(0,IOFILE1)
            NRLIST(no_nr)=NRLIST_HIST(no_nr,IOFILE1)
            NRLIST2(no_nr)=NRLIST_HIST(no_nr,IOFILE1)
          ENDDO
          NIYLIST(0)=NIYLIST_HIST(0,IOFILE1)
          CALL ASSERT(NIYLIST(0).LE.16,'>>Increase NIYLIST size',
     '      ERROR,*9999)
          DO nony=1,NIYLIST(0)
            NIYLIST(nony)=NIYLIST_HIST(nony,IOFILE1)
          ENDDO

          NIQLIST(0)=NIQLIST_HIST(0,IOFILE1)
          CALL ASSERT(NIQLIST(0).LE.1,'>>Increase NIQLIST size',
     '      ERROR,*9999)
          NIQLIST(1)=NIQLIST_HIST(1,IOFILE1)

          NIQSLIST(0)=NIQSLIST_HIST(0,IOFILE1)
          CALL ASSERT(NIQSLIST(0).LE.1,'>>Increase NIQSLIST size',
     '      ERROR,*9999)
          NIQSLIST(1)=NIQLIST_HIST(1,IOFILE1)

C***  Open the OUTPUT history file
          CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '      NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',FILEFORMAT,
     '      OUTFILE,'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)


          DO notime=1,NUMTIMEDATA !loop over all the time steps
C***  Read YP for each time series
            CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '        NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '        TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '        FILEFORMAT,INFILE1,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,
     '        YQDATA,YQSDATA,ERROR,*9999)

            REF_VAL=0.0d0
            DO no_ref_loc=1,NP_REF_LOC(0)
              no_nr=NP_INTERFACE(NP_REF_LOC(no_ref_loc),0)
              CALL ASSERT(no_nr.EQ.1,
     '          '>>Reference node is in more than one region',
     '          ERROR,*9999)
              nr=NP_INTERFACE(NP_REF_LOC(no_ref_loc),1)
              ny=NYNP(1,1,1,NP_REF_LOC(no_ref_loc),0,1,nr)
              REF_VAL=REF_VAL+YP(ny,1,nx)
            ENDDO
            REF_VAL=REF_VAL/DBLE(NP_REF_LOC(0))
            DO no_nr=1,NRLIST(0)
              nr=NRLIST(no_nr)
              DO nc=1,NCT(nr,nx)
                DO nony=1,NYNR(0,0,nc,nr,nx)
                  ny=NYNR(nony,0,nc,nr,nx)
                  YP(ny,1,nx)=YP(ny,1,nx)-REF_VAL
                ENDDO !no_nynr
              ENDDO !nc
            ENDDO !no_nrlist
C***    Write out new solution

            YPDATA=.TRUE.
            YQDATA=.FALSE.

            CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '        NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '        TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',
     '        FILEFORMAT,OUTFILE,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,
     '        YQDATA,YQSDATA,ERROR,*9999)
          ENDDO !nottime

C*** Close both history files
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '      NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '      INFILE1,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)
          CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '      NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '      OUTFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

        ELSE !.not.HISTORY (i.e. SIGNAL)

          TIME=0.D0
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '      'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '      'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

          IF(FILEFORMAT.EQ.'ASCII') THEN
            NUMTIMEDATA=0
            ENDFILE=.FALSE.
            DO WHILE(.NOT.ENDFILE)
              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
            ENDDO
          ENDIF

          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '      'RESET',ENDFILE,.TRUE.,ERROR,*9999)

          SIGNAL_HOWGENERATED(IOFILE2)=SIGNAL_HOWGENERATED(IOFILE1)
          SIGNAL_HEADER(IOFILE2)=SIGNAL_HEADER(IOFILE1)
          SIGNAL_NUMREGIONS(IOFILE2)=SIGNAL_NUMREGIONS(IOFILE1)
          SIGNAL_ELEMLOC(IOFILE2)=SIGNAL_ELEMLOC(IOFILE1)
          DO nr=1,SIGNAL_NUMREGIONS(IOFILE2)
            SIGNAL_NUMELEC(nr,IOFILE2)=SIGNAL_NUMELEC(nr,IOFILE1)
            SIGNAL_REGNAME(nr,IOFILE2)=SIGNAL_REGNAME(nr,IOFILE1)
            SIGNAL_REGTYPE(nr,IOFILE2)=SIGNAL_REGTYPE(nr,IOFILE1)
            IF(SIGNAL_ELEMLOC(IOFILE2).GT.0) THEN
              SIGNAL_NUMXI(nr,IOFILE2)=SIGNAL_NUMXI(nr,IOFILE1)
            ENDIF
          ENDDO !nr

          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '      'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C*** Write the electrode data (from the first input file)

          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '      'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Loop over the times

          DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
            DO notime=1,NUMTIMEDATA

              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C***      Calculate the new signal

              REF_VAL=0.0d0
              DO nd=1,NP_REF_LOC(0)
                REF_VAL=REF_VAL+ZD(NJT+1,NP_REF_LOC(nd))
              ENDDO
              REF_VAL=REF_VAL/DBLE(NP_REF_LOC(0))
              DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE1)
                nd=NDDATA(noelec,nr)
                ZD(NJT+1,nd)=ZD(NJT+1,nd)-REF_VAL
              ENDDO

              CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,
     '          SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,
     '          OUTFILE,'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

            ENDDO !notime
          ENDDO !nr

          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE1,
     '      ' ',ENDFILE,.TRUE.,ERROR,*9999)
          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE1,
     '      ' ',ENDFILE,.TRUE.,ERROR,*9999)

        ENDIF !HISTORY

      ENDIF

      CALL EXITS('APREFE')
      RETURN
 9999 CALL ERRORS('APREFE',ERROR)
C*** Close both history files correctly
      IF(HISTORY) THEN
        CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '    NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '    NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '    YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '    INFILE1,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR1,*9998)
      ELSE
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE1,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)
      ENDIF
 9998 IF(HISTORY) THEN
        CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '    NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '    NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '    YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '    OUTFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR1,*9997)
      ELSE
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,OUTFILE,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)
      ENDIF
 9997 CALL EXITS('APREFE')
      RETURN 1
      END


