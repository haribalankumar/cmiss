      SUBROUTINE COMBSIG(LD,NBJ,NDDATA,NRLIST,WD,XID,ZD,STRING,
     '  ERROR,*)
C#### Subroutine: COMBSIG
C###  Description:
C###    <HTML>
C###      Combined two sets of signal files into one master file.
C###
C###      <OL>
C###         <LI> TIME1 & TIME2 are the start times of each signal file
C###         <LI> LENGTH is the number of time steps to combine
C###         <LI>
C###      </OL>
C###
C###      Stores the second signals electrode information in
C###      WD2,ZD2,XID2 before appending the information to the original
C###      arrays.
C###    </HTML>

C***  Created Leo Cheng 5-MAY-1998

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'sign00.cmn'

!     Parameter List
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     '  NRLIST(0:NRM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,LENGTH,LD2(NDM),
     '  N3CO,nd,ni,nj,nodata,noelec,
     '  nr,NUMTIMEDATA,NUMTIMEDATA1,notime
      REAL*8 SIGNALMAX(9),SIGNALMIN(9),
     '  TIME,TIME1,TIME2,
     '  WD2(NJM,NDM),XID2(NIM,NDM),ZD2(NJM,NDM)
      LOGICAL ALL_REGIONS,CBBREV,ENDFILE
      CHARACTER ERROR_DUMMY*255,FILEFORMAT*6,
     '  INFILE1*(MXCH),INFILE2*(MXCH),OUTFILE*(MXCH)

!     Functions
      INTEGER IFROMC
      REAL*8 RFROMC

      CALL ENTERS('COMBSIG',*9999)
      CALL STRING_TRIM(FILE00,IBEG1,IEND1)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C#### Command: FEM combine signal
C###  Parameter:      <infile1 FILENAME[$current]>
C###    Name of the first input signal file.
C###  Parameter:      <infile2 FILENAME[$current]>
C###    Name of the second input signal file.
C###  Parameter:      <outfile FILENAME[combined]>
C###    Name of the combined output signal file.
C###  Parameter:      <(time1  #[time1 0.0]>
C###  Parameter:      <(time2  #[time2 0.0]>
C###    Defines the bounds of the time interval for which the signal file
C###    is to be combined.
C###  Parameter:      <(length #[length 10]>
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the signal file is stored as a binary or
C###    ascii file.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined. The all
C###    value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not
C###    specified will be skipped.
C###  Description:
C###    Combines two signal files with different electrode locations.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//
     '    '<infile1 FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//
     '    '<infile2 FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(4)=BLANK(1:15)//
     '    '<outfile FILENAME[combined]>'
        OP_STRING(5)=BLANK(1:15)//'<time1 #[time1 0.0]>'
        OP_STRING(6)=BLANK(1:15)//'<time2 #[time2 0.0]>'
        OP_STRING(7)=BLANK(1:15)//'<length  #[length 10]>'
        OP_STRING(8)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(9)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','COMBSIG',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'INFILE1',5,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
         INFILE1=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          INFILE1=FILE00(IBEG1:IEND1)
        ENDIF
        IF(CBBREV(CO,'INFILE2',5,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          INFILE2=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          INFILE2=FILE00(IBEG1:IEND1)
        ENDIF

        IF(CBBREV(CO,'OUTFILE',5,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          OUTFILE=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          OUTFILE='combined'
        ENDIF

        CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)
        CALL STRING_TRIM(INFILE1,IBEG,IEND)
        IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
     '    INFILE1(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
          ERROR='>>Input and output file names must be different'
          GOTO 9999
        ENDIF
        CALL STRING_TRIM(INFILE2,IBEG,IEND)
        IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
     '    INFILE2(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
          ERROR='>>Input and output file names must be different'
          GOTO 9999
        ENDIF

        IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF

        IF(CBBREV(CO,'TIME1',5,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
         TIME1=RFROMC(CO(N3CO+1)(IBEG1:IEND1))
        ELSE
          TIME1=0.D0
        ENDIF

        IF(CBBREV(CO,'TIME2',5,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
         TIME2=RFROMC(CO(N3CO+1)(IBEG1:IEND1))
        ELSE
          TIME2=0.D0
        ENDIF

        IF(CBBREV(CO,'LENGTH',3,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          LENGTH=IFROMC(CO(N3CO+1)(IBEG1:IEND1))
        ELSE
          LENGTH=10
        ENDIF


        TIME=0.D0
C*** Open up the two input signal files and reset
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE2,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C*** Read in first & second electode information
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD2,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD2,XID2,ZD2,'READ',
     '    FILEFORMAT,INFILE2,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Check the two input files are compatible
        CALL ASSERT(SIGNAL_HOWGENERATED(IOFILE1).EQ.
     '    SIGNAL_HOWGENERATED(IOFILE2),
     '    'Not both generated by CMISS',ERROR,*9999)
        CALL ASSERT(SIGNAL_NUMREGIONS(IOFILE1).EQ.
     '    SIGNAL_NUMREGIONS(IOFILE2),
     '    '>>Number of regions differ',ERROR,*9999)
        CALL ASSERT(SIGNAL_ELEMLOC(IOFILE1).EQ.
     '    SIGNAL_ELEMLOC(IOFILE2),
     '    '>> Elements not stored in both files',ERROR,*9999)
        DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
          CALL ASSERT(SIGNAL_NUMXI(nr,IOFILE1).EQ.
     '      SIGNAL_NUMXI(nr,IOFILE2),
     '      '>> Number of xi dirns not equal',ERROR,*9999)
        ENDDO


C*** Create new electrode information
        NDT=0
        DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
          DO nodata=1,SIGNAL_NUMELEC(nr,IOFILE2)
            nd=NDDATA(nodata,nr)
C*** Electrode Numbers
            NDDATA(nodata+SIGNAL_NUMELEC(nr,IOFILE1),nr)=
     '        SIGNAL_NUMELEC(nr,IOFILE1)+nodata

C*** Electrode positions & weights
            DO nj=1,NJT
              ZD(nj,nd+SIGNAL_NUMELEC(nr,IOFILE1))=ZD2(nj,nd)
              WD(nj,nd+SIGNAL_NUMELEC(nr,IOFILE1))=WD2(nj,nd)
            ENDDO
            WD(NJT+1,nd+SIGNAL_NUMELEC(nr,IOFILE1))=WD2(NJT+1,nd)

C*** Elements & Xi locations
            LD(nd+SIGNAL_NUMELEC(nr,IOFILE1))=LD2(nd)
            DO ni=1,SIGNAL_NUMXI(nr,IOFILE1)
              XID(ni,nd+SIGNAL_NUMELEC(nr,IOFILE1))=XID2(ni,nd)
            ENDDO

          ENDDO !nodata
          NDT=NDT+SIGNAL_NUMELEC(nr,IOFILE1)+SIGNAL_NUMELEC(nr,IOFILE2)
        ENDDO !nr
        CALL ASSERT(NDM.GE.NDT,'>>Increase NDM',ERROR,*9999)

C***    Reset both files
C        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
C     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
C     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)
C        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
C     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE2,
C     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)

          SIGNAL_HOWGENERATED(IOFILE3)=SIGNAL_HOWGENERATED(IOFILE1)
          SIGNAL_HEADER(IOFILE3)=SIGNAL_HEADER(IOFILE1)
          SIGNAL_NUMREGIONS(IOFILE3)=SIGNAL_NUMREGIONS(IOFILE1)
          SIGNAL_ELEMLOC(IOFILE3)=SIGNAL_ELEMLOC(IOFILE1)
          DO nr=1,SIGNAL_NUMREGIONS(IOFILE3)
            SIGNAL_NUMELEC(nr,IOFILE3)=
     '        SIGNAL_NUMELEC(nr,IOFILE1)+SIGNAL_NUMELEC(nr,IOFILE2)
            SIGNAL_REGNAME(nr,IOFILE3)=SIGNAL_REGNAME(nr,IOFILE1)
            SIGNAL_REGTYPE(nr,IOFILE3)=SIGNAL_REGTYPE(nr,IOFILE1)
            IF(SIGNAL_ELEMLOC(IOFILE3).GT.0) THEN
              SIGNAL_NUMXI(nr,IOFILE3)=SIGNAL_NUMXI(nr,IOFILE1)
            ENDIF
          ENDDO !nr
        CALL IOSIGN(IOFILE3,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C*** Write the combined electrode data
        CALL IOSIGN(IOFILE3,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

        IF(FILEFORMAT.EQ.'ASCII') THEN
          NUMTIMEDATA=0
          ENDFILE=.FALSE.
          DO WHILE(.NOT.ENDFILE)
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
          ENDDO

          NUMTIMEDATA1=0
          ENDFILE=.FALSE.
          DO WHILE(.NOT.ENDFILE)
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            IF(.NOT.ENDFILE) NUMTIMEDATA1=NUMTIMEDATA1+1
          ENDDO
        ENDIF

C This is only a crude check as you cannot tell how far into the
C file you want to start. If the length exceeds the time data
C in the file, it will just output the same information.
        CALL ASSERT(NUMTIMEDATA.GE.LENGTH,'>> LENGTH too long for file1'
     '    ,ERROR,*9999)
        CALL ASSERT(NUMTIMEDATA1.GE.LENGTH,
     '    '>> LENGTH too long for file2',ERROR,*9999)

C***    Reset both files
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE2,
     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)


C*** Get each file to the correct location
C          NUMTIMEDATA=0
          ENDFILE=.FALSE.
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME1,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '      'SIGNAL_DATA',ENDFILE,.FALSE.,ERROR,*9999)
          CALL ASSERT(.NOT.ENDFILE,
     '      '>>Could not find specified time in FILE1',
     '      ERROR,*9999)

C          NUMTIMEDATA1=0
          ENDFILE=.FALSE.
          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,1,SIGNALMAX,
     '      SIGNALMIN,TIME2,WD,XID,ZD,'READ',FILEFORMAT,INFILE2,
     '      'SIGNAL_DATA',ENDFILE,.FALSE.,ERROR,*9999)
          CALL ASSERT(.NOT.ENDFILE,
     '      '>>Could not find specified time in FILE2',
     '      ERROR,*9999)

C*** Write the currect time set
          DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
            DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE2) !2nd lot of elects
              nd=NDDATA(noelec+SIGNAL_NUMELEC(nr,IOFILE1),nr)
              ZD(NJT+1,nd+SIGNAL_NUMREGIONS(IOFILE1)-1)
     '          =ZD(NJT+2,noelec)
            ENDDO
          ENDDO

          TIME=TIME1 !using the first signal to set the timing
          CALL IOSIGN(IOFILE3,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '      'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

          DO notime=2,length
C***      Read the two sets of signal data
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '        SIGNALMIN,TIME2,WD,XID,ZD,'READ',FILEFORMAT,INFILE2,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

            DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
              DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE2) !2nd lot of elects
                nd=NDDATA(noelec+SIGNAL_NUMELEC(nr,IOFILE1),nr)
                ZD(NJT+1,nd+SIGNAL_NUMREGIONS(IOFILE1)-1)
     '            =ZD(NJT+2,noelec)
              ENDDO
            ENDDO
            CALL IOSIGN(IOFILE3,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
          ENDDO

C*** Close signal files
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE1,
     '      ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9998)
          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE2,
     '      ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9997)
          CALL IOSIGN(IOFILE3,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,OUTFILE,
     '      ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9997)

      ENDIF ! '?'

      CALL EXITS('COMBSIG')
      RETURN
 9999 CALL ERRORS('COMBSIG',ERROR)
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE1,
     '  ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9998)
 9998 CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE2,
     '  ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9997)
 9997 CALL EXITS('COMBSIG')
      RETURN 1
      END


