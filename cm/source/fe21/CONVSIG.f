      SUBROUTINE CONVSIG(LD,NBJ,NDDATA,WD,XID,ZD,STRING,ERROR,*)

C#### Subroutine:CONVSIG
C###  Description:
C###    Converts a CMISS signal from ascii to binary format or
C###    vice versa.
C***  Created By    : Leo Cheng 30 November 1997

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'sign00.cmn'

!     Parameter List
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,N3CO,nr,
     '  NUMTIMEDATA,NUMTIMEDATA1,notime,SAMPLE
      REAL*8  SIGNALMIN(9),SIGNALMAX(9),TIME
      LOGICAL CBBREV,ENDFILE
      CHARACTER FORMAT_FROM*6,FORMAT_TO*6,INFILE*(MXCH),OUTFILE*(MXCH),
     '  ERROR_DUMMY*255
!     Functions
      INTEGER IFROMC

      CALL ENTERS('CONVSIG',*9999)
      CALL STRING_TRIM(FILE00,IBEG,IEND)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG1,IEND1)
C---------------------------------------------------------------------

C#### Command: FEM convert signal
C###  Parameter:      <from (ascii/binary)[ascii]>
C###    Specify whether the signal file is stored as a binary or
C###    ascii file.
C###  Parameter:      <infile FILENAME[$current]>
C###    The name of the input file with the format specified above
C###  Parameter:      <outfile FILENAME[$current]>
C###    The name of the output file with the alternative format
C###  Parameter:      <samples #[all]>
C###    The number of time steps to convert.
C###  Description:
C###    Converts a signal file between the two file formats,
C###    binary or ascii. The ASCII file will have a file
C###    extension of .sign and the BINARY file will have an
C###    extension of .binsig. SAMPLES specifies the number of
C###    time samples which are written out - useful to write out
C###    out 1 sample to examine information about a .binsig file.
C###
        OP_STRING(1)=STRING(1:IEND1)
        OP_STRING(2)=BLANK(1:15)//
     '    '<from (ascii/binary)[ascii]>'
        OP_STRING(3)=BLANK(1:15)//
     '    '<infile FILENAME['//FILE00(IBEG:IEND)//']>'
        OP_STRING(4)=BLANK(1:15)//
     '    '<outfile FILENAME['//FILE00(IBEG:IEND)//']>'
        OP_STRING(5)=BLANK(1:15)//
     '    '<samples #[all]>'

        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CONVSIG',ERROR,*9999)
      ELSE


        CALL ASSERT(NJM.GT.NJT+1,'Increase NJM to NJT+1',ERROR,*9999)


C*** I/O Files
        IF(CBBREV(CO,'INFILE',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          INFILE=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          INFILE=FILE00(IBEG:IEND)
        ENDIF

        IF(CBBREV(CO,'OUTFILE',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          OUTFILE=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          CALL STRING_TRIM(INFILE,IBEG1,IEND1)
          OUTFILE=INFILE(IBEG1:IEND1)
        ENDIF

C*** Set file format
        IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
          FORMAT_FROM='BINARY'
          FORMAT_TO='ASCII'
        ELSE
          FORMAT_FROM='ASCII'
          FORMAT_TO='BINARY'
        ENDIF

C*** Set number of sample output
        IF(CBBREV(CO,'SAMPLES',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          SAMPLE=IFROMC(CO(N3CO+1)(IBEG1:IEND1))
        ELSE
          SAMPLE=-1 !sample of -1 indicates output all
        ENDIF


C LKC 4_APR-1998 Initialise Time
        TIME=0.D0

C*** Open up the input signal files,read elect. data and reset
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT_FROM,INFILE,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT_FROM,INFILE,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT_FROM,INFILE,
     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)

        IF(FORMAT_FROM.EQ.'ASCII') THEN
          NUMTIMEDATA=0
          ENDFILE=.FALSE.
          DO WHILE(.NOT.ENDFILE)
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT_FROM,INFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
          ENDDO
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT_FROM,INFILE,
     '      'RESET',ENDFILE,.TRUE.,ERROR,*9999)
        ENDIF

C LKC 2-JUN-1998 Change to warning
C        CALL ASSERT(NUMTIMEDATA.GT.0,
C     '    '>>No Signal ?, NUMTIMEDATA.LE.0',ERROR,*9999)
        IF(NUMTIMEDATA.LE.0)THEN
          OP_STRING(1)='WARNING: No signal, NUMTIMEDATA <= 0'
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        ENDIF

C*** Open up the ouput signal file
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

C*** Create and write electrode data to the output file
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FORMAT_TO,OUTFILE,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FORMAT_TO,OUTFILE,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Convert the signal data
        IF(SAMPLE.LT.0) THEN !output all samples
          DO notime=1,NUMTIMEDATA
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT_FROM,INFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FORMAT_TO,OUTFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
          ENDDO
        ELSE
          SAMPLE=MIN(SAMPLE,NUMTIMEDATA)
          DO notime=1,SAMPLE

            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT_FROM,INFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FORMAT_TO,OUTFILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
          ENDDO !sample
        ENDIF

        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FORMAT_FROM,INFILE,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FORMAT_TO,OUTFILE,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)
      ENDIF

      CALL EXITS('CONVSIG')
      RETURN
 9999 CALL ERRORS('CONVSIG',ERROR)
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FORMAT_FROM,INFILE,
     '  ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9998)
 9998 CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FORMAT_TO,OUTFILE,
     '  ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9997)
 9997 CALL EXITS('CONVSIG')
      RETURN 1
      END


