      SUBROUTINE ADDSIG(LD,NBJ,NDDATA,WD,XID,ZD,STRING,ERROR,*)

C#### Subroutine: ADDSIG
C###  Description:
C###    ADDSIG creates a third signal file from a linear combination
C###    of two additional signal files i.e. signal3=a*signal1+b*signal2.

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
      INTEGER IBEG,IBEG1,IEND,IEND1,N3CO,nd,noelec,notime,nr,
     '  NUMTIMEDATA,NUMTIMEDATA1
      REAL*8 A_COEFF,B_COEFF,RFROMC,SIGNALMAX(9),SIGNALMIN(9),TIME
      CHARACTER ERROR1*255,FILEFORMAT*6,INFILE1*(MXCH),INFILE2*(MXCH),
     '  OUTFILE*(MXCH)
      LOGICAL CBBREV,ENDFILE

      CALL ENTERS('ADDSIG',*9999)

      CALL STRING_TRIM(FILE00,IBEG1,IEND1)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM add signal
C###  Parameter:      <infile1 FILENAME[$current]>
C###    The name of the first signal input file
C###  Parameter:      <a_coefficient=#[1.0]>
C###    The weighting coefficient for the first signal file
C###  Parameter:      <infile2 FILENAME[$current]>
C###    The name of the second signal input file
C###  Parameter:      <b_coefficient=#[1.0]>
C###    The weighting coefficient for the second signal file
C###  Parameter:      <outfile FILENAME[$current_new]>
C###    The output signal file
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Description:
C###    Creates a signal file with the specified outfile filename from
C###    a linear combination of the two specified signal files i.e.
C###    outfile=a_coefficient*infile1+b_coefficient*infile2.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<infile1 FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<a_coeffient=#[1.0]>'
        OP_STRING(4)=BLANK(1:15)//'<infile2 FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(5)=BLANK(1:15)//'<b_coeffient=#[1.0]>'
        OP_STRING(6)=BLANK(1:15)//'<outfile FILENAME['
     '    //FILE00(IBEG1:IEND1)//'_new]>'
        OP_STRING(7)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','ADDSIG',ERROR,*9999)
      ELSE
C LKC 24-AUG-1997 added assert line
        CALL ASSERT(NJM.GE.NJT+2,'>>Increase NJM to NJT+2',
     '        ERROR,*9999)
        IF(CBBREV(CO,'INFILE1',7,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          INFILE1=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          INFILE1=FILE00(IBEG1:IEND1)
        ENDIF

        IF(CBBREV(CO,'A_COEFFICIENT',1,noco+1,NTCO,N3CO)) THEN
          IF(N3CO.LT.NTCO) THEN
            A_COEFF=RFROMC(CO(N3CO+1))
          ELSE
            ERROR='>>No a_coefficient value supplied'
            GOTO 9999
          ENDIF
        ELSE
          A_COEFF=1.0d0
        ENDIF

        IF(CBBREV(CO,'INFILE2',7,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          INFILE2=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          INFILE2=FILE00(IBEG1:IEND1)
        ENDIF

        IF(CBBREV(CO,'B_COEFFICIENT',1,noco+1,NTCO,N3CO)) THEN
          IF(N3CO.LT.NTCO) THEN
            B_COEFF=RFROMC(CO(N3CO+1))
          ELSE
            ERROR='>>No b_coefficient value supplied'
            GOTO 9999
          ENDIF
        ELSE
          B_COEFF=1.0d0
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

C LKC 8-JAN-98 Initialise time to some value
        TIME=0.D0

C*** Open up the two input signal files and reset then to the
C*** beginning of signal data
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE2,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C LKC 9-SEP-97  assert for numtimedata moved below

        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE2,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

        CALL ASSERT(SIGNAL_NUMREGIONS(IOFILE1).EQ.
     '    SIGNAL_NUMREGIONS(IOFILE2),
     '    '>>Number of regions is not equal in input signal files',
     '    ERROR,*9999)
        DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
          CALL ASSERT(SIGNAL_NUMELEC(nr,IOFILE1).EQ.
     '      SIGNAL_NUMELEC(nr,IOFILE2),
     '      '>>Number of region electrodes is not equal in input '
     '      //'signal files',ERROR,*9999)
        ENDDO

C LKC 9-SEP-97  added to check for ascii time numtimedata
C***    IF ASCII need to calculate NUMTIMEDATA

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

        CALL ASSERT(NUMTIMEDATA.EQ.NUMTIMEDATA1,
     '    '>>Number of times is not equal in the input signal files',
     '    ERROR,*9999)

C***    Reset both files
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)

        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE2,
     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C*** Open up the ouput signal file

          SIGNAL_HOWGENERATED(IOFILE3)=SIGNAL_HOWGENERATED(IOFILE1)
          SIGNAL_HEADER(IOFILE3)=SIGNAL_HEADER(IOFILE1)
          SIGNAL_NUMREGIONS(IOFILE3)=SIGNAL_NUMREGIONS(IOFILE1)
          SIGNAL_ELEMLOC(IOFILE3)=SIGNAL_ELEMLOC(IOFILE1)
          DO nr=1,SIGNAL_NUMREGIONS(IOFILE3)
            SIGNAL_NUMELEC(nr,IOFILE3)=SIGNAL_NUMELEC(nr,IOFILE1)
            SIGNAL_REGNAME(nr,IOFILE3)=SIGNAL_REGNAME(nr,IOFILE1)
            SIGNAL_REGTYPE(nr,IOFILE3)=SIGNAL_REGTYPE(nr,IOFILE1)
            IF(SIGNAL_ELEMLOC(IOFILE3).GT.0) THEN
              SIGNAL_NUMXI(nr,IOFILE3)=SIGNAL_NUMXI(nr,IOFILE1)
            ENDIF
          ENDDO !nr

        CALL IOSIGN(IOFILE3,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C*** Write the electrode data (from the first input file)

        CALL IOSIGN(IOFILE3,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Loop over the times

        DO notime=1,NUMTIMEDATA1

C***      Read the two sets of signal data

          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE1,
     '      'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILE2,
     '      'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C***      Calculate the new signal

          DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
            DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE1)
              nd=NDDATA(noelec,nr)
              ZD(NJT+1,nd)=A_COEFF*ZD(NJT+1,nd)+B_COEFF*ZD(NJT+2,nd)
            ENDDO
          ENDDO !nr

C***      Write the new signal

          CALL IOSIGN(IOFILE3,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,OUTFILE,
     '      'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

        ENDDO !notime

C***    Close the files

        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE1,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE2,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE3,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,OUTFILE,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)

      ENDIF

      CALL EXITS('ADDSIG')
      RETURN
 9999 CALL ERRORS('ADDSIG',ERROR)
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE1,
     '  ' ',ENDFILE,.TRUE.,ERROR1,*9998)
 9998 CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILE2,
     '  ' ',ENDFILE,.TRUE.,ERROR1,*9997)
 9997 CALL IOSIGN(IOFILE3,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,OUTFILE,
     '  ' ',ENDFILE,.TRUE.,ERROR1,*9996)
 9996 CALL EXITS('ADDSIG')
      RETURN 1
      END



