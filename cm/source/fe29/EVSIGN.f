      SUBROUTINE EVSIGN(LD,NBJ,NDDATA,WD,XID,ZD,STRING,ERROR,*)

C#### Subroutine: EVSIGN
C###  Description:
C###    EVSIGN interpolates a signal file.
C###    Currently only downsamples signal in the time domain by
C###    dropping  time samples.
C**** Create Leo Cheng, July 2003

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'sign00.cmn'

      
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM)

!     Local Variables
      INTEGER FACTOR,IBEG,IEND,IBEG1,IEND1,N3CO,notime,nr,NUMTIMEDATA,
     &  NUMTIMEDATA1 
      REAL*8 SIGNALMIN(9),SIGNALMAX(9),TIME
      LOGICAL ENDFILE
      CHARACTER FORMAT*6,INFILE*(MXCH),INTERP*(4),OUTFILE*(MXCH)

!     Function
      LOGICAL CBBREV
      INTEGER IFROMC
      CALL ENTERS('EVSIGN',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG1,IEND1)      
        CALL STRING_TRIM(FILE00,IBEG,IEND)      

C---------------------------------------------------------------------

C#### Command: FEM evalulate signal
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the signal file is stored as a binary or
C###    ascii file.
C###  Parameter:      <infile FILENAME[$current]>
C###    The name of the input file.
C###  Parameter:      <outfile FILENAME[${current}_sample]>
C###    The name of the resampled output file.
C###  Parameter:      <(space/time) [time]>
C###    Whether to interpolate in space or time.
C###  Parameter:      <factor #[1]>
C###    Indicate the number of time samples to skip/interpolate
C###    The default value of 1 indicates to take every step.
C###    Negative value indicates to downsample the signal.        
C###  Description:
C###    Interpolates a signal file in the temporal or spatial domain.
C###    Currently only implemented in the temporal domain for
C###    downsampling.

        OP_STRING(1)=STRING(1:IEND1)
        OP_STRING(2)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(3)=BLANK(1:15)//
     &    '<infile FILENAME['//FILE00(IBEG:IEND)//']>'
        OP_STRING(4)=BLANK(1:15)//
     &    '<outfile FILENAME['//FILE00(IBEG:IEND)//'_sample]>'
        OP_STRING(5)=BLANK(1:15)//'<(space/time) [time]>'
        OP_STRING(6)=BLANK(1:15)//'<factor #[1]>'

        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CONVSIG',ERROR,*9999)
      ELSE

        IF(CBBREV(CO,'INFILE',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          INFILE=CO(N3CO+1)(IBEG:IEND)
        ELSE
          CALL STRING_TRIM(FILE00,IBEG,IEND)
          INFILE=FILE00(IBEG:IEND)
        ENDIF
        
        IF(CBBREV(CO,'OUTFILE',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          OUTFILE=CO(N3CO+1)(IBEG:IEND)
        ELSE
          CALL STRING_TRIM(FILE00,IBEG,IEND)
          OUTFILE=FILE00(IBEG:IEND)//'_sample'
        ENDIF

C*** Set file format
        IF(CBBREV(CO,'BINARY',3,noco+1,NTCO,N3CO)) THEN
          FORMAT='BINARY'
        ELSE
          FORMAT='ASCII'
        ENDIF

        IF(CBBREV(CO,'FACTOR',3,noco+1,NTCO,N3CO)) THEN
            FACTOR=IFROMC(CO(N3CO+1))
        ELSE
          FACTOR=1
        ENDIF

        IF(FACTOR.GT.1) THEN
          ERROR='>> EVSIGN not implemented for upsampling in time'
          GOTO 9999
        ENDIF
        
C*** Specify the domain to interpolate in
C*** Eventually will implement the BOTH option.        
        IF(CBBREV(CO,'SPACE',3,noco+1,NTCO,N3CO)) THEN
          INTERP='SPAC'
        ELSE
          INTERP='TIME'
        ENDIF
        
C*** Open up the input signal files,read elect. data and reset
        TIME=0.D0
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     &    SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT,INFILE,
     &    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     &    SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT,INFILE,
     &    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     &    SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT,INFILE,
     &    'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C*** Calculate the number of timesamples for ascii file
        IF(FORMAT.EQ.'ASCII') THEN
          NUMTIMEDATA=0
          ENDFILE=.FALSE.
          DO WHILE(.NOT.ENDFILE)
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     &        SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT,INFILE,
     &        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
          ENDDO
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     &      SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT,INFILE,
     &      'RESET',ENDFILE,.TRUE.,ERROR,*9999)
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
     &    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FORMAT,OUTFILE,
     &    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     &    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FORMAT,OUTFILE,
     &    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)


C*** Do the temporal interpolation/sampling
        IF(INTERP(1:4).EQ.'TIME') THEN
          IF(FACTOR.LT.1) THEN !downsample
            FACTOR=-FACTOR

            DO notime=1,NUMTIMEDATA
            
              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA1,0,
     '          SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',FORMAT,INFILE,
     &          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              IF(MOD(notime,FACTOR).EQ.0) THEN
                CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,0,
     &            SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FORMAT,
     &            OUTFILE,'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              ENDIF
              
            ENDDO !notime
          ELSE
            ERROR='>> Not implemented for upsampling in time'
            GOTO 9999
          ENDIF

        ENDIF !INTERP
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     &    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FORMAT,INFILE,
     &    ' ',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     &      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FORMAT,OUTFILE,
     &    ' ',ENDFILE,.TRUE.,ERROR,*9999)

      ENDIF

      CALL EXITS('EVSIGN')
      RETURN
 9999 CALL ERRORS('EVSIGN',ERROR)
      CALL EXITS('EVSIGN')
      RETURN
      END
        

