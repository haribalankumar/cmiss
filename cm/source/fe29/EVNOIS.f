      SUBROUTINE EVNOIS(NBJ,LD,LIST,NDDATA,WD,XID,ZD,
     '  STRING,ERROR,*)

C#### Subroutine: EVNOIS
C###  Description:
C##       Evaluates the noise level in a set of signals. Assumes
C###      that the master is the clean signal and the comparefile
C###      is the noisy version of that file.

C*** Created by Leo Cheng 7 August 2000

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'sign00.cmn'

!     Parameter List
      INTEGER LD(NDM),LIST(0:NLISTM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IBEG,IBEG1,ICOL,IEND,IEND1,N3CO,
     '  noelec,notime,nr,NUMTIMEDATA,NUMTIMEDATA1
      REAL*8 CLEAN_F,NOISE_F,SIGNALMAX(9),SIGNALMIN(9),
     '  TIME
      CHARACTER COMPAREFNAME*(MXCH),FILEFORMAT*6,
     '  MASTERFNAME*(MXCH)
      LOGICAL CBBREV,ENDFILE,INLIST



      CALL ENTERS('EVNOIS',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------
C#### Command: FEM evaluate noise
C###  Parameter: <masterfile FILENAME>
C###    The master signal
C###  Parameter: <comparefile FILENAME>
C###    The signal file to compare with
C###  Parameter:        <electrodes (#s/all)[all]>
C###    Specify the electrode numbers to be evaluated.
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Whether the files are ascii (.ipsign) or binary (.binsig)
        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<masterfile FILENAME['//
     '    FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<comparefile FILENAME['//
     '    FILE00(IBEG1:IEND1)//']>'
        OP_STRING(4)=BLANK(1:15)//'<electrodes (#s/all)[all]>'
        OP_STRING(5)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe29','doc','EVNOIS',ERROR,*9999)
      ELSE

C***    Set File format
        IF(CBBREV(CO,'BINARY',3,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF

C***    Set master file
        IF(CBBREV(CO,'MASTERFILE',3,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          MASTERFNAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          MASTERFNAME=FILE00(IBEG:IEND)
        ENDIF

C***    Set compare file
        IF(CBBREV(CO,'COMPAREFILE',3,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          COMPAREFNAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          COMPAREFNAME=FILE00(IBEG:IEND)
        ENDIF
        CALL ASSERT(NJM.GE.NJT+2,'Increase NJM to NJT+2',ERROR,*9999)

        TIME=0.0d0 !initialise time

C***    Open master signal file
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFNAME,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
C***    Read the electrode data
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFNAME,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)
C***    IF ASCII need to calculate NUMTIMEDATA
        IF(FILEFORMAT.EQ.'ASCII') THEN
          NUMTIMEDATA=0
          ENDFILE=.FALSE.
          DO WHILE(.NOT.ENDFILE)
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFNAME,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
          ENDDO
        ENDIF
C***    Reset the master file
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFNAME,
     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C***    Open compare signal file
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFNAME,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
C***    Read the electrode data
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFNAME,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)
C***    IF ASCII need to calculate NUMTIMEDATA
        IF(FILEFORMAT.EQ.'ASCII') THEN
          NUMTIMEDATA1=0
          ENDFILE=.FALSE.
          DO WHILE(.NOT.ENDFILE)
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFNAME,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            IF(.NOT.ENDFILE) NUMTIMEDATA1=NUMTIMEDATA1+1
          ENDDO
        ENDIF
C***    Reset the compare file
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFNAME,
     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)


C*** Check that the 2 files are compatible (#electrodes, #timesteps)

        CALL ASSERT(
     '    SIGNAL_NUMREGIONS(IOFILE1).EQ.SIGNAL_NUMREGIONS(IOFILE2),
     '    'Number of regions not the same',ERROR,*9999)

        DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
          CALL ASSERT(
     '      SIGNAL_NUMELEC(nr,IOFILE1).EQ.SIGNAL_NUMELEC(nr,IOFILE2),
     '      'Number of electrodes not the same',ERROR,*9999)
        ENDDO !nr

        CALL ASSERT(NUMTIMEDATA.EQ.NUMTIMEDATA1,
     '    'Number of timesteps not the same',ERROR,*9999)
        nr=1

C*** Set subset of electrodes
        IF(CBBREV(CO,'ELECTRODES',3,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NDM,LIST(0),LIST(1),ERROR,*9999)
        ELSE !default to all
          LIST(0)=SIGNAL_NUMELEC(nr,IOFILE1)
          DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE1)
            LIST(noelec)=noelec
          ENDDO

        ENDIF
        CALL ASSERT(
     '    LIST(0).GT.0,'There are no electrodes to compare',ERROR,*9999)

C*** Do the noise calculations

        DO nr=1,SIGNAL_NUMREGIONS(IOFILE1)
          NOISE_F=0.D0
          CLEAN_F=0.D0
          DO notime=1,NUMTIMEDATA !Loop over the times

            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,MASTERFNAME,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,1,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,COMPAREFNAME,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Loop over the electrodes
            DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE1)
              IF(INLIST(noelec,LIST(1),LIST(0),icol)) THEN
                NOISE_F=NOISE_F+DABS(ZD(NJT+1,noelec)-ZD(NJT+2,noelec))
                CLEAN_F=CLEAN_F+DABS(ZD(NJT+1,noelec))
              ENDIF
            ENDDO !noelec
          ENDDO !notime

          WRITE(OP_STRING(1),'('' '')')
          WRITE(OP_STRING(2),'('' Region Number    = '',I4)') nr
          WRITE(OP_STRING(3),'('' # Time Steps     = '',I4)')
     '      NUMTIMEDATA
          WRITE(OP_STRING(4),'('' # Electrodes     = '',I4)')
     '      LIST(0)

          WRITE(OP_STRING(5),'('' '')')
          WRITE(OP_STRING(6),'('' Avg ABS noise    = '',F12.5)')
     '      NOISE_F/(NUMTIMEDATA*LIST(0))
          WRITE(OP_STRING(7),'('' % Noise          = '',F12.2)')
     '      NOISE_F/CLEAN_F*100.D0
          WRITE(OP_STRING(8),'('' SNR (db)         = '',F12.2)')
     '      20.D0*DLOG10(CLEAN_F/NOISE_F)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ENDDO ! nr

C*** Close the signal files
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,MASTERFNAME,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,COMPAREFNAME,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)

      ENDIF

      CALL EXITS('EVNOIS')
      RETURN
 9999 CALL ERRORS('EVNOIS',ERROR)
      CALL EXITS('EVNOIS')

      RETURN 1
      END


