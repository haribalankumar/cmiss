      SUBROUTINE OPSIGN(LD,NBJ,NDDATA,WD,XID,ZD,
     '  FILEFORMAT,SIGFNAME,ERROR,*)

C#### Subroutine: OPSIGN
C###  Description:
C###    Output signal information

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'sign00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'mach00.inc'

!     Parameter List
      CHARACTER SIGFNAME*100,ERROR*(*),FILEFORMAT*6
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM)

!     Local Variables
      INTEGER notime,NUMTIMEDATA,noelec,noelec_end,noelec_start,
     '  nonr,noreject,nr,NREJECT(256)
      REAL*8 DUMMY(1),SIGNALMAX(9),SIGNALMIN(9),TIME,DELTAT
      LOGICAL ENDFILE
      CHARACTER ERROR_DUMMY*255

      CALL ENTERS('OPSIGN',*9999)

      CALL ASSERT(NJM.GE.NJT+2,'>>Increase NJM',ERROR,*9999)


C LKC 2-MAY-1999 Initialise time & dummy
      TIME=0.D0
      DUMMY(1)=0.D0

C*** Open up the two signal files
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGFNAME,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C*** Read in the electrode data
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGFNAME,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C***    IF ASCII need to calculate NUMTIMEDATA
        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
          NUMTIMEDATA=0
          ENDFILE=.FALSE.
          DO WHILE(.NOT.ENDFILE)
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGFNAME,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
          ENDDO

        ENDIF
        CALL ASSERT(NUMTIMEDATA.GE.1,
     '    '>> No timedata information in MAST signal',ERROR,*9999)


C*** Reset both files
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGFNAME,
     '    'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C*** Read in the first time
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,1,
     '    SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '    FILEFORMAT,SIGFNAME,
     '    'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
        DELTAT=TIME
        DO notime=2,NUMTIMEDATA
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,1,
     '      SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ',
     '      FILEFORMAT,SIGFNAME,
     '      'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
        ENDDO

        WRITE(OP_STRING(1),'('' '')')
        WRITE(OP_STRING(2),'('' Number of time samples '',I12)')
     '     NUMTIMEDATA
        WRITE(OP_STRING(3),'('' Start time             '',E12.5)')
     '    DELTAT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(NUMTIMEDATA.EQ.1) THEN !only 1 time step
          DELTAT=0
          WRITE(OP_STRING(1),'('' Final time             '',E12.5)')
     '      TIME
          WRITE(OP_STRING(2),'('' Time step              '',E12.5)')
     '      DELTAT
          WRITE(OP_STRING(3),'('' Frequency              '',E12.5)')
     '      0
        ELSE
          DELTAT=(TIME-DELTAT)/DBLE((NUMTIMEDATA-1))

          WRITE(OP_STRING(1),'('' Final time             '',E12.5)')
     '      TIME
          WRITE(OP_STRING(2),'('' Time step              '',E12.5)')
     '      DELTAT
          WRITE(OP_STRING(3),'('' Frequency              '',E12.5)')
     '      1.D0/DELTAT
        ENDIF
          WRITE(OP_STRING(4),'('' Total number of regions'',I12)')
     '      SIGNAL_NUMREGIONS(IOFILE1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C*** Collect electrode information
        DO nonr=1,SIGNAL_NUMREGIONS(IOFILE1)
          noreject=0
C!!! LKC 30-AUG-1999 Needs to change when data has region dependency
          nr=nonr

C LKC 15-SEP-1999 Modifications for multiple regions
C assumes that the electrodes are sequentially sorted
C
C          DO noelec=1,SIGNAL_NUMELEC(nr,IOFILE1)

          IF(nr.EQ.1) THEN
            noelec_start=1
            noelec_end=SIGNAL_NUMELEC(nr,IOFILE1)
          ELSE
            noelec_start=1+SIGNAL_NUMELEC(nr-1,IOFILE1)
            noelec_end=SIGNAL_NUMELEC(nr,IOFILE1)+
     '        SIGNAL_NUMELEC(nr-1,IOFILE1)
          ENDIF

          DO noelec=noelec_start,noelec_end
            IF(ABS(WD(NJT+1,noelec)).LT.ZERO_TOL) THEN
              noreject=noreject+1
              CALL ASSERT(noreject.LE.50,
     '          '>>Increase NREJECT array',ERROR,*9999)
              NREJECT(noreject)=noelec
            ENDIF
          ENDDO

          WRITE(OP_STRING(1),'('' '')')
          WRITE(OP_STRING(2),'('' Region '',I8)') nr
          WRITE(OP_STRING(3),'('' #Electrodes is          '',I5)')
     '      NDDATA(0,nr)
          WRITE(OP_STRING(4),'('' The Electrodes are :    '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL WRITE_LONG(INTTYPE,1,1,IOFI,NDDATA(0,nr),10,10,
     '      NDDATA(1,nr),DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)



          IF(noreject.EQ.0) THEN
            WRITE(OP_STRING(1),
     '        '('' There are no rejected electrodes. '')')
            WRITE(OP_STRING(2),'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING(1),'('' #Rejected electrodes is '',I5,'
     '        //''', Electrodes are:'')') noreject
            WRITE(OP_STRING(2),'('' '')')
            WRITE(OP_STRING(3),'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

            CALL WRITE_LONG(INTTYPE,1,1,IOFI,noreject,10,10,
     '        NREJECT,DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
          ENDIF

        ENDDO


        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,SIGFNAME,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)

      CALL EXITS('OPSIGN')

      RETURN
 9999 CALL ERRORS('OPSIGN',ERROR)
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,SIGFNAME,
     '  ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9998)

 9998 CALL EXITS('OPSIGN')
      RETURN 1
      END


