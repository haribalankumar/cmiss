      SUBROUTINE IPTIME(NTIME_INTERP,NTIME_POINTS,TIME_VALUES,
     '  TIME_VARIABLE_NAMES,ERROR,*)

C#### Subroutine: IPTIME
C###  Description:
C###    Inputs time variable information. At the moment this is used
C###    to input time dependent boundary conditions.
C**** Created by Martin Buist, 15 November 1999.

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'time_variable.cmn'

!     Parameter list
      INTEGER NTIME_INTERP(NTIMEVARSM),NTIME_POINTS(NTIMEVARSM)
      REAL*8 TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM)
      CHARACTER ERROR*(*),TIME_VARIABLE_NAMES(NTIMEVARSM)*(*)
!     Local variables
      INTEGER i,IBEG,ICHAR,IEND,INFO,j,NOQUES,NUMTIMES
      CHARACTER CHAR1*4
      LOGICAL FILEIP

      CALL ENTERS('IPTIME',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      INFO=0
      ICHAR=999

      IF(IOTYPE.NE.3) THEN

 802    FORMAT='($,'' Enter the name of the time '
     '    //'variable [EXIT]: '',A32)'
        CDEFLT(1)='EXIT'
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

        IF(CDATA(1).NE.'EXIT') THEN !not default exit
          NTIMEVARST=NTIMEVARST+1
          CALL ASSERT(NTIMEVARST.LE.NTIMEVARSM,
     '      '>>Increase NTIMEVARSM in *.ippara',ERROR,*9999)
          CALL STRING_TRIM(CDATA(1),IBEG,IEND)
          TIME_VARIABLE_NAMES(NTIMEVARST)=CDATA(1)(IBEG:IEND)

          DO i=1,NTIMEVARST-1
            CALL ASSERT(TIME_VARIABLE_NAMES(i).NE.
     '        TIME_VARIABLE_NAMES(NTIMEVARST),
     '        '>>Time variable names must be unique',ERROR,*9999)
          ENDDO !i

C PM 26-JUL-01
          FORMAT='(/$,'' Do you want this time variable to be'
     '      //' periodic [N]? '',A)'
          ADEFLT(1)='N'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '      IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '      ERROR,*9999)
          PERIODIC=ADATA(1)

          IDEFLT(1)=2
          FORMAT='(/$,'' Enter the number of time points to '
     '      //'be set [2]: '',I4)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          NUMTIMES=IDATA(1)
          IF(NUMTIMES.GT.NTIMEPOINTST) NTIMEPOINTST=NUMTIMES
          CALL ASSERT(NTIMEPOINTST.LE.NTIMEPOINTSM,
     '      '>>Increase NTIMEPOINTSM in *.ippara',ERROR,*9999)
          NTIME_POINTS(NTIMEVARST)=NUMTIMES
          CALL ASSERT(NTIME_POINTS(NTIMEVARST).GE.2,
     '      '>>At least 2 time points must be set',ERROR,*9999)

C PM 26-JUL-01
          IF((PERIODIC.EQ.'N').OR.(PERIODIC.EQ.'n')) THEN
C LKC 29-JAN-2000 Incorrect precision
C            RDEFLT(1)=0.0
            RDEFLT(1)=0.D0
            FORMAT='($,'' Enter the time variable value before '
     '        //'time point 1 [0.0]: '',E12.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            TIME_VALUES(2,0,NTIMEVARST)=RDATA(1)
            TIME_VALUES(1,0,NTIMEVARST)=-RMAX
          ENDIF

          DO i=1,NUMTIMES
            WRITE(CHAR1,'(I4)') i
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
C LKC 29-JAN-2000 Incorrect precision
C            RDEFLT(1)=0.0
            RDEFLT(1)=0.D0
            FORMAT='($,'' Enter the time for time point '
     '          //CHAR1(IBEG:IEND)//' [0.0]: '',E12.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            TIME_VALUES(1,i,NTIMEVARST)=RDATA(1)

C LKC 29-JAN-2000 Incorrect precision
C            RDEFLT(1)=0.0
            RDEFLT(1)=0.D0
            FORMAT='($,'' Enter the time variable value at time point '
     '          //CHAR1(IBEG:IEND)//' [0.0]: '',E12.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            TIME_VALUES(2,i,NTIMEVARST)=RDATA(1)
          ENDDO !i

          WRITE(CHAR1,'(I4)') NUMTIMES
          CALL STRING_TRIM(CHAR1,IBEG,IEND)

C PM 26-JUL-01
          IF((PERIODIC.EQ.'N').OR.(PERIODIC.EQ.'n')) THEN
C LKC 29-JAN-2000 Incorrect precision
C            RDEFLT(1)=0.0
            RDEFLT(1)=0.D0
            FORMAT='($,'' Enter the time variable value after '
     '        //'time point '//CHAR1(IBEG:IEND)//' [0.0]: '',E12.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            TIME_VALUES(2,NUMTIMES+1,NTIMEVARST)=RDATA(1)
            TIME_VALUES(1,NUMTIMES+1,NTIMEVARST)=RMAX
          ENDIF

          IDEFLT(1)=2
          FORMAT='(/'' The type of interpolation is [2]: '''//
     '      '/''   (1) Constant           '''//
     '      '/''   (2) Linear Lagrange    '''//
     '      '/''  *(3) Quadratic Lagrange '''//
     '      '/''  *(4) Cubic Lagrange     '''//
     '      '/''  *(5) Cubic Hermite      '''//
     '      '/$,''    '',I1)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          NTIME_INTERP(NTIMEVARST)=IDATA(1)

          GOTO 802
        ENDIF

      ELSE IF(IOTYPE.EQ.3) THEN

        DO i=1,NTIMEVARST
          CALL STRING_TRIM(TIME_VARIABLE_NAMES(i),IBEG,IEND)
          CDATA(1)=TIME_VARIABLE_NAMES(i)(IBEG:IEND)
          FORMAT='($,'' Enter the name of the time '
     '      //'variable [EXIT]: '',A32)'
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

C PM 26-JUL-01
          FORMAT='(/$,'' Do you want this time variable to be'
     '      //' periodic [N]? '',A)'
          ADEFLT(1)='N'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '      IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '      ERROR,*9999)
          PERIODIC=ADATA(1)

          IDATA(1)=NTIME_POINTS(i)
          FORMAT='(/$,'' Enter the number of time points to '
     '      //'be set [2]: '',I4)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

C PM 26-JUL-01
          IF((PERIODIC.EQ.'N').OR.(PERIODIC.EQ.'n')) THEN
            RDATA(1)=TIME_VALUES(2,0,i)
            FORMAT='($,'' Enter the time variable value before '
     '        //'time point 1 [0.0]: '',E12.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          ENDIF

          DO j=1,NTIME_POINTS(i)
            WRITE(CHAR1,'(I4)') j
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            RDATA(1)=TIME_VALUES(1,j,i)
            FORMAT='($,'' Enter the time for time point '
     '        //CHAR1(IBEG:IEND)//' [0.0]: '',E12.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

            RDATA(1)=TIME_VALUES(2,j,i)
            FORMAT='($,'' Enter the time variable value at time point '
     '        //CHAR1(IBEG:IEND)//' [0.0]: '',E12.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          ENDDO !j
          WRITE(CHAR1,'(I4)') NTIME_POINTS(i)
          CALL STRING_TRIM(CHAR1,IBEG,IEND)

C PM 26-JUL-01
          IF((PERIODIC.EQ.'N').OR.(PERIODIC.EQ.'n')) THEN
            RDATA(1)=TIME_VALUES(2,NTIME_POINTS(i)+1,i)
            FORMAT='($,'' Enter the time variable value after '
     '        //'time point '//CHAR1(IBEG:IEND)//' [0.0]: '',E12.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          ENDIF

          IDATA(1)=NTIME_INTERP(i)
          FORMAT='(/'' The type of interpolation is [2]: '''//
     '      '/''   (1) Constant           '''//
     '      '/''   (2) Linear Lagrange    '''//
     '      '/''  *(3) Quadratic Lagrange '''//
     '      '/''  *(4) Cubic Lagrange     '''//
     '      '/''  *(5) Cubic Hermite      '''//
     '      '/$,''    '',I1)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        ENDDO !i

        CDATA(1)='EXIT'
        FORMAT='($,'' Enter the name of the time '
     '    //'variable [EXIT]: '',A32)'
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      ENDIF !iotype

      CALL EXITS('IPTIME')
      RETURN
 9999 CALL ERRORS('IPTIME',ERROR)
      CALL EXITS('IPTIME')
      RETURN 1
      END


