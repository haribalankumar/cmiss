      SUBROUTINE EVTIME(NTIME_INTERP,NTIME_POINTS,NVARINDEX,EVALTIME,
     '  TIME_VALUES,VARIABLE_VALUE,STRING,TIME_VARIABLE_NAMES,OUTPUT,
     '  YP,ERROR,*)

C#### Subroutine: EVTIME
C###  Description:
C###    EVTIME evaluates a time variable value given a time.
C**** Created by Martin Buist, 18 November 1999.

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NTIME_INTERP(NTIMEVARSM),NTIME_POINTS(NTIMEVARSM),
     '  NVARINDEX
      REAL*8 EVALTIME,TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM)
     '  ,VARIABLE_VALUE,YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH),
     '  TIME_VARIABLE_NAMES(NTIMEVARSM)*(*)
      LOGICAL OUTPUT
!     Local Variables
      INTEGER i,IBEG,IBEG1,IEND,IEND1,IFROMC,INDEX,N3CO
      REAL*8 ALPHA,RFROMC
      CHARACTER FILE*(MXCH),VARIABLE_NAME*255
      LOGICAL CBBREV,OPFILE

      CALL ENTERS('EVTIME',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate time_variable<;FILENAME> variable <name>
C###  Parameter:   <time #[0.0]>
C###    Specify the time to evaluate the time variable at.
C###  Description:
C###    Evauluate time variable value at a given time. Note that
C###    the variable name is not an optional parameter.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> variable <name>'
        OP_STRING(2)=BLANK(1:15)//'time #[0.0]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVTIME',ERROR,*9999)
      ELSE
        CALL ASSERT(CALL_TIME,'>>Must define time first',ERROR,*9999)
        IF(OUTPUT) THEN
          IF(NTCOQU(noco).GT.0) THEN
            FILE=COQU(noco,1)
            CALL STRING_TRIM(FILE,IBEG,IEND)
            OPFILE=.TRUE.
            IOFI=IOFILE1
            CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.evtime','NEW',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
          ELSE
            OPFILE=.FALSE.
            IOFI=IOOP
          ENDIF

          IF(CBBREV(CO,'VARIABLE',3,noco+1,NTCO,N3CO)) THEN
            VARIABLE_NAME=CO(N3CO+1)
            CALL STRING_TRIM(VARIABLE_NAME,IBEG1,IEND1)
            NVARINDEX=0
            DO i=1,NTIMEVARST
              CALL STRING_TRIM(TIME_VARIABLE_NAMES(i),IBEG,IEND)
              IF((IEND1-IBEG1).EQ.(IEND-IBEG)) THEN
                IF(VARIABLE_NAME(IBEG1:IEND1).EQ.
     '            TIME_VARIABLE_NAMES(i)(IBEG:IEND)) NVARINDEX=i
              ENDIF
            ENDDO !i
            CALL ASSERT(NVARINDEX.NE.0,'>>Time variable not found',
     '        ERROR,*9999)
          ELSE
            ERROR='>>A variable name must be specified'
            GOTO 9999
          ENDIF

          IF(CBBREV(CO,'TIME',2,noco+1,NTCO,N3CO)) THEN
            EVALTIME=RFROMC(CO(N3CO+1))
          ELSE
            EVALTIME=0.0d0
          ENDIF
        ENDIF

        IF(CBBREV(CO,'YP',1,noco+1,NTCO,N3CO)) THEN
          INDEX=IFROMC(CO(N3CO+1))
        ELSE
          INDEX=1
        ENDIF

        VARIABLE_VALUE=0.0d0
        !Before first time point
        IF(EVALTIME.LT.TIME_VALUES(1,1,NVARINDEX)) THEN
          VARIABLE_VALUE=TIME_VALUES(2,0,NVARINDEX)
        !After last time point
        ELSE IF(EVALTIME.GT.TIME_VALUES(1,NTIME_POINTS(NVARINDEX),
     '    NVARINDEX)) THEN
          VARIABLE_VALUE=TIME_VALUES(2,(NTIME_POINTS(NVARINDEX)+1),
     '      NVARINDEX)
        ELSE
          IF(NTIME_INTERP(NVARINDEX).EQ.1) THEN
            DO i=1,NTIME_POINTS(NVARINDEX)-1
              IF((EVALTIME.GE.TIME_VALUES(1,i,NVARINDEX)).AND.
     '          (EVALTIME.LE.TIME_VALUES(1,i+1,NVARINDEX)))
     '          VARIABLE_VALUE=TIME_VALUES(2,i,NVARINDEX)
            ENDDO !i
          ELSE IF(NTIME_INTERP(NVARINDEX).EQ.2) THEN
            DO i=1,NTIME_POINTS(NVARINDEX)-1
              IF((EVALTIME.GE.TIME_VALUES(1,i,NVARINDEX)).AND.
     '          (EVALTIME.LE.TIME_VALUES(1,i+1,NVARINDEX))) THEN
                ALPHA=(TIME_VALUES(1,i+1,NVARINDEX)-EVALTIME)/
     '            (TIME_VALUES(1,i+1,NVARINDEX)-
     '            TIME_VALUES(1,i,NVARINDEX))
                VARIABLE_VALUE=(ALPHA*TIME_VALUES(2,i,NVARINDEX))+
     '            ((1.0d0-ALPHA)*TIME_VALUES(2,i+1,NVARINDEX))
              ENDIF
            ENDDO !i
          ENDIF
        ENDIF

        YP(INDEX,8,1)=VARIABLE_VALUE

        IF(OUTPUT) THEN
          WRITE(OP_STRING,'(/'' Time variable name: '',A)')
     '      VARIABLE_NAME(IBEG:IEND)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Value at time '',E12.4,'' is: '',E12.4)')
     '      EVALTIME,VARIABLE_VALUE
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          IF(OPFILE) THEN
            CALL CLOSEF(IOFI,ERROR,*9999)
            IOFI=IOOP
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('EVTIME')
      RETURN
 9999 CALL ERRORS('EVTIME',ERROR)
      CALL EXITS('EVTIME')
      RETURN 1
      END


