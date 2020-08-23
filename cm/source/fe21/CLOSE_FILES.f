      SUBROUTINE CLOSE_FILES(NHQ,NPNY,NQNY,NRLIST,NRLIST2,NXLIST,
     '  NYNR,NYQNR,YP,YQ,YQS,STRING,ERROR,*)

C#### Subroutine: CLOSE_FILES
C###  Description:
C###    CLOSE_FILES closes files.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NHQ(NRM,NXM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NQNY(2,NYQM,0:NRCM,NXM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NXLIST(0:NXM),NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER CLOCAT,IBEG,IBEG1,IEND,IEND1,IFROMC,IPOS,IUNIT,N3CO,
     '  na,NIQLIST(0:1),NIQSLIST(0:1),NIYLIST(0:1),NUMTIMEDATA,nx,nxc
      REAL*8 TIME,YPMAX(1),YPMIN(1)
      CHARACTER FILE*100,FILEFORMAT*6,TYPE*8
      LOGICAL ABBREV,CBBREV,ENDFILE,YPDATA,YQDATA,YQSDATA

      CALL ENTERS('CLOSE_FILES',*9999)

      IF(CO(noco+1).EQ.'?'.OR.CO(noco+2).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        IF(CO(noco+1).NE.'?') THEN !? is second command option
C         delete first command option from string
          CALL STRING_TRIM(CO(noco+1),IBEG1,IEND1)
          IPOS=CLOCAT(CO(noco+1)(IBEG1:IEND1),STRING)
          STRING=STRING(1:IPOS-1)//STRING(IPOS+IEND1-IBEG1+1:IEND)
        ENDIF
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM close [history]
C###  Parameter:      <unit #[20]>
C###    Specifies the unit number of the file to close
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Closes an ASCII/BINARY history file with unit number UNIT.

        OP_STRING(1)=STRING(1:IEND)//' [history]'
        OP_STRING(2)=BLANK(1:15)//'<unit #[20]>'
        OP_STRING(3)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CLOSE',ERROR,*9999)
      ELSE
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(ABBREV(CO(noco+1),'HISTORY',1)) THEN
          TYPE='HISTORY'
        ELSE
          TYPE='HISTORY'
        ENDIF

        IF(TYPE(1:8).EQ.'HISTORY') THEN

          IF(CBBREV(CO,'UNIT',1,noco+1,NTCO,N3CO)) THEN
            IUNIT=IFROMC(CO(N3CO+1))
          ELSE
            IUNIT=20
          ENDIF

          IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF

C LKC 4-APR-98 Initialise FILE
          FILE=' '
C LKC 12-APR-98 Initialise NIYLIST & YPDATA,YQDATA - these should not matter
          NIYLIST(0)=0
          NIYLIST(1)=0
          NIQLIST(0)=0
          NIQSLIST(0)=0
          YPDATA=.FALSE.
          YQDATA=.FALSE.
          YQSDATA=.FALSE.
          na=1

          CALL IOHIST(IUNIT,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,nx,
     '      NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,
     '      YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,FILE,' ',ENDFILE,
     '      .TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)

        ENDIF

      ENDIF

      CALL EXITS('CLOSE_FILES')
      RETURN
 9999 CALL ERRORS('CLOSE_FILES',ERROR)
      CALL EXITS('CLOSE_FILES')
      RETURN 1
      END


