      SUBROUTINE OPEN_FILES(LD,NBJ,NDDATA,NHQ,NPNY,NQNY,NRLIST,NRLIST2,
     '  NXLIST,NYNR,NYQNR,WD,XID,YP,YQ,YQS,ZD,STRING,ERROR,*)

C#### Subroutine: OPEN_FILES
C###  Description:
C###    OPEN_FILES opens files.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'sign00.cmn'
!     Parameter List
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     '  NHQ(NRM,NXM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NQNY(2,NYQM,0:NRCM,NXM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NXLIST(0:NXM),NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER CLOCAT,i,IBEG,IBEG1,IEND,IEND1,IFROMC,IPOS,IUNIT,N3CO,
     '  na,niq,niqs,niy,NTVARLIST,nx,nxc
      REAL*8 SIGNALMAX(9),SIGNALMIN(9),TIME,YPMAX(1),YPMIN(1)
      CHARACTER COMMAND*5,FILE*100,FILEFORMAT*6,TYPE*8
      CHARACTER*10 VARLIST(3)
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,ENDFILE,YPDATA,YQDATA,YQSDATA

      CALL ENTERS('OPEN_FILES',*9999)

 1    IF(CO(noco+1).EQ.'?'.OR.CO(noco+2).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        IF(CO(noco+1).NE.'?') THEN !? is second command option
C         delete first command option from string
          CALL STRING_TRIM(CO(noco+1),IBEG1,IEND1)
          IPOS=CLOCAT(CO(noco+1)(IBEG1:IEND1),STRING)
          STRING=STRING(1:IPOS-1)//STRING(IPOS+IEND1-IBEG1+1:IEND)
        ENDIF
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM open history<;FILENAME[$current]>
C###  Description:
C###    Opens a history file.
C###  Parameter: <(read/write)[read]>
C###    Specify whether the history file is being opened for reading
C###    or writing.
C###  Parameter: <unit #[20]>
C###    Specify the unit number of the history file to open.
C###  Parameter: <variables VARLIST[yp,yq,yqs]>
C###    Specify the list of variable types that the history file will
C###    contain. The list of available variable types are: yp,yq and
C###    yqs.
C###  Parameter: <niylist (#s/all)[all]>
C###    Specify the list of niy numbers of the yp data to store in the
C###    history file.
C###  Parameter: <niqlist (#s/all)[all]>
C###    Specify the list of niq numbers of the yq data to store in the
C###    history file.
C###  Parameter: <niqslist (#s/all)[all]>
C###    Specify the list of niqs numbers of the yqs data to store in
C###    the history file.
C###  Parameter: <(ascii/binary)[ascii]>
C###    Specify whether the history file to close is an ascii or
C###    binary file.
C###  Parameter: <region (#s/all)[all]>
C###  Parameter: <class #[1]>

        OP_STRING(1)=STRING(1:IEND)//' history<;FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<(read/write)[read]>'
        OP_STRING(3)=BLANK(1:15)//'<unit #[20]>'
        OP_STRING(4)=BLANK(1:15)//'<variables VARLIST[yp,yq,yqs]>'
        OP_STRING(5)=BLANK(1:15)//'<niylist (#s/all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<niqlist (#s/all)[all]>'
        OP_STRING(7)=BLANK(1:15)//'<niqslist (#s/all)[all]>'
        OP_STRING(8)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(9)=BLANK(1:15)//'<region (#s/all)[all]>'
        OP_STRING(10)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C#### Command: FEM open signal<;FILENAME[$current]>
C###  Description:
C###    Opens a signal file.
C###  Parameter: <(read/write)[read]>
C###    Specify whether the signal file is being opened for reading
C###    or writing.
C###  Parameter: <unit #[20]>
C###    Specify the unit number of the signal file to open.
C###  Parameter: <(ascii/binary)[ascii]>
C###    Specify whether the history file to close is an ascii or
C###    binary file.
C###  Parameter: <region (#s/all)[all]>

        OP_STRING(1)=STRING(1:IEND)//' signal<;FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<(read/write)[read]>'
        OP_STRING(3)=BLANK(1:15)//'<unit #[20]>'
        OP_STRING(4)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','OPEN',ERROR,*9999)
      ELSE
        IF(ABBREV(CO(noco+1),'HISTORY',1)) THEN
          TYPE='HISTORY'
        ELSE IF(ABBREV(CO(noco+1),'SIGNAL',1)) THEN
          TYPE='SIGNAL'
        ELSE
          ERROR='>>Invalid open command'
          GOTO 9999
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '    ERROR,*9999)
        CALL CHECKF(1,noco+1,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(CBBREV(CO,'WRITE',1,noco+1,NTCO,N3CO)) THEN
          COMMAND='WRITE'
        ELSE
          COMMAND='READ'
        ENDIF

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

        IF(TYPE(1:8).EQ.'HISTORY') THEN


          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)

          IF(CBBREV(CO,'VARIABLES',1,noco+1,NTCO,N3CO)) THEN
            YPDATA=.FALSE.
            YQDATA=.FALSE.
            YQSDATA=.FALSE.
            CALL PARSSL(CO(N3CO+1),3,NTVARLIST,VARLIST,ERROR,*9999)
            DO i=1,NTVARLIST
              IF(ABBREV(VARLIST(i),'YP',2)) THEN
                YPDATA=.TRUE.
              ELSE IF(ABBREV(VARLIST(i),'YQS',3)) THEN
                YQSDATA=.TRUE.
              ELSE IF(ABBREV(VARLIST(i),'YQ',2)) THEN
                YQDATA=.TRUE.
              ENDIF
            ENDDO
          ELSE
            YPDATA=.TRUE.
            YQDATA=.TRUE.
            YQSDATA=.TRUE.
          ENDIF

          IF(YPDATA) THEN
            IF(CBBREV(CO,'NIYLIST',3,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),20,NIYLIST_HIST(0,IUNIT),
     '          NIYLIST_HIST(1,IUNIT),ERROR,*9999)
            ELSE
              NIYLIST_HIST(0,IUNIT)=NIYM
              CALL ASSERT(NIYM.LE.99,
     '          '>>Increase NIYLIST_HIST in hist00.cmn',ERROR,*9999)
              DO niy=1,NIYM
                NIYLIST_HIST(niy,IUNIT)=niy
              ENDDO
            ENDIF
          ENDIF
          IF(YQDATA) THEN
            IF(CBBREV(CO,'NIQLIST',4,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),20,NIQLIST_HIST(0,IUNIT),
     '          NIQLIST_HIST(1,IUNIT),ERROR,*9999)
            ELSE
              CALL ASSERT(NIQM.LE.99,
     '          '>>Increase NIQLIST_HIST in hist00.cmn',ERROR,*9999)
              NIQLIST_HIST(0,IUNIT)=NIQM
              DO niq=1,NIQM
                NIQLIST_HIST(niq,IUNIT)=niq
              ENDDO
            ENDIF
          ENDIF
          IF(YQSDATA) THEN
            IF(CBBREV(CO,'NIQSLIST',4,noco+1,NTCO,N3CO)) THEN
C *** DPN 14 February 2000 - change max number of integers to match
C     hard-coded NIQSLIST_HIST(0:500,99) in hist00.cmn
              CALL PARSIL(CO(N3CO+1),500,NIQSLIST_HIST(0,IUNIT),
     '          NIQSLIST_HIST(1,IUNIT),ERROR,*9999)
            ELSE
              CALL ASSERT(NIQSM.LE.99,
     '          '>>Increase NIQSLIST_HIST in hist00.cmn',ERROR,*9999)
              NIQSLIST_HIST(0,IUNIT)=NIQSM
              DO niqs=1,NIQSM
                NIQSLIST_HIST(niqs,IUNIT)=niqs
              ENDDO
            ENDIF
          ENDIF

          na=1
          CALL IOHIST(IUNIT,na,NHQ(1,nx),NIQLIST_HIST(0,IUNIT),
     '      NIQSLIST_HIST(0,IUNIT),NIYLIST_HIST(0,IUNIT),
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '      NUMTIMEDATA_HIST(IUNIT),nx,NYNR(0,0,1,0,nx),
     '      NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,YPMIN,
     '      YQ(1,1,1,nx),YQS,COMMAND,FILEFORMAT,FILE,'OPEN',ENDFILE,
     '      .TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)

        ELSE IF(TYPE(1:6).EQ.'SIGNAL') THEN

          CALL IOSIGN(IUNIT,LD,NBJ,NDDATA,NUMTIMEDATA_SIG(IUNIT),
     '      0,SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,COMMAND,FILEFORMAT,
     '      FILE,'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

          CALL IOSIGN(IUNIT,LD,NBJ,NDDATA,NUMTIMEDATA_SIG(IUNIT),
     '      0,SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,COMMAND,FILEFORMAT,
     '      FILE,'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

        ENDIF

      ENDIF

      CALL EXITS('OPEN_FILES')
      RETURN
 9999 CALL ERRORS('OPEN_FILES',ERROR)
      CALL EXITS('OPEN_FILES')
      RETURN 1
      END


