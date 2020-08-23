      SUBROUTINE IOHIST(FILEID,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '  NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,TIME,YP,YPMAX,
     '  YPMIN,YQ,YQS,COMMAND,FILEFORMAT,FILENAME,SUBCOMMAND,ENDFILE,
     '  NEXTTIME,YPDATA,YQDATA,YQSDATA,ERROR,*)

C#### Subroutine: IOHIST
C###  Description:
C###    <HTML><P>
C###    IOHIST reads and writes history files. COMMAND determines
C###    whether data is to be read from ('READ') or written to
C###    ('WRITE') IUNIT. FILEFORMAT indicates whether or not the file
C###    is an 'ASCII' or a 'BINARY' file.
C###    </P><P>
C###    Available commands are:
C###    <UL>
C###    <LI> 'READ' - reading from a history file
C###         NIQLIST, NIQSLIST, and NIYLIST should list
C###         variables required from the file.
C###    <LI> 'WRITE' - writting to a history file
C###    <LI> 'CLOSE' - close a history file
C###    </UL>
C###    Available subcommands are:
C###    <UL>
C###    <LI> For 'READ'
C###         <UL>
C###         <LI> 'OPEN' - open a history file for reading.
C###         <LI> 'RESET' - reset an open history file back to the
C###              beginning of the time series data.
C###         <LI> 'TIME_DATA' - read the time data. If NEXTTIME is
C###              .TRUE. the next time data will be read. If NEXTTIME
C###              is .FALSE. the time read will be the next time in the
C###              file that is greater than or equal to the time
C###              specified by TIME.
C###         </UL>
C###    <LI> For 'WRITE'
C###         <UL>
C###         <LI> 'OPEN' - open a history file for writing
C###         <LI> 'TIME_DATA' - write the time data given
C###         </UL>
C###    </UL></P><P>
C###    ENDFILE is set when an attempt to read a 'TIME_DATA' record/tag
C###    fails because an end of file occured.
C###    </P></HTML>

C#### Comment: BINARY HISTORY FILE FORMAT
C###  Description:
C###    <HTML>
C###    For time series (history) data there are two tags.
C###    <P>
C###    TAG INDEX = 1 : Time series header information.
C###    </P>
C###    <UL>
C###    <LI>An integer to indicate the number of regions stored in the
C###        file.
C###    <LI>(number of regions)xintegers to indicate the region numbers.
C###    <LI>Integer to indicate whether or not the time series data
C###        does not contain (=0) or does contain (=1) YP data
C###    <LI>Integer to indicate whether or not the time series data
C###        does not contain (=0) or does contain (=1) YQ data
C###    <LI>If the version is later than 2 for ASCII or 1.1.0 for binary
C###        then an integer to indicate whether or not the time series
C###        data does not contain (=0) or does contain (=1) YQS data
C###    <LI>If the file contains YP then the following information is
C###        placed in the file.
C###        <UL>
C###        <LI>An integer indicating the total number of ny values
C###            stored in the file.
C###        <LI>(number of regions)xintegers to indicate the number of
C###            nc's for each region that are stored in the file.
C###        <LI>For each nc of each region an integer to indicate the
C###            number of ny's contained in the file for that nc of the
C###            region.
C###        <LI>For each nc of each region an integers indicating the
C###            number of ny's in that nc of that nr and then the list
C###            and then a sequence of integers indicating the list of
C###            ny numbers.
C###        <LI>An integer to specify the number of NIY indices stored
C###            in the file.
C###        <LI>(number of NIY values)xintegers to indicate the NIY
C###            indices that are stored in the file.
C###        <LI>An integer to specify the number of NPNY values stored
C###            in the file for each ny.
C###        <LI>A sequence of integers for the the NPNY values stored
C###            as the the number of NPNY values (see above) for each
C###            ny. Note: the sequence of ny's is given by the list
C###            of ny's (see above) for each nc of each nr.
C###        </UL>
C###    <LI>If the file contains YQ then the following information is
C###        placed in the file.
C###        <UL>
C###        <LI>An integer to specify the number of nyq's contained in
C###            the file.
C###        <LI>(number of regions)xintegers to indicate the number of
C###            nc's for each region that are stored in the file.
C###        <LI>(number of regions)xinteger to specify the number of
C###            nyq's specified in each region.
C###        <LI>An integer to specify the number of NIQ indices stored
C###            in the file.
C###        <LI>(number of NIQ values)xintegers to indicate the NIQ
C###            indices that are stored in the file.
C###        <LI>An integer to specify the na index stored in the file.
C###        <LI>An integer to specify the number of NQNY values stored
C###            in the file for each nyq.
C###        <LI>A sequence of integers for the the NQNY values stored
C###            as the the number of NQNY values (see above) for each
C###            nyq.
C###        </UL>
C###    </UL>
C###    <LI>If the file contains YQS then the following information is
C###        placed in the file.
C###        <UL>
C###        <LI>An integer to specify the number of nq's contained in
C###            the file.
C###        <LI>(number of regions)xinteger to specify the start nq
C###            number in each region.
C###        <LI>(number of regions)xinteger to specify the stop nq
C###            number in each region.
C###        <LI>An integer to specify the number of niqs indices stored
C###            in the file.
C###        <LI>(number of niqs values)xintegers to indicate the NIQS
C###            indices that are stored in the file.
C###        </UL>
C###    </UL>
C###    TAG INDEX = 2 : Time series data
C###    <UL>
C###    <LI>double precision number for the time.
C###    <LI>An integer to indicate whether or not this time contains
C###        (=1) or does not contain (=0) any YP values
C###    <LI>An integer to indicate whether or not this time contains
C###        (=1) or does not contain (=0) any YQ values
C###    <LI>If the version is later than 2 for ASCII or 1.1.0 for binary
C###        then an integer to indicate whether or not this time
C###        contains (=1) or does not contain (=0) any YQS values
C###    <LI>If this time contains YP data (number of ny values stored x
C###        number of niy values stored)xdouble precision numbers for
C###        the YP time series data. Note: the sequence of values is
C###        the list of niys for each ny (given by the list of ny's -
C###        see above) for nc of each nr.
C###    <LI>If this time contains YQ data (number of nyq values stored)x
C###        (number of niq indicies stored)x double precision numbers
C###         for the YQ time series data.
C###    <LI>If this time contains YQS data (number of niqs indicies
C###        stored)x(number of nq values stored)xdouble precision
C###        numbers for the YQS time series data.
C###    </UL>
C###    Note: The number of times stored is given by the number of tags
C###    minus 1.
C###    </HTML>


      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'binf00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'mach00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'seqf00.cmn'
!     Parameter List
      INTEGER FILEID,na,NHQ(NRM),NIQLIST(0:*),NIQSLIST(0:*),
     '  NIYLIST(0:*),NPNY(0:6,NYM,0:NRCM),NQNY(2,NYQM,0:NRCM),
     '  NRLIST(0:NRM),NRLIST2(0:NRM),NUMTIMEDATA,nx,
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM)
      REAL*8 TIME,YP(NYM,NIYM),YPMAX(*),YPMIN(*),YQ(NYQM,NIQM,NAM),
     '  YQS(NIQSM,NQM)
      CHARACTER COMMAND*(*),ERROR*(*),FILEFORMAT*(*),FILENAME*(*),
     '  SUBCOMMAND*(*)
      LOGICAL ENDFILE,NEXTTIME,YPDATA,YQDATA,YQSDATA
!     Local Variables
      INTEGER CERROR(50),COMMAPOS,ERR,FILETYPE,IBEG,IDUMMY,IEND,na1(1),
     '  nc,niq,niqs,niqpos,niqspos,niy,niypos,no_niqlist,no_niqslist,
     '  no_niylist,no_npny,noo_npny,noo_nqny,no_nqnr,no_nqny,no_nrlist,
     '  no_nynr,nq,nr,nrr,nrfile,nrpos,NUM_NPNY,NUM_NQNY,NUMSKIPBYTES,
     '  ny,nyq_index,OLDIOTYPE,VERSION(3),TOTAL,DUMMY
      REAL*8 DUMMYS(500),FILETIM(1),FILETIME,TIME2(1)
      CHARACTER LINE*132,STRING*30,FMT*500
      LOGICAL INLIST,ISBINFILEOPEN,ISENDBINFILE

      PARAMETER(NUM_NPNY=7)
      PARAMETER(NUM_NQNY=2)

      CALL ENTERS('IOHIST',*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Command : '',A,'
     '    //'/'' Subcommand : '',A,'
     '    //'/'' Fileformat : '',A,'
     '    //'/'' Filename : '',A,'
     '    //'/'' FileID : '',I5)') COMMAND,SUBCOMMAND,FILEFORMAT,
     '    FILENAME,FILEID
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IF(COMMAND.EQ.'READ') THEN

        IF(FILEFORMAT.EQ.'ASCII') THEN

C*** Read ASCII Open
          IF(SUBCOMMAND.EQ.'OPEN') THEN

            OLDIOTYPE=IOTYPE
            IOTYPE=2 !Read
            VERSION(1)=3 !11/3/99
            CALL OPEN_SEQ_FILE(VERSION(1),FILEID,FILENAME,
     '        'hist','OLD',.FALSE.,ERROR,*9999)
            IOTYPE=OLDIOTYPE

C cpb 29/3/96 There is an ambiguity when reading a history file that
C does not have the same number of nc's or a different ny/niy/nr setup
C than the current setup. You need to keep track of these numbers in
C order to read the file but do not want to overwrite the current setup.

C cpb 14/4/97
C The following is implemented: For the # of nc's and npny's the file is
C read so long as the number in the file is <= the number in
C the current setup. For the the list quantities (ie. nr list, niy list,
C ny list) the open statement will read in and set the HISTORY lists
C (e.g. NRLIST_HIST). The mapping between what is in the file and what
C gets read into what location will be determined by the supplied
C NRLIST and NRLIST2 lists. For example if the file contains two
C regions and you only want to read in one region, say the second
C region in the file into the third current region, then the supplied
C NRLIST would be NRLIST(0)=1 to indicate the number of regions to
C read. NRLIST(1)=2 to indicate that the second region in the file is
C to be read. Note that if regions in the file are not contained in
C in the NRLIST list they are skipped. The supplied NRLIST2 would be
C NRLIST2(0)=1 (must be the same as NRLIST(0)) to indicate that there
C is one region to be mapped and NRLIST2(1)=3 to indicate that the
C region specified in NRLIST(1) is to be read into the third region.

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C             READ(FILEID,'('' Number of regions stored: '',I2)')
C     '        NRLIST_HIST(0,FILEID)
           FMT='('' Number of regions stored: '',I2)'
            READ(FILEID,FMT) NRLIST_HIST(0,FILEID)


C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C            READ(FILEID,'('' Regions are:'',10(1X,I2),'
C     '        //':/(13X,10(1X,I2)))') (NRLIST_HIST(no_nrlist,FILEID),
C     '        no_nrlist=1,NRLIST_HIST(0,FILEID))
            FMT='('' Regions are:'',10(1X,I2),:/(13X,10(1X,I2)))'
            READ(FILEID,FMT) (NRLIST_HIST(no_nrlist,FILEID),
     '        no_nrlist=1,NRLIST_HIST(0,FILEID))


C LKC 14-SEP-1999 Need to check if the regions are valid
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(.NOT.INLIST(nr,NRLIST_HIST(1,FILEID),
     '          NRLIST_HIST(0,FILEID),IDUMMY)) THEN

C LKC 20-SEP-1999 may be better to have a warning rather than assert
C                ERROR='Region requested is not stored in history file'
C                GOTO 9999

                IEND=0
                CALL APPENDC(IEND,' WARNING: Region ',OP_STRING(1))
                CALL APPENDI(IEND,nr,OP_STRING(1))
                CALL APPENDC(IEND,
     '            ' requested, not stored in history file',OP_STRING(1))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO

            READ(FILEID,'(A)') LINE
            IF(SEQFILEVERSION(FILEID).LE.2) THEN
              IF(LINE(20:23).EQ.'both') THEN
                YPDATA=.TRUE.
                YQDATA=.TRUE.
              ELSE IF(LINE(20:23).EQ.'only') THEN
                IF(LINE(25:26).EQ.'YP') THEN
                  YPDATA=.TRUE.
                  YQDATA=.FALSE.
                ELSE IF(LINE(25:26).EQ.'YQ') THEN
                  YPDATA=.FALSE.
                  YQDATA=.TRUE.
                ELSE
                  ERROR='>>Invalid file format'
                  GOTO 9999
                ENDIF
              ELSE
                ERROR='>>Invalid file format'
                GOTO 9999
              ENDIF
              YQSDATA=.FALSE.
            ELSE
              YPDATA=.FALSE.
              YQDATA=.FALSE.
              YQSDATA=.FALSE.
              CALL STRING_TRIM(LINE,IBEG,IEND)
              IBEG=35
              COMMAPOS=INDEX(LINE(IBEG:IEND),',')
              DO WHILE(COMMAPOS.NE.0)
                IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YP') THEN
                  YPDATA=.TRUE.
                  IBEG=IBEG+3
                ELSE IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YQ') THEN
                  YQDATA=.TRUE.
                  IBEG=IBEG+3
                ELSE IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YQS') THEN
                  YQSDATA=.TRUE.
                  IBEG=IBEG+4
                ELSE
                  ERROR='>>Invalid file format'
                  GOTO 9999
                ENDIF
                COMMAPOS=INDEX(LINE(IBEG:IEND),',')
              ENDDO
              IF(LINE(IBEG:IEND).EQ.'YP') THEN
                YPDATA=.TRUE.
              ELSE IF(LINE(IBEG:IEND).EQ.'YQ') THEN
                YQDATA=.TRUE.
              ELSE IF(LINE(IBEG:IEND).EQ.'YQS') THEN
                YQSDATA=.TRUE.
              ELSE
                ERROR='>>Invalid file format'
                GOTO 9999
              ENDIF
            ENDIF
            IF(YPDATA) THEN
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Total number of ny values stored: '','
C     '          //'I10)') TOTALNY(FILEID)
              FMT='('' Total number of ny values stored: '',I10)'
              READ(FILEID,FMT) TOTALNY(FILEID)

              CALL ASSERT(TOTALNY(FILEID).LE.NYM,
     '          '>>Increase NYM',ERROR,*9999)

              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Region '',I2,'' : Number of '
C     '            //'nc.s stored = '',I1)') nrr,NCT_HIST(nr,FILEID)
C                READ(FILEID,'(13X,''Number of ny.s for each nc :'','
C     '            //'3(1X,I10),:/(41X,3(1X,I10)))')
C     '            (NYNR_HIST(nc,nr,FILEID),nc=1,NCT_HIST(nr,FILEID))
                FMT='('' Region '',I2,
     '            '' : Number of nc.s stored = '',I1)'
                READ(FILEID,FMT) nrr,NCT_HIST(nr,FILEID)
                FMT='(13X,''Number of ny.s for each nc :'','
     '            //'3(1X,I10),:/(41X,3(1X,I10)))'
                READ(FILEID,FMT)
     '            (NYNR_HIST(nc,nr,FILEID),nc=1,NCT_HIST(nr,FILEID))

              ENDDO !no_nrlist


C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Number of niy indicies stored = '',I2)')
C     '          NIYLIST_HIST(0,FILEID)
              FMT='('' Number of niy indicies stored = '',I2)'
              READ(FILEID,FMT) NIYLIST_HIST(0,FILEID)

              CALL ASSERT(NIYLIST_HIST(0,FILEID).LE.500,
     '          '>>Increase dimension of DUMMYS',ERROR,*9999)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Indicies are :'',20(1X,I2))')
C     '          (NIYLIST_HIST(no_niylist,FILEID),
C     '          no_niylist=1,NIYLIST_HIST(0,FILEID))
              FMT='('' Indicies are :'',20(1X,I2))'
              READ(FILEID,FMT)
     '          (NIYLIST_HIST(no_niylist,FILEID),
     '          no_niylist=1,NIYLIST_HIST(0,FILEID))



              DO no_niylist=1,NIYLIST_HIST(0,FILEID)
                niy=NIYLIST_HIST(no_niylist,FILEID)
                CALL ASSERT(niy.LE.NIYM,'>>Increase NIYM',ERROR,*9999)
              ENDDO

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Number of NPNY indicies stored = '','
C     '          //'I2)') NUM_NPNY_HIST(FILEID)
              FMT='('' Number of NPNY indicies stored = '',I2)'
              READ(FILEID,FMT) NUM_NPNY_HIST(FILEID)

              CALL ASSERT(NUM_NPNY_HIST(FILEID).LE.NUM_NPNY,'>>File '
     '          //'NUM_NPNY is different than current NUM_NPNY',
     '          ERROR,*9999)
              IF(NUM_NPNY_HIST(FILEID).NE.NUM_NPNY) THEN
                WRITE(OP_STRING,'('' >>Warning: File NUM_NPNY '
     '            //'is different than the current NUM_NPNY'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              READ(FILEID,'(A)') LINE

              DO no_npny=0,NUM_NPNY_HIST(FILEID)-1
                TOTAL=0
                DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                  nr=NRLIST_HIST(no_nrlist,FILEID)
                  DO nc=1,NCT_HIST(NR,FILEID)
                    TOTAL=TOTAL+NYNR_HIST(nc,nr,FILEID)
                  ENDDO
                ENDDO
                IF(TOTAL.LE.10) THEN


C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                  READ(FILEID,'('' NPNY index '',I1,'' :'',10(1X,I5))')
C     '              noo_npny,(IDUMMY,no_nynr=1,TOTAL)
                  FMT='('' NPNY index '',I1,'' :'',10(1X,I5))'
                  READ(FILEID,FMT) noo_npny,(IDUMMY,no_nynr=1,TOTAL)

                ELSE


C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C
C                  READ(FILEID,'('' NPNY index '',I1,'' :'',10(1X,I5),'
C     '              //':/(15X,10(1X,I5)))') noo_npny,
C     '              (((IDUMMY,no_nynr=1,NYNR_HIST(nc,
C     '              NRLIST_HIST(no_nrlist,FILEID),FILEID)),
C     '              nc=1,NCT_HIST(NRLIST_HIST(no_nrlist,FILEID),
C     '              FILEID)),no_nrlist=1,NRLIST_HIST(0,FILEID))

                  FMT='('' NPNY index '',I1,'' :'',10(1X,I5),'
     '              //':/(15X,10(1X,I5)))'
                  READ(FILEID,FMT) noo_npny, (((IDUMMY,no_nynr=1,
     '              NYNR_HIST(nc,NRLIST_HIST(no_nrlist,FILEID),FILEID)),
     '              nc=1,NCT_HIST(NRLIST_HIST(no_nrlist,FILEID),FILEID))
     '              ,no_nrlist=1,NRLIST_HIST(0,FILEID))

                ENDIF
              ENDDO !no_npny
            ENDIF
            IF(YQDATA) THEN
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Total number of nyq values stored:'
C     '          //' '',I10)') TOTALNQ(FILEID)

              FMT='('' Total number of nyq values stored: '',I10)'
              READ(FILEID,FMT) TOTALNQ(FILEID)

              CALL ASSERT(TOTALNQ(FILEID).LE.NYQM,
     '          '>>Increase NYQM',ERROR,*9999)

              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Region '',I2,'' : Number of '
C     '            //'nc values stored = '',I10))') nrr,nc

                FMT='('' Region '',I2,'' : Number of '
     '            //'nc values stored = '',I10))'
                READ(FILEID,FMT) nrr,nc
                CALL ASSERT(nc.EQ.1,'>>Number of nc values > 1',
     '            ERROR,*9999)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Region '',I2,'' : Number of '
C     '            //'nyq.s stored = '',I10))') nrr,NYQ_HIST(nr,FILEID)

                FMT='('' Region '',I2,'' : Number of '
     '            //'nyq.s stored = '',I10))'
                READ(FILEID,FMT) nrr,NYQ_HIST(nr,FILEID)

              ENDDO !no_nrlist (nr)


C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Number of niq indicies stored = '',I2)')
C     '          NIQLIST_HIST(0,FILEID)

              FMT='('' Number of niq indicies stored = '',I2)'
              READ(FILEID,FMT) NIQLIST_HIST(0,FILEID)

              CALL ASSERT(NIQLIST_HIST(0,FILEID).LE.500,
     '          '>>Increase dimension of DUMMYS',ERROR,*9999)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Indicies are :'',20(1X,I2))')
C     '          (NIQLIST_HIST(no_niqlist,FILEID),no_niqlist=1,
C     '          NIQLIST_HIST(0,FILEID))

              FMT='('' Indicies are :'',20(1X,I2))'
              READ(FILEID,FMT) (NIQLIST_HIST(no_niqlist,FILEID),
     '          no_niqlist=1,NIQLIST_HIST(0,FILEID))

              DO no_niqlist=1,NIQLIST_HIST(0,FILEID)
                niq=NIQLIST_HIST(no_niqlist,FILEID)
                CALL ASSERT(niq.LE.NIQM,'>>Increase NIQM',ERROR,*9999)
              ENDDO

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Value of na index stored = '',I2)') na
              FMT='('' Value of na index stored = '',I2)'
              READ(FILEID,FMT) na

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Number of NQNY indicies stored = '','
C     '          //'I2)') NUM_NQNY_HIST(FILEID)
              FMT='('' Number of NQNY indicies stored = '',I2)'
              READ(FILEID,FMT) NUM_NQNY_HIST(FILEID)

              CALL ASSERT(NUM_NQNY_HIST(FILEID).LE.NUM_NQNY,'>>File '
     '          //'NUM_NQNY is different than current NUM_NQNY',
     '          ERROR,*9999)
              IF(NUM_NQNY_HIST(FILEID).NE.NUM_NQNY) THEN
                WRITE(OP_STRING,'('' >>Warning: File NUM_NQNY '
     '            //'is different than the current NUM_NQNY'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              READ(FILEID,'(A)') LINE

              DO no_nqny=1,NUM_NQNY_HIST(FILEID)
                IF(TOTALNQ(FILEID).LE.10) THEN


C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                  READ(FILEID,'('' NQNY index '',I1,'' :'',10(1X,I5))')
C     '              noo_nqny,(IDUMMY,no_nqnr=1,TOTALNQ(FILEID))
                  FMT='('' NQNY index '',I1,'' :'',10(1X,I5))'
                  READ(FILEID,FMT)
     '              noo_nqny,(IDUMMY,no_nqnr=1,TOTALNQ(FILEID))

                ELSE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                  READ(FILEID,'('' NQNY index '',I1,'' :'',10(1X,I5),'
C    '              //':/(15X,10(1X,I5)))') noo_nqny,((IDUMMY,no_nqnr=1,
C     '              NYQ_HIST(NRLIST_HIST(no_nrlist,FILEID),FILEID)),
C     '              no_nrlist=1,NRLIST_HIST(0,FILEID))

                  FMT='('' NQNY index '',I1,'' :'',10(1X,I5),'
     '              //':/(15X,10(1X,I5)))'
                  READ(FILEID,FMT) noo_nqny,((IDUMMY,no_nqnr=1,
     '              NYQ_HIST(NRLIST_HIST(no_nrlist,FILEID),FILEID)),
     '              no_nrlist=1,NRLIST_HIST(0,FILEID))
                ENDIF
              ENDDO !no_nqny
            ENDIF
            IF(YQSDATA) THEN
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Total number of nq values stored = '''
C     '          //',I10)') TOTALNQS(FILEID)
              FMT='('' Total number of nq values stored = '',I10)'
              READ(FILEID,FMT) TOTALNQS(FILEID)
              CALL ASSERT(TOTALNQS(FILEID).LE.NQM,'>>Increase NQM',
     '          ERROR,*9999)

              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Region '',I2,'' : Start '
C     '            //'nq number = '',I10))') nrr,NQRS_HIST(1,nr,FILEID)

                FMT='('' Region '',I2,'' : Start nq number = '',I10)'
                READ(FILEID,FMT) nrr,NQRS_HIST(1,nr,FILEID)


C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Region '',I2,'' : Stop '
C     '            //'nq number = '',I10))') nrr,NQRS_HIST(2,nr,FILEID)

                FMT='('' Region '',I2,'' : Stop nq number = '',I10))'
                READ(FILEID,FMT) nrr,NQRS_HIST(2,nr,FILEID)
              ENDDO !no_nrlist (nr)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C             READ(FILEID,'('' Number of niqs indicies stored = '',I2)')
C     '          NIQSLIST_HIST(0,FILEID)
              FMT='('' Number of niqs indicies stored = '',I2)'
              READ(FILEID,FMT) NIQSLIST_HIST(0,FILEID)
              CALL ASSERT(NIQSLIST_HIST(0,FILEID).LE.500,
     '          '>>Increase dimension of DUMMYS',ERROR,*9999)
              IF(NIQSLIST_HIST(0,FILEID).LE.15) THEN



C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Indicies are :'',15(1X,I3))')
C     '            (NIQSLIST_HIST(no_niqslist,FILEID),no_niqslist=1,
C     '            NIQSLIST_HIST(0,FILEID))
                FMT='('' Indicies are :'',15(1X,I3))'
                READ(FILEID,FMT)
     '            (NIQSLIST_HIST(no_niqslist,FILEID),no_niqslist=1,
     '            NIQSLIST_HIST(0,FILEID))

              ELSE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Indicies are :'',15(1X,I3),'
C     '            //':/(15X,15(1X,I3)))') (NIQSLIST_HIST(no_niqslist,
C     '            FILEID),no_niqslist=1,NIQSLIST_HIST(0,FILEID))
                FMT='('' Indicies are :'',15(1X,I3),'
     '            //':/(15X,15(1X,I3)))'
                READ(FILEID,FMT) (NIQSLIST_HIST(no_niqslist,
     '            FILEID),no_niqslist=1,NIQSLIST_HIST(0,FILEID))
              ENDIF

              DO no_niqslist=1,NIQSLIST_HIST(0,FILEID)
                niqs=NIQSLIST_HIST(no_niqslist,FILEID)
                CALL ASSERT(niqs.LE.NIQSM,'>>Increase NIQSM',
     '            ERROR,*9999)
              ENDDO

            ENDIF


C*** Read ASCII Reset
          ELSE IF(SUBCOMMAND.EQ.'RESET') THEN

            CALL CLOSEF(FILEID,ERROR,*9999)
            OLDIOTYPE=IOTYPE
            IOTYPE=2 !Read
            VERSION(1)=3 !11/3/99
            CALL OPEN_SEQ_FILE(VERSION(1),FILEID,FILENAME,
     '        'hist','OLD',.FALSE.,ERROR,*9999)
            IOTYPE=OLDIOTYPE
            READ(FILEID,'(A)') LINE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C            READ(FILEID,'('' Regions are:'',10(1X,I2),'
C     '        //':/(13X,10(1X,I2)))') (IDUMMY,no_nrlist=1,
C     '        NRLIST_HIST(0,FILEID))

            FMT='('' Regions are:'',10(1X,I2),:/(13X,10(1X,I2)))'
            READ(FILEID,FMT) (IDUMMY,no_nrlist=1,NRLIST_HIST(0,FILEID))

            READ(FILEID,'(A)') LINE
            IF(SEQFILEVERSION(FILEID).LE.2) THEN
              IF(LINE(20:23).EQ.'both') THEN
                YPDATA=.TRUE.
                YQDATA=.TRUE.
              ELSE IF(LINE(20:23).EQ.'only') THEN
                IF(LINE(25:26).EQ.'YP') THEN
                  YPDATA=.TRUE.
                  YQDATA=.FALSE.
                ELSE IF(LINE(25:26).EQ.'YQ') THEN
                  YPDATA=.FALSE.
                  YQDATA=.TRUE.
                ELSE
                  ERROR='>>Invalid file format'
                  GOTO 9999
                ENDIF
              ELSE
                ERROR='>>Invalid file format'
                GOTO 9999
              ENDIF
              YQSDATA=.FALSE.
            ELSE
              YPDATA=.FALSE.
              YQDATA=.FALSE.
              YQSDATA=.FALSE.
              CALL STRING_TRIM(LINE,IBEG,IEND)
              IBEG=35
              COMMAPOS=INDEX(LINE(IBEG:IEND),',')
              DO WHILE(COMMAPOS.NE.0)
                IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YP') THEN
                  YPDATA=.TRUE.
                  IBEG=IBEG+3
                ELSE IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YQ') THEN
                  YQDATA=.TRUE.
                  IBEG=IBEG+3
                ELSE IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YQS') THEN
                  YQSDATA=.TRUE.
                  IBEG=IBEG+4
                ELSE
                  ERROR='>>Invalid file format'
                  GOTO 9999
                ENDIF
                COMMAPOS=INDEX(LINE(IBEG:IEND),',')
              ENDDO
              IF(LINE(IBEG:IEND).EQ.'YP') THEN
                YPDATA=.TRUE.
              ELSE IF(LINE(IBEG:IEND).EQ.'YQ') THEN
                YQDATA=.TRUE.
              ELSE IF(LINE(IBEG:IEND).EQ.'YQS') THEN
                YQSDATA=.TRUE.
              ELSE
                ERROR='>>Invalid file format'
                GOTO 9999
              ENDIF
            ENDIF
            IF(YPDATA) THEN
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)
                READ(FILEID,'(A)') LINE
C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'(13X,''Number of ny.s for each nc :'','
C     '            //'3(1X,I10),:/(41X,3(1X,I10)))')
C     '            (NYNR_HIST(nc,nr,FILEID),nc=1,NCT_HIST(nr,FILEID))

                FMT='(13X,''Number of ny.s for each nc :'','
     '            //'3(1X,I10),:/(41X,3(1X,I10)))'
                READ(FILEID,FMT)
     '            (NYNR_HIST(nc,nr,FILEID),nc=1,NCT_HIST(nr,FILEID))
              ENDDO !nr
              READ(FILEID,'(A)') LINE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Indicies are :'',20(1X,I2))')
C     '          (IDUMMY,no_niylist=1,NIYLIST_HIST(0,FILEID))
              FMT='('' Indicies are :'',20(1X,I2))'
              READ(FILEID,FMT)
     '          (IDUMMY,no_niylist=1,NIYLIST_HIST(0,FILEID))
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
              DO no_npny=0,NUM_NPNY_HIST(FILEID)-1
                TOTAL=0
                DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                  nr=NRLIST_HIST(no_nrlist,FILEID)
                  DO nc=1,NCT_HIST(NR,FILEID)
                    TOTAL=TOTAL+NYNR_HIST(nc,nr,FILEID)
                  ENDDO
                ENDDO

                IF(TOTAL.LE.10) THEN

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C
C                  READ(FILEID,'('' NPNY index '',I1,'' :'',10(1X,I5))')
C     '              noo_npny,(IDUMMY,no_nynr=1,TOTAL)
C                ELSE
C                  READ(FILEID,'('' NPNY index '',I1,'' :'',10(1X,I5),'
C     '              //':/(15X,10(1X,I5)))') noo_npny,
C     '              (((IDUMMY,no_nynr=1,NYNR_HIST(nc,
C     '              NRLIST_HIST(no_nrlist,FILEID),FILEID)),
C     '              nc=1,NCT_HIST(NRLIST_HIST(no_nrlist,FILEID),
C     '              FILEID)),no_nrlist=1,NRLIST_HIST(0,FILEID))

                  FMT='('' NPNY index '',I1,'' :'',10(1X,I5))'
                  READ(FILEID,FMT) noo_npny,(IDUMMY,no_nynr=1,TOTAL)
                ELSE
                  FMT='('' NPNY index '',I1,'' :'',10(1X,I5),'
     '              //':/(15X,10(1X,I5)))'
                  READ(FILEID,FMT) noo_npny,
     '              (((IDUMMY,no_nynr=1,NYNR_HIST(nc,
     '              NRLIST_HIST(no_nrlist,FILEID),FILEID)),
     '              nc=1,NCT_HIST(NRLIST_HIST(no_nrlist,FILEID),
     '              FILEID)),no_nrlist=1,NRLIST_HIST(0,FILEID))
                ENDIF
              ENDDO !no_npny
            ENDIF
            IF(YQDATA) THEN
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
              READ(FILEID,'('' Total number of nyq values stored:'
     '          //' '',I10)') TOTALNQ(FILEID)

              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Region '',I2,'' : Number of '
C     '            //'nc values stored = '',I10))') nrr,nc
                FMT='('' Region '',I2,'' : Number of '
     '            //'nc values stored = '',I10))'
                READ(FILEID,FMT) nrr,nc
                CALL ASSERT(nc.EQ.1,'>>Number of nc values > 1',
     '            ERROR,*9999)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Region '',I2,'' : Number of '
C     '            //'nyq.s stored = '',I10))') nrr,NYQ_HIST(nr,FILEID)
                FMT='('' Region '',I2,'' : Number of '
     '            //'nyq.s stored = '',I10))'
                READ(FILEID,FMT) nrr,NYQ_HIST(nr,FILEID)
              ENDDO !no_nrlist (nr)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Number of niq indicies stored = '',I2)')
C     '          NIQLIST_HIST(0,FILEID)
C              READ(FILEID,'('' Indicies are :'',20(1X,I2))')
C     '          (NIQLIST_HIST(no_niqlist,FILEID),no_niqlist=1,
C     '          NIQLIST_HIST(0,FILEID))
C              READ(FILEID,'('' Value of na index stored = '',I2)') na
C              READ(FILEID,'('' Number of NQNY indicies stored = '','
C     '          //'I2)') NUM_NQNY_HIST(FILEID)

              FMT='('' Number of niq indicies stored = '',I2)'
              READ(FILEID,FMT) NIQLIST_HIST(0,FILEID)
              FMT='('' Indicies are :'',20(1X,I2))'
              READ(FILEID,FMT)
     '          (NIQLIST_HIST(no_niqlist,FILEID),no_niqlist=1,
     '          NIQLIST_HIST(0,FILEID))
              FMT='('' Value of na index stored = '',I2)'
              READ(FILEID,FMT) na
              FMT='('' Number of NQNY indicies stored = '',I2)'
              READ(FILEID,FMT) NUM_NQNY_HIST(FILEID)

              CALL ASSERT(NUM_NQNY_HIST(FILEID).LE.NUM_NQNY,'>>File '
     '          //'NUM_NQNY is different than current NUM_NQNY',
     '          ERROR,*9999)
              IF(NUM_NQNY_HIST(FILEID).NE.NUM_NQNY) THEN
                WRITE(OP_STRING,'('' >>Warning: File NUM_NQNY '
     '            //'is different than the current NUM_NQNY'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              READ(FILEID,'(A)') LINE

              DO no_nqny=1,NUM_NQNY_HIST(FILEID)
                IF(TOTALNQ(FILEID).LE.10) THEN

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                  READ(FILEID,'('' NQNY index '',I1,'' :'',10(1X,I5))')
C     '              noo_nqny,(IDUMMY,no_nqnr=1,TOTALNQ(FILEID))
C                ELSE
C                  READ(FILEID,'('' NQNY index '',I1,'' :'',10(1X,I5),'
C    '              //':/(15X,10(1X,I5)))') noo_nqny,((IDUMMY,no_nqnr=1,
C     '              NYQ_HIST(NRLIST_HIST(no_nrlist,FILEID),FILEID)),
C     '              no_nrlist=1,NRLIST_HIST(0,FILEID))

                  FMT='('' NQNY index '',I1,'' :'',10(1X,I5))'
                  READ(FILEID,FMT)
     '              noo_nqny,(IDUMMY,no_nqnr=1,TOTALNQ(FILEID))
                ELSE
                  FMT='('' NQNY index '',I1,'' :'',10(1X,I5),'
     '              //':/(15X,10(1X,I5)))'
                  READ(FILEID,FMT) noo_nqny,((IDUMMY,no_nqnr=1,
     '              NYQ_HIST(NRLIST_HIST(no_nrlist,FILEID),FILEID)),
     '              no_nrlist=1,NRLIST_HIST(0,FILEID))

                ENDIF
              ENDDO !no_nqny
            ENDIF
            IF(YQSDATA) THEN
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE
              READ(FILEID,'(A)') LINE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(FILEID,'('' Total number of nq values stored = '''
C     '          //',I10)') TOTALNQS(FILEID)
              FMT='('' Total number of nq values stored = '',I10)'
              READ(FILEID,FMT) TOTALNQS(FILEID)

              CALL ASSERT(TOTALNQS(FILEID).LE.NQM,'>>Increase NQM',
     '          ERROR,*9999)

              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Region '',I2,'' : Start '
C     '            //'nq number = '',I10))') nrr,NQRS_HIST(1,nr,FILEID)
C                READ(FILEID,'('' Region '',I2,'' : Stop '
C     '            //'nq number = '',I10))') nrr,NQRS_HIST(2,nr,FILEID)

                FMT='('' Region '',I2,'' : Start nq number = '',I10))'
                READ(FILEID,FMT) nrr,NQRS_HIST(1,nr,FILEID)

                FMT='('' Region '',I2,'' : Stop nq number = '',I10))'
                READ(FILEID,FMT) nrr,NQRS_HIST(2,nr,FILEID)

              ENDDO !no_nrlist (nr)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C             READ(FILEID,'('' Number of niqs indicies stored = '',I2)')
C     '          NIQSLIST_HIST(0,FILEID)
              FMT='('' Number of niqs indicies stored = '',I2)'
              READ(FILEID,FMT) NIQSLIST_HIST(0,FILEID)
              CALL ASSERT(NIQSLIST_HIST(0,FILEID).LE.500,
     '          '>>Increase dimension of DUMMYS',ERROR,*9999)
              IF(NIQSLIST_HIST(0,FILEID).LE.15) THEN

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(FILEID,'('' Indicies are :'',15(1X,I3))')
C     '            (NIQSLIST_HIST(no_niqslist,FILEID),no_niqslist=1,
C     '            NIQSLIST_HIST(0,FILEID))
C              ELSE
C                READ(FILEID,'('' Indicies are :'',15(1X,I3),'
C     '            //':/(15X,15(1X,I3)))') (NIQSLIST_HIST(no_niqslist,
C     '            FILEID),no_niqslist=1,NIQSLIST_HIST(0,FILEID))

                FMT='('' Indicies are :'',15(1X,I3))'
                READ(FILEID,FMT)
     '            (NIQSLIST_HIST(no_niqslist,FILEID),no_niqslist=1,
     '            NIQSLIST_HIST(0,FILEID))
              ELSE
                FMT='('' Indicies are :'',15(1X,I3),:/(15X,15(1X,I3)))'
                READ(FILEID,FMT) (NIQSLIST_HIST(no_niqslist,FILEID),
     '            no_niqslist=1,NIQSLIST_HIST(0,FILEID))
              ENDIF
            ENDIF

C*** Read ASCII Timedata
          ELSE IF(SUBCOMMAND.EQ.'TIME_DATA') THEN

            CALL ASSERT(NRLIST(0).EQ.NRLIST2(0),
     '        '>>Invalid NRLIST/NRLIST2 setup',ERROR,*9999)

            ENDFILE=.FALSE.
            IF(TOTALNY(FILEID).GT.0) THEN
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
                  nr=NRLIST2(nrpos)
                  CALL ASSERT(NCT_HIST(nrfile,FILEID).LE.NCT(nr,nx),
     '             '>>File NCT(nr,nx) is > than NCT(nr,nx)',ERROR,*9999)
                  IF(NCT_HIST(nrfile,FILEID).NE.NCT(nr,nx)) THEN
                    WRITE(OP_STRING,'('' >>Warning: File NCT is '
     '                //'different than current NCT'')')
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                  DO nc=1,NCT_HIST(nrfile,FILEID)
                    CALL ASSERT(NYNR_HIST(nc,nrfile,FILEID).LE.
     '                NYNR(0,0,nc,nr),
     '               '>>File NYNR(0,0,nc,nr) is > than NYNR(0,0,nc,nr)',
     '                ERROR,*9999)
                    IF(NYNR_HIST(nc,nrfile,FILEID).NE.NYNR(0,0,nc,nr))
     '                THEN
                      WRITE(OP_STRING,'('' >>Warning: File '
     '                  //'NYNR(0,0,nc,nr) is different than current '
     '                  //'NYNR'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDDO !nc
                ENDIF
              ENDDO !no_nrlist
            ENDIF
            IF(TOTALNQ(FILEID).GT.0) THEN
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
                  nr=NRLIST2(nrpos)
                  CALL ASSERT(NYQ_HIST(nrfile,FILEID).LE.
     '              NYQNR(0,0,1,nr),
     '              '>>File NYQNR(0,0,1,nr) is > than NYQNR(0,0,1,nr)',
     '              ERROR,*9999)
                  IF(NYQ_HIST(nrfile,FILEID).NE.NYQNR(0,0,1,nr)) THEN
                    WRITE(OP_STRING,'('' >>Warning: File '
     '                //'NYQNR(0,0,1,nr) is different than current '
     '                //'NYQNR'')')
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDDO !no_nrlist
            ENDIF
            IF(TOTALNQS(FILEID).GT.0) THEN
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
                  nr=NRLIST2(nrpos)
                  CALL ASSERT(NQRS_HIST(1,nrfile,FILEID).EQ.NQR(1,nr),
     '              '>>File NQR(1,nr) is different than current '
     '              //'NQR(1,nr)',ERROR,*9999)
                  CALL ASSERT(NQRS_HIST(2,nrfile,FILEID).EQ.NQR(2,nr),
     '              '>>File NQR(2,nr) is different than current '
     '              //'NQR(2,nr)',ERROR,*9999)
                ENDIF
              ENDDO !no_nrlist
            ENDIF

 10         READ(FILEID,'(A)',END=9997) LINE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C            READ(FILEID,'('' Time = '',D11.4)') FILETIME
            FMT='('' Time = '',D11.4)'
            READ(FILEID,FMT) FILETIME

C cpb 14/4/97 Implementing reading specified time
            IF(NEXTTIME) THEN
              TIME=FILETIME
            ELSE
              IF(FILETIME.LT.TIME) THEN
                READ(FILEID,'(A)') LINE
                IF(SEQFILEVERSION(FILEID).LE.2) THEN
                  IF(LINE(23:26).EQ.'both') THEN
                    YPDATA=.TRUE.
                    YQDATA=.TRUE.
                  ELSE IF(LINE(23:26).EQ.'only') THEN
                    IF(LINE(28:29).EQ.'YP') THEN
                      YPDATA=.TRUE.
                      YQDATA=.FALSE.
                    ELSE IF(LINE(28:29).EQ.'YQ') THEN
                      YPDATA=.FALSE.
                      YQDATA=.TRUE.
                    ELSE
                      ERROR='>>Invalid file format'
                      GOTO 9999
                    ENDIF
                  ELSE
                    ERROR='>>Invalid file format'
                    GOTO 9999
                  ENDIF
                  YQSDATA=.FALSE.
                ELSE
                  YPDATA=.FALSE.
                  YQDATA=.FALSE.
                  YQSDATA=.FALSE.
                  CALL STRING_TRIM(LINE,IBEG,IEND)
                  IBEG=42
                  COMMAPOS=INDEX(LINE(IBEG:IEND),',')
                  DO WHILE(COMMAPOS.NE.0)
                    IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YP') THEN
                      YPDATA=.TRUE.
                      IBEG=IBEG+3
                    ELSE IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YQ') THEN
                      YQDATA=.TRUE.
                      IBEG=IBEG+3
                    ELSE IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YQS') THEN
                      YQSDATA=.TRUE.
                      IBEG=IBEG+4
                    ELSE
                      ERROR='>>Invalid file format'
                      GOTO 9999
                    ENDIF
                    COMMAPOS=INDEX(LINE(IBEG:IEND),',')
                  ENDDO
                  IF(LINE(IBEG:IEND).EQ.'YP') THEN
                    YPDATA=.TRUE.
                  ELSE IF(LINE(IBEG:IEND).EQ.'YQ') THEN
                    YQDATA=.TRUE.
                  ELSE IF(LINE(IBEG:IEND).EQ.'YQS') THEN
                    YQSDATA=.TRUE.
                  ELSE
                    ERROR='>>Invalid file format'
                    GOTO 9999
                  ENDIF
                ENDIF
                READ(FILEID,'(A)') LINE
                IF(YPDATA) THEN
                  DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                    nr=NRLIST_HIST(no_nrlist,FILEID)
                    DO nc=1,NCT_HIST(nr,FILEID)
                      DO no_nynr=1,NYNR_HIST(nc,nr,FILEID)
                        IF(NIYLIST_HIST(0,FILEID).LE.4) THEN


C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                          READ(FILEID,'('' YP('',I7,'',niy):'','
C     '                      //'4(1X,D12.5))') DUMMY,(DUMMY,
C     '                      no_niylist=1,NIYLIST_HIST(0,FILEID))
                          FMT='('' YP('',I7,'',niy):'','
     '                      //'4(1X,D12.5))'
                          READ(FILEID,FMT) DUMMY,(DUMMY,
     '                      no_niylist=1,NIYLIST_HIST(0,FILEID))
                        ELSE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                          READ(FILEID,'('' YP('',I7,'',niy):'','
C     '                      //'4(1X,D12.5),:/(17X,4(1X,D12.5)))')
C     '                      DUMMY,(DUMMY,no_niylist=1,
C     '                      NIYLIST_HIST(0,FILEID))
                          FMT='('' YP('',I7,'',niy):'','
     '                      //'4(1X,D12.5),:/(17X,4(1X,D12.5)))'
                          READ(FILEID,FMT)
     '                      DUMMY,(DUMMY,no_niylist=1,
     '                      NIYLIST_HIST(0,FILEID))

                        ENDIF
                      ENDDO !no_nynr
                    ENDDO !nc
                  ENDDO !no_nrlist
                ENDIF
                IF(YQDATA) THEN
                  DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                    nr=NRLIST_HIST(no_nrlist,FILEID)
                    DO no_nqnr=1,NYQ_HIST(nr,FILEID)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                      READ(FILEID,'('' YQ('',I7,'',niq):'','
C     '                  //'4(1X,D12.5),:/(17X,4(1X,D12.5)))') DUMMY,
C     '                  (DUMMY,niq=1,NIQLIST_HIST(0,FILEID))
                      FMT='('' YQ('',I7,'',niq):'','
     '                  //'4(1X,D12.5),:/(17X,4(1X,D12.5)))'
                      READ(FILEID,FMT) DUMMY,
     '                  (DUMMY,niq=1,NIQLIST_HIST(0,FILEID))

                    ENDDO
                  ENDDO !no_nrlist (nr)
                ENDIF
                IF(YQSDATA) THEN
                  DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                    nr=NRLIST_HIST(no_nrlist,FILEID)
                    DO nq=NQRS_HIST(1,nr,FILEID),NQRS_HIST(2,nr,FILEID)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                      READ(FILEID,'('' YQS(niqs,'',I7,''):'','
C     '                  //'4(1X,D12.5),:/(17X,4(1X,D12.5)))') DUMMY,
C     '                  (DUMMY,niqs=1,NIQSLIST_HIST(0,FILEID))
                      FMT='('' YQS(niqs,'',I7,''):'','
     '                  //'4(1X,D12.5),:/(17X,4(1X,D12.5)))'
                      READ(FILEID,FMT) DUMMY,
     '                  (DUMMY,niqs=1,NIQSLIST_HIST(0,FILEID))

                    ENDDO !nq
                  ENDDO !no_nrlist (nr)
                ENDIF
                GOTO 10
              ELSE
                TIME=FILETIME
              ENDIF
            ENDIF

            READ(FILEID,'(A)') LINE
            IF(SEQFILEVERSION(FILEID).LE.2) THEN
              IF(LINE(23:26).EQ.'both') THEN
                YPDATA=.TRUE.
                YQDATA=.TRUE.
              ELSE IF(LINE(23:26).EQ.'only') THEN
                IF(LINE(28:29).EQ.'YP') THEN
                  YPDATA=.TRUE.
                  YQDATA=.FALSE.
                ELSE IF(LINE(28:29).EQ.'YQ') THEN
                  YPDATA=.FALSE.
                  YQDATA=.TRUE.
                ELSE
                  ERROR='>>Invalid file format'
                  GOTO 9999
                ENDIF
              ELSE
                ERROR='>>Invalid file format'
                GOTO 9999
              ENDIF
              YQSDATA=.FALSE.
            ELSE
              YPDATA=.FALSE.
              YQDATA=.FALSE.
              YQSDATA=.FALSE.
              CALL STRING_TRIM(LINE,IBEG,IEND)
              IBEG=42
              COMMAPOS=INDEX(LINE(IBEG:IEND),',')
              DO WHILE(COMMAPOS.NE.0)
                IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YP') THEN
                  YPDATA=.TRUE.
                  IBEG=IBEG+3
                ELSE IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YQ') THEN
                  YQDATA=.TRUE.
                  IBEG=IBEG+3
                ELSE IF(LINE(IBEG:IBEG+COMMAPOS-2).EQ.'YQS') THEN
                  YQSDATA=.TRUE.
                  IBEG=IBEG+4
                ELSE
                  ERROR='>>Invalid file format'
                  GOTO 9999
                ENDIF
                COMMAPOS=INDEX(LINE(IBEG:IEND),',')
              ENDDO
              IF(LINE(IBEG:IEND).EQ.'YP') THEN
                YPDATA=.TRUE.
              ELSE IF(LINE(IBEG:IEND).EQ.'YQ') THEN
                YQDATA=.TRUE.
              ELSE IF(LINE(IBEG:IEND).EQ.'YQS') THEN
                YQSDATA=.TRUE.
              ELSE
                ERROR='>>Invalid file format'
                GOTO 9999
              ENDIF
            ENDIF
            IF(YPDATA) THEN
              DO no_niylist=1,NIYLIST(0)
                niy=NIYLIST(no_niylist)
                YPMIN(niy)=RMAX
                YPMAX(niy)=-RMAX
              ENDDO !niy
              READ(FILEID,'(A)') LINE
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
C                 Read region
                  nr=NRLIST2(nrpos)
                  DO nc=1,NCT_HIST(nrfile,FILEID)
                    DO no_nynr=1,NYNR_HIST(nc,nrfile,FILEID)
                      ny=NYNR(no_nynr,0,nc,nr)
                      IF(NIYLIST_HIST(0,FILEID).LE.4) THEN

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                        READ(FILEID,'('' YP('',I7,'',niy):'','
C     '                    //'4(1X,D12.5))') DUMMY,(DUMMYS(no_niylist),
C     '                    no_niylist=1,NIYLIST_HIST(0,FILEID))
                        FMT='('' YP('',I7,'',niy):'','
     '                    //'4(1X,D12.5))'
                        READ(FILEID,FMT) DUMMY,(DUMMYS(no_niylist),
     '                    no_niylist=1,NIYLIST_HIST(0,FILEID))
                      ELSE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                        READ(FILEID,'('' YP('',I7,'',niy):'','
C     '                    //'4(1X,D12.5),:/(17X,5(1X,D12.5)))') DUMMY,
C     '                    (DUMMYS(no_niylist),no_niylist=1,
C     '                    NIYLIST_HIST(0,FILEID))
                        FMT='('' YP('',I7,'',niy):'','
     '                    //'4(1X,D12.5),:/(17X,5(1X,D12.5)))'
                        READ(FILEID,FMT) DUMMY,
     '                    (DUMMYS(no_niylist),no_niylist=1,
     '                    NIYLIST_HIST(0,FILEID))
                      ENDIF
                      DO no_niylist=1,NIYLIST_HIST(0,FILEID)
                        niy=NIYLIST_HIST(no_niylist,FILEID)
                        IF(INLIST(niy,NIYLIST(1),NIYLIST(0),niypos))
     '                    THEN !Store into YP
                          YP(ny,niy)=DUMMYS(no_niylist)
                          IF(YP(ny,niy).GT.YPMAX(niy))
     '                      YPMAX(niy)=YP(ny,niy)
                          IF(YP(ny,niy).LT.YPMIN(niy))
     '                      YPMIN(niy)=YP(ny,niy)
                        ENDIF
                      ENDDO !niy
                    ENDDO !no_nynr (ny)
                  ENDDO !nc
                ELSE !Skip region
                  DO nc=1,NCT_HIST(nrfile,FILEID)
                    DO no_nynr=1,NYNR_HIST(nc,nrfile,FILEID)
                      IF(NIYLIST_HIST(0,FILEID).LE.4) THEN

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                        READ(FILEID,'('' YP('',I7,'',niy):'','
C     '                    //'4(1X,D12.5))') DUMMY,(DUMMY,
C     '                    no_niylist=1,NIYLIST_HIST(0,FILEID))
                        FMT='('' YP('',I7,'',niy):'',4(1X,D12.5))'
                        READ(FILEID,FMT) DUMMY,(DUMMY,
     '                    no_niylist=1,NIYLIST_HIST(0,FILEID))
                      ELSE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                        READ(FILEID,'('' YP('',I7,'',niy):'','
C     '                    //'4(1X,D12.5),:/(17X,4(1X,D12.5)))')
C     '                    DUMMY,(DUMMY,no_niylist=1,
C     '                    NIYLIST_HIST(0,FILEID))
                        FMT='('' YP('',I7,'',niy):'','
     '                    //'4(1X,D12.5),:/(17X,4(1X,D12.5)))'
                        READ(FILEID,FMT)
     '                    DUMMY,(DUMMY,no_niylist=1,
     '                    NIYLIST_HIST(0,FILEID))
                      ENDIF
                    ENDDO !no_nynr
                  ENDDO !nc
                ENDIF
              ENDDO !no_nrlist
            ENDIF
            IF(YQDATA) THEN
              READ(FILEID,'(A)') LINE
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
C                 Read region
                  nr=NRLIST2(nrpos)
                  DO no_nqnr=1,NYQ_HIST(nr,FILEID)
                    nyq_index=NYQNR(no_nqnr,0,1,nr)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                    READ(FILEID,'('' YQ('',I7,'',niq):'',4(1X,D12.5),'
C     '                //':/(17X,4(1X,D12.5)))') DUMMY,
C     '                (DUMMYS(no_niqlist),no_niqlist=1,
C     '                NIQLIST_HIST(0,FILEID))
                    FMT='('' YQ('',I7,'',niq):'',4(1X,D12.5),'
     '                //':/(17X,4(1X,D12.5)))'
                    READ(FILEID,FMT) DUMMY,
     '                (DUMMYS(no_niqlist),no_niqlist=1,
     '                NIQLIST_HIST(0,FILEID))

                    DO no_niqlist=1,NIQLIST_HIST(0,FILEID)
                      niq=NIQLIST_HIST(no_niqlist,FILEID)
                      IF(INLIST(niq,NIQLIST(1),NIQLIST(0),niqpos)) THEN
                        YQ(nyq_index,niq,na)=DUMMYS(no_niqlist)
                      ENDIF
                    ENDDO !no_niqlist (niq)
                  ENDDO
                ELSE !Skip region
                  DO no_nqnr=1,NYQ_HIST(nr,FILEID)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                    READ(FILEID,'('' YQ('',I7,'',niq):'','
C     '                //'4(1X,D12.5),:/(17X,4(1X,D12.5)))') DUMMY,
C     '                (DUMMY,niq=1,NIQLIST_HIST(0,FILEID))
                    FMT='('' YQ('',I7,'',niq):'','
     '                //'4(1X,D12.5),:/(17X,4(1X,D12.5)))'
                    READ(FILEID,FMT) DUMMY,
     '                (DUMMY,niq=1,NIQLIST_HIST(0,FILEID))
                  ENDDO
                ENDIF
              ENDDO !no_nrlist (nr)
            ENDIF
            IF(YQSDATA) THEN
              READ(FILEID,'(A)') LINE
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
C                 Read region
                  nr=NRLIST2(nrpos)
                  DO nq=NQRS_HIST(1,nr,FILEID),NQRS_HIST(2,nr,FILEID)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                   READ(FILEID,'('' YQS(niqs,'',I7,''):'',4(1X,D12.5),'
C     '                //':/(17X,4(1X,D12.5)))') DUMMY,
C     '                (DUMMYS(no_niqslist),no_niqslist=1,
C     '                NIQSLIST_HIST(0,FILEID))
                    FMT='('' YQS(niqs,'',I7,''):'',4(1X,D12.5),'
     '                //':/(17X,4(1X,D12.5)))'
                    READ(FILEID,FMT) DUMMY,
     '                (DUMMYS(no_niqslist),no_niqslist=1,
     '                NIQSLIST_HIST(0,FILEID))

                    DO no_niqslist=1,NIQSLIST_HIST(0,FILEID)
                      niqs=NIQSLIST_HIST(no_niqslist,FILEID)
                      IF(INLIST(niqs,NIQSLIST(1),NIQSLIST(0),niqspos))
     '                  THEN
                        YQS(niqs,nq)=DUMMYS(no_niqslist)
                      ENDIF
                    ENDDO !no_niqlist (niqs)
                  ENDDO
                ELSE !Skip region
                  DO nq=NQRS_HIST(1,nr,FILEID),NQRS_HIST(2,nr,FILEID)


C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
                    FMT='('' YQS(niqs,'',I7,''):'','
     '                //'4(1X,D12.5),:/(17X,4(1X,D12.5)))'
                    READ(FILEID,FMT) DUMMY,
     '                (DUMMY,niqs=1,NIQSLIST_HIST(0,FILEID))

                  ENDDO
                ENDIF
              ENDDO !no_nrlist (nr)
            ENDIF
          ELSE
            ERROR='>>Invalid subcommand'
            GOTO 9999
          ENDIF

        ELSE IF(FILEFORMAT.EQ.'BINARY') THEN

C*** Read Binary Open
          IF(SUBCOMMAND.EQ.'OPEN') THEN

C***        Open the history file
            FILETYPE=2 !Binary history file
            VERSION(1)=1 !version 11/3/99
            VERSION(2)=2
            VERSION(3)=0
            CALL OPEN_BIN_FILE(FILEID,FILETYPE,NUMBERTAGS,VERSION,
     '        'READ','his',FILENAME,ERROR,*9999)

C***        Read the header information tag.
            CALL READ_BIN_TAG_HEADER(FILEID,ERROR,*9999)
            CALL ASSERT(TAGINDEX.EQ.1,
     '        '>>Time series header tag not found',ERROR,*9999)

C cpb 13/11/95 There is an ambiguity when reading a history file that
C does not have the same number of nc's or a different ny/niy/nr setup
C than the current setup. You need to keep track of these numbers in
C order to read the file but do not want to overwrite the current setup.

C cpb 14/4/97
C The following is implemented: For the # of nc's and npny's the file is
C read so long as the number in the file is <= the number in
C the current setup. For the the list quantities (ie. nr list, niy list,
C ny list) the open statement will read in and set the HISTORY lists
C (e.g. NRLIST_HIST). The mapping between what is in the file and what
C gets read into what location will be determined by the supplied lists.
C NRLIST would be NRLIST(0)=1 to indicate the number of regions to
C read. NRLIST(1)=2 to indicate that the second region in the file is
C to be read. Note that if regions in the file are not contained in
C in the NRLIST list they are skipped. The supplied NRLIST2 would be
C NRLIST2(0)=1 (must be the same as NRLIST(0)) to indicate that there
C is one region to be mapped and NRLIST2(1)=3 to indicate that the
C region specified in NRLIST(1) is to be read into the third region.

C***        Read the number of regions in the file
            CALL BINREADFILE(FILEID,INTTYPE,1,NRLIST_HIST(0,FILEID),
     '        REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C***        Read the number of regions and the region list
            CALL BINREADFILE(FILEID,INTTYPE,NRLIST_HIST(0,FILEID),
     '        NRLIST_HIST(1,FILEID),REAL4DATA,REAL8DATA,CHARDATA,
     '        LOGDATA,SINTDATA,ERROR,*9999)


C LKC 14-SEP-1999 Need to check if the regions are valid
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(.NOT.INLIST(nr,NRLIST_HIST(1,FILEID),
     '          NRLIST_HIST(0,FILEID),IDUMMY)) THEN

C LKC 20-SEP-1999 may be better to have a warning rather than assert
C                ERROR='Region requested is not stored in history file'
C                GOTO 9999

                IEND=0
                CALL APPENDC(IEND,' WARNING: Region ',OP_STRING(1))
                CALL APPENDI(IEND,nr,OP_STRING(1))
                CALL APPENDC(IEND,
     '            ' requested, not stored in history file',OP_STRING(1))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)


              ENDIF
            ENDDO


C***        Read in whether YP data present or not
            CALL BINREADFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            IF(INTDATA(1).EQ.1) THEN
              YPDATA=.TRUE.
            ELSE
              YPDATA=.FALSE.
            ENDIF
C***        Read in whether YQ data present or not
            CALL BINREADFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            IF(INTDATA(1).EQ.1) THEN
              YQDATA=.TRUE.
            ELSE
              YQDATA=.FALSE.
            ENDIF
            IF(BINFILEVERSION(1,FILEID).GE.1.AND.
     '        BINFILEVERSION(2,FILEID).GE.2.AND.
     '        BINFILEVERSION(3,FILEID).GE.0) THEN
C***          Read in whether YQS data present or not
              CALL BINREADFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              IF(INTDATA(1).EQ.1) THEN
                YQSDATA=.TRUE.
              ELSE
                YQSDATA=.FALSE.
              ENDIF
            ELSE
              YQSDATA=.FALSE.
            ENDIF
            IF(YPDATA) THEN
C***          Read in the total number of ny's stored
              CALL BINREADFILE(FILEID,INTTYPE,1,TOTALNY(FILEID),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)
              CALL ASSERT(TOTALNY(FILEID).LE.NYM,
     '          '>>Increase NYM',ERROR,*9999)

C***          Read in the total number of nc's for each region
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)
                CALL BINREADFILE(FILEID,INTTYPE,1,NCT_HIST(nr,FILEID),
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
              ENDDO !no_nrlist (nr)
C***          Read in the list of ny's for each nc for each region
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)
                DO nc=1,NCT_HIST(nr,FILEID)
C***              Read in the number of ny's for each nc for each region
                  CALL BINREADFILE(FILEID,INTTYPE,1,
     '              NYNR_HIST(nc,nr,FILEID),REAL4DATA,REAL8DATA,
     '              CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***              Skip the list of ny's for each nc for each region
                  NUMSKIPBYTES=NYNR_HIST(nc,nr,FILEID)*INTSIZE
                  CALL BINSKIPFILE(FILEID,NUMSKIPBYTES,ERROR,*9999)
                ENDDO !nc
              ENDDO !no_nrlist (nr)
C***          Read in the number of NIY indices stored in the file
              CALL BINREADFILE(FILEID,INTTYPE,1,NIYLIST_HIST(0,FILEID),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
C***          Read in the number of NIY index list stored in the file
              CALL BINREADFILE(FILEID,INTTYPE,NIYLIST_HIST(0,FILEID),
     '          NIYLIST_HIST(1,FILEID),REAL4DATA,REAL8DATA,CHARDATA,
     '          LOGDATA,SINTDATA,ERROR,*9999)

              DO no_niylist=1,NIYLIST_HIST(0,FILEID)
                niy=NIYLIST_HIST(no_niylist,FILEID)
                CALL ASSERT(niy.LE.NIYM,'>>Increase NIYM',
     '            ERROR,*9999)
              ENDDO

C***          Read in the number of NPNY values (second dimension of
C***          NPNY).
              CALL BINREADFILE(FILEID,INTTYPE,1,NUM_NPNY_HIST(FILEID),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
              CALL ASSERT(NUM_NPNY_HIST(FILEID).LE.NUM_NPNY,'>>File '
     '          //'NUM_NPNY is different than current NUM_NPNY',
     '          ERROR,*9999)
              IF(NUM_NPNY_HIST(FILEID).NE.NUM_NPNY) THEN
                WRITE(OP_STRING,'('' >>Warning: File NUM_NPNY '
     '            //'is different than the current NUM_NPNY'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
C***          Skip in NPNY
              NUMSKIPBYTES=0
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)
                DO nc=1,NCT_HIST(nr,FILEID)
                  NUMSKIPBYTES=NUMSKIPBYTES+NYNR_HIST(nc,nr,FILEID)*
     '              NUM_NPNY_HIST(FILEID)*INTSIZE
                ENDDO !nc
              ENDDO !no_nrlist (nr)
              CALL BINSKIPFILE(FILEID,NUMSKIPBYTES,ERROR,*9999)
            ENDIF
            IF(YQDATA) THEN
C*** Read in the total number of nyq's stored
              CALL BINREADFILE(FILEID,INTTYPE,1,TOTALNQ(FILEID),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)
              CALL ASSERT(TOTALNQ(FILEID).LE.NYQM,'>>Increase NYQM',
     '          ERROR,*9999)

              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)
C*** Read in the number of nc's for the current region
                CALL BINREADFILE(FILEID,INTTYPE,1,na1,
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
              ENDDO

              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)
C*** Read in the number of nyq's for the current region
                CALL BINREADFILE(FILEID,INTTYPE,1,NYQ_HIST(nr,FILEID),
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
              ENDDO

C*** Read in the number of niq indices
              CALL BINREADFILE(FILEID,INTTYPE,1,NIQLIST_HIST(0,FILEID),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
C*** Read in the niq index values
              CALL BINREADFILE(FILEID,INTTYPE,NIQLIST_HIST(0,FILEID),
     '          NIQLIST_HIST(1,FILEID),REAL4DATA,REAL8DATA,CHARDATA,
     '          LOGDATA,SINTDATA,ERROR,*9999)

              DO no_niqlist=1,NIQLIST_HIST(0,FILEID)
                niq=NIQLIST_HIST(no_niqlist,FILEID)
                CALL ASSERT(niq.LE.NIQM,'>>Increase NIQM',
     '            ERROR,*9999)
              ENDDO

C*** Read in the multigrid level index
              CALL BINREADFILE(FILEID,INTTYPE,1,na1,
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
              na=na1(1)

C***          Write out the number of NQNY indicies for each nyq
              CALL BINREADFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              NUM_NQNY_HIST(FILEID)=INTDATA(1)
              CALL ASSERT(NUM_NQNY_HIST(FILEID).LE.NUM_NQNY,'>>File '
     '          //'NUM_NQNY is different than current NUM_NQNY',
     '          ERROR,*9999)
              IF(NUM_NQNY_HIST(FILEID).NE.NUM_NQNY) THEN
                WRITE(OP_STRING,'('' >>Warning: File NUM_NQNY '
     '            //'is different than the current NUM_NQNY'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF

C***          Skip in NPNY
              NUMSKIPBYTES=0
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)
                NUMSKIPBYTES=NUMSKIPBYTES+NYQ_HIST(nr,FILEID)*
     '            NUM_NQNY_HIST(FILEID)*INTSIZE
              ENDDO !no_nrlist (nr)
              CALL BINSKIPFILE(FILEID,NUMSKIPBYTES,ERROR,*9999)
            ENDIF
            IF(YQSDATA) THEN
C*** Read in the total number of nq's stored
              CALL BINREADFILE(FILEID,INTTYPE,1,TOTALNQS(FILEID),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)
              CALL ASSERT(TOTALNQS(FILEID).LE.NQM,'>>Increase NQM',
     '          ERROR,*9999)

              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)
C*** Read in the start nq for the current region
                CALL BINREADFILE(FILEID,INTTYPE,1,
     '            NQRS_HIST(1,nr,FILEID),REAL4DATA,REAL8DATA,
     '            CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              ENDDO
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nr=NRLIST_HIST(no_nrlist,FILEID)
C*** Read in the stop nq for the current region
                CALL BINREADFILE(FILEID,INTTYPE,1,
     '            NQRS_HIST(2,nr,FILEID),REAL4DATA,REAL8DATA,
     '            CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              ENDDO

C*** Read in the number of niq indices
              CALL BINREADFILE(FILEID,INTTYPE,1,NIQSLIST_HIST(0,FILEID),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
C*** Read in the niq index values
              CALL BINREADFILE(FILEID,INTTYPE,NIQSLIST_HIST(0,FILEID),
     '          NIQSLIST_HIST(1,FILEID),REAL4DATA,REAL8DATA,CHARDATA,
     '          LOGDATA,SINTDATA,ERROR,*9999)

              DO no_niqslist=1,NIQSLIST_HIST(0,FILEID)
                niqs=NIQSLIST_HIST(no_niqslist,FILEID)
                CALL ASSERT(niqs.LE.NIQSM,'>>Increase NIQSM',
     '            ERROR,*9999)
              ENDDO

            ENDIF

            NUMTIMEDATA=NUMBERTAGS-1

C*** Read Binary Reset
          ELSE IF(SUBCOMMAND.EQ.'RESET') THEN

C           Reset file to beginning of time series (history) data

            CALL BINSETFILE(FILEID,0,ERROR,*9999)
            CALL BINSKIPCMHEADER(FILEID,3,ERROR,*9999)
C***        Read time series header tag
            CALL READ_BIN_TAG_HEADER(FILEID,ERROR,*9999)
            CALL ASSERT(TAGINDEX.EQ.1,
     '        '>>Time series header tag not found',ERROR,*9999)
            CALL BINSKIPFILE(FILEID,NUMTAGBYTES,ERROR,*9999)

C*** Read Binary Timedata
          ELSE IF(SUBCOMMAND.EQ.'TIME_DATA') THEN

            CALL ASSERT(NRLIST(0).EQ.NRLIST2(0),
     '        '>>Invalid NRLIST/NRLIST2 setup',ERROR,*9999)

            ENDFILE=.FALSE.
            IF(TOTALNY(FILEID).GT.0) THEN
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
                  nr=NRLIST2(nrpos)
                  CALL ASSERT(NCT_HIST(nrfile,FILEID).LE.NCT(nr,nx),
     '              '>>File NCT(nr,nx) is > than NCT(nr,nx)',
     '              ERROR,*9999)
                  IF(NCT_HIST(nrfile,FILEID).NE.NCT(nr,nx)) THEN
                    WRITE(OP_STRING,'('' >>Warning: File NCT is '
     '                //'different than current NCT'')')
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                  DO nc=1,NCT_HIST(nrfile,FILEID)
                    CALL ASSERT(NYNR_HIST(nc,nrfile,FILEID).LE.
     '                NYNR(0,0,nc,nr),
     '               '>>File NYNR(0,0,nc,nr) is > than NYNR(0,0,nc,nr)',
     '                ERROR,*9999)
                    IF(NYNR_HIST(nc,nrfile,FILEID).NE.
     '                NYNR(0,0,nc,nr)) THEN
                      WRITE(OP_STRING,'('' >>Warning: File '
     '                  //'NYNR(0,0,nc,nr) is different than current '
     '                  //'NYNR'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDDO !nc
                ENDIF
              ENDDO !no_nrlist
            ENDIF
            IF(TOTALNQ(FILEID).GT.0) THEN
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
                  nr=NRLIST2(nrpos)
                  CALL ASSERT(NYQ_HIST(nrfile,FILEID).LE.
     '              NYQNR(0,0,1,nr),
     '              '>>File NYQNR(0,0,1,nr) is > than NYQNR(0,0,1,nr)',
     '              ERROR,*9999)
                  IF(NYQ_HIST(nrfile,FILEID).NE.NYQNR(0,0,1,nr)) THEN
                    WRITE(OP_STRING,'('' >>Warning: File '
     '                //'NYQNR(0,0,1,nr) is different than current '
     '                //'NYQNR'')')
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDDO !no_nrlist
            ENDIF
            IF(TOTALNQS(FILEID).GT.0) THEN
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
                  nr=NRLIST2(nrpos)
                  CALL ASSERT(NQRS_HIST(1,nrfile,FILEID).EQ.NQR(1,nr),
     '              '>>File NQR(1,nr) is different than current '
     '              //'NQR(1,nr)',ERROR,*9999)
                  CALL ASSERT(NQRS_HIST(2,nrfile,FILEID).EQ.NQR(2,nr),
     '              '>>File NQR(2,nr) is different than current '
     '              //'NQR(2,nr)',ERROR,*9999)
                ENDIF
              ENDDO !no_nrlist
            ENDIF

C***        Read the tag header information
 20         CALL READ_BIN_TAG_HEADER(FILEID,ERROR,*9999)
            CALL ASSERT(TAGINDEX.EQ.2,
     '        '>>Time series data tag not found',ERROR,*9999)
C***        Read in the time
            FILETIM(1)=FILETIME
            CALL BINREADFILE(FILEID,DPTYPE,1,INTDATA,REAL4DATA,
     '        FILETIM,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            FILETIME=FILETIM(1)
C cpb 14/4/97 Implementing reading specified time
            IF(.NOT.NEXTTIME) THEN
              IF(FILETIME.LT.TIME) THEN
                CALL BINSKIPFILE(FILEID,NUMTAGBYTES-DPSIZE,ERROR,*9999)
                IF(ISENDBINFILE(FILEID)) GOTO 9997
                GOTO 20
              ELSE
                TIME=FILETIME
              ENDIF
            ELSE
              TIME=FILETIME
            ENDIF

C***        Read in whether YP data present or not
            CALL BINREADFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            IF(INTDATA(1).EQ.1) THEN
              YPDATA=.TRUE.
            ELSE
              YPDATA=.FALSE.
            ENDIF
C***        Read in whether YQ data present or not
            CALL BINREADFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            IF(INTDATA(1).EQ.1) THEN
              YQDATA=.TRUE.
            ELSE
              YQDATA=.FALSE.
            ENDIF
            IF(BINFILEVERSION(1,FILEID).GE.1.AND.
     '        BINFILEVERSION(2,FILEID).GE.2.AND.
     '        BINFILEVERSION(3,FILEID).GE.0) THEN
C***          Read in whether YQS data present or not
              CALL BINREADFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              IF(INTDATA(1).EQ.1) THEN
                YQSDATA=.TRUE.
              ELSE
                YQSDATA=.FALSE.
              ENDIF
            ELSE
              YQSDATA=.FALSE.
            ENDIF
            IF(YPDATA) THEN
C***          Read in the YP time data
              DO no_niylist=1,NIYLIST(0)
                niy=NIYLIST(no_niylist)
                YPMIN(niy)=RMAX
                YPMAX(niy)=-RMAX
              ENDDO !no_niylist (niy)
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
                  nr=NRLIST2(nrpos)
                  DO nc=1,NCT_HIST(nrfile,FILEID)
                    DO no_nynr=1,NYNR_HIST(nc,nrfile,FILEID)
                      ny=NYNR(no_nynr,0,nc,nr)
                      DO no_niylist=1,NIYLIST_HIST(0,FILEID)
                        niy=NIYLIST_HIST(no_niylist,FILEID)
                        IF(INLIST(niy,NIYLIST(1),NIYLIST(0),niypos))
     '                    THEN !Read in value
                          CALL BINREADFILE(FILEID,DPTYPE,1,INTDATA,
     '                      REAL4DATA,YP(ny,niy),CHARDATA,LOGDATA,
     '                      SINTDATA,ERROR,*9999)
                          IF(YP(ny,niy).GT.YPMAX(niy))
     '                      YPMAX(niy)=YP(ny,niy)
                          IF(YP(ny,niy).LT.YPMIN(niy))
     '                      YPMIN(niy)=YP(ny,niy)
                        ELSE !Skip value
                          CALL BINSKIPFILE(FILEID,DPSIZE,ERROR,*9999)
                        ENDIF
                      ENDDO !no_niylist (niy)
                    ENDDO !no_nynr (ny)
                  ENDDO !nc
                ELSE ! Skip region
                  DO nc=1,NCT_HIST(nrfile,FILEID)
                    NUMSKIPBYTES=NYNR_HIST(nc,nrfile,FILEID)*
     '                NIYLIST_HIST(0,FILEID)*DPSIZE
                    CALL BINSKIPFILE(FILEID,NUMSKIPBYTES,ERROR,*9999)
                  ENDDO !nc
                ENDIF
              ENDDO !no_nrlist (nr)
            ENDIF
            IF(YQDATA) THEN
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
                  nr=NRLIST2(nrpos)
                  DO no_nqnr=1,NYQ_HIST(nrfile,FILEID)
                    nyq_index=NYQNR(no_nqnr,0,1,nr)
                    DO no_niqlist=1,NIQLIST_HIST(0,FILEID)
                      niq=NIQLIST_HIST(no_niqlist,FILEID)
                      IF(INLIST(niq,NIQLIST(1),NIQLIST(0),niqpos)) THEN
                        CALL BINREADFILE(FILEID,DPTYPE,1,INTDATA,
     '                    REAL4DATA,YQ(nyq_index,niq,na),CHARDATA,
     '                    LOGDATA,SINTDATA,ERROR,*9999)
                      ELSE !Skip value
                        CALL BINSKIPFILE(FILEID,DPSIZE,ERROR,*9999)
                      ENDIF
                    ENDDO !niq
                  ENDDO !no_nqnr (nyq)
                ELSE !Skip region
                  NUMSKIPBYTES=NYQ_HIST(nrfile,FILEID)*
     '              NIQLIST_HIST(0,FILEID)*DPSIZE
                  CALL BINSKIPFILE(FILEID,NUMSKIPBYTES,ERROR,*9999)
                ENDIF
              ENDDO !no_nrlist (nr)
            ENDIF
            IF(YQSDATA) THEN
              DO no_nrlist=1,NRLIST_HIST(0,FILEID)
                nrfile=NRLIST_HIST(no_nrlist,FILEID)
                IF(INLIST(nrfile,NRLIST(1),NRLIST(0),nrpos)) THEN
                  nr=NRLIST2(nrpos)
                  DO nq=NQRS_HIST(1,nr,FILEID),NQRS_HIST(2,nr,FILEID)
                    DO no_niqslist=1,NIQSLIST_HIST(0,FILEID)
                      niqs=NIQSLIST_HIST(no_niqslist,FILEID)
                      IF(INLIST(niqs,NIQSLIST(1),NIQSLIST(0),niqspos))
     '                  THEN
                        CALL BINREADFILE(FILEID,DPTYPE,1,INTDATA,
     '                    REAL4DATA,YQS(niqs,nq),CHARDATA,
     '                    LOGDATA,SINTDATA,ERROR,*9999)
                      ELSE !Skip value
                        CALL BINSKIPFILE(FILEID,DPSIZE,ERROR,*9999)
                      ENDIF
                    ENDDO !niqs
                  ENDDO !nq
                ELSE !Skip region
                  NUMSKIPBYTES=(NQRS_HIST(2,nrfile,FILEID)-
     '              NQRS_HIST(1,nrfile,FILEID)+1)*
     '              NIQSLIST_HIST(0,FILEID)*DPSIZE
                  CALL BINSKIPFILE(FILEID,NUMSKIPBYTES,ERROR,*9999)
                ENDIF
              ENDDO !no_nrlist (nr)
            ENDIF
          ELSE
            ERROR='>>Invalid subcommand'
            GOTO 9999
          ENDIF

        ELSE
          ERROR='>>Invalid file format'
          GOTO 9999
        ENDIF

C       Check that the data requested is actually present.
C       It should be enough just to check this for the OPEN case to save
C       checking for each set of time data.
        IF(SUBCOMMAND.EQ.'OPEN') THEN

          IF(NIYLIST(0).NE.0) THEN
            CALL ASSERT(YPDATA,'No yp variables in history file',
     '        ERROR,*9999)
            DO no_niylist=1,NIYLIST(0)
              niy=NIYLIST(no_niylist)
              IF(.NOT.INLIST(niy,NIYLIST_HIST(1,FILEID),
     '          NIYLIST_HIST(0,FILEID),IDUMMY)) THEN
                IEND=0
                CALL APPENDC(IEND,'yp variable ',ERROR)
                CALL APPENDI(IEND,niy,ERROR)
                CALL APPENDC(IEND,' not in history file',ERROR)
                CALL FLAG_ERROR(0,ERROR(:IEND))
                GOTO 9998
              ENDIF !.not.inlist
            ENDDO !niy
          ENDIF !NIYLIST(0).NE.0

          IF(NIQLIST(0).NE.0) THEN
            CALL ASSERT(YQDATA,'No yq variables in history file',
     '        ERROR,*9999)
            DO no_niqlist=1,NIQLIST(0)
              niq=NIQLIST(no_niqlist)
              IF(.NOT.INLIST(niq,NIQLIST_HIST(1,FILEID),
     '          NIQLIST_HIST(0,FILEID),IDUMMY)) THEN
                IEND=0
                CALL APPENDC(IEND,'yq variable ',ERROR)
                CALL APPENDI(IEND,niq,ERROR)
                CALL APPENDC(IEND,' not in history file',ERROR)
                CALL FLAG_ERROR(0,ERROR(:IEND))
                GOTO 9998
              ENDIF !.not.inlist
            ENDDO !niq
          ENDIF !NIQLIST(0).NE.0

          IF(NIQSLIST(0).NE.0) THEN
            CALL ASSERT(YQSDATA,'No yqs variables in history file',
     '        ERROR,*9999)
            DO no_niqslist=1,NIQSLIST(0)
              niqs=NIQSLIST(no_niqslist)
              IF(.NOT.INLIST(niqs,NIQSLIST_HIST(1,FILEID),
     '          NIQSLIST_HIST(0,FILEID),IDUMMY)) THEN
                IEND=0
                CALL APPENDC(IEND,'yqs variable ',ERROR)
                CALL APPENDI(IEND,niqs,ERROR)
                CALL APPENDC(IEND,' not in history file',ERROR)
                CALL FLAG_ERROR(0,ERROR(:IEND))
                GOTO 9998
              ENDIF !.not.inlist
            ENDDO !niqs
          ENDIF !NIQSLIST(0).NE.0

        ENDIF !SUBCOMMAND.EQ.'OPEN'

        ENDFILE=.FALSE.

      ELSE IF(COMMAND.EQ.'WRITE') THEN

        IF(FILEFORMAT.EQ.'ASCII') THEN

C*** Write ASCII Open
          IF(SUBCOMMAND.EQ.'OPEN') THEN

            OLDIOTYPE=IOTYPE
            IOTYPE=3 !Write
            VERSION(1)=3 !26/10/95
            CALL OPEN_SEQ_FILE(VERSION(1),FILEID,FILENAME,
     '        'hist','NEW',.FALSE.,ERROR,*9999)
            IOTYPE=OLDIOTYPE

            WRITE(FILEID,'('' Number of regions stored: '',I2)')
     '        NRLIST(0)
            WRITE(FILEID,'('' Regions are:'',10(1X,I2),'
     '        //':/(13X,10(1X,I2)))') (NRLIST(no_nrlist),no_nrlist=1,
     '        NRLIST(0))
            NRLIST_HIST(0,FILEID)=NRLIST(0)
            DO no_nrlist=1,NRLIST(0)
              NRLIST_HIST(no_nrlist,FILEID)=NRLIST(no_nrlist)
            ENDDO !no_nrlist (nr)
C            IF(YPDATA.AND.YQDATA) THEN
C              WRITE(FILEID,'('' The file contains both YP and YQ '
C     '          //'data'')')
C            ELSE IF(YPDATA) THEN
C              WRITE(FILEID,'('' The file contains only YP data'')')
C            ELSE IF(YQDATA) THEN
C              WRITE(FILEID,'('' The file contains only YQ data'')')
C            ENDIF
            IF(YPDATA) THEN
              IF(YQDATA.AND.YQSDATA) THEN
                STRING='YP,YQ,YQS'
              ELSE IF(YQDATA) THEN
                STRING='YP,YQ'
              ELSE
C LKC 2-JUN-1999 this does not appear correct
C                STRING='YP,YQS'
                STRING='YP'
              ENDIF
            ELSE IF(YQDATA) THEN
              IF(YQSDATA) THEN
                STRING='YQ,YQS'
              ELSE
                STRING='YQ'
              ENDIF
            ELSE
              STRING='YQS'
            ENDIF
            CALL STRING_TRIM(STRING,IBEG,IEND)
            WRITE(FILEID,'('' The file contains the variables: '',A)')
     '        STRING(IBEG:IEND)
            IF(YPDATA) THEN
              WRITE(FILEID,'(/,'' YP Data:'')')
              TOTALNY(FILEID)=0
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO nc=1,NCT(nr,nx)
                  TOTALNY(FILEID)=TOTALNY(FILEID)+NYNR(0,0,nc,nr)
                ENDDO !nc
              ENDDO !no_nrlist (nr)
              WRITE(FILEID,'(/,'' Total number of ny values stored: '','
     '          //'I10)') TOTALNY(FILEID)
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                WRITE(FILEID,'('' Region '',I2,'' : Number of '
     '            //'nc.s stored = '',I1)') nr,NCT(nr,nx)
                WRITE(FILEID,'(13X,''Number of ny.s for each nc :'','
     '            //'3(1X,I10),:/(41X,3(1X,I10)))')
     '            (NYNR(0,0,nc,nr),nc=1,NCT(nr,nx))
              ENDDO !no_nrlist (nr)
              WRITE(FILEID,'('' Number of niy indicies stored = '',I2)')
     '          NIYLIST(0)
              WRITE(FILEID,'('' Indicies are :'',20(1X,I2))')
     '          (NIYLIST(no_niylist),no_niylist=1,NIYLIST(0))
C cpb 17/7/96 Keep a copy of the niy indicies that the file was opened
C with. This list should be FILEID dependent.
              NIYLIST_HIST(0,FILEID)=NIYLIST(0)
              DO no_niylist=1,NIYLIST(0)
                NIYLIST_HIST(no_niylist,FILEID)= NIYLIST(no_niylist)
              ENDDO !no_niylist

              DO no_niylist=1,NIYLIST_HIST(0,FILEID)
                niy=NIYLIST_HIST(no_niylist,FILEID)
                CALL ASSERT(niy.LE.NIYM,'>>Increase NIYM',
     '            ERROR,*9999)
              ENDDO

              WRITE(FILEID,'('' Number of NPNY indicies stored = '','
     '          //'I2)') NUM_NPNY
              WRITE(FILEID,'('' NPNY values:'')')
              DO no_npny=0,NUM_NPNY-1
                IF(TOTALNY(FILEID).LE.10) THEN
                  WRITE(FILEID,'('' NPNY index '',I1,'' :'',10(1X,I5))')
     '              no_npny,(((NPNY(no_npny,NYNR(no_nynr,0,nc,
     '              NRLIST(no_nrlist)),0),
     '              no_nynr=1,NYNR(0,0,nc,NRLIST(no_nrlist))),
     '              nc=1,NCT(NRLIST(no_nrlist),nx)),
     '              no_nrlist=1,NRLIST(0))
                ELSE
                  WRITE(FILEID,'('' NPNY index '',I1,'' :'',10(1X,I5),'
     '              //':/(15X,10(1X,I5)))') no_npny,
     '              (((NPNY(no_npny,
     '              NYNR(no_nynr,0,nc,NRLIST(no_nrlist)),0),
     '              no_nynr=1,NYNR(0,0,nc,NRLIST(no_nrlist))),
     '              nc=1,NCT(NRLIST(no_nrlist),nx)),
     '              no_nrlist=1,NRLIST(0))
                ENDIF
              ENDDO !no_npny
            ENDIF
            IF(YQDATA) THEN
              WRITE(FILEID,'(/,'' YQ Data:'')')

              TOTALNQ(FILEID)=0
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                TOTALNQ(FILEID)=TOTALNQ(FILEID)+(NHQ(nr)*
     '            (NQR(2,nr)-NQR(1,nr)+1))
              ENDDO !no_nrlist (nr)
              CALL ASSERT(TOTALNQ(FILEID).EQ.NYQT(1,1,nx),
     '          '>>Inconsistent number of nyqs',ERROR,*9999)
              WRITE(FILEID,'(/,'' Total number of nyq values stored:'
     '          //' '',I10)') TOTALNQ(FILEID)

              nc=1 !for all grid problems at present
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                WRITE(FILEID,'('' Region '',I2,'' : Number of '
     '            //'nc values stored = '',I10))') nr,nc
                WRITE(FILEID,'('' Region '',I2,'' : Number of '
     '            //'nyq.s stored = '',I10))') nr,
     '            NHQ(nr)*(NQR(2,nr)-NQR(1,nr)+1)
              ENDDO !no_nrlist (nr)

              WRITE(FILEID,'('' Number of niq indicies stored = '',I2)')
     '          NIQLIST(0)
              WRITE(FILEID,'('' Indicies are :'',20(1X,I2))')
     '          (NIQLIST(no_niqlist),no_niqlist=1,NIQLIST(0))
              WRITE(FILEID,'('' Value of na index stored = '',I2)') na

              NIQLIST_HIST(0,FILEID)=NIQLIST(0)
              DO no_niqlist=1,NIQLIST(0)
                NIQLIST_HIST(no_niqlist,FILEID)=NIQLIST(no_niqlist)
              ENDDO !no_niqlist

              DO no_niqlist=1,NIQLIST_HIST(0,FILEID)
                niq=NIQLIST_HIST(no_niqlist,FILEID)
                CALL ASSERT(niq.LE.NIQM,'>>Increase NIQM',
     '            ERROR,*9999)
              ENDDO

              WRITE(FILEID,'('' Number of NQNY indicies '
     '          //'stored = '',I2)') NUM_NQNY

              WRITE(FILEID,'('' NQNY values:'')')

              DO no_nqny=1,NUM_NQNY
                IF(TOTALNQ(FILEID).LE.10) THEN
                  WRITE(FILEID,'('' NQNY index '',I1,'' :'',10(1X,I5))')
     '              no_nqny,((NQNY(no_nqny,NYQNR(no_nqnr,1,1,
     '              NRLIST(no_nrlist)),1),no_nqnr=1,NYQNR(0,0,1,
     '              NRLIST(no_nrlist))),no_nrlist=1,NRLIST(0))
                ELSE
                  WRITE(FILEID,'('' NQNY index '',I1,'' :'','
     '              //'10(1X,I5),:/(15X,10(1X,I5)))')
     '              no_nqny,((NQNY(no_nqny,
     '              NYQNR(no_nqnr,1,1,NRLIST(no_nrlist)),1),
     '              no_nqnr=1,NYQNR(0,0,1,NRLIST(no_nrlist))),
     '              no_nrlist=1,NRLIST(0))
                ENDIF
              ENDDO
            ENDIF
            IF(YQSDATA) THEN
              WRITE(FILEID,'(/,'' YQS Data:'')')

              TOTALNQS(FILEID)=0
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                TOTALNQS(FILEID)=TOTALNQS(FILEID)+(NQR(2,nr)-
     '            NQR(1,nr)+1)
              ENDDO !no_nrlist (nr)
              WRITE(FILEID,'(/,'' Total number of nq values stored = '''
     '          //',I10)') TOTALNQS(FILEID)

              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                WRITE(FILEID,'('' Region '',I2,'' : Start '
     '            //'nq number = '',I10))') nr,NQR(1,nr)
                WRITE(FILEID,'('' Region '',I2,'' : Stop '
     '            //'nq number = '',I10))') nr,NQR(2,nr)
              ENDDO !no_nrlist (nr)

              WRITE(FILEID,'('' Number of niqs indicies stored = '','
     '          //'I2)') NIQSLIST(0)
              WRITE(FILEID,'('' Indicies are :'',15(1X,I3),'
     '          //':/(15x,15(1X,I3)))') (NIQSLIST(no_niqslist),
     '          no_niqslist=1,NIQSLIST(0))

              NIQSLIST_HIST(0,FILEID)=NIQSLIST(0)
              DO no_niqslist=1,NIQSLIST(0)
                NIQSLIST_HIST(no_niqslist,FILEID)=NIQSLIST(no_niqslist)
              ENDDO !no_niqslist

              DO no_niqslist=1,NIQSLIST_HIST(0,FILEID)
                niqs=NIQSLIST_HIST(no_niqslist,FILEID)
                CALL ASSERT(niqs.LE.NIQSM,'>>Increase NIQSM',
     '            ERROR,*9999)
              ENDDO

            ENDIF

            NUMTIMEDATA=0

C*** Write ASCII Timedata
          ELSE IF(SUBCOMMAND.EQ.'TIME_DATA') THEN

            WRITE(FILEID,'(/,'' Time = '',D11.4)') TIME
C            IF(YPDATA.AND.YQDATA) THEN
C              WRITE(FILEID,'('' Time series contains both YP and YQ '
C     '          //'data'')')
C            ELSE IF(YPDATA) THEN
C              WRITE(FILEID,'('' Time series contains only YP data'')')
C            ELSE IF(YQDATA) THEN
C              WRITE(FILEID,'('' Time series contains only YQ data'')')
C            ENDIF
            IF(YPDATA) THEN
              IF(YQDATA.AND.YQSDATA) THEN
                STRING='YP,YQ,YQS'
              ELSE IF(YQDATA) THEN
                STRING='YP,YQ'
              ELSE
C 2-JUN-1999 LKC This is not correct?
C                STRING='YP,YQS'
                STRING='YP'
              ENDIF
            ELSE IF(YQDATA) THEN
              IF(YQSDATA) THEN
                STRING='YQ,YQS'
              ELSE
                STRING='YQ'
              ENDIF
            ELSE
              STRING='YQS'
            ENDIF
            CALL STRING_TRIM(STRING,IBEG,IEND)
            WRITE(FILEID,'('' The time series contains the '
     '        //'variables: '',A)') STRING(IBEG:IEND)
            IF(YPDATA) THEN
              WRITE(FILEID,'('' YP Data:'')')
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO nc=1,NCT(nr,nx)
                  DO no_nynr=1,NYNR(0,0,nc,nr)
                    ny=NYNR(no_nynr,0,nc,nr)
                    WRITE(FILEID,'('' YP('',I7,'',niy):'',4(1X,D12.5),'
     '                //':/(17X,4(1X,D12.5)))') ny,(YP(ny,
     '                NIYLIST(no_niylist)),no_niylist=1,NIYLIST(0))
                  ENDDO !no_nynr
                ENDDO !nc
              ENDDO !nr
            ENDIF
            IF(YQDATA) THEN
              WRITE(FILEID,'('' YQ Data:'')')
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO no_nqnr=1,NYQNR(0,0,1,nr)
                  nyq_index=NYQNR(no_nqnr,0,1,nr)
                  WRITE(FILEID,'('' YQ('',I7,'',niq):'',4(1X,D12.5),'
     '              //':/(17X,4(1X,D12.5)))') nyq_index,(YQ(nyq_index,
     '              NIQLIST(no_niqlist),na),no_niqlist=1,NIQLIST(0))
                ENDDO !no_nqnr (nyq)
              ENDDO !no_nrlist (nr)
            ENDIF
            IF(YQSDATA) THEN
              WRITE(FILEID,'('' YQS Data:'')')
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO nq=NQR(1,nr),NQR(2,nr)
                  WRITE(FILEID,'('' YQS(niqs,'',I7,''):'',4(1X,D12.5),'
     '              //':/(17X,4(1X,D12.5)))') nq,(YQS(NIQSLIST(
     '              no_niqslist),nq),no_niqslist=1,NIQSLIST(0))
                ENDDO !nq
              ENDDO !no_nrlist (nr)
            ENDIF

            NUMTIMEDATA=NUMTIMEDATA+1

          ELSE
            ERROR='>>Not implemented yet'
            GOTO 9999
          ENDIF

        ELSE IF(FILEFORMAT.EQ.'BINARY') THEN

C*** Write Binary Open
          IF(SUBCOMMAND.EQ.'OPEN') THEN

C***        Open the history file for writing
            FILETYPE=2 !Binary history file
            VERSION(1)=1 !Version 11/3/99
            VERSION(2)=2
            VERSION(3)=0
            NUMBERTAGS=1 !Just header tag for the moment.
            CALL OPEN_BIN_FILE(FILEID,FILETYPE,NUMBERTAGS,VERSION,
     '        'WRITE','his',FILENAME,ERROR,*9999)
C***        Write the header information tag. (Time tags added later).
            TAGINDEX=1
            TAGHEADER='History header information'
            NUMBERSUBTAGS=0
C            NUMTAGBYTES=(NRLIST(0)+1)*INTSIZE+2*INTSIZE
            NUMTAGBYTES=(NRLIST(0)+1)*INTSIZE+3*INTSIZE
            IF(YPDATA) THEN
              NUMTAGBYTES=NUMTAGBYTES+2*INTSIZE !total_ny and NUM_NPNY
C                                                not in the loop
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                NUMTAGBYTES=NUMTAGBYTES+INTSIZE !NCT
                DO nc=1,NCT(nr,nx)
                  NUMTAGBYTES=NUMTAGBYTES+INTSIZE+NYNR(0,0,nc,nr)*
     '              INTSIZE !NYNR(0,0,nc,nr)+ny list
                  DO no_nynr=1,NYNR(0,0,nc,nr)
                    NUMTAGBYTES=NUMTAGBYTES+NUM_NPNY*INTSIZE !NPNY
                  ENDDO !no_nynr
                ENDDO !nc
              ENDDO !no_nrlist (nr)
              NUMTAGBYTES=NUMTAGBYTES+(NIYLIST(0)+1)*INTSIZE !niy #+list
            ENDIF
            IF(YQDATA) THEN
              NUMTAGBYTES=NUMTAGBYTES+INTSIZE !total nyq
              NUMTAGBYTES=NUMTAGBYTES+NRLIST(0)*INTSIZE !nc
              NUMTAGBYTES=NUMTAGBYTES+NRLIST(0)*INTSIZE !nyq per region
              NUMTAGBYTES=NUMTAGBYTES+(NIQLIST(0)+1)*INTSIZE !niq
              NUMTAGBYTES=NUMTAGBYTES+INTSIZE !na
              NUMTAGBYTES=NUMTAGBYTES+INTSIZE !num_nqny
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO no_nqnr=1,NYQNR(0,1,1,nr)
                  NUMTAGBYTES=NUMTAGBYTES+INTSIZE*NUM_NQNY !nqny indices
                ENDDO !no_nqnr (nyq)
              ENDDO !no_nrlist (nr)
            ENDIF
            IF(YQSDATA) THEN
              NUMTAGBYTES=NUMTAGBYTES+INTSIZE !total nq
              NUMTAGBYTES=NUMTAGBYTES+NRLIST(0)*INTSIZE !nq per region
              NUMTAGBYTES=NUMTAGBYTES+NRLIST(0)*INTSIZE !start nq
              NUMTAGBYTES=NUMTAGBYTES+NRLIST(0)*INTSIZE !stop nq
              NUMTAGBYTES=NUMTAGBYTES+(NIQSLIST(0)+1)*INTSIZE !niqs
            ENDIF

            CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)

C***        Write the number of regions and the region list
            CALL BINWRITEFILE(FILEID,INTTYPE,NRLIST(0)+1,NRLIST,
     '        REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '        *9999)
            NRLIST_HIST(0,FILEID)=NRLIST(0)
            DO no_nrlist=1,NRLIST(0)
              NRLIST_HIST(no_nrlist,FILEID)=NRLIST(no_nrlist)
            ENDDO !no_nrlist (nr)
C***        Write out if YP data present
            IF(YPDATA) THEN
              INTDATA(1)=1
            ELSE
              INTDATA(1)=0
            ENDIF
            CALL BINWRITEFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***        Write out if YQ data present
            IF(YQDATA) THEN
              INTDATA(1)=1
            ELSE
              INTDATA(1)=0
            ENDIF
            CALL BINWRITEFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***        Write out if YQS data present
            IF(YQSDATA) THEN
              INTDATA(1)=1
            ELSE
              INTDATA(1)=0
            ENDIF
            CALL BINWRITEFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            IF(YPDATA) THEN
C***          Write out the total number of ny's stored
              TOTALNY(FILEID)=0
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO nc=1,NCT(nr,nx)
                  TOTALNY(FILEID)=TOTALNY(FILEID)+NYNR(0,0,nc,nr)
                ENDDO !nc
              ENDDO !no_nrlist (nr)
              CALL BINWRITEFILE(FILEID,INTTYPE,1,TOTALNY(FILEID),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)
C***          Write out the total number of nc's for each region
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                CALL BINWRITEFILE(FILEID,INTTYPE,1,NCT(nr,nx),REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              ENDDO !no_nrlist (nr)
C***          Write out the list of ny's for each nc for each region
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO nc=1,NCT(nr,nx)
C***              Write out the # of ny's for each nc for each region
                  CALL BINWRITEFILE(FILEID,INTTYPE,1,NYNR(0,0,nc,nr),
     '              REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '              ERROR,*9999)
C***              Write out the list of ny's for each nc for each region
                  CALL BINWRITEFILE(FILEID,INTTYPE,NYNR(0,0,nc,nr),
     '              NYNR(1,0,nc,nr),REAL4DATA,REAL8DATA,CHARDATA,
     '              LOGDATA,SINTDATA,ERROR,*9999)
                ENDDO !nc
              ENDDO !no_nrlist (nr)
C***          Write out the number of NIY indices stored in the file
              CALL BINWRITEFILE(FILEID,INTTYPE,1,NIYLIST(0),REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***          Write out the number of NIY index list stored in the file
              CALL BINWRITEFILE(FILEID,INTTYPE,NIYLIST(0),NIYLIST(1),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)
C cpb 17/7/96 Keep a copy of the niy indicies that the file was opened
C with. This list should be FILEID dependent.
              NIYLIST_HIST(0,FILEID)=NIYLIST(0)
              DO no_niylist=1,NIYLIST(0)
                NIYLIST_HIST(no_niylist,FILEID)= NIYLIST(no_niylist)
              ENDDO !no_niylist

              DO no_niylist=1,NIYLIST_HIST(0,FILEID)
                niy=NIYLIST_HIST(no_niylist,FILEID)
                CALL ASSERT(niy.LE.NIYM,'>>Increase NIYM',
     '            ERROR,*9999)
              ENDDO

C***          Write out the number of NPNY values (second dimension of
C***          NPNY).
              INTDATA(1)=NUM_NPNY
              CALL BINWRITEFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***          Write out NPNY
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO nc=1,NCT(nr,nx)
                  DO no_nynr=1,NYNR(0,0,nc,nr)
                    ny=NYNR(no_nynr,0,nc,nr)
                    CALL BINWRITEFILE(FILEID,INTTYPE,NUM_NPNY,
     '                NPNY(0,ny,0),REAL4DATA,REAL8DATA,CHARDATA,
     '                LOGDATA,SINTDATA,ERROR,*9999)
                  ENDDO !no_nynr (ny)
                ENDDO !nc
              ENDDO !no_nrlist (nr)
            ENDIF
            IF(YQDATA) THEN
              TOTALNQ(FILEID)=0
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                TOTALNQ(FILEID)=TOTALNQ(FILEID)+(NHQ(nr)*
     '            (NQR(2,nr)-NQR(1,nr)+1))
              ENDDO !no_nrlist (nr)
              CALL ASSERT(TOTALNQ(FILEID).EQ.NYQT(1,1,nx),
     '          '>>Inconsistent number of nyqs',ERROR,*9999)

C***          Write out the total number of nyq's stored
              CALL BINWRITEFILE(FILEID,INTTYPE,1,TOTALNQ(FILEID),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)

              nc=1 !for all grid problems at present
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
C***            Number of nc's per region
                na1(1)=nc
                CALL BINWRITEFILE(FILEID,INTTYPE,1,na1,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              ENDDO !no_nrlist (nr)

              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
C***            Number of nyq's per region
                na1(1)=NHQ(nr)*(NQR(2,nr)-NQR(1,nr)+1)
                CALL BINWRITEFILE(FILEID,INTTYPE,1,na1,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              ENDDO !no_nrlist (nr)

C***          Write out the number of NIQ indices stored in the file
              CALL BINWRITEFILE(FILEID,INTTYPE,1,NIQLIST(0),REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***          Write out the number of NIQ index list stored in the file
              CALL BINWRITEFILE(FILEID,INTTYPE,NIQLIST(0),NIQLIST(1),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)
C***          Write out the value of the na index
              na1(1)=na
              CALL BINWRITEFILE(FILEID,INTTYPE,1,na1,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

              NIQLIST_HIST(0,FILEID)=NIQLIST(0)
              DO no_niqlist=1,NIQLIST(0)
                NIQLIST_HIST(no_niqlist,FILEID)=NIQLIST(no_niqlist)
              ENDDO !no_niqlist

              DO no_niqlist=1,NIQLIST_HIST(0,FILEID)
                niq=NIQLIST_HIST(no_niqlist,FILEID)
                CALL ASSERT(niq.LE.NIQM,'>>Increase NIQM',
     '            ERROR,*9999)
              ENDDO

C***          Write out the number of NQNY indicies for each nyq
              INTDATA(1)=NUM_NQNY
              CALL BINWRITEFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C***          Write out NPNY
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO no_nqnr=1,NYQNR(0,1,1,nr)
                  ny=NYQNR(no_nqnr,1,1,nr)
                  DO nq=1,NUM_NQNY
                    CALL BINWRITEFILE(FILEID,INTTYPE,1,NQNY(nq,ny,1),
     '                REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                ERROR,*9999)
                  ENDDO !nq
                ENDDO !no_nqnr (nyq)
              ENDDO !no_nrlist (nr)

            ENDIF
            IF(YQSDATA) THEN
              TOTALNQS(FILEID)=0
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                TOTALNQS(FILEID)=TOTALNQS(FILEID)+NQR(2,nr)-NQR(1,nr)+1
              ENDDO !no_nrlist (nr)

C***          Write out the total number of nq's stored
              CALL BINWRITEFILE(FILEID,INTTYPE,1,TOTALNQS(FILEID),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)

              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
C***            start nq per region
                NQRS_HIST(1,nr,FILEID)=NQR(1,nr)
                CALL BINWRITEFILE(FILEID,INTTYPE,1,NQR(1,nr),
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
              ENDDO !no_nrlist (nr)

              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
C***            stop nq per region
                NQRS_HIST(2,nr,FILEID)=NQR(2,nr)
                CALL BINWRITEFILE(FILEID,INTTYPE,1,NQR(2,nr),REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              ENDDO !no_nrlist (nr)

C***          Write out the number of NIQS indices stored in the file
              CALL BINWRITEFILE(FILEID,INTTYPE,1,NIQSLIST(0),REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***          Write out the NIQS index list stored in the file
              CALL BINWRITEFILE(FILEID,INTTYPE,NIQSLIST(0),NIQSLIST(1),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)

              NIQSLIST_HIST(0,FILEID)=NIQSLIST(0)
              DO no_niqslist=1,NIQSLIST(0)
                NIQSLIST_HIST(no_niqslist,FILEID)=NIQSLIST(no_niqslist)
              ENDDO !no_niqslist

              DO no_niqslist=1,NIQSLIST_HIST(0,FILEID)
                niqs=NIQSLIST_HIST(no_niqslist,FILEID)
                CALL ASSERT(niqs.LE.NIQSM,'>>Increase NIQSM',
     '            ERROR,*9999)
              ENDDO

            ENDIF
            NUMTIMEDATA=0

C*** Write Binary Append
          ELSE IF(SUBCOMMAND.EQ.'APPEND') THEN

            ERROR='>>Not implemented yet'
            GOTO 9999

C*** Write Binary Timedata
          ELSE IF(SUBCOMMAND.EQ.'TIME_DATA') THEN

C***        Set the number of tags in the file
            NUMTIMEDATA=NUMTIMEDATA+1 !Adjust the number of times
            NUMBERTAGS=NUMTIMEDATA+1 !Adjust the number of tags
            CALL BINRESETNUMTAG(FILEID,NUMBERTAGS,ERROR,*9999)
C***        Write the tag header information
            TAGINDEX=2
            TAGHEADER='Time series data'
            NUMBERSUBTAGS=0
            NUMTAGBYTES=DPSIZE+3*INTSIZE

            IF(YPDATA) NUMTAGBYTES=NUMTAGBYTES+TOTALNY(FILEID)*
     '        NIYLIST_HIST(0,FILEID)*DPSIZE
            IF(YQDATA) NUMTAGBYTES=NUMTAGBYTES+TOTALNQ(FILEID)*
     '        NIQLIST_HIST(0,FILEID)*DPSIZE
            IF(YQSDATA) NUMTAGBYTES=NUMTAGBYTES+TOTALNQS(FILEID)*
     '        NIQSLIST_HIST(0,FILEID)*DPSIZE
            CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
C***        Write the time
            TIME2(1)=TIME
            CALL BINWRITEFILE(FILEID,DPTYPE,1,INTDATA,REAL4DATA,
     '        TIME2,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            TIME=TIME2(1)
C***        Write whether YP data present or not
            IF(YPDATA) THEN
              INTDATA(1)=1
            ELSE
              INTDATA(1)=0
            ENDIF
            CALL BINWRITEFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***        Write whether YQ data present or not
            IF(YQDATA) THEN
              INTDATA(1)=1
            ELSE
              INTDATA(1)=0
            ENDIF
            CALL BINWRITEFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***        Write whether YQS data present or not
            IF(YQSDATA) THEN
              INTDATA(1)=1
            ELSE
              INTDATA(1)=0
            ENDIF
            CALL BINWRITEFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            IF(YPDATA) THEN
C***          Write the YP time data
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO nc=1,NCT(nr,nx)
                  DO no_nynr=1,NYNR(0,0,nc,nr)
                    ny=NYNR(no_nynr,0,nc,nr)
                    DO no_niylist=1,NIYLIST(0)
                      niy=NIYLIST(no_niylist)
                      CALL BINWRITEFILE(FILEID,DPTYPE,1,INTDATA,
     '                  REAL4DATA,YP(ny,niy),CHARDATA,LOGDATA,
     '                  SINTDATA,ERROR,*9999)
                    ENDDO !no_niylist (niy)
                  ENDDO !no_nynr (ny)
                ENDDO !nc
              ENDDO !no_nrlist (nr)
            ENDIF
            IF(YQDATA) THEN
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO no_nqnr=1,NYQNR(0,0,1,nr)
                  nyq_index=NYQNR(no_nqnr,0,1,nr)
                  DO no_niqlist=1,NIQLIST(0)
                    niq=NIQLIST(no_niqlist)
                    CALL BINWRITEFILE(FILEID,DPTYPE,1,INTDATA,REAL4DATA,
     '                YQ(NYQNR(nyq_index,0,1,nr),niq,na),CHARDATA,
     '                LOGDATA,SINTDATA,ERROR,*9999)
                  ENDDO !no_niqlist (niq)
                ENDDO !no_nqnr (nyq)
              ENDDO !no_nrlist (nr)
            ENDIF
            IF(YQSDATA) THEN
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO no_niqslist=1,NIQSLIST(0)
                    niqs=NIQSLIST(no_niqslist)
                    CALL BINWRITEFILE(FILEID,DPTYPE,1,INTDATA,
     '                REAL4DATA,YQS(niqs,nq),CHARDATA,LOGDATA,
     '                SINTDATA,ERROR,*9999)
                  ENDDO !no_niqslist (niqs)
                ENDDO !nq
              ENDDO !no_nrlist (nr)
            ENDIF
          ELSE
            ERROR='>>Invalid subcommand'
            GOTO 9999
          ENDIF

        ELSE
          ERROR='>>Invalid file format'
          GOTO 9999
        ENDIF

      ELSE IF(COMMAND.EQ.'CLOSE') THEN

        IF(FILEFORMAT.EQ.'ASCII') THEN

          CALL CLOSEF(FILEID,ERROR,*9999)

        ELSE IF(FILEFORMAT.EQ.'BINARY') THEN

          CALL BINCLOSEFILE(FILEID,ERROR,*9999)

        ELSE
          ERROR='>>Invalid file format'
          GOTO 9999
        ENDIF

      ELSE
        ERROR='>>Invalid command'
        GOTO 9999
      ENDIF

      CALL EXITS('IOHIST')
      RETURN
 9997 ENDFILE=.TRUE.
      CALL EXITS('IOHIST')
      RETURN
 9998 ERROR=' '
 9999 CALL ERRORS('IOHIST',ERROR)
      CALL EXITS('IOHIST')
      IF(ISBINFILEOPEN(FILEID)) CALL BINARYCLOSEFILE(FILEID,ERR,CERROR)
      RETURN 1
      END


