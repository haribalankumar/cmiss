      SUBROUTINE IPLEAD(ERROR,*)

C#### Subroutine: IPLEAD
C###  Description:
C###    IPLEAD inputs electrocardiographic leads.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'lead00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,ICHAR,IEND,INFO,nlead,nleadelec,NOQUES
      CHARACTER CHAR2*2,CHAR3*2
      LOGICAL FILEIP

      CALL ENTERS('IPLEAD',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='(/$,'' The number of leads is [1]: '',I2)'
      IF(IOTYPE.EQ.3) IDATA(1)=NUMLEADS
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,NUMLEADMX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NUMLEADS=IDATA(1)

      DO nlead=1,NUMLEADS
        WRITE(CHAR2,'(I2)') nlead
        CDEFLT(1)='Lead'//CHAR2
        CALL STRING_TRIM(CDEFLT(1),IBEG,IEND)
        FORMAT='(/$,'' Enter the title for lead '//CHAR2
     '    //' ['//CDEFLT(1)(IBEG:IEND)//']: '',A15)'
        IF(IOTYPE.EQ.3) CDATA(1)=LEADTITLE(nlead)
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,15,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) LEADTITLE(nlead)=CDATA(1)(1:15)
        FORMAT='($,'' Enter the number of electrodes in lead '//CHAR2
     '      //' [2]: '',I2)'
        IDEFLT(1)=2
        IF(IOTYPE.EQ.3) IDATA(1)=LEADELECS(0,nlead)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NUMLEADELECMX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) LEADELECS(0,nlead)=IDATA(1)
        WRITE(CHAR3,'(I2)') LEADELECS(0,nlead)
        FORMAT='($,'' Enter the '//CHAR3//' electrode numbers in '
     '    //'lead '//CHAR2//': '','//CHAR3//'I5)'
        IF(IOTYPE.EQ.3) THEN
          DO nleadelec=1,LEADELECS(0,nlead)
            IDATA(nleadelec)=LEADELECS(nleadelec,nlead)
          ENDDO !nleadelec
        ELSE
          DO nleadelec=1,LEADELECS(0,nlead)
            IDEFLT(nleadelec)=1
          ENDDO !nleadelec
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    LEADELECS(0,nlead),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '    IDATA,IDEFLT,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO nleadelec=1,LEADELECS(0,nlead)
            LEADELECS(nleadelec,nlead)=IDATA(nleadelec)
          ENDDO !nleadelec
        ENDIF
        FORMAT='($,'' Enter the additive potential constant in '
     '    //'lead '//CHAR2//' [0.0]: '',D12.4)'
        IF(IOTYPE.EQ.3) THEN
            RDATA(1)=LEADCOUP(0,nlead)
        ELSE
            RDEFLT(1)=0.0d0
        ENDIF
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) LEADCOUP(0,nlead)=RDATA(1)
        FORMAT='($,'' Enter the '//CHAR3//' electrode coefficents in '
     '    //'lead '//CHAR2//': '','//CHAR3//'F7.3)'
        IF(IOTYPE.EQ.3) THEN
          DO nleadelec=1,LEADELECS(0,nlead)
            RDATA(nleadelec)=LEADCOUP(nleadelec,nlead)
          ENDDO !nleadelec
        ELSE
          DO nleadelec=1,LEADELECS(0,nlead)
            RDEFLT(nleadelec)=1.0d0
          ENDDO !nleadelec
        ENDIF
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    LEADELECS(0,nlead),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '    IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '    RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO nleadelec=1,LEADELECS(0,nlead)
            LEADCOUP(nleadelec,nlead)=RDATA(nleadelec)
          ENDDO !nleadelec
        ENDIF
      ENDDO !nlead

C LKC 24-JUL-98
C
C      FORMAT='(/$,'' Enter the filename for the electrode signal '
C     '  //'file [current]: '',A30)'
C      CDEFLT(1)=FILE00
C      IF(IOTYPE.EQ.3) CDATA(1)=SIGFNAME
C      CALL GINOUT(IOTYPE,2,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '  ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
C     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C      IF(IOTYPE.NE.3) SIGFNAME=CDATA(1)(1:30)

      CALL EXITS('IPLEAD')
      RETURN
 9999 CALL ERRORS('IPLEAD',ERROR)
      CALL EXITS('IPLEAD')
      RETURN 1
      END


