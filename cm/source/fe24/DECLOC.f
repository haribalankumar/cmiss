      SUBROUTINE DECLOC(ISEG,ISCLOC,CSEG,STRING,ERROR,*)

C#### Subroutine: DECLOC
C###  Description:
C###    DECLOC defines clock.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cloc00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISCLOC(NWM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,ICHAR,IEND,IEND1,IEND2,
     '  iw,IWK(6),INDEX,INDEX_OLD,
     '  INFO,INSTAT,ISTEMP,NOQUES,NTIW
      REAL*8 XWC1,XWC2,YWC1,YWC2
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL CALCU,FILEIP,FILIO,GENER,MOUSE

      CALL ENTERS('DECLOC',*9999)
      ICHAR=999
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define clock;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define position and size of clock face. Clock parameters can be
C###    read from or written to the file FILENAME.ipcloc, with $current
C###    giving the current default file.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C#### Command: FEM define clock;m
C###  Parameter:      <on WS>[1]
C###  Description:
C###    Define position and size of clock face with mouse on specified
C###    workstation.

        OP_STRING(1)=STRING(1:IEND)//';m'
        OP_STRING(2)=BLANK(1:15)//'<on WS>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DECLOC',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS('DLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(MOUSE) CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        FILEIP=.FALSE.
        NOQUES=0

        IF(FILIO) THEN
          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipcloc',STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
          FORMAT='($,'' Enter coords of clock centre [0.0,0.0]: '','
     '      //'2E11.4)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=X_CLOCK(1)
            RDATA(2)=X_CLOCK(2)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,
     '      INFO,ERROR,*9999)
c         CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,2,RDATA,RZERO,
c    '      -RMAX,RMAX,INFO,ERROR,*9999)
          IF(iotype.NE.3) THEN
            X_CLOCK(1)=RDATA(1)
            X_CLOCK(2)=RDATA(2)
          ENDIF
          FORMAT='($,'' Enter clock radius [0.1]: '',E11.4)'
          RDEFLT(1)=0.1D0
          IF(IOTYPE.EQ.3) RDATA(1)=R_CLOCK
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,
     '      INFO,ERROR,*9999)
c         CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,
c    '      0.0D0,RMAX,INFO,ERROR,*9999)
          IF(iotype.NE.3) R_CLOCK=RDATA(1)
          CALL CLOSEF(IFILE,ERROR,*9999)

        ELSE IF(MOUSE) THEN
          iw=IWK(1)
          CALL ACWK(iw,0,ERROR,*9999)
          CALL OPEN_SEGMENT(ISTEMP,ISEG,iw,'Temp',INDEX,INDEX_OLD,
     '      0,1,CSEG,ERROR,*9999)
          INSTAT=1
          WRITE(OP_STRING,'('' >>Locate centre of clock'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL LOCATOR(INSTAT,0.0D0,XWC1,0.0D0,YWC1,
     '      ERROR,*9999)
          IF(INSTAT.EQ.1) THEN
            X_CLOCK(1)=XWC1
            X_CLOCK(2)=YWC1
            WRITE(OP_STRING,'('' >>Locate radius of clock'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(INSTAT,XWC1,XWC2,YWC1,YWC2,
     '        ERROR,*9999)
            R_CLOCK=DSQRT((XWC2-XWC1)**2+(YWC2-YWC1)**2)
          ENDIF
          CALL CLOSE_SEGMENT(ISTEMP,iw,ERROR,*9999)
          CALL SGCLOC(INDEX,ISCLOC(iw),ISEG,iw,CSEG,0.0d0,ERROR,*9999)
          CALL DAWK(iw,0,ERROR,*9999)

        ENDIF
      ENDIF

      CALL EXITS('DECLOC')
      RETURN
 9999 CALL ERRORS('DECLOC',ERROR)
      CALL EXITS('DECLOC')
      RETURN 1
      END


