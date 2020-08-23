      SUBROUTINE IPGROW(ERROR,*)

C#### Subroutine: IPGROW
C###  Description:
C###    IPGROW defines growth parameters.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'grow00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp60.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES
      CHARACTER CHAR10*10,CHAR9*9
      LOGICAL FILEIP

      CALL ENTERS('IPGROW',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='('' Specify type of growth law [1]:'''//
     '  '/''   (1) Density = k* strain energy density'''//
     '  '/''   (2) Density = constant + k* strain energy density'''//
     '  '/''   (3) Carter bone law'''//
     '  '/''   (4) Huiskes bone law'''//
     '  '/''   (5) Tree growth #1'''//
     '  '/''   (6) Tree growth #2'''//
     '  '/''   (7) Tree growth #3'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP60
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP60=IDATA(1)

      IF(KTYP60.EQ.1) THEN
        GROW1=0.0D0
        FORMAT='($,'' Specify growth proportionality constant [0]: '','
     '    //'D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=GROW2
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) GROW2=RDATA(1)

      ELSE IF(KTYP60.EQ.2) THEN
        FORMAT='($,'' Specify growth density constant [0]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=GROW1
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) GROW1=RDATA(1)

        FORMAT='($,'' Specify growth proportionality constant [0]: '','
     '    //'D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=GROW2
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) GROW2=RDATA(1)

      ELSE IF(KTYP60.EQ.3) THEN !Carter growth law
        FORMAT='($,'' Specify # loading cycles [1]: '',I10)'
        IF(IOTYPE.EQ.3) IDATA(1)=NT_LOADING_CYCLES
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NT_LOADING_CYCLES=IDATA(1)

        RDEFLT(1)=100.0D0
        CHAR9='100Kg/m^3'
        FORMAT='($,'' Specify minimum density ['//CHAR9//']: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=DENSITY_MIN
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) DENSITY_MIN=RDATA(1)

        RDEFLT(1)=1800.0D0
        CHAR10='1800Kg/m^3'
        FORMAT='($,'' Specify maximum density ['//CHAR10//']: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=DENSITY_MAX
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) DENSITY_MAX=RDATA(1)

      ELSE IF(KTYP60.EQ.4) THEN !Huiskes growth law

      ELSE IF(KTYP60.EQ.5) THEN !Tree growth law #1
        FORMAT='($,'' Specify growth length change constant [0]: '','
     '    //'D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=GROW1
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) GROW1=RDATA(1)

        FORMAT='($,'' Specify growth proportionality constant [0]: '','
     '    //'D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=GROW2
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) GROW2=RDATA(1)

      ELSE IF(KTYP60.EQ.6) THEN !Tree growth law #2

      ELSE IF(KTYP60.EQ.7) THEN !Tree growth law #3
      ENDIF

      CALL EXITS('IPGROW')
      RETURN
 9999 CALL ERRORS('IPGROW',ERROR)
      CALL EXITS('IPGROW')
      RETURN 1
      END


