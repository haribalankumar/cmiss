      SUBROUTINE DEGAUS(NBJ,NRLIST,YG,STRING,ERROR,*)

C#### Subroutine: DEGAUS
C###  Description:
C###    DEGAUS defines Gauss point positions.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NRLIST(0:NRM)
      REAL*8 YG(NIYGM,NGM,NEM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IFROMC,IPFILE,
     '  N3CO,no_nrlist,nr,
     '  NT_YG
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,CBBREV,FILIO,FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DEGAUS',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define gauss;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Gauss point values are read from or written to a file
C###    FILENAME.ippres.  This is for reading or writing a Gauss point
C###    pressure file for VSAero.
C###  Parameter:      <number VARIABLE#[1]>
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<number VARIABLE#[1]>'
C       OP_STRING(2)=BLANK(1:15)//'<pressure>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
C       OP_STRING(4)=BLANK(1:15)//'<scale FACTOR#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEGAUS',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 24-Jan-1990
        CALL PARSE_QUALIFIERS('CDLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
C         IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
C           FACTOR=RFROMC(CO(N3CO+1))
C         ELSE
C           FACTOR=1.0D0
C         ENDIF
          IF(CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
            NT_YG=IFROMC(CO(N3CO+1))
          ELSE
            NT_YG=1
          ENDIF

        ENDIF

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'gaus',STATUS,
     '        ERR,ERROR,*9999)
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              CALL IPGAUS(NBJ,nr,NT_YG,YG,ERROR,*9999)
            ENDDO
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF
      ENDIF

C MHT 04-05-01 no path to 9998
C 9998 CALL EXITS('DEGAUS')
C      RETURN

      CALL EXITS('DEGAUS')
      RETURN
 9999 CALL ERRORS('DEGAUS',ERROR)
      CALL EXITS('DEGAUS')
      RETURN 1
      END


