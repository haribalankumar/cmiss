      SUBROUTINE DEREFI(IBT,NBJ,NEELEM,NELIST,NKJE,NPF,NPNE,NRLIST,NVJE,
     '  NEERR,PG,RG,SE,WG,XA,XE,XG,XP,STRING,ERROR,*)

C#### Subroutine: DEREFI
C###  Description:
C###    DEREFI defines refinement parameters.
C**** Created by Carey Stevens 1 October 1997

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 NEERR(NEM,3),PG(NSM,NUM,NGM,NBM),
     '  RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DEREFI',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define refinement;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Parameter:       <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Description:
C###    Define refinement parameters.
C###    Parameters can be read from or written to the file
C###    FILENAME.iprefi, with $current specifing the current
C###    default file.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:13)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEREFI',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 30-Sept-97
        CALL PARSE_QUALIFIERS('LPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          IF(NRLIST(0).GT.1) THEN
            ALL_REGIONS=.TRUE.
          ELSE
            ALL_REGIONS=.FALSE.
          ENDIF
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'refi',STATUS,
     '        ERR,ERROR,*9999)
            CALL IPREFI(IBT,NBJ,NEELEM,NELIST,NKJE,NPF,NPNE,NRLIST,NVJE,
     '        NEERR,PG,RG,SE,WG,XA,XE,XG,XP,ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF
      ENDIF
      CALL_DEREFI=.TRUE.

      CALL EXITS('DEREFI')
      RETURN
 9999 CALL ERRORS('DEREFI',ERROR)
      CALL EXITS('DEREFI')
      RETURN 1
      END


