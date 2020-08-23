      SUBROUTINE DECELL(CELL_ICQS_SPATIAL,CELL_ICQS_VALUE,
     '  CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,IBT,IDO,INP,NEELEM,NELIST,
     '  NENQ,NPNE,NPNODE,NQET,NQNE,NQS,NQXI,NRLIST,NXLIST,
     '  CE,CELL_RCQS_VALUE,CELL_YQS_VALUE,CELL_ICQS_NAMES,
     '  CELL_RCQS_NAMES,CELL_YQS_NAMES,CP,CQ,XE,STRING,ERROR,*)

C#### Subroutine: DECELL
C###  Description:
C###    DECELL reads in cell model parameters from a *.ipcell file
C###    which has been created using the CMGUI cellular modelling
C###    interface
C**** Created by Martin Buist, July 1998

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER CELL_ICQS_SPATIAL(NQIM,NQVM),CELL_ICQS_VALUE(NQIM,NQVM),
     '  CELL_RCQS_SPATIAL(NQRM,NQVM),CELL_YQS_SPATIAL(NIQSM,NQVM),
     '  IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NENQ(0:8,NQM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),
     '  NQNE(NEQM,NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),
     '  NXLIST(0:NXM)
      REAL*8 CE(NMM,NEM,NXM),CELL_RCQS_VALUE(NQRM,NQVM),
     '  CELL_YQS_VALUE(NIQSM,NQVM),CP(NMM,NPM,NXM),CQ(NMM,NQM,NXM),
     '  XE(NSM,NJM)
      CHARACTER CELL_ICQS_NAMES(NQIM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH),
     '  ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,
     '  IEND2,IPFILE,nr,nx,nxc
      LOGICAL ALL_REGIONS,CALCU,FILIO,GENER,MOUSE
      CHARACTER FILE*(MXCH),STATUS*3

      CALL ENTERS('DECELL',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C--------------------------------------------------------------------

C#### Command: FEM define cell;l/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    The material parameters are read from the file FILENAME.ipcell
C###    in the directory specified by PATH with $current specifing the
C###    current default file.

        OP_STRING(1)=STRING(1:IEND)//';l/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C--------------------------------------------------------------------

C#### Command: FEM define cell;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Parameter:      <region #[1]>
C###    Specify the region number to use.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    The material parameters are read from the file FILENAME.ipcell
C###    in the directory specified by PATH with $current specifing the
C###    current default file.
C###    Note that prompted input is for simple models only.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C--------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DECELL',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRID.EQ.1,'Must set USE_GRID to 1',ERROR,*9999)
        CALL ASSERT(USE_CELL.EQ.1,'Must set USE_CELL to 1',ERROR,*9999)
        CALL ASSERT(CALL_EQUA,'>>Must define equations first',
     '    ERROR,*9999)
C        CALL PARSE_QUALIFIERS('LR',noco,1,CO,COQU,CALCU,FILIO,GENER,
C     '    MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        nr=NRLIST(1)
        IF((ITYP3(nr,nx).LE.4).AND.(ITYP19(nr,nx).EQ.1)) THEN
          IPFILE=1
          CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,CALCU,FILIO,
     '      GENER,MOUSE,STATUS,ERROR,*1)
          IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

          IF(FILIO) THEN
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'cell',STATUS,
     '        ERR,ERROR,*9999)

            CALL IPCELL_PROMPT(IBT,IDO,INP,NEELEM,NELIST,NENQ,NPNE,
     '        NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,CE(1,1,nx),CP(1,1,nx),
     '        CQ(1,1,nx),XE,ERROR,*9999)

            CALL CLOSEF(IFILE,ERROR,*9999)
          ENDIF !filio
          CALL_CELL=.TRUE.
        ELSE IF((ITYP19(nr,nx).EQ.1).AND.(ITYP3(nr,
     '      nx).EQ.10).AND.(KTYP33.EQ.6)) THEN
          ! We have a user defined electrical model, using a model
          ! provided by a cellML description.
          IPFILE=1
          CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,CALCU,FILIO,
     '      GENER,MOUSE,STATUS,ERROR,*1)
          IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

          IF(FILIO) THEN
            IF(IOTYPE.EQ.3) THEN
              CALL OPEN_SEQ_FILE(IPFILE,IFILE,FILE,'cell',STATUS,
     '          ALL_REGIONS,ERROR,*9999)
              CALL IPCELL_WRITE(CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '          CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,nr,CELL_RCQS_VALUE,
     '          CELL_YQS_VALUE,CELL_ICQS_NAMES,CELL_RCQS_NAMES,
     '          CELL_YQS_NAMES,CQ(1,1,nx),ERROR,*9999)
              CALL CLOSEF(IFILE,ERROR,*9999)
            ELSE
              CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'cell',STATUS,
     '          ERR,ERROR,*9999)
              CALL IPCELL_CELLML(CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '          CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,nr,CELL_RCQS_VALUE,
     '          CELL_YQS_VALUE,CELL_ICQS_NAMES,CELL_RCQS_NAMES,
     '          CELL_YQS_NAMES,CQ(1,1,nx),ERROR,*9999)
              CALL CLOSEF(IFILE,ERROR,*9999)
              CALL_CELL=.TRUE.
            ENDIF
          ENDIF !filio
        ELSE
          CALL PARSE_QUALIFIERS('LRW',noco,1,CO,COQU,CALCU,FILIO,GENER,
     '      MOUSE,STATUS,ERROR,*1)
          IPFILE=1
          IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          IF(FILIO) THEN
            CALL OPEN_SEQ_FILE(IPFILE,IFILE,FILE,'cell',STATUS,
     '        ALL_REGIONS,ERROR,*9999)
            IF(IOTYPE.EQ.3) THEN
              CALL IPCELL_WRITE(CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '          CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,nr,CELL_RCQS_VALUE,
     '          CELL_YQS_VALUE,CELL_ICQS_NAMES,CELL_RCQS_NAMES,
     '          CELL_YQS_NAMES,CQ(1,1,nx),ERROR,*9999)
            ELSE
              CALL IPCELL_READ(CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '          CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,nr,CELL_RCQS_VALUE,
     '          CELL_YQS_VALUE,CELL_ICQS_NAMES,CELL_RCQS_NAMES,
     '          CELL_YQS_NAMES,CQ(1,1,nx),ERROR,*9999)
              CALL_CELL=.TRUE.
            ENDIF
            CALL CLOSEF(IFILE,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('DECELL')
      RETURN
 9999 CALL ERRORS('DECELL',ERROR)
      CALL EXITS('DECELL')
      RETURN 1
      END


