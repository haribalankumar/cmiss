      SUBROUTINE DECORN(INP,NBJ,NKJE,NPF,NPNE,NRLIST,NVHP,NVJE,NW,
     '  PG,SE,XA,XE,XP,STRING,ERROR,*)

C#### Subroutine: DECORN
C###  Description:
C###    DECORN defines the corner nodes for a BEM problem.

C###    MPN 25-Nov-94: NVJP is used differenlty now so this routine
C###    needs to be rewritten

C###    NEW 31-Oct-94
C###    If NVJP(nc=0) > 1 Then we have a corner node.
C###    Considering the list of surrounding ne's. If the current ne is
C###    the smallest in the list then nc=2. If it is the next highest
C###    then nc=3, and if it is the highest (3D only) then nc=4.
C###    Can avoid using NVJE for this case -> NVJE to be removed.

!OLD
!C**** 7-7-92.  Sets up NVJE.
!C#### NVJE is the mapping between nn,nb,ne,nr and the nc value. Its
!C###  value is 2 unless the nn in ne is at a corner, in which case it
!C###  can take values 2,3,or 4 depending on whether ne is the 1st,2nd
!C###  or 3rd corner element.

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
C*** NRLIST needs to be passed
      INTEGER INP(NNM,NIM,NBFM),NBJ(NJM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),
     '  NRLIST(0:NRM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM)
      REAL*8 PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IPFILE,nr,nx
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,
     '  GENER,MOUSE

      CALL ENTERS('DECORN',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define corners;c/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define corner details for boundary elements. Corner details
C###    are read from or written to the file FILENAME (with extension
C###    .ipcorn) in the directory specified by PATH.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';c/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define corners;m
C###  Description:
C###    Define corner details using the mouse.
C###  Parameter:      <on WS#[1]>
C###    Specify the workstation number on which to define corners.

        OP_STRING(1)=STRING(1:IEND)//';m'
        OP_STRING(2)=BLANK(1:15)//'<on WS#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DECORN',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 24-Jan-1990
        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
C old        IF(MOUSE) CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        CALL ASSERT(NPT(1).GT.0,'>>no nodes are defined',
     '    ERROR,*9999)
        CALL ASSERT(NET(1).GT.0,'>>no elements are defined',
     '    ERROR,*9999)
        CALL ASSERT(CALL_EQUA,'>>Problem type not defined',ERROR,*9999)
        !Need to define equation first so nvne and nvjp arrays set up.

        nx=1 ! MPN 14Jun2000  may need generalising?

C For the case when some corners have already been defined
        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'corn',STATUS,
     '        ERR,ERROR,*9999)
            CALL IPCORN(INP,NBJ,NKJE,NPF,NPNE,
     '        nr,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),
     '        PG,SE,XA,XE,XP,ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('DECORN')
      RETURN
 9999 CALL ERRORS('DECORN',ERROR)
      CALL EXITS('DECORN')
      RETURN 1
      END


