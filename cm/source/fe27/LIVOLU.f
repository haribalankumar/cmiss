      SUBROUTINE LIVOLU(NBJ,NEELEM,NELIST,NKJE,NPF,NPNE,NRLIST,NVJE,NW,
     '  PG,RG,SE,WG,XA,XE,XG,XN,XP,STRING,ERROR,*)

C#### Subroutine: LIVOLU
C###  Description:
C###    LIVOLU lists the volume enclosed by a selected group of
C###    elements. Currently only implemented for closed groups of
C###    2d elements in 3d space.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM)
      REAL*8 PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,nx
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,OPFILE

      CALL ENTERS('LIVOLU',*9999)

      OPFILE=.FALSE.

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list volume<;FILENAME>
C###  Parameter:         <element (#s/all)[all]>
C###    Specifies the elements numbers to include. The 'all' command
C###    prompts all currently defined elements to be included.
C###  Parameter:         <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###  Description:
C###    LIVOLU lists the volume enclosed by a selected group of
C###    elements. Currently only implemented for closed groups of
C###    2d elements in 3d space.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C----------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIVOLU',ERROR,*9999)
      ELSE
        CALL ASSERT(NJT.EQ.3,'>> Only implemented for 3d',ERROR,*9999)
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opvolu','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        nx=1 ! MPN 14Jun2000  may need generalising?

        CALL OPVOLU(NBJ,NELIST,NKJE,NPF,NPNE,NRLIST,NVJE,NW(1,1,nx),
     '    PG,RG,SE,WG,XA,XE,XG,XN,XP,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIVOLU')
      RETURN
 9999 CALL ERRORS('LIVOLU',ERROR)
      CALL EXITS('LIVOLU')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
      ENDIF
      RETURN 1
      END


