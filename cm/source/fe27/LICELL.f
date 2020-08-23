      SUBROUTINE LICELL(CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '  CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,NEELEM,NPNODE,NRLIST,NXLIST,
     '  CE,CELL_RCQS_VALUE,CELL_YQS_VALUE,CELL_ICQS_NAMES,
     '  CELL_RCQS_NAMES,CELL_YQS_NAMES,CP,CQ,STRING,ERROR,*)

C#### Subroutine: LICELL
C###  Description:
C###    LICELL lists cellular material parameters
C**** Created   LKC 25-MAR-1998
C *** Re-Written - DPN June 1999

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER CELL_ICQS_VALUE(NQIM,NQVM),CELL_ICQS_SPATIAL(NQIM,NQVM),
     '  CELL_RCQS_SPATIAL(NQRM,NQVM),CELL_YQS_SPATIAL(NIQSM,NQVM),
     '  NEELEM(0:NE_R_M,0:NRM),NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NXLIST(0:NXM)
      REAL*8 CE(NMM,NEM,NXM),CELL_RCQS_VALUE(NQRM,NQVM),
     '  CELL_YQS_VALUE(NIQSM,NQVM),CP(NMM,NPM,NXM),CQ(NMM,NQM,NXM)
      CHARACTER CELL_ICQS_NAMES(NQIM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH),
     '  ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,nr,nx,nxc
      LOGICAL ALL_REGIONS,OPFILE
      CHARACTER FILE*100

      CALL ENTERS('LICELL',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C-------------------------------------------------------------------

C#### Command: FEM list cell;<FILENAME>
C###  Parameter:      <region #[1]>
C###    Specify the region number to use.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Lists cell parameters. Cell parameters can be listed to the
C###    file FILENAME.opcell.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C-------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LICELL',ERROR,*9999)
      ELSE
        CALL ASSERT(CALL_CELL,'>>Must define cell first',ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
        nr=NRLIST(1)

        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opcell','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF
        IF((ITYP3(nr,nx).LE.4).AND.(ITYP19(nr,nx).EQ.1)) THEN
          CALL OPCELL_PROMPT(NEELEM,NPNODE,nr,nx,CE(1,1,nx),CP(1,1,nx),
     '      CQ(1,1,nx),ERROR,*9999)
        ELSE
          CALL OPCELL(CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '      CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,CELL_RCQS_VALUE,
     '      CELL_YQS_VALUE,CELL_ICQS_NAMES,CELL_RCQS_NAMES,
     '      CELL_YQS_NAMES,ERROR,*9999)
        ENDIF
        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LICELL')
      RETURN
 9999 CALL ERRORS('LICELL',ERROR)
      CALL EXITS('LICELL')
      RETURN 1
      END


