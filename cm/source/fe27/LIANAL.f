      SUBROUTINE LIANAL(NPNODE,NRLIST,NXLIST,NYNP,STRING,YP,ERROR,*)

C#### Subroutine: LIANAL
C###  Description:
C###    LIANAL lists analytic formula parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NRLIST(0:NRM),NXLIST(0:NXM),NPNODE(0:NP_R_M,0:NRM),
     '   NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,no_nrlist,nr,nx,nxc,N3CO
      CHARACTER FILE*100
CC AJPs - 191297
      LOGICAL ACTIVATION,ALL_REGIONS,CBBREV,OPFILE
CC AJPe

      CALL ENTERS('LIANAL',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list analytic<;FILENAME>
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specifies the class number (of solve type) of the elements to
C###    list
C###  Description:
C###    Lists analytic parameters for the analytic solution.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
CC AJPs - 191297 - RGB

C#### Command: FEM list analytic<;FILENAME>
C###  Parameter:    activation
C###  Parameter:    <region (#s/all)[1]>
C###  Parameter:    <class #[1]>
C###  Description:
C###    Lists analytic activation parameters

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'activation'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
CC AJPe

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIANAL',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opanal','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF
CC AJPs 23/10/97 - 191297 - RGB
        IF(CBBREV(CO,'ACTIVATION',1,noco+1,NTCO,N3CO)) THEN
          ACTIVATION=.TRUE.
        ELSE
          ACTIVATION=.FALSE.
        ENDIF
CC AJPe
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
CC AJPs
          CALL OPANAL(NPNODE,nr,nx,NYNP,YP,ACTIVATION,ERROR,*9999)
CC AJPe
         ENDDO !no_nrlist (nr)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIANAL')
      RETURN
 9999 CALL ERRORS('LIANAL',ERROR)
      CALL EXITS('LIANAL')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
      ENDIF
      RETURN 1
      END


