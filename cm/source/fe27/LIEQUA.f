      SUBROUTINE LIEQUA(NBH,NEELEM,NELIST,NHE,NHP,NKH,
     '  NPL,NPNODE,NRLIST,NVHE,NVHP,NW,NXLIST,STRING,ERROR,*)

C#### Subroutine: LIEQUA
C###  Description:
C**** LIEQUA lists equation parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPL(5,0:3,NLM),NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),NW(NEM,3,NXM),
     '  NXLIST(0:NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,no_nrlist,nr,nx,nxc
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,CBBREV,FULL,OPFILE

      CALL ENTERS('LIEQUA',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list equation<;FILENAME>
C###  Description:
C###    List equation type, solution method, and dependent variable
C###    interpolation to the screen or to file FILENAME.opequa if
C###    qualifier is present.
C###  Parameter:      <full>
C###    Provide extended output including node version information.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to list. The "all" keyword indicates
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to list.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<full>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIEQUA',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opequa','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'FULL',2,noco+1,NTCO,N3CO)) THEN
          FULL=.TRUE.
        ELSE
          FULL=.FALSE.
        ENDIF

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          IF(ITYP2(nr,nx).GT.0) THEN
            CALL OPEQUA(NBH,NEELEM,NELIST,NHE(1,nx),NHP(1,nr,nx),
     '        NKH(1,1,1,nr),
     '        NPL,NPNODE,nr,NVHE,NVHP(1,1,1,nr),NW(1,1,nx),nx,nxc,FULL,
     '        ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'('' >>Problem type not defined for '
     '        //'region '',I2)') nr
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !no_nrlist

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIEQUA')
      RETURN
 9999 CALL ERRORS('LIEQUA',ERROR)
      CALL EXITS('LIEQUA')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
      ENDIF
      RETURN 1
      END


