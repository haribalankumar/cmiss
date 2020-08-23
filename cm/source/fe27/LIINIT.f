      SUBROUTINE LIINIT(ITHRES,NBH,NBJ,NEELEM,NHE,NHP,NKH,NODENVCB,
     '  NPNODE,NRLIST,NVCB,NVHP,NW,NWQ,NXLIST,NYNE,NYNP,
     '  AQ,THRES,YG,YP,YQ,STRING,FIX,ERROR,*)

C#### Subroutine: LIINIT
C###  Description:
C###    LIINIT lists initial and boundary conditions.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ITHRES(3,NGM,NEM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NODENVCB(NVCBM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NVCB(-1:3,NVCBM),
     '  NVHP(NHM,NPM,NCM,0:NRM),
     '  NW(NEM,3,NXM),NWQ(8,0:NQM,NAM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 AQ(NMAQM,NQM),THRES(3,NGM,NEM),YG(NIYGM,NGM,NEM),
     '  YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,nolist,nr,nx,nxc
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,CBBREV,FULL,OPFILE

      CALL ENTERS('LIINIT',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list initial<;FILENAME>
C###  Description:
C###    List initial and boundary conditions to the screen or to file
C###    FILENAME.opinit if qualifier is present.
C###  Parameter:      <full>
C###    Provide extended output.
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
        CALL DOCUM('fe27','doc','LIINIT',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opinit','NEW',
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

        DO nolist=1,NRLIST(0)
          nr=NRLIST(nolist)
          CALL ASSERT(ITYP1(nr,nx).GT.0,
     '      '>>No problem type defined',ERROR,*9999)
        ENDDO

        CALL ASSERT(CALL_INIT,
     '    '>>No initial conditions defined',ERROR,*9999)

        DO nolist=1,NRLIST(0)
          nr=NRLIST(nolist)
          CALL OPINIT(ITHRES,NBH,NBJ,NEELEM(0,nr),
     '      NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NODENVCB,
     '      NPNODE(0,nr),nr,NVCB,NVHP(1,1,1,nr),NW(1,1,nx),NWQ,nx,nxc,
     '      NYNE,NYNP,AQ,THRES,
     '      YG,YP(1,1,nx),YQ(1,1,1,nx),FIX(1,1,nx),FULL,ERROR,*9999)
C SMAR009 19/01/99 removed NVCNODE,
        ENDDO

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIINIT')
      RETURN
 9999 CALL ERRORS('LIINIT',ERROR)
      CALL EXITS('LIINIT')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
      ENDIF
      RETURN 1
      END


