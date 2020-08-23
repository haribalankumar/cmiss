      SUBROUTINE LIMATE(ILPIN,LIST,NBJ,NEELEM,NMBIN,NPNODE,
     '  NRLIST,NW,NXLIST,CE,CGE,CIN,CP,CQ,YG,STRING,ERROR,*)

C#### Subroutine: LIMATE
C###  Description:
C###    LIMATE lists material parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ILPIN(NMM,NRM,NXM),LIST(0:NLISTM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NMBIN(NMM,NRM,NXM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),NW(NEM,3,NXM),NXLIST(0:NXM)
      REAL*8 CE(NMM,NEM,NXM),CGE(NMM,NGM,NEM,NXM),
     '  CIN(NMM,0:NGM,NNEPM),CP(NMM,NPM,NXM),
     '  CQ(NMM,NQM,NXM),YG(NIYGM,NGM,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,nolist,nr,nx,nxc
      CHARACTER FILE*100
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,CONSTIT,OPFILE

      CALL ENTERS('LIMATE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list material<;FILENAME>
C###  Description:
C###    List the material parameters of the solution domain to the
C###    screen or to file FILENAME.opmate if qualifier is present.
C###  Parameter:    <(entered/constitutive)[entered]>
C###    For finite elasticity problem, specify either the entered or
C###    (calculated) constitutive parameters.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to list. The "all" keyword indicates
C###    all currently defined regions.
C###  Parameter:      <using (fit/optimisation/solve)[solve]>
C###    Specify the problem type to list.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of problem type) to list.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<(entered/constitutive)[entered]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)
     '    //'<using (fit/optimisation/solve)[solve]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)


C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIMATE',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opmate','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL ASSERT(CALL_MATE,'>>Material properties not defined',
     '    ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'SOLVE',2)) THEN
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this solve class',ERROR,*9999)
          ELSE IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this fit class',ERROR,*9999)
          ELSE IF(ABBREV(CO(N3CO+1),'OPTIMISE',2)) THEN
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this optimisation class',
     '        ERROR,*9999)
          ELSE
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this solve class',ERROR,*9999)
          ENDIF
        ELSE
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,
     '      '>>No nx defined for this solve class',ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'CONSTITUTIVE',1,noco+1,NTCO,N3CO)) THEN
          CONSTIT=.TRUE.
        ELSE
          CONSTIT=.FALSE.
        ENDIF

C LKC 26-MAR-1998 Split off into opcell
C        IF(CBBREV(CO,'CELL',2,noco+1,NTCO,N3CO)) THEN
C          CELL=.TRUE.
C        ELSE
C          CELL=.FALSE.
C        ENDIF

        DO nolist=1,NRLIST(0)
          nr=NRLIST(nolist)
          CALL OPMATE(ILPIN(1,nr,nx),LIST,NBJ,
     '      NEELEM,NMBIN(1,nr,nx),NPNODE,nr,NW(1,1,nx),nx,nxc,
     '      CE(1,1,nx),CGE(1,1,1,nx),
     '      CIN,CP(1,1,nx),CQ(1,1,nx),YG,CONSTIT,
     '      ERROR,*9999)
        ENDDO

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIMATE')
      RETURN
 9999 CALL ERRORS('LIMATE',ERROR)
      CALL EXITS('LIMATE')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
      ENDIF
      RETURN 1
      END


