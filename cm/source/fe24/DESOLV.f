      SUBROUTINE DESOLV(IBT,IDO,INP,NAN,NAQ,NBH,NBJ,NEELEM,NELIST,
     '  NENP,NHE,NKB,NKHE,NKJE,NLL,NNB,NNF,NNL,NONY,
     '  NP_INTERFACE,NPF,NPL,NPNE,NPNODE,NPNY,NRE,NRLIST,NVHE,NVHP,NVJE,
     '  NW,NWP,NWQ,NXI,NXLIST,NXQ,NYNE,NYNO,NYNP,NYNR,NYNY,NYQNR,
     '  AQ,CONY,CYNO,CYNY,SE,SP,XA,XE,XP,YP,STRING,FIX,ERROR,*)

C#### Subroutine: DESOLV
C###  Description:
C###    DESOLV defines solution parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'gen000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp90.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NAQ(NQM,NAM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NHE(NEM,NXM),NKB(2,2,2,NNM,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NP_INTERFACE(0:NPM,0:3),
     '  NPF(9,NFM),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NWP(NPM,2),
     '  NWQ(8,0:NQM,NAM),NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),NYNY(0:NYYM,NYM,NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 AQ(NMAQM,NQM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),CYNY(0:NYYM,NYM,NRM,NXM),
     '  SE(NSM,NBFM,NEM),
     '  SP(NKM,NBFM,NPM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE,N3CO,nc,
     '  nonr,no_coup_nrlist,no_nrlist,no_nynr,nr,nrc,nx,nxc,ny
      CHARACTER FILE*(MXCH),FILE_NAME*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,CBBREV,FILIO,FIRST_TIME,
     '  GENER,LINE_NAME,MOUSE,SEND_GLOBAL

      CALL ENTERS('DESOLV',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define solve;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define the details of the solution procedure and the level of
C###    solver output desired. These are read from or written to the file
C###    FILENAME.ipsolv in the directory PATH.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <(separate/coupled)[separate]>
C###    Specify whether to solve the problem as a collection of separate
C###    regions or to couple the regions together. The 'coupled' option
C###    sets up the coupled solution mapping arrays whereas the default
C###    'separate' option does not.
C###  Parameter:      <nosend_global>
C###    For parallel processing problems, signal the solver to not
C###    re-send the global arrays and common blocks to the remote
C###    processes, as the user may decide that they will not change
C###    (saves time). Note the default is to (re)send global arrays and
C###    common blocks.
C###  Parameter:      <name FILENAME>
C###    Specify the filename for output of time-dependent problems.
C###    The default is the filename specified in the .ipsolv file.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<(separate/coupled)[separate]>'
        OP_STRING(5)=BLANK(1:15)//'<nosend_global>'
        OP_STRING(6)=BLANK(1:15)//'<name FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DESOLV',ERROR,*9999)
      ELSE
        IPFILE=3 !is input file version number on 2-Mar-1994

        CALL PARSE_QUALIFIERS('DLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        CALL ASSERT(CALL_NODE,'>>Nodes not defined',ERROR,*9999)
        CALL ASSERT(CALL_ELEM,'>>Elements not defined',ERROR,*9999)
        CALL ASSERT(CALL_EQUA,'>>Problem type not defined',ERROR,*9999)
!new AJP 12-1-94  Need to define initial conditions first in every case?
! DPN 6-5-99 NO, for cellular modelling only need to define CELL not
!            initials.
        IF (ITYP2(NRLIST(1),NXLIST(1)).NE.12) THEN
          CALL ASSERT(CALL_INIT,'>>Define b.c.s first',ERROR,*9999)
        ELSE
          CALL ASSERT(CALL_CELL,'>>Define cell first',ERROR,*9999)
          CALL_MATE=.TRUE. !don't need to define a material for a
                           !single cell
          CALL_INIT=CALL_CELL !need to define CELL but not b.c.'s
        ENDIF
c        IF(ITYP5(NRLIST(1)).EQ.3) THEN !modal analysis
c          CALL ASSERT(CALL_INIT,'>>Define b.c.s first',ERROR,*9999)
c        ENDIF
! end new

C CPB 21/11/94 Adding NX_LOC
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'NOSEND_GLOBAL',2,noco+1,NTCO,N3CO)) THEN
          SEND_GLOBAL=.FALSE.
        ELSE
          SEND_GLOBAL=.TRUE.
        ENDIF
        IF(CBBREV(CO,'NAME',1,noco+1,NTCO,N3CO)) THEN
          LINE_NAME=.TRUE.
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          FILE_NAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          LINE_NAME=.FALSE.
        ENDIF

C****   New for coupled AJP/CPB 23/2/95
C JWF 24/06/02   Modifed for Contact Mechanics.
C Problem is also considered coupled if contact option was
C selected in (ipequ). ie KTYP5G(nr).GE.1

        IF(CBBREV(CO,'COUPLED',2,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT((KTYP90.GT.0).OR.(KTYP5G(1).GE.1),'
     '      >> Define coupling first',ERROR,*9999)
          CALL ASSERT(NRLIST(0).GT.1,'>>Must couple > 1 regions',
     '      ERROR,*9999)
          IS_COUPLED(nx)=.TRUE.

C#### IS_COUPLED determines whether whether a coupled problem is
C###    solved as coupled or as a collection of separate regions.
          DO nonr=1,NRLIST(0) !hold on to list of coupled regions
            COUP_NRLIST(nonr,nx)=NRLIST(nonr)
          ENDDO
          COUP_NRLIST(0,nx)=NRLIST(0)
C****     Need to sort the coupled region list to ensure the master
C****     region (lower region #) mappings are calculated before
C****     slave regions

          CALL ISORT(COUP_NRLIST(0,nx),COUP_NRLIST(1,nx))

C****     Calculate list of ny's for coupled region nr=0
          DO nc=1,NCM
            DO nrc=0,2
              NYNR(0,nrc,nc,0,nx)=0
            ENDDO !nrc
          ENDDO !nc
          DO no_coup_nrlist=1,COUP_NRLIST(0,nx)
            nr=COUP_NRLIST(no_coup_nrlist,nx)
            DO nc=1,NCT(nr,nx)
              DO nrc=0,2
                DO no_nynr=1,NYNR(0,nrc,nc,nr,nx)
                  ny=NYNR(no_nynr,nrc,nc,nr,nx)
                  NYNR(0,nrc,nc,0,nx)=NYNR(0,nrc,nc,0,nx)+1
                  IF(NYNR(0,nrc,nc,0,nx).LE.NY_R_M) THEN
                    NYNR(NYNR(0,nrc,nc,0,nx),nrc,nc,0,nx)=ny
                  ENDIF
                ENDDO !no_nynr
                CALL ASSERT(NYNR(0,nrc,nc,0,nx).LE.NY_R_M,
     '            '>>Increase NY_R_M',ERROR,*9999)
              ENDDO !nrc
            ENDDO !nc
          ENDDO !no_coup_nrlist

C***      Check to see if there are any BEM regions in the coupled
C***      problem

C#### COUPLED_BEM is true if any of the regions in the coupled problem
C###    use the Boundary Element Method and is false otherwise.

          COUPLED_BEM(nx)=.FALSE.
          DO no_coup_nrlist=1,COUP_NRLIST(0,nx)
            nr=COUP_NRLIST(no_coup_nrlist,nx)
            IF(ITYP4(nr,nx).EQ.2.OR.ITYP4(nr,nx).EQ.3) THEN
              COUPLED_BEM(nx)=.TRUE.
            ENDIF
          ENDDO !no_coup_nrlist

        ELSE
          IS_COUPLED(nx)=.FALSE.
          COUPLED_BEM(nx)=.FALSE.
          COUP_NRLIST(0,nx)=0
        ENDIF

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'solv',STATUS,
     '        ERR,ERROR,*9999)

            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)

              CALL IPSOLV(IBT,IDO,INP,NAN,NAQ,NBH,NBJ,NEELEM,NELIST,
     '          NENP,NHE(1,nx),NKB,NKHE,NKJE,NLL,NNB,NNF,NNL,
     '          NONY(0,1,1,0,nx),NP_INTERFACE,NPF,NPL,NPNE,NPNODE,
     '          NPNY(0,1,0,nx),nr,NRE,NVHE,NVHP,NVJE,NW(1,1,nx),NWP,NWQ,
     '          nx,NXI,NXQ,NYNE,NYNO(0,1,1,0,nx),NYNP,NYNR(0,0,1,0,nx),
     '          NYNY(0,1,1,nx),NYQNR(0,0,1,0,nx),AQ,
     '          CONY(0,1,1,0,nx),CYNO(0,1,1,0,nx),CYNY(0,1,1,nx),SE,SP,
     '          XA,XE,XP,YP,FIX(1,1,nx),LINE_NAME,SEND_GLOBAL,FILE_NAME,
     '          ERROR,*9999)
            ENDDO
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF
        FIRSTS(nx)=.TRUE.
        CALL_SOLV=.TRUE.
      ENDIF

      CALL EXITS('DESOLV')
      RETURN
 9999 CALL ERRORS('DESOLV',ERROR)
      CALL EXITS('DESOLV')
      RETURN 1
      END


