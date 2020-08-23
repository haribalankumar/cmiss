      SUBROUTINE LIGAUS(NBH,NBJ,NEELEM,NELIST,NGLIST,NHE,NHP,NKH,
     '  NKHE,NKJE,NPF,NPNE,NPNODE,NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,
     '  NYNP,CURVCORRECT,PG,RG,SE,XA,XE,XG,XIG,XP,YG,YP,ZA,ZE,ZG,ZP,
     '  STRING,ERROR,*)

C#### Subroutine: LIGAUS
C###  Description:
C###    LIGAUS lists element Gauss point array XG(nj,nu).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NGLIST(0:NGM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,MIN0,N3CO,nb,no_nglist,no_nrlist,nr,
     '  NUMYGCMPTS,nx,nxc
      CHARACTER FILE*100,TYPE*8
      LOGICAL ALL_REGIONS,CBBREV,OPFILE

      CALL ENTERS('LIGAUS',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list gauss<;FILENAME>
C###  Description:
C###    Lists Gauss point values of geometric variables (default) or
C###    pressures in specified elements to screen or file
C###    FILENAME.opgaus if qualifier present.
C###  Parameter:    <(xg/yg/xg&yg/pressure/solution)[xg&yg]>
C###    Specify what type of output is required. XG is position
C###    information, YG is the gauss point solution array and
C###    pressure/solution are both calculated prpoerties.
C###  Parameter:    <element (#s/all)[all]>
C###    Specify the list of elements in which to list the data. The
C###    'all' option will list the data for all currently defined
C###    elements in the specified regions.
C###  Parameter:    <basis (#/geom)[geom]>
C###    Specify the basis to use to output the information. By
C###    default the basis which defines the geometry is used.
C###  Parameter:    <at (GAUSS_PT#s/all)[all]>
C###    This parameter allows you to specify a subset of the gauss
C###    points within each element. The default is 'all' so all
C###    gauss points are listed.
C###  Parameter:    <numygcmpts #[1]>
C###    Specify the number of YG components to write out.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//
     '    '<(xg/yg/xg&yg/pressure/solution)[xg&yg]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<basis (#/geom)[geom]>'
        OP_STRING(5)=BLANK(1:15)//'<at (GAUSS_PT#s/all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<numygcmpts #[1]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list gauss<;FILENAME> groups
C###  Description:
C###    This command lists all Gauss point groups. Group information
C###    is written to the file FILENAME (with extension .opgaus)

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> groups'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIGAUS',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opgaus','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        IF(CBBREV(CO,'PRESSURE',1,noco+1,NTCO,N3CO)) THEN
          TYPE(1:8)='PRESSURE'
        ELSE IF(CBBREV(CO,'SOLUTION',1,noco+1,NTCO,N3CO)) THEN
          TYPE(1:8)='SOLUTION'
        ELSE IF(CBBREV(CO,'XG',1,noco+1,NTCO,N3CO)) THEN
          TYPE(1:2)='XG'
        ELSE IF(CBBREV(CO,'YG',1,noco+1,NTCO,N3CO)) THEN
          TYPE(1:2)='YG'
        ELSE IF(CBBREV(CO,'GROUPS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='GROUPS'
        ELSE
          TYPE(1:5)='XG&YG'
        ENDIF

        IF(CBBREV(CO,'NUMYGCMPTS',1,noco+1,NTCO,N3CO)) THEN
          NUMYGCMPTS=IFROMC(CO(N3CO+1))
C news MPN 1Jun2000: first index of YG is dimensioned to NIYGM not NJM
          CALL ASSERT(NUMYGCMPTS.GT.0.AND.NUMYGCMPTS.LE.NIYGM,
     '      '>>ERROR: Must satisfy {0<NUMYGCMPTS<=NIYGM}',ERROR,*9999)
        ELSE
C news MPN 1Jun2000: changed default to min(NIYGM,3)
          NUMYGCMPTS=MIN0(NIYGM,3)
C old          NUMYGCMPTS=3
        ENDIF

        IF(CBBREV(CO,'BASIS',2,noco+1,NTCO,N3CO)) THEN
          nb=IFROMC(CO(N3CO+1))
        ELSE
          nb=0
        ENDIF

        IF(TYPE(1:6).NE.'GROUPS') THEN
          IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NGM,NGLIST(0),NGLIST(1),ERROR,*9999)
          ELSE
            IF(nb.EQ.0) THEN
              NGLIST(0)=NGT(NBJ(1,NEELEM(1,NRLIST(1))))
            ELSE
              NGLIST(0)=NGT(nb)
            ENDIF
            DO no_nglist=1,NGLIST(0)
              NGLIST(no_nglist)=no_nglist
            ENDDO
          ENDIF
        ENDIF

        IF(TYPE(1:6).EQ.'GROUPS') THEN
          CALL OPGAUSG(ERROR,*9999)
        ELSE
          DO no_nrlist=1,NRLIST(0) !region list
            nr=NRLIST(no_nrlist) !region#

            IF(TYPE(1:8).EQ.'SOLUTION') THEN
              CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr,nx),NKH(1,1,1,nr),
     '          NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,
     '          ERROR,*9999)
            ENDIF !type

            CALL OPGAUSX(nb,NBH,NBJ,NELIST,NGLIST,NHE,NKHE,NKJE,
     '        NPF,NPNE,nr,NUMYGCMPTS,NVHE,NVJE,NW,nx,
     '        CURVCORRECT,PG,RG,SE,XA,XE,XG,XIG,XP,YG,ZA,ZE,ZG,ZP,TYPE,
     '        ERROR,*9999)
          ENDDO
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIGAUS')
      RETURN
 9999 CALL ERRORS('LIGAUS',ERROR)
      CALL EXITS('LIGAUS')
      RETURN 1
      END


