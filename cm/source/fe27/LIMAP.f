      SUBROUTINE LIMAP(IBT,IDO,INP,NBH,NBJ,
     '  NEELEM,NENP,NHE,NHP,NKB,NKHE,NKH,NKJ,NKJE,NLL,NNB,
     '  NNF,NNL,NONY,NPF,NPL,NPNE,NPNODE,NPNY,
     '  NRLIST,NVHE,NVHP,NVJE,NVJP,NWP,NXI,NXLIST,
     '  NYNE,NYNO,NYNP,NYNR,NYNY,CONY,CYNO,CYNY,SE,SP,XA,XE,XP,
     '  STRING,FIX,ERROR,*)

C#### Subroutine: LIMAP
C###  Description:
C###    LIMAP lists mapping arrays.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKB(2,2,2,NNM,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NNB(4,4,4,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  NWP(NPM,2),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),NYNY(0:NYYM,NYM,NRM,NXM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),CYNY(0:NYYM,NYM,NRM,NXM),
     '  SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,nc,NCLIST(0:4),no_nc,no_nrc,
     '  nonr,nr,nrc,NRCLIST(0:3),nx,nxc
      CHARACTER FILE*100,NXTYPE*8,TYPE*8,VERSION_TYPE*11
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,OPFILE

      CALL ENTERS('LIMAP',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list map<;FILENAME> [mesh]
C###  Description:
C###    List the mappings between parameters nk,nv,nh,np,nc,nr and mesh
C###    degrees of freedom ny.  That is arrays NYNR, NYNP, NYNE (if
C###    applicable), and NPNY.
C###  Parameter:    <nrc (#s/all)[all]>
C###    Specify the nrc numbers of mappings to list.  0 corresponds to
C###    global variable numbers, 1 corresponds to local equations /
C###    rows, and 2 corresponds to local variables / columns.
C###  Parameter:    <nc (#s/all)[all]>
C###    Specify the nc (dependent variable derivative type) values of
C###    the dependent variable matrices. If the option 'all' is
C###    specified then nc values for that dependent variable matrix are
C###    used.
C###  Parameter:    <coupled/region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped. 'coupled' specifies the coupled region.
C###  Parameter:    <using (solve/fit/optimise)[solve]>
C###    Specify the problem type to list.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> [mesh]'
        OP_STRING(2)=BLANK(1:15)//'<nrc (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<nc (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<coupled/region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<using (solve/fit/optimise)'
     '    //'[solve]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list map<;FILENAME> solution
C###  Description:
C###    List the mappings between mesh degrees of freedom ny and
C###    solution degress of freedom.  That is arrays NONY, CONY, CYNO,
C###    and NYNO.
C###  Parameter:    <nrc (#s/all)[all]>
C###    Specify the nrc numbers to list.  1 corresponds to local
C###    equations / rows, and 2 corresponds to local variables / columns.
C###  Parameter:    <nc (#s/all)[all]>
C###    Specify the nc (dependent variable derivative type) values of
C###    the dependent variable matrices. If the option 'all' is
C###    specified then nc values for that dependent variable matrix are
C###    used.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped. 'coupled' specifies the coupled region.
C###  Parameter:    <using (solve/fit/optimise)[solve]>
C###    Specify the problem type to list.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> solution'
        OP_STRING(2)=BLANK(1:15)//'<nrc (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<nc (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<coupled/region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<using (solve/fit/optimise)'
     '    //'[solve]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list map<;FILENAME> versions [independent]
C###  Description:
C###    List the number of versions at each node and the versions used
C###    by each element for independent variables.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped. 'coupled' specifies the coupled region.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> versions '
     '    //'[independent]'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list map<;FILENAME> versions dependent
C###  Description:
C###    List the number of versions at each node and the versions used
C###    by each element for dependent variables.
C###  Parameter:    <nc (#s/all)[all]>
C###    Specify the nc (dependent variable derivative type) values of
C###    the dependent variable matrices. If the option 'all' is
C###    specified then nc values for that dependent variable matrix are
C###    used.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped. 'coupled' specifies the coupled region.
C###  Parameter:    <using (solve/fit/optimise)[solve]>
C###    Specify the problem type to list.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> versions '
     '    //'dependent'
        OP_STRING(2)=BLANK(1:15)//'<nc (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<coupled/region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<using (solve/fit/optimise)'
     '    //'[solve]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list map<;FILENAME> material
C###  Description:
C###    (not implemented)
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped. 'coupled' specifies the coupled region.
C###  Parameter:    <using (solve/fit/optimise)[solve]>
C###    Specify the problem type to list.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> material'
        OP_STRING(2)=BLANK(1:15)//'<coupled/region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<using (solve/fit/optimise)'
     '    //'[solve]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list map<;FILENAME> lines
C###  Description:
C###    (not implemented)
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped. 'coupled' specifies the coupled region.
C###  Parameter:    <using (solve/fit/optimise)[solve]>
C###    Specify the problem type to list.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> lines'
        OP_STRING(2)=BLANK(1:15)//'<coupled/region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<using (solve/fit/optimise)'
     '    //'[solve]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIMAP',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFILE1,'DISK',FILE(IBEG:IEND)//'.opmap','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

C MPN 13Aug2003: moved below
C old        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'MESH',2,noco+1,NTCO,N3CO)) THEN
          TYPE='MESH'
        ELSE IF(CBBREV(CO,'SOLUTION',2,noco+1,NTCO,N3CO)) THEN
          TYPE='SOLUTION'
        ELSE IF(CBBREV(CO,'LINES',2,noco+1,NTCO,N3CO)) THEN
          TYPE='LINE'
        ELSE IF(CBBREV(CO,'MATERIAL',2,noco+1,NTCO,N3CO)) THEN
          TYPE='MATERIAL'
        ELSE IF(CBBREV(CO,'VERSIONS',2,noco+1,NTCO,N3CO)) THEN
          TYPE='VERSIONS'
        ELSE
          TYPE='MESH'
        ENDIF

        IF(TYPE(1:4).EQ.'MESH') THEN
          IF(CBBREV(CO,'NRC',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO),'ALL',2)) THEN
              NRCLIST(0)=3
              DO nrc=0,2
                NRCLIST(nrc+1)=nrc
              ENDDO
            ELSE
              CALL PARSIL(CO(N3CO+1),3,NRCLIST(0),NRCLIST(1),ERROR,
     '          *9999)
              DO no_nrc=1,NRCLIST(0)
                nrc=NRCLIST(no_nrc)
                CALL ASSERT(nrc.GT.0.AND.nrc.LT.3,
     '            '>>Invalid nrc number',
     '            ERROR,*9999)
              ENDDO
            ENDIF
          ELSE
            NRCLIST(0)=3
            DO nrc=0,2
              NRCLIST(nrc+1)=nrc
            ENDDO
          ENDIF
          IF(CBBREV(CO,'NC',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'ALL',2)) THEN
              NCLIST(0)=NCM
              DO nc=1,NCM
                NCLIST(nc)=nc
              ENDDO
            ELSE
              CALL PARSIL(CO(N3CO+1),4,NCLIST(0),NCLIST(1),ERROR,*9999)
              DO no_nc=1,NCLIST(0)
                nc=NCLIST(no_nc)
                CALL ASSERT(nc.GT.0.AND.nc.LE.NCM,
     '            '>>Invalid nc number',ERROR,*9999)
              ENDDO !no_nc
            ENDIF
          ELSE
            NCLIST(0)=NCM
            DO nc=1,NCM
              NCLIST(nc)=nc
            ENDDO
          ENDIF
          IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'SOLVE',2)) THEN
              CALL ASSERT(CALL_EQUA,'>>No equation defined',ERROR,*9999)
              NXTYPE='SOLVE'
            ELSE IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
              CALL ASSERT(CALL_FIT,'>>No fit defined',ERROR,*9999)
              NXTYPE='FIT'
              NCLIST(0)=1
              NCLIST(1)=1
            ELSE IF(ABBREV(CO(N3CO+1),'OPTIMISE',2)) THEN
              IF(KTYP8.EQ.6) THEN !data fitting by optimisation
                CALL ASSERT(CALL_FIT,'>>No fit defined',ERROR,*9999)
                NCLIST(0)=1
                NCLIST(1)=1
                NRCLIST(0)=1
                NRCLIST(1)=0
              ELSE
                CALL ASSERT(CALL_OPTI,'>>No optimisation defined',
     '            ERROR,*9999)
              ENDIF
              NXTYPE='OPTIMISE'
            ELSE
              CALL ASSERT(CALL_EQUA,'>>No equation defined',ERROR,*9999)
              NXTYPE='SOLVE'
            ENDIF
          ELSE
            CALL ASSERT(CALL_EQUA,'>>No equation defined',ERROR,*9999)
            NXTYPE='SOLVE'
          ENDIF
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
        ELSE IF(TYPE(1:8).EQ.'SOLUTION') THEN
          IF(CBBREV(CO,'NRC',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO),'ALL',2)) THEN
              NRCLIST(0)=2
              DO nrc=1,2
                NRCLIST(nrc)=nrc
              ENDDO
            ELSE
              CALL PARSIL(CO(N3CO+1),3,NRCLIST(0),NRCLIST(1),ERROR,
     '          *9999)
              DO no_nrc=1,NRCLIST(0)
                nrc=NRCLIST(no_nrc)
                CALL ASSERT(nrc.GT.0.AND.nrc.LT.3,
     '            '>>Invalid nrc number',
     '            ERROR,*9999)
              ENDDO
            ENDIF
          ELSE
            NRCLIST(0)=2
            DO nrc=1,2
              NRCLIST(nrc)=nrc
            ENDDO
          ENDIF
          IF(CBBREV(CO,'NC',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'ALL',2)) THEN
              NCLIST(0)=NCM
              DO nc=1,NCM
                NCLIST(nc)=nc
              ENDDO
            ELSE
              CALL PARSIL(CO(N3CO+1),4,NCLIST(0),NCLIST(1),ERROR,*9999)
              DO no_nc=1,NCLIST(0)
                nc=NCLIST(no_nc)
                CALL ASSERT(nc.GT.0.AND.nc.LE.NCM,
     '            '>>Invalid nc number',ERROR,*9999)
              ENDDO !no_nc
            ENDIF
          ELSE
            NCLIST(0)=NCM
            DO nc=1,NCM
              NCLIST(nc)=nc
            ENDDO
          ENDIF
          IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'SOLVE',2)) THEN
              CALL ASSERT(CALL_SOLV,'>>No solution defined',ERROR,*9999)
              NXTYPE='SOLVE'
            ELSE IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
              CALL ASSERT(CALL_FIT,'>>No fit defined',ERROR,*9999)
              NXTYPE='FIT'
              NCLIST(0)=1
              NCLIST(1)=1
            ELSE IF(ABBREV(CO(N3CO+1),'OPTIMISE',2)) THEN
              CALL ASSERT(CALL_OPTI,'>>No optimisation defined',
     '          ERROR,*9999)
              IF(KTYP8.EQ.6) THEN !data fitting by optimisation
                NCLIST(0)=1
                NCLIST(1)=1
                NRCLIST(0)=1
                NRCLIST(1)=0 !MPN and AJP global variables
              ENDIF
              NXTYPE='OPTIMISE'
            ELSE
              CALL ASSERT(CALL_SOLV,'>>No solution defined',ERROR,*9999)
              NXTYPE='SOLVE'
            ENDIF
          ELSE
            CALL ASSERT(CALL_SOLV,'>>No solution defined',ERROR,*9999)
            NXTYPE='SOLVE'
          ENDIF
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
        ELSE IF(TYPE(1:8).EQ.'MATERIAL') THEN
          CALL ASSERT(CALL_OPTI,'>>No optimisation defined',ERROR,*9999)
          NXTYPE='OPTIMISE'
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
        ELSE IF(TYPE(1:4).EQ.'LINES') THEN
          CALL ASSERT(CALL_OPTI,'>>No optimisation defined',ERROR,*9999)
          NXTYPE='OPTIMISE'
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
        ELSE IF(TYPE(1:8).EQ.'VERSIONS') THEN
          CALL ASSERT(CALL_NODE,'>>Nodes not defined',ERROR,*9999)
          IF(CBBREV(CO,'DEPENDENT',2,noco+1,NTCO,N3CO)) THEN
            VERSION_TYPE='DEPENDENT'
          ELSE
            VERSION_TYPE='INDEPENDENT'
          ENDIF
          IF(VERSION_TYPE(1:9).EQ.'DEPENDENT') THEN
            IF(CBBREV(CO,'NC',2,noco+1,NTCO,N3CO)) THEN
              IF(ABBREV(CO(N3CO+1),'ALL',2)) THEN
                NCLIST(0)=NCM
                DO nc=1,NCM
                  NCLIST(nc)=nc
                ENDDO
              ELSE
                CALL PARSIL(CO(N3CO+1),4,NCLIST(0),NCLIST(1),ERROR,
     '            *9999)
                DO no_nc=1,NCLIST(0)
                  nc=NCLIST(no_nc)
                  CALL ASSERT(nc.GT.0.AND.nc.LE.NCM,
     '              '>>Invalid nc number',ERROR,*9999)
                ENDDO !no_nc
              ENDIF
            ELSE
              NCLIST(0)=NCM
              DO nc=1,NCM
                NCLIST(nc)=nc
              ENDDO
            ENDIF
            CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
            nxc=NXLIST(1)
          ELSE
            NCLIST(0)=1
            NCLIST(1)=0
          ENDIF
        ENDIF

        IF(.NOT.(TYPE(1:8).EQ.'VERSIONS'.AND.VERSION_TYPE(1:11).NE.
     '    'INDEPENDENT')) THEN
          IF(NXTYPE(1:5).EQ.'SOLVE') THEN
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this solve class',ERROR,*9999)
          ELSE IF(NXTYPE(1:3).EQ.'FIT') THEN
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this fit class',ERROR,*9999)
          ELSE IF(NXTYPE(1:8).EQ.'OPTIMISE') THEN
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this optimisation class',
     '        ERROR,*9999)
          ENDIF
        ELSE
          nx=1
        ENDIF

C new MPN 13Aug2003: adding option to interrogate mappings for coupled problems
        IF(CBBREV(CO,'COUPLED',2,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(IS_COUPLED(nx),
     &      '>>Define class (eg solve) for coupled problem',ERROR,*9999)
          NRLIST(0)=1
          NRLIST(1)=0 !the coupled region
        ELSE
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     &      ERROR,*9999)
        ENDIF

        DO nonr=1,NRLIST(0)
          nr=NRLIST(nonr)
          CALL OPMAP(IBT,IDO,INP,NBH,NBJ,
     '      NCLIST,NEELEM,NENP,NHE(1,nx),
     '      NHP(1,nr,nx),NKB,NKHE,NKH(1,1,1,nr),
     '      NKJ,NKJE,NLL,NNB,NNF,NNL,
     '      NONY(0,1,1,0,nx),NPF,NPL,NPNE,NPNODE,NPNY(0,1,0,nx),nr,
     '      NRCLIST,NVHE,NVHP(1,1,1,nr),NVJE,NVJP,NWP,nx,NXI,
     '      NYNE,NYNO(0,1,1,0,nx),NYNP,NYNR(0,0,1,0,nx),NYNY(0,1,1,nx),
     '      CONY(0,1,1,0,nx),CYNO(0,1,1,0,nx),CYNY(0,1,1,nx),SE,SP,XA,
     '      XE,XP,NXTYPE,TYPE,FIX(1,1,nx),ERROR,*9999)
        ENDDO !nonr

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIMAP')
      RETURN
 9999 CALL ERRORS('LIMAP',ERROR)
      CALL EXITS('LIMAP')
      RETURN 1
      END


