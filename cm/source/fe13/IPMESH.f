      SUBROUTINE IPMESH(IBT,IDO,INP,LD,MIN_ORDER,NBJ,NBJF,NEELEM,NEL,
     '  NELIST,NELIST2,NENFVC,NENP,NEP,NFF,NFFACE,NFVC,NKJE,NKEF,NKJ,
     '  NLF,NLL,NLLINE,NNB,NNF,NNL,NODENVC,NODENVCB,NORD,NP_ATTACH,
     &  NP_INTERFACE,NPF,NPL,NPLIST,NPNE,NPNF,NPNODE,NPQ,NQLIST,
     &  nr_coronary,NRE,Nrefine,nr_host,NRLIST,NVCNODE,NVJE,NVJF,NVJL,
     &  NVJP,NXI,NXQ,TERMINAL_ORDER,BBM,CE,DF,DL,PG,RG,MRNA,RNA,
     &  SE,SF,Spread,VC,VC_INIT,WG,XA,XAB,XE,XG,XID,XIP,XNFV,XP,XQ,ZA,
     &  ZD,GROUP_NAME,INPUT_PATH,ADD_SUPER,ARTERIES,CALCU,LUNG_TOTAL,
     &  MAKE_GROUP,ORDER,REDUCE,VEINS,ERROR,*)
      
C#### Subroutine: IPMESH
C###  Description:
C###    IPMESH defines mesh parameters for specialized problems.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'ptr01.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'parameters.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),MIN_ORDER,NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NELIST2(0:NEM),NENFVC(0:NFVCM,NFVM),NENP(NPM,0:NEPM,0:NRM),
     '  NEP(NPM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NFVC(2,0:NFVCM,NVCM),
     '  NKJE(NKM,NNM,NJM,NEM),NKEF(0:4,16,6,NBFM),NKJ(NJM,NPM),
     '  NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),NNB(4,4,4,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NODENVC(NVCM),NODENVCB(NVCBM),
     '  NORD(5,NE_R_M),NPF(9,NFM),NPL(5,0:3,NLM),NPLIST(0:NPM),
     &  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     &  NPNODE(0:NP_R_M,0:NRM),NPQ(NQM),NP_ATTACH,
     &  NP_INTERFACE(0:NPM,0:3),NQLIST(0:NQM),nr_coronary,NRE(NEM),
     &  Nrefine,nr_host,NRLIST(0:NRM),NVCNODE(2,NP_R_M),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),
     &  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  NXQ(-NIM:NIM,0:4,0:NQM,NAM),TERMINAL_ORDER
      REAL*8 BBM(2,NEM),CE(NMM,NEM),DF(NFM),DL(3,NLM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),MRNA(3),RNA(3),
     &  SE(NSM,NBFM,NEM),SF(NSM,NBFM),Spread,VC(0:NVCM),VC_INIT(2,NVCM),
     &  WG(NGM,NBM),XA(NAM,NJM,NEM),XAB(NORM,NEM),XE(NSM,NJM),
     &  XG(NJM,NUM),XID(NIM,NDM),XIP(NIM,NPM),XNFV(-(NJM+1):NJM,NFVM),
     &  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM)
      LOGICAL ADD_SUPER,ARTERIES,CALCU,LUNG_TOTAL,MAKE_GROUP,ORDER,
     &  REDUCE,VEINS
      CHARACTER*(*) INPUT_PATH
      CHARACTER GROUP_NAME*30,ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,ICHAR,IEND,INFO,nb_face,noelem,noelem_start,
     &  no_nrlist,NOQUES,nr
      INTEGER*4 NSTORE_PTR,RSTORE_PTR,XSTORE_PTR
      LOGICAL FILEIP
      CHARACTER STRING*255

      CALL ENTERS('IPMESH',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(ARTERIES.OR.VEINS.OR.LUNG_TOTAL)THEN
        JTYP14=4               !stochastic fractal tree
      ELSE

        IF(NJT.EQ.2)THEN
          FORMAT='('' Enter mesh type [1]:'''//
     '      '/''   (1) Regular'''//
     '      '/''   (2) Circular mesh'''//
     '      '/''   (3) Regular fractal tree'''//
     '      '/''   (4) Stochastic fractal tree'''//
     '      '/''   (5) '''//
     '      '/''   (6) Coronary tree'''//
     '      '/''   (7) Grid point'''//
     '      '/''   (8) Purkinje fibre tree'''//
     '      '/''   (9) Arbitrary 2D domain mesh'''//
     '      '/''  (10) Surface coronary arteries'''//
     '      '/''  (11) Voronoi mesh'''//
     '      '/$,''    '',I2)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=JTYP14
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,11,LDATA,
     '      LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) JTYP14=IDATA(1)

        ELSE IF(NJT.EQ.3)THEN
          FORMAT='('' Enter mesh type [1]:'''//
     '      '/''   (1) Regular'''//
     '      '/''   (2) Eccentric spheres mesh'''//
     '      '/''   (3) Regular fractal tree'''//
     '      '/''   (4) Stochastic fractal tree'''//
     '      '/''   (5) Cylindrical mesh (unclosed)'''//
     '      '/''   (6) Coronary tree'''//
     '      '/''   (7) Grid point'''//
     '      '/''   (8) Purkinje fibre tree'''//
     '      '/''   (9) DNA'''//
     '      '/''  (10) Surface coronary arteries '''//
     '      '/''  (11) Voronoi mesh'''//
     '      '/$,''    '',I2)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=JTYP14
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,11,LDATA,
     '      LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) JTYP14=IDATA(1)
        ENDIF
      ENDIF

      IF(JTYP14.EQ.1) THEN      !Regular mesh
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL INIT_NJ_LOC(NJT,NJL_GEOM,nr,ERROR,*9999)
C GMH 14/2/97 Destroy this region and recreate -
C             necessary for cmgui locking
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_DESTROY(nr,ERROR,*9999)
          ENDIF
          CALL IPMESH1(NBJ,NEELEM,NENP,NKJE,NKJ,NP_INTERFACE,
     '      NPNE,NPNODE,nr,NRE,NVJE,NVJP,SE,XP,CALCU,ERROR,*9999)
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_CREATE(nr,ERROR,*9999)
          ENDIF
        ENDDO !no_nrlist

      ELSE IF(JTYP14.EQ.2) THEN !eccentric spheres mesh
C GMH 14/2/97 Region updating inside call
        CALL IPMESH3(IBT,NBJ,NEELEM,NENP,NKJE,NKJ,
     '    NP_INTERFACE,NPNE,NPNODE,NRE,NVJE,NVJP,SE,XP,ERROR,*9999)

      ELSE IF(JTYP14.EQ.3.OR.JTYP14.EQ.4) THEN !Fractal tree

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL INIT_NJ_LOC(NJT,NJL_GEOM,nr,ERROR,*9999)
C GMH 14/2/97 Destroy this region and recreate -
C             necessary for cmgui locking
c stupid if adding a mesh
          IF(.NOT.ADD)THEN
            IF(NPNODE(0,nr).GT.0) THEN
              CALL REGION_DESTROY(nr,ERROR,*9999)
            ENDIF
          ENDIF
          noelem_start=NEELEM(0,nr) !records highest elem # before new created
          NSTORE_PTR=0
          RSTORE_PTR=0
          XSTORE_PTR=0
          CALL ALLOCATE_MEMORY(NE_R_M*3,1,INTTYPE,NSTORE_PTR,.TRUE.,
     &      ERROR,*9999)
          CALL ALLOCATE_MEMORY(NE_R_M*3,1,DPTYPE,RSTORE_PTR,.TRUE.,
     &      ERROR,*9999)
          CALL ALLOCATE_MEMORY(NE_R_M*3,1,DPTYPE,XSTORE_PTR,.TRUE.,
     &      ERROR,*9999)
          
          CALL IPMESH2(MIN_ORDER,NBJ,NEELEM,NELIST,NELIST2,NENP,NKJE,
     '      NKJ,NLL,NNB,NORD,NP_ATTACH,NP_INTERFACE,NPL,NPLIST,
     &      NPNE,NPNODE,nr,NRE,Nrefine,nr_host,%VAL(NSTORE_PTR),NVJE,
     &      NVJP,NXI,TERMINAL_ORDER,BBM,CE,%VAL(RSTORE_PTR),SE,Spread,
     &      XAB,XP,%VAL(XSTORE_PTR),ADD_SUPER,ARTERIES,LUNG_TOTAL,ORDER,
     &      REDUCE,VEINS,ERROR,*9999)

          CALL FREE_MEMORY(NSTORE_PTR,ERROR,*9999)
          CALL FREE_MEMORY(RSTORE_PTR,ERROR,*9999)

          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_CREATE(nr,ERROR,*9999)
          ENDIF
        ENDDO


      ELSE IF(JTYP14.EQ.5) THEN !cylindrical mesh
C GMH 14/2/97 Region updating inside call
        CALL IPMESH5(NBJ,NEELEM,NENP,NKJE,NKJ,NP_INTERFACE,
     '    NPNE,NPNODE,NRE,NVJE,NVJP,SE,XP,ERROR,*9999)

      ELSE IF(JTYP14.EQ.6) THEN !coronary mesh
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL INIT_NJ_LOC(NJT,NJL_GEOM,nr,ERROR,*9999)
C GMH 14/2/97 Destroy this region and recreate -
C             necessary for cmgui locking
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_DESTROY(nr,ERROR,*9999)
          ENDIF
          IF(NPE_PTR.EQ.0) THEN
             CALL ALLOCATE_MEMORY((NPM+1)*NEM*NRM,0,
     '            INTTYPE,NPE_PTR,.TRUE.,ERROR,*9999)
          ENDIF
          CALL IPMESH6(IBT,IDO,INP,NBJ,NEELEM,NENP,NEP,NKJE,
     '         NKJ,%VAL(NPE_PTR),NPF,NP_INTERFACE,NPNE,NPNODE,nr,
     '         nr_coronary,NRE,NXI,NVJE,NVJP,SE,XA,XE,XIP,XP,
     '         INPUT_PATH,ERROR,*9999)
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_CREATE(nr,ERROR,*9999)
          ENDIF
        ENDDO !no_nrlist

      ELSE IF(JTYP14.EQ.7) THEN !grid point mesh
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL INIT_NJ_LOC(NJT,NJL_GEOM,nr,ERROR,*9999)
C GMH 14/2/97 Destroy this region and recreate -
C             necessary for cmgui locking
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_DESTROY(nr,ERROR,*9999)
          ENDIF
          CALL IPMESH7(NBJ,NEELEM,NENP,NKJE,NKJ,NPNE,NPNODE,
     '      NPQ,NQLIST,nr,NRE,NVJE,NVJP,NXQ,SE,XP,XQ,ERROR,*9999)
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_CREATE(nr,ERROR,*9999)
          ENDIF
        ENDDO !no_nrlist

      ELSE IF(JTYP14.EQ.8) THEN !Purkinje fibre mesh
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL INIT_NJ_LOC(NJT,NJL_GEOM,nr,ERROR,*9999)
C GMH 14/2/97 Destroy this region and recreate -
C             necessary for cmgui locking
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_DESTROY(nr,ERROR,*9999)
          ENDIF
          CALL IPMESH8(IBT,IDO,INP,NBJ,NEELEM,NKJE,NKJ,NPF,
     '      NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,SE,XA,XE,XP,
     '      ERROR,*9999)
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_CREATE(nr,ERROR,*9999)
          ENDIF
        ENDDO !no_nrlist

      ELSE IF(JTYP14.EQ.9) THEN      !DNA mesh
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL INIT_NJ_LOC(NJT,NJL_GEOM,nr,ERROR,*9999)
          CALL IPMESH9(NBJ,NEELEM,NENP,NKJE,NKJ,NP_INTERFACE,
     '      NPNE,NPNODE,nr,NRE,NVJE,NVJP,MRNA,RNA,SE,XP,ERROR,*9999)
        ENDDO !no_nrlist
      ELSE IF(JTYP14.EQ.10) THEN !Coronary surface mesh
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL INIT_NJ_LOC(NJT,NJL_GEOM,nr,ERROR,*9999)
C GMH 14/2/97 Destroy this region and recreate -
C             necessary for cmgui locking
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_DESTROY(nr,ERROR,*9999)
          ENDIF
          IF(NPE_PTR.EQ.0) THEN
             CALL ALLOCATE_MEMORY((NPM+1)*NEM*NRM,0,
     '            INTTYPE,NPE_PTR,.TRUE.,ERROR,*9999)
          ENDIF
          CALL IPMESH10(IBT,IDO,INP,LD,NBJ,NEELEM,NENP,NEP,NKJE,
     '      %VAL(NPE_PTR),NPF,NP_INTERFACE,NPNE,NPNODE,nr,nr_coronary,
     '         NRE,NVJE,SE,XA,XE,XID,XIP,XP,ZD,ERROR,*9999)
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_CREATE(nr,ERROR,*9999)
          ENDIF
        ENDDO !no_nrlist

      ELSEIF(JTYP14.EQ.11) THEN !Voronoi mesh
        CALL ASSERT(NRLIST(0).EQ.1,'>>Can only have one Voronoi region',
     '    ERROR,*9999)
        nr=NRLIST(1)
        FORMAT='('' The Delaunay elements are to be [1]:'''//
     '    '/''   (1) Retained'''//
     '    '/''   (2) Discarded'''//
     '    '/$,''    '',I2)'
        IDEFLT(1)=1
        DEL_RET=1
        IF(IOTYPE.EQ.3) IDATA(1)=DEL_RET
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,LDATA,
     '    LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) DEL_RET=IDATA(1)
        IF(DEL_RET.EQ.2)THEN
          FORMAT='('' Enter the face basis function [1]: '',I3)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBFM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) nb_face=IDATA(1)
        ENDIF
C        NF_LIST_PTR=0
C        CALL ALLOCATE_MEMORY(NFM+1,0,INTTYPE,NF_LIST_PTR,.TRUE.,ERROR,
C     '    *9999)
        CALL IPMESH11(IBT,NBJ,nb_face,NEELEM,NENFVC,NENP,
     '    NFVC,NKJ,NKJE,NODENVC,NODENVCB,NP_INTERFACE,NPLIST,NPNE,
     '    NPNODE,nr,NRE,NVCNODE,NVJE,NVJP,NXI,SE,VC,VC_INIT,XNFV,XP,ZA,
     '    ERROR,*9999)
C        CALL FREE_MEMORY(NF_LIST_PTR,ERROR,*9999)
      ELSE
        CALL FLAG_ERROR(0,'Invalid mesh type')
        GOTO 9998

      ENDIF !jtyp14

C CPB Adding average arc-length scaling so must determine NXI before SE
      IF((ITYP5(1,1).EQ.2.AND.((ITYP2(1,1).EQ.5.AND.ITYP3(1,1).EQ.2)
     '  .OR.(ITYP2(1,1).EQ.11))).OR.(JTYP14.EQ.3.OR.JTYP14.EQ.4))THEN
C rgb dont calcuate if using voronoi cells
      ELSEIF(JTYP14.NE.11.OR.(JTYP14.EQ.11.AND.DEL_RET.EQ.2)) THEN
C     only call NENXI, LINCAL and FACCAL if not already done in IPMESH2
        CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
        CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '    NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,DL,SE,XP,ERROR,*9999)
        CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,NKJE,NKEF,
     '    NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,NVJF,DF,PG,RG,SE,
     '    SF,WG,XA,XE,XG,XP,ERROR,*9999)
      ENDIF !ITYP3

C GMH 22/7/96 Calculate the global list of nodes
C Not for Voronoi nodes as they should have already been inputted rgb
      IF(JTYP14.NE.11) THEN
        CALL CALC_GLOBAL_NODES(NPNODE,ERROR,*9999)
      ENDIF

      IF(MAKE_GROUP)THEN
        IF(JTYP14.EQ.3.OR.JTYP14.EQ.4) THEN !Fractal tree
          CALL STRING_TRIM(GROUP_NAME,IBEG,IEND)
          STRING=GROUP_NAME(IBEG:IEND)
          i=0
          DO noelem=noelem_start+1,NEELEM(0,nr)
            i=i+1
            NELIST(i)=NEELEM(noelem,nr)
          ENDDO !noelem
          NELIST(0)=i
          CALL GRELEM_SUB(NELIST,STRING,.TRUE.,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,
     &      '('' Sorry, grouping not implemented for this mesh type'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL_BASE=.TRUE.
      CALL_ELEM=.TRUE.
      CALL_MESH=.TRUE.
      CALL_NODE=.TRUE.

      CALL EXITS('IPMESH')
      RETURN
 9998 ERROR=' '
 9999 CALL ERRORS('IPMESH',ERROR)
      CALL EXITS('IPMESH')
      RETURN 1
      END


