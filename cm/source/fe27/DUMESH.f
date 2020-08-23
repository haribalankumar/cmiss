      SUBROUTINE DUMESH(IBT,IDO,INP,NBJ,NBJF,NEELEM,NEL,NELIST,
     '  NENP,NFF,NFFACE,NKJE,NKEF,NKJ,NLF,NLL,NLLINE,NNB,NNF,NNL,
     '  NPF,NP_INTERFACE,NPL,NPNE,NPNF,NPNODE,NRE,NRLIST,
     '  NVJE,NVJF,NVJL,NVJP,NXI,DF,DL,PG,RG,SE,SF,STRING,WG,XA,
     '  XE,XG,XP,ERROR,*)

C#### Subroutine: DUMESH
C###  Description:
C###    DUMESH duplicates elements.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NBJF(NJM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),
     '  NFFACE(0:NF_R_M,NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NKEF(0:4,16,6,NBFM),NKJ(NJM,NPM),
     '  NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),NNB(4,4,4,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NRLIST(0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),
     '  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 DF(NFM),DL(3,NLM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG6,IEND,IEND1,IEND6,IFROMC,N3CO,
     '  nb,ne,ne_new,nj,nk,nn,noelem,nonode,
     '  np,np_max,np_new,np_start,
     '  nr,ns,NEELEM_0_nr,NEELEM_0_0,
     '  NPNODE_0_nr,NPNODE_0_0,nv
      CHARACTER CHAR6*6
      LOGICAL ALL_REGIONS,CBBREV,GEN_ELEMENT,GEN_NODE

      CALL ENTERS('DUMESH',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
C        CHAR6=CFROMI(NPT(1)+1,'(I6)')
        WRITE(CHAR6,'(I6)') NPT(1)+1
        CALL STRING_TRIM(CHAR6,IBEG6,IEND6)

C---------------------------------------------------------------------

C#### Command: FEM duplicate mesh
C###  Parameter:         <(gen_element/gen_node/gen_both)[gen_both]>
C###    Specifies whether to generate new node and element files
C###    or both.
C###  Parameter:         <element (#s/all)[all]>
C###    Specifies the elements numbers to include. The 'all' command
C###    prompts all currently defined elements to be included.
C###  Parameter:         <node #[1]>
C###    Specifies the node numbers to be included
C###  Parameter:         <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Duplicates the current mesh - producing either nodes, elements
C###    (linked to original nodes) or both (a standalone mesh).  The
C###    node number given is the position from where new nodes are
C###    numbered.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)
     '    //'<(gen_element/gen_node/gen_both)[gen_both]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)
     '    //'<node #['//CHAR6(IBEG6:IEND6)//']>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DUMESH',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        nr=NRLIST(1)
        IF(CBBREV(CO,'NODE',1,noco+1,NTCO,N3CO)) THEN
          np_start=IFROMC(CO(N3CO+1))
        ELSE
          np_start=NPT(nr)+1
        ENDIF

        GEN_NODE=.TRUE.
        GEN_ELEMENT=.TRUE.
        IF(CBBREV(CO,'GEN_ELEMENT',5,noco+1,NTCO,N3CO)) THEN
          GEN_NODE=.FALSE.
        ELSE IF(CBBREV(CO,'GEN_NODE',5,noco+1,NTCO,N3CO)) THEN
          GEN_ELEMENT=.FALSE.
        ENDIF

        IF(GEN_NODE) THEN
          NPNODE_0_nr=NPNODE(0,nr)
          NPNODE_0_0 =NPNODE(0,0)
          np_max=0
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            np_new=np+np_start-NPNODE(1,nr)
            IF(nonode+NPNODE(0,nr).LE.NP_R_M) THEN
              NPNODE(nonode+NPNODE(0,nr),nr)=np_new
            ENDIF
            NPNODE_0_nr =NPNODE_0_nr +1 !#nodes in region nr
            NPNODE_0_0  =NPNODE_0_0  +1 !#nodes in all regions
            IF(np_new.LE.NPM) THEN
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np_new,.FALSE.,ERROR,*9999)
              DO nj=1,NJT
                NVJP(nj,np_new)=NVJP(nj,np)
                NKJ(nj,np_new)=NKJ(nj,np)
                DO nv=1,NVJP(nj,np_new)
                  DO nk=1,NKJ(nj,np_new)
                    XP(nk,nv,nj,np_new)=XP(nk,nv,nj,np)
                  ENDDO !nk
                ENDDO !nv
              ENDDO
            ENDIF !le.npm
            IF(np_new.GT.np_max) np_max=np_new
          ENDDO
        ENDIF !gen_node

        IF(GEN_ELEMENT) THEN
          NEELEM_0_nr=NEELEM(0,nr)
          NEELEM_0_0 =NEELEM(0,0)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            ne_new=ne+NET(nr)
            NEELEM(noelem+NEELEM(0,nr),nr)=ne_new
            NEELEM_0_nr =NEELEM_0_nr +1 !#elements in region nr
            NEELEM_0_0  =NEELEM_0_0  +1 !#elements in all regions
            nb=NBJ(1,ne)
            DO nn=1,NNT(nb)
              IF(GEN_NODE) THEN
                NPNE(nn,nb,ne_new)=NPNE(nn,nb,ne)+np_start-NPNODE(1,nr)
              ELSE
                NPNE(nn,nb,ne_new)=NPNE(nn,nb,ne)
              ENDIF
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                NVJE(nn,nb,nj,ne_new)=NVJE(nn,nb,nj,ne)
              ENDDO
            ENDDO
            DO ns=1,NST(nb)+NAT(nb)
              SE(ns,nb,ne_new)= SE(ns,nb,ne)
            ENDDO
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              nb=NBJ(nj,ne)
              NBJ(nj,ne_new)=nb
              DO nn=1,NNT(nb)
                DO nk=1,NKT(nn,nb)
                  NKJE(nk,nn,nj,ne_new)=NKJE(nk,nn,nj,ne)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF !gen_element

        IF(GEN_NODE) THEN
          NPT(nr)= np_max !highest node# in region nr
          NPT(0) = NPT(nr) !highest node# in all regions
          NPNODE(0,nr)=NPNODE_0_nr !#nodes in region nr
          NPNODE(0,0) =NPNODE_0_0 !#nodes in all regions
C         Update interface info
          CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
        ENDIF

        IF(GEN_ELEMENT) THEN
          CALL ASSERT(2*NET(nr).LE.NEM,'>>Increase NEM',
     '      ERROR,*9999)

C??? CS 3/2/98 Need to duplicate NRE as well.
          DO noelem=NET(nr)+1,2*NET(nr)
            NRE(noelem)=NRE(noelem-NET(nr))
          ENDDO
          NET(nr)= 2*NET(nr) !highest element# in region nr
          NET(0) =   NET(nr) !highest element# in all regions
          NEELEM(0,nr)=NEELEM_0_nr !#elements in region nr
          NEELEM(0,0) =NEELEM_0_0 !#elements in all regions

          CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)

          CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '      NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '      DL,SE,XP,ERROR,*9999)
          CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '      NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,NVJF,
     '      DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,ERROR,*9999)
        ENDIF !gen_element

        CALL_EQUA=.FALSE. !to ensure that NVJP    is redefined
        CALL_INIT=.FALSE. !to ensure that NHP     is redefined
        CALL_MATE=.FALSE. !to ensure that CP etc are redefined
      ENDIF
      CALL EXITS('DUMESH')

      RETURN
 9999 CALL ERRORS('DUMESH',ERROR)
      CALL EXITS('DUMESH')
      RETURN 1
      END


