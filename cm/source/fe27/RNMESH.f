      SUBROUTINE RNMESH(IBT,IDO,INP,NBJ,NBJF,NEELEM,
     '  NEL,NENP,NFF,NFFACE,NKJE,NKEF,NKJ,NLF,NLL,NLLINE,NNB,NNF,
     '  NNL,NPF,NPL,NPNE,NPNF,NPNODE,NRE,NRLIST,NVJE,NVJF,NVJL,
     '  NVJP,NXI,DF,DL,PG,RG,SE,SF,STRING,WG,XA,XE,XG,XP,
     '  ERROR,*)

C#### Subroutine: RNMESH
C###  Description:
C###    RNMESH renumbers mesh by adding a fixed offset to the node
C###    and element numbers.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!    Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NBJF(NJM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),
     '  NFFACE(0:NF_R_M,NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKEF(0:4,16,6,NBFM),
     '  NKJ(NJM,NPM),NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NPF(9,NFM),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     '  NPNODE(0:NP_R_M,0:NRM),NRE(NEM),NRLIST(0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),
     '  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 DF(NFM),DL(3,NLM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,n,N3CO,nb,ne,NE_NEW,NE_OFFSET,ne_old,njj,
     '  njl,nj,nk,nn,noelem,nonode,np,NP_NEW,NP_OFFSET,np_old,nr,ns,nv
      LOGICAL ALL_REGIONS,CBBREV

      CALL ENTERS('RNMESH',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM renumber mesh
C###  Parameter:        <node_offset #[0]>
C###    Specifies the 'node_offset' number which will be added to
C###    each node
C###  Parameter:        <element_offset #[0]>
C###    Specifies the 'element_offset' number which will be added to
C###    each element
C###  Parameter:        <region #[1]>
C###    Specify the element file region numbers to be defined.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<node_offset #[0]>'
        OP_STRING(3)=BLANK(1:15)//'<element_offset #[0]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','RNMESH',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
C CPB 24/11/93 NOTE: Regions not handled properly. Interconnections
C between regions etc are not maintained.
        nr=NRLIST(1)
        IF(CBBREV(CO,'NODE_OFFSET',1,noco+1,NTCO,N3CO)) THEN
          NP_OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          NP_OFFSET=0
        ENDIF
        IF(CBBREV(CO,'ELEMENT_OFFSET',1,noco+1,NTCO,N3CO)) THEN
          NE_OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          NE_OFFSET=0
        ENDIF
C CPB 25/11/95 - Initially shift all nodes and elements beyond the
C current set of nodes and elements (+their offsets) and then write
C them back into the correct locations to avoid overwriting arrays.
        CALL ASSERT(2*NPT(nr)+NP_OFFSET.LE.NPM,
     '    '>>Increase NPM s.b > 2*NPT+NODE_OFFSET',ERROR,*9999)
        CALL ASSERT(2*NET(nr)+NE_OFFSET.LE.NEM,
     '    '>>Increase NEM s.b > 2*NET+ELEMENT_OFFSET',ERROR,*9999)

C Shift the nodes and elements beyond the the current set of nodes and
C elements + their offsets and zero the old node and element locations

        IF(NP_OFFSET.NE.0) THEN
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            NP_NEW=np+NP_OFFSET+NPT(nr)
!           CALL NODE_CHANGE(np,.TRUE.,ERROR,*9999)
!           CALL NODE_CHANGE(NP_NEW,.TRUE.,ERROR,*9999)
            NPNODE(nonode,nr)=NP_NEW
            DO njl=1,3
              DO njj=1,NJ_LOC(njl,0,nr)
                nj=NJ_LOC(njl,njj,nr)
                DO nv=1,NVJP(nj,np)
                  DO nk=1,NKJ(nj,np)
                    XP(nk,nv,nj,NP_NEW)=XP(nk,nv,nj,np)
                    XP(nk,nv,nj,np)=0.0D0
                  ENDDO !nk
                ENDDO !nv
                NKJ(nj,NP_NEW)=NKJ(nj,np)
                NKJ(nj,np)=0
                NVJP(nj,NP_NEW)=NVJP(nj,np)
                NVJP(nj,np)=0
              ENDDO !njj
            ENDDO !njl
          ENDDO !nonode (np)
        ENDIF !no_offset

        IF(NE_OFFSET.NE.0) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NE_NEW=ne+NE_OFFSET+NET(nr)
            NEELEM(noelem,nr)=NE_NEW
            DO nb=1,NBFT
              DO nn=1,NNT(nb)
                NPNE(nn,nb,NE_NEW)=NPNE(nn,nb,ne) !adjusted later for np
                NPNE(nn,nb,ne)=0
              ENDDO !nn
              DO ns=1,NST(nb)+NAT(nb)
                SE(ns,nb,NE_NEW)=SE(ns,nb,ne)
                SE(ns,nb,ne)=0.0D0
              ENDDO !ns
            ENDDO !nb
            DO njl=1,3
              DO njj=1,NJ_LOC(njl,0,nr)
                nj=NJ_LOC(njl,njj,nr)
                nb=NBJ(nj,ne)
                NBJ(nj,NE_NEW)=nb
                NBJ(nj,ne)=0
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    NKJE(nk,nn,nj,NE_NEW)=NKJE(nk,nn,nj,ne)
                    NKJE(nk,nn,nj,ne)=0
                  ENDDO !nk
                ENDDO !nn
              ENDDO !njj
            ENDDO !njl
          ENDDO !noelem
        ENDIF !ne_offset

C Shift the nodes and elements back to their correct locations

        IF(NP_OFFSET.NE.0) THEN
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            NP_NEW=np-NPT(nr)
!           CALL NODE_CHANGE(NP_NEW,.TRUE.,ERROR,*9999)
            NPNODE(nonode,nr)=NP_NEW
            DO njl=1,3
              DO njj=1,NJ_LOC(njl,0,nr)
                nj=NJ_LOC(njl,njj,nr)
                DO nv=1,NVJP(nj,np)
                  DO nk=1,NKJ(nj,np)
                    XP(nk,nv,nj,NP_NEW)=XP(nk,nv,nj,np)
                  ENDDO
                ENDDO
                NKJ(nj,NP_NEW)=NKJ(nj,np)
                NVJP(nj,NP_NEW)=NVJP(nj,np)
              ENDDO !njj
            ENDDO !njl
          ENDDO
        ENDIF !np_offset

        IF(NE_OFFSET.NE.0) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NE_NEW=ne-NET(nr)
            NEELEM(noelem,nr)=NE_NEW
            DO nb=1,NBFT
              DO nn=1,NNT(nb)
                NPNE(nn,nb,NE_NEW)=NPNE(nn,nb,ne) !adjusted later for np
              ENDDO
              DO ns=1,NST(nb)+NAT(nb)
                SE(ns,nb,NE_NEW)=SE(ns,nb,ne)
              ENDDO
            ENDDO
            DO njl=1,3
              DO njj=1,NJ_LOC(njl,0,nr)
                nj=NJ_LOC(njl,njj,nr)
                nb=NBJ(nj,ne)
                NBJ(nj,NE_NEW)=nb
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    NKJE(nk,nn,nj,NE_NEW)=NKJE(nk,nn,nj,ne)
                  ENDDO
                ENDDO !nn
              ENDDO !njj
            ENDDO !njl
          ENDDO !noelem
        ENDIF !ne_offset

C Zero all the arrays above the shifted nodes and elements i.e.
C workspace arrays

        IF(NP_OFFSET.NE.0) THEN
          DO np=NPT(nr)+NP_OFFSET+1,2*NPT(nr)+NP_OFFSET
!         CALL NODE_CHANGE(np,.TRUE.,ERROR,*9999)
            DO njl=1,3
              DO njj=1,NJ_LOC(njl,0,nr)
                nj=NJ_LOC(njl,njj,nr)
                DO nv=1,NVJP(nj,np)
                  DO nk=1,NKJ(nj,np)
                    XP(nk,nv,nj,np)=0.0D0
                  ENDDO !nk
                ENDDO !nv
                NKJ(nj,np)=0
                NVJP(nj,np)=0
              ENDDO !njj
            ENDDO !njl
          ENDDO !np
        ENDIF !np_offset

        IF(NE_OFFSET.NE.0) THEN
          DO ne=NET(nr)+NE_OFFSET+1,2*NET(nr)+NE_OFFSET
            DO nb=1,NBFT
              DO nn=1,NNT(nb)
                NPNE(nn,nb,ne)=0
C KAT 2001-05-24: Why set this array for an element with no basis function?
C                DO nk=1,NKT(nn,nb)
C                  NKE(nk,nn,nb,ne)=0
C                ENDDO
              ENDDO
              DO ns=1,NST(nb)+NAT(nb)
                SE(ns,nb,ne)=0.0D0
              ENDDO
            ENDDO
            DO njl=1,3
              DO njj=1,NJ_LOC(njl,0,nr)
                nj=NJ_LOC(njl,njj,nr)
                NBJ(nj,ne)=0
              ENDDO !njj
            ENDDO !njl
          ENDDO !ne
        ENDIF !ne_offset

C Adjust the highest node numbers
        NPT(nr)= NPT(nr)+NP_OFFSET !highest node# in region nr
        NPT(0) = NPT(nr)           !highest node# in all regions

        NET(nr)= NET(nr)+NE_OFFSET !highest element# in region nr
        NET(0) = NET(nr)           !highest element# in all regions

C Update other arrays (added PJH 5may97)
        IF(NP_OFFSET.NE.0) THEN
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr) !new node#
            np_old=np-NP_OFFSET  !old node#
C!!! only works for non-overlapping sequences
            NENP(np,0,nr)=NENP(np_old,0,nr)
            NENP(np_old,0,nr)=0
            DO n=1,NENP(np,0,nr)
              NENP(np,n,nr)=NENP(np_old,n,nr) !adjusted later for ne
              NENP(np_old,n,nr)=0
            ENDDO !n
          ENDDO !nonode (np)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr) !new element#
            DO nb=1,NBFT
              DO nn=1,NNT(nb)
                NPNE(nn,nb,ne)=NPNE(nn,nb,ne)+NP_OFFSET
              ENDDO !nn
            ENDDO !nb
          ENDDO !noelem (ne)
        ENDIF !np_offset

        IF(NE_OFFSET.NE.0) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr) !new element#
            ne_old=ne-NE_OFFSET  !old element#
            NRE(ne)=NRE(ne_old)
            NRE(ne_old)=0
C!!! only works for non-overlapping sequences
            DO njl=1,3
              DO njj=1,NJ_LOC(njl,0,nr)
                nj=NJ_LOC(njl,njj,nr)
                DO nb=1,NBFT
                  DO nn=1,NNT(nb)
                    NVJE(nn,nb,nj,ne)=NVJE(nn,nb,nj,ne_old)
                    NVJE(nn,nb,nj,ne_old)=0
                  ENDDO !nn
                ENDDO !nb
              ENDDO !njj
            ENDDO !njl
          ENDDO !noelem (ne)
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr) !new node#
            DO n=1,NENP(np,0,nr)
              NENP(np,n,nr)=NENP(np,n,nr)+NE_OFFSET
            ENDDO !n
          ENDDO !nonode (np)
        ENDIF !ne_offset

C Recalculate the connectivity lines and faces
CC       .. average arc-length scaling means must calc NXI before SE
        CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)

        CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '    NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '    DL,SE,XP,ERROR,*9999)
        CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,NKJE,NKEF,
     '    NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,NVJF,
     '    DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,ERROR,*9999)

        CALL_BASE=.TRUE.
        CALL_ELEM=.TRUE.
        CALL_MESH=.TRUE.
        CALL_NODE=.TRUE.
        CALL_EQUA=.FALSE. !to ensure that NVJP    is redefined
        CALL_INIT=.FALSE. !to ensure that NHP     is redefined
        CALL_MATE=.FALSE. !to ensure that CP etc are redefined
      ENDIF

      CALL EXITS('RNMESH')
      RETURN
 9999 CALL ERRORS('RNMESH',ERROR)
      CALL EXITS('RNMESH')
      RETURN 1
      END


