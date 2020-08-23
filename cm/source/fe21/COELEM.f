      SUBROUTINE COELEM(IBT,IDO,INP,NBH,NBJ,NEELEM,NEL,NELIST,
     '  NENP,NHE,NKHE,NKJ,NKJE,NLL,NLLINE,NNL,NPL,NPLIST,NPNE,NPNODE,
     '  NVJE,NVJL,NW,DL,SE,XP,STRING,ERROR,*)

C#### Subroutine: COELEM
C###  Description:
C###    COELEM copies elements.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NHE(NEM,NXM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NLLINE(0:NL_R_M,0:NRM),NNL(0:4,12,NBFM),NPL(5,0:3,NLM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM),NW(NEM,3,NXM)
      REAL*8 DL(3,NLM),SE(NSM,NBFM,NEM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG3,IBEG4,IEND,IEND3,IEND4,IFROMC,N,N3CO,
     '  nb,nc,ne,NEOLD,NENEW,nh,nhx,nj,nk,nl,
     '  noelem,nolist,nn,np,NPTOT,nr,ns,nx
      CHARACTER CHAR3*3,CHAR4*4
      LOGICAL BELONGS,CBBREV,INLIST

      CALL ENTERS('COELEM',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(CHAR3(1:3),'(I3)') NET(1)
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
        NPTOT=0
        DO np=1,NPT(1) !and find whether np belongs to current elements
          BELONGS=.FALSE.
          DO ne=1,NET(1)
            IF(INLIST(np,NPNE(1,1,ne),NNT(1),N)) BELONGS=.TRUE.
          ENDDO
          IF(BELONGS) NPTOT=NPTOT+1
        ENDDO
        WRITE(CHAR4(1:4),'(I4)') NPTOT
        CALL STRING_TRIM(CHAR4,IBEG4,IEND4)
C---------------------------------------------------------------------
C#### Command: FEM copy elements
C###  Description:
C###    Copies elements to new element numbers
C###  Parameter: <from ELEMENT#s[1..0]>
C###  Parameter: <incr NODE#_INCR[0]>
C---------------------------------------------------------------------

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)
     '    //'<from ELEMENT#s[1..'//CHAR3(IBEG3:IEND3)//']>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<incr NODE#_INCR['//CHAR4(IBEG4:IEND4)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','COELEM',ERROR,*9999)
      ELSE
        nx=1 !temporary
        nc=1 !Temporary AJP 18-12-91

        IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NET(1),NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          NELIST(0)=0
          DO nr=1,NRT
            DO noelem=NELIST(0)+1,NELIST(0)+NEELEM(0,nr)
              NELIST(noelem)=NEELEM(noelem,nr)
            ENDDO
            NELIST(0)=NELIST(0)+NEELEM(0,nr)
          ENDDO
        ENDIF
        IF(CBBREV(CO,'INCREMENT',1,noco+1,NTCO,N3CO)) THEN
          NPTOT=IFROMC(CO(N3CO+1))
        ELSE
          NPTOT=0
          DO np=1,NPT(1) !& find whether np belongs to current elements
            BELONGS=.FALSE.
            DO ne=1,NET(1)
              IF(INLIST(np,NPNE(1,1,ne),NNT(1),N)) BELONGS=.TRUE.
            ENDDO
            IF(BELONGS) NPTOT=NPTOT+1
          ENDDO
        ENDIF
        CALL ASSERT(NET(1)+NELIST(0).LE.NEM,'>>NEM too small',
     '    ERROR,*9999)

        DO nolist=1,NELIST(0)
          NEOLD=NELIST(nolist)
          NENEW=NET(1)+nolist
          NEELEM(NEELEM(0,1)+nolist,1)=NENEW
          DO nj=1,NJ_LOC(0,0,0)
            nb=NBJ(nj,NEOLD)
            NBJ(nj,NENEW)=nb
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                NKJE(nk,nn,nj,NENEW)=NKJE(nk,nn,nj,NEOLD)
              ENDDO !nk
            ENDDO !nn
          ENDDO !nj
          NHE(NENEW,nx)=NHE(NEOLD,nx)
          DO nhx=1,NHE(NENEW,nx)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh,nc,NEOLD)
            NBH(nh,nc,NENEW)=nb
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
              NKHE(nk,nn,nh,NENEW)=NKHE(nk,nn,nh,NEOLD)
              ENDDO !nk
            ENDDO !nn
          ENDDO !nh
          NW(NENEW,1,nx)=NW(NEOLD,1,nx)
          DO nb=1,NBFT
            DO nn=1,NNT(nb)
              NPNE(nn,nb,NENEW)=NPNE(nn,nb,NEOLD)+NPTOT
C KAT 23Feb01: now handled by NKJE,NKHE
C              DO nk=1,NKT(nn,nb)
C                NKE(nk,nn,nb,NENEW)=NKE(nk,nn,nb,NEOLD)
C              ENDDO
            ENDDO
            DO ns=1,NST(nb)+NAT(nb)
              SE(ns,nb,NENEW)=SE(ns,nb,NEOLD)
            ENDDO
          ENDDO
          nb=NBJ(1,NENEW)
        ENDDO
        NET(1)=NET(1)+NELIST(0)
        NEELEM(0,1)=NEELEM(0,1)+NELIST(0)

C       set up basis fn, #derivs and mapping arrays for geometry
        CALL GLOBALJ(NBJ,NEELEM,NKJ,NPLIST,NPNE,NPNODE,
     '    ERROR,*9999)

C       define global line segments
        CALL LINSEG(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,NLL,NLLINE,
     '    NNL,NPL,NPNE,NPNODE,NVJE,NVJL,XP,ERROR,*9999)
        DO nl=1,NLT
          CALL ARCSCA(IDO,0,0,0,NBJ,NEL(0,nl),nl,
     '      NPL(1,0,nl),NPNE,NVJL(1,1,nl),DL,1.0D-6,XP,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('COELEM')
      RETURN
 9999 CALL ERRORS('COELEM',ERROR)
      CALL EXITS('COELEM')
      RETURN 1
      END


