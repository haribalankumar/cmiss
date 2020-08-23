      SUBROUTINE CHBASE(IBT,IDO,INP,NAN,NBJ,NBJF,NDET,NEELEM,NEL,NENP,
     '  NFF,NFFACE,NGAP,NKB,NKJE,NKEF,NLF,NLL,NLLINE,NNB,NNF,
     '  NNL,NPF,NPL,NPNE,NPNF,NPNODE,NRE,NSB,NVJE,NVJF,NVJL,
     '  DET,DF,DL,PG,RG,SE,SF,WG,XA,XE,XG,XIG,XP,STRING,ERROR,*)

C#### Subroutine: CHBASE
C###  Description:
C###   CHBASE changes basis parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
!     Parameter List

      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NDET(NBFM,0:NNM),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NKB(2,2,2,NNM,NBFM),
     '  NKJE(NKM,NNM,NJM,NEM),NKEF(0:4,16,6,NBFM),
     '  NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),NNB(4,4,4,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NSB(NKM,NNM,NBFM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM)
      REAL*8 DET(NBFM,0:NNM,NGM,6),DF(NFM),DL(3,NLM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,nb,ne,nj,nk,nn,noelem,nr,ns
      CHARACTER FILE*100
      LOGICAL ABBREV,CBBREV,UPDATE

      CALL ENTERS('CHBASE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM change bases
C###  Parameter:      <number BASIS#[1]>
C###    Defines the basis number to be changed
C###  Parameter:      <update>
C###  Description:
C###    Change an existing basis number to give a new definition of the
C###    basis.  Update recalculates elements, lines and faces etc.

        OP_STRING(1)=STRING(1:IEND)//' <number BASIS#[1]>'
        OP_STRING(2)=BLANK(1:15)//'<update>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHBASE',ERROR,*9999)
      ELSE
        IF(ABBREV(CO(noco+1),'NUMBER',1)) THEN
          nb=IFROMC(CO(noco+2))
        ELSE
          nb=1
        ENDIF
        IF(CBBREV(CO,'UPDATE',1,noco+1,NTCO,N3CO)) THEN
          UPDATE=.TRUE.
        ELSE
          UPDATE=.FALSE.
        ENDIF

        IOTYPE=1
        FILE='NEW'
        CALL STRING_TRIM(FILE,IBEG,IEND)
        CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.base','NEW',
     '    'DIRECT','FORMATTED',132,ERROR,*9999)
        CALL IPBASE(IBT,IDO,INP,NAN,nb,NDET,NGAP,NKB,NNB,NSB,
     '    DET,PG,WG,XIG,ERROR,*9999)
        CALL CLOSEF(IFILE,ERROR,*9999)

        IF(UPDATE) THEN
          DO nr=1,NRT
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nj=1,NJ_LOC(NJL_GEOM,0,0)
                IF(NBJ(nj,ne).EQ.nb) THEN
                  DO nn=1,NNT(nb)
                    DO nk=1,NKT(nn,nb)
                      NKJE(nk,nn,nj,ne)=nk
                    ENDDO
                  ENDDO
                ENDIF !nb
              ENDDO !nj
              DO ns=1,NST(nb)+NAT(nb)
                SE(ns,nb,ne)=1.d0
              ENDDO
            ENDDO
          ENDDO
          CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '      NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '      DL,SE,XP,ERROR,*9999)
          CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '      NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,NVJF,
     '      DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('CHBASE')
      RETURN
 9999 CALL ERRORS('CHBASE',ERROR)
      CALL EXITS('CHBASE')
      RETURN 1
      END


