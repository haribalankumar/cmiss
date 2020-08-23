      SUBROUTINE DEFACE(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,NKJE,
     '  NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,
     '  NVJF,DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,STRING,ERROR,*)

C#### Subroutine: DEFACE
C###  Description:
C###    DEFACE defines faces.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),
     '  NFFACE(0:NF_R_M,NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NKEF(0:4,16,6,NBFM),NLF(4,NFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM)
      REAL*8 DF(NFM),PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  SF(NSM,NBFM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IEND
      CHARACTER STATUS*3
      LOGICAL CALCU,FILIO,GENER,MOUSE

      CALL ENTERS('DEFACE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define face;c
C###  Description:
C###    This command calculates the element faces within all regions
C###    for a finite element mesh.

        OP_STRING(1)=STRING(1:IEND)//';c'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEFACE',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS('C',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)

        IF(CALCU) THEN
          CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '      NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,NVJF,
     '      DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,ERROR,*9999)
          CALL_FACE=.TRUE.
        ENDIF
        CALL_FACE=.TRUE.
      ENDIF

      CALL EXITS('DEFACE')
      RETURN
 9999 CALL ERRORS('DEFACE',ERROR)
      CALL EXITS('DEFACE')
      RETURN 1
      END


