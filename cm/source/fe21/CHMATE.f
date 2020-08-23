      SUBROUTINE CHMATE(CELL_ICQS_SPATIAL,CELL_ICQS_VALUE,
     '  CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,GRNGLIST,
     '  IBT,ICQS,ICQS_SPATIAL,IDO,
     '  IICQS_SPATIAL,IRCQS_SPATIAL,ILPIN,ILTIN,INP,
     '  NBJ,NEELEM,NELIST,NENQ,NGLIST,NMBIN,NPLIST,NPNE,NPNODE,NQLIST,
     '  NQNE,NQS,NQXI,NRLIST,NW,NXI,NXLIST,CE,CELL_CP,CELL_RCQS_VALUE,
     '  CELL_YQS_VALUE,CGE,CIN,CP,CQ,RCQS,RCQS_SPATIAL,XE,XIG,
     '  YG,YP,YQS,FIX,STRING,ERROR,*)

C#### Subroutine: CHMATE
C###  Description:
C###    CHMATE changes the material parameter input type or value

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER CELL_ICQS_SPATIAL(NQIM,NQVM),CELL_ICQS_VALUE(NQIM,NQVM),
     '  CELL_RCQS_SPATIAL(NQRM,NQVM),CELL_YQS_SPATIAL(NIQSM,NQVM),
     '  GRNGLIST(0:NEGM),IBT(3,NIM,NBFM),
     '  ICQS(NQIM),ICQS_SPATIAL(NQISVM,NQM),
     '  IDO(NKM,NNM,0:NIM,NBFM),
     '  IICQS_SPATIAL(0:NQISVM,NQVM),IRCQS_SPATIAL(0:NQRSVM,NQVM),
     '  ILPIN(NMM,NRM,NXM),ILTIN(NRM,NXM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NENQ(0:8,NQM),NGLIST(0:NGM),
     '  NMBIN(NMM,NRM,NXM),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NQLIST(0:NQM),NQNE(NEQM,NQEM),
     '  NQS(NEQM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),
     '  NW(NEM,3),NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM)
      REAL*8 CE(NMM,NEM,NXM),CELL_CP(NMQM,NPM),
     '  CELL_RCQS_VALUE(NQRM,NQVM),
     '  CELL_YQS_VALUE(NIQSM,NQVM),CGE(NMM,NGM,NEM,NXM),
     '  CIN(NMM,0:NGM,NNEPM),CP(NMM,NPM,NXM),
     '  CQ(NMM,NQM,NXM),RCQS(NQRM),
     '  RCQS_SPATIAL(NQRSVM,NQM),XE(NSM,NJM),XIG(NIM,NGM,NBM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),YQS(NIQSM,NQM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,il,INT_TEMP(1),
     '  N3CO,nb,ne,ne_mod,ng,ng_mod,noelem,nr,NTRL,nx,nxc
      REAL*8 max,percent,REAL_TEMP(1)
      LOGICAL ALL_REGIONS,CBBREV,DECREASE,GAUSS

      CALL ENTERS('CHMATE',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM change material parameter #[1] gauss
C###  Description:
C###    Specify that the material parameter be defined by Gauss
C###    points
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//' parameter #[1]'
        OP_STRING(2)=BLANK(1:15)//'<gauss>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change material parameter # decrease #

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHMATE',ERROR,*9999)
      ELSE

        OP_STRING(1)='Warning: Not fully implemented'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'GAUSS',2,noco+1,NTCO,N3CO)) THEN
          GAUSS=.TRUE.
        ELSE
          GAUSS=.FALSE.
        ENDIF

        IF(CBBREV(CO,'DECREASE',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),1,NTRL,REAL_TEMP,ERROR,*9999)
          percent=REAL_TEMP(1)
          DECREASE=.TRUE.
        ELSE
          DECREASE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'PARAMETER',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),1,NTRL,INT_TEMP,ERROR,*9999)
          il=INT_TEMP(1)
        ENDIF

        IF(GAUSS) THEN
          ILPIN(il,nr,nx)=4

          DO noelem=1,NEELEM(0,nr)
            ne = NEELEM(noelem,nr)
            nb=NBJ(1,ne)
            DO ng=1,NGT(nb)
              CGE(il,ng,ne,nx)=CE(il,ne,nx)
            ENDDO
          ENDDO
        ELSEIF(DECREASE) THEN
           max= 0.0d0
           DO noelem=1,NEELEM(0,nr)
            ne = NEELEM(noelem,nr)
            nb=NBJ(1,ne)
            DO ng=1,NGT(nb)
              IF(YG(il,ng,ne).GT.max) THEN
                ng_mod=ng
                ne_mod=ne
                max=YG(il,ng,ne)
              ENDIF
            ENDDO
          ENDDO
          CGE(il,ng_mod,ne_mod,nx)=CGE(il,ng_mod,ne_mod,nx)*
     '      (1.0d0-percent/100d0)
        ENDIF
      ENDIF

      CALL EXITS('CHMATE')
      RETURN
 9999 CALL ERRORS('CHMATE',ERROR)
      CALL EXITS('CHMATE')
      RETURN 1
      END


