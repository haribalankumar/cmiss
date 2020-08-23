      SUBROUTINE DEMATE(CELL_ICQS_SPATIAL,CELL_ICQS_VALUE,
     '  CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,GRNGLIST,IBT,
     '  ICQS_SPATIAL,IDO,IICQS_SPATIAL,IRCQS_SPATIAL,ILPIN,ILTIN,INP,
     '  NBJ,NEELEM,NELIST,NENP,NENQ,NGLIST,NMBIN,NNB,NPLIST,NPNE,NPNODE,
     '  NQET,NQLIST,NQNE,NQS,NQXI,NRLIST,NW,NXI,NXLIST,CE,CELL_CP,
     '  CELL_RCQS_VALUE,CELL_YQS_VALUE,CGE,CIN,CP,CQ,RCQS_SPATIAL,
     '  XE,XIG,YG,YP,YQS,FIX,CELL_ICQS_NAMES,CELL_RCQS_NAMES,
     '  CELL_YQS_NAMES,STRING,ERROR,*)

C#### Subroutine: DEMATE
C###  Description:
C###    DEMATE defines materials with prompted input or from
C###    filename.ipmate (or filename.ipmatc if the cell qualifier is
C###    used).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER CELL_ICQS_SPATIAL(NQIM,NQVM),CELL_ICQS_VALUE(NQIM,NQVM),
     '  CELL_RCQS_SPATIAL(NQRM,NQVM),CELL_YQS_SPATIAL(NIQSM,NQVM),
     '  GRNGLIST(0:NEGM),IBT(3,NIM,NBFM),
     '  ICQS_SPATIAL(NQISVM,NQM),IDO(NKM,NNM,0:NIM,NBFM),
     '  IICQS_SPATIAL(0:NQISVM,NQVM),IRCQS_SPATIAL(0:NQRSVM,NQVM),
     '  ILPIN(NMM,NRM,NXM),ILTIN(NRM,NXM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),NGLIST(0:NGM),
     '  NMBIN(NMM,NRM,NXM),NNB(4,4,4,NBFM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),
     '  NQLIST(0:NQM),NQNE(NEQM,NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),
     '  NRLIST(0:NRM),NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NXLIST(0:NXM)
      REAL*8 CE(NMM,NEM,NXM),CELL_CP(NMQM,NPM),
     '  CELL_RCQS_VALUE(NQRM,NQVM),
     '  CELL_YQS_VALUE(NIQSM,NQVM),CGE(NMM,NGM,NEM,NXM),
     '  CIN(NMM,0:NGM,NNEPM),CP(NMM,NPM,NXM),
     '  CQ(NMM,NQM,NXM),RCQS_SPATIAL(NQRSVM,NQM),XE(NSM,NJM),
     '  XIG(NIM,NGM,NBM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),
     '  YQS(NIQSM,NQM)
      CHARACTER CELL_ICQS_NAMES(NQIM,NQVM)*(*),CELL_RCQS_NAMES(NQIM,
     '  NQVM)*(*),CELL_YQS_NAMES(NQIM,NQVM)*(*),ERROR*(*),STRING*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER CPBASIS,ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IFROMC,IPFILE,N3CO,nmq,no_nrlist,nq,nqsv,nqv,nr,nx,nxc,POINTS
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,CBBREV,CELL,FILIO,FIRST_TIME,GENER,
     '  MOUSE,COUPLING

      CALL ENTERS('DEMATE',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define material;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define the material properties of the solution domain. For
C###    example, isotropy, orthotropy, or anisotropy can be specified
C###    for an elasticity problem. If the problem is a heat conduction
C###    problem, the values of the coefficients are entered here, ie the
C###    source term and the values of the conductivity in the x,y and z
C###    directions are defined.  If the problem is a threshold
C###    activation problem, the wavespeeds in the fibre direction and in
C###    the transverse direction can be specified. The material
C###    parameters are read from or written to the file FILENAME.ipmate
C###    in the directory PATH.
C###  Parameter:      <cell>
C###    Specify that the cell material parameters are to be defined.
C###    If this is the call the material parameters are read from
C###    or writen to the file FILENAME.ipmatc
C###  Parameter:      <coupling>
C###    Specify that the coupling material parameters are to be defined.
C###  Parameter:      <constant/elements/nodes/grid/gauss>
C###    Specify if the material parameters are to be defined as spatially
C###    constant, at elements, at nodes, at grid points or at gauss points.
C###  Parameter:      <basis #[1]>
C###    Specify the basis type number if the material parameters are defined
C###    at node points.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<cell>'
        OP_STRING(3)=BLANK(1:15)//'<coupling>'
        OP_STRING(4)=BLANK(1:15)//'<(constant/elements/nodes/
     '                                  grid_pts/gauss_pts)[constant]>'
        OP_STRING(5)=BLANK(1:15)//'<basis #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEMATE',ERROR,*9999)
      ELSE
        IPFILE=2 !is input file version number on 7-Dec-1995
C GMH 7/12/95 1->2 Method of asking for element based values changed
        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        CALL ASSERT(CALL_EQUA,'>>Problem type not defined',ERROR,*9999)

        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'CELL',2,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_CELL,'>>Must define cell first',ERROR,*9999)
          CELL=.TRUE.
        ELSE
          CELL=.FALSE.
        ENDIF
        IF(CBBREV(CO,'COUPLING',3,noco+1,NTCO,N3CO)) THEN
          COUPLING=.TRUE.
        ELSE
          COUPLING=.FALSE.
        ENDIF

        IF(CBBREV(CO,'CONSTANT',3,noco+1,NTCO,N3CO)) THEN
          POINTS=1
        ELSEIF(CBBREV(CO,'ELEMENTS',2,noco+1,NTCO,N3CO))THEN
          POINTS=2
        ELSEIF(CBBREV(CO,'NODES',2,noco+1,NTCO,N3CO))THEN
          POINTS=3
        ELSEIF(CBBREV(CO,'GRID_PTS',2,noco+1,NTCO,N3CO))THEN
          POINTS=4
        ELSEIF(CBBREV(CO,'GAUSS_PTS',2,noco+1,NTCO,N3CO))THEN
          POINTS=5
        ELSE
          POINTS=1
        ENDIF

        IF(CBBREV(CO,'BASIS',3,noco+1,NTCO,N3CO)) THEN
          CPBASIS=IFROMC(CO(N3CO+1))
        ELSE
          CPBASIS=1
        ENDIF

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            IF(CELL) THEN
              CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'matc',
     '          STATUS,ERR,ERROR,*9999)
            ELSE
              CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'mate',
     '          STATUS,ERR,ERROR,*9999)
            ENDIF

            nmq=0
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)

              IF(CELL) THEN
                CALL IPMAT3_CELL(CELL_ICQS_SPATIAL,CELL_ICQS_VALUE,
     '            CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,IBT,
     '            ICQS_SPATIAL,IDO,IICQS_SPATIAL,IRCQS_SPATIAL,INP,
     '            NEELEM,NELIST,nmq,NPNE,NPNODE,NQET,NQLIST,NQNE,NQS,
     '            NQXI,nr,nx,CELL_CP,CELL_RCQS_VALUE,CELL_YQS_VALUE,
     '            RCQS_SPATIAL,XE,YQS,CELL_ICQS_NAMES,
     '            CELL_RCQS_NAMES,CELL_YQS_NAMES,ERROR,*9999)
C *** DPN 28 March 2000 - check for spatially varying Cm and/or Am
                DO nqv=1,CELL_NUM_VARIANTS
                  DO nqsv=1,IRCQS_SPATIAL(0,nqv)
                    IF(IRCQS_SPATIAL(nqsv,nqv).EQ.1) THEN
                      !Cm is spatially varying, so need to overwrite the
                      !default values from IPCELL with the spatially
                      !varying values for this variant
                      DO nq=1,NQT
                        IF (ICQS_SPATIAL(1,nq).EQ.nqv) CQ(1,nq,nx)=
     '                    RCQS_SPATIAL(nqsv,nq)
                      ENDDO !nq
                    ELSEIF(IRCQS_SPATIAL(nqsv,nqv).EQ.2) THEN
                      !Am is spatially varying, so need to overwrite the
                      !default values from IPCELL with the spatially
                      !varying values for this variant
                      DO nq=1,NQT
                        IF (ICQS_SPATIAL(1,nq).EQ.nqv) CQ(2,nq,nx)=
     '                    RCQS_SPATIAL(nqsv,nq)
                      ENDDO !nq
                    ENDIF
                  ENDDO !nqsv
                ENDDO !nqv
              ELSEIF(COUPLING)THEN
                CALL IPMATE_COUP(nr,nx,ERROR,*9999)
              ELSE
               write(*,*) 'calling ipmate', POINTS
                CALL IPMATE(CPBASIS,GRNGLIST,IBT,ICQS_SPATIAL,
     '            IDO,ILPIN(1,nr,nx),ILTIN(nr,nx),INP,IRCQS_SPATIAL,NBJ,
     '            NEELEM,NELIST,NENP,NENQ,NGLIST,NMBIN(1,nr,nx),NNB,
     '            NPLIST,NPNE,NPNODE,NQET,NQNE,NQS,NQXI,nr,NW(1,1,nx),
     '            nx,NXI,POINTS,CE(1,1,nx),CELL_RCQS_VALUE
     '            ,CGE(1,1,1,nx),CIN,CP(1,1,nx),CQ(1,1,nx),
     '            RCQS_SPATIAL,XE,
     '            XIG,YG,YP(1,1,nx),ALL_REGIONS,
     '            FIX(1,1,nx),ERROR,*9999)
              ENDIF

            ENDDO !no_nrlist
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF

        IF(IOTYPE.NE.3) THEN
          DO nr=1,NRT
            ASSEMBLE_GLOBAL(nr,nx)=.FALSE.
          ENDDO !nr
          CALL_UPGAUS_EIK=.FALSE.
        ENDIF
        IF(CELL) THEN
          CALL_CELL_MATE=.TRUE.
        ELSE
          CALL_MATE=.TRUE.
        ENDIF
      ENDIF

      CALL EXITS('DEMATE')
      RETURN
 9999 CALL ERRORS('DEMATE',ERROR)
      CALL EXITS('DEMATE')
      RETURN 1
      END


