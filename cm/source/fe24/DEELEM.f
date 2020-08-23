      SUBROUTINE DEELEM(IBT,IDO,INP,ISEG,ISELNO,ISLINE,ISLINO,
     '  ISNONO,MXI,NBJ,NBJF,NEELEM,NEL,NELIST,NENP,NFF,NFFACE,
     '  NKJE,NKEF,NKJ,NLF,NLL,NLLINE,NLLIST,NNB,NNF,NNL,NPF,NPL,NPLIST1,
     '  NPLIST2,NPNE,NPNF,NPNODE,NRE,NRLIST,NUNK,NVJE,NVJF,NVJL,NVJP,
     '  NXI,DF,DL,PG,RG,SE,SF,WG,XA,XE,XG,XP,ZC,ZP,CSEG,STRING,
     '  ERROR,*)

      
C#### Subroutine: DEELEM
C###  Description:
C###    DEELEM defines element parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'fbgr00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'linc00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'parameters.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISELNO(NWM,NEM),ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),
     '  ISNONO(NWM,NPM),MXI(2,NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKEF(0:4,16,6,NBFM),
     '  NKJ(NJM,NPM),NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),
     '  NLLIST(0:NLM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NPF(9,NFM),NPL(5,0:3,NLM),NPLIST1(0:NPM),NPLIST2(0:NPM),
     '  NPNE(NNM,NBFM,NEM),
     '  NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NRLIST(0:NRM),NUNK(NKM,NJM,NPM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),
     '  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 DF(NFM),DL(3,NLM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),ZC(NJM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,i,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IFROMC,INDEX,
     &  INDEX_OLD,INDEX_POLYLINE,INSTAT,IPFILE,IPICK,ISCIRC(27),ISEGM,
     &  iw,IWK(6),n1list,n2list,N3CO,nb,nbasis,NB1,nbb,NBLIST(0:4),nb_f,
     '  ne,NE_OFFSET,NETOLD,nf,nj,njj,njtype,nk,nl,nn,noelem,
     &  noelem_start,nolist,nonode,np,NP1,NP2,NP3,NP4,NP5,NP6,NP7,NP8,
     &  NP_OFFSET,no_nrlist,nr,ns,nss,NTIW,nx
      REAL*8 DX(3),RFROMC,SUM,X(3),XCENT,YCENT,Z(3),ZCENT
      CHARACTER FILE*(MXCH),GROUP_NAME*255,STATUS*3,TYPE*8
      LOGICAL ALL_REGIONS,BASIS_OVERRIDE,BOUNDARY,BOX,CALCU,CBBREV,
     &  FILIO,FIRST_TIME,GENER,MAKE_GROUP,MOUSE,ONE_DIMENSION,ASCII,
     &  BINARY,ENDFILE,TREE

      CALL ENTERS('DEELEM',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define elements;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    This command defines the number of elements, the number of
C###    coordinates for each element, the basis functions used for
C###    each geometric/fibre/field variable, and the node numbers
C###    defining the element.
C###    If scale factors are calculated from arc lengths, then this is
C###    done iteratively.  The maximum number of iterations in this
C###    process and the desired maximum error in scale factors and be
C###    specified by iterate and error respectively.
C###  Parameter:      <(geometry/fibre/field)[geometry]>
C###    Specify the type of element field to be defined. The file
C###    extensions are .ipelem for geometry, .ipelfb for fibres and
C###    .ipelfd for fields.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions. Any
C###    region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:      <(binary/ascii)[ascii]>
C###    Specify whether the element file is stored as a binary or
C###    ascii file.
C###  Parameter:      <iterate MAXIMUM[50]>
C###    Maximum number of iterations to be used in calculating arc
C###    length scale factors.
C###  Parameter:      <error MAXIMUM[CONVERG_TOL]>
C###    The maximum error allowable in calculating arc length scale
C###    factors.
C###  Parameter:    <elements (GROUP/#s/all)[all]>
C###    Specify either element group, element numbers or all elements
C###    to be exported.
C###  Parameter:    <basis (#)[1]>
C###    Overrides the basis function number in the file.
C###  Parameter:      <as GROUP_NAME>
C###    Allows you to specify a group name for the elements        
      

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<(geometry/fibre/field)[geometry]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<(binary/ascii)[ascii]>'
        OP_STRING(5)=BLANK(1:15)//'<1dimension>'
        WRITE(OP_STRING(6),'(15X,''<iterate MAXIMUM['',I3,'']>'')')
     '    LINC_ITMAX
        WRITE(OP_STRING(7),'(15X,''<error MAXIMUM['',D7.1,'']>'')')
     '    LINC_TOL
        OP_STRING(6)=BLANK(1:15)//'<as GROUP_NAME>'

        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define elements;m<;FILENAME[default]><;(PATH/example)[default]>
C###  Description:
C###    The elements are defined with the mouse using the specified
C###    basis function.  Use the left mouse button to 'pick' a node.
C###    Use the central mouse button to quit when you have finished
C###    defining elements.
C###  Parameter:      <on WS#[1]>
C###    WS is workstation window number.
C###  Parameter:      <basis #[1]>
C###    Define a basis number for the elements to be defined as.
C###  Parameter:      <rgb=RGB[blue]>
C###    Define the colour of the element lines that are to be drawn in
C###    in the window.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the element region numbers to be defined. The all value
C###    specifies all currently defined regions. Any region numbers
C###    contained in the file but are not specified will be skipped.

        OP_STRING(1)=STRING(1:IEND)//';m'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<on WS#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<basis #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[blue]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define elements;c
C###  Description:
C###    This command does not appear to do anything at the moment.
C###  Parameter:      <basis #s[1,1,1]>
C###
C###  Parameter:      <between NODES_1 and NODES_2>
C###
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the element region numbers to be defined. The all value
C###    specifies all currently defined regions. Any region numbers
C###    contained in the file but are not specified will be skipped.


        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(3)=BLANK(1:15)//'<basis #s[1,1,1]>'
        OP_STRING(4)=BLANK(1:15)//'<between NODES_1 and NODES_2>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEELEM',ERROR,*9999)
      ELSE

        CALL ASSERT(NBT.GT.0,'>>No basis functions are defined',
     '    ERROR,*9999)

        IPFILE=2 !is input file version number on 21May93
        CALL PARSE_QUALIFIERS('CDLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(MOUSE) THEN
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
          CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        ENDIF

        TREE=.FALSE. !set default for option on versions
        
        IF(CBBREV(CO,'1DIMENSION',2,noco+1,NTCO,N3CO)) THEN
          ONE_DIMENSION=.TRUE.
        ELSE
          ONE_DIMENSION=.FALSE.
        ENDIF
        IF(CBBREV(CO,'BASIS',3,noco+1,NTCO,N3CO)) THEN
          BASIS_OVERRIDE=.TRUE.
          nbasis=IFROMC(CO(N3CO+1))
        ELSE
          BASIS_OVERRIDE=.FALSE.
          nbasis=0
        ENDIF

        IF(CBBREV(CO,'OFFSET_NODES',8,noco+1,NTCO,N3CO)) THEN
          NP_OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          NP_OFFSET=0
        ENDIF

        IF(CBBREV(CO,'OFFSET_ELEMENTS',8,noco+1,NTCO,N3CO)) THEN
          NE_OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          NE_OFFSET=0
        ENDIF

        IF(CBBREV(CO,'ITERATE',2,noco+1,NTCO,N3CO)) THEN
          LINC_ITMAX=IFROMC(CO(N3CO+1))
        ENDIF
        IF(CBBREV(CO,'ERROR',2,noco+1,NTCO,N3CO)) THEN
          LINC_TOL=DABS(RFROMC(CO(N3CO+1)))
        ENDIF

C KSB 28 June 2004 Group elements 'as'
        IF(CBBREV(CO,'AS',2,noco+1,NTCO,N3CO)) THEN
          MAKE_GROUP=.TRUE.
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          GROUP_NAME=CO(N3CO+1)(IBEG:IEND)
          noelem_start=NEELEM(0,NRLIST(1))
        ELSE
          MAKE_GROUP=.FALSE.
        ENDIF
        
        BOX=.FALSE.

        IF(MOUSE) THEN
          CALL ASSERT(CALL_NODE,'>>Define nodes first',ERROR,*9999)
          IF(CBBREV(CO,'BASIS',1,noco+1,NTCO,N3CO)) THEN
            nb=IFROMC(CO(N3CO+1))
          ELSE
            nb=1
          ENDIF
          IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
            INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
          ELSE
            INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','RED')
          ENDIF

        ELSE IF(CALCU) THEN
          CALL ASSERT(CALL_NODE,'>>Define nodes first',ERROR,*9999)
          IF(CBBREV(CO,'BASIS',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),4,NBLIST(0),NBLIST(1),ERROR,*9999)
          ELSE
            NBLIST(0)=3
            DO nolist=1,NBLIST(0)
              NBLIST(nolist)=1
            ENDDO
          ENDIF
          NB1=NBLIST(1)
          IF(CBBREV(CO,'BETWEEN',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NPT(1),NPLIST1(0),NPLIST1(1),ERROR,
     '        *9999)
          ELSE
            IF(NIT(1).NE.1) THEN
              CO(noco+1)='?'
              GO TO 1
            ENDIF
          ENDIF
          IF(CBBREV(CO,'AND',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NPT(1),NPLIST2(0),NPLIST2(1),ERROR,
     '        *9999)
          ELSE
            IF(NIT(1).NE.1) THEN
              CO(noco+1)='?'
              GO TO 1
            ENDIF
          ENDIF
        ELSE IF(FILIO) THEN
          BINARY=.FALSE.
          ASCII=.TRUE.
          IF(CBBREV(CO,'BINARY',3,noco+1,NTCO,N3CO)) THEN
            BINARY=.TRUE.
            ASCII=.FALSE.

C           ..Check to see if binary is ok
            IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.4.OR.MOUSE.OR.CALCU) THEN
              ERROR='Binary elements only for r/w '
              GOTO 9999
            ENDIF
          ELSEIF(CBBREV(CO,'GEOMETRY',2,noco+1,NTCO,N3CO)) THEN
            CALL ASSERT(CALL_NODE,'>>Define nodes first',ERROR,*9999)
            TYPE='GEOMETRY'
          ELSE IF(CBBREV(CO,'FIBRE',3,noco+1,NTCO,N3CO)) THEN
            CALL ASSERT(CALL_FIBR,'>>Define fibres first',ERROR,*9999)
            TYPE='FIBRE'
            IPFILE=1 !version as of 10/3/97
          ELSE IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
            CALL ASSERT(CALL_FIEL,'>>Define field first',ERROR,*9999)
            TYPE='FIELD'
            IPFILE=1 !version as of 10/3/97
            IF(CBBREV(CO,'TREE',4,noco+1,NTCO,N3CO)) THEN
              TREE=.TRUE.
            ELSE
              TREE=.FALSE.
            ENDIF

          ELSE
            CALL ASSERT(CALL_NODE,'>>Define nodes first',ERROR,*9999)
            TYPE='GEOMETRY'
          ENDIF
        ENDIF

        IF(FILIO) THEN
          IF(ASCII) THEN
            FIRST_TIME=.TRUE.
            BACKUP=.FALSE.
            DO WHILE(FIRST_TIME.OR.BACKUP)
              FIRST_TIME=.FALSE.
              IF(TYPE(1:8).EQ.'GEOMETRY') THEN
                CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'elem',STATUS,
     '            ERR,ERROR,*9999)
                njtype=NJL_GEOM
              ELSE IF(TYPE(1:5).EQ.'FIBRE') THEN
                CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'elfb',STATUS,
     '            ERR,ERROR,*9999)
                njtype=NJL_FIBR
              ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
                CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'elfd',STATUS,
     '            ERR,ERROR,*9999)
                njtype=NJL_FIEL
              ENDIF

              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)

                CALL IPELEM(IBT,nbasis,NBJ,NEELEM,NELIST,NENP,NE_OFFSET,
     &            NKJE,NPNE,NPNODE,NP_OFFSET,nr,NRLIST,NRE,NVJE,NVJP,
     &            BASIS_OVERRIDE,TYPE,TREE,ERROR,*9999)

C LKC 15-APR-1999 breaking up the asserts
C Note : NET(nr) will always be less than NET(0) so don't need to check
C                CALL ASSERT(NET(nr).LE.NEM.AND.NLT.LE.NLM,
C     '            '>>NEM or NLM too small',ERROR,*9999)
C                CALL ASSERT(NET(0).LE.NEM.AND.NLT.LE.NLM,
C     '            '>>NEM or NLM too small',ERROR,*9999)

                CALL ASSERT(NET(0).LE.NEM,
     '            '>>Increase NEM ',ERROR,*9999)
                WRITE(OP_STRING,'(''>>Increase NLM '',I6)') NLT
                CALL ASSERT(NLT.LE.NLM,
     '            '>>Increase NLM',ERROR,*9999)

                IF(NPNODE(0,nr).GT.0) THEN !Nodes have been defined

C               Give warning if #node derivs doesn't match NKT(0,nb)
                  DO njj=1,NJ_LOC(njtype,0,nr)
                    nj=NJ_LOC(njtype,njj,nr)
                    nb=NBJ(nj,NEELEM(1,nr))
                    IF(NKT(0,nb).GT.NKJ(nj,NPNODE(1,nr))) THEN
                      WRITE(OP_STRING,'('' >>Warning: May need to '
     '                  //'update node derivatives'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDDO !nj
                ENDIF
              ENDDO !nr
              CALL CLOSEF(IFILE,ERROR,*9999)
              IF(BACKUP) IOTYPE=2
            ENDDO
            CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

            IF(IOTYPE.NE.3) THEN
              IF(NPT(0).GT.0.AND.TYPE(1:8).EQ.'GEOMETRY') THEN
                IF(ONE_DIMENSION)THEN
                  nb=NBJ(1,NEELEM(1,nr))
                  CALL NENXI_1D(nb,NEELEM,NENP,NPNE,nr,NXI,ERROR,*9999)
                  CALL LINCAL_1D(NBJ,NEELEM,NKJE,nr,NRE,NVJE,
     '              SE,ERROR,*9999)
                ELSE
                  CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
                  CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '              NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '              DL,SE,XP,ERROR,*9999)
                  CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '              NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,
     '              NRE,NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,XE,
     '              XG,XP,ERROR,*9999)
                ENDIF !ITYP3

              ELSE !fibre or field
C             Set up NBJF
                DO nf=1,NFT
                  ne=NPF(6,nf)
                  nr=NRE(ne) !perhaps should check both adjacent elements
                  DO njj=1,NJ_LOC(njtype,0,nr)
                    nj=NJ_LOC(njtype,njj,nr)
                    nb=NBJ(nj,ne) !is basis # of parent elem for nj
                    IF(nb.NE.0) THEN
                      IF(NIT(nb).EQ.2) THEN !face is the element
                        NBJF(nj,nf)=nb
                      ELSE !select face basis based on parent element
                        CALL FIND_FACE_BASIS(IBT,nb,nb_f,NPF(8,nf),NNF,
     '                    ERROR,*9999)
                        IF(nb_f.EQ.0) THEN
                          WRITE(ERROR,'('' >>No face basis found for '
     '                      //'face '',I5)') nf
                          GOTO 9999
                        ENDIF
C                     Have found basis nb that matches the face basis
                        NBJF(nj,nf)=nb_f
                      ENDIF !NIT=2
                    ELSE !nb=0 so no basis defined for this nj
                      NBJF(nj,nf)=0
                    ENDIF
                  ENDDO !njj
                ENDDO !nf
              ENDIF !geometry
              CALL GLOBALJ(NBJ,NEELEM,NKJ,NPLIST1,NPNE,NPNODE,
     '          ERROR,*9999)
C           23-6-92 AJP.  Added to ensure NKJ is set correctly even
C           when the derivatives do not need updating. This must be
C           done after call to LINCAL above.
            ENDIF !IOTYPE.NE.3
          ELSEIF(BINARY) THEN

            IF(IOTYPE.EQ.2) THEN !Read the file
              CALL IOELEM(IOFILE1,NBJ,NEELEM,NKJE,NPNE,NRE,NRLIST,NVJE,
     '          'READ','BINARY',FILE00,'OPEN',ENDFILE,ERROR,*9999)
              CALL IOELEM(IOFILE1,NBJ,NEELEM,NKJE,NPNE,NRE,NRLIST,NVJE,
     '          'READ','BINARY',FILE00,'ELEMDATA',ENDFILE,ERROR,*9999)
              CALL IOELEM(IOFILE1,NBJ,NEELEM,NKJE,NPNE,NRE,NRLIST,NVJE,
     '          'CLOSE','BINARY',FILE00,'SEEYA',ENDFILE,ERROR,*9999)
            ELSEIF(IOTYPE.EQ.3) THEN !Write the file
              CALL IOELEM(IOFILE1,NBJ,NEELEM,NKJE,NPNE,NRE,NRLIST,NVJE,
     '          'WRITE','BINARY',FILE00,'OPEN',ENDFILE,ERROR,*9999)
              CALL IOELEM(IOFILE1,NBJ,NEELEM,NKJE,NPNE,NRE,NRLIST,NVJE,
     '          'WRITE','BINARY',FILE00,'ELEMDATA',ENDFILE,ERROR,*9999)
              CALL IOELEM(IOFILE1,NBJ,NEELEM,NKJE,NPNE,NRE,NRLIST,NVJE,
     '          'CLOSE','BINARY',FILE00,'SEEYA',ENDFILE,ERROR,*9999)
            ELSE
              ERROR='> Incorrect option in binary ioelem'
              GOTO 9999
            ENDIF

C           ..Face, line calculations

            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)

              CALL ASSERT(NET(nr).LE.NEM.AND.NLT.LE.NLM,
     '          '>>NEM or NLM too small',ERROR,*9999)
              CALL ASSERT(NET(0).LE.NEM.AND.NLT.LE.NLM,
     '          '>>NEM or NLM too small',ERROR,*9999)

              IF(NPNODE(0,nr).GT.0) THEN !Nodes have been defined

C               Give warning if #node derivs doesn't match NKT(0,nb)
                DO njtype=1,3
                  DO njj=1,NJ_LOC(njtype,0,nr)
                    nj=NJ_LOC(njtype,njj,nr)
                    nb=NBJ(nj,NEELEM(1,nr))
                    IF(NKT(0,nb).GT.NKJ(nj,NPNODE(1,nr))) THEN
                      WRITE(OP_STRING,'('' >>Warning: May need to '
     '                  //'update node derivatives'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDDO !nj
                ENDDO !njtype
              ENDIF
            ENDDO !nr
            CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

            IF(IOTYPE.NE.3) THEN
              IF(ONE_DIMENSION)THEN
                nb=NBJ(1,NEELEM(1,nr))
                CALL NENXI_1D(nb,NEELEM,NENP,NPNE,nr,NXI,ERROR,*9999)
              ELSE
C               Determine NXI before SE for average arc-length scaling
                CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
C               Determine element connectivity
                IF((USE_VORONOI.EQ.0).OR.nr.EQ.1)THEN
C               Calculate line & face info if node set complete
                  CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '              NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '              DL,SE,XP,ERROR,*9999)
                  CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '              NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,
     '              NRE,NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,XE,
     '              XG,XP,ERROR,*9999)
                ELSE
                  DO nr=1,NRLIST(0)
                    DO noelem=1,NEELEM(0,nr)
                      ne=NEELEM(noelem,nr)
                      DO nb=1,NBFT
                        DO ns=1,NST(nb)+NAT(nb)
                          SE(ns,nb,ne)=1.d0
                        ENDDO !ns
                      ENDDO !nb
                    ENDDO !noelem
                  ENDDO !nr
                ENDIF

C             Set up NBJF
                DO nf=1,NFT
                  ne=NPF(6,nf)
                  nr=NRE(ne)
                  DO njtype=2,3 !fibre/field
                    DO njj=1,NJ_LOC(njtype,0,nr)
                      nj=NJ_LOC(njtype,njj,nr)
                      nb=NBJ(nj,ne) !is basis # of parent elem for nj
                      IF(nb.NE.0) THEN
                        IF(NIT(nb).EQ.2) THEN !face is the element
                          NBJF(nj,nf)=nb
                        ELSE !select face basis based on parent elem
                          CALL FIND_FACE_BASIS(IBT,nb,nb_f,NPF(8,nf),
     '                      NNF,ERROR,*9999)
                          IF(nb_f.EQ.0) THEN
                            WRITE(ERROR,'('' >>No face basis found '
     '                        //'for face '',I5)') nf
                            GOTO 9999
                          ENDIF
C                        Have found basis nb that matches the face basis
                          NBJF(nj,nf)=nb_f
                        ENDIF !NIT=2
                      ELSE !nb=0 so no basis defined for this nj
                        NBJF(nj,nf)=0
                      ENDIF
                    ENDDO !njj
                  ENDDO !njtype
                ENDDO !nf
              ENDIF
              CALL GLOBALJ(NBJ,NEELEM,NKJ,NPLIST1,NPNE,NPNODE,
     '          ERROR,*9999)
            ENDIF !IOTYPE.NE.3

          ELSE
            ERROR='Incorrect option in file i/o'
            GOTO 9999
          ENDIF

        ELSE IF(MOUSE) THEN

C LKC 12-JUL-1999 need USE_GRAPHICS for ZC now
          CALL ASSERT(USE_GRAPHICS.EQ.1,'>> Set USE_GRAPHICS to 1',
     '      ERROR,*9999)

          nr=1 !Temporary PJH 28Feb95
          CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'elem',STATUS,
     '      ERR,ERROR,*9999)
          iw=IWK(1)
          NP1=NPNODE(1,nr)
          CALL ASSERT(ISNONO(iw,NP1).GT.0,
     '      '>>Node segments not defined',ERROR,*9999)
          CALL ACWK(iw,0,ERROR,*9999)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(ISEG(ISELNO(1,ne)).EQ.1) THEN
              CALL VISIB(iw,ISEG,ISNONO(iw,np1),'VISIBLE',ERROR,*9999)
            ENDIF
          ENDDO !noelem
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            CALL DETECT(iw,ISEG,ISNONO(iw,np),'DETECTABLE',ERROR,*9999)
          ENDDO
C cpb 14/3/94 Clear any existing elements unless doing a define;add
          CALL DAWK(iw,0,ERROR,*9999)
          IF(NET(nr).GT.0) THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO i=1,NTIW
                iw=IWK(i)
                IF(ISEG(ISELNO(iw,ne)).GT.0) THEN
                  CALL ACWK(iw,1,ERROR,*9999)
                  CALL DELETE_SEGMENT(ISELNO(iw,ne),ISEG,iw,ERROR,*9999)
                  CALL DAWK(iw,1,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO !noelem
          ENDIF
          DO I=1,27
            ISCIRC(I)=0
          ENDDO !I
          CALL ACWK(iw,0,ERROR,*9999)
          NET(nr)=0
          NEELEM(0,nr)=0
          NETOLD=NET(nr)
          ne=0
          noelem=0
          INSTAT=1
          DO WHILE(INSTAT.EQ.1)
C ***       Create element segment
            WRITE(OP_STRING,'('' >>Pick '',I2,'
     '        //''' nodes on '',I2,'':'')')
     '        NNT(nb),iw
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            noelem=noelem+1
            ne=ne+1
            CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
            XCENT=0.0D0
            YCENT=0.0D0
            ZCENT=0.0D0
            DO nn=1,NNT(nb)
              CALL PICK(iw,'REQUEST',INSTAT,ISEGM,IPICK,ERROR,*9999)
              IF(INSTAT.EQ.1) THEN
                np=IFROMC(CSEG(ISEGM)(53:57))
                WRITE(OP_STRING,'('' Node '',I5)') np
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                NPNE(nn,nb,ne)=np
                DO nj=1,NJT
                  X(nj)=XP(1,1,nj,np)
                ENDDO
                CALL OPEN_SEGMENT(ISCIRC(nn),ISEG,iw,'NODECIRC',INDEX,
     '            INDEX_OLD,1,1,CSEG,ERROR,*9999)
                CALL XZ(ITYP10(1),X,Z)
                CALL CIRCLE('FILL AREA',iw,0.007D0*DBLE(DIAG),Z,10,
     '            ERROR,*9999)
                CALL CLOSE_SEGMENT(ISCIRC(nn),iw,ERROR,*9999)
                XCENT=XCENT+Z(1)
                YCENT=YCENT+Z(2)
C TVK 16July1999: Fix fem de elem;m bug for 2D
                IF(NJT.EQ.3) THEN
                   ZCENT=ZCENT+Z(3)
                ELSE
                  ZCENT=0
                ENDIF
              ELSE
                NEELEM(0,nr)=noelem-1
                GO TO 210
              ENDIF
            ENDDO !nn
 210        IF(INSTAT.EQ.1) THEN
              XCENT=XCENT/DBLE(NNT(nb))
              YCENT=YCENT/DBLE(NNT(nb))
              ZCENT=ZCENT/DBLE(NNT(nb))
              ZC(1,ne)=XCENT
              ZC(2,ne)=YCENT
              ZC(3,ne)=ZCENT
              NRE(ne)=nr
              NEELEM(noelem,nr)=ne
              NEELEM(0,nr)=NEELEM(0,nr)+1
              DO nj=1,NJT
                NBJ(nj,ne)=nb
              ENDDO !nj

C LKC 12-APR-98 Initialise Boundary to ????
              BOUNDARY=.FALSE.
              CALL SGELEM(INDEX,ISEG,ISELNO(iw,ne),iw,MXI(1,ne),
     '          NBJ(1,ne),ne,NLL(1,ne),NPL,nr,BOUNDARY,BOX,
     '          CSEG,DL,XP,ZC(1,ne),ERROR,*9999)
              DO nn=1,NNT(nb)
                CALL DELETE_SEGMENT(ISCIRC(nn),ISEG,iw,ERROR,*9999)
              ENDDO

C ***         Define deriv#s, version#s & scaling factors
              DO nn=1,NNT(nb)
C                DO nk=1,NKT(nn,nb)
C                  NKE(nk,nn,nb,ne)=nk
C                ENDDO
                DO nj=1,NJT
                  NVJE(nn,nb,nj,ne)=1
                  DO nk=1,NKT(nn,nb)
                    NKJE(nk,nn,nj,ne)=nk
                  ENDDO
                ENDDO
              ENDDO
              DO ns=1,NST(nb)+NAT(nb)
                SE(ns,nb,ne)=1.0D0
              ENDDO
            ELSE
              DO nn=1,NNT(nb)
                IF(ISCIRC(nn).NE.0) THEN
                  IF(ISEG(ISCIRC(nn)).GT.0) THEN
                    CALL DELETE_SEGMENT(ISCIRC(nn),ISEG,iw,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDDO !while instat=1

          NET(nr)=NEELEM(0,nr)
          IF(NET(nr).GT.NET(0)) NET(0)=NET(nr)
          IF(NEELEM(0,nr).GT.NEELEM(0,0)) NEELEM(0,0)=NEELEM(0,nr)

          IF(NKT(0,nb).GT.1.AND.NIT(nb).EQ.1) THEN
C ***       Calculate derivatives for cubic Hermite basis
            DO ne=NETOLD+1,NET(nr)
              NP1=NPNE(1,nb,ne)
              NP2=NPNE(2,nb,ne)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
              SUM=0.0D0
              DO nj=1,NJT
                DX(nj)=XP(1,1,nj,NP2)-XP(1,1,nj,NP1)
                SUM=SUM+DX(nj)**2
              ENDDO
              SUM=DSQRT(SUM)
              DO nj=1,NJT
C KAT 1Sep98:  IF not necessary
C                IF(JTYP2.EQ.0.OR.JTYP2.EQ.2.OR.JTYP2.EQ.3) THEN
                  XP(2,1,nj,NP1)=DX(nj)/SUM
                  XP(2,1,nj,NP2)=DX(nj)/SUM
C                ENDIF
              ENDDO
            ENDDO
          ENDIF !nkt,nit

C ***     Define other basis function information
          DO nbb=1,NBFT
            IF(nbb.NE.nb.AND.NNT(nbb).EQ.NNT(nb)) THEN
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                DO nn=1,NNT(nbb)
                  NPNE(nn,nbb,ne)=NPNE(nn,nb,ne)
                ENDDO
                DO nss=1,NST(nbb)+NAT(nbb)
                  SE(nss,nbb,ne)=1.0D0
                ENDDO
              ENDDO
            ENDIF
          ENDDO !nbb

C         Calculate elements surrounding nodes
          CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

C         Calculate line & face info if node set complete
          CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '      NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '      DL,SE,XP,ERROR,*9999)
          CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '      NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,
     '      NRE,NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,XE,
     '      XG,XP,ERROR,*9999)

C         Set up basis fn, #derivs and mapping arrays for geometry
          CALL GLOBALJ(NBJ,NEELEM,NKJ,NPLIST1,NPNE,NPNODE,
     '      ERROR,*9999)

          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            CALL DETECT(iw,ISEG,ISNONO(iw,np),'UNDETECTABLE',ERROR,
     '        *9999)
          ENDDO
          CALL DAWK(iw,0,ERROR,*9999)
          CALL ASSERT(NET(nr).LE.NEM.AND.NLT.LE.NLM,
     '      '>>NEM or NLM too small',ERROR,*9999)
          NLLIST(0)=NLT
          DO nolist=1,NLLIST(0)
            NLLIST(nolist)=nolist
          ENDDO
          NTLINE=NTLINE+1
          nx=1 !temporary
          DO iw=1,2*NJT-3
            CALL ACWK(iw,1,ERROR,*9999)
            CALL SGLINE(INDEX,ISEG,ISLINE(iw,NTLINE),ISLINO(iw),iw,
     '        NLLIST,NTLINE,NPL,nr,nx,CSEG,'UNDEFORMED',DL,'SOLID',
     '        .TRUE.,XP,ZP,ERROR,*9999)
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO
          IF(FBGRAF) THEN
            DO iw=5,6
              CALL ACWK(iw,1,ERROR,*9999)
              CALL SGLINE(INDEX,ISEG,ISLINE(iw,1),ISLINO(iw),iw,NLLIST,
     '          1,NPL,nr,nx,CSEG,'UNDEFORMED',DL,'SOLID',.TRUE.,
     '          XP,ZP,ERROR,*9999)
              CALL DAWK(iw,1,ERROR,*9999)
            ENDDO
          ENDIF

C         Write ipelem file (iotype=3 is set in PARSE_QUALIFIERS)
          CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'elem',
     '      STATUS,ERR,ERROR,*9999)
          TYPE='GEOMETRY'
          CALL IPELEM(IBT,nbasis,NBJ,NEELEM,NELIST,NENP,NE_OFFSET,NKJE,
     &      NPNE,NPNODE,NP_OFFSET,nr,NRLIST,NRE,NVJE,NVJP,
     &      BASIS_OVERRIDE,TYPE,TREE,ERROR,*9999)
          CALL CLOSEF(IFILE,ERROR,*9999)

          CALL_LINE=.TRUE.

        ELSE IF(CALCU) THEN
          nr=1 !Temporary PJH 28Feb95
          IF(NIT(NB1).EQ.1) THEN !calc 1D elems along a line of nodes
            DO np=1,NPT(nr)-1
              ne=np
              NPNE(1,NB1,ne)=np
              NPNE(2,NB1,ne)=np+1
              NEELEM(ne,nr)=ne
              DO nj=1,NJT
                nb=NBLIST(nj)
                NBJ(nj,ne)=nb
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    NKJE(nk,nn,nj,ne)=nk
                  ENDDO
                ENDDO
              ENDDO
              DO nolist=1,NBLIST(0)
                nb=NBLIST(nolist)
C KAT 23Feb01: now handled by NKJE
C                DO nn=1,NNT(nb)
C                  DO nk=1,NKT(nn,nb)
C                    NKE(nk,nn,nb,ne)=nk
C                  ENDDO
C                ENDDO
                DO ns=1,NST(nb)+NAT(nb)
                  SE(ns,nb,ne)=1.0D0
                ENDDO
              ENDDO
              IF(NKT(0,NB1).GT.1) THEN !Calc derivs for cubic Herm basis
                NP1=NPNE(1,NB1,ne)
                NP2=NPNE(2,NB1,ne)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
                SUM=0.0D0
                DO nj=1,NJT
                  DX(nj)=XP(1,1,nj,NP2)-XP(1,1,nj,NP1)
                  SUM=SUM+DX(nj)**2
                ENDDO
                SUM=DSQRT(SUM)
                DO nj=1,NJT
                  XP(2,1,nj,NP1)=DX(nj)/SUM
                  XP(2,1,nj,NP2)=DX(nj)/SUM
                ENDDO
              ENDIF
            ENDDO
            NET(nr)=NPT(nr)-1
            NEELEM(0,nr)=NET(nr)

          ELSE IF(NIT(NB1).EQ.3) THEN !calculate 3D elements
            ne=0
            DO n1list=1,NPLIST1(0) !loops over 1st layer of nodes
              NP1=NPLIST1(n1list)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' NP1='',I5)') NP1
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              NP2=0
              DO nl=1,NLT !to find node NP2 at Xi(1)=1 along from NP1
                IF(NPL(1,0,nl).EQ.1.AND.NPL(2,1,nl).EQ.NP1) THEN
                  NP2=NPL(3,1,nl)
                  GO TO 20
                ENDIF
              ENDDO
 20           NP3=0
              DO nl=1,NLT !to find node NP3 at Xi(2)=1 along from NP1
                IF(NPL(1,0,nl).EQ.2.AND.NPL(2,1,nl).EQ.NP1) THEN
                  NP3=NPL(3,1,nl)
                  GO TO 21
                ENDIF
              ENDDO
 21           IF(NP2.GT.0.AND.NP3.GT.0) THEN
                ne=ne+1
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ne='',I5)') ne
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                DO nolist=1,NBLIST(0)
                  nb=NBLIST(nolist)
                  NPNE(1,nb,ne)=NP1
                  NPNE(2,nb,ne)=NP2
                  NPNE(3,nb,ne)=NP3
                ENDDO
                NP4=0
                DO nl=1,NLT !to find node NP4 at Xi(2)=1 along from NP2
                  IF(NPL(1,0,nl).EQ.2.AND.NPL(2,1,nl).EQ.NP2) THEN
                    NP4=NPL(3,1,nl)
                    GO TO 22
                  ENDIF
                ENDDO
 22             IF(NP4.GT.0) THEN
                  DO nolist=1,NBLIST(0)
                    nb=NBLIST(nolist)
                    NPNE(4,nb,ne)=NP4
                  ENDDO
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' np1..4='',4I6)')
     '                NP1,NP2,NP3,NP4
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSE
                  ERROR=' >>Node for nn=4 not found'
                  GO TO 9999
                ENDIF
              ENDIF
            ENDDO

            ne=0
            DO n2list=1,NPLIST2(0) !loops over 2nd layer of nodes
              NP5=NPLIST2(n2list)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' NP5='',I5)') NP5
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              NP6=0
              DO nl=1,NLT !to find node NP6 at Xi(1)=1 along from NP5
                IF(NPL(1,0,nl).EQ.1.AND.NPL(2,1,nl).EQ.NP5) THEN
                  NP6=NPL(3,1,nl)
                  GO TO 30
                ENDIF
              ENDDO
 30           NP7=0
              DO nl=1,NLT !to find node NP7 at Xi(2)=1 along from NP5
                IF(NPL(1,0,nl).EQ.2.AND.NPL(2,1,nl).EQ.NP5) THEN
                  NP7=NPL(3,1,nl)
                  GO TO 31
                ENDIF
              ENDDO
 31           IF(NP6.GT.0.AND.NP7.GT.0) THEN
                ne=ne+1
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ne='',I5)') ne
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                DO nolist=1,NBLIST(0)
                  nb=NBLIST(nolist)
                  NPNE(5,nb,ne)=NP5
                  NPNE(6,nb,ne)=NP6
                  NPNE(7,nb,ne)=NP7
                ENDDO
                NP8=0
                DO nl=1,NLT !to find node NP8 at Xi(2)=1 along from NP6
                  IF(NPL(1,0,nl).EQ.2.AND.NPL(2,1,nl).EQ.NP6) THEN
                    NP8=NPL(3,1,nl)
                    GO TO 32
                  ENDIF
                ENDDO
 32             IF(NP8.GT.0) THEN
                  DO nolist=1,NBLIST(0)
                    nb=NBLIST(nolist)
                    NPNE(8,nb,ne)=NP8
                  ENDDO
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' NP5..8='',4I6)')
     '                NP5,NP6,NP7,NP8
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  NEELEM(ne,1)=ne
                  DO nj=1,NJT
                    nb=NBLIST(nj)
                    NBJ(nj,ne)=nb
                    DO nn=1,NNT(nb)
                      DO nk=1,NKT(nn,nb)
                        NKJE(nk,nn,nj,ne)=nk
                      ENDDO
                    ENDDO
                  ENDDO
                  DO nolist=1,NBLIST(0)
                    nb=NBLIST(nolist)
C KAT 23Feb01: now handled by NKJE
C                    DO nn=1,NNT(nb)
C                      DO nk=1,NKT(nn,nb)
C                        NKE(nk,nn,nb,ne)=nk
C                      ENDDO
C                    ENDDO
                    DO ns=1,NST(nb)+NAT(nb)
                      SE(ns,nb,ne)=1.0D0
                    ENDDO
                  ENDDO
                ELSE
                  ERROR=' >>Node for nn=8 not found'
                  GO TO 9999
                ENDIF
              ENDIF
            ENDDO
            NET(nr)=ne
            NEELEM(0,nr)=NET(nr)

          ENDIF

C         set up basis fn, #derivs and mapping arrays for geometry
          CALL GLOBALJ(NBJ,NEELEM,NKJ,NPLIST1,NPNE,NPNODE,
     '      ERROR,*9999)

C         define global line segments
          CALL LINSEG(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '      NLLINE,NNL,NPL,NPNE,NPNODE,NVJE,NVJL,XP,ERROR,*9999)
          DO nl=1,NLT
            CALL ARCSCA(IDO,0,0,0,NBJ,NEL(0,nl),nl,
     '        NPL(1,0,nl),NPNE,NVJL(1,1,nl),DL,1.0D-6,XP,ERROR,*9999)
          ENDDO

        ENDIF
        IF(FILIO) THEN
          IF(BINARY) THEN
            DO no_nrlist=1,NRLIST(0)
              IF(NJ_LOC(NJL_GEOM,0,nr).GT.0) THEN
                CALL_ELEM=.TRUE.
              ELSEIF(NJ_LOC(NJL_FIBR,0,nr).GT.0) THEN
                CALL_ELFB=.TRUE.
              ELSEIF(NJ_LOC(NJL_FIEL,0,nr).GT.0) THEN
                CALL_ELFD=.TRUE.
              ENDIF
            ENDDO
            CALL_ELEM=.TRUE.
          ELSE
            IF(TYPE(1:8).EQ.'GEOMETRY') THEN
              CALL_ELEM=.TRUE.
            ELSE IF(TYPE(1:5).EQ.'FIBRE') THEN
              CALL_ELFB=.TRUE.
            ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
              CALL_ELFD=.TRUE.
            ENDIF
          ENDIF
        ELSE
          CALL_ELEM=.TRUE.
        ENDIF

C KAT 4May99: Can only calculate NUNK if NENP is set up
C CS 14/9/98 new added NUNK
        CALL CALC_NUNK(IDO,NBJ,NENP,NKJE,NKJ,NPNE,NPNODE,
     '    NRLIST,NUNK,ERROR,*9999)

C       CS 6/7/99 not using for now - very slow, need a better way
C       CS 18/5/98 new
C        IF(JTYP2C.EQ.1) THEN
C          DO no_nrlist=1,NRLIST(0)
C            nr=NRLIST(no_nrlist)
C            CALL HANGING_NODE_DETECT(IBT,IDO,INP,NBJ,NEELEM,NELIST,
C     '        NENP,NKE,NPF,NPLIST1,
C     '        NPNE,NPNODE,nr,NVJE,NWP,SE,XA,XE,XP,ERROR,
C     '        *9999)
C          ENDDO ! no_nrlist
C        ENDIF ! JTYP2C
C KSB Adding elem group:
        IF(MAKE_GROUP) THEN
          STRING=GROUP_NAME
          i=0
          DO noelem=noelem_start+1,NEELEM(0,nr)
            i=i+1
            NELIST(i)=NEELEM(noelem,nr)
          ENDDO !noelem
          NELIST(0)=i
          CALL GRELEM_SUB(NELIST,STRING,.TRUE.,ERROR,*9999)          
        ENDIF !MAKE_GROUP
        
      ENDIF

      CALL EXITS('DEELEM')
      RETURN
 9999 CALL ERRORS('DEELEM',ERROR)
      CALL EXITS('DEELEM')
      RETURN 1
      END


