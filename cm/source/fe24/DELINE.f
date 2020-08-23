      SUBROUTINE DELINE(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,
     '  NLL,NLLINE,NNL,NPINTER,NPL,NPNE,NPNODE,NRE,
     '  NRLIST,NVJE,NVJL,DL,SE,XP,STRING,ERROR,*)

C#### Subroutine: DELINE
C###  Description:
C###    DELINE defines lines.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),NNL(0:4,12,NBFM),
     '  NPINTER(3,0:20),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRE(NEM),
     '  NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM)
      REAL*8 DL(3,NLM),SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IPFILE,IWK(6),JDER,KOUNT,N3CO,
     '  nb,ni,nl,NL1,NL2,NL3,nn,nointer,NP1,NP2,NP3,NTIW
      CHARACTER FILE*(MXCH),STATUS*3,TYPE*8
      LOGICAL ALL_REGIONS,CALCU,CBBREV,FILIO,FIRST_TIME,GENER,
     '  INTERPOLATED_NODES,MOUSE

      CALL ENTERS('DELINE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define line;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    The line scaling information is read from or written to a file
C###    FILENAME.ipline.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define line;c
C###  Description:
C###    This command calculates the lines to connect the nodes, and
C###    which elements they belong to.
C###  Parameter:      <derivatives>
C###    Will update nodal derivatives during scale factor calculation.
C###  Parameter:      <interpolated_nodes>
C###    Update line derivatives for global node derivatives and bicubic
C###    Note: this is designed for constructing 3D heart mesh from
C###    2D surfaces & need to interpolate nodes in midwall which use
C###    read in global line derivs rather than calculated arclengths
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'<derivatives>'
        OP_STRING(3)=BLANK(1:15)//'<interpolated_nodes>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define line;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]> mapping
C###  Description:
C###    Defines any non-standard line mappings, to be read from a
C###    file.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
     '    //' mapping'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DELINE',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 24-Jan-1990
        CALL PARSE_QUALIFIERS('CDLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(MOUSE) CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(FILIO) THEN
          CALL ASSERT(CALL_LINE,'>>no lines defined',ERROR,*9999)
          IF(CBBREV(CO,'MAPPING',1,noco+1,NTCO,N3CO)) THEN
            TYPE='MAPPING'
            CALL ASSERT(JTYP2B.EQ.1,'>>Non-standard mapping not set',
     '        ERROR,*9999)
            WRITE(OP_STRING,'('' >>WARNING: Line mappings may not be '
     '        //'maintained in all places.'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE
            TYPE='STANDARD'
          ENDIF

        ELSE IF(CALCU) THEN
          IF(CBBREV(CO,'DERIVATIVES',1,noco+1,NTCO,N3CO)) THEN
            JDER=1 !to update nodal derivs during scale factor calc
          ELSE
            JDER=0 !to not change nodal derivs during scale factor calc
          ENDIF
          IF(CBBREV(CO,'INTERPOLATED_NODES',1,noco+1,NTCO,N3CO)) THEN
            INTERPOLATED_NODES=.TRUE.
          ELSE
            INTERPOLATED_NODES=.FALSE.
          ENDIF
        ENDIF

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'line',
     '        STATUS,ERR,ERROR,*9999)
            CALL IPLINE(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,
     '        NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRLIST,NVJE,NVJL,
     '        DL,SE,XP,TYPE,ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
          IF(TYPE(1:7).EQ.'MAPPING') THEN
C           Update the DL and SE array for those lines that have been
C           mapped
            DO nl=1,NLT
              IF(NPL(4,0,nl).LT.0) THEN
                DL(1,nl)=-1.0D0*DL(2,ABS(NPL(4,0,nl)))
                DL(2,nl)=-1.0D0*DL(1,ABS(NPL(4,0,nl)))
                DL(3,nl)=DL(3,ABS(NPL(4,0,nl)))
              ENDIF
            ENDDO
            DO nb=1,NBFT
              CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,
     '          DL,SE,ERROR,*9999)
            ENDDO
          ENDIF

        ELSE IF(CALCU) THEN
          IF(INTERPOLATED_NODES) THEN
C ***       Update line derivatives for global node derivs & bicubic
C ***       Note: this is designed for constructing 3D heart mesh from
C ***       2D surfaces & need to interpolate nodes in midwall which use
C ***       read in global line derivs rather than calculated arclengths
            CALL LINSEG(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '        NLLINE,NNL,NPL,NPNE,NPNODE,NVJE,NVJL,XP,ERROR,*9999)
            nb=NBJ(1,NEELEM(1,NRLIST(1)))

C CPB 8/10/94 swapping over nbi=4/5
C            IF(NBI(nb).EQ.3.AND.NKT(0,nb).EQ.4) THEN
            IF(NBI(nb).EQ.3.AND.NKT(0,nb).EQ.5) THEN
              DO nointer=1,NPINTER(1,0) !loops over interpolated nodes
                WRITE(OP_STRING,'('' Interpolated node #'',I2)')
     '            nointer
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                NP1=NPINTER(1,nointer) !is interpolated node
                NP2=NPINTER(2,nointer) !is adjacent node
                NP3=NPINTER(3,nointer) !is adjacent node
                WRITE(OP_STRING,'('' np1='',I3,'' np2='',I3,'
     '            //''' np3='',I3)') NP1,NP2,NP3
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                DO ni=1,2 !loops over 2 Xi directions
                  DO nn=1,2 !loops over node position on line
                    IF(DOP) THEN
                      WRITE(OP_STRING,
     '                  '('' ni='',I1,'' nn='',I1)') ni,nn
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    NL1=0
                    NL2=0
                    NL3=0
                    KOUNT=0
                    DO nl=1,NLT !loops over global lines
                      IF(NPL(1,0,nl).eq.ni) THEN !line in Xi(ni) direction
                        IF(NPL(1+nn,1,nl).EQ.NP1) THEN
                          NL1=nl
                          KOUNT=KOUNT+1
                        ELSE IF(NPL(1+nn,1,nl).EQ.NP2) THEN
                          NL2=nl
                        ELSE IF(NPL(1+nn,1,nl).EQ.NP3) THEN
                          NL3=nl
                        ENDIF
                      ENDIF
                    ENDDO
                    IF(KOUNT.GT.1) THEN
                      WRITE(OP_STRING,
     '                  '('' >>Warning! Examine '
     '                  //'line scaling factors at node '',I4)') NP1
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                    WRITE(OP_STRING,'('' nl1='',I3,'' nl2='',I3,'
     '                //''' nl3='',I3)') NL1,NL2,NL3
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    IF(nl1.NE.0) THEN
                      IF(nl2.NE.0.and.nl3.NE.0) THEN
                        DL(1,NL1)=0.5D0*(DL(1,NL2)+DL(1,NL3))
                        DL(2,NL1)=0.5D0*(DL(2,NL2)+DL(2,NL3))
                      ELSE IF(nl2.NE.0) THEN
                        DL(1,NL1)=DL(1,NL2)
                        DL(2,NL1)=DL(2,NL2)
                      ELSE IF(nl3.NE.0) THEN
                        DL(1,NL1)=DL(1,NL3)
                        DL(2,NL1)=DL(2,NL3)
                      ELSE
                        DL(1,NL1)=0.0D0
                        DL(2,NL1)=0.0D0
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

          CALL LINCAL(IBT,IDO,INP,JDER,NBJ,NEELEM,NEL,NENP,
     '      NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '      DL,SE,XP,ERROR,*9999)

        ELSE IF(MOUSE) THEN
        ENDIF
        CALL_LINE=.TRUE.
      ENDIF

      CALL EXITS('DELINE')
      RETURN
 9999 CALL ERRORS('DELINE',ERROR)
      CALL EXITS('DELINE')
      RETURN 1
      END


