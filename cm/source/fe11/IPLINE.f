      SUBROUTINE IPLINE(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,
     '  NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRLIST,NVJE,NVJL,
     '  DL,SE,XP,TYPE,ERROR,*)

C#### Subroutine: IPLINE
C###  Description:
C###    IPLINE inputs scaling factors for NBI=2 and 3 cases

C**** NBI(nb) is type of scaling factor input for basis nb.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NENP(NPM,0:NEPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),
     '  NNL(0:4,12,NBFM),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM)
      REAL*8 DL(3,NLM),SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER IBEG0,IBEG1,IBEG2,IBEG3,ICHAR,IEND0,IEND1,IEND2,IEND3,
     '  INFO,nb,NB_NJ,ne,NETOT,ni,nj,njj,njj1,njj2,nk,NL,NL_TOT,nn,
     '  noelem,noline,nonrlist,NOQUES,NP1,NP2,nr,ns
      CHARACTER CHAR0*2,CHAR1*4,CHAR2*2,CHAR3*4,CHAR4*4
      LOGICAL CONTINUE,CUBIC,FILEIP,UPDATE_DL,FOUND

      CALL ENTERS('IPLINE',*9999)

      ICHAR=999
      NOQUES=0
C CPB 18/9/93 Adding non-standard line mappings
      IF(TYPE(1:7).EQ.'MAPPING') THEN !non-standard line mapping
        CONTINUE=.TRUE.
        IF(IOTYPE.NE.3) THEN
          DO nl=1,NLM
            NPL(4,0,nl)=0
          ENDDO
        ENDIF
        nl=0
        DO WHILE(CONTINUE)
          FORMAT='(/$,'' Enter the line number to be mapped'
     '      //' [EXIT]: '',I4)'
          IF(IOTYPE.EQ.3) THEN
            CONTINUE=.FALSE.
            FOUND=.FALSE.
            IDATA(1)=0
C           Try and find a line that has been mapped to write out
            DO WHILE(nl.LT.NLT.AND.(.NOT.FOUND))
              nl=NL+1
              IF(NPL(4,0,nl).LT.0) THEN
                CONTINUE=.TRUE.
                FOUND=.TRUE.
                IDATA(1)=NL
              ENDIF
            ENDDO
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NLT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c         CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,1,NLT,
c    '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(IDATA(1).EQ.0) THEN
              CONTINUE=.FALSE.
            ELSE IF(IDATA(1).NE.0) THEN
              CONTINUE=.TRUE.
              nl=IDATA(1)
            ENDIF
          ENDIF
          IF(CONTINUE) THEN
            WRITE(CHAR1,'(I4)') nl
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            FORMAT='($,'' Enter the line number to map line '//
     '        CHAR1(IBEG1:IEND1)//' from ['//CHAR1(IBEG1:IEND1)
     '        //']: '',I4)'
            IDEFLT(1)=NL
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=ABS(NPL(4,0,nl))
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        NLT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c           CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,NLT,
c    '        INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              NPL(4,0,IDATA(1))=NL
              NPL(4,0,nl)=-1*IDATA(1)
            ENDIF
          ENDIF
        ENDDO
      ELSE IF(TYPE(1:8).EQ.'STANDARD') THEN
        FILEIP=.FALSE.
        NOQUES=0
        DO nb=1,NBFT
          CONTINUE=.FALSE.
          DO nonrlist=1,NRLIST(0)
            nr=NRLIST(nonrlist)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(NBI(nb).EQ.2.OR.NBI(nb).EQ.3) THEN
C               Read or write element/global scale factors SE/DL
                DO njj1=1,3 !geometry, fibres, field
                  DO njj2=1,NJ_LOC(njj1,0,nr)
                    nj=NJ_LOC(njj1,njj2,nr)
                    IF(nb.EQ.NBJ(nj,ne)) CONTINUE=.TRUE.
                  ENDDO !njj2
                ENDDO !njj1
              ELSE
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
C                 to check NB is geometric var basis function
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  IF(nb.EQ.NBJ(nj,ne)) CONTINUE=.TRUE.
                ENDDO !nj
              ENDIF
            ENDDO !ne
          ENDDO !nr
          IF(CONTINUE) THEN
            IF(NBI(nb).EQ.1) THEN
C             Unit scale factors
              IF(NNT(nb).GT.0) THEN
                DO nonrlist=1,NRLIST(0)
                  nr=NRLIST(nonrlist)
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    DO ns=1,NST(nb)+NAT(nb)
                      SE(ns,nb,ne)=1.0D0
                    ENDDO !ns
                  ENDDO !ne
                ENDDO !nr
C10/3/89    ELSE IF(IBT(1,1,nb).EQ.5) THEN
C             CALL SPLINE(IBT(1,1,nb),nb,NEELEM,VE,ERROR,*9999)
              ELSE
                DO nonrlist=1,NRLIST(0)
                  nr=NRLIST(nonrlist)
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    DO ns=1,NST(nb)+NAT(nb)
                      SE(ns,nb,ne)=1.0D0
                    ENDDO !ns
                  ENDDO !ne
                ENDDO !nr
              ENDIF

            ELSE IF(NBI(nb).EQ.2) THEN
C             Read or write element scale factors SE
              IF(IOTYPE.GT.0) THEN
                WRITE(CHAR0,'(I2)') nb
                CALL STRING_TRIM(CHAR0,IBEG0,IEND0)
                DO nonrlist=1,NRLIST(0)
                  nr=NRLIST(nonrlist)
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    WRITE(CHAR1,'(I4)') ne
                    CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
                    ns=0
                    DO nn=1,NNT(nb)
                      WRITE(CHAR2,'(I2)') nn
                      CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
                      DO nk=1,NKT(nn,nb)
                        ns=ns+1
                        WRITE(CHAR3,'(I2)') nk
                        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
                        FORMAT='($,'' Basis '//CHAR0(IBEG0:IEND0)
     '                    //'; element '//CHAR1(IBEG1:IEND1)
     '                    //'; vertex '//CHAR2(IBEG2:IEND2)
     '                    //': scale factor for nk='//
     '                    CHAR3(IBEG3:IEND3)//' is '',D25.16)'
                        RDEFLT(1)=SE(ns,nb,ne)
                        IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
                        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '                    FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,
     '                    CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
     '                    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
                        IF(IOTYPE.NE.3) SE(ns,nb,ne)=RDATA(1)
                      ENDDO !nk
C                      nss=ns
C          DO nk=1,NKT(nn,nb)
C            nss=nss+1
C            RDEFLT(nk)=SE(nss,nb,ne)
C            IF(IOTYPE.EQ.3) RDATA(nk)=RDEFLT(nk)
C          ENDDO
C          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,
C     '                  FILEIP,FORMAT,NKT(nn,nb),ADATA,ADEFLT,CDATA,
C     '                  CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
C     '                  RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C                      DO nk=1,NKT(nn,nb)
C                        ns=ns+1 !needs to be incremented for all IOTYPE
C                        IF(IOTYPE.NE.3) SE(ns,nb,ne)=RDATA(nk)
C                      ENDDO
                    ENDDO !nn
                  ENDDO !noelem
                ENDDO !nonrlist (nr)
              ENDIF

C             Calculate DL from SE (restricted cases handled at present)
              DO nl=1,NLT
C CPB 10/12/92 Find out if DL needs to be changed
                UPDATE_DL=.FALSE.
                DO nj=1,NJT
                  NETOT=NPL(4,2,nl)
                  DO noelem=1,NETOT
                    ne=NPL(4+noelem,2,nl)
                    NB_NJ=NBJ(nj,ne)
C                   If any of the elements around the current line
C                   uses the current basis as their geometric basis
C                   then update DL
                    IF(NB_NJ.EQ.NB) UPDATE_DL=.TRUE.
                  ENDDO
                ENDDO
                IF(UPDATE_DL) THEN
                  NI =NPL(1,0,nl)
                  NP1=NPL(2,1,nl)
                  NP2=NPL(3,1,nl)
                  NE =NEL(1,nl)
                  DO nn=1,NNT(nb)
C!!! ni to nk conversion not done!  BAD luck if ni=3!
                    IF(NPNE(nn,nb,ne).EQ.NP1) THEN
                      DL(1,nl)=SE(1+NI+(NN-1)*NKT(0,nb),nb,ne)
                    ELSE IF(NPNE(nn,nb,ne).EQ.NP2) THEN
                      DL(2,nl)=SE(1+NI+(NN-1)*NKT(0,nb),nb,ne)
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO

            ELSE IF(NBI(nb).EQ.3) THEN
C             Read or write global scale factors DL(1,nl) & DL(2,nl)
C             First do line topology
              CALL LINSEG(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '          NLLINE,NNL,NPL,NPNE,NPNODE,NVJE,NVJL,XP,ERROR,*9999)
              IDEFLT(1)=NLT
              WRITE(CHAR1,'(I4)') NLT
              IF(IOTYPE.EQ.3) THEN
                NL_TOT=NLT
                IDATA(1)=NL_TOT
              ENDIF
              FORMAT='($,'' Number of lines ['//CHAR1(1:4)//']: '',I4)'
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          0,NLM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '          *9999)
              IF(IOTYPE.NE.3) NL_TOT=IDATA(1)

              DO noline=1,NL_TOT
                IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.3) THEN
C                 prompted input or write file
C!!!!!May need updating?
                  nl=noline
                  CUBIC=.FALSE.
                  DO nj=1,NJT
                    IF(NPL(1,nj,nl).EQ.4) CUBIC=.TRUE.
                  ENDDO
                  IF(CUBIC) THEN
                    CONTINUE=.TRUE.
                    IDEFLT(1)=NPL(2,1,nl)
                    IDEFLT(2)=NPL(3,1,nl)
                    RDEFLT(1)=DL(1,nl)
                    RDEFLT(2)=DL(2,nl)
                    WRITE(CHAR2,'(I1)') NPL(1,0,nl)
                    WRITE(CHAR3,'(I4)') IDEFLT(1)
                    WRITE(CHAR4,'(I4)') IDEFLT(2)
                  ENDIF
                ELSE
                  CONTINUE=.TRUE.
                  CHAR2='?'
                  CHAR3='   ?'
                  CHAR4='   ?'
                ENDIF
                IF(CONTINUE) THEN
                  WRITE(CHAR1,'(I4)') noline
                  FORMAT='($,'' '//CHAR1(1:4)//') Arc in Xi('//
     '              CHAR2(1:1)//') direction has nodes ['//
     '              CHAR3(1:4)//','//CHAR4(1:4)//']: '',2I5)'
                  IF(IOTYPE.EQ.3) THEN
                    IDATA(1)=IDEFLT(1)
                    IDATA(2)=IDEFLT(2)
                  ENDIF
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,1,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '              RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) THEN
                    NP1=IDATA(1)
                    NP2=IDATA(2)
                  ENDIF

                  FORMAT='($,7X,''Scale factors are: '',2(D11.4,1X))'
                  IF(IOTYPE.EQ.3) THEN
                    RDATA(1)=RDEFLT(1)
                    RDATA(2)=RDEFLT(2)
                  ENDIF
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '              RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) THEN
                    DO nl=1,NLT
                      IF(NPL(2,1,nl).EQ.NP1.AND.NPL(3,1,nl).EQ.NP2) THEN
                        DL(1,nl)=RDATA(1)
                        DL(2,nl)=RDATA(2)
                      ENDIF
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO
C             Calculate SE from DL
              CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,
     '          DL,SE,ERROR,*9999)

c cpb 8/10/94 Changing around nbi(nb) = 4 and 5

            ELSE IF(NBI(nb).EQ.4) THEN
C             Calculate DL derivatives from angle change
              DO nl=1,NLT
C CPB 10/12/92 Find out if DL needs to be changed
                UPDATE_DL=.FALSE.
                DO nj=1,NJT
                  NETOT=NPL(4,2,nl)
                  DO noelem=1,NETOT
                    ne=NPL(4+noelem,2,nl)
                    NB_NJ=NBJ(nj,ne)
C                   If any of the elements around the current line
C                   uses the current basis as their geometric basis
C                   then update DL
                    IF(NB_NJ.EQ.NB) UPDATE_DL=.TRUE.
                  ENDDO
                ENDDO
                IF(UPDATE_DL) THEN
                  CALL ANGSCA(NPL(1,0,nl),DL(1,nl),XP,ERROR,*9999)
                ENDIF
              ENDDO

C             Calculate SE from DL
              CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,
     '          DL,SE,ERROR,*9999)

            ELSE IF(NBI(nb).GE.5.AND.NBI(nb).LE.7) THEN
C             Calculate DL from average arclength
              DO nl=1,NLT
C CPB 10/12/92 Find out if DL needs to be changed
                UPDATE_DL=.FALSE.
                DO nj=1,NJT
                  NETOT=NPL(4,2,nl)
                  DO noelem=1,NETOT
                    ne=NPL(4+noelem,2,nl)
                    NB_NJ=NBJ(nj,ne)
C                   If any of the elements around the current line
C                   uses the current basis as their geometric basis
C                   then update DL
                    IF(NB_NJ.EQ.NB) UPDATE_DL=.TRUE.
                  ENDDO
                ENDDO
                IF(UPDATE_DL) THEN
                  CALL ARCSCA(IDO,0,0,0,NBJ,NEL(0,nl),NL,
     '              NPL(1,0,nl),NPNE,NVJL,DL,1.0D-6,XP,ERROR,*9999)
                ENDIF
              ENDDO

C             Calculate SE from DL
              CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,
     '          DL,SE,ERROR,*9999)

            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('IPLINE')
      RETURN
 9999 CALL ERRORS('IPLINE',ERROR)
      CALL EXITS('IPLINE')
      RETURN 1
      END


