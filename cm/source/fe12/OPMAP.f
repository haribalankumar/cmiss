      SUBROUTINE OPMAP(IBT,IDO,INP,NBH,NBJ,
     '  NCLIST,NEELEM,NENP,NHE,NHP,NKB,NKHE,NKH,
     '  NKJ,NKJE,NLL,NNB,NNF,NNL,NONY,NPF,NPL,NPNE,NPNODE,NPNY,nr,
     '  NRCLIST,NVHE,NVHP,NVJE,NVJP,NWP,nx,
     '  NXI,NYNE,NYNO,NYNP,NYNR,NYNY,CONY,CYNO,CYNY,SE,SP,XA,XE,XP,
     '  NXTYPE,TYPE,FIX,ERROR,*)

C#### Subroutine: OPMAP
C###  Description:
C###    Outputs mapping arrays for region nr.

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NCLIST(0:4),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NHE(NEM),NHP(NPM),NKB(2,2,2,NNM,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NKH(NHM,NPM,NCM),NKJ(NJM,NPM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NNB(4,4,4,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NONY(0:NOYM,NYM,NRCM,0:NRM),
     '  NPF(9,NFM),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),nr,NRCLIST(0:3),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJP(NJM,NPM),NWP(NPM,2),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NYNY(0:NYYM,NYM,NRM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  CYNY(0:NYYM,NYM,NRM),SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),NXTYPE*(*),TYPE*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER i,na,nb,nc,ne,nh,nhj,nhx,nj,njj1,njj2,nk,NK_TOT,nn,no,
     '  noelem,nonc,nonode,nonrc,no_nynr,noy,np,nrc,nrcc,nv,ny,ny_type,
     '  nyo
      CHARACTER CHAR1,CHAR2,CFROMI,NC_STR(4)*25,NRC_STR(0:2)*25,
     '  NY_STR(2)*4
      DATA NC_STR  /'   (LHS stiffness matrix)',
     '              '   (RHS matrix terms)    ',
     '              '   (Damping, etc. matrix)',
     '              '   (Mass matrix terms)   '/
      DATA NRC_STR /'   (global variables)    ',
     '              '   (local equations|rows)',
     '              '   (local variables|cols)'/
      DATA NY_STR /'NYNP','NYNE'/

      CALL ENTERS('OPMAP',*9999)

      WRITE(OP_STRING,'(/'' Region '',I1,'':'')') nr
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(TYPE(1:4).EQ.'LINE') THEN
        WRITE(OP_STRING,'(/'' Line mapping arrays:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(TYPE(1:8).EQ.'MATERIAL') THEN
        WRITE(OP_STRING,'(/'' Material mapping arrays:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(TYPE(1:4).EQ.'MESH') THEN
        WRITE(OP_STRING,'(/'' Mesh mapping arrays:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nonc=1,NCLIST(0)
          nc=NCLIST(nonc)
          WRITE(OP_STRING,'(/'' nc='',I1,A25)') nc,NC_STR(nc)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nonrc=1,NRCLIST(0)
            nrc=NRCLIST(nonrc)
            WRITE(OP_STRING,'(/'' nrc='',I1,A)')
     '        nrc,NRC_STR(nrc)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(/'' NYNR(0,'',I1,'','',I1,'',nr)='',I5)')
     '        nrc,nc,NYNR(0,nrc,nc,nr)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C cpb 1/2/00 Adding long write
C            WRITE(OP_STRING,'('' NYNR(1..,'',I1,'','',I1,'',nr): '','
C     '        //'10(I5,X),/:(19X,10(I5,X)))') nrc,nc,
C     '        (NYNR(ny,nrc,nc,nr),ny=1,NYNR(0,nrc,nc,nr))
C            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            CHAR1=CFROMI(nrc,'(I1)')
            CHAR2=CFROMI(nc,'(I1)')
            CALL WRITE_LONG(INTTYPE,1,1,IOFI,NYNR(0,nrc,nc,nr),10,10,
     '        NYNR(1,nrc,nc,nr),%VAL(0),'('' NYNR(1..'//CHAR1//
     '        ','//CHAR2//',nr): '',10(I5,X))','(18X,10(I5,X))',
     '        ERROR,*9999)
            WRITE(OP_STRING,'(/)')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(NXTYPE(1:5).EQ.'SOLVE') THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nhx=1,NHP(np)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,nc)
                    IF(nrc.EQ.1)THEN
                      NK_TOT=MAX(NKH(nh,np,1),NKH(nh,np,2))
                    ELSE
                      NK_TOT=NKH(nh,np,nc)
                    ENDIF
                    DO nk=1,NK_TOT
                      ny=NYNP(nk,nv,nh,np,nrc,nc,nr)
                      WRITE(OP_STRING,'('' nc='',I1,'',nrc='',I1,'
     '                  //''',np='',I4,'',nh='',I1,'',nv='',I2,'
     '                  //''',nk='',I1,'
     '                  //''': NYNP='',I5,'
     '                  //''' NPNY(0..,ny,nrc): '',I1,1X,I1,1X,'
     '                  //'I2,1X,I1,1X,I4,1X,I1,1X,I1)')
     '                  nc,nrc,np,nh,nv,nk,ny,(NPNY(i,ny,nrc),i=0,6)
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nh
              ENDDO !nonode (np)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                DO nhx=1,NHE(ne)
                  nh=NH_LOC(nhx,nx)
                  IF(nrc.EQ.1) THEN
C!!! Use the L.H.S. (nc=1) basis to determine the # of equations
                    nb=NBH(nh,1,ne)
                  ELSE
                    nb=NBH(nh,nc,ne)
                  ENDIF
                  DO na=1,NAT(nb)
                    ny=NYNE(na,nh,nrc,nc,ne)
                    WRITE(OP_STRING,'('' nc='',I1,'',nrc='',I1,'
     '                //''',ne='',I4,'',nh='',I1,'',na='',I2,5X,'
     '                //''': NYNE='',I5,'
     '                //''' NPNY(0..,ny,nrc): '',I1,1X,I1,1X,'
     '                //'I2,1X,I1,1X,I4,1X,I1)') nc,nrc,ne,nh,na,ny,
     '                (NPNY(i,ny,nrc),i=0,5)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDDO !na
                ENDDO !nh
              ENDDO !noelem (ne)

            ELSE IF(NXTYPE(1:3).EQ.'FIT') THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO njj1=1,NUM_FIT(0)
                  DO nhj=1,NUM_FIT(njj1)
                    nhx=NLH_FIT(nhj,3,njj1)
                    nh=NH_LOC(nhx,nx)
                    DO nv=1,NVHP(nh,np,1)
                      DO nk=1,NKH(nh,np,1)
                        ny=NYNP(nk,nv,nh,np,nrc,nc,nr)
                        WRITE(OP_STRING,'('' nc='',I1,'',nrc='',I1,'
     '                    //''',np='',I5,'',nh='',I1,'',nv='',I1,'
     '                    //''',nk='',I1,'
     '                    //''': NYNP='',I5,'
     '                    //''' NPNY(0..,ny,nrc): '',I1,1X,I1,1X,'
     '                    //'I1,1X,I1,1X,I5,1X,I1,1X,I1)') nc,nrc,np,nh,
     '                    nv,nk,ny,(NPNY(i,ny,nrc),i=0,6)
                        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                      ENDDO !nk
                    ENDDO !nv
                  ENDDO !nhj
                ENDDO !njj1
              ENDDO !nonode (np)
            ELSE IF(NXTYPE(1:8).EQ.'OPTIMISE') THEN
              IF(KTYP8.EQ.6) THEN !data fitting by optimisation
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  DO njj1=1,NUM_FIT(0)
                    DO nhj=1,NUM_FIT(njj1)
                      nj=NLH_FIT(nhj,1,njj1)
                      nhx=NLH_FIT(nhj,3,njj1)
                      nh=NH_LOC(nhx,nx)
                      DO nv=1,NVJP(nj,np)
                        DO nk=1,NKJ(nj,np)
                          ny=NYNP(nk,nv,nh,np,nrc,nc,nr)
                          WRITE(OP_STRING,'('' nc='',I1,'',nrc='',I1,'
     '                      //''',np='',I4,'',nh='',I1,'',nv='',I2,'
     '                      //''',nk='',I1,'': NYNP='',I5,'
     '                      //''' NPNY(0..,ny,nrc): '',I1,1X,I1,1X,'
     '                      //'I2,1X,I1,1X,I4,1X,I1,1X,I1)') nc,nrc,np,
     '                      nh,nv,nk,ny,(NPNY(i,ny,nrc),i=0,6)
                          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                        ENDDO !nk
                      ENDDO !nv
                    ENDDO !nhj
                  ENDDO !njj1
                ENDDO !nonode (np)
              ELSE
C*** Other optimisation cases
              ENDIF !ktyp8
            ENDIF !nxtype
          ENDDO !nrc
        ENDDO !nonc (nc)
      ELSE IF(TYPE(1:8).EQ.'SOLUTION') THEN
        WRITE(OP_STRING,'(/'' Solution mapping arrays:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Mesh <--> Solution arrays'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nonc=1,NCLIST(0)
          nc=NCLIST(nonc)
          WRITE(OP_STRING,'(/'' nc='',I1,A25)') nc,NC_STR(nc)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nonrc=1,NRCLIST(0)
            nrc=NRCLIST(nonrc)
            IF(nrc.EQ.2) THEN
              nrcc=0  !for accessing NYNR etc arrays with nrc
            ELSE
              nrcc=nrc
            ENDIF
            WRITE(OP_STRING,'(/'' nrc='',I1,A)')
     '        nrc,NRC_STR(nrcc)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C cpb 12/3/95 Changing the loops for the solution arrays. Must loop
C over NYNR for the particular region being listed and not loop
C over node, variables, derivatives etc. as this does not pick up
C multiple ny's at an interface node.

            IF(NXTYPE(1:5).EQ.'SOLVE') THEN

C cpb 12/3/95 old
C              DO nonode=1,NPNODE(0,nr)
C                np=NPNODE(nonode,nr)
C                DO nh=1,NHP(np)
C                  DO nv=1,NVHP(nh,np,nc)
C                    IF(nrc.EQ.1) THEN
C                      NK_TOT=MAX(NKH(nh,np,1),NKH(nh,np,2))
C                    ELSE
C                      NK_TOT=NKH(nh,np,nc)
C                    ENDIF
C                    DO nk=1,NK_TOT
C                      ny=NYNP(nk,nv,nh,np,nrcc,nc,nr)
C cpb 12/3/95 new
              DO no_nynr=1,NYNR(0,nrcc,nc,nr)
                ny=NYNR(no_nynr,nrcc,nc,nr)
                ny_type=NPNY(0,ny,nrcc)
                IF(NONY(0,ny,nrc,nr).EQ.0) THEN
                  WRITE(OP_STRING,'(1X,A,''='',I5,'
     '              //''' NONY(0,ny)  ='',I6)')
     '              NY_STR(ny_type),ny,NONY(0,ny,nrc,nr)
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ELSE IF(NONY(0,ny,nrc,nr).GT.0) THEN
                  WRITE(OP_STRING,'(1X,A,''='',I5,'
     '              //''' NONY(0,ny)  ='',I6,'
     '              //''' CONY(0,ny)='',F10.4)')
     '              NY_STR(ny_type),ny,NONY(0,ny,nrc,nr),
     '              CONY(0,ny,nrc,nr)
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(12X,''NONY(1..,ny)='','
     '              //'9(1X,I5))') (NONY(noy,ny,nrc,nr),noy=1,
     '              NONY(0,ny,nrc,nr))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(12X,''CONY(1..,ny)='','
     '              //'9(1X,F10.4))') (CONY(noy,ny,nrc,nr),noy=1,
     '              NONY(0,ny,nrc,nr))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO !no_nynr

C 24/2/97 LC archived section : cpb 12/3/95 old

            ELSE IF(NXTYPE(1:3).EQ.'FIT') THEN
              DO njj1=1,NUM_FIT(0)
                WRITE(OP_STRING,'(/'' Fit variable : '',I1)') njj1
                CALL GLOBALF(IBT,IDO,INP,NBH,NBJ,
     '            NENP,njj1,NKB,NKH,NKHE,NKJE,NLL,NNB,NNF,NNL,NONY,NPF,
     '            NPL,NPNE,NPNODE,NPNY,nr,NVHE,NVHP,NVJE,NWP,
     '            nx,NXI,NYNE,NYNO,NYNP,
     '            NYNR,NYNY,CONY,CYNO,CYNY,SE,SP,XA,XE,XP,FIX,ERROR,
     '            *9999)
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  DO nhj=1,NUM_FIT(njj1)
                    nhx=NLH_FIT(nhj,3,njj1)
                    nh=NH_LOC(nhx,nx)
                    WRITE(OP_STRING,'('' Variable '',I1,'', nh='',I1)')
     '                nhj,nh
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    DO nv=1,NVHP(nh,np,1)
                      DO nk=1,NKH(nh,np,1)
                        ny=NYNP(nk,nv,nh,np,nrcc,1,nr)
                        IF(NONY(0,ny,nrc,nr).EQ.0) THEN
                          WRITE(OP_STRING,'('' NYNP='',I5,'
     '                      //''' NONY(0,ny)  ='',I6)')
     '                      ny,NONY(0,ny,nrc,nr)
                          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                        ELSE IF(NONY(0,ny,nrc,nr).GT.0) THEN
                          WRITE(OP_STRING,'('' NYNP='',I5,'
     '                      //''' NONY(0,ny)  ='',I6,'
     '                      //''' CONY(0,ny)='',F10.4)')
     '                      ny,NONY(0,ny,nrc,nr),CONY(0,ny,nrc,nr)
                          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                          WRITE(OP_STRING,'(12X,''NONY(1..,ny)='','
     '                      //'9(1X,I5))') (NONY(noy,ny,nrc,nr),noy=1,
     '                      NONY(0,ny,nrc,nr))
                          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                          WRITE(OP_STRING,'(12X,''CONY(1..,ny)='','
     '                      //'9(1X,F10.4))') (CONY(noy,ny,nrc,nr),
     '                      noy=1,NONY(0,ny,nrc,nr))
                          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDDO !nk
                    ENDDO !nv
                  ENDDO !nonode (np)
                ENDDO !nhj
              ENDDO !njj1
            ELSE IF(NXTYPE(1:8).EQ.'OPTIMISE') THEN
              IF(KTYP8.EQ.6) THEN !data fitting by optimisation
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  DO njj1=1,NUM_FIT(0)
                    DO nhj=1,NUM_FIT(njj1)
                      nj=NLH_FIT(nhj,1,njj1)
                      nhx=NLH_FIT(nhj,3,njj1)
                      nh=NH_LOC(nhx,nx)
                      DO nv=1,NVJP(nj,np)
                        DO nk=1,NKJ(nj,np)
                          ny=NYNP(nk,nv,nh,np,nrcc,nc,nr)
                          IF(NONY(0,ny,nrc,nr).EQ.0) THEN
                            WRITE(OP_STRING,'('' NYNP='',I5,'
     '                        //''' NONY(0,ny)  ='',I6)')
     '                        ny,NONY(0,ny,nrc,nr)
                            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                          ELSE IF(NONY(0,ny,nrc,nr).GT.0) THEN
                            WRITE(OP_STRING,'('' NYNP='',I5,'
     '                        //''' NONY(0,ny)  ='',I6,'
     '                        //''' CONY(0,ny)='',F10.4)')
     '                        ny,NONY(0,ny,nrc,nr),CONY(0,ny,nrc,nr)
                            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                            WRITE(OP_STRING,'(12X,''NONY(1..,ny)='','
     '                        //'9(1X,I5))') (NONY(noy,ny,nrc,nr),noy=1,
     '                        NONY(0,ny,nrc,nr))
                            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                            WRITE(OP_STRING,'(12X,''CONY(1..,ny)='','
     '                        //'9(1X,F10.4))')
     '                        (CONY(noy,ny,nrc,nr),noy=1,
     '                        NONY(0,ny,nrc,nr))
                            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                          ENDIF
                        ENDDO !nk
                      ENDDO !nv
                    ENDDO !nhj
                  ENDDO !njj1
                ENDDO !nonode (np)
              ELSE
C*** Other optimisation cases
              ENDIF !ktyp8
            ENDIF !nxtype
          ENDDO !nonrc (nrc)
        ENDDO !nonc (nc)
        WRITE(OP_STRING,'(/'' Solution <--> Mesh arrays'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(NXTYPE(1:3).EQ.'FIT') THEN
          DO njj1=1,NUM_FIT(0)
            WRITE(OP_STRING,'(/'' Fit variable : '',I1)') njj1
            CALL GLOBALF(IBT,IDO,INP,NBH,NBJ,
     '        NENP,njj1,NKB,NKH,NKHE,NKJE,NLL,NNB,NNF,NNL,NONY,NPF,
     '        NPL,NPNE,NPNODE,NPNY,nr,NVHE,NVHP,NVJE,
     '        NWP,nx,NXI,NYNE,NYNO,NYNP,
     '        NYNR,NYNY,CONY,CYNO,CYNY,SE,SP,XA,XE,XP,FIX,ERROR,*9999)
            DO nonrc=1,NRCLIST(0)
              nrc=NRCLIST(nonrc)
              IF(nrc.EQ.2) THEN
                nrcc=0 !for accessing NYNR etc.
              ELSE
                nrcc=nrc
              ENDIF
              WRITE(OP_STRING,'(/'' nrc='',I1,'', NOT('',I1,'',1,'
     '          //''',I1,'','',I1,'')='',I5,A)')
     '          nrc,nrc,nr,nx,NOT(nrc,1,nr,nx),NRC_STR(nrcc)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO no=1,NOT(nrc,1,nr,nx)
                WRITE(OP_STRING,'('' no='',I5,'' NYNO(0,no)  ='',I6,'
     '            //''' CYNO(0,no)='',F10.4)')
     '            no,NYNO(0,no,nrc,nr),CYNO(0,no,nrc,nr)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                IF(NYNO(0,no,nrc,nr).GT.0) THEN
                  WRITE(OP_STRING,'(10X,''NYNO(1..,no)='',9(1X,I5))')
     '              (NYNO(nyo,no,nrc,nr),nyo=1,NYNO(0,no,nrc,nr))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(10X,''CYNO(1..,no)='','
     '              //'9(1X,F10.4))') (CYNO(nyo,no,nrc,nr),
     '              nyo=1,NYNO(0,no,nrc,nr))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO !no
            ENDDO !nonrc (nrc)
          ENDDO !njj1
        ELSE
          DO nonrc=1,NRCLIST(0)
            nrc=NRCLIST(nonrc)
            IF(nrc.EQ.2) THEN
              nrcc=0 !for accessing NYNR etc.
            ELSE
              nrcc=nrc
            ENDIF
            WRITE(OP_STRING,'(/'' nrc='',I1,'', NOT('',I1,'',1,'',I1,'
     '        //''','',I1,'')='',I5,A)')
     '        nrc,nrc,nr,nx,NOT(nrc,1,nr,nx),NRC_STR(nrcc)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO no=1,NOT(nrc,1,nr,nx)
              WRITE(OP_STRING,'('' no='',I5,'' NYNO(0,no)  ='',I6,'
     '          //''' CYNO(0,no)='',F10.4)')
     '          no,NYNO(0,no,nrc,nr),CYNO(0,no,nrc,nr)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              IF(NYNO(0,no,nrc,nr).GT.0) THEN
                WRITE(OP_STRING,'(10X,''NYNO(1..,no)='',6(4X,I7),'
     '            //'/(23X,6(4X,I7)))')
     '            (NYNO(nyo,no,nrc,nr),nyo=1,NYNO(0,no,nrc,nr))
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(10X,''CYNO(1..,no)='',6(1X,F10.4),'
     '            //'/(23X,6(1X,F10.4)))')
     '            (CYNO(nyo,no,nrc,nr),nyo=1,NYNO(0,no,nrc,nr))
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !no
          ENDDO !nonrc (nrc)
        ENDIF
      ELSE IF(TYPE(1:8).EQ.'VERSIONS') THEN
        WRITE(OP_STRING,'(/'' Version mapping arrays:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nonc=1,NCLIST(0)
          nc=NCLIST(nonc)
          IF(nc.EQ.0) THEN
            WRITE(OP_STRING,'(/'' Number of versions for '
     '        //'geometric variables at nodes:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              WRITE(OP_STRING,'('' NVJP(nj=1..,np='',I4,''):'',12I2)')
     '          np,(NVJP(NJ_LOC(NJL_GEOM,njj2,nr),np),
     '          njj2=1,NJ_LOC(NJL_GEOM,0,nr))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !np
            WRITE(OP_STRING,'(/'' Version numbers for geometric '
     '        //'variables in elements:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nb=1,NBFT
                IF(NIT(nb).EQ.NJT.AND.NNT(nb).GT.0) THEN
                  DO njj1=1,3 !geom/fibres/field
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      WRITE(OP_STRING,'('' NVJE(nn=1..,nb='',I2,'
     '                  //''',nj='',I2,'',ne='',I4,''):'',32I2)')
     '                  nb,nj,ne,(NVJE(nn,nb,nj,ne),nn=1,NNT(nb))
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    ENDDO !njj2
                  ENDDO !njj1
                ENDIF
              ENDDO !nb
            ENDDO !ne

          ELSE IF(CALL_EQUA) THEN
            WRITE(OP_STRING,'(/'' nc='',I1,A25)') nc,NC_STR(nc)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(/'' Number of versions for '
     '        //'dependent variables at nodes:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              WRITE(OP_STRING,'('' NVHP(nh=1..,np='',I4,'',nc='',I1,'
     '          //'''):'',12I2)') np,nc,
     '         (NVHP(NH_LOC(nhx,nx),np,nc),nhx=1,NHP(np))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO  !np
            WRITE(OP_STRING,'(/'' Version numbers for dependent '
     '        //'variables in elements:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nb=1,NBFT
                IF(NIT(nb).EQ.NJT.AND.NNT(nb).GT.0) THEN
                  DO nhx=1,NHE(NE)
                    nh=NH_LOC(nhx,nx)
                    WRITE(OP_STRING,'('' NVHE(nn=1..,nb='',I2,'
     '                //''',nh='',I2,'',ne='',I4,''):'',32I2)')
     '                nb,nh,ne,(NVHE(nn,nb,nh,ne),nn=1,NNT(nb))
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDDO  !nh
                ENDIF
              ENDDO  !nb
            ENDDO  !ne
          ENDIF
        ENDDO !nonc (nc)
      ENDIF

      CALL EXITS('OPMAP')
      RETURN
 9999 CALL ERRORS('OPMAP',ERROR)
      CALL EXITS('OPMAP')
      RETURN 1
      END


