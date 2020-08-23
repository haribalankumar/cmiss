      SUBROUTINE DIVH1(IBT,IDO,IDRN,INP,NBH,NBJ,ne,NEELEM,
     '  NHE,NKJE,NKJ,NPF,NPLIST,NPNE,NPNODE,NPSTART,nr,
     '  NRE,NUNK,NVJE,NVJP,NW,nx,
     '  CE,SE,SP,XA,XE,XII,XP,NOCROSS,ERROR,*)

C#### Subroutine: DIVH1
C###  Description:
C###    DIVH1 subdivides Hermite tensor product element ne in the
C###    Xi(IDRN) direction at Xi(IDRN)=XII.

C**** IDR2 & IDR3 are the two orthogonal directions.
C**** NPSTART is first node number for new nodes
C**** NNIP(i1,i2,i3) records the element vertex number nn for each
C****   set of node position indices I1,I2,I3 where these are in the
C****   IDRN,IDR2 & IDR3 directions, respectively.
C**** IV(nn,idrn) are vertices of new element corresponding to newly
C****   created nodes
C**** DSDXI(ni) are derivatives of arc-length wrt Xi

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),IDRN,
     '  INP(NNM,NIM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),ne,
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),
     '  NPF(9,NFM),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPSTART,nr,NRE(NEM),NUNK(NKM,NJM,NPM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NW(NEM,3),nx
      REAL*8 CE(NMM,NEM),SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XII,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL NOCROSS
!     Local Variables
      INTEGER i1,i2,i3,IDR(5),IDR2,IDR3,IMAP(3,4,4,4),nb,nb1,NBC1,nc,
     '  ne_NEW,nh,nhx,NIT1,nj,njj1,njj2,nk,nm,nn,nn_ne,nn_nenew,
     '  NNIP(4,4,4),NODE,nonode,np,NPNE1(64,99),NPNE2(64,99),
     '  NPTNEW,NPTOLD,NUMI1,NUMI2,NUMI3,NVJE1(64,99,3*NJ_LOC_MX)
      REAL*8 REFINE_GETXIPOS,XI(3)
      LOGICAL ATCOLLAPSE(4,4,4),EXISTS

      DATA IDR/1,2,3,1,2/

      CALL ENTERS('DIVH1',*9999)

      nc=1 ! temporary cpb 22/11/94
      NPTNEW=NPSTART-1
      ne_NEW=NET(0)+1
      CALL ASSERT(ne_NEW.LE.NEM,'>>Too many elements -Increase NEM',
     '  ERROR,*9999)

      nb=NBJ(1,ne)
      NIT1=NIT(nb)
      NBC1=NBC(nb)
      IF(NIT1.EQ.1) THEN
        IDR2=0
        IDR3=0
      ELSEIF(NIT1.EQ.2) THEN
        IDR3=0
        IF(IDRN.EQ.1) THEN
          IDR2=2
        ELSE
          IDR2=1
        ENDIF
      ELSE
        IDR2=IDR(IDRN+1)
        IDR3=IDR(IDRN+2)
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' ne='',I5,'', nb='',I2,'', NIT(nb)='','
     '    //'I2,'', NNT(nb)='',I2,'', IDRN='',I1,'', IDR2='',I1,'
     '    //''', IDR3='',I1)') ne,nb,NIT(nb),NNT(nb),IDRN,IDR2,IDR3
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' XII='',D12.5)') XII
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DO nb1=1,NBFT
        DO nn=1,64
          NPNE1(nn,nb1)=0
          NPNE2(nn,nb1)=0
        ENDDO !nn
      ENDDO !nb1
      DO njj1=1,3
        DO njj2=1,NJ_LOC(njj1,0,nr)
          nj=NJ_LOC(njj1,njj2,nr)
          DO nb1=1,NBFT
            DO nn=1,NNT(nb1)
              NVJE1(nn,nb1,nj)=NVJE(nn,nb1,nj,ne)
            ENDDO !nn
          ENDDO !nb1
        ENDDO !njj2
      ENDDO !njj1

      CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '  NRE(ne),NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)

      DO njj1=1,3
        DO njj2=1,NJ_LOC(njj1,0,nr)
          nj=NJ_LOC(njj1,njj2,nr)
          NBJ(nj,ne_NEW)=NBJ(nj,ne)
        ENDDO !njj2
      ENDDO !njj1

C KAT 8Oct99: Interpolation may have different nodes for different nb
C             so we must loop over nb.

C      DO njj1=1,3
C        DO njj2=1,NJ_LOC(njj1,0,nr)
C          nj=NJ_LOC(njj1,njj2,nr)
C          nb=NBJ(nj,ne)
      DO nb=1,NBFT
        nb1=nb
        IF(NIT(nb).EQ.NIT1.AND.NBC(nb).EQ.NBC1) THEN

C KAT 29Aug00: MIPSpro f77 7.3.1.1m optimized code requires CALC_NNIP to
C         be able to access dummy var nb2 even though it isn't used, so
C         we can't give it a null pointer.  A constant 0 is provided
C         instead.
          CALL CALC_NNIP(0,IBT,IDRN,IDR2,IDR3,IMAP,INP,nb,0,NNIP,
     '      %VAL(0),0,NUMI1,NUMI2,NUMI3,1,ATCOLLAPSE,.FALSE.,.FALSE.,
     '      ERROR,*9999)

C***  Loop over Xi directions
          DO i3=1,NUMI3
            IF(NUMI3.EQ.1) THEN
              XI(3)=0.0d0
            ELSE
              XI(IDR3)=DBLE(i3-1)/DBLE(NUMI3-1)
            ENDIF
            DO i2=1,NUMI2
              IF(NUMI2.EQ.1) THEN
                XI(2)=0.0d0
              ELSE
                XI(IDR2)=DBLE(i2-1)/DBLE(NUMI2-1)
              ENDIF

              DO i1=1,NUMI1-1

C*** Get the xi location of the refine point
                XI(IDRN)=REFINE_GETXIPOS(i1,IBT,IDRN,nb,XII)

C*** Get the node number of the refine point
                CALL REFINE_FINDNODE(IBT,IDO,INP,nb,NBJ,ne,NKJ,NODE,
     '            NPLIST,NPNODE,NPTNEW,NPTOLD,nr,NUNK,NVJP,XE,XI,XP,
     '            EXISTS,ERROR,*9999)
                CALL NODE_CHANGE(NODE,.FALSE.,ERROR,*9999)

C*** Calculate the element connectivity array for the node
                IF(IBT(1,IDRN,nb).EQ.1.AND.IBT(2,IDRN,nb).EQ.2)
     '            THEN !Quadratic Lagrange
                  IF(XII.EQ.0.5d0) THEN
                    IF(i1.EQ.1) THEN
                      nn_ne=NNIP(2,i2,i3)
                      nn_nenew=0
                    ELSE
                      nn_ne=0
                      nn_nenew=NNIP(2,i2,i3)
                    ENDIF
                  ELSE IF(XII.GT.0.5d0) THEN
                    IF(i1.EQ.1) THEN
                      nn_ne=NNIP(3,i2,i3)
                    ELSE
                      nn_ne=0
                    ENDIF
                    nn_nenew=NNIP(i1,i2,i3)
                  ELSE
                    nn_ne=NNIP(i1+1,i2,i3)
                    IF(i1.EQ.2) THEN
                      nn_nenew=NNIP(1,i2,i3)
                    ELSE
                      nn_nenew=0
                    ENDIF
                  ENDIF
                ELSE IF(IBT(1,IDRN,nb).EQ.1.AND.IBT(2,IDRN,nb).EQ.3)
     '              THEN !Cubic Lagrange
                  IF(XII.EQ.(1.0d0/3.0d0)) THEN
                    IF(i1.EQ.3) THEN
                      nn_ne=0
                      nn_nenew=NNIP(2,i2,i3)
                    ELSE
                      nn_ne=NNIP(i1+1,i2,i3)
                      nn_nenew=0
                    ENDIF
                  ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                    IF(i1.EQ.1) THEN
                      nn_ne=NNIP(3,i2,i3)
                      nn_nenew=0
                    ELSE
                      nn_ne=0
                      nn_nenew=NNIP(i1,i2,i3)
                    ENDIF
                  ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                    IF(i1.EQ.1) THEN
                      nn_ne=NNIP(4,i2,i3)
                    ELSE
                      nn_ne=0
                    ENDIF
                    nn_nenew=NNIP(i1,i2,i3)
                  ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                    IF(i1.LE.2) THEN
                      nn_ne=NNIP(i1+2,i2,i3)
                    ELSE
                      nn_ne=0
                    ENDIF
                    IF(i1.GE.2) THEN
                      nn_nenew=NNIP(i1-1,i2,i3)
                    ELSE
                      nn_nenew=0
                    ENDIF
                  ELSE
                    nn_ne=NNIP(i1+1,i2,i3)
                    IF(i1.EQ.3) THEN
                      nn_nenew=NNIP(1,i2,i3)
                    ELSE
                      nn_nenew=0
                    ENDIF
                  ENDIF
                ELSE
                  nn_ne=NNIP(2,i2,i3)
                  nn_nenew=NNIP(1,i2,i3)
                ENDIF

C*** Put the node in the element
                IF(nn_ne.NE.0) NPNE1(nn_ne,nb1)=NODE
                IF(nn_nenew.NE.0) NPNE2(nn_nenew,nb1)=NODE

C*** Set the node properties
                DO njj1=1,3
                  DO njj2=1,NJ_LOC(njj1,0,nr)
                    nj=NJ_LOC(njj1,njj2,nr)
                    IF(NBJ(nj,ne).EQ.nb) THEN
                      CALL REFINE_SETNODE(i1,i2,i3,IBT,IDO,IDRN,INP,NBJ,
     '                  ne,ne_NEW,nj,NKJ,nn_ne,nn_nenew,NNIP,NODE,NPNE,
     '                  nr,NUMI1,NUNK,NVJE,
     '                  NVJP,SE,SP,XE,XI,XII,XP,
     '                  EXISTS,NOCROSS,ERROR,*9999)
                    ENDIF !nb
                  ENDDO !njj2
                ENDDO !njj1

              ENDDO !i1

C*** Fill in the old nodes into the old and new elements and set the
C*** local node properties for the existing nodes (versions and
C*** derivatives).
              IF(IBT(1,IDRN,nb).EQ.1.AND.IBT(2,IDRN,nb).EQ.2)
     '          THEN !Quadratic Lagrange
                IF(XII.EQ.0.5d0) THEN
                  NPNE1(NNIP(1,i2,i3),nb1)=NPNE(NNIP(1,i2,i3),nb,ne)
                  NPNE1(NNIP(3,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
                  NPNE2(NNIP(1,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
                  NPNE2(NNIP(3,i2,i3),nb1)=NPNE(NNIP(3,i2,i3),nb,ne)
C              DO nb1=1,NBFT
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      IF(NBJ(nj,ne).EQ.nb) THEN
                        NVJE(NNIP(1,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(2,i2,i3),nb1,nj)
                        NVJE(NNIP(3,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(3,i2,i3),nb1,nj)
                      ENDIF !nb
                    ENDDO !njj2
                  ENDDO !njj1
C              ENDDO !nb1
                ELSE IF(XII.GT.0.5d0) THEN
                  NPNE1(NNIP(1,i2,i3),nb1)=NPNE(NNIP(1,i2,i3),nb,ne)
                  NPNE1(NNIP(2,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
                  NPNE2(NNIP(3,i2,i3),nb1)=NPNE(NNIP(3,i2,i3),nb,ne)
C              DO nb1=1,NBFT
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      IF(NBJ(nj,ne).EQ.nb) THEN
                        NVJE(NNIP(3,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(3,i2,i3),nb1,nj)
                      ENDIF !nb
                    ENDDO !njj2
                  ENDDO !njj1
C              ENDDO !nb1
                ELSE
                  NPNE1(NNIP(1,i2,i3),nb1)=NPNE(NNIP(1,i2,i3),nb,ne)
                  NPNE2(NNIP(2,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
                  NPNE2(NNIP(3,i2,i3),nb1)=NPNE(NNIP(3,i2,i3),nb,ne)
C              DO nb1=1,NBFT
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      IF(NBJ(nj,ne).EQ.nb) THEN
                        NVJE(NNIP(2,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(2,i2,i3),nb1,nj)
                        NVJE(NNIP(3,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(3,i2,i3),nb1,nj)
                      ENDIF !nb
                    ENDDO !njj2
                  ENDDO !njj1
C              ENDDO !nb1
                ENDIF
              ELSE IF(IBT(1,IDRN,nb).EQ.1.AND.
     '            IBT(2,IDRN,nb).EQ.3) THEN !Cubic Lagrange
                IF(XII.EQ.(1.0d0/3.0d0)) THEN
                  NPNE1(NNIP(1,i2,i3),nb1)=NPNE(NNIP(1,i2,i3),nb,ne)
                  NPNE1(NNIP(4,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
                  NPNE2(NNIP(1,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
                  NPNE2(NNIP(3,i2,i3),nb1)=NPNE(NNIP(3,i2,i3),nb,ne)
                  NPNE2(NNIP(4,i2,i3),nb1)=NPNE(NNIP(4,i2,i3),nb,ne)
C              DO nb1=1,NBFT
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      IF(NBJ(nj,ne).EQ.nb) THEN
                        NVJE(NNIP(1,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(2,i2,i3),nb1,nj)
                        NVJE(NNIP(3,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(3,i2,i3),nb1,nj)
                        NVJE(NNIP(4,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(4,i2,i3),nb1,nj)
                      ENDIF !nb
                    ENDDO !njj2
                  ENDDO !njj1
C              ENDDO !nb1
                ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                  NPNE1(NNIP(1,i2,i3),nb1)=NPNE(NNIP(1,i2,i3),nb,ne)
                  NPNE1(NNIP(2,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
                  NPNE1(NNIP(4,i2,i3),nb1)=NPNE(NNIP(3,i2,i3),nb,ne)
                  NPNE2(NNIP(1,i2,i3),nb1)=NPNE(NNIP(3,i2,i3),nb,ne)
                  NPNE2(NNIP(4,i2,i3),nb1)=NPNE(NNIP(4,i2,i3),nb,ne)
C              DO nb1=1,NBFT
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      IF(NBJ(nj,ne).EQ.nb) THEN
                        NVJE(NNIP(1,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(3,i2,i3),nb1,nj)
                        NVJE(NNIP(4,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(4,i2,i3),nb1,nj)
                      ENDIF !nb
                    ENDDO !njj2
                  ENDDO !njj1
C              ENDDO !nb1
                ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                  NPNE1(NNIP(1,i2,i3),nb1)=NPNE(NNIP(1,i2,i3),nb,ne)
                  NPNE1(NNIP(2,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
                  NPNE1(NNIP(3,i2,i3),nb1)=NPNE(NNIP(3,i2,i3),nb,ne)
                  NPNE2(NNIP(4,i2,i3),nb1)=NPNE(NNIP(4,i2,i3),nb,ne)
C              DO nb1=1,NBFT
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      IF(NBJ(nj,ne).EQ.nb) THEN
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP(4,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(4,i2,i3),nb1,nj)
                      ENDIF !nb
                    ENDDO !njj2
                  ENDDO !njj1
C              ENDDO !nb1
                ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                  NPNE1(NNIP(1,i2,i3),nb1)=NPNE(NNIP(1,i2,i3),nb,ne)
                  NPNE1(NNIP(2,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
                  NPNE2(NNIP(3,i2,i3),nb1)=NPNE(NNIP(3,i2,i3),nb,ne)
                  NPNE2(NNIP(4,i2,i3),nb1)=NPNE(NNIP(4,i2,i3),nb,ne)
C              DO nb1=1,NBFT
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      IF(NBJ(nj,ne).EQ.nb) THEN
                        NVJE(NNIP(3,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(3,i2,i3),nb1,nj)
                        NVJE(NNIP(4,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(4,i2,i3),nb1,nj)
                      ENDIF !nb
                    ENDDO !njj2
                  ENDDO !njj1
C              ENDDO !nb1
                ELSE
                  NPNE1(NNIP(1,i2,i3),nb1)=NPNE(NNIP(1,i2,i3),nb,ne)
                  NPNE2(NNIP(2,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
                  NPNE2(NNIP(3,i2,i3),nb1)=NPNE(NNIP(3,i2,i3),nb,ne)
                  NPNE2(NNIP(4,i2,i3),nb1)=NPNE(NNIP(4,i2,i3),nb,ne)
C              DO nb1=1,NBFT
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      IF(NBJ(nj,ne).EQ.nb) THEN
                        NVJE(NNIP(2,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(2,i2,i3),nb1,nj)
                        NVJE(NNIP(3,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(3,i2,i3),nb1,nj)
                        NVJE(NNIP(4,i2,i3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP(4,i2,i3),nb1,nj)
                      ENDIF !nb
                    ENDDO !njj2
                  ENDDO !njj1
C              ENDDO !nb1
                ENDIF
              ELSE
                NPNE1(NNIP(1,i2,i3),nb1)=NPNE(NNIP(1,i2,i3),nb,ne)
                NPNE2(NNIP(2,i2,i3),nb1)=NPNE(NNIP(2,i2,i3),nb,ne)
C            DO nb1=1,NBFT
                DO njj1=1,3
                  DO njj2=1,NJ_LOC(njj1,0,nr)
                    nj=NJ_LOC(njj1,njj2,nr)
                    IF(NBJ(nj,ne).EQ.nb) THEN
                      NVJE(NNIP(2,i2,i3),nb1,nj,ne_NEW)=
     '                  NVJE1(NNIP(2,i2,i3),nb1,nj)
                    ENDIF !nb
                  ENDDO !njj2
                ENDDO !njj1
C            ENDDO !nb1
              ENDIF
            ENDDO !i2
          ENDDO !i3

C        ENDDO !njj2
C      ENDDO !njj1

        ENDIF !NIT
      ENDDO !nb

      NPSTART=NPT(0)+1

      DO nb1=1,NBFT
        DO nn=1,NNT(nb1)
C KAT 22Jul99: Interpolation may have different nodes for different nj
          NPNE(nn,nb1,ne)=NPNE1(nn,nb1)
          NPNE(nn,nb1,ne_NEW)=NPNE2(nn,nb1)
        ENDDO !nn
      ENDDO !nb1

      DO njj1=1,3
        DO njj2=1,NJ_LOC(njj1,0,nr)
          nj=NJ_LOC(njj1,njj2,nr)
          nb1=NBJ(nj,ne)
          DO nn=1,NNT(nb1)
            DO nk=1,NKT(nn,nb1)
              NKJE(nk,nn,nj,ne_NEW)=NKJE(nk,nn,nj,ne)
            ENDDO !nk
          ENDDO !nn
        ENDDO !njj2
      ENDDO !njj1

      CALL REFINE_SETSCALEFACTORS(IBT,IDO,IDRN,ne,ne_NEW,NIT1,SE,XII,
     '  ERROR,*9999)

      NHE(ne_NEW)=NHE(ne)
      NW(ne_NEW,1)=NW(ne,1)
      NW(ne_NEW,2)=NW(ne,2)
      NRE(ne_NEW)=NRE(ne)
      DO nm=1,NMM
        CE(nm,ne_NEW)=CE(nm,ne)
      ENDDO !nm

      IF(CALL_EQUA.OR.CALL_FIT.OR.CALL_OPTI) THEN
        DO nhx=1,NHE(ne)
          nh=NH_LOC(nhx,nx)
          DO nc=1,NCM
            NBH(nh,nc,ne_NEW)=NBH(nh,nc,ne)
          ENDDO !nc
        ENDDO !nhx
      ENDIF !call_equa etc

      NET(0)=ne_NEW   !is new highest element # in mesh
      NET(nr)=ne_NEW  !is new highest element # in region nr
      CALL ASSERT(NET(nr).LE.NEM,'>>Increase NEM',ERROR,*9999)
      NEELEM(0,nr)=NEELEM(0,nr)+1 !is new #elements in region nr
      CALL ASSERT(NEELEM(0,nr).LE.NE_R_M,'>>Increase NE_R_M',
     '  ERROR,*9999)
      NEELEM(NEELEM(0,nr),nr)=ne_NEW !is new element number

!     Update element list in each group
      CALL UPGREL(ne,ne_NEW,ERROR,*9999)

! Check that the global nodes of new element ne_NEW are in the
! node list for the current region
      nb1=NBJ(1,ne_NEW)
      DO nn=1,NNT(nb1)
        IF(NPNODE(0,nr).LE.NP_R_M) THEN
          np=NPNE(nn,nb1,ne_NEW)
          EXISTS=.FALSE.
          DO nonode=1,NPNODE(0,nr)
            IF(NPNODE(nonode,nr).eq.np) EXISTS=.TRUE.
          ENDDO
          IF(.NOT.EXISTS) THEN
            NPNODE(0,nr)=NPNODE(0,nr)+1
            IF(NPNODE(0,nr).LE.NP_R_M) NPNODE(NPNODE(0,nr),nr)=np
          ENDIF
        ENDIF
      ENDDO
      CALL ASSERT(NPNODE(0,nr).LE.NP_R_M,'>>Increase NP_R_M',
     '  ERROR,*9999)

C*** Diagnostic output
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' ne='',I5,'', ne_new='',I5)') ne,ne_new
        DO nb=1,NBFT
          WRITE(OP_STRING,'('' NPNE(nn,'',I2,'',ne    ):'','
     '      //'10(1X,I5),/:(20X,10(1X,I5)))') nb,(NPNE(nn,nb,ne),
     '      nn=1,NNT(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NPNE(nn,'',I2,'',ne_new):'','
     '      //'10(1X,I5),/:(20X,10(1X,I5)))') nb,(NPNE(nn,nb,ne_new),
     '      nn=1,NNT(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nb
CC$      call mp_unsetlock()
      ENDIF !dop

      CALL EXITS('DIVH1')
      RETURN
 9999 CALL ERRORS('DIVH1',ERROR)
      CALL EXITS('DIVH1')
      RETURN 1
      END


