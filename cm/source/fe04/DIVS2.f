      SUBROUTINE DIVS2(IBT,IDO,IDRN,INP,NBH,NBJ,ne,NEELEM,
     '  NHE,NKJE,NKJ,NPF,NPLIST,NPNE,NPNODE,NPSTART,nr,
     '  NRE,NUNK,NVJE,NVJP,NW,nx,
     '  CE,SE,SP,XA,XE,XII,XP,NOCROSS,ERROR,*)

C#### Subroutine: DIVS2
C###  Description:
C###    DIVS2 subdivides the sector element ne in the Xi(IDRN)
C###    direction at Xi(IDRN)=XII.

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
     '  NPNODE(0:NP_R_M,0:NRM),NPSTART,
     '  nr,NRE(NEM),NUNK(NKM,NJM,NPM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NW(NEM,3),nx
      REAL*8 CE(NMM,NEM),SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XII,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL NOCROSS
!     Local Variables
      INTEGER CASE,COLLAPSED_XI,i1,i2,ii1,ii2,i3,ii3,IBT_CHECK(2,3),
     '  IDR(5),IDR2,IDR3,IMAP(3,4,4,4),MATCH,nb,nb1,nbb,NBC_CHECK,nc,
     '  ne_NEW,nh,nhx,ni,NITB,nj,njj1,njj2,nk,nm,nn,nn_ne,nn_nenew,
     '  NNIP1(4,4,4),NNIP2(4,4,4),NODE,nonode,np,NPNE1(64),NPNE2(64),
     '  NPTNEW,NPTOLD,NUM_COLLAPSED,NUMI1,NUMI2,NUMI3,PERP_XI,
     '  NVJE1(64,99,3*NJ_LOC_MX)
      REAL*8 REFINE_GETXIPOS,XI(3)
      LOGICAL ATCOLLAPSE(4,4,4),EXISTS,FOUND

      DATA IDR/1,2,3,1,2/

      CALL ENTERS('DIVS2',*9999)

C!!! It is assumed that the same basis function is used in all geometric
C!!! directions.

C*** With sector elements there are three cases to consider when
C*** refining. The first case comes about if we are refining in a
C*** direction perpendicular to a collapsed Xi direction (as defined
C*** by IBT(3,ni,nb)) In this case we need to create another element
C*** of the same type as the current element.
C*** The second case occurs when we are refining in the collpased
C*** direction. In this case we need to create a new element that is
C*** the uncollapsed version of the current element.
C*** The third case occurs when we are refining in a direction that
C*** is neither collapsed nor perpendicular to a collapsed xi
C*** direction. This case is the standard refine as per divh1.

C cpb 30/11/96 I will separately code each case rather than trying to
C be clever and think of a neat way to handle everything.

      IF(NET(0)+1.LE.NEM) THEN
        nc=1 ! temporary cpb 22/11/94
C       LLOOSE_TOL=DSQRT(LOOSE_TOL)
        NPTNEW=NPSTART-1
        ne_NEW=NET(0)+1
        nb=NBJ(1,ne)
        NITB=NIT(nb)
        IF(NITB.EQ.2) THEN
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
        NUM_COLLAPSED=0
        DO ni=1,NITB
          IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
            NUM_COLLAPSED=NUM_COLLAPSED+1
            COLLAPSED_XI=ni
            PERP_XI=IBT(3,ni,nb)
          ENDIF
        ENDDO !ni
        CALL ASSERT(NUM_COLLAPSED.GT.0,'>>Should use DIVH1',ERROR,*9999)
        IF(IDRN.EQ.PERP_XI) THEN
          CASE=2
        ELSE IF(IBT(1,IDRN,nb).EQ.5.OR.IBT(1,IDRN,nb).EQ.6) THEN
          CASE=1
        ELSE
          CASE=3
        ENDIF

        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' ne='',I5,'', nb='',I2,'', NIT(nb)='','
     '      //'I2,'', NNT(nb)='',I2,'', IDRN='',I1,'', IDR2='',I1,'
     '      //''', IDR3='',I1)') ne,nb,NIT(nb),NNT(nb),IDRN,IDR2,IDR3
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' XII='',D12.5)') XII
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' CASE='',I1)') CASE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        DO nn=1,64
          NPNE1(nn)=0
          NPNE2(nn)=0
        ENDDO !nn
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
     '    NRE(ne),NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,
     '    ERROR,*9999)

        IF(CASE.EQ.1) THEN !Case 1

          CALL CALC_NNIP(COLLAPSED_XI,IBT,IDRN,IDR2,IDR3,IMAP,INP,
     '      nb,nb,NNIP1,NNIP2,NUM_COLLAPSED,NUMI1,NUMI2,NUMI3,PERP_XI,
     '      ATCOLLAPSE,.FALSE.,.TRUE.,ERROR,*9999)

          DO njj1=1,3
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              NBJ(nj,ne_NEW)=NBJ(nj,ne)
            ENDDO !njj2
          ENDDO !njj1

C*** Loop over Xi directions
          DO i3=1,NUMI3
            IF(NUMI3.EQ.1) THEN
              XI(3)=0.0d0
            ELSE
              XI(IDR3)=DBLE(i3-1)/DBLE(NUMI3-1)
            ENDIF
            DO i2=1,NUMI2
              XI(IDR2)=DBLE(i2-1)/DBLE(NUMI2-1)
              DO i1=1,NUMI1-1
                ii1=IMAP(1,i1,i2,i3)
                ii2=IMAP(2,i1,i2,i3)
                ii3=IMAP(3,i1,i2,i3)
                IF(.NOT.ATCOLLAPSE(i1,i2,i3)) THEN

C*** Get the xi position
                  XI(IDRN)=REFINE_GETXIPOS(ii1,IBT,IDRN,nb,XII)

C*** Get the node number of the refine point
                  CALL REFINE_FINDNODE(IBT,IDO,INP,nb,NBJ,ne,NKJ,NODE,
     '              NPLIST,NPNODE,NPTNEW,NPTOLD,nr,NUNK,NVJP,XE,XI,XP,
     '              EXISTS,ERROR,*9999)

C*** Calculate the element arrays for the node
                  IF((IBT(1,IDRN,nb).EQ.1.OR.IBT(1,IDRN,nb).EQ.5.OR.
     '              IBT(1,IDRN,nb).EQ.6).AND.IBT(2,IDRN,nb).EQ.2)
     '              THEN !Quadratic Lagrange
                    IF(XII.EQ.0.5d0) THEN
                      IF(ii1.EQ.1) THEN
                        nn_ne=NNIP1(2,ii2,ii3)
                        nn_nenew=0
                      ELSE
                        nn_ne=0
                        nn_nenew=NNIP1(2,ii2,ii3)
                      ENDIF
                    ELSE IF(XII.GT.0.5d0) THEN
                      IF(ii1.EQ.1) THEN
                        nn_ne=NNIP1(3,ii2,ii3)
                      ELSE
                        nn_ne=0
                      ENDIF
                      nn_nenew=NNIP1(ii1,ii2,ii3)
                    ELSE
                      nn_ne=NNIP1(ii1+1,ii2,ii3)
                      IF(ii1.EQ.2) THEN
                        nn_nenew=NNIP1(1,ii2,ii3)
                      ELSE
                        nn_nenew=0
                      ENDIF
                    ENDIF
                  ELSE IF((IBT(1,IDRN,nb).EQ.1.OR.
     '                IBT(1,IDRN,nb).EQ.5.OR.IBT(1,IDRN,nb).EQ.6).AND.
     '                IBT(2,IDRN,nb).EQ.3) THEN !Cubic Lagrange
                    IF(XII.EQ.(1.0d0/3.0d0)) THEN
                      IF(ii1.EQ.3) THEN
                        nn_ne=0
                        nn_nenew=NNIP1(2,ii2,ii3)
                      ELSE
                        nn_ne=NNIP1(ii1+1,ii2,ii3)
                        nn_nenew=0
                      ENDIF
                    ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                      IF(ii1.EQ.1) THEN
                        nn_ne=NNIP1(3,ii2,ii3)
                        nn_nenew=0
                      ELSE
                        nn_ne=0
                        nn_nenew=NNIP1(ii1,ii2,ii3)
                      ENDIF
                    ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                      IF(ii1.EQ.1) THEN
                        nn_ne=NNIP1(4,ii2,ii3)
                      ELSE
                        nn_ne=0
                      ENDIF
                      nn_nenew=NNIP1(ii1,ii2,ii3)
                    ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                      IF(ii1.LE.2) THEN
                        nn_ne=NNIP1(ii1+2,ii2,ii3)
                      ELSE
                        nn_ne=0
                      ENDIF
                      IF(ii1.GE.2) THEN
                        nn_nenew=NNIP1(ii1-1,ii2,ii3)
                      ELSE
                        nn_nenew=0
                      ENDIF
                    ELSE
                      nn_ne=NNIP1(ii1+1,ii2,ii3)
                      IF(ii1.EQ.3) THEN
                        nn_nenew=NNIP1(1,ii2,ii3)
                      ELSE
                        nn_nenew=0
                      ENDIF
                    ENDIF
                  ELSE
                    nn_ne=NNIP1(2,ii2,ii3)
                    nn_nenew=NNIP1(1,ii2,ii3)
                  ENDIF

C*** Put the node in the element
                IF(nn_ne.NE.0) NPNE1(nn_ne)=NODE
                IF(nn_nenew.NE.0) NPNE2(nn_nenew)=NODE

C*** Set the node properties
C KAT 22Jul99: Interpolation may have different nodes for different nj.
C              This nj loop should include CALC_NNIP (see DIVH1).
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      CALL REFINE_SETNODE(ii1,ii2,ii3,IBT,IDO,IDRN,
     '                  INP,NBJ,ne,ne_new,nj,NKJ,nn_ne,nn_nenew,
     '                  NNIP1,NODE,NPNE,nr,NUMI1,
     '                  NUNK,NVJE,NVJP,SE,SP,
     '                  XE,XI,XII,XP,EXISTS,NOCROSS,
     '                  ERROR,*9999)
                    ENDDO !njj2
                  ENDDO !njj1

                ENDIF

              ENDDO !i1

C*** Fill in the old nodes into the old and new elements and set the
C*** local node properties for the exisiting nodes (versions and
C*** derivatives).
              IF((IBT(1,IDRN,nb).EQ.1.OR.IBT(1,IDRN,nb).EQ.5.OR.
     '          IBT(1,IDRN,nb).EQ.6).AND.IBT(2,IDRN,nb).EQ.2)
     '          THEN !Quadratic Lagrange
                IF(XII.EQ.0.5d0) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(3,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(1,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(1,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE IF(XII.GT.0.5d0) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(2,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ENDIF
              ELSE IF((IBT(1,IDRN,nb).EQ.1.OR.
     '            IBT(1,IDRN,nb).EQ.5.OR.IBT(1,IDRN,nb).EQ.6).AND.
     '            IBT(2,IDRN,nb).EQ.3) THEN !Cubic Lagrange
                IF(XII.EQ.(1.0d0/3.0d0)) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(4,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(1,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(1,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(4,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(4,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(1,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(1,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(4,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(4,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(4,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(2,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(4,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ENDIF
              ELSE
                NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                NPNE2(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                DO nb1=1,NBFT
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      NVJE(NNIP1(2,ii2,ii3),nb1,nj,ne_NEW)=
     '                  NVJE1(NNIP1(2,ii2,ii3),nb1,nj)
                    ENDDO !njj2
                  ENDDO !njj1
                ENDDO !nb1
              ENDIF
            ENDDO !i2
          ENDDO !i3

C*** Update scaling factors in current element & new element
          CALL REFINE_SETSCALEFACTORS(IBT,IDO,IDRN,ne,ne_new,NITB,SE,
     '      XII,ERROR,*9999)

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' NPNE1(nn):'',10(1X,I5),'
     '        //'/:(11X,10(1X,I5)))') (NPNE1(nn),nn=1,NNT(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NPNE2(nn):'',10(1X,I5),'
     '        //'/:(11X,10(1X,I5)))') (NPNE2(nn),nn=1,NNT(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

          DO nb1=1,NBFT
            DO nn=1,NNT(nb1)
              NPNE(nn,nb1,ne)=NPNE1(nn)
              NPNE(nn,nb1,ne_NEW)=NPNE2(nn)
C              DO nk=1,NKT(nn,nb1)
C                NKE(nk,nn,nb1,ne_NEW)=NKE(nk,nn,nb1,ne)
C              ENDDO !nk
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

C*** Form other arrays for new element ne_NEW
          NHE(ne_NEW)=NHE(ne)
          NW(ne_NEW,1)=NW(ne,1)
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

          NPSTART=NPT(0)+1
          NET(0)=ne_NEW !is new highest element # in mesh
          NET(nr)=ne_NEW !is new highest element # in region nr
          CALL ASSERT(NET(nr).LE.NEM,
     '      '>>Too many elements - increase NEM',ERROR,*9999)
          NEELEM(0,nr)=NEELEM(0,nr)+1 !is new #elements in region nr
          NEELEM(NEELEM(0,nr),nr)=ne_NEW !is new element number

C*** Update element list in each group
          CALL UPGREL(ne,ne_NEW,ERROR,*9999)

C*** Check that the global nodes of new element ne_NEW are in the
C*** node list for the current region
          DO nn=1,NNT(nb)
            IF(NPNODE(0,nr).LE.NP_R_M) THEN
              np=NPNE(nn,nb,ne_NEW)
              FOUND=.FALSE.
              nonode=1
              DO WHILE(.NOT.FOUND.AND.nonode.LE.NPNODE(0,nr))
                IF(NPNODE(nonode,nr).eq.np) THEN
                  FOUND=.TRUE.
                ELSE
                  nonode=nonode+1
                ENDIF
              ENDDO
              IF(.NOT.FOUND) THEN
                NPNODE(0,nr)=NPNODE(0,nr)+1
                IF(NPNODE(0,nr).LE.NP_R_M) NPNODE(NPNODE(0,nr),nr)=np
              ENDIF
            ENDIF
          ENDDO !nn
          CALL ASSERT(NPNODE(0,nr).LE.NP_R_M,'>>Increase NP_R_M',
     '      ERROR,*9999)

        ELSE IF(CASE.EQ.2) THEN !Case 2

C*** Need to check that the uncollpased basis function has been set
C*** up.
          DO ni=1,NITB
            IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
              IF(IBT(2,ni,nb).EQ.4) THEN !Hermite
                IBT_CHECK(1,ni)=2
                IBT_CHECK(2,ni)=1
              ELSE
                IBT_CHECK(1,ni)=1
                IBT_CHECK(2,ni)=IBT(2,ni,nb)
              ENDIF
            ELSE
              IBT_CHECK(1,ni)=IBT(1,ni,nb)
              IF(IBT(1,ni,nb).EQ.2) THEN
                IBT_CHECK(2,ni)=1
              ELSE
                IBT_CHECK(2,ni)=IBT(2,ni,nb)
              ENDIF
            ENDIF
          ENDDO !ni
          IF(NBC(nb).EQ.2) THEN
            NBC_CHECK=1
          ELSE IF(NBC(nb).EQ.6) THEN
            NBC_CHECK=5
          ELSE
            ERROR='>>Invalid NBC'
            GOTO 9999
          ENDIF
          FOUND=.FALSE.
          nbb=1
C MLB July 2001 - This should definitely be NBFT not NBT. IBT is
C         dimensioned to NBFT and looping to NBT causes BEM sector
C         refinement to crash.
C          DO WHILE(.NOT.FOUND.AND.nbb.LE.NBT)
          DO WHILE(.NOT.FOUND.AND.nbb.LE.NBFT)
            MATCH=0
            DO ni=1,NITB
              IF(IBT(1,ni,nbb).EQ.IBT_CHECK(1,ni).AND.IBT(2,ni,nbb).EQ.
     '          IBT_CHECK(2,ni).AND.NBC(nbb).EQ.NBC_CHECK) MATCH=MATCH+1
            ENDDO !ni
            IF(MATCH.EQ.NITB) THEN
              FOUND=.TRUE.
            ELSE
              nbb=nbb+1
            ENDIF
          ENDDO ! not found
          CALL ASSERT(FOUND,'>>Need to set up uncollapsed basis',
     '      ERROR,*9999)

          CALL CALC_NNIP(COLLAPSED_XI,IBT,IDRN,IDR2,IDR3,IMAP,INP,
     '      nb,nbb,NNIP1,NNIP2,NUM_COLLAPSED,NUMI1,NUMI2,NUMI3,PERP_XI,
     '      ATCOLLAPSE,.TRUE.,.TRUE.,ERROR,*9999)

          DO njj1=1,3
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              NBJ(nj,ne_NEW)=nbb
            ENDDO !njj2
          ENDDO !njj1

C*** Initialise NKJE for the new element. Set those elements that need
C*** updating later on.
          DO njj1=1,3
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              nb1=NBJ(nj,ne)
              DO nn=1,NNT(nb1)
                DO nk=1,NKT(nn,nb1)
                  NKJE(nk,nn,nj,ne_NEW)=1
                ENDDO !nk
              ENDDO !nn
            ENDDO !njj2
          ENDDO !njj1
C          DO nb1=1,NBFT
C            DO nn=1,NNT(nb1)
C              DO nk=1,NKT(nn,nb1)
C                NKE(nk,nn,nb1,ne_NEW)=1
C              ENDDO !nk
C            ENDDO !nn
C          ENDDO !nb1

C*** Loop over Xi directions
          DO i3=1,NUMI3
            IF(NUMI3.EQ.1) THEN
              XI(3)=0.0d0
            ELSE
              XI(IDR3)=DBLE(i3-1)/DBLE(NUMI3-1)
            ENDIF
            DO i2=1,NUMI2
              XI(IDR2)=DBLE(i2-1)/DBLE(NUMI2-1)
              DO i1=1,NUMI1-1

C*** Get the xi position
                XI(IDRN)=REFINE_GETXIPOS(i1,IBT,IDRN,nb,XII)

C*** Get the node number of the refine point
                CALL REFINE_FINDNODE(IBT,IDO,INP,nb,NBJ,ne,NKJ,NODE,
     '            NPLIST,NPNODE,NPTNEW,NPTOLD,nr,NUNK,NVJP,XE,XI,XP,
     '            EXISTS,ERROR,*9999)

C*** Calculate the element arrays for the node
                IF(IBT(1,COLLAPSED_XI,nb).EQ.5) THEN
                  IF(IBT(1,IDRN,nb).EQ.1.AND.
     '              IBT(2,IDRN,nb).EQ.2) THEN !Quadratic Lagrange
                    IF(XII.EQ.0.5d0) THEN
                      IF(i1.EQ.1) THEN
                        nn_ne=NNIP1(2,i2,i3)
                        nn_nenew=0
                      ELSE
                        nn_ne=0
                        nn_nenew=NNIP2(2,i2,i3)
                      ENDIF
                    ELSE IF(XII.GT.0.5d0) THEN
                      IF(i1.EQ.1) THEN
                        nn_ne=NNIP1(3,i2,i3)
                      ELSE
                        nn_ne=0
                      ENDIF
                      nn_nenew=NNIP2(i1,i2,i3)
                    ELSE
                      nn_ne=NNIP1(i1+1,i2,i3)
                      IF(i1.EQ.2) THEN
                        nn_nenew=NNIP2(1,i2,i3)
                      ELSE
                        nn_nenew=0
                      ENDIF
                    ENDIF
                  ELSE IF(IBT(1,IDRN,nb).EQ.1.AND.
     '              IBT(2,IDRN,nb).EQ.3) THEN !Cubic Lagrange
                    IF(XII.EQ.(1.0d0/3.0d0)) THEN
                      IF(i1.EQ.3) THEN
                        nn_ne=0
                        nn_nenew=NNIP2(2,i2,i3)
                      ELSE
                        nn_ne=NNIP1(i1+1,i2,i3)
                        nn_nenew=0
                      ENDIF
                    ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                      IF(i1.EQ.1) THEN
                        nn_ne=NNIP1(3,i2,i3)
                        nn_nenew=0
                      ELSE
                        nn_ne=0
                        nn_nenew=NNIP2(i1,i2,i3)
                      ENDIF
                    ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                      IF(i1.EQ.1) THEN
                        nn_ne=NNIP1(4,i2,i3)
                      ELSE
                        nn_ne=0
                      ENDIF
                      nn_nenew=NNIP2(i1,i2,i3)
                    ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                      IF(i1.LE.2) THEN
                        nn_ne=NNIP1(i1+2,i2,i3)
                      ELSE
                        nn_ne=0
                      ENDIF
                      IF(i1.GE.2) THEN
                        nn_nenew=NNIP2(i1-1,i2,i3)
                      ELSE
                        nn_nenew=0
                      ENDIF
                    ELSE
                      nn_ne=NNIP1(i1+1,i2,i3)
                      IF(i1.EQ.3) THEN
                        nn_nenew=NNIP2(1,i2,i3)
                      ELSE
                        nn_nenew=0
                      ENDIF
                    ENDIF
                  ELSE
                    nn_ne=NNIP1(i1+1,i2,i3)
                    nn_nenew=NNIP2(i1,i2,i3)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        DO nk=1,NKT(nn_ne,nb)
                          NKJE(nk,nn_nenew,nj,ne_NEW)=
     '                      NKJE(nk,nn_ne,nj,ne)
                        ENDDO !nk
                      ENDDO !njj2
                    ENDDO !njj1
C                    DO nk=1,NKT(nn_ne,nb)
C                      NKE(nk,nn_nenew,nbb,ne_NEW)=NKE(nk,nn_ne,nb,ne)
C                    ENDDO !nk
                  ENDIF
                ELSE !IBT(1,COLLAPSED_XI,nb)=6
                  IF(IBT(1,IDRN,nb).EQ.1.AND.
     '              IBT(2,IDRN,nb).EQ.2) THEN !Quadratic Lagrange
                    IF(XII.EQ.0.5d0) THEN
                      IF(i1.EQ.1) THEN
                        nn_ne=0
                        nn_nenew=NNIP2(2,i2,i3)
                      ELSE
                        nn_ne=NNIP1(2,i2,i3)
                        nn_nenew=0
                      ENDIF
                    ELSE IF(XII.GT.0.5d0) THEN
                      nn_ne=NNIP1(i1,i2,i3)
                      IF(i1.EQ.1) THEN
                        nn_nenew=NNIP2(3,i2,i3)
                      ELSE
                        nn_nenew=0
                      ENDIF
                    ELSE
                      IF(i1.EQ.2) THEN
                        nn_ne=NNIP1(1,i2,i3)
                      ELSE
                        nn_ne=0
                      ENDIF
                      nn_nenew=NNIP2(i1+1,i2,i3)
                    ENDIF
                  ELSE IF(IBT(1,IDRN,nb).EQ.1.AND.
     '              IBT(2,IDRN,nb).EQ.3) THEN !Cubic Lagrange
                    IF(XII.EQ.(1.0d0/3.0d0)) THEN
                      IF(i1.EQ.3) THEN
                        nn_ne=NNIP1(2,i2,i3)
                        nn_nenew=0
                      ELSE
                        nn_ne=0
                        nn_nenew=NNIP2(i1+1,i2,i3)
                      ENDIF
                    ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                      IF(i1.EQ.1) THEN
                        nn_ne=0
                        nn_nenew=NNIP2(3,i2,i3)
                      ELSE
                        nn_ne=NNIP1(i1,i2,i3)
                        nn_nenew=0
                      ENDIF
                    ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                      nn_ne=NNIP1(i1,i2,i3)
                      IF(i1.EQ.1) THEN
                        nn_nenew=NNIP2(4,i2,i3)
                      ELSE
                        nn_nenew=0
                      ENDIF
                    ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                      IF(i1.GE.2) THEN
                        nn_ne=NNIP1(i1-1,i2,i3)
                      ELSE
                        nn_ne=0
                      ENDIF
                      IF(i1.LE.2) THEN
                        nn_nenew=NNIP2(i1+2,i2,i3)
                      ELSE
                        nn_nenew=0
                      ENDIF
                    ELSE
                      IF(i1.EQ.3) THEN
                        nn_ne=NNIP1(1,i2,i3)
                      ELSE
                        nn_ne=0
                      ENDIF
                      nn_nenew=NNIP2(i1+1,i2,i3)
                    ENDIF
                  ELSE
                    nn_ne=NNIP1(i1,i2,i3)
                    nn_nenew=NNIP2(i1+1,i2,i3)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        DO nk=1,NKT(nn_ne,nb)
                          NKJE(nk,nn_nenew,nj,ne_NEW)=
     '                      NKJE(nk,nn_ne,nj,ne)
                        ENDDO !nk
                      ENDDO !njj2
                    ENDDO !njj1
C                    DO nk=1,NKT(nn_ne,nb)
C                      NKE(nk,nn_nenew,nbb,ne_NEW)=NKE(nk,nn_ne,nb,ne)
C                    ENDDO !nk
                  ENDIF
                ENDIF

C*** Put the node in the element
                IF(nn_ne.NE.0) NPNE1(nn_ne)=NODE
                IF(nn_nenew.NE.0) NPNE2(nn_nenew)=NODE

C*** Set the node properties
C KAT 22Jul99: Interpolation may have different nodes for different nj.
C              This nj loop should include CALC_NNIP (see DIVH1).
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      CALL REFINE_SETNODE(i1,i2,i3,IBT,IDO,IDRN,
     '                  INP,NBJ,ne,ne_new,nj,NKJ,nn_ne,nn_nenew,
     '                  NNIP1,NODE,NPNE,nr,NUMI1,NUNK,
     '                  NVJE,NVJP,SE,SP,XE,XI,XII,XP,EXISTS,NOCROSS,
     '                  ERROR,*9999)
                    ENDDO !njj2
                  ENDDO !njj1

              ENDDO !i1

C*** Fill in the old nodes into the old and new elements and set the
C*** local node properties for the exisiting nodes (versions and
C*** derivatives).
              IF(NUM_COLLAPSED.EQ.2) THEN
                ii2=1
                ii3=1
              ELSE
                IF(COLLAPSED_XI.EQ.IDR2) THEN
                  ii2=1
                  ii3=i3
                ELSE
                  ii3=1
                  ii2=i2
                ENDIF
              ENDIF
              IF(IBT(1,COLLAPSED_XI,nb).EQ.5) THEN
                IF(IBT(1,IDRN,nb).EQ.1.AND.
     '            IBT(2,IDRN,nb).EQ.2) THEN !Quadratic Lagrange
                  IF(XII.EQ.0.5d0) THEN
                    NPNE1(NNIP1(3,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE2(NNIP2(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,i2,i3),nb,nj)
                        NVJE(NNIP2(3,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE IF(XII.GT.0.5d0) THEN
                    NPNE1(NNIP1(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(3,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE
                    NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE2(NNIP2(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(2,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,i2,i3),nb,nj)
                        NVJE(NNIP2(3,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDIF
                ELSE IF(IBT(1,IDRN,nb).EQ.1.AND.
     '            IBT(2,IDRN,nb).EQ.3) THEN !Cubic Lagrange
                  IF(XII.EQ.(1.0d0/3.0d0)) THEN
                    NPNE1(NNIP1(4,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE2(NNIP2(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    NPNE2(NNIP2(4,i2,i3))=NPNE(NNIP1(4,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,i2,i3),nb,nj)
                        NVJE(NNIP2(3,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,i2,i3),nb,nj)
                        NVJE(NNIP2(4,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                    NPNE1(NNIP1(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE1(NNIP1(4,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    NPNE2(NNIP2(4,i2,i3))=NPNE(NNIP1(4,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,i2,i3),nb,nj)
                        NVJE(NNIP2(4,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                    NPNE1(NNIP1(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE1(NNIP1(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(4,i2,i3))=NPNE(NNIP1(4,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(4,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                    NPNE1(NNIP1(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    NPNE2(NNIP2(4,i2,i3))=NPNE(NNIP1(4,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(3,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,i2,i3),nb,nj)
                        NVJE(NNIP2(4,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE
                    NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE2(NNIP2(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    NPNE2(NNIP2(4,i2,i3))=NPNE(NNIP1(4,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(2,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,i2,i3),nb,nj)
                        NVJE(NNIP2(3,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,i2,i3),nb,nj)
                        NVJE(NNIP2(4,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDIF
                ELSE
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE2(NNIP2(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      NVJE(NNIP2(2,i2,i3),nbb,nj,ne_NEW)=
     '                  NVJE1(NNIP1(2,i2,i3),nb,nj)
                      DO nk=1,NKT(NNIP1(2,i2,i3),nb)
                        NKJE(nk,NNIP2(2,i2,i3),nj,ne_NEW)=
     '                    NKJE(nk,NNIP1(2,i2,i3),nj,ne)
                      ENDDO !nk
                    ENDDO !njj2
                  ENDDO !njj1
C                  DO nk=1,NKT(NNIP1(2,i2,i3),nb)
C                    NKE(nk,NNIP2(2,i2,i3),nbb,ne_NEW)=
C     '                NKE(nk,NNIP1(2,i2,i3),nb,ne)
C                  ENDDO !nk
                ENDIF
              ELSE !IBT(1,COLLAPSED_XI,nb)=6
                IF(IBT(1,IDRN,nb).EQ.1.AND.
     '            IBT(2,IDRN,nb).EQ.2) THEN !Quadratic Lagrange
                  IF(XII.EQ.0.5d0) THEN
                    NPNE1(NNIP1(1,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE1(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(1,i2,i3),nb,ne)
                    NPNE2(NNIP2(3,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(1,i2,i3),nb,nj)
                        NVJE(NNIP2(3,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE IF(XII.GT.0.5d0) THEN
                    NPNE1(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(1,i2,i3),nb,ne)
                    NPNE2(NNIP2(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(1,i2,i3),nb,nj)
                        NVJE(NNIP2(2,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE
                    NPNE1(NNIP1(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE1(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(1,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(1,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDIF
                ELSE IF(IBT(1,IDRN,nb).EQ.1.AND.
     '            IBT(2,IDRN,nb).EQ.3) THEN !Cubic Lagrange
                  IF(XII.EQ.(1.0d0/3.0d0)) THEN
                    NPNE1(NNIP1(1,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE1(NNIP1(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    NPNE1(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(1,i2,i3),nb,ne)
                    NPNE2(NNIP2(4,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(1,i2,i3),nb,nj)
                        NVJE(NNIP2(4,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                    NPNE1(NNIP1(1,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    NPNE1(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(1,i2,i3),nb,ne)
                    NPNE2(NNIP2(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE2(NNIP2(4,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(1,i2,i3),nb,nj)
                        NVJE(NNIP2(2,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,i2,i3),nb,nj)
                        NVJE(NNIP2(4,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                    NPNE1(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                    NPNE2(NNIP1(1,i2,i3))=NPNE(NNIP1(1,i2,i3),nb,ne)
                    NPNE2(NNIP2(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE2(NNIP2(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(1,i2,i3),nb,nj)
                        NVJE(NNIP2(2,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,i2,i3),nb,nj)
                        NVJE(NNIP2(3,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                    NPNE1(NNIP1(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    NPNE1(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(1,i2,i3),nb,ne)
                    NPNE2(NNIP2(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(1,i2,i3),nb,nj)
                        NVJE(NNIP2(2,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ELSE
                    NPNE1(NNIP1(2,i2,i3))=NPNE(NNIP1(2,i2,i3),nb,ne)
                    NPNE1(NNIP1(3,i2,i3))=NPNE(NNIP1(3,i2,i3),nb,ne)
                    NPNE1(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                    NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(1,i2,i3),nb,ne)
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                    NVJE1(NNIP1(1,i2,i3),nb,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDIF
                ELSE
                  NPNE1(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP2(1,i2,i3))=NPNE(NNIP1(1,i2,i3),nb,ne)
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      NVJE(NNIP2(1,i2,i3),nbb,nj,ne_NEW)=
     '                  NVJE1(NNIP1(1,i2,i3),nb,nj)
                      DO nk=1,NKT(NNIP1(1,i2,i3),nb)
                        NKJE(nk,NNIP2(1,i2,i3),nj,ne_NEW)=
     '                    NKJE(nk,NNIP1(1,i2,i3),nj,ne)
                      ENDDO !nk
                    ENDDO !njj2
                  ENDDO !njj1
C                  DO nk=1,NKT(NNIP1(1,i2,i3),nb)
C                    NKE(nk,NNIP2(1,i2,i3),nbb,ne_NEW)=
C     '                NKE(nk,NNIP1(1,i2,i3),nb,ne)
C                  ENDDO !nk
                ENDIF
              ENDIF
            ENDDO !i2
          ENDDO !i3

C*** Update scaling factors in current element & new element
          CALL REFINE_SETSCALEFACTORS(IBT,IDO,IDRN,ne,ne_new,NITB,SE,
     '      XII,ERROR,*9999)

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' NPNE1(nn):'',10(1X,I5),'
     '        //'/:(11X,10(1X,I5)))') (NPNE1(nn),nn=1,NNT(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NPNE2(nn):'',10(1X,I5),'
     '        //'/:(11X,10(1X,I5)))') (NPNE2(nn),nn=1,NNT(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

          DO nb1=1,NBFT
            DO nn=1,NNT(nb1)
              NPNE(nn,nb1,ne)=NPNE1(nn)
              NPNE(nn,nb1,ne_NEW)=NPNE2(nn)
            ENDDO !nn
          ENDDO !nb1

C*** Form other arrays for new element ne_NEW
          NHE(ne_NEW)=NHE(ne)
          NW(ne_NEW,1)=NW(ne,1)
          NRE(ne_NEW)=NRE(ne)
          DO nm=1,NMM
            CE(nm,ne_NEW)=CE(nm,ne)
          ENDDO !nm

          IF(CALL_EQUA.OR.CALL_FIT.OR.CALL_OPTI) THEN
            DO nhx=1,NHE(ne)
              nh=NH_LOC(nhx,nx)
              DO nc=1,NCM
                NBH(nh,nc,ne_NEW)=nbb
              ENDDO !nc
            ENDDO !nhx
          ENDIF !call_equa etc

          NPSTART=NPT(0)+1
          NET(0)=ne_NEW !is new highest element # in mesh
          NET(nr)=ne_NEW !is new highest element # in region nr
          CALL ASSERT(NET(nr).LE.NEM,
     '      '>>Too many elements - increase NEM',ERROR,*9999)
          NEELEM(0,nr)=NEELEM(0,nr)+1 !is new #elements in region nr
          NEELEM(NEELEM(0,nr),nr)=ne_NEW !is new element number

C*** Update element list in each group
          CALL UPGREL(ne,ne_NEW,ERROR,*9999)

C*** Check that the global nodes of new element ne_NEW are in the
C*** node list for the current region
          DO nn=1,NNT(nbb)
            IF(NPNODE(0,nr).LE.NP_R_M) THEN
              np=NPNE(nn,nbb,ne_NEW)
              FOUND=.FALSE.
              nonode=1
              DO WHILE(.NOT.FOUND.AND.nonode.LE.NPNODE(0,nr))
                IF(NPNODE(nonode,nr).eq.np) THEN
                  FOUND=.TRUE.
                ELSE
                  nonode=nonode+1
                ENDIF
              ENDDO
              IF(.NOT.FOUND) THEN
                NPNODE(0,nr)=NPNODE(0,nr)+1
                IF(NPNODE(0,nr).LE.NP_R_M) NPNODE(NPNODE(0,nr),nr)=np
              ENDIF
            ENDIF
          ENDDO !nn
          CALL ASSERT(NPNODE(0,nr).LE.NP_R_M,'>>Increase NP_R_M',
     '      ERROR,*9999)

        ELSE !Case=3

          CALL CALC_NNIP(COLLAPSED_XI,IBT,IDRN,IDR2,IDR3,IMAP,INP,
     '      nb,nb,NNIP1,NNIP2,NUM_COLLAPSED,NUMI1,NUMI2,NUMI3,PERP_XI,
     '      ATCOLLAPSE,.FALSE.,.TRUE.,ERROR,*9999)

          DO njj1=1,3
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              NBJ(nj,ne_NEW)=NBJ(nj,ne)
            ENDDO !njj2
          ENDDO !njj1

C*** Loop over Xi directions
          DO i3=1,NUMI3
            IF(NUMI3.EQ.1) THEN
              XI(3)=0.0d0
            ELSE
              XI(IDR3)=DBLE(i3-1)/DBLE(NUMI3-1)
            ENDIF
            DO i2=1,NUMI2
              XI(IDR2)=DBLE(i2-1)/DBLE(NUMI2-1)
              DO i1=1,NUMI1-1
                ii1=IMAP(1,i1,i2,i3)
                ii2=IMAP(2,i1,i2,i3)
                ii3=IMAP(3,i1,i2,i3)

C*** Get the xi position
                XI(IDRN)=REFINE_GETXIPOS(ii1,IBT,IDRN,nb,XII)

C*** Get the node number of the refine point
                CALL REFINE_FINDNODE(IBT,IDO,INP,nb,NBJ,ne,NKJ,NODE,
     '            NPLIST,NPNODE,NPTNEW,NPTOLD,nr,NUNK,NVJP,XE,XI,XP,
     '            EXISTS,ERROR,*9999)

C*** Calculate the element arrays for the node
                IF((IBT(1,IDRN,nb).EQ.1.OR.IBT(1,IDRN,nb).EQ.5.OR.
     '            IBT(1,IDRN,nb).EQ.6).AND.IBT(2,IDRN,nb).EQ.2)
     '            THEN !Quadratic Lagrange
                  IF(XII.EQ.0.5d0) THEN
                    IF(ii1.EQ.1) THEN
                      nn_ne=NNIP1(2,ii2,ii3)
                      nn_nenew=0
                    ELSE
                      nn_ne=0
                      nn_nenew=NNIP1(2,ii2,ii3)
                    ENDIF
                  ELSE IF(XII.GT.0.5d0) THEN
                    IF(ii1.EQ.1) THEN
                      nn_ne=NNIP1(3,ii2,ii3)
                    ELSE
                      nn_ne=0
                    ENDIF
                    nn_nenew=NNIP1(ii1,ii2,ii3)
                  ELSE
                    nn_ne=NNIP1(ii1+1,ii2,ii3)
                    IF(ii1.EQ.2) THEN
                      nn_nenew=NNIP1(1,ii2,ii3)
                    ELSE
                      nn_nenew=0
                    ENDIF
                  ENDIF
                ELSE IF((IBT(1,IDRN,nb).EQ.1.OR.
     '              IBT(1,IDRN,nb).EQ.5.OR.IBT(1,IDRN,nb).EQ.6).AND.
     '              IBT(2,IDRN,nb).EQ.3) THEN !Cubic Lagrange
                  IF(XII.EQ.(1.0d0/3.0d0)) THEN
                    IF(ii1.EQ.3) THEN
                      nn_ne=0
                      nn_nenew=NNIP1(2,ii2,ii3)
                    ELSE
                      nn_ne=NNIP1(ii1+1,ii2,ii3)
                      nn_nenew=0
                    ENDIF
                  ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                    IF(ii1.EQ.1) THEN
                      nn_ne=NNIP1(3,ii2,ii3)
                      nn_nenew=0
                    ELSE
                      nn_ne=0
                      nn_nenew=NNIP1(ii1,ii2,ii3)
                    ENDIF
                  ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                    IF(ii1.EQ.1) THEN
                      nn_ne=NNIP1(4,ii2,ii3)
                    ELSE
                      nn_ne=0
                    ENDIF
                    nn_nenew=NNIP1(ii1,ii2,ii3)
                  ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                    IF(ii1.LE.2) THEN
                      nn_ne=NNIP1(ii1+2,ii2,ii3)
                    ELSE
                      nn_ne=0
                    ENDIF
                    IF(ii1.GE.2) THEN
                      nn_nenew=NNIP1(ii1-1,ii2,ii3)
                    ELSE
                      nn_nenew=0
                    ENDIF
                  ELSE
                    nn_ne=NNIP1(ii1+1,ii2,ii3)
                    IF(ii1.EQ.3) THEN
                      nn_nenew=NNIP1(1,ii2,ii3)
                    ELSE
                      nn_nenew=0
                    ENDIF
                  ENDIF
                ELSE
                  nn_ne=NNIP1(2,ii2,ii3)
                  nn_nenew=NNIP1(1,ii2,ii3)
                ENDIF

C*** Put the node in the element
                IF(nn_ne.NE.0) NPNE1(nn_ne)=NODE
                IF(nn_nenew.NE.0) NPNE2(nn_nenew)=NODE

C*** Set the node properties
C KAT 22Jul99: Interpolation may have different nodes for different nj.
C              This nj loop should include CALC_NNIP (see DIVH1).
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      CALL REFINE_SETNODE(ii1,ii2,ii3,IBT,IDO,IDRN,
     '                  INP,NBJ,ne,ne_new,nj,NKJ,nn_ne,nn_nenew,
     '                  NNIP1,NODE,NPNE,nr,NUMI1,
     '                  NUNK,NVJE,NVJP,SE,SP,
     '                  XE,XI,XII,XP,EXISTS,NOCROSS,
     '                  ERROR,*9999)
                    ENDDO !njj2
                  ENDDO !njj1

              ENDDO !i1

C*** Fill in the old nodes into the old and new elements and set the
C*** local node properties for the exisiting nodes (versions and
C*** derivatives).
              IF((IBT(1,IDRN,nb).EQ.1.OR.IBT(1,IDRN,nb).EQ.5.OR.
     '          IBT(1,IDRN,nb).EQ.6).AND.IBT(2,IDRN,nb).EQ.2)
     '          THEN !Quadratic Lagrange
                IF(XII.EQ.0.5d0) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(3,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(1,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(1,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE IF(XII.GT.0.5d0) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(2,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ENDIF
              ELSE IF((IBT(1,IDRN,nb).EQ.1.OR.
     '            IBT(1,IDRN,nb).EQ.5.OR.IBT(1,IDRN,nb).EQ.6).AND.
     '            IBT(2,IDRN,nb).EQ.3) THEN !Cubic Lagrange
                IF(XII.EQ.(1.0d0/3.0d0)) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(4,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(1,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(1,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(4,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(4,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(1,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(1,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(4,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(4,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE1(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(4,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ELSE
                  NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(3,ii2,ii3))=NPNE(NNIP1(3,ii2,ii3),nb,ne)
                  NPNE2(NNIP1(4,ii2,ii3))=NPNE(NNIP1(4,ii2,ii3),nb,ne)
                  DO nb1=1,NBFT
                    DO njj1=1,3
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NVJE(NNIP1(2,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(2,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(3,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(3,ii2,ii3),nb1,nj)
                        NVJE(NNIP1(4,ii2,ii3),nb1,nj,ne_NEW)=
     '                    NVJE1(NNIP1(4,ii2,ii3),nb1,nj)
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDDO !nb1
                ENDIF
              ELSE
                NPNE1(NNIP1(1,ii2,ii3))=NPNE(NNIP1(1,ii2,ii3),nb,ne)
                NPNE2(NNIP1(2,ii2,ii3))=NPNE(NNIP1(2,ii2,ii3),nb,ne)
                DO nb1=1,NBFT
                  DO njj1=1,3
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      NVJE(NNIP1(2,ii2,ii3),nb1,nj,ne_NEW)=
     '                  NVJE1(NNIP1(2,ii2,ii3),nb1,nj)
                    ENDDO !njj2
                  ENDDO !njj1
                ENDDO !nb1
              ENDIF
            ENDDO !i2
          ENDDO !i3

C*** Update scaling factors in current element & new element
          CALL REFINE_SETSCALEFACTORS(IBT,IDO,IDRN,ne,ne_new,NITB,SE,
     '      XII,ERROR,*9999)

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' NPNE1(nn):'',10(1X,I5),'
     '        //'/:(11X,10(1X,I5)))') (NPNE1(nn),nn=1,NNT(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NPNE2(nn):'',10(1X,I5),'
     '        //'/:(11X,10(1X,I5)))') (NPNE2(nn),nn=1,NNT(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

          DO nb1=1,NBFT
            DO nn=1,NNT(nb1)
              NPNE(nn,nb1,ne)=NPNE1(nn)
              NPNE(nn,nb1,ne_NEW)=NPNE2(nn)
C              DO nk=1,NKT(nn,nb1)
C                NKE(nk,nn,nb1,ne_NEW)=NKE(nk,nn,nb1,ne)
C              ENDDO !nk
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

C*** Form other arrays for new element ne_NEW
          NHE(ne_NEW)=NHE(ne)
          NW(ne_NEW,1)=NW(ne,1)
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

          NPSTART=NPT(0)+1
          NET(0)=ne_NEW !is new highest element # in mesh
          NET(nr)=ne_NEW !is new highest element # in region nr
          CALL ASSERT(NET(nr).LE.NEM,
     '      '>>Too many elements - increase NEM',ERROR,*9999)
          NEELEM(0,nr)=NEELEM(0,nr)+1 !is new #elements in region nr
          NEELEM(NEELEM(0,nr),nr)=ne_NEW !is new element number

C*** Update element list in each group
          CALL UPGREL(ne,ne_NEW,ERROR,*9999)

C*** Check that the global nodes of new element ne_NEW are in the
C*** node list for the current region
          DO nn=1,NNT(nb)
            IF(NPNODE(0,nr).LE.NP_R_M) THEN
              np=NPNE(nn,nb,ne_NEW)
              FOUND=.FALSE.
              nonode=1
              DO WHILE(.NOT.FOUND.AND.nonode.LE.NPNODE(0,nr))
                IF(NPNODE(nonode,nr).eq.np) THEN
                  FOUND=.TRUE.
                ELSE
                  nonode=nonode+1
                ENDIF
              ENDDO
              IF(.NOT.FOUND) THEN
                NPNODE(0,nr)=NPNODE(0,nr)+1
                IF(NPNODE(0,nr).LE.NP_R_M) NPNODE(NPNODE(0,nr),nr)=np
              ENDIF
            ENDIF
          ENDDO !nn
          CALL ASSERT(NPNODE(0,nr).LE.NP_R_M,'>>Increase NP_R_M',
     '      ERROR,*9999)
        ENDIF
      ENDIF
      NET(0)=ne_NEW   !is new highest element # in mesh

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

      CALL EXITS('DIVS2')
      RETURN
 9999 CALL ERRORS('DIVS2',ERROR)
      CALL EXITS('DIVS2')
      RETURN 1
      END


