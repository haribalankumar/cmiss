      SUBROUTINE DIVS1(IBT,IDO,IDRN,INP,NBH,NBJ,ne,
     '  NEELEM,NHE,NKJE,NKJ,NPF,NPNE,NPNODE,
     '  NPSTART,nr,NRE,NVJE,NVJP,NW,nx,
     '  CE,SE,XA,XE,XII,XP,NOCROSS,ERROR,*)

C#### Subroutine: DIVS1
C###  Description:
C###    DIVS1 subdivides special Hermite simplex element ne in the
C###    Xi(IDRN) direction at Xi(IDRN)=XII. Currently only for 2d
C###    elements.  If the element is refined in the xi1 direction,
C###    then a new Hermite simplex elements are obtained. If xi2 is
C###    the direction of refinement, then a bicubic hermite element is
C###    created.
C**** NPSTART is first node number for new nodes
C**** DSDXI(ni) are derivatives of arc-length wrt Xi

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),IDRN,
     '  INP(NNM,NIM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),ne,
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NPSTART,
     '  nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJP(NJM,NPM),NW(NEM,3),nx
      REAL*8 CE(NMM,NEM),SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XII,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL NOCROSS
!     Local Variables
      INTEGER i3,IK,k,nb,NB_HERMITE,nb1,nc,ne_NEW,nh,nhx,
     '  ni,nj,nk,nm,nn,NNTB,NODE,NODE1,noelem,nonode,np,
     '  NPNE1(27),NPNE2(27),NPTNEW,NPTOLD,NP_TEST,nr2,ns,nu
      REAL*8 AA,DIFF,DSDXI(3),G1,G3,PXI,R,RC,RR,RRC,SLX,SMX,
     '  SUM,X(10,6),XI(3),XS(4)
      LOGICAL EXIST,FOUND

      CALL ENTERS('DIVS1',*9999)

      nc=1 ! temporary cpb 23/11/94
      NPTNEW=NPSTART-1
      ne_NEW=NET(0)+1
      CALL ASSERT(ne_NEW.LE.NEM,'>>Too many elements -Increase NEM',
     '  ERROR,*9999)
      nb1=NBJ(1,ne)
      NNTB=NNT(nb1)
      IF(IDRN.EQ.1) THEN !Refining in xi1 direction
        XI(1)=XII
        IF(NKT(1,nb1).EQ.1)THEN !Apex at node 1
          XI(2)=1.0D0
        ELSE IF(NKT(3,nb1).EQ.1)THEN !Apex at node 3
          XI(2)=0.0D0
        ENDIF

        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)

        DO nj=1,NJT !to find coords of proposed new node position
          nb=NBJ(nj,ne)
          XS(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '               nb,1,XI,XE(1,nj))
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,FMT='('' XI: '',3E12.3)') (XI(ni),ni=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,FMT='('' XS: '',3E12.3)') (XS(nj),nj=1,NJT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        DO np=1,NPTNEW
C ***     Search for existing nodes
          IK=0
          DO nj=1,NJT
            DIFF= DABS(XS(nj)-XP(1,1,nj,np))
            IF((DIFF.LT.(1.d-3)* DABS(XS(nj)).OR.DIFF.LT.1.d-6).
     '        OR.(ITYP10(1).eq.2.and.nj.EQ.2.AND.( DABS(XP(1,1,1,np)).
     '        LT.1.d-6.OR. DABS(DIFF-2.0D0*PI).LT.1.d-6)).
     '        OR.(ITYP10(1).eq.3.and.nj.EQ.2.AND.( DABS(XP(1,1,1,np)).
     '        LT.1.d-6.OR. DABS(DIFF-2.0D0*PI).LT.1.d-6)).
     '        OR.(ITYP10(1).eq.3.and.nj.EQ.3.AND.( DABS(XP(1,1,1,np)).
     '        LT.1.d-6.OR. DABS(DIFF-2.0D0*PI).LT.1.d-6)).
     '        OR.(ITYP10(1).eq.4.and.nj.EQ.3.AND.( DABS(XP(1,1,2,np)).
     '        LT.1.d-6.OR. DABS(XP(1,1,2,np)-PI).LT.1.d-6.
     '        OR. DABS(DIFF-2.d0*PI).LT.1.d-5))) THEN
                IK=IK+1
            ENDIF
          ENDDO
          IF(IK.EQ.NJT) THEN
            EXIST=.TRUE.
            NODE=np
            GO TO 41
          ELSE
            EXIST=.FALSE.
          ENDIF
        ENDDO

 41     IF(.NOT.EXIST) THEN
          NPTNEW=NPTNEW+1
          CALL ASSERT(NPTNEW.LE.NPM,
     '      '>>Too many nodes - increase NPM',ERROR,*9999)
          NPNODE(0,nr)=NPNODE(0,nr)+1
          NPNODE(NPNODE(0,nr),nr)=NPTNEW
          NODE=NPTNEW
          NPT(nr)=NPTNEW
          NPT( 0)=NPTNEW

C ***     Calculate nodal coords & derivs wrt Xi-directions
          DO nj=1,NJ_LOC(0,0,nr)
            nb=NBJ(nj,ne)
!AJP 12/7/96 Check on nb added because if a fibre field has been
!defined in another region nj_loc(0,0) may exceed what is needed for
!the current region and nb may not be defined.
            IF(nb.GT.0) THEN
              DO nu=1,NUT(nb)
                X(nu,nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,nu,XI,XE(1,nj))
              ENDDO
            ELSE
              X(1,nj)=0.0d0 !since this is used later
            ENDIF !nb
          ENDDO

C ***     If an angle is 2*pi set it to zero
          IF(ITYP10(1).GE.2) THEN
            IF(DABS(X(1,2)-2.0D0*PI).LT.1.d-5) X(1,2)=0.0D0
            IF(ITYP10(1).GE.3) THEN
              IF(DABS(X(1,3)-2.0D0*PI).LT.1.d-5) X(1,3)=0.0D0
            ENDIF
          ENDIF

          IF(DOP) THEN
            DO nj=1,NJT
              WRITE(OP_STRING,'('' EXIST='',L1,'' NODE='',I4,'
     '          //''' X(nk,'',i1,''): '',7E12.3)')
     '          EXIST,NODE,nj,(X(nk,nj),nk=1,NKT(0,NBJ(nj,ne)))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF

          NPTOLD=NPNODE(NPNODE(0,nr)-1,nr)
          DO nj=1,NJ_LOC(0,0,nr)
            NKJ(nj,NPTNEW)=NKJ(nj,NPTOLD)
            NVJP(nj,NPTNEW)=NVJP(nj,NPTOLD)
          ENDDO

          IF(NBI(nb1).EQ.3) THEN !global line derivatives
            WRITE(OP_STRING,'('' Not implemented in DIVS1'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C cpb 5/12/96 re-adding arc-length scaling
c cpb 8/10/94 changing around nbi(nb) = 4 and 5
C          ELSE IF(NBI(nb1).EQ.4) THEN !arc-length derivatives
          ELSE IF(NBI(nb1).GE.5.AND.NBI(nb1).LE.7)
     '        THEN !arc-length or ave. arc-length derivatives
C ***       Calculate derivatives of arc-length wrt Xi
            DO ni=1,NIT(nb1)
              nu=1+ni*(1+ni)/2
              SUM=0.0D0
              IF(ITYP10(1).EQ.1) THEN
                DO nj=1,NJT
                  SUM=SUM+X(nu,nj)**2
                ENDDO
              ELSE IF(ITYP10(1).EQ.2) THEN
                R=X(1,1)
                RR=R*R
                SUM=SUM+X(nu,1)**2+RR*X(nu,2)**2
                IF(NJT.GT.2) SUM=SUM+X(nu,3)**2
              ELSE IF(ITYP10(1).EQ.3) THEN
                R=X(1,1)
                RR=R*R
                RC=R*DCOS(X(1,3))
                RRC=RC*RC
                SUM=SUM+X(nu,1)**2+RRC*X(nu,2)**2+RR*X(nu,3)**2
              ELSE IF(ITYP10(1).EQ.4) THEN
                AA=FOCUS*FOCUS
                SLX=DSINH(X(1,1))
                SMX=DSIN(X(1,2))
                G1=AA*(SLX*SLX+SMX*SMX)
                G3=AA* SLX*SLX*SMX*SMX
                SUM=SUM+G1*(X(nu,1)**2+X(nu,2)**2)
                IF(NJT.GT.2) SUM=SUM+G3*X(nu,3)**2
              ENDIF
              DSDXI(ni)=DSQRT(SUM)
            ENDDO
            IF(DOP) THEN
              WRITE(OP_STRING,'('' DSDXI(ni): '',3E12.3)')
     '          (DSDXI(ni),ni=1,NIT(nb1))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF

C cpb 5/12/96 re-adding arc-length scaling
c cpb 8/10/94 changing around nbi(nb) = 4 and 5
C          IF(NBI(nb1).EQ.3.OR.NBI(nb1).EQ.4) THEN
          IF(NBI(nb1).EQ.3.OR.(NBI(nb1).GE.5.AND.NBI(nb1).LE.7)) THEN
C ***       Store nodal derivatives wrt arc-length
            DO nj=1,NJ_LOC(0,0,nr)
              XP(1,1,nj,NODE)=X(1,nj)
              nb=NBJ(nj,ne)
!AJP 12/7/96 Check on nb added becuase if a fibre field has been
!defined in another region nj_loc(0,0,nr) may exceed what is needed for
!the current region and nb may not be defined.
              IF(nb.GT.0) THEN
                DO nk=2,NKT(0,nb)
                  IF(nk.EQ.2) THEN
                    IF(DABS(DSDXI(1)).GT.1.d-6) THEN
                      XP(2,1,nj,NODE)=X(2,nj)/DSDXI(1)
                    ELSE
                      WRITE(OP_STRING,
     '                  '('' Warning: Arc length deriv is zero'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(nk.EQ.3) THEN
                    IF(DABS(DSDXI(2)).GT.1.d-6) THEN
                      XP(3,1,nj,NODE)=X(4,nj)/DSDXI(2)
                    ELSE
                      WRITE(OP_STRING,
     '                  '('' Warning: Arc length deriv is zero'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(nk.EQ.4) THEN
                    IF(NOCROSS) THEN
                      XP(4,1,nj,NODE)=0.0D0
                    ELSE
                      IF(DABS(DSDXI(1)).GT.1.d-6.AND.DABS(DSDXI(2))
     '                   .GT.1.d-6) THEN
                        XP(4,1,nj,NODE)=X(6,nj)/(DSDXI(1)*DSDXI(2))
                      ELSE
                        WRITE(OP_STRING,
     '                     '('' Warning: Arc length deriv is zero'')')
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF !nb
            ENDDO
          ELSE
            DO nj=1,NJ_LOC(0,0,nr)
              XP(1,1,nj,NODE)=X(1,nj)
            ENDDO
          ENDIF

        ELSE !Node exits, but possibly only in another region
          FOUND=.FALSE.
          nonode=1
          DO WHILE((nonode.LE.NPNODE(0,nr)).AND.(.NOT.FOUND))
            NP_TEST=NPNODE(nonode,nr)
            IF(np_test.eq.np)THEN !Node np is in current region
              FOUND=.TRUE.
            ELSE
              nonode=nonode+1
            ENDIF
          ENDDO
          IF(.NOT.FOUND)THEN !Update arrays to include np
            NPNODE(0,nr)=NPNODE(0,nr)+1
            NPNODE(NPNODE(0,nr),nr)=np
            IF(np.GT.NPT(nr))THEN
              NPT(nr)=np
            ENDIF
! New AJP 21-1-94
            !Node exists but in another region.
            !Update current region arrays.
            !Firstly find other region number of node np
            FOUND=.FALSE.
            nr2=1
            DO WHILE ((nr2.LE.NRT).AND.(.NOT.FOUND))
              IF(nr2.ne.nr)THEN
                nonode=1
                DO WHILE ((nonode.LE.NPNODE(0,nr2)).AND.
     '                    (.NOT.FOUND))
                  NP_TEST=NPNODE(nonode,nr2)
                  IF(np.EQ.NP_TEST)THEN
                    FOUND=.TRUE.
                  ELSE
                    nonode=nonode+1
                  ENDIF
                ENDDO
                IF(.NOT.FOUND)nr2=nr2+1
              ELSE
                nr2=nr2+1
              ENDIF
            ENDDO
            IF(nr2.GT.NRT)THEN
              WRITE(OP_STRING,
     '          '('' Cant find node in other regions ?'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              GOTO 9999
            ENDIF
!End New
          ENDIF
        ENDIF !.not.exist
        IF(NKT(1,nb1).EQ.1)THEN !Apex at node 1
          NPNE1(1)=NPNE(1,nb1,ne)
          NPNE1(2)=NPNE(2,nb1,ne)
          NPNE1(3)=NODE
          NPNE2(1)=NPNE(1,nb1,ne)
          NPNE2(2)=NODE
          NPNE2(3)=NPNE(3,nb1,ne)
        ELSE IF(NKT(3,nb1).EQ.1)THEN !Apex at node 3
          NPNE1(1)=NPNE(1,nb1,ne)
          NPNE1(2)=NODE
          NPNE1(3)=NPNE(3,nb1,ne)
          NPNE2(1)=NODE
          NPNE2(2)=NPNE(2,nb1,ne)
          NPNE2(3)=NPNE(3,nb1,ne)
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'('' NPNE1: '',27I3)') (NPNE1(nn),nn=1,NNTB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NPNE2: '',27I3)') (NPNE2(nn),nn=1,NNTB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO nb=1,NBFT
          IF(NNT(nb).LE.3)THEN
!           Same type of basis function, or a basis function with
!           no more nodes than the simplex
            DO nn=1,NNT(nb)
              NPNE(nn,nb,ne)=NPNE1(nn)
              NPNE(nn,nb,ne_NEW)=NPNE2(nn)
            ENDDO
          ELSE
!           Possibly a Lagrange product basis function
            IF(NKT(1,nb1).EQ.1)THEN
              NPNE(1,nb,ne)=NPNE1(1)
              NPNE(1,nb,ne_NEW)=NPNE2(1)
              NPNE(2,nb,ne)=NPNE1(1)
              NPNE(2,nb,ne_NEW)=NPNE2(1)
              NPNE(3,nb,ne)=NPNE1(2)
              NPNE(3,nb,ne_NEW)=NPNE2(2)
              NPNE(4,nb,ne)=NPNE1(3)
              NPNE(4,nb,ne_NEW)=NPNE2(3)
            ELSEIF(NKT(3,nb1).EQ.1)THEN
              NPNE(1,nb,ne)=NPNE1(1)
              NPNE(1,nb,ne_NEW)=NPNE2(1)
              NPNE(2,nb,ne)=NPNE1(2)
              NPNE(2,nb,ne_NEW)=NPNE2(2)
              NPNE(3,nb,ne)=NPNE1(3)
              NPNE(3,nb,ne_NEW)=NPNE2(3)
              NPNE(4,nb,ne)=NPNE1(3)
              NPNE(4,nb,ne_NEW)=NPNE2(3)
            ENDIF
            IF(NNT(nb).GT.4)THEN
              WRITE(OP_STRING,
     '          '('' Warning: Cant update NPNE correctly in DIVS1'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nn=5,NNT(nb)
                NPNE(nn,nb,ne)=NPNE(4,nb,ne)
                NPNE(nn,nb,ne_NEW)=NPNE(4,nb,ne_NEW)
              ENDDO
            ENDIF
          ENDIF
        ENDDO

C ***   Update scaling factors in current element & new element
        IF(NKT(1,nb1).EQ.1)THEN !Apex at node 1
        DO nb=1,NBFT
C ***       Xi(1) derivs
            SE( 3,nb,ne_NEW)=0.25D0*(SE( 3,nb,ne)+SE( 7,nb,ne))
            SE( 7,nb,ne_NEW)=0.50D0* SE( 7,nb,ne)
            SE( 7,nb,ne)    =0.25D0*(SE( 3,nb,ne)+SE( 7,nb,ne))
            SE( 3,nb,ne)    =0.50D0* SE( 3,nb,ne)
C ***       Xi(2) derivs
            SE( 4,nb,ne_NEW)=0.50D0*(SE( 4,nb,ne)+SE( 8,nb,ne))
            SE( 8,nb,ne_NEW)=      SE( 8,nb,ne)
            SE( 8,nb,ne)    =0.50D0*(SE( 4,nb,ne)+SE( 8,nb,ne))
            SE( 4,nb,ne)    =      SE( 4,nb,ne)

C ***       Xi(1),Xi(2) derivs
            SE( 5,nb,ne_NEW)=SE( 3,nb,ne_NEW)*SE( 4,nb,ne_NEW)
            SE( 9,nb,ne_NEW)=SE( 7,nb,ne_NEW)*SE( 8,nb,ne_NEW)
            SE( 5,nb,ne)  =SE( 3,nb,ne)  *SE( 4,nb,ne)
            SE( 9,nb,ne)  =SE( 7,nb,ne)  *SE( 8,nb,ne)
            ns=0
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                ns=ns+1
                SE(ns,nb,ne_NEW)=SE(ns,nb,ne)
              ENDDO
            ENDDO
            IF(DOP) THEN
              DO ns=1,NST(nb)
                WRITE(OP_STRING,'('' ne='',I4,'' nb='',i1,'' ns='','
     '            //'I2,'//' '' SE='',E12.3)') ne,nb,ns,SE(ns,nb,ne)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO
        ELSE IF(NKT(3,nb1).EQ.3)THEN !Apex at node 3
          DO nb=1,NBFT
C ***       Xi(1) derivs
            SE( 2,nb,ne_NEW)=0.25D0*(SE( 2,nb,ne)+SE( 6,nb,ne))
            SE( 6,nb,ne_NEW)=0.50D0* SE( 6,nb,ne)
            SE( 6,nb,ne)    =0.25D0*(SE( 2,nb,ne)+SE( 6,nb,ne))
            SE( 2,nb,ne)    =0.50D0* SE( 2,nb,ne)
C ***       Xi(2) derivs
            SE( 3,nb,ne_NEW)=0.50D0*(SE( 3,nb,ne)+SE( 7,nb,ne))
            SE( 7,nb,ne_NEW)=      SE( 7,nb,ne)
            SE( 7,nb,ne)    =0.50D0*(SE( 3,nb,ne)+SE( 7,nb,ne))
            SE( 3,nb,ne)    =      SE( 3,nb,ne)

C ***       Xi(1),Xi(2) derivs
            SE( 4,nb,ne_NEW)=SE( 2,nb,ne_NEW)*SE( 3,nb,ne_NEW)
            SE( 8,nb,ne_NEW)=SE( 6,nb,ne_NEW)*SE( 7,nb,ne_NEW)
            SE( 4,nb,ne)  =SE( 2,nb,ne)  *SE( 3,nb,ne)
            SE( 8,nb,ne)  =SE( 6,nb,ne)  *SE( 7,nb,ne)
            ns=0
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                ns=ns+1
                SE(ns,nb,ne_NEW)=SE(ns,nb,ne)
              ENDDO
            ENDDO
            IF(DOP) THEN
              DO ns=1,NST(nb)
                WRITE(OP_STRING,'('' ne='',I4,'' nb='',i1,'' ns='','
     '            //'I2,'//' '' SE='',E12.3)') ne,nb,ns,SE(ns,nb,ne)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO
        ENDIF

      ELSE IF(IDRN.EQ.2)THEN !Refining in xi2 direction
!       !Firstly check that a bicubic hermite basis function has been
!       !set up.
        nb=1
        FOUND=.FALSE.
        DO WHILE((.NOT.FOUND).AND.(nb.LE.NBT))
          k=0
          DO ni=1,NIT(nb1)
            IF((IBT(1,ni,nb).EQ.2).AND.(IBT(2,ni,nb).EQ.1))k=k+1
          ENDDO
          IF(k.EQ.NIT(nb1))THEN
            FOUND=.TRUE.
          ELSE
            nb=nb+1
          ENDIF
        ENDDO
        CALL ASSERT(nb.LE.NBT,'Need to set up a bicubic basis function',
     '    ERROR,*9999)
        NB_HERMITE=nb
        DO i3=1,2 !Need to create 2 new nodes
          XI(1)=i3-1
          XI(2)=XII
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)

          DO nj=1,NJT !to find coords of proposed new node position
            nb=NBJ(nj,ne)
            XS(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                 nb,1,XI,XE(1,nj))
          ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,FMT='('' XI: '',3E12.3)') (XI(ni),ni=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,FMT='('' XS: '',3E12.3)') (XS(nj),nj=1,NJT)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          DO np=1,NPTNEW
C ***       Search for existing nodes
            IK=0
            DO nj=1,NJT
              DIFF= DABS(XS(nj)-XP(1,1,nj,np))
              IF((DIFF.LT.(1.d-3)* DABS(XS(nj)).OR.DIFF.LT.1.d-6).
     '          OR.(ITYP10(1).eq.2.and.nj.EQ.2.AND.(DABS(XP(1,1,1,np)).
     '          LT.1.d-6.OR. DABS(DIFF-2.0D0*PI).LT.1.d-6)).
     '          OR.(ITYP10(1).eq.3.and.nj.EQ.2.AND.(DABS(XP(1,1,1,np)).
     '          LT.1.d-6.OR. DABS(DIFF-2.0D0*PI).LT.1.d-6)).
     '          OR.(ITYP10(1).eq.3.and.nj.EQ.3.AND.(DABS(XP(1,1,1,np)).
     '          LT.1.d-6.OR. DABS(DIFF-2.0D0*PI).LT.1.d-6)).
     '          OR.(ITYP10(1).eq.4.and.nj.EQ.3.AND.(DABS(XP(1,1,2,np)).
     '          LT.1.d-6.OR. DABS(XP(1,1,2,np)-PI).LT.1.d-6.
     '          OR. DABS(DIFF-2.0D0*PI).LT.1.d-5))) THEN
                  IK=IK+1
              ENDIF
            ENDDO
            IF(IK.EQ.NJT) THEN
              EXIST=.TRUE.
              NODE=np
              GO TO 42
            ELSE
              EXIST=.FALSE.
            ENDIF
          ENDDO

 42       IF(.NOT.EXIST) THEN
            NPTNEW=NPTNEW+1
            CALL ASSERT(NPTNEW.LE.NPM,
     '        '>>Too many nodes - increase NPM',ERROR,*9999)
            NPNODE(0,nr)=NPNODE(0,nr)+1
            NPNODE(NPNODE(0,nr),nr)=NPTNEW
            NODE=NPTNEW
            NPT(nr)=NPTNEW
            NPT( 0)=NPTNEW

C ***       Calculate nodal coords & derivs wrt Xi-directions
            DO nj=1,NJ_LOC(0,0,nr)
              nb=NBJ(nj,ne)
!AJP 12/7/96 Check on nb added becuase if a fibre field has been
!defined in another region nj_loc(0,0,nr) may exceed what is needed for
!the current region and nb may not be defined.
              IF(nb.GT.0) THEN
                DO nu=1,NUT(nb)
                  X(nu,nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,nu,XI,XE(1,nj))
                ENDDO
              ELSE
                X(1,nj)=0.0d0 !since this is used later
              ENDIF !nb
            ENDDO

C ***       If an angle is 2*pi set it to zero
            IF(ITYP10(1).GE.2) THEN
              IF(DABS(X(1,2)-2.0D0*PI).LT.1.d-5) X(1,2)=0.0D0
              IF(ITYP10(1).GE.3) THEN
                IF(DABS(X(1,3)-2.0D0*PI).LT.1.d-5) X(1,3)=0.0D0
              ENDIF
            ENDIF

            IF(DOP) THEN
              DO nj=1,NJT
                WRITE(OP_STRING,'('' EXIST='',L1,'' NODE='',I4,'
     '            //''' X(nk,'',i1,''): '',7E12.3)')
     '            EXIST,NODE,nj,(X(nk,nj),nk=1,NKT(0,NBJ(nj,ne)))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF

            NPTOLD=NPNODE(NPNODE(0,nr)-1,nr)
            DO nj=1,NJ_LOC(0,0,nr)
              NKJ(nj,NPTNEW)=NKJ(nj,NPTOLD)
              NVJP(nj,NPTNEW)=NVJP(nj,NPTOLD)
            ENDDO

            IF(NBI(nb1).EQ.3) THEN !global line derivatives
              WRITE(OP_STRING,'('' Not implemented in DIVS1'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C cpb 5/12/96 re-adding arc-length scaling
c cpb 8/10/94 changing around nbi(nb) = 4 and 5
C            ELSE IF(NBI(nb1).EQ.4) THEN !arc-length derivatives
            ELSE IF(NBI(nb1).GE.5.AND.NBI(nb1).LE.7)
     '          THEN !arc-length or ave. arc-length derivatives
C ***         Calculate derivatives of arc-length wrt Xi
              DO ni=1,NIT(nb1)
                nu=1+ni*(1+ni)/2
                SUM=0.0D0
                IF(ITYP10(1).EQ.1) THEN
                  DO nj=1,NJT
                    SUM=SUM+X(nu,nj)**2
                  ENDDO
                ELSE IF(ITYP10(1).EQ.2) THEN
                  R=X(1,1)
                  RR=R*R
                  SUM=SUM+X(nu,1)**2+RR*X(nu,2)**2
                  IF(NJT.GT.2) SUM=SUM+X(nu,3)**2
                ELSE IF(ITYP10(1).EQ.3) THEN
                  R=X(1,1)
                  RR=R*R
                  RC=R*DCOS(X(1,3))
                  RRC=RC*RC
                  SUM=SUM+X(nu,1)**2+RRC*X(nu,2)**2+RR*X(nu,3)**2
                ELSE IF(ITYP10(1).EQ.4) THEN
                  AA=FOCUS*FOCUS
                  SLX=DSINH(X(1,1))
                  SMX=DSIN(X(1,2))
                  G1=AA*(SLX*SLX+SMX*SMX)
                  G3=AA* SLX*SLX*SMX*SMX
                  SUM=SUM+G1*(X(nu,1)**2+X(nu,2)**2)
                  IF(NJT.GT.2) SUM=SUM+G3*X(nu,3)**2
                ENDIF
                DSDXI(ni)=DSQRT(SUM)
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' DSDXI(ni): '',3E12.3)')
     '            (DSDXI(ni),ni=1,NIT(nb1))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF

C cpb 5/12/96 re-adding arc-length scaling
c cpb 8/10/94 changing around nbi(nb) = 4 and 5
C            IF(NBI(nb1).EQ.3.OR.NBI(nb1).EQ.4) THEN
            IF(NBI(nb1).EQ.3.OR.(NBI(nb1).GE.5.AND.NBI(nb1).LE.7)) THEN
C ***         Store nodal derivatives wrt arc-length
              DO nj=1,NJ_LOC(0,0,nr)
                XP(1,1,nj,NODE)=X(1,nj)
                nb=NBJ(nj,ne)
!AJP 12/7/96 Check on nb added becuase if a fibre field has been
!defined in another region nj_loc(0,0,nr) may exceed what is needed for
!the current region and nb may not be defined.
                IF(nb.GT.0) THEN
                  DO nk=2,NKT(0,nb)
                    IF(nk.EQ.2) THEN
                      IF(DABS(DSDXI(1)).GT.1.d-6) THEN
                        XP(2,1,nj,NODE)=X(2,nj)/DSDXI(1)
                      ELSE
                        WRITE(OP_STRING,
     '                    '('' Warning: Arc length deriv is zero'')')
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                    ELSE IF(nk.EQ.3) THEN
                      IF(DABS(DSDXI(2)).GT.1.d-6) THEN
                        XP(3,1,nj,NODE)=X(4,nj)/DSDXI(2)
                      ELSE
                        WRITE(OP_STRING,
     '                    '('' Warning: Arc length deriv is zero'')')
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                    ELSE IF(nk.EQ.4) THEN
                      IF(NOCROSS) THEN
                        XP(4,1,nj,NODE)=0.0D0
                      ELSE
                        IF(DABS(DSDXI(1)).GT.1.d-6.AND.DABS(DSDXI(2))
     '                     .GT.1.d-6) THEN
                          XP(4,1,nj,NODE)=X(6,nj)/(DSDXI(1)*DSDXI(2))
                        ELSE
                          WRITE(OP_STRING,
     '                      '('' Warning: Arc length deriv is zero'')')
                          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF !nb
              ENDDO
            ELSE
              DO nj=1,NJ_LOC(0,0,nr)
                XP(1,1,nj,NODE)=X(1,nj)
              ENDDO
            ENDIF

          ELSE !Node exits, but possibly only in another region
            FOUND=.FALSE.
            nonode=1
            DO WHILE((nonode.LE.NPNODE(0,nr)).AND.(.NOT.FOUND))
              NP_TEST=NPNODE(nonode,nr)
              IF(np_test.eq.np)THEN !Node np is in current region
                FOUND=.TRUE.
              ELSE
                nonode=nonode+1
              ENDIF
            ENDDO
            IF(.NOT.FOUND)THEN !Update arrays to include np
              NPNODE(0,nr)=NPNODE(0,nr)+1
              NPNODE(NPNODE(0,nr),nr)=np
              IF(np.GT.NPT(nr))THEN
                NPT(nr)=np
              ENDIF
            ENDIF
          ENDIF !.not.exist
          IF(i3.EQ.1)NODE1=NODE !Store first node
        ENDDO !End of loop over both directions
        IF(NKT(1,nb1).EQ.1)THEN !Apex at node 1
          NPNE1(1)=NPNE(1,nb1,ne) !Simplex
          NPNE1(2)=NODE1
          NPNE1(3)=NODE
          NPNE2(1)=NODE1 !Bicubic hermite
          NPNE2(2)=NODE
          NPNE2(3)=NPNE(2,nb1,ne)
          NPNE2(4)=NPNE(3,nb1,ne)
        ELSE IF(NKT(3,nb1).EQ.1)THEN !Apex at node 3
          NPNE1(1)=NPNE(1,nb1,ne) !Bicubic hermite
          NPNE1(2)=NPNE(2,nb1,ne)
          NPNE1(3)=NODE1
          NPNE1(4)=NODE
          NPNE2(1)=NODE1 !Simplex
          NPNE2(2)=NODE
          NPNE2(3)=NPNE(3,nb1,ne)
        ENDIF

        IF(DOP) THEN
          WRITE(OP_STRING,'('' NPNE1: '',27I3)') (NPNE1(nn),nn=1,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NPNE2: '',27I3)') (NPNE2(nn),nn=1,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
          DO nb=1,NBFT
          IF(NNT(nb).LE.3)THEN
!           Same type of basis function, or a basis function with
!           no more nodes than the simplex
            DO nn=1,NNT(nb)
              NPNE(nn,nb,ne)=NPNE1(nn)
              NPNE(nn,nb,ne_NEW)=NPNE2(nn)
            ENDDO
          ELSE
!           Possibly a Lagrange product basis function
            IF(NKT(1,nb1).EQ.1)THEN
              NPNE(1,nb,ne)=NPNE1(1)
              NPNE(1,nb,ne_NEW)=NPNE2(1)
              NPNE(2,nb,ne)=NPNE1(1)
              NPNE(2,nb,ne_NEW)=NPNE2(2)
              NPNE(3,nb,ne)=NPNE1(2)
              NPNE(3,nb,ne_NEW)=NPNE2(3)
              NPNE(4,nb,ne)=NPNE1(3)
              NPNE(4,nb,ne_NEW)=NPNE2(4)
            ELSEIF(NKT(3,nb1).EQ.1)THEN
              NPNE(1,nb,ne)=NPNE1(1)
              NPNE(1,nb,ne_NEW)=NPNE2(1)
              NPNE(2,nb,ne)=NPNE1(2)
              NPNE(2,nb,ne_NEW)=NPNE2(2)
              NPNE(3,nb,ne)=NPNE1(3)
              NPNE(3,nb,ne_NEW)=NPNE2(3)
              NPNE(4,nb,ne)=NPNE1(4)
              NPNE(4,nb,ne_NEW)=NPNE2(3)
            ENDIF
            IF(NNT(nb).GT.4)THEN
              WRITE(OP_STRING,
     '          '('' Warning: Cant update NPNE correctly in DIVS1'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nn=5,NNT(nb)
                NPNE(nn,nb,ne)=NPNE(4,nb,ne)
                NPNE(nn,nb,ne_NEW)=NPNE(4,nb,ne_NEW)
              ENDDO
            ENDIF
          ENDIF
        ENDDO

C ***   Update scaling factors in current element & new element
        IF(NKT(1,nb1).EQ.1)THEN !Apex at node 1
!         !ne_NEW is the bicubic hermite element.  This code relies
!         !on the user having set up a bicubic hermite element with
!         !the correct nodes in it (i.e. a collapsed node element).
          DO nb=1,NBFT
C ***       Xi(1) derivs
            SE( 2,nb,ne_NEW)=0.50D0*(SE( 2,nb,ne)+SE( 10,nb,ne))
            SE( 6,nb,ne_NEW)=0.50D0*(SE( 6,nb,ne)+SE( 14,nb,ne))
            SE(10,nb,ne_NEW)=     SE(10,nb,ne)
            SE(14,nb,ne_NEW)=     SE(14,nb,ne)
            SE( 3,nb,ne)    =       SE( 3,nb,ne)
            SE( 7,nb,ne)    =       SE( 7,nb,ne)

C ***       Xi(2) derivs
            SE( 3,nb,ne_NEW)=0.25D0*(SE( 3,nb,ne)+SE( 11,nb,ne))
            SE( 7,nb,ne_NEW)=0.25D0*(SE( 7,nb,ne)+SE( 15,nb,ne))
            SE(11,nb,ne_NEW)=0.50D0* SE(11,nb,ne)
            SE(15,nb,ne_NEW)=0.50D0* SE(15,nb,ne)
            SE( 4,nb,ne)    =0.50D0* SE( 4,nb,ne)
            SE( 8,nb,ne)    =0.50D0* SE( 8,nb,ne)

C ***       Xi(1),Xi(2) derivs
            SE( 4,nb,ne_NEW)=SE( 2,nb,ne_NEW)*SE( 3,nb,ne_NEW)
            SE( 8,nb,ne_NEW)=SE( 6,nb,ne_NEW)*SE( 7,nb,ne_NEW)
            SE(12,nb,ne_NEW)=SE(10,nb,ne_NEW)*SE(11,nb,ne_NEW)
            SE(16,nb,ne_NEW)=SE(14,nb,ne_NEW)*SE(15,nb,ne_NEW)
            SE( 5,nb,ne)  =SE( 3,nb,ne)  *SE( 4,nb,ne)
            SE( 9,nb,ne)  =SE( 7,nb,ne)  *SE( 8,nb,ne)
            ns=0
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                ns=ns+1
                SE(ns,nb,ne_NEW)=SE(ns,nb,ne)
              ENDDO
            ENDDO
            IF(DOP) THEN
              DO ns=1,NST(nb)
                WRITE(OP_STRING,'('' ne='',I4,'' nb='',i1,'' ns='','
     '            //'I2,'' SE='',E12.3)') ne,nb,ns,SE(ns,nb,ne)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO
        ELSE IF(NKT(3,nb1).EQ.3)THEN !Apex at node 3
!         !ne_NEW is the hermite simplex element.  This code relies
!         !on the user having set up a bicubic hermite element with
!         !the correct nodes in it (i.e. a collapsed node element).
          DO nb=1,NBFT
C ***       Xi(1) derivs
            SE( 3,nb,ne_NEW)  =        SE( 3,nb,ne)
            SE( 7,nb,ne_NEW)  =        SE( 7,nb,ne)
            SE(10,nb,ne)      =0.50D0*(SE( 2,nb,ne)+SE( 10,nb,ne))
            SE(14,nb,ne)      =0.50D0*(SE( 6,nb,ne)+SE( 14,nb,ne))
            SE( 2,nb,ne)      =        SE( 2,nb,ne)
            SE( 6,nb,ne)      =        SE( 6,nb,ne)

C ***       Xi(2) derivs
            SE( 4,nb,ne_NEW)  =0.50D0* SE( 4,nb,ne)
            SE( 8,nb,ne_NEW)  =0.50D0* SE( 8,nb,ne)
            SE( 3,nb,ne)      =0.25D0*(SE( 3,nb,ne)+SE( 11,nb,ne))
            SE( 7,nb,ne)      =0.25D0*(SE( 7,nb,ne)+SE( 15,nb,ne))
            SE(11,nb,ne)      =0.50D0* SE(11,nb,ne)
            SE(15,nb,ne)      =0.50D0* SE(15,nb,ne)

C ***       Xi(1),Xi(2) derivs
            SE( 5,nb,ne_NEW)  =SE( 3,nb,ne_NEW)  *SE( 4,nb,ne_NEW)
            SE( 9,nb,ne_NEW)  =SE( 7,nb,ne_NEW)  *SE( 8,nb,ne_NEW)
            SE( 4,nb,ne)=SE( 2,nb,ne)*SE( 3,nb,ne)
            SE( 8,nb,ne)=SE( 6,nb,ne)*SE( 7,nb,ne)
            SE(12,nb,ne)=SE(10,nb,ne)*SE(11,nb,ne)
            SE(16,nb,ne)=SE(14,nb,ne)*SE(15,nb,ne)

            ns=0
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                ns=ns+1
                SE(ns,nb,ne_NEW)=SE(ns,nb,ne)
              ENDDO
            ENDDO
            IF(DOP) THEN
              DO ns=1,NST(nb)
                WRITE(OP_STRING,'('' ne='',I4,'' nb='',i1,'' ns='','
     '            //'I2,'' SE='',E12.3)') ne,nb,ns,SE(ns,nb,ne)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDIF

C *** Form other arrays for new element
      NHE(ne_NEW)=NHE(ne)
      NW(ne_NEW,1)=NW(ne,1)
      NRE(ne_NEW)=NRE(ne)
      DO nm=1,NMM
        CE(nm,ne_NEW)=CE(nm,ne)
      ENDDO
      IF(IDRN.EQ.1) THEN
        DO nj=1,NJ_LOC(0,0,nr)
          NBJ(nj,ne_NEW)=NBJ(nj,ne)
        ENDDO
        IF(CALL_EQUA.OR.CALL_FIT.OR.CALL_OPTI) THEN
          DO nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
            DO nc=1,NCM
              NBH(nh,nc,ne_NEW)=NBH(nh,nc,ne)
            ENDDO !nc
          ENDDO !nhx
        ENDIF
        DO nb=1,NBFT
          DO nn=1,NNT(nb)
            DO nj=1,NJ_LOC(0,0,nr)
              NVJE(nn,nb,nj,ne_NEW)=NVJE(nn,nb,nj,ne)
            ENDDO
          ENDDO
        ENDDO
C KAT 23Feb01: now handled by NKJE below
C        DO nb=1,NBFT
C          DO nn=1,NNT(nb)
C            DO nk=1,NKT(nn,nb)
C              NKE(nk,nn,nb,ne_NEW)=NKE(nk,nn,nb,ne)
C            ENDDO
C          ENDDO
C        ENDDO
      ELSE IF(IDRN.EQ.2)THEN
        IF(NKT(1,nb1).EQ.1)THEN !Apex at node 1
          DO nj=1,NJ_LOC(0,0,nr)
            NBJ(nj,ne_NEW)=NB_HERMITE !????? Very specialised
          ENDDO
          IF(CALL_EQUA.OR.CALL_FIT.OR.CALL_OPTI) THEN
            DO nhx=1,NHE(ne)
              nh=NH_LOC(nhx,nx)
              DO nc=1,NCM
                NBH(nh,nc,ne_NEW)=NB_HERMITE !????? Very specialised
              ENDDO !nc
            ENDDO !nhx
          ENDIF
          DO nb=1,NBFT
            DO nn=1,NNT(nb)
              DO nj=1,NJ_LOC(0,0,nr)
                NVJE(nn,nb,nj,ne_NEW)=NVJE(nn,nb,nj,ne)
              ENDDO
            ENDDO
          ENDDO
C KAT 23Feb01: now handled by NKJE below
C          DO nb=1,NBFT
C            DO nn=1,NNT(nb)
C              DO nk=1,NKT(nn,nb)
C                NKE(nk,nn,nb,ne_NEW)=NKE(nk,nn,nb,ne)
C              ENDDO
C            ENDDO
C          ENDDO
        ELSE IF(NKT(3,nb1).EQ.1)THEN !Apex at node 3
          DO nj=1,NJ_LOC(0,0,nr)
            NBJ(nj,ne_NEW)=NBJ(nj,ne)
            NBJ(nj,ne)=NB_HERMITE
          ENDDO
          IF(CALL_EQUA.OR.CALL_FIT.OR.CALL_OPTI) THEN
            DO nhx=1,NHE(ne)
              nh=NH_LOC(nhx,nx)
              DO nc=1,NCM
                NBH(nh,nc,ne_NEW)=NBH(nh,nc,ne)
                NBH(nh,nc,ne)=NB_HERMITE
              ENDDO !nc
            ENDDO !nhx
          ENDIF
          DO nb=1,NBFT
            DO nn=1,NNT(nb)
              DO nj=1,NJ_LOC(0,0,nr)
                NVJE(nn,nb,nj,ne_NEW)=NVJE(nn,nb,nj,ne)
              ENDDO
            ENDDO
          ENDDO
C KAT 23Feb01: now handled by NKJE below
C          DO nb=1,NBFT
C            DO nn=1,NNT(nb)
C              DO nk=1,NKT(nn,nb)
C                NKE(nk,nn,nb,ne_NEW)=NKE(nk,nn,nb,ne)
C              ENDDO
C            ENDDO
C          ENDDO
        ENDIF
      ENDIF

      DO nj=1,NJ_LOC(0,0,nr)
        nb1=NBJ(nj,ne)
        DO nn=1,NNT(nb1)
          DO nk=1,NKT(nn,nb1)
            NKJE(nk,nn,nj,ne_NEW)=NKJE(nk,nn,nj,ne)
          ENDDO !nk
        ENDDO !nn
      ENDDO !nj

      NPSTART=NPT(0)+1
      NET(0)=ne_NEW   !is new highest element # in mesh
      NET(nr)=ne_NEW  !is new highest element # in region nr
      CALL ASSERT(NET(nr).LE.NEM,
     '  '>>Too many elements - increase NEM',ERROR,*9999)
      NEELEM(0,nr)=NEELEM(0,nr)+1 !is new #elements in region nr
      NEELEM(NEELEM(0,nr),nr)=ne_NEW !is new element number

!     !Update element list in each group
      CALL UPGREL(ne,ne_NEW,ERROR,*9999)

! Check that the global nodes of new element ne_NEW are in the
! node list for the current region
      nb=NBJ(1,ne_NEW)
      DO nn=1,NNT(nb)
        np=NPNE(nn,nb,ne_NEW)
        EXIST=.FALSE.
        DO nonode=1,NPNODE(0,nr)
          IF(NPNODE(nonode,nr).eq.np) EXIST=.TRUE.
        ENDDO
        IF(.NOT.EXIST) THEN
          NPNODE(0,nr)=NPNODE(0,nr)+1
          NPNODE(NPNODE(0,nr),nr)=np
        ENDIF
      ENDDO

C *** Diagnostic output
      IF(DOP) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nb=1,NBFT
            WRITE(OP_STRING,'('' NPNE(nn,'',i1,'','',I4,'
     '        //'''): '',20I4)') nb,ne,(NPNE(nn,nb,ne),nn=1,NNT(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
        DO np=1,NPT(nr)
          DO nj=1,NJT
            WRITE(OP_STRING,'('' XP(nk,1,'',i1,'','',I4,''): '','
     '        //'10E11.4)') nj,np,(XP(nk,1,nj,np),nk=1,NKJ(nj,np))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('DIVS1')
      RETURN
 9999 CALL ERRORS('DIVS1',ERROR)
      CALL EXITS('DIVS1')
      RETURN 1
      END

c cpb 24/11/96 Adding sector refinement

