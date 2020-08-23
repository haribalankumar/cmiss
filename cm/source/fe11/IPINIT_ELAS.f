      SUBROUTINE IPINIT_ELAS(NAN,NBH,NBHF,NEELEM,NEL,NENP,NFF,
     '  NHE,NHP,NKEF,NKH,NKHE,NLL,NNF,NNL,NPF,NPL,NPLIST,NPNE,
     '  NPNODE,nr,NVHE,NVHP,NW,nx,NYNE,NYNP,DF,DL,PG,SE,WG,XP,
     '  YP,ZA,ZP,FIX,NOFIX,ERROR,*)

C#### Subroutine: IPINIT_ELAS
C###  Description:
C###    Inputs initial conditions and boundary conditions for elasticity
C###    problems.

C For linear elasticity
C**** YP(ny1v,1) contains essential b.c.s,defined by FIX(ny1v,1).
C**** YP(ny2v,1)    "     natural           "        FIX(ny2v,1).
C For finite elasticity
C**** FIX(ny,1) defines prescribed dep var/force bdry conds
C**** FIX(ny,2) defines prescribed incremental dep var/force bdry conds
C**** FIX(ny,3) defines prescribed dep var/force initial conds
C****  YP(ny,1) has current equilibrium solution
C****  YP(ny,2) has prescribed dep var/force increms set by FIX(ny,1)
C****  YP(ny,3) has prescribed initial equilibrium solution
C****  YP(ny,4) has current set of equilibrium equation residuals
C**** This subroutine initialises/inputs FIX(ny,1), FIX(ny,2)
C****   and FIX(ny,3), YP(ny,1), YP(ny,2) and YP(ny,3)

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='IPINIT_ELAS')
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),NHE(NEM),NHP(NPM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),
     '  NLL(12,NEM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NPF(9,NFM),NPL(5,0:3,NLM),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 
     '  DF(NFM),DL(3,NLM),PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM),NOFIX
!     Local Variables
      INTEGER INFO,iy,kf,kl,loop,LOOPT,
     '  n,N1,N1ELEM,N1NODE,N2ELEM,N2NODE,
     '  na,nb,nc,ne,NELIST(0:NEM),nf,nh,NHT,nhx,nj,nk,NKHF(NKM,NNM,NHM),
     '  NKHMAX,nl,nn,noelem,nonode,NOQUES,np,NPNF(NNM,NBFM),nt,NUMGT1,
     '  NUMFIXED(3),nv,NVHF(NNM,NBFM,NHM),NVHMAX,ny,ny_first
      REAL*8 COEFF(4),fn,SF(NSM,NBFM),ZP_GROUP(NKM,NVM)
      CHARACTER CHAR1*100,CHAR2*100
      LOGICAL CONNECTED,FILEIP,FIX_GROUP(NKM,NVM),
     '  GROUPED_FACE,GROUPED_LINE,INLIST,INTEGRATED,ROTATE,TRANSLATE

      CALL ENTERS(ROUTINENAME,*9999)

      FILEIP=.FALSE.
      NOQUES=0

      IF(.NOT.NOFIX.AND.IOTYPE.NE.3) THEN
        ! Initialise FIX
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          DO nc=1,NCT(nr,nx)
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,nc)
                DO nk=1,NKH(nh,np,nc)
C KAT 2003-04-14: FIX does not have an nc index so the global ny seems
C the appropriate ny to use.  What was the logic behind also setting ny'
C s for
C the local row and column numbers?  This was already turned off if
C KTYP5G(nr).EQ.1 (Contact coupled initialisation).
C                       DO nrc=0,2
                  ny=NYNP(nk,nv,nh,np,0,nc,nr)
                  FIX(ny,1)=.FALSE.
                  ! Don't know what this means for finite elasticity,
                  ! but, for linear, SOLVE1 uses this for mixed boundary
                  ! conditions.
                  FIX(ny,2)=.FALSE.
                  IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
                    FIX(ny,3)=.TRUE. !always specify i.c.s
C KAT 2003-04-14: We certainly shouldn't have to initialize YP for all
C NIYM elements.  Why initialize it at all?  It is set in ZPYP calls
C below.
C                         DO loop=1,NIYM
C                           YP(ny,loop)=0.0d0
C                         ENDDO !loop
C                       ENDDO !nrc
                  ENDIF ! finite
                ENDDO !nk
              ENDDO !nv
            ENDDO !nh
          ENDDO !nc
        ENDDO !nonode
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nc=1,NCT(nr,nx)
            DO nhx=1,NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
              DO na=1,NAT(NBH(nh,nc,ne))
                ny=NYNE(na,nh,0,nc,ne)
                FIX(ny,1)=.FALSE.
                IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
                  FIX(ny,2)=.FALSE.
                  FIX(ny,3)=.TRUE. !always specify initial conditions
                ENDIF ! finite
              ENDDO !na
            ENDDO !nh
          ENDDO !nc
          IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
            ! Initialise NW
            NW(ne,1)=1
          ENDIF ! finite
        ENDDO !noelem
      ENDIF

      IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
        FORMAT='('' Specify option [1]:'''//
     '    '/''   (1) Initial displacements all zero'''//
     '    '/''   (2) Read in initial conditions    '''//
     '    '/''   (3) Restart from previous solution'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=2 !to write out current solution
          KTYP5=2
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IONE,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP5=IDATA(1)
      ELSE ! linear
        KTYP5=0 !no initial conditions
      ENDIF

      IF(KTYP5.EQ.0) THEN !no initial conditions
        LOOPT=2
      ELSEIF(KTYP5.EQ.1) THEN !Initial displacements all zero
        LOOPT=2
C       Initialise ZP,ZA
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nc=1,NCT(nr,nx)
            IF(nc.EQ.1) THEN ! dependent var
              DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,nhx,nr)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,nc)
                  DO nk=1,NKH(nh,np,nc)
                    ZP(nk,nv,nh,np,nc)=XP(nk,nv,nj,np)
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nh
            ELSE  !not reading in initial conditions
              DO nhx=1,NHP(np)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,nc)
                  DO nk=1,NKH(nh,np,nc)
                    ZP(nk,nv,nh,np,nc)=0.0d0
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nh
            ENDIF !nc
          ENDDO !nc
        ENDDO !nonode (np)
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nc=1,NCT(nr,nx)
            DO nhx=1,NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
              DO na=1,NAT(NBH(nh,nc,ne))
                ZA(na,nh,nc,ne)=0.0d0
              ENDDO !na
            ENDDO !nh
          ENDDO !nc
        ENDDO !noelem (ne)

        ! Set current solution
        CALL ZPYP(1,NBH,NEELEM,NHE,NHP,NKH,
     '    NPNODE,nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)

      ELSE IF(KTYP5.EQ.2) THEN !Read in initial conditions
        LOOPT=3
C TVK 12/01/2000 Set logical var to calculate correct equilm press
        IF(KTYP52(nr).EQ.4) THEN
          EQUIM_PRESSURE=.TRUE.
        ENDIF
      ELSE IF(KTYP5.EQ.3) THEN !Restart from previous solution
        LOOPT=2
        RESTAR=.TRUE.
      ENDIF

      DO loop=2,LOOPT
        IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
          IF(loop.EQ.2) THEN
            iy=2 ! boundary conditions (increments)
          ELSE ! loop=3
            iy=1 !initial conditions/current solution
          ENDIF
        ELSE ! linear, loop=2
          iy=1 ! boundary conditions
        ENDIF
        IF(IOTYPE.EQ.3) THEN
          CALL YPZP(iy,NBH,NEELEM,NHE,NHP,NKH,
     '      NPNODE,nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        ENDIF

        DO nc=1,2  !displacement/force
          IF(loop.EQ.2) THEN !boundary conditions
            IF(nc.EQ.1) THEN !dependent vars
              IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
                FORMAT='(/'' Dependent variable boundary conditions:'')'
              ELSE ! linear elasticity
                FORMAT='(/'' Displacement boundary conditions:'')'
              ENDIF
            ELSE !forces
              FORMAT='(/'' Force boundary conditions:'')'
            ENDIF
          ELSE ! loop=3, finite, initial conditions
            IF(nc.EQ.1) THEN          !dependent vars
              FORMAT='(/'' Dependent variable initial conditions:'')'
            ELSE                      !forces
              FORMAT='(/'' Force initial conditions:'')'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,0,
     '      IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

          IF(ITYP1(nr,nx).EQ.4 ! linear elasticity
     '      .AND.nc.EQ.2) THEN ! forces
            FORMAT=
     '        '(/'' Specify whether the force values are [1]:'''//
     '        '/''   (1) Weighted integrals of force fields'''//
     '        '/''   (2) Tractions to be interpolated'//
     '        ' and integrated approximately'''//
     '        '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IONE,
     '        1,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '        *9999)
            INTEGRATED=IDATA(1).EQ.1
          ELSE
            INTEGRATED=.TRUE.
          ENDIF ! linear, forces

          IF(IOTYPE.NE.3) THEN
            ! Initialize ZP (all nc) and ZA (for nc=1)
            IF(ITYP4(nr,nx).EQ.2) THEN !Boundary elements
              !MLB FLux variables are ny's with niy=1 in BEM
              !    but niy=2 in FEM. Need to initialise FIX
              !    correctly for BEM.
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nhx=1,NHP(np)
                  nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
                  DO nv=1,NVHP(nh,np,nc)
                    DO nk=1,NKH(nh,np,nc)
                      ny=NYNP(nk,nv,nh,np,0,nc,nr)
                      FIX(ny,1)=.FALSE. !nc=1,2
                    ENDDO       !nk
                  ENDDO         !nv
                ENDDO           !nh
              ENDDO             !nonode (np)
            ENDIF
            IF(loop.EQ.3.AND.nc.EQ.1) THEN
              ! Read in dependent var initial conds
C             Set ZP to XP for nodes not initialised
C MPN 8-Dec-94: this seems to do it for all nodes. may need fixing
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,nhx,nr)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,nc)
                    DO nk=1,NKH(nh,np,nc)
                      ZP(nk,nv,nh,np,nc)=XP(nk,nv,nj,np)
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nh
              ENDDO !nonode (np)
            ELSE  !not dependent var initial conditions
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nhx=1,NHP(np)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,nc)
                    DO nk=1,NKH(nh,np,nc)
                      ZP(nk,nv,nh,np,nc)=0.0d0
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nh
              ENDDO !nonode (np)
            ENDIF !loop & nc
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nhx=1,NH_LOC(0,nx)
                nh=NH_LOC(nhx,nx)
                DO na=1,NAT(NBH(nh,nc,ne))
                  ZA(na,nh,nc,ne)=0.0d0
                ENDDO !na
              ENDDO !nh
            ENDDO !noelem (ne)
          ENDIF !iotype.ne.3

          IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
            NHT=NHT50(KTYP52(nr),KTYP51(nr))
          ELSE ! linear
            NHT=NHP(np)
          ENDIF
          DO nhx=1,NHT
            nh=NH_LOC(nhx,nx)
            IF(nc.EQ.1.OR.nhx.LE.NJT) THEN
              WRITE(CHAR1,'(I1)') nhx
              IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
                FORMAT='(/'' Dependent variable/equation number '
     '            //CHAR1(1:1)//' : '')'
              ELSE              ! linear
                IF(nc.EQ.1) THEN !dependent vars
                  FORMAT='(/'' Displacement component '
     '              //CHAR1(1:1)//' : '')'
                ELSE            !forces
                  FORMAT='(/'' Force component '//CHAR1(1:1)//' : '')'
                ENDIF
              ENDIF             ! finite/linear
              CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,IMIN,
     '          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '          ERROR,*9999)
              nb=NBH(nh,nc,NEELEM(1,nr))
              IF(IOTYPE.EQ.3) THEN
                N1NODE=1
                IF(loop.EQ.3) nonode=0 !ini node loop for writing ics
              ENDIF
              IF(NNT(nb).EQ.0) THEN ! no nodes for interpolation
                NPLIST(0)=0 ! don't request nodes.
              ELSE
                NPLIST(0)=-1    ! so that we start the loop
              ENDIF
              DO WHILE(NPLIST(0).NE.0)
                FORMAT='($,'' Enter node #s/name [EXIT]: '',I5)'
                IF(IOTYPE.EQ.3) THEN
                  IF(loop.EQ.2) THEN !boundary conditions
                    DO nonode=N1NODE,NPNODE(0,nr)
                      np=NPNODE(nonode,nr)
                      N2NODE=nonode
                      DO nv=1,NVHP(nh,np,nc)
                        DO nk=1,NKH(nh,np,nc)
                          IF(FIX(NYNP(nk,nv,nh,np,0,nc,nr),1)) GOTO 6200
                        ENDDO !nk
                      ENDDO !nv
                    ENDDO !nonode (np)
                    N2NODE=0
 6200               IF(N2NODE.EQ.0) THEN
                      IDATA(1)=0
                    ELSE
                      IDATA(1)=np
                      N1NODE=N2NODE+1
                    ENDIF
                  ELSE IF(loop.EQ.3) THEN !initial conditions
                    nonode=nonode+1
                    IF(nonode.LE.NPNODE(0,nr)) THEN
                      np=NPNODE(nonode,nr)
                      IDATA(1)=np
                    ELSE
                      IDATA(1)=0 !to terminate node loop
                    ENDIF
                  ENDIF !loop
                  IDATA(0)=1 !write out one node at a time
                ENDIF !iotype=3

 6500           CDATA(1)='NODES' !for use with group input
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,
     '            NOQUES,FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '            0,IDATA,IZERO,0,NPT(nr),
     '            LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                IF(IDATA(1).EQ.0) THEN ! default exit
                  NPLIST(0)=0
                ELSE
                  NPLIST(0)=IDATA(0)
                  DO n=1,IDATA(0)
                    NPLIST(n)=IDATA(n)
                    np=IDATA(n)
                    IF(iotype.NE.3) THEN
                      IF(.NOT.INLIST(np,NPNODE(1,nr),
     '                  NPNODE(0,nr),N1)) THEN
                        WRITE(OP_STRING,'('' Node '',I5,'' does not '
     '                    //'belong to the current region'')') np
                        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                        GOTO 6500
                      ENDIF
                    ENDIF !iotype.NE.3
                  ENDDO !n
                ENDIF !IDATA(1)

                ! Define bdry condition for maximum number of versions and
                ! derivatives in group.
                NVHMAX=0
                NKHMAX=0
                DO n=1,NPLIST(0)
                  np=NPLIST(n)
                  IF(NVHP(nh,np,nc).GT.NVHMAX) NVHMAX=NVHP(nh,np,nc)
                  IF(NKH(nh,np,nc).GT.NKHMAX) NKHMAX=NKH(nh,np,nc)
                ENDDO !np
                DO nv=1,NVHMAX
                  IF(NVHMAX.GT.1) THEN
                    WRITE(CHAR1,'(I2)') nv
                    FORMAT='('' For version number '//CHAR1(1:2)
     '                //':'')'
                    CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,
     '                IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '                RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                  ENDIF
                  DO nk=1,NKHMAX
                    ny=NYNP(nk,nv,nh,np,0,nc,nr)
C KAT 2003-04-10: If we are setting displacement boundary conditions and
C there is more than one derivative then we need to determine which
C derivatives the user wants to fix.  With finite elasticity a question
C has historically been asked for each derivative even if there is only
C one derivative and also for natural boundary conditions.  For linear
C elasticity, we don't ask the unnecessary questions, which gives some
C backward compatibility with files from when these questions were not
C asked.
                    IF(loop.EQ.2 ! boundary conditions
     '                .AND.(ITYP1(nr,nx).EQ.5 ! finite elasticity
     '                .OR.(nc.EQ.1.AND.NKHMAX.GT.1))) THEN
                      ! Ask the user
                      IF(nc.EQ.1) THEN !dependent vars
                        IF(nk.EQ.1) THEN
                          IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
                            FORMAT=
     '                        '($,'' Do you want to prescribe the'//
     '                        ' dependent variable [Y]? '',A)'
                          ELSE ! linear
                            FORMAT=
     '                        '($,'' Do you want to prescribe the'//
     '                        ' displacement [Y]? '',A)'
                          ENDIF ! finite/linear
                          ADEFLT(1)='Y'
                        ELSE IF(nk.GT.1) THEN
                          WRITE(CHAR1,'(I1)') nk
                          FORMAT='($,'' Do you want to prescribe'
     '                      //' derivative number '//CHAR1(1:1)
     '                      //' [N]? '',A)'
                          ADEFLT(1)='N'
                        ENDIF
                      ELSE ! nc=2 forces
                        IF(nk.EQ.1) THEN
                          FORMAT=
     '                      '($,'' Do you want to prescribe the'//
     '                      ' force component [Y]? '',A)'
                          ADEFLT(1)='Y'
                        ELSE IF(nk.GT.1) THEN
                          WRITE(CHAR1,'(I1)') nk
                          FORMAT=
     '                      '($,'' Do you want to prescribe moment'
     '                      //' number '//CHAR1(1:1)//' [N]? '',A)'
                          ADEFLT(1)='N'
                        ENDIF
                      ENDIF
                      IF(IOTYPE.EQ.3) THEN
                        IF(FIX(ny,1)) THEN
                          ADATA(1)='Y'
                        ELSE
                          ADATA(1)='N'
                        ENDIF
                      ENDIF
                      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,
     '                  NOQUES,FILEIP,FORMAT,1,
     '                  ADATA,ADEFLT,CDATA,CDEFLT,
     '                  0,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '                  RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
                      FIX_GROUP(nk,nv)=ADATA(1).EQ.'Y'
                    ELSE
                      ! No need to ask user.
                      FIX_GROUP(nk,nv)=.TRUE.
                    ENDIF
                    IF(FIX_GROUP(nk,nv)) THEN
                      !initial conditions or specified bc
                      IF(loop.EQ.2) THEN       !incremental bcs
                        RDEFLT(1)=0.0d0
                        IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
                          FORMAT=
     '                      '($,'' The increment is [0.0]: '',G25.17)'
                        ELSEIF(nc.EQ.1) THEN ! linear, dependent vars
                          IF(nk.EQ.1) THEN
                            FORMAT=
     '                        '($,'' The displacement is [0.0]: '','//
     '                        'G12.5)'
                          ELSE
                            WRITE(FORMAT,'(A,I1,A)')
     '                        '($,'' The value of displacement '//
     '                        'derivative ',nk,' is [0.0]:'',G12.5)'
                          ENDIF
                        ELSE IF(nc.EQ.2) THEN !linear, forces
                          IF(nk.EQ.1) THEN
                            FORMAT='($,'' The force is [0.0]: '','
     '                        //'G12.5)'
                          ELSEIF(INTEGRATED) THEN !weighted integrals
                            WRITE(FORMAT,'(A,I1,A)')
     '                        '($,'' The value of moment number ',
     '                        nk,' is [0.0]:'',G12.5)'
                          ELSE !tractions to be interpolated
                            WRITE(FORMAT,'(A,I1,A)')
     '                        '($,'' The value of force '//
     '                        'derivative ',nk,' is [0.0]:'',G12.5)'
                          ENDIF
                        ENDIF
                      ELSE ! loop=3, initial conditions
                        RDEFLT(1)=ZP(nk,nv,nh,np,nc)
C                        CHAR1=CFROMR(RDEFLT(1),'(D12.5)')
                        WRITE(CHAR1,'(D12.5)') RDEFLT(1)
                        IF(nk.EQ.1) THEN
                          FORMAT='($,'' The init. value of the '
     '                      //'dependent variable is ['//CHAR1(1:12)
     '                      //']: '',G25.17)'
                        ELSE IF(nk.GT.1) THEN
                          WRITE(CHAR2,'(I1)') nk
                          FORMAT='($,'' The init. value of derivative '
     '                      //'number '//CHAR2(1:1)//' is ['
     '                      //CHAR1(1:12)//']: '',G25.17)'
                        ENDIF
                      ENDIF ! loop
                      IF(IOTYPE.EQ.3) RDATA(1)=ZP(nk,nv,nh,np,nc)
                      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '                  FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '                  0,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '                  RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
                      IF(IOTYPE.NE.3) THEN
                        ZP_GROUP(nk,nv)=RDATA(1)
                      ENDIF ! IOTYPE==3
                    ENDIF ! specified
                  ENDDO !nk
                ENDDO !nv
                IF(IOTYPE.NE.3) THEN
C                 Apply bdry conditions to all nodes in group
                  DO n=1,NPLIST(0)
                    np=NPLIST(n)
                    DO nv=1,NVHP(nh,np,nc)
                      DO nk=1,NKH(nh,np,nc)
                        ny=NYNP(nk,nv,nh,np,0,nc,nr)
                        IF(loop.EQ.2) THEN ! boundary conditions
                          FIX(ny,1)=FIX_GROUP(nk,nv)
                          IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
C??? What is the purpose of FIX(ny,2)?
                            FIX(ny,2)=FIX(ny,1)
                          ENDIF
                        ENDIF ! loop==2
                        IF(FIX_GROUP(nk,nv)) THEN
C NEW CS 29/10/97
                          IF(nc.EQ.1.OR.INTEGRATED) THEN
                            ! dependent var or integrated forces
                            ZP(nk,nv,nh,np,nc)=ZP_GROUP(nk,nv)
                          ELSE ! to be interpolated
C                       Compute integrated force values
C                       Integrate over lines/faces connected to
C                       node that are in only one element (ie external)
                            IF(nk.GT.1) THEN
                              CALL FLAG_ERROR(-1,
     '                          'Traction integration'//
     '                          ' not implemented for derivatives')
                              ERROR=' '
                              GOTO 9999
                            ENDIF
                            DO noelem=1,NENP(np,0,nr)
                              ne=NENP(np,noelem,nr)
                              IF(NIT(nb).EQ.2) THEN
C                             loop over local lines of current element
                                DO kl=1,NLE(nb)
                                  nl=NLL(kl,ne)
                                  IF(nl.NE.0) THEN
                                    nn=0
                                    CONNECTED=.FALSE.
                                    DO WHILE ((nn.LT.NNL(0,kl,nb)).AND.
     '                                (.NOT.CONNECTED))
                                      nn=nn+1
                                      IF((np.EQ.NPL(nn+1,1,nl)).AND.
     '                                  (NEL(0,nl).EQ.1)) THEN
C!!! Is this connected to the correct version?
C!!! What prevents the contribution from being multiplied by the number
C!!! of versions?
                                        CONNECTED=.TRUE.
                                      ENDIF
                                    ENDDO !while
                                    IF(CONNECTED) THEN
                                      nt=2
C!!! Shouldn't the COEFF depend on the derivative number?
                                      CALL CALC_CONTRIB_COEFF(COEFF,
     '                                  nb,nb,nt,
     '                                  NPL(1,0,nl),PG,WG,ERROR,*9999)
C!!! Not valid if the Jacobian varies through the line?
                                      fn=COEFF(nn)*DL(3,nl)*ZP_GROUP(nk,
     '                                  nv)
C!!! What happens to the contributions to other equations with non-zero
C!!! weights on this line?
                                      ZP(nk,nv,nh,np,nc)=
     '                                  ZP(nk,nv,nh,np,nc)+fn
                                    ENDIF !CONNECTED
                                  ENDIF !nl.NE.0
                                ENDDO !kl
                              ELSE IF(NIT(nb).EQ.3) THEN
C                             loop over local faces of current element
                                DO kf=1,NFE(nb)
                                  nf=NFF(kf,ne)
                                  CALL CALC_FACE_INFORMATION_DEP(
     '                              NBH(1,1,ne),NBHF(1,1,nf),kf,NHE(ne),
     '                              NKHE(1,1,1,ne),NKEF,NKHF,NNF,
     '                              NPNE(1,1,ne),NPNF,NVHE(1,1,1,ne),
     '                              NVHF,nx,SE(1,1,ne),SF,ERROR,*9999)
                                  nn=0
                                  CONNECTED=.FALSE.
                                  DO WHILE ((nn.LT.NNF(0,kf,nb)).AND.
     '                              (.NOT.CONNECTED))
                                    nn=nn+1
                                    IF((np.EQ.NPNE(NNF(nn+1,kf,nb),
     '                                nb,ne)).AND.
     '                                (NPF(5,nf).EQ.1)) THEN
C!!! Is this connected to the correct version?
C!!! What prevents the contribution from being multiplied by the number
C!!! of versions?
                                      CONNECTED=.TRUE.
                                    ENDIF
                                  ENDDO !while
                                  IF(CONNECTED) THEN
C!!! Shouldn't the COEFF depend on the derivative number?
                                    CALL  CALC_CONTRIB_COEFF(COEFF,nb,
     '                                NBHF(NH_LOC(1,nx),1,nf),
     '                                NNF(0,kf,nb),NPL(1,0,1),PG,WG,
     '                                ERROR,*9999)
C!!! Not valid if the Jacobian varies through the face?
                                    fn=COEFF(nn)*DF(nf)*ZP_GROUP(nk,nv)
C!!! What happens to the contributions to other equations with non-zero
C!!! weights on this line?
                                    ZP(nk,nv,nh,np,nc)=
     '                                ZP(nk,nv,nh,np,nc)+fn
                                  ENDIF !CONNECTED
                                ENDDO !kf
                              ELSE ! NIT=1
                                CALL FLAG_ERROR(0,
     '                            'Traction integration'//
     '                            ' not implemented for 1D')
                                ERROR=' '
                                GOTO 9999
                              ENDIF !NIT(nb)
                            ENDDO !noelem
                          ENDIF !INTEGRATED/not
                        ENDIF !FIX_GROUP
                      ENDDO !nk
                    ENDDO !nv
                  ENDDO !np
                ENDIF !IO_TYPE!=3
              ENDDO ! NPLIST(0)!=0

              IF(nc.EQ.1.AND.NAT(nb).GT.0) THEN !aux vars for nc=1 only
                IF(IOTYPE.EQ.3) THEN
                  N1ELEM=1
                  IF(loop.EQ.3) noelem=0 !ini elem loop for writing ics
                ENDIF
 6720           FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
                IF(IOTYPE.EQ.3) THEN
                  IF(loop.EQ.2) THEN       !incremental bcs
                    DO noelem=N1ELEM,NEELEM(0,nr)
                      ne=NEELEM(noelem,nr)
                      N2ELEM=noelem
                      DO na=1,NAT(NBH(nh,nc,ne))
                        ny=NYNE(na,nh,0,nc,ne)
                        IF(FIX(ny,1)) GOTO 6740
                      ENDDO !na
                    ENDDO !noelem (ne)
                    N2ELEM=0
 6740               IF(N2ELEM.EQ.0) THEN
                      IDATA(1)=0
                    ELSE
                      IDATA(1)=ne
                      N1ELEM=N2ELEM+1
                    ENDIF
                  ELSE IF(loop.EQ.3) THEN  !initial conditions
                    noelem=noelem+1
                    IF(noelem.LE.NEELEM(0,nr)) THEN
                      ne=NEELEM(noelem,nr)
                      IDATA(1)=ne
                    ELSE
                      IDATA(1)=0 !to terminate element loop
                    ENDIF
                  ENDIF
                ENDIF
 6760           CDATA(1)='ELEMENTS' !for use with group input
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IZERO,
     '            0,NET(nr),
     '            LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                IF(IDATA(1).NE.0) THEN !not default exit
                  NELIST(0)=IDATA(0)
                  DO n=1,IDATA(0)
                    NELIST(n)=IDATA(n)
                    ne=IDATA(n)
                    IF(.NOT.INLIST(ne,NEELEM(1,nr),
     '                NEELEM(0,nr),N1)) THEN
                      WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '                  //'in the current region'')') ne
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                      GOTO 6760
                    ENDIF
                  ENDDO !n

C                 Define bdry condition for first element in group
                  ne=NELIST(1) !rest of group filled at end of nh loop
                  DO na=1,NAT(NBH(nh,nc,ne))
                    ny=NYNE(na,nh,0,nc,ne)
                    WRITE(CHAR1,'(I1)') na
                    IF(loop.EQ.2) THEN       !incremental bcs
                      FORMAT='($,'' Do you want to prescribe auxiliary'
     '                  //' variable/rhs number '//CHAR1(1:1)
     '                  //' [N]? '',A)'
                      IF(IOTYPE.EQ.3) THEN
                        IF(FIX(ny,1)) THEN
                          ADATA(1)='Y'
                        ELSE
                          ADATA(1)='N'
                        ENDIF
                      ENDIF
                      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,
     '                  FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,0,
     '                  IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '                  RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
                      IF(IOTYPE.NE.3) THEN
                        IF(ADATA(1).EQ.'Y') THEN
                          FIX(ny,1)=.TRUE.
                          FIX(ny,2)=.TRUE.
                        ELSE IF(ADATA(1).EQ.'N') THEN
                          FIX(ny,1)=.FALSE.
                          FIX(ny,2)=.FALSE.
                        ENDIF
                      ENDIF
                      IF(ADATA(1).EQ.'Y') THEN
                        FORMAT='($,'' The increment is [0.0]: '','
     '                    //'G25.17)'
                        IF(IOTYPE.EQ.3) RDATA(1)=ZA(na,nh,nc,ne)
                        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,
     '                    0,0,NOQUES,
     '                    FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '                    0,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '                    RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
                        IF(IOTYPE.NE.3) ZA(na,nh,nc,ne)=RDATA(1)
                        IF(KTYP57(nr).EQ.2) THEN
C                         pressure bc's entered as aux vars: set up NW
                          IF(NAN(3,na,NBH(nh,nc,ne)).EQ.-1) THEN
C                           Xi3=0 face pressure bc
                            IF(NW(ne,1).EQ.1) NW(ne,1)=2
                            IF(NW(ne,1).EQ.3) NW(ne,1)=4
                          ELSE IF(NAN(3,na,NBH(nh,nc,ne)).EQ.-2) THEN
C                           Xi3=1 face pressure bc
                            IF(NW(ne,1).EQ.1) NW(ne,1)=3
                            IF(NW(ne,1).EQ.2) NW(ne,1)=4
                          ENDIF
                        ENDIF
                      ENDIF
                    ELSE IF(loop.EQ.3) THEN !initial conditions
                      FORMAT='($,'' The initial value of auxiliary'//
     '                  ' variable number '//CHAR1(1:1)
     '                  //' is [0.0]: '',G25.17)'
                      IF(IOTYPE.EQ.3) RDATA(1)=ZA(na,nh,nc,ne)
                      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '                  FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,
     '                  IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RZERO,
     '                  -RMAX,RMAX,INFO,ERROR,*9999)
                      IF(IOTYPE.NE.3) ZA(na,nh,nc,ne)=RDATA(1)
                    ENDIF
                  ENDDO !na

                  IF(IOTYPE.NE.3) THEN
C                   Apply bdry conditions to rest of elements group
                    DO n=2,NELIST(0)
                      ne=NELIST(n)
                      NW(ne,1)=NW(NELIST(1),1)
                      DO na=1,NAT(NBH(nh,nc,ne))
                        ny=NYNE(na,nh,0,nc,ne)
                        ny_first=NYNE(na,nh,0,nc,NELIST(1))
                        ZA(na,nh,nc,ne)=ZA(na,nh,nc,NELIST(1))
                        FIX(ny,1)=FIX(ny_first,1)
                        FIX(ny,2)=FIX(ny_first,2)
                      ENDDO !na
                    ENDDO !n
                  ENDIF

                  GO TO 6720 !for more elements
                ENDIF !idata(1).ne.0
              ENDIF !nc=1 & nat(nb)>0

            ENDIF !nc=1 or nhx<=NJT
          ENDDO !nhx (nh)
        ENDDO !nc

        IF(ITYP1(nr,nx).EQ.4) THEN ! linear elasticity
          FORMAT='(/$,'' Are there any uniformly distributed loads '
     '      //'[NO]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='N'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,0,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
          IF((IOTYPE.NE.3).AND.(ADATA(1).EQ.'Y')) THEN
            nc=2 !force
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              IF(nc.EQ.1.OR.nh.LE.NJT) THEN
                WRITE(CHAR1,'(I1)') nhx
                FORMAT='(/'' Force component '//CHAR1(1:1)//' : '')'
                CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,0,IDATA,IDEFLT,
     '            IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
                nb=NBH(nh,nc,NEELEM(1,nr))
                IF(NNT(nb).GT.0) THEN
                  IF(IOTYPE.EQ.3) THEN
                    N1NODE=1
                  ENDIF

 7100             FORMAT='($,'' Enter node #s/name [EXIT]: '',I5)'
                  IF(IOTYPE.EQ.3) THEN
                    DO nonode=N1NODE,NPNODE(0,nr)
                      np=NPNODE(nonode,nr)
                      N2NODE=nonode
                      DO nv=1,NVHP(nh,np,nc)
                        DO nk=1,NKH(nh,np,nc)
                          IF(FIX(NYNP(nk,nv,nh,np,0,nc,nr),1)) GOTO 7200
                        ENDDO !nk
                      ENDDO !nv
                    ENDDO !nonode (np)
                    N2NODE=0
 7200               IF(N2NODE.EQ.0) THEN
                      IDATA(1)=0
                    ELSE
                      IDATA(1)=NP
                      N1NODE=N2NODE+1
                    ENDIF
                  ENDIF !iotype=3

 7500             CDATA(1)='NODES' !for use with group input
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,
     '              NOQUES,FILEIP,
     '              FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,
     '              IZERO,0,NPT(nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '              RMAX,INFO,ERROR,*9999)
                  IF(IDATA(1).NE.0) THEN !not default exit
                    IF(iotype.NE.3) THEN
                      NPLIST(0)=IDATA(0)
                      DO n=1,IDATA(0)
                        NPLIST(n)=IDATA(n)
                        np=IDATA(n)
                        IF(.NOT.INLIST(np,NPNODE(1,nr),
     '                    NPNODE(0,nr),N1)) THEN
                          WRITE(OP_STRING,'('' Node '',I5,'' does not '
     '                      //'belong to the current region'')') np
                          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                          GOTO 7500
                        ENDIF
                      ENDDO !n
                    ENDIF !iotype.NE.3

C               Define bdry condition for first node in group
                    np=NPLIST(1) !rest of group filled later
                    DO nv=1,NVHP(nh,np,nc)
                      IF(NVHP(nh,np,nc).GT.1) THEN
                        WRITE(CHAR1,'(I2)') nv
                        FORMAT='('' For version number '//CHAR1(1:2)
     '                    //':'')'
                        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,
     '                    FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
     '                    0,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '                    RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                      ENDIF
                      DO nk=1,NKH(nh,np,nc)
                        ny=NYNP(nk,nv,nh,np,0,nc,nr)
                        FIX(ny,1)=.TRUE. !ny=ny1v for nc=1, ny=ny2v for nc=2
C PJH 27May99       FIX(ny,2)=.TRUE.
                        IF(nk.EQ.1) THEN
                          FORMAT='($,'' The force is [0.0]: '',G12.5)'
                        ELSE
                          WRITE(FORMAT,'(A,I1,A)')
     '                      '($,'' The value of force derivative ',nk,
     '                      ' is [0.0]: '',G12.5)'
                        ENDIF
                        IF(IOTYPE.EQ.3) RDATA(1)=ZP(nk,nv,nh,np,nc)
                        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,
     '                    0,0,NOQUES,
     '                    FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '                    0,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '                    RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
                        IF(IOTYPE.NE.3) THEN
C                     Compute integrated force values
C                     Integrate over lines connected to node that
C                     are in only one element (ie external lines)
C                     AND and nodes on line are in group
                          DO n=1,NPLIST(0)
                            np=NPLIST(n)
                            ny=NYNP(nk,nv,nh,np,0,nc,nr)
                            FIX(ny,1)=.TRUE.
                            DO noelem=1,NENP(np,0,nr)
                              ne=NENP(np,noelem,nr)
                              IF(NIT(nb).EQ.2) THEN
C                             loop over local lines of current element
                                DO kl=1,NLE(nb)
                                  nl=NLL(kl,ne)
                                  IF(nl.NE.0) THEN
                                    nn=0
                                    GROUPED_LINE=.TRUE. !Check if in group
                                    DO WHILE ((nn.LT.NNL(0,kl,nb)).AND.
     '                                GROUPED_LINE) !loop over line nodes
                                      nn=nn+1
                                      IF(.NOT.INLIST(NPL(nn+1,1,nl),
     '                                  NPLIST(1),NPLIST(0),N1)) THEN
                                        GROUPED_LINE=.FALSE.
                                      ENDIF
                                    ENDDO !while
                                    IF(GROUPED_LINE) THEN
                                      nn=1
                                      DO WHILE(NPL(nn+1,1,nl).NE.np)
                                        nn=nn+1
                                      ENDDO
                                      nt=2
                                      CALL  CALC_CONTRIB_COEFF(COEFF,nb,
     '                                  nb,nt,NPL(1,0,nl),PG,WG,ERROR,
     '                                  *9999)
                                      fn=COEFF(nn)*DL(3,nl)*RDATA(1)
                                      ZP(nk,nv,nh,np,nc)=
     '                                  ZP(nk,nv,nh,np,nc)+fn
                                    ENDIF !GROUPED_LINE
                                  ENDIF !nl.NE.0
                                ENDDO !kl
                              ELSE IF(NIT(nb).EQ.3) THEN
C                           loop over local faces of current element
                                DO kf=1,NFE(nb)
                                  nf=NFF(kf,ne)
                                  CALL CALC_FACE_INFORMATION_DEP(
     '                              NBH(1,1,ne),NBHF(1,1,nf),kf,NHE(ne),
     '                              NKHE(1,1,1,ne),NKEF,NKHF,NNF,
     '                              NPNE(1,1,ne),NPNF,NVHE(1,1,1,ne),
     '                              NVHF,nx,SE(1,1,ne),SF,ERROR,*9999)
                                  nn=0
                                  GROUPED_FACE=.TRUE.
                                  DO WHILE ((nn.LT.NNF(0,kf,nb)).AND.
     '                              GROUPED_FACE)
                                    nn=nn+1
                                    IF(.NOT.INLIST(NPNE(NNF(nn+1,kf,nb),
     '                                nb,ne),
     '                                NPLIST(1),NPLIST(0),N1)) THEN
                                      GROUPED_FACE=.FALSE.
                                    ENDIF
                                  ENDDO !while
                                  IF(GROUPED_FACE) THEN
                                    nn=1
                                    DO WHILE(NPL(nn+1,1,nl).NE.np)
                                      nn=nn+1
                                    ENDDO
                                    CALL  CALC_CONTRIB_COEFF(COEFF,nb,
     '                                NBHF(NH_LOC(1,nx),1,nf),
     '                                NNF(0,kf,nb),NPL(1,0,1),PG,WG,
     '                                ERROR,*9999)
                                    fn=COEFF(nn)*DF(nf)*RDATA(1)
                                    ZP(nk,nv,nh,np,nc)=
     '                                ZP(nk,nv,nh,np,nc)+fn
                                  ENDIF !GROUPED_FACE
                                ENDDO !kf
                              ENDIF !NIT(nb)
                            ENDDO !noelem
                          ENDDO !n
                        ENDIF !IO_TYPE
                      ENDDO !nk
                    ENDDO !nv
                    GO TO 7100
                  ENDIF !idata(1).NE.0
                ENDIF !nnt(nb)>0
              ENDIF !nc.EQ.1.OR.nh.LE.NJT
            ENDDO !nh
          ENDIF !UDL
        ENDIF ! linear

      IF(loop.EQ.2) THEN ! boundary conditions
C cpb 23/3/96 Adding check for ability to translate or rotate the
C mesh and hence obtain a singular matrix. In order for the mesh not
C to translate displacement must be specified for at least one node
C in all directions. In order for the mesh not to rotate the
C displacement must be specified for at least two nodes for at least
C the number of directions minus 1.
        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
          NUMFIXED(nhx)=0
        ENDDO !nj
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
            nh=NH_LOC(nhx,nx)
            ny=NYNP(1,1,nh,np,0,1,nr)
            IF(FIX(ny,1)) NUMFIXED(nhx)=NUMFIXED(nhx)+1
          ENDDO !nhx
        ENDDO !nonode
        TRANSLATE=.FALSE.
        NUMGT1=0
        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
          nh=NH_LOC(nhx,nx)
          IF(NUMFIXED(nhx).EQ.0) TRANSLATE=.TRUE.
          IF(NUMFIXED(nh).GE.2) NUMGT1=NUMGT1+1
        ENDDO
C!!! Fixing derivatives aids in prevention of rotation too!
        ROTATE=NUMGT1.LT.(NH_LOC(0,nx)-1).AND..FALSE.
        IF(TRANSLATE) THEN
          WRITE(OP_STRING,'('' >>WARNING: The boundary conditions '
     '      //'permit rigid body translation'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(ROTATE) THEN
          WRITE(OP_STRING,'('' >>WARNING: The boundary conditions '
     '      //'permit rigid body rotation'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF !boundary conditions

        IF(IOTYPE.NE.3) THEN
          CALL ZPYP(iy,NBH,NEELEM,NHE,NHP,NKH,
     '      NPNODE,nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        ENDIF !IOTYPE!=3
      ENDDO !loop

      IF(ITYP1(nr,nx).EQ.5) THEN ! finite elasticity
         write(*,*) 'IOTYPE inside IPINIT_ELAS',IOTYPE
        IF(IOTYPE.EQ.3) THEN
          RDATA(1)=b(1)
          RDATA(2)=b(2)
          RDATA(3)=b(3)
          RDEFLT(1)=b(1)
          RDEFLT(2)=b(2)
          RDEFLT(3)=b(3)
        ELSE
          RDATA(1)=0
          RDATA(2)=0
          RDATA(3)=0
          RDEFLT(1)=0
          RDEFLT(2)=0
          RDEFLT(3)=0
C VJ 6Mar2005: Adding b_inc to contain the incremented body force values when using
C increment option in fem solve
          b_inc(1)=0.0d0
          b_inc(2)=0.0d0
          b_inc(3)=0.0d0  
        ENDIF
        FORMAT='($,'' Specify the gravity vector '
     '    //'components [0,0,0]: '',3G25.17)'
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '    FILEIP,FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,
     '    0,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          b(1)=RDATA(1)
          b(2)=RDATA(2)
          b(3)=RDATA(3)
        ENDIF
      ENDIF !finite

      write(*,*) 'RDATA and b'
      write(*,*) 'RDATA(1), RDATA(2),RDATA(3)'
      write(*,*) RDATA(1), RDATA(2),RDATA(3)
      write(*,*) 'b(1), b(2), b(3)'
      write(*,*) b(1), b(2), b(3)


C JWF 1/5/03 Inertia question - Defines KTYP5I(nr).
C Inertia term is for general finite elasticity 
C problems but is currently only optional in contact
C problems.       
      IF (KTYP5G(nr).GE.1) THEN
        FORMAT='($,'' Is this body inertia dependent [N]? '',A)'
        ADEFLT(1)='N'
        IF(IOTYPE.EQ.3) THEN
          IF (KTYP5I(nr).GT.0) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,0,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN 
          KTYP5I(nr)=1 ! set inertia option
          IF(IOTYPE.EQ.3) THEN  
            RDATA(1)=T_inc
            RDEFLT(1)=T_inc   
          ELSE 
            RDATA(1)=1.0d0
            RDEFLT(1)=1.0d0       
          ENDIF 
          FORMAT='($,'' Specify time increment '
     '    //'[1.0]?: '',D20.10)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '      FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '      0,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            T_inc=RDATA(1)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


