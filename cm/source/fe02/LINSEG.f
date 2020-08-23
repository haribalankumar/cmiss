      SUBROUTINE LINSEG(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '  NLLINE,NNL,NPL,NPNE,NPNODE,NVJE,NVJL,XP,ERROR,*)

C#### Subroutine: LINSEG
C###  Description:
C###    LINSEG defines the various parameters associated with
C###    global line segments nl=1,NLT, and element line segments
C###    NAE=1,NLE(nb).

C#### Variable: NVJL(nn,nj,nl)
C###  Type: INTEGER
C###  Set_up: LINSEG
C###  Description:
C###    NVJL(nn,nj,nl) is the global version number used for
C###    geometric variable nj on node nn of line nl.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NLLINE(0:NL_R_M,0:NRM),NNL(0:4,12,NBFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR,i1,i2,IBMAX,IBTYP,IBTYP1,IBTYP2,INPN,INPOS,JE,
     '  KOUNT,KOUNT1,l,MJ(5,2),mn,N,n1,n1elem,n2,NAE,
     '  nb,nbb,ne,nee,ni,ni1,ni2,ni3,NIA,nicollapse,NITB,nj,nj_nb,njj,
     '  njtype,nk,NKEL(2),nl,NLTOT,nn,nnn,NNEL(0:4),NNTB,noelem,nonee,
     '  nonode,np,np1,np2,NPLL(4),nr,NM(3),LIADJ,NUMCOLLAPSED
      LOGICAL COLLAPSED,COLLAPSEDNODE,CONT,CORNER,DERIV,FOUND,FOUND2,
     '  ISATCOLLAPSE,SECTOR

      DATA MJ/1,2,1,0,0,
     '        1,2,3,1,2/

      CALL ENTERS('LINSEG',*9999)

      NLT=0
      DO nr=1,NRT
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(/'' Region: '',I2)') nr
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        NLLINE(0,nr)=0 !initializes #lines in nr
C***    Loop over all nodes np
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' Node: '',I5)') np
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(/'' NENP(np,0..,nr): '',10I6)')
     '        (NENP(np,noelem,nr),noelem=0,NENP(np,0,nr))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
C***      Loop over the elements surrounding node np
          DO noelem=1,NENP(np,0,nr)
            ne=NENP(np,noelem,nr)
C ***       Find nb type for highest degree geometric interpolant
            IBMAX=1
            nj_nb=1
            DO nj=2,NJ_LOC(NJL_GEOM,0,nr)
              nb=NBJ(nj,ne)
              DO ni=1,NIT(nb)
                IF(IBT(1,ni,nb).EQ.1.AND.IBT(2,ni,nb).GT.IBMAX) THEN
                  IBMAX=IBT(2,ni,nb)
                  nj_nb=nj
                ELSE IF(IBT(1,ni,nb).EQ.2) THEN
                  IF(IBT(2,ni,nb).EQ.1) THEN
                    IBMAX=3
                    nj_nb=nj
                  ELSE IF(IBMAX.LT.2) THEN
                    IBMAX=2
                    nj_nb=nj
                  ENDIF
                ELSE IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
                  IF(IBT(2,ni,nb).EQ.4) THEN
                    IBMAX=2
                    nj_nb=nj
                  ELSE IF(IBT(2,ni,nb).GT.IBMAX) THEN
                    IBMAX=IBT(2,ni,nb)
                    nj_nb=nj
                  ENDIF
                ENDIF
              ENDDO !ni
            ENDDO !nj
            nb=NBJ(nj_nb,ne)
            nn=1
            DO WHILE(nn.LE.NNT(nb).AND.NPNE(nn,nb,ne).NE.NP)
              nn=nn+1
            ENDDO
            IF(nn.GT.NNT(nb)) THEN
              ERROR='>>Could not find local node'
              GOTO 9999
            ENDIF
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(/'' Element: '',I5,'', Local node:'','
     '          //'I3)') ne,nn
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' nb number of highest geometric'
     '          //' interpolant is '',I2)') nb
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            NITB=NIT(nb)
            NNTB=NNT(nb)
C***        Calculate whether node np is corner of current element ne
            CORNER=.TRUE.
            SECTOR=.FALSE.
            NUMCOLLAPSED=0
            DO ni=1,NITB
              IBTYP=IBT(2,ni,nb)
              INPOS=INP(nn,ni,nb)
              IF(IBT(1,ni,nb).EQ.1.OR.IBT(1,ni,nb).EQ.3) THEN
C               Lagrange or Simplex
                IF(IBTYP.EQ.2) THEN !quadratic
                  IF(INPOS.EQ.2) CORNER=.FALSE.
                ELSE IF(IBTYP.EQ.3) THEN !cubic
                  IF(INPOS.EQ.2.OR.INPOS.EQ.3) CORNER=.FALSE.
                ENDIF
              ELSE IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
C               Sector
                IF(IBTYP.LT.4) THEN !Lagrange
                  IF(IBTYP.EQ.2) THEN !quadratic
                    IF(INPOS.EQ.2) CORNER=.FALSE.
                  ELSE IF(IBTYP.EQ.3) THEN !cubic
                    IF(INPOS.EQ.2.OR.INPOS.EQ.3) CORNER=.FALSE.
                  ENDIF
                ENDIF
                SECTOR=.TRUE.
                NUMCOLLAPSED=NUMCOLLAPSED+1
              ENDIF
            ENDDO !ni
            IF(CORNER) THEN
              COLLAPSEDNODE=ISATCOLLAPSE(IBT(1,1,nb),INP(1,1,nb),nb,nn)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' SECTOR='',L1,'', COLLAPSEDNODE='','
     '            //'L1)') SECTOR,COLLAPSEDNODE
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              DO mn=1,NNTB
                IF(mn.NE.nn) THEN
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,'(/'' nn='',I3,'', mn='',I3)')
     '                nn,mn
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF
C***              Calculate whether MN is corner of element ne
                  CORNER=.TRUE.
                  DO ni=1,NITB
                    IBTYP=IBT(2,ni,nb)
                    INPOS=INP(mn,ni,nb)
                    IF(IBT(1,ni,nb).EQ.1.OR.IBT(1,ni,nb).EQ.3) THEN
C                     Lagrange or Simplex
                      IF(IBTYP.EQ.2) THEN !quadratic
                        IF(INPOS.EQ.2) CORNER=.FALSE.
                      ELSE IF(IBTYP.EQ.3) THEN !cubic
                        IF(INPOS.EQ.2.OR.INPOS.EQ.3) CORNER=.FALSE.
                      ENDIF
                    ELSE IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
C                     Sector
                      IF(IBTYP.LT.4) THEN !Lagrange
                        IF(IBTYP.EQ.2) THEN !quadratic
                          IF(INPOS.EQ.2) CORNER=.FALSE.
                        ELSE IF(IBTYP.EQ.3) THEN !cubic
                          IF(INPOS.EQ.2.OR.INPOS.EQ.3) CORNER=.FALSE.
                        ENDIF
                      ENDIF
                    ENDIF
                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'('' IBTYP='',I2,'' INPOS='',I2,'
     '                  //''' CORNER='',L1)') IBTYP,INPOS,CORNER
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF
                  ENDDO !ni
C***              If mn is a corner then try and find a line between
C***              mn and nn.
                  KOUNT=0
                  NIA=1
                  IF(CORNER) THEN
                    IF(COLLAPSEDNODE.AND.SECTOR) THEN
                      IF(NITB.EQ.2.OR.NUMCOLLAPSED.EQ.2) THEN
                        KOUNT=NITB-1
                        DO ni=1,NITB
                          IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6)
     '                      NIA=IBT(3,ni,nb)
                        ENDDO !ni
                      ELSE
                        DO ni=1,NITB
                          IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6)
     '                      THEN
                            ni2=IBT(3,ni,nb)
                            IF((ni.EQ.1.AND.ni2.EQ.2).OR.
     '                        (ni.EQ.2.AND.ni2.EQ.1)) THEN
                              ni3=3
                            ELSE IF((ni.EQ.1.AND.ni2.EQ.3).OR.
     '                          (ni.EQ.3.AND.ni2.EQ.1)) THEN
                              ni3=2
                            ELSE
                              ni3=1
                            ENDIF
                            IF(INP(nn,ni3,nb).EQ.INP(mn,ni3,nb)) THEN
                              KOUNT=NITB-1
                              NIA=IBT(3,ni,nb)
                            ELSE
                              IF(INP(nn,ni,nb).EQ.INP(mn,ni,nb).AND.
     '                          INP(nn,ni2,nb).EQ.INP(mn,ni2,nb)) THEN
                                KOUNT=NITB-1
                                NIA=ni3
                              ELSE
                                KOUNT=0
                                NIA=0
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDDO !ni
                      ENDIF
                    ELSE
                      DO ni=1,NITB
                        i1=INP(nn,ni,nb)
                        i2=INP(mn,ni,nb)
                        IF(i1.EQ.i2) THEN
C                       KOUNT is the number of 'lines' orthogonal to
C                       the proposed line between nn and mn
                          KOUNT=KOUNT+1
                        ELSE
                          NIA=ni !the xi direction of the line
                        ENDIF
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          WRITE(OP_STRING,
     '                      '('' I1='',I2,'' I2='',I2,'' KOUNT='','
     '                      //'I2,'' NIA='',I1)') i1,i2,KOUNT,NIA
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF
                      ENDDO !ni
                    ENDIF
                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'('' FINAL KOUNT='',I2,'' NIA='','
     '                  //'I2)') KOUNT,NIA
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF
                    IF(NIA.NE.0) THEN
                      IBTYP1=IBT(1,NIA,nb)
                      IBTYP2=IBT(2,NIA,nb)
                      INPOS=INP(nn,NIA,nb)
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' IBTYP1='',I2,'' IBTYP2='','
     '                    //'I2,'' INPOS='',I2)') IBTYP1,IBTYP2,INPOS
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF
                    ENDIF
                  ENDIF
                  IF(CORNER.AND.
     '              ((IBTYP1.NE.3.AND.KOUNT.EQ.(NIT(nb)-1))
     '              .OR.(IBTYP1.EQ.3.AND.KOUNT.LT.NIT(nb)))) THEN
C***                Find the nodes and derivatives along the line
                    IF(IBTYP1.EQ.1.OR.IBTYP1.EQ.3.OR.
     '                (IBTYP1.EQ.5.AND.IBTYP2.LE.3).OR.
     '                (IBTYP1.EQ.6.AND.IBTYP2.LE.3)) THEN
                      IF(IBTYP2.EQ.1) THEN
C                       Linear Lagrange
                        DERIV=.FALSE.
                        NNEL(0)=2
                        IF(INPOS.EQ.1) THEN
                          NNEL(1)=nn
                          NNEL(2)=mn
                        ELSE IF(INPOS.EQ.2) THEN
                          NNEL(1)=mn
                          NNEL(2)=nn
                        ENDIF
                      ELSE IF(IBTYP2.EQ.2) THEN
C                       Quadratic Lagrange
                        DERIV=.FALSE.
                        NNEL(0)=3
                        IF(INPOS.EQ.1) THEN
                          NNEL(1)=nn
                          NNEL(3)=mn
                        ELSE IF(INPOS.EQ.3) THEN
                          NNEL(1)=mn
                          NNEL(3)=nn
                        ENDIF
C                       Find midside node NP2
                        IF(NITB.EQ.1) THEN
C                         1D quadratic Lagrange
                          DO n=1,NNTB
                            IF(INP(n,NIA,nb).EQ.2) NNEL(2)=n
                          ENDDO !N
                        ELSE IF(NITB.GT.1.AND.IBTYP1.EQ.3) THEN
C                         quadratic simplex (note use of default INP)
                          IF(nn.EQ.1.AND.mn.EQ.3.OR
     '                      .nn.EQ.3.AND.mn.EQ.1) THEN
                            NNEL(2)=2
                          ELSE IF(nn.EQ.1.AND.mn.EQ.6.
     '                        or.nn.EQ.6.AND.mn.EQ.1) THEN
                            NNEL(2)=4
                          ELSE IF(nn.EQ.3.AND.mn.EQ.6.
     '                        or.nn.EQ.6.AND.mn.EQ.3) THEN
                            NNEL(2)=5
                          ENDIF
                        ELSE
C                         normal quadratic Lagrange
                          ni2=MJ(NIA+1,NITB-1)
                          ni3=MJ(NIA+2,NITB-1)
                          n=1
                          FOUND=.FALSE.
                          IF(COLLAPSEDNODE) THEN
                            nnn=mn
                          ELSE
                            nnn=nn
                          ENDIF
                          DO WHILE(n.LE.NNTB.AND..NOT.FOUND)
                            IF(INP(n,NIA,nb).EQ.2) THEN
                              IF(INP(n,ni2,nb).EQ.INP(nnn,ni2,nb)) THEN
! PJH 24May97 reorganise IF to avoid ni3=0 when NITB=2
                                IF(NITB.EQ.2) THEN
                                  NNEL(2)=n
                                  FOUND=.TRUE.
                                ELSE IF(NITB.EQ.3.AND.INP(n,ni3,nb)
     '                            .EQ.INP(nnn,ni3,nb)) THEN
                                  NNEL(2)=n
                                  FOUND=.TRUE.
                                ENDIF !nitb
                              ENDIF !inp
                            ENDIF !inp
                            n=n+1
                          ENDDO !n
                        ENDIF
                      ELSE IF(IBTYP2.EQ.3) THEN
C                       Cubic Lagrange
                        DERIV=.FALSE.
                        NNEL(0)=4
                        IF(INPOS.EQ.1) THEN
                          NNEL(1)=nn
                          NNEL(4)=mn
                        ELSE IF(INPOS.EQ.4) THEN
                          NNEL(1)=mn
                          NNEL(4)=nn
                        ENDIF
                        IF(COLLAPSEDNODE) THEN
                          nnn=mn
                        ELSE
                          nnn=nn
                        ENDIF
                        DO N=1,NNTB
                          INPN=INP(N,NIA,nb)
                          IF(INPN.EQ.2.OR.INPN.EQ.3) THEN
                            IF(NITB.EQ.1) THEN
                              NNEL(INPN)=N
                            ELSE IF(NITB.EQ.2) THEN
                              ni2=MJ(NIA+1,NITB-1)
                              IF(INP(N,ni2,nb).EQ.INP(nnn,ni2,nb)) THEN
                                NNEL(INPN)=N
                              ENDIF
                            ELSE IF(NITB.EQ.3) THEN
                              ni2=MJ(NIA+1,NITB-1)
                              ni3=MJ(NIA+2,NITB-1)
                              IF(INP(N,ni2,nb).EQ.INP(nnn,ni2,nb).AND.
     '                          (INP(N,ni3,nb).EQ.INP(nnn,ni3,nb))) THEN
                                NNEL(INPN)=N
                              ENDIF
                            ENDIF
                          ENDIF !INPN
                        ENDDO !n
                      ELSE IF((IBTYP1.EQ.3).AND.(IBTYP2.EQ.4))THEN
C                       Special Hermite simplex interpolation AJP 1-6-93
                        DERIV=.TRUE.
                        NNEL(0)=2
                        IF(NIA.EQ.1)THEN
C                         Always standard hermite in first direction
                          IF(INPOS.EQ.1) THEN
                            NNEL(1)=nn
                            NNEL(2)=mn
                          ELSE IF(INPOS.EQ.2) THEN
                            NNEL(1)=mn
                            NNEL(2)=nn
                          ENDIF
                          nk=1
                          FOUND=.FALSE.
                          DO WHILE(nk.LE.NKT(nn,nb).AND..NOT.FOUND)
                            IF(IDO(nk,nn,NIA,nb).EQ.2) THEN
                              NKEL(1)=nk
                              NKEL(2)=nk
                              FOUND=.TRUE.
                            ENDIF
                            nk=nk+1
                          ENDDO !nk
                        ELSE IF(NIA.EQ.2)THEN
C                         Special interpolation in this direction
C                         Note: only implemented for 2d at the moment
                          IF(INPOS.EQ.1) THEN
                            NNEL(1)=nn
                            NNEL(2)=mn
                          ELSE IF(INPOS.EQ.2) THEN
                            NNEL(1)=mn
                            NNEL(2)=nn
                          ENDIF
                          IF(NKT(1,nb).EQ.1)THEN !Apex at node 1
                            nk=1
                            FOUND=.FALSE.
                            DO WHILE(nk.LE.NKT(0,nb).AND..NOT.FOUND)
                              IF(IDO(nk,nn,NIA,nb).EQ.2) THEN
                                NKEL(1)=0
                                NKEL(2)=nk
                                FOUND=.TRUE.
                              ENDIF
                              nk=nk+1
                            ENDDO
                          ELSE IF(NKT(3,nb).EQ.1) THEN !Apex at node 3
                            nk=1
                            FOUND=.FALSE.
                            DO WHILE(nk.LE.NKT(0,nb).AND..NOT.FOUND)
                              IF(IDO(nk,nn,NIA,nb).EQ.2) THEN
                                NKEL(1)=nk
                                NKEL(2)=0
                              ENDIF
                              FOUND=.TRUE.
                              nk=nk+1
                            ENDDO
                          ENDIF
                        ENDIF !End of ni choice
                      ENDIF
                    ELSE IF(IBTYP1.EQ.2.OR.(IBTYP1.EQ.5.AND.IBTYP2
     '                .EQ.4).OR.(IBTYP1.EQ.6.AND.IBTYP2.EQ.4)) THEN
                      DERIV=.TRUE.
                      NNEL(0)=2
                      IF(INPOS.EQ.1) THEN
                        NNEL(1)=nn
                        NNEL(2)=mn
                      ELSE IF(INPOS.EQ.2) THEN
                        NNEL(1)=mn
                        NNEL(2)=nn
                      ENDIF
                      IF(IBTYP2.EQ.1.OR.IBTYP2.EQ.4) THEN
C                       Cubic Hermite or Hermite Sector in the collapsed
C                       xi direction
                        nk=1
                        FOUND=.FALSE.
                        DO WHILE(nk.LE.NKT(nn,nb).AND..NOT.FOUND)
                          IF(IDO(nk,nn,NIA,nb).EQ.2) THEN
                            NKEL(1)=nk
                            NKEL(2)=nk
                            FOUND=.TRUE.
                          ENDIF
                          nk=nk+1
                        ENDDO
                      ELSE IF(IBTYP2.EQ.2.OR.IBTYP2.EQ.3) THEN
C                       Quadratic Hermite
                        IF(COLLAPSEDNODE) THEN
                          nnn=mn
                        ELSE
                          nnn=nn
                        ENDIF
                        IF(IBTYP2.EQ.2) THEN
                          NKEL(1)=0
                          nk=1
                          FOUND=.FALSE.
                          DO WHILE(nk.LE.NKT(nnn,nb).AND..NOT.FOUND)
                            IF(IDO(nk,nnn,NIA,nb).EQ.2) THEN
                              NKEL(2)=nk
                              FOUND=.TRUE.
                            ENDIF
                            nk=nk+1
                          ENDDO !nk
                        ELSE IF(IBTYP2.EQ.3) THEN
                          NKEL(2)=0
                          nk=1
                          FOUND=.FALSE.
                          DO WHILE(nk.LE.NKT(nnn,nb).AND..NOT.FOUND)
                            IF(IDO(nk,nnn,NIA,nb).EQ.2) THEN
                              NKEL(1)=nk
                              FOUND=.TRUE.
                            ENDIF
                            nk=nk+1
                          ENDDO !nk
                        ENDIF
                      ENDIF
                    ENDIF
C                   Calculate values for NPL from local node and deriv,
C                   and see if the line is collapsed.
                    IF(DERIV) THEN
                      DO n=1,2
                        IF(NKEL(n).EQ.0) THEN
                          NPLL(2+n)=0
                        ELSE
                          NPLL(2+n)=NKJE(NKEL(n),NNEL(n),nj_nb,ne)
                        ENDIF
                      ENDDO
                    ELSE
                      DO n=NNEL(0)+1,4
                        NPLL(n)=0
                      ENDDO
                    ENDIF !DERIV
                    COLLAPSED=.TRUE.
                    DO n=1,NNEL(0)
                      np2=NPNE(NNEL(n),nb,ne)
                      NPLL(n)=np2
                      IF(n.EQ.1) THEN
                        np1=np2
                      ELSE IF(COLLAPSED) THEN
                        IF(np2.NE.np1) THEN
                          COLLAPSED=.FALSE.
                        ELSE IF(ITYP10(nr).GE.2.AND
     '                    .XP(1,1,NJ_LOC(NJL_GEOM,1,nr),np2).NE.0.0d0)
     '                    THEN !not rect. Cart.
C                         Only collapsed if r=0 (cylindrical or spherical),
C                         or lambda=0 or mu=0 (prolate or oblate).
                          IF(ITYP10(nr).LE.3) THEN !cyl. or sph.
                            COLLAPSED=.FALSE.
                          ELSE IF(XP(1,1,NJ_LOC(NJL_GEOM,2,nr),np2)
     '                      .NE.0.0d0) THEN
                            COLLAPSED=.FALSE.
                          ENDIF
                        ENDIF
                      ENDIF !COLLAPSED
                    ENDDO
                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'('' NP1='',I5,'' NP2='',I5,'
     '                  //''' NP3='',I5,'' NP4='',I5)')
     '                  NPLL(1),NPLL(2),NPLL(3),NPLL(4)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF

                    IF(COLLAPSED) THEN
                      FOUND=.TRUE. !to avoid following sections

                    ELSE
C***                  Check to see if the line has already been defined
                      IF(NIT(nb).EQ.1) THEN
                        FOUND=.FALSE.
                        nl=1
                        DO WHILE(nl.LE.NLT.AND.nl.LE.NLM.AND..NOT.FOUND)
                          IF(NEL(1,nl).EQ.ne) THEN
                            FOUND=.TRUE.
                          ELSE
                            nl=nl+1
                          ENDIF
                        ENDDO !nl
                      ELSE
C!!! KAT 2001-12-19 Slow
                        FOUND=.FALSE.
                        nl=0
                        NLTOT=MIN(NLT,NLM)
                        DO WHILE(nl.LT.NLTOT.AND..NOT.FOUND)
                          nl=nl+1
C!!! KAT 2001-12-19 This only works for similar interpolation for each nj
                          IF(NPL(2,1,nl).EQ.NPLL(1)
     '                      .AND.NPL(3,1,nl).EQ.NPLL(2)
     '                      .AND.NPL(4,1,nl).EQ.NPLL(3)
     '                      .AND.NPL(5,1,nl).EQ.NPLL(4)) THEN
                            FOUND=.TRUE.
                          ENDIF
C KAT 2001-12-19:         This could be used to check for match in
C                         either direction but the meaning of NPL(2:3,0,nl)
C                         would have to be considered.
C                          IF(DERIV) THEN
C                            IF(NPL(1,1,nl).GT.3) THEN !hermite
C                              IF((NPL(2,1,nl).EQ.NPLL(1)
C     '                          .AND.NPL(3,1,nl).EQ.NPLL(2)
C     '                          .AND.NPL(4,1,nl).EQ.NPLL(3)
C     '                          .AND.NPL(5,1,nl).EQ.NPLL(4))
C     '                          .OR.(NPL(2,1,nl).EQ.NPLL(2)
C     '                          .AND.NPL(3,1,nl).EQ.NPLL(1)
C     '                          .AND.NPL(4,1,nl).EQ.NPLL(4)
C     '                          .AND.NPL(5,1,nl).EQ.NPLL(3))) THEN
C                                FOUND=.TRUE.
C                              ENDIF
C                            ENDIF !hermite
C                          ELSE !not DERIV
C                            IF(NPL(1,1,nl).EQ.1) THEN !linear
C                              IF((NPL(2,1,nl).EQ.NPLL(1)
C     '                          .AND.NPL(3,1,nl).EQ.NPLL(2)
C     '                          .AND.NPL(4,1,nl).EQ.NPLL(3)
C     '                          .AND.NPL(5,1,nl).EQ.NPLL(4))
C     '                          .OR.(NPL(2,1,nl).EQ.NPLL(2)
C     '                          .AND.NPL(3,1,nl).EQ.NPLL(1)
C     '                          .AND.0.EQ.NPLL(4)
C     '                          .AND.0.EQ.NPLL(3))) THEN
C                                FOUND=.TRUE.
C                              ENDIF
C                            ELSEIF(NPL(1,1,nl).EQ.1) THEN !quadratic
C                              IF((NPL(2,1,nl).EQ.NPLL(1)
C     '                          .AND.NPL(3,1,nl).EQ.NPLL(2)
C     '                          .AND.NPL(4,1,nl).EQ.NPLL(3)
C     '                          .AND.NPL(5,1,nl).EQ.NPLL(4))
C     '                          .OR.(NPL(2,1,nl).EQ.NPLL(3)
C     '                          .AND.NPL(3,1,nl).EQ.NPLL(2)
C     '                          .AND.NPL(4,1,nl).EQ.NPLL(1)
C     '                          .AND.0.EQ.NPLL(4))) THEN
C                                FOUND=.TRUE.
C                              ENDIF
C                            ELSEIF(NPL(1,1,nl).EQ.1) THEN !cubic
C                              IF((NPL(2,1,nl).EQ.NPLL(1)
C     '                          .AND.NPL(3,1,nl).EQ.NPLL(2)
C     '                          .AND.NPL(4,1,nl).EQ.NPLL(3)
C     '                          .AND.NPL(5,1,nl).EQ.NPLL(4))
C     '                          .OR.(NPL(2,1,nl).EQ.NPLL(4)
C     '                          .AND.NPL(3,1,nl).EQ.NPLL(3)
C     '                          .AND.NPL(4,1,nl).EQ.NPLL(2)
C     '                          .AND.NPL(5,1,nl).EQ.NPLL(1))) THEN
C                                FOUND=.TRUE.
C                              ENDIF
C                            ENDIF !hermite
C                          ENDIF !DERIV
                        ENDDO !nl
                      ENDIF
C MLB 16 Sept 1997
C Check to see if the line is defined in the current region -
C needed for FEM/BEM coupling
                      IF(FOUND) THEN
                        FOUND2=.FALSE.

C LKC 7-NOV-97 added assert
                        CALL ASSERT(NL_R_M.GE.NLLINE(0,nr),
     '                    '>>Increase NL_R_M',ERROR,*9999)
C!!! KAT 2001-12-19 Slow
                        DO ni1=1,NLLINE(0,nr)
                          IF(NLLINE(ni1,nr).EQ.nl) FOUND2=.TRUE.
                        ENDDO
                        IF(.NOT.FOUND2) THEN
                          NLLINE(0,nr)=NLLINE(0,nr)+1
                          IF(NLLINE(0,nr).LE.NL_R_M) THEN
                            NLLINE(NLLINE(0,nr),nr)=nl !rec global line#
                          ENDIF
                          CALL ASSERT(NEL(0,nl)+1.LE.NELM,
     '                      '>>NELM too small',ERROR,*9999)
                          NEL(0,nl)=NEL(0,nl)+1
                          NEL(NEL(0,nl),nl)=ne
                        ENDIF
                      ENDIF
C end of check
                    ENDIF !COLLAPSED

                    IF(.NOT.FOUND) THEN
                      NLT=NLT+1
                      IF(NLT.LE.NLM) THEN
C***                    Set up NPL
                        DO i1=1,5
                          DO nj=0,3
                            NPL(i1,nj,NLT)=0
                          ENDDO
                        ENDDO
                        NPL(1,0,NLT)=NIA
                        DO n=1,4
                          NPL(1+n,1,NLT)=NPLL(n)
                        ENDDO
C                       Versions
C                       Fibre/Field version are done here also as no
C                       other routine has been written yet.
                        DO njtype=1,3
                          DO njj=1,NJ_LOC(njtype,0,nr)
                            nj=NJ_LOC(njtype,njj,nr)
                            nbb=NBJ(nj,ne)
                            IF(nbb.NE.0) THEN
                              DO n=1,NNEL(0)
                                NVJL(n,nj,NLT)=NVJE(NNEL(n),nbb,nj,ne)
                              ENDDO !n
                            ENDIF !nbb
                          ENDDO !njj
                        ENDDO !njtype
                        NLLINE(0,nr)=NLLINE(0,nr)+1 !incr #lines in nr
                        IF(NLLINE(0,nr).LE.NL_R_M) THEN
                          NLLINE(NLLINE(0,nr),nr)=NLT !rec global line#
                        ENDIF
                        DO nj=1,NJT
                          nbb=NBJ(nj,ne)
                          IF(IBT(1,NIA,nbb).EQ.1) THEN !Lagrange
                            NPL(1,nj,NLT)=IBT(2,NIA,nbb)
                          ELSE IF(IBT(1,NIA,nbb).EQ.2) THEN !Hermite
                            IF(IBT(2,NIA,nbb).EQ.1) THEN !cubic
                              NPL(1,nj,NLT)=4
                            ELSE IF(IBT(2,NIA,nbb).EQ.2) THEN !quad 1
                              NPL(1,nj,NLT)=6
                            ELSE IF(IBT(2,NIA,nbb).EQ.3) THEN !quad 2
                              NPL(1,nj,NLT)=7
                            ENDIF
                          ELSE IF(IBT(1,NIA,nbb).EQ.3) THEN !Simplex
                            IF(IBT(2,NIA,nbb).EQ.1) THEN !lin Lagrange
                              NPL(1,nj,NLT)=1
                            ELSE IF(IBT(2,NIA,nbb).EQ.2) THEN !quad L.
                              NPL(1,nj,NLT)=2
                            ELSE
                              IF(NIA.EQ.1) THEN !Hermite
                                NPL(1,nj,NLT)=4
                              ELSE IF(NIA.EQ.2)THEN !Special hermite
                                IF(NKT(1,nb).EQ.1)THEN !Apex node 1
                                  NPL(1,nj,NLT)=6
                                ELSE IF(NKT(3,nb).EQ.1)THEN !Apex node 3
                                  NPL(1,nj,NLT)=7
                                ENDIF
                              ELSE
                                NPL(1,nj,NLT)=IBT(2,NIA,nbb)
                              ENDIF
                            ENDIF
                          ELSE IF(IBT(1,NIA,nbb).EQ.5.OR.
     '                      IBT(1,NIA,nbb).EQ.6) THEN !Sector
                            NPL(1,nj,NLT)=IBT(2,NIA,nbb)
                          ELSE
                            NPL(1,nj,NLT)=1
                          ENDIF
                        ENDDO !nj
C***                    Set up NEL
                        NEL(0,NLT)=1
                        NEL(1,NLT)=ne
                        IF(NIT(nb).GT.1) THEN
                          JE=1
                          DO nonee=1,NENP(np,0,nr)
                            nee=NENP(np,nonee,nr)
                            IF(nee.NE.ne) THEN
                              nbb=NBJ(1,nee)
                              FOUND=.FALSE.
                              nnn=1
                              DO WHILE(nnn.LE.NNT(nbb).AND..NOT.FOUND)
                                IF((NPNE(nnn,nbb,nee).NE.NPNE(nn,nb,ne))
     '                            .AND.(NPNE(nnn,nbb,nee).EQ.
     '                            NPNE(mn,nb,ne))) THEN
                                  JE=JE+1
                                  NEL(0,NLT)=NEL(0,NLT)+1
                                  IF(JE.LE.NELM) THEN
                                    NEL(JE,NLT)=nee
                                  ENDIF
                                  FOUND=.TRUE.
                                ELSE
                                  nnn=nnn+1
                                ENDIF
                              ENDDO !nn
                            ENDIF
                          ENDDO !nonee (nee)
                        ENDIF
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          DO nj=0,NJ_LOC(NJL_GEOM,0,nr)
                            WRITE(OP_STRING,'('' NPL(1..5,nj='',I1,''
     '                        //'','',I5,''):'','
     '                        //'10I5,/,18X,10I5)') nj,NLT,
     '                        (NPL(L,nj,NLT),L=1,5)
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDDO
                          WRITE(OP_STRING,'('' NEL(0,'',I5,''):'','
     '                      //'I5)') NLT,NEL(0,NLT)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          WRITE(OP_STRING,'('' NEL(1..,'',I5,''):'','
     '                      //'20I5)') NLT,(NEL(L,NLT),L=1,NEL(0,NLT))
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF !DOP
                      ENDIF !NLT<NLM
                    ENDIF !new line
                  ENDIF
                ENDIF
              ENDDO !mn
            ENDIF !corner
          ENDDO !noelem
        ENDDO !nonode (np)
        CALL ASSERT(NLLINE(0,nr).LE.NL_R_M,'>>Increase NL_R_M',ERROR,
     '    *9999)
      ENDDO !nr
      IF(NLT.GT.NLM) THEN
        CALL FLAG_ERROR(0,' ')
        CALL WRITE_CHAR(IOER,'Increase NLM.  Try ',ERR)
        CALL WRITE_INT(IOER,NLT,ERR)
        CALL WRITE_CHAR(IOER,NEWLINE,ERR)
        GOTO 9998
      ENDIF

      DO nb=1,NBFT
C GMH 19-4-96 Ensure NNL is initialised
        DO ni1=0,4
          DO ni2=1,12
            NNL(ni1,ni2,nb)=0
          ENDDO !ni2
        ENDDO !ni1
        NITB=NIT(nb)
        IF(IBT(1,NITB,nb).EQ.9) THEN !time tensor product basis
          NITB=NITB-1                !spatial dimension
        ENDIF
C***    Define the NNL array (element node #s) & NLE array (tot # edges)
        SECTOR=.FALSE.
        DO ni=1,NITB
          IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) SECTOR=.TRUE.
        ENDDO !ni
        IF((IBT(1,1,nb).EQ.1.OR.IBT(1,1,nb).EQ.2).AND..NOT.SECTOR) THEN

C 21/2/97 LC archived section :  Lagrange/Hermite
C                   cpb 7/8/95 Rewritting this code

          IF(NITB.EQ.1) THEN
            DO n1=1,IBT(2,1,nb)+1
              DO nn=1,NNT(nb)
                IF(INP(nn,1,nb).EQ.n1) THEN
                  NNL(n1,1,nb)=nn
                  NNL(0,1,nb)=NNL(0,1,nb)+1
                ENDIF
              ENDDO !nn
            ENDDO !n1
            NLE(nb)=1
          ELSE IF(NITB.EQ.2) THEN
C           First find the max extents of this basis function
            DO ni=1,2
              NM(ni)=0
              DO nn=1,NNT(nb)
                IF(INP(nn,ni,nb).GT.NM(ni)) NM(ni)=INP(nn,ni,nb)
              ENDDO !nn
            ENDDO !ni
            NAE=1
            DO ni1=1,2
              ni2=MJ(ni1+1,1)
C             We are looking for lines going in the NI1 direction,
C             From the starting point NI1=0.
              DO n1=1,NM(ni2),NM(ni2)-1
                KOUNT=0
                DO nn=1,NNT(nb)
                  IF(INP(nn,ni2,nb).EQ.n1) THEN
                    KOUNT=KOUNT+1
                    NNL(KOUNT,NAE,nb)=nn
                    NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                  ENDIF
                ENDDO !nn
                IF(KOUNT.GT.1) THEN
                  IF(KOUNT.LT.4) NNL(KOUNT+1,NAE,nb)=0
                  NAE=NAE+1
                ENDIF
              ENDDO !n1
            ENDDO !ni1
            NLE(nb)=NAE-1
          ELSE
C           First find the max extents of this basis function
            DO ni=1,3
              NM(ni)=0
              DO nn=1,NNT(nb)
                IF(INP(nn,ni,nb).GT.NM(ni)) NM(ni)=INP(nn,ni,nb)
              ENDDO !nn
            ENDDO !ni
            NAE=1
            DO ni1=1,3
              ni2=MJ(ni1+1,2)
              ni3=MJ(ni1+2,2)
C             We are looking for lines going in the NI1 direction,
C             From the starting point NI1=0.
              DO n2=1,NM(ni3),NM(ni3)-1
                DO n1=1,NM(ni2),NM(ni2)-1
                  KOUNT=0
                  DO nn=1,NNT(nb)
                    IF(INP(nn,ni2,nb).EQ.n1.AND.
     '                INP(nn,ni3,nb).EQ.n2) THEN
                      KOUNT=KOUNT+1
                      NNL(KOUNT,NAE,nb)=nn
                      NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                    ENDIF
                  ENDDO !nn
                  IF(KOUNT.GT.1) THEN
                    IF(KOUNT.LT.4) NNL(KOUNT+1,NAE,nb)=0
                    NAE=NAE+1
                  ENDIF
                ENDDO !n1
              ENDDO !n2
            ENDDO !ni1
            NLE(nb)=NAE-1
          ENDIF
        ELSE IF(IBT(1,1,nb).EQ.3) THEN !Simplex basis 2D
          IF(NIT(nb).EQ.2) THEN
            NNL(0,1,nb)=2
            NNL(1,1,nb)=1
            NNL(2,1,nb)=2
            NNL(0,2,nb)=2
            NNL(1,2,nb)=1
            NNL(2,2,nb)=3
            NNL(0,3,nb)=2
            NNL(1,3,nb)=2
            NNL(2,3,nb)=3
            NLE(nb)=3
          ELSEIF(NIT(nb).EQ.3) THEN !Simplex basis 3D
            NNL(0,1,nb) = 2
            NNL(1,1,nb) = 1
            NNL(2,1,nb) = 2
            NNL(0,2,nb) = 2
            NNL(1,2,nb) = 1
            NNL(2,2,nb) = 3
            NNL(0,3,nb) = 2
            NNL(1,3,nb) = 1
            NNL(2,3,nb) = 4
            NNL(0,4,nb) = 2
            NNL(1,4,nb) = 2
            NNL(2,4,nb) = 3
            NNL(0,5,nb) = 2
            NNL(1,5,nb) = 2
            NNL(2,5,nb) = 4
            NNL(0,6,nb) = 2
            NNL(1,6,nb) = 3
            NNL(2,6,nb) = 4
            NLE(nb) = 6 ! Lines
          ENDIF
        ELSE IF(IBT(1,1,nb).EQ.4) THEN
        ELSE IF(SECTOR) THEN !Sector
c cpb 7/8/95 Generalising sector elements
          IF(NITB.EQ.2) THEN
C           First find the max extents of this basis function
            DO ni=1,2
              NM(ni)=0
              DO nn=1,NNT(nb)
                IF(INP(nn,ni,nb).GT.NM(ni)) NM(ni)=INP(nn,ni,nb)
              ENDDO !nn
            ENDDO !ni
            NAE=1
            DO ni1=1,2
              ni2=MJ(ni1+1,1)
C             We are looking for lines going in the NI1 direction,
C             From the starting point NI1=0.
              DO n1=1,NM(ni2),NM(ni2)-1
                KOUNT=0
                KOUNT1=0
                DO nn=1,NNT(nb)
                  IF(IBT(1,ni2,nb).EQ.5) THEN
                    IF(INP(nn,ni2,nb).EQ.n1.OR.INP(nn,ni1,nb).EQ.1) THEN
                      KOUNT=KOUNT+1
                      NNL(KOUNT,NAE,nb)=nn
                      NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                    ENDIF
                  ELSE IF(IBT(1,ni2,nb).EQ.6) THEN
                    IF(INP(nn,ni2,nb).EQ.n1) THEN
                      KOUNT=KOUNT+1
                      NNL(KOUNT,NAE,nb)=nn
                      NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                    ELSE IF(INP(nn,ni1,nb).EQ.NM(ni1)) THEN
                      IF(ni2.EQ.2) THEN
                        KOUNT1=KOUNT1+1
                        NNL(NM(ni1),NAE,nb)=nn
                        NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                      ELSE
                        KOUNT=KOUNT+1
                        NNL(KOUNT,NAE,nb)=nn
                        NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                      ENDIF
                    ENDIF
                  ELSE IF(INP(nn,ni2,nb).EQ.n1) THEN
                    KOUNT=KOUNT+1
                    NNL(KOUNT,NAE,nb)=nn
                    NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                  ENDIF
                ENDDO !nn
                IF((KOUNT+KOUNT1).GT.1) THEN
                  IF((KOUNT+KOUNT1).LT.4) NNL(KOUNT+KOUNT1+1,NAE,nb)=0
                  NAE=NAE+1
                ENDIF
              ENDDO !n1
            ENDDO !ni1
            NLE(nb)=NAE-1
          ELSE
C           First find the max extents of this basis function
            DO ni=1,3
              NM(ni)=0
              DO nn=1,NNT(nb)
                IF(INP(nn,ni,nb).GT.NM(ni)) NM(ni)=INP(nn,ni,nb)
              ENDDO !nn
            ENDDO !ni
            NAE=1
            DO ni1=1,3
              ni2=MJ(ni1+1,2)
              ni3=MJ(ni1+2,2)
C             We are looking for lines going in the NI1 direction,
C             From the starting point NI1=0.
              NUMCOLLAPSED=0
              IF(IBT(1,ni2,nb).EQ.5.OR.IBT(1,ni2,nb).EQ.6) THEN
                NUMCOLLAPSED=1
                nicollapse=IBT(3,ni2,nb)
              ENDIF
              IF(IBT(1,ni3,nb).EQ.5.OR.IBT(1,ni3,nb).EQ.6) THEN
                NUMCOLLAPSED=NUMCOLLAPSED+1
                nicollapse=IBT(3,ni3,nb)
              ENDIF
              DO n2=1,NM(ni3),NM(ni3)-1
                DO n1=1,NM(ni2),NM(ni2)-1
                  KOUNT=0
                  KOUNT1=0
                  DO nn=1,NNT(nb)
                    IF(NUMCOLLAPSED.EQ.1) THEN
                      IF(IBT(1,ni2,nb).EQ.5) THEN
                        IF((INP(nn,ni2,nb).EQ.n1.OR.
     '                    (INP(nn,nicollapse,nb).EQ.1.AND.
     '                    nicollapse.EQ.ni1)).AND.
     '                    INP(nn,ni3,nb).EQ.n2) THEN
                          KOUNT=KOUNT+1
                          NNL(KOUNT,NAE,nb)=nn
                          NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                        ENDIF
                      ELSE IF(IBT(1,ni2,nb).EQ.6) THEN
                        IF(INP(nn,ni2,nb).EQ.n1.AND.
     '                    INP(nn,ni3,nb).EQ.n2) THEN
                          KOUNT=KOUNT+1
                          NNL(KOUNT,NAE,nb)=nn
                          NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                        ELSE IF((INP(nn,nicollapse,nb).EQ.NM(ni1).AND.
     '                    nicollapse.EQ.ni1).AND.
     '                    INP(nn,ni3,nb).EQ.n2) THEN
                          IF(nicollapse.LT.ni2) THEN
                            KOUNT1=KOUNT1+1
                            NNL(NM(ni1),NAE,nb)=nn
                            NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                          ELSE
                            KOUNT=KOUNT+1
                            NNL(KOUNT,NAE,nb)=nn
                            NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                          ENDIF
                        ENDIF
                      ELSE IF(IBT(1,ni3,nb).EQ.5) THEN
                        IF((INP(nn,ni3,nb).EQ.n2.OR.
     '                    (INP(nn,nicollapse,nb).EQ.1.AND.
     '                    nicollapse.EQ.ni1)).AND.
     '                    INP(nn,ni2,nb).EQ.n1) THEN
                          KOUNT=KOUNT+1
                          NNL(KOUNT,NAE,nb)=nn
                          NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                        ENDIF
                      ELSE IF(IBT(1,ni3,nb).EQ.6) THEN
                        IF(INP(nn,ni3,nb).EQ.n2.AND.
     '                    INP(nn,ni2,nb).EQ.n1) THEN
                          KOUNT=KOUNT+1
                          NNL(KOUNT,NAE,nb)=nn
                          NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                        ELSE IF((INP(nn,nicollapse,nb).EQ.NM(ni1).AND.
     '                    nicollapse.EQ.ni1).AND.
     '                    INP(nn,ni2,nb).EQ.n1) THEN
                          IF(nicollapse.LT.ni3) THEN
                            KOUNT1=KOUNT1+1
                            NNL(NM(ni1),NAE,nb)=nn
                            NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                          ELSE
                            KOUNT=KOUNT+1
                            NNL(KOUNT,NAE,nb)=nn
                            NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                          ENDIF
                        ENDIF
                      ENDIF
                    ELSE IF(NUMCOLLAPSED.EQ.2) THEN
                      IF(IBT(1,ni2,nb).EQ.5) THEN
                        IF((INP(nn,ni2,nb).EQ.n1.OR.
     '                    INP(nn,nicollapse,nb).EQ.1).AND.
     '                    (INP(nn,ni3,nb).EQ.n2.OR.
     '                    INP(nn,nicollapse,nb).EQ.1))
     '                    THEN
                          KOUNT=KOUNT+1
                          NNL(KOUNT,NAE,nb)=nn
                          NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                        ENDIF
                      ELSE
                        IF(INP(nn,ni2,nb).EQ.n1.AND.
     '                    INP(nn,ni3,nb).EQ.n2) THEN
                          KOUNT=KOUNT+1
                          NNL(KOUNT,NAE,nb)=nn
                          NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                        ELSE IF(INP(nn,nicollapse,nb).EQ.NM(ni1).AND.
     '                    INP(nn,nicollapse,nb).EQ.NM(ni1)) THEN
                          IF(nicollapse.LT.3) THEN
                            KOUNT1=KOUNT1+1
                            NNL(NM(ni1),NAE,nb)=nn
                            NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                          ELSE
                            KOUNT=KOUNT+1
                            NNL(KOUNT,NAE,nb)=nn
                            NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                          ENDIF
                        ENDIF
                      ENDIF
                    ELSE
                      IF(INP(nn,ni2,nb).EQ.n1.AND.
     '                  INP(nn,ni3,nb).EQ.n2) THEN
                        KOUNT=KOUNT+1
                        NNL(KOUNT,NAE,nb)=nn
                        NNL(0,NAE,nb)=NNL(0,NAE,nb)+1
                      ENDIF
                    ENDIF
                  ENDDO !nn
                  IF((KOUNT+KOUNT1).GT.1) THEN
                    IF((KOUNT+KOUNT1).LT.4) NNL(KOUNT+KOUNT1+1,NAE,nb)=0
                    NAE=NAE+1
                  ENDIF
                ENDDO !n1
              ENDDO !n2
            ENDDO !ni1
            NLE(nb)=NAE-1
          ENDIF
        ENDIF
      ENDDO !nb

      DO nr=1,NRT
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          DO NAE=1,12
            NLL(NAE,ne)=0
          ENDDO
          DO NAE=1,NLE(nb)
            DO nl=1,NLT
C             Check whether all global nodes of current element edge
C             NAE of current element ne match global nodes of current
C             global line nl
              CONT=.TRUE.
              DO n1=1,4
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  DO nj=0,NJ_LOC(NJL_GEOM,0,nr)
                    WRITE(OP_STRING,'(/'' ne='',I5,'' NAE='','
     '                //'I2,'' nl='',I5,'' N1='',I1,'' NNL='','
     '                //'I5,'' nj='',I1,'' NPL(N1+1,nj,nl)='','
     '                //'I5)') ne,NAE,nl, n1,NNL(n1,NAE,nb),
     '                nj,NPL(n1+1,nj,nl)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDDO
CC$                call mp_unsetlock()
                ENDIF
                IF(NNL(n1,NAE,nb).GT.0) THEN
                  IF((NPNE(NNL(n1,NAE,nb),nb,ne).NE.NPL(n1+1,1,nl)))
     '              THEN
C KAT 27Oct98: This doesn't work properly and its not needed at the
C              moment anyway.
CC New CS 15/11/97 Also check in opposite direction if non-standard
CC line mappings are not used.
C                    IF(JTYP2B.NE.1.AND..NOT.SECTOR) THEN
C                      IF(NPNE(NNL(n1,NAE,nb),nb,ne).NE.
C     '                  NPL(NNL(0,NAE,nb)+2-n1,1,nl)) THEN
C                        CONT=.FALSE.
C                      ENDIF
C                    ELSE
                    CONT=.FALSE.
C                    ENDIF
                  ENDIF
                ENDIF !node exixts
              ENDDO !n1
              IF(CONT) THEN
                IF(NIT(nb).EQ.1) THEN
                  n1elem=1
                  FOUND=.FALSE.
                  DO WHILE(n1elem.LE.(noelem-1).AND..NOT.FOUND)
                    nee=NEELEM(n1elem,nr)
                    IF(nl.EQ.NLL(1,nee)) THEN
                      FOUND=.TRUE.
                    ELSE
                      n1elem=n1elem+1
                    ENDIF
                  ENDDO
                  IF(.NOT.FOUND) NLL(NAE,ne)=nl
                ELSE
                  NLL(NAE,ne)=nl
                ENDIF
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,
     '              '(/'' For element ne='',I4,'' local arc NAE='',I2,'
     '              //''' is global arc NLL(NAE,ne)='',I4)') ne,NAE,nl
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                GOTO 90
              ENDIF
            ENDDO !nl
 90         CONTINUE
          ENDDO !NAE
        ENDDO !noelem
      ENDDO !nr

      DO nl=1,NLT
        CALL ASSERT(NEL(0,nl).LE.NELM,'>>NEL cannot hold all elements:'
     '    //' increase NELM',ERROR,*9999)
      ENDDO

      DO nl=1,NLT
        DO nn=1,2
C!!! KAT 2001-12-19 Slow
          NPL(1+nn,0,nl)=LIADJ(nn,nl,NPL,NVJL)
        ENDDO
      ENDDO

      CALL EXITS('LINSEG')
      RETURN
 9998 ERROR=' '
 9999 CALL ERRORS('LINSEG',ERROR)
      CALL EXITS('LINSEG')
      RETURN 1
      END


