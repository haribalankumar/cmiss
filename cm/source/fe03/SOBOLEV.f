      SUBROUTINE SOBOLEV(IBT,IDO,INP,NBJ,NEELEM,NKJE,NLL,NONL,
     '  NONY,NPF,NPNE,nr,NRE,NVJE,NYNO,NYNP,
     '  DL,JACOBIAN,GRADIENT,PAOPTI,
     '  PG,RG,SE,SOB_VALUE,WG,WU,XA,XE,XG,XIG,XP,
     '  CALCJAC,ERROR,*)

C#### Subroutine: SOBOLEV
C###  Description:
C###    SOBOLEV evaluates the Sobolev value of the mesh and the
C###    jacobian/gradient of the value wrt to the optimisation
C###    parameters for data fitting by optimisation.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NLL(12,NEM),NONL(NLM),NONY(0:NOYM,NYM,NRCM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),nr,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NYNO(0:NYOM,NOOPM,NRCM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 DL(3,NLM),DXIX(3,3),GL(3,3),GRADIENT(*),GU(3,3),
     '  JACOBIAN(NREM,*),PAOPTI(*),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SOB_VALUE,XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),WG(NGM,NBM),
     '  WU(0:NUM+1,NEM)
      LOGICAL CALCJAC
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,nb,nc,ne,ng,nj,nk,nl,nn,noelem,no2,
     '  nol,noopti,noy,np,
     '  nrc,ns,nu,NUTOT,nv,ny
      INTEGER MAP1(4,2),MAP2(3,2),MAP3(3,2)
      REAL*8 PXI,PXIG(10),SFACTOR,SUM1,SUM2,TEMP_VAL

      DATA MAP1 /1,1,2,2,3,4,3,4/ ! MAP1(nn,nk-1) gives local line # for
                                  ! local node nn and deriv nk (bicubic)
      DATA MAP2 /0,3,3,0,1,2/     ! MAP3 for apex node 1 hermite-simplex
      DATA MAP3 /1,1,0,2,3,0/     ! MAP2 for apex node 3 hermite-simplex

      CALL ENTERS('SOBOLEV',*9999)

      nc=1 ! temporary
      nrc=2 ! temporary AJP 30-11-94
      nv=1 ! temporary cpb 22/11/94

      G1SCALING=.FALSE. ! temporary

      IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.1) THEN
        NUTOT=3
      ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.2) THEN
        NUTOT=6
      ELSE
        ERROR='NIT=3 not implemented'
        GOTO 9999
      ENDIF

      SOB_VALUE=0.d0
      IF(CALCJAC) THEN
        IF(KTYP29.EQ.1.OR.KTYP1B.EQ.2) THEN
          DO noopti=1,NTOPTI
            JACOBIAN(NT_RES,noopti)=0.d0
          ENDDO
        ENDIF
      ENDIF
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        WU(NUM+1,ne)=0.d0
        DO nj=1,NJT
          nb=NBJ(nj,ne)
          DO nu=2,NUTOT
            SFACTOR=WU(nu,ne)*WU(0,ne)
            SUM1=0.d0

C cpb 13/5/98 Adding the Jacobian to the Sobolev integration
            DO ng=1,NGT(nb)
C LKC 17-MAY-1998 Dimensions have not been dropped correctly
C              CALL XEXG(NBJ(1,ne),ng,NRE,PG,XE,XG,ERROR,*9999)
C              CALL XGMG(0,NIT(nb),nb,NRE,DXIX,GL,GU,RG(ng),XG,
              CALL XEXG(NBJ(1,ne),ng,NRE(ne),PG,XE,XG,ERROR,*9999)
              CALL XGMG(0,NIT(nb),nb,NRE(ne),DXIX,GL,GU,RG(ng),XG,
     '          ERROR,*9999)
              PXIG(ng)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '          nu,XIG(1,ng,nb),XE(1,nj))
              SUM1=SUM1+PXIG(ng)*PXIG(ng)*WG(ng,nb)*RG(ng)
            ENDDO
            SUM1=SUM1*SFACTOR
            IF(CALCJAC) THEN
              ns=0
              DO nn=1,NNT(nb)
                np=NPNE(nn,nb,ne)
                DO nk=1,NKT(nn,nb)
                  ns=ns+1
                  SUM2=0.d0
                  DO ng=1,NGT(nb)
                    SUM2=SUM2+PXIG(ng)*PG(ns,nu,ng,nb)*WG(ng,nb)*RG(ng)
                  ENDDO
                  SUM2=2.d0*SUM2*SFACTOR
                  ny=NYNP(nk,nv,nj,np,nrc,nc,nr)
C                   noopti=NONY(1,ny,2)
C                  IF(noopti.GT.0) THEN ! dof is in the fit
C                     IF(KTYP29.EQ.1) THEN
C                      JACOBIAN(NT_RES,noopti)=JACOBIAN(NT_RES,noopti)+
C     '                  SUM2*SE(ns,nb,ne)
C                     ELSE IF(KTYP29.EQ.2) THEN
C                      GRADIENT(noopti)=GRADIENT(noopti)+SUM2*
C     '                  SE(ns,nb,ne)
C                     ENDIF
C                  ENDIF
C GMH 21/3/96 add for angle as derivative
                  DO noy=1,NONY(0,ny,2)
                    noopti=NONY(noy,ny,2)
C This code is taken from RESFUN (with minor modifications)
                    IF((nk.EQ.2.OR.nk.EQ.3).AND.KTYP1B.EQ.2) THEN !derivative as angle
                      IF(NYNO(0,noopti,2).EQ.3) THEN !two angles
                        IF(noy.EQ.1) THEN !del/del theta
                          no2=NONY(2,ny,2) !get phi
                          IF(NYNO(1,noopti,2).EQ.ny) THEN !x
                            TEMP_VAL=
     '                        -SUM2*SE(ns,nb,ne)
     '                        *DSIN(PAOPTI(no2))*
     '                        DSIN(PAOPTI(noopti))
                          ELSEIF(NYNO(2,noopti,2).EQ.ny) THEN !y
                            TEMP_VAL=
     '                        SUM2*SE(ns,nb,ne)*
     '                        DSIN(PAOPTI(no2))*
     '                        DCOS(PAOPTI(noopti))
                          ELSEIF(NYNO(3,noopti,2).EQ.ny) THEN !z
                            TEMP_VAL=0.0d0
                          ELSE
                            ERROR='>>Invalid component'
                            GOTO 9999
                          ENDIF
                        ELSEIF(noy.EQ.2) THEN !del/del phi
                          no2=NONY(1,ny,2) !get theta
                          IF(NYNO(1,noopti,2).EQ.ny) THEN !x
                            TEMP_VAL=
     '                        SUM2*SE(ns,nb,ne)*
     '                        DCOS(PAOPTI(noopti))*
     '                        DCOS(PAOPTI(no2))
                          ELSEIF(NYNO(2,noopti,2).EQ.ny) THEN !y
                            TEMP_VAL=
     '                        SUM2*SE(ns,nb,ne)*
     '                        DCOS(PAOPTI(noopti))*
     '                        DSIN(PAOPTI(no2))
                          ELSEIF(NYNO(3,noopti,2).EQ.ny) THEN !z
                            TEMP_VAL=
     '                        -SUM2*SE(ns,nb,ne)*
     '                        DSIN(PAOPTI(noopti))
                          ELSE
                            ERROR='>>Invalid component'
                            GOTO 9999
                          ENDIF
                        ELSE
                          ERROR='>>Only two angles '
     '                      //'for three components'
                          GOTO 9999
                        ENDIF
                      ELSEIF(NYNO(0,noopti,2).EQ.2) THEN !one angle
                        IF(noy.EQ.1) THEN !del/del theta
                          IF(NYNO(1,noopti,2).EQ.ny) THEN !equiv of x
                            TEMP_VAL=
     '                        -SUM2*SE(ns,nb,ne)*
     '                        DSIN(PAOPTI(noopti))
                          ELSEIF(NYNO(2,noopti,2).EQ.ny) THEN !y
                            TEMP_VAL=
     '                        SUM2*SE(ns,nb,ne)*
     '                        DCOS(PAOPTI(noopti))
                          ELSE
                            ERROR='>>Invalid component'
                            GOTO 9999
                          ENDIF
                        ELSE
                          ERROR='>>Only one angle '
     '                      //'for two components'
                          GOTO 9999
                        ENDIF
                      ELSE
                        ERROR='>>Must be two '
     '                    //'or three components'
                        GOTO 9999
                      ENDIF !first angle
                    ELSE
                      TEMP_VAL=SUM2*SE(ns,nb,ne)
                    ENDIF
                    IF(KTYP29.EQ.1) THEN
                      JACOBIAN(NT_RES,noopti)=JACOBIAN(NT_RES,noopti)+
     '                  TEMP_VAL
                    ELSE IF(KTYP29.EQ.2) THEN
                      GRADIENT(noopti)=GRADIENT(noopti)+TEMP_VAL
                    ENDIF
                  ENDDO !noy
                  IF(G1SCALING) THEN
                    IF(nk.GT.1) THEN ! Could have a line in the fit
                                     ! associated with it
                      IF(NIT(nb).EQ.1) THEN ! 1D Element
                        nl=NLL(1,ne)
                        IF(nl.GT.0) THEN
                          noopti=NONL(nl)
                          IF(noopti.GT.0) THEN ! Line is in the optimisation
                            IF(KTYP29.EQ.1) THEN
                              JACOBIAN(NT_RES,noopti)=
     '                          JACOBIAN(NT_RES,noopti)+SUM2*
     '                          XP(nk,nv,nj,np)
                            ELSE IF(KTYP29.EQ.2) THEN
                              GRADIENT(noopti)=GRADIENT(noopti)+SUM2*
     '                          XP(nk,nv,nj,np)
                            ENDIF
                          ENDIF
                        ENDIF !nl.GT.0
                      ELSE IF(NIT(nb).EQ.2) THEN ! 2D Element
                        IF(IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2)
     '                    THEN ! Bicubic-Hermite
                          IF(nk.eq.2.or.nk.EQ.3) THEN
                            nol=MAP1(nn,nk-1)
                            nl=NLL(nol,ne)
                            IF(nl.GT.0) THEN
                              noopti=NONL(nl)
                              IF(noopti.GT.0) THEN ! line is optimised
                                IF(KTYP29.EQ.1) THEN
                                  JACOBIAN(NT_RES,noopti)=
     '                              JACOBIAN(NT_RES,noopti)+SUM2*
     '                              XP(nk,nv,nj,np)
                                ELSE IF(KTYP29.EQ.2) THEN
                                  GRADIENT(noopti)=
     '                              GRADIENT(noopti)+SUM2*XP(nk,nv,nj,
     '                              np)
                                ENDIF
                              ENDIF
                            ENDIF !nl.GT.0
                          ELSE IF(nk.EQ.4) THEN
                            DO i=2,3
                              nol=MAP1(nn,nk-i)
                              nl=NLL(nol,ne)
                              IF(nl.GT.0) THEN
                                noopti=NONL(nl)
                                IF(noopti.GT.0) THEN ! line is optimised
                                  IF(KTYP29.EQ.1) THEN
                                    JACOBIAN(NT_RES,noopti)=
     '                                JACOBIAN(NT_RES,noopti)+SUM2*
     '                                XP(nk,nv,nj,np)*DL(3,nl)
                                  ELSE IF(KTYP29.EQ.2) THEN
                                    GRADIENT(noopti)=GRADIENT(noopti)
     '                                +SUM2*XP(nk,nv,nj,np)*DL(3,nl)
                                  ENDIF
                                ENDIF
                              ENDIF
                            ENDDO
                          ENDIF
C LKC 4-DEC-97 Previously not smoothing for hermite-simplex elems
C       ELSE IF(IBT(1,1,nb).EQ.3.AND.IBT(1,2,nb).EQ.3)

                        ELSE IF(IBT(1,1,nb).EQ.3.AND.IBT(1,2,nb).EQ.4)
     '                    THEN ! Hermite-Simplex
                          IF(nk.eq.2.or.nk.EQ.3) THEN
                            IF(NKT(1,nb).EQ.1) THEN ! Apex at node 1
                              nol=MAP2(nn,nk-1)
                            ELSE                    ! Apex at node 3
                              nol=MAP3(nn,nk-1)
                            ENDIF
                            nl=NLL(nol,ne)
                            IF(nl.GT.0) THEN
                              noopti=NONL(nl)
                              IF(noopti.GT.0) THEN ! line is optimised
                                IF(KTYP29.EQ.1) THEN
                                  JACOBIAN(NT_RES,noopti)=
     '                              JACOBIAN(NT_RES,noopti)+SUM2*
     '                              XP(nk,nv,nj,np)
                                ELSE IF(KTYP29.EQ.2) THEN
                                  GRADIENT(noopti)=
     '                              GRADIENT(noopti)+SUM2*XP(nk,nv,nj,
     '                              np)
                                ENDIF
                              ENDIF
                            ENDIF !nl.GT.0
                          ELSE IF(nk.EQ.4) THEN
                            DO i=2,3
                              IF(NKT(1,nb).EQ.1) THEN ! Apex at node 1
                                nol=MAP2(nn,nk-i)
                              ELSE                    ! Apex at node 3
                                nol=MAP3(nn,nk-i)
                              ENDIF
                              nl=NLL(nol,ne)
                              IF(nl.GT.0) THEN
                                noopti=NONL(nl)
                                IF(noopti.GT.0) THEN ! line is optimised
                                  IF(KTYP29.EQ.1) THEN
                                    JACOBIAN(NT_RES,noopti)=
     '                                JACOBIAN(NT_RES,noopti)+SUM2*
     '                                XP(nk,nv,nj,np)*DL(3,nl)
                                  ELSE IF(KTYP29.EQ.2) THEN
                                    GRADIENT(noopti)=GRADIENT(noopti)
     '                                +SUM2*XP(nk,nv,nj,np)*DL(3,nl)
                                  ENDIF
                                ENDIF
                              ENDIF
                            ENDDO
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
            SOB_VALUE=SOB_VALUE+SUM1
            WU(NUM+1,ne)=WU(NUM+1,ne)+SUM1
          ENDDO
        ENDDO
      ENDDO

C      DO noelem=1,NEELEM(0,nr)
C        ne=NEELEM(noelem,nr)
C        DO nj=1,NJT
C          nb=NBJ(nj,ne)
C          DO nu=1,NUTOT
C            ns1=0
C            DO nn1=1,NNT(nb)
C              NP1=NPNE(nn1,nb,ne)
C              DO nk1=1,NKT(nn1,nb)
C                ns1=ns1+1
C                ns2=0
C                DO nn2=1,NNT(nb)
C                  NP2=NPNE(nn2,nb,ne)
C                  DO nk2=1,NKT(nn2,nb)
C                    ns2=ns2+1
C                    SUM1=PGG(ns1,ns2,nu+1,nb)*WU(nu,ne)*WU(0,ne)
C                    SOB_VALUE=SOB_VALUE+SUM1*
C     '                XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '                XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                    IF(CALCJAC) THEN
C                      ny=NYNP(nk1,nv,nj,np1,nc,nr)
C                      noopti=NONY(1,ny,2)
C                      IF(noopti.GT.0) THEN ! Geometric param. is optimised
C                        IF(KTYP29.EQ.1) THEN
C                           JACOBIAN(NT_RES,noopti)=
C     '                     JACOBIAN(NT_RES,noopti)+SUM1*
C     '                     SE(ns1,nb,ne)*
C     '                     XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                        ELSE IF(KTYP29.EQ.2) THEN
C                           GRADIENT(noopti)=
C     '                       GRADIENT(noopti)+SUM1*
C     '                       SE(ns1,nb,ne)*
C     '                       XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                        ENDIF
C                      ENDIF
C                      ny=NYNP(nk2,nv,nj,np2,nc,nr)
C                      noopti=NONY(1,ny,2)
C                      IF(noopti.GT.0) THEN ! Geometric param. is optimised
C                        IF(KTYP29.EQ.1) THEN
C                         JACOBIAN(NT_RES,noopti)=
C     '                     JACOBIAN(NT_RES,noopti)+SUM1*
C     '                     XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '                     SE(ns2,nb,ne)
C                        ELSE IF(KTYP29.EQ.2) THEN
C                           GRADIENT(noopti)=
C     '                       GRADIENT(noopti)+SUM1*
C     '                       XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '                       SE(ns2,nb,ne)
C                        ENDIF
C                      ENDIF
C                      IF(NIT(nb).EQ.1) THEN ! 1D Element
C                        IF(nk1.GT.1) THEN
C                          noopti=NONL(NLL(1,ne))
C                          IF(noopti.GT.0) THEN ! Line is optimised
C                            IF(KTYP29.EQ.1) THEN
C                  JACOBIAN(NT_RES,noopti)=
C     '        JACOBIAN(NT_RES,noopti)+SUM1*
C     '              XP(nk1,nv,nj,NP1)*
C     '        XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                            ELSE IF(KTYP29.EQ.2) THEN
C                  GRADIENT(noopti)=
C     '        GRADIENT(noopti)+SUM1*
C     '        XP(nk1,nv,nj,NP1)*
C     '        XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                            ENDIF
C                          ENDIF
C                        ENDIF
C                        IF(nk2.GT.1) THEN
C                          noopti=NONL(NLL(1,ne))
C                          IF(noopti.GT.0) THEN ! Line is optimised
C                            IF(KTYP29.EQ.1) THEN
C                  JACOBIAN(NT_RES,noopti)=
C     '        JACOBIAN(NT_RES,noopti)+SUM1*
C     '        XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '        XP(nk2,nv,nj,NP2)
C                            ELSE IF(KTYP29.EQ.2) THEN
C                  GRADIENT(noopti)=
C     '        GRADIENT(noopti)+SUM1*
C     '        XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '        XP(nk2,nv,nj,NP2)
C                            ENDIF
C                          ENDIF
C                        ENDIF
C                      ELSE IF(NIT(nb).EQ.2) THEN ! 2D Element
CC cpb 4/2/94 the following code for a 2d element was implemented quickly
CC and is not efficient for handling bihermite and hermite-simplex elements.
C
CC NOTE: non-standard line mappings are not maintained here!
C
C                        IF(IBT(1,1,nb).EQ.2.AND.
C     '                    IBT(1,2,nb).EQ.2) THEN ! BiHermite
C                          IF(nk1.EQ.2.OR.nk1.EQ.3) THEN
C                            nol=MAP1(nn1,nk1-1)
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              GRADIENT(noopti)=
C     '          GRADIENT(noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ENDIF
C                            ENDIF
C                          ELSE IF(nk1.EQ.4) THEN
C                            nol=MAP1(nn1,nk1-3)
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              GRADIENT(noopti)=
C     '          GRADIENT(noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ENDIF
C                            ENDIF
C                            nol=MAP1(nn1,nk1-2)
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              GRADIENT(noopti)=
C     '          GRADIENT(noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ENDIF
C                            ENDIF
C                          ENDIF
C                          IF(nk2.EQ.2.OR.nk2.EQ.3) THEN
C                            nol=MAP1(nn2,nk2-1)
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              GRADIENT(noopti)=
C     '          GRADIENT(noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ENDIF
C                            ENDIF
C                          ELSE IF(nk2.EQ.4) THEN
C                            nol=MAP1(nn2,nk2-3)
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              GRADIENT(noopti)=
C     '          GRADIENT(noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ENDIF
C                            ENDIF
C                            nol=MAP1(nn2,nk2-2)
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              GRADIENT(noopti)=
C     '          GRADIENT(noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ENDIF
C                            ENDIF
C                          ENDIF
C                        ELSE IF(IBT(1,1,nb).EQ.3.AND.
C     '                    IBT(1,2,nb).EQ.3) THEN ! Hermite-Simplex
C                          IF(nk1.EQ.2.OR.nk1.EQ.3) THEN
C                            IF(NKT(1,nb).EQ.1) THEN ! Apex at node 1
C                              nol=MAP2(nn1,nk1-1)
C                            ELSE                    ! Apex at node 3
C                              nol=MAP3(nn1,nk1-1)
C                            ENDIF
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              GRADIENT(noopti)=
C     '          GRADIENT(noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ENDIF
C                            ENDIF
C                          ELSE IF(nk1.EQ.4) THEN
C                            IF(NKT(1,nb).EQ.1) THEN ! Apex at node 1
C                              nol=MAP2(nn1,nk1-3)
C                            ELSE                    ! Apex at node 3
C                              nol=MAP3(nn1,nk1-3)
C                            ENDIF
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ENDIF
C                            ENDIF
C                            IF(NKT(1,nb).EQ.1) THEN ! Apex at node 1
C                              nol=MAP2(nn1,nk1-2)
C                            ELSE                    ! Apex at node 3
C                              nol=MAP3(nn1,nk1-2)
C                            ENDIF
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              GRADIENT(noopti)=
C     '          GRADIENT(noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*
C     '          XP(nk2,nv,nj,NP2)*SE(ns2,nb,ne)
C                              ENDIF
C                            ENDIF
C                          ENDIF
C                          IF(nk2.EQ.2.OR.nk2.EQ.3) THEN
C                            IF(NKT(1,nb).EQ.1) THEN ! Apex at node 1
C                              nol=MAP2(nn2,nk2-1)
C                            ELSE                    ! Apex at node 3
C                              nol=MAP3(nn2,nk2-1)
C                            ENDIF
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ENDIF
C                            ENDIF
C                          ELSE IF(nk2.EQ.4) THEN
C                            IF(NKT(1,nb).EQ.1) THEN ! Apex at node 1
C                              nol=MAP2(nn2,nk2-3)
C                            ELSE                    ! Apex at node 3
C                              nol=MAP3(nn2,nk2-3)
C                            ENDIF
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              GRADIENT(noopti)=
C     '          GRADIENT(noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ENDIF
C                            ENDIF
C                            IF(NKT(1,nb).EQ.1) THEN ! Apex at node 1
C                              nol=MAP2(nn2,nk2-2)
C                            ELSE                    ! Apex at node 3
C                              nol=MAP3(nn2,nk2-2)
C                            ENDIF
C                            noopti=NONL(NLL(nol,ne))
C                            IF(noopti.GT.0) THEN ! line is optimised
C                              IF(KTYP29.EQ.1) THEN
C              JACOBIAN(NT_RES,noopti)=
C     '          JACOBIAN(NT_RES,noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ELSE IF(KTYP29.EQ.2) THEN
C              GRADIENT(noopti)=
C     '          GRADIENT(noopti)+SUM1*
C     '          XP(nk1,nv,nj,NP1)*SE(ns1,nb,ne)*
C     '          XP(nk2,nv,nj,NP2)
C                              ENDIF
C                            ENDIF
C                          ENDIF
C                        ENDIF
C                      ELSE
C                        ERROR='NIT=3 not implemented'
C                        GOTO 9999
C                      ENDIF
C                    ENDIF
C                  ENDDO
C                ENDDO
C              ENDDO
C            ENDDO
C          ENDDO
C        ENDDO
C      ENDDO

      CALCJAC=.FALSE.
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(''Sobolev Value = '',D12.5)') SOB_VALUE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('SOBOLEV')
      RETURN
 9999 CALL ERRORS('SOBOLEV',ERROR)
      CALL EXITS('SOBOLEV')
      RETURN 1
      END


