      SUBROUTINE GLOBALO(IDO,INP,ISIZE_MFI,ISIZE_PHI,LD,LDR,NBH,NBJ,
     '  NEELEM,NENP,NKB,NKHE,NKH,
     '  NLNO,NMNO,NNB,NONL,NONM,NONY,NPL,NPLIST4,NPNE,NPNODE,NPNY,
     '  nr,NVHE,NVHP,nx_opt,nx_sol,NXI,NYNE,NYNO,NYNP,NYNR,NYNY,
     '  PAOPTY,CONY,CYNO,CYNY,PAOPTI,PMAX,PMIN,XP,FIX,ERROR,*)

C#### Subroutine: GLOBALO
C###  Description:
C###    GLOBALO calculates global mapping from nodal, element and line
C###    etc. arrays for region nr for optimisation problems ie NLNO,
C###    NONL,NMNO,NONM,NONY,NYNO,CONY,CYNO.  It calculates and sets up
C###    the optimisation variables, residuals,and constraints.

      IMPLICIT NONE
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISIZE_MFI(3,NSSM),ISIZE_PHI(2),
     '  LD(NDM),LDR(0:NDM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKB(2,2,2,NNM,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NKH(NHM,NPM,NCM),NLNO(NOPM),
     '  NMNO(1:2,0:NOPM),NNB(4,4,4,NBFM),NONL(NLM),NONM(NMM,NPM),
     '  NONY(0:NOYM,NYM,NRCM),NPL(5,0:3,NLM),NPLIST4(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),
     '  nr,NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),nx_opt,nx_sol,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),NYNY(0:NYYM,NYM),
     '  PAOPTY(NOPM)
      REAL*8 CONY(0:NOYM,NYM,NRCM),CYNO(0:NYOM,NOOPM,NRCM),
     '  CYNY(0:NYYM,NYM),PAOPTI(*),PMAX(*),PMIN(*),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER GETNYR,nc,nd,ndr,ne,nh,nj,nj2,nk,
     '  nl,nm,nmlist,no,no2,nocont,noelem,nonode,no_nynr,no_max,
     '  nores,noy,NOYT,np,nss,nrc,nv,nv2,ny,nyo,nyy(2),ny1,ny2
      REAL*8 COY,RATIO
      LOGICAL DERIVIN1,DERIVIN2,DONE,FOUND,
     '  INCR_NOPM,INCR_NOYM,INCR_NYOM,NODEIN1,NODEIN2

      CALL ENTERS('GLOBALO',*9999)

C GMH 19/3/96 Initialise mapping arrays (from globalh)
      DO nrc=1,2
        DO no=1,NOM
          DO nyo=0,NYOM
            NYNO(nyo,no,nrc)=0
            CYNO(nyo,no,nrc)=0.0d0
          ENDDO !nyo
        ENDDO !no
        NOT(nrc,1,nr,nx_opt)=0
      ENDDO !nrc
      DO nc=1,NCT(nr,nx_opt) !GK, GQ (and GD) variables
        DO no_nynr=1,NYNR(0,0,nc,nr,nx_opt) !Loop over variables in nrr
          ny=NYNR(no_nynr,0,nc,nr,nx_opt) !variable #
          ny2=GETNYR(1,NPNY,nr,1,0,ny,NYNE,NYNP) !equiv rhs row #
          nyy(1)=ny2
          nyy(2)=ny
          DO nrc=1,2
            DO noy=0,NOYM
              NONY(noy,nyy(nrc),nrc)=0
              CONY(noy,nyy(nrc),nrc)=0.0d0
            ENDDO !noy
          ENDDO !nrc
        ENDDO !no_nynr (ny)
      ENDDO !nc

CC AJPs Moved inside data fitting loop.  Others do not use NH. 191297
C      CALL ASSERT(NHM.EQ.NJM,'>>NHM must be equal to NJM',ERROR,*9999)
CC AJPe
      IF(KTYP26.EQ.1) THEN      !material parameter optimisation
C news AJP 17/3/96
        NTOPTI=0
        DO nmlist=1,NMNO(1,0) !list of material params in fit
          nm=NMNO(1,nmlist) !is material param number (setup in IPOPTI)
C news HS 26/4/04
C added proper setup of NTOPTI for grid coupling problems
          IF(ILP(nm,1,nr,nx_sol).EQ.1.OR.      !material param constant
     &      (KTYP54(nr).EQ.3.AND.KTYP3B.EQ.2)) THEN! then grid coupling
C newe
C          IF(ILP(nm,1,nr,nx_sol).EQ.1) THEN !material param constant
            NTOPTI=NTOPTI+1
            PAOPTY(NTOPTI)=3
            NMNO(2,nmlist)=nmlist
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NONM(nm,ne)=NTOPTI
            ENDDO !noelem (ne)
          ELSE IF(ILP(nm,1,nr,nx_sol).EQ.2) THEN !piecewise constant
            NMNO(2,nmlist)=NEELEM(1,nr)
            NONM(nm,1)=NEELEM(1,nr)
            NTOPTI=NTOPTI+1
            PAOPTY(NTOPTI)=3
            no_max=NMNO(1,0)*NEELEM(0,nr)
            CALL ASSERT(no_max.LE.NOPM,'>>Increase NOPM',ERROR,*9999)
            DO noelem=2,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NMNO(1,NMNO(1,0)*(noelem-1)+nmlist)=nm
              NMNO(2,NMNO(1,0)*(noelem-1)+nmlist)=ne
              NTOPTI=NTOPTI+1
              PAOPTY(NTOPTI)=3
              NONM(nm,noelem)=NTOPTI
              PMIN(NMNO(1,0)*(noelem-1)+nmlist)=PMIN(nmlist)
              PMAX(NMNO(1,0)*(noelem-1)+nmlist)=PMAX(nmlist)
              PAOPTI(NMNO(1,0)*(noelem-1)+nmlist)=PAOPTI(nmlist)
            ENDDO !noelem
          ELSE IF(ILP(nm,1,nr,nx_sol).EQ.3) THEN !piecewise linear
            NMNO(2,nmlist)=NPNODE(1,nr)
            NONM(nm,1)=NPNODE(1,nr)
            NTOPTI=NTOPTI+1
            PAOPTY(NTOPTI)=3
            no_max=NMNO(1,0)*NPNODE(0,nr)
            CALL ASSERT(no_max.LE.NOPM,'>>Increase NOPM',ERROR,*9999)
            DO nonode=2,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              NMNO(1,NMNO(1,0)*(nonode-1)+nmlist)=nm
              NMNO(2,NMNO(1,0)*(nonode-1)+nmlist)=np
              NTOPTI=NTOPTI+1
              PAOPTY(NTOPTI)=3
              NONM(nm,nonode)=NTOPTI
              PMIN(NMNO(1,0)*(nonode-1)+nmlist)=PMIN(nmlist)
              PMAX(NMNO(1,0)*(nonode-1)+nmlist)=PMAX(nmlist)
              PAOPTI(NMNO(1,0)*(nonode-1)+nmlist)=PAOPTI(nmlist)
            ENDDO !nonode
          ELSE IF(ILP(nm,1,nr,nx_sol).EQ.4) THEN !Gauss points
C !!!       needs completing
          ENDIF !ilp
        ENDDO !nmlist
C newe AJP 17/3/96
        IF(KTYP27.EQ.5) THEN
          ndr=0
          DO nd=1,NDT
            IF(LD(nd).NE.0) THEN
              ndr=ndr+1
              LDR(ndr)=nd
            ENDIF
          ENDDO
          LDR(0)=ndr
        ENDIF
        IF(KTYP27.EQ.6) THEN
          LDR(0)=NTOPTI
        ENDIF

      ELSE IF(KTYP26.EQ.2) THEN !geometric parameter optimisation
        IF(KTYP27.EQ.1) THEN
        ELSE IF(KTYP27.EQ.2) THEN
        ELSE IF(KTYP27.EQ.3) THEN
        ELSE IF(KTYP27.EQ.4) THEN
        ELSE IF(KTYP27.EQ.5) THEN !Data fitting
CC AJPs Moved from above 191297
          CALL ASSERT(NHM.EQ.NJM,
     '      '>>NHM must be equal to NJM',ERROR,*9999)
CC AJPe

C*** Calculate the optimisation variables

          INCR_NOPM=.FALSE.
          INCR_NOYM=.FALSE.
          INCR_NYOM=.FALSE.
          no=0
          DO no_nynr=1,NYNR(0,0,1,nr,nx_opt)
            ny=NYNR(no_nynr,0,1,nr,nx_opt)
            CONY(0,ny,2)=0.0d0
            IF(.NOT.FIX(ny,1,nx_opt)) THEN !mesh dof is in the fit
              DONE=.FALSE.

              IF(JTYP2A.EQ.1) THEN !treat versions as coincident
C MPN 3Jul2003: dropping down of indices was all wrong here!!!
                CALL GETEQVNONY(IDO,INP,NBH,NBJ,NENP(1,0,nr),NKB,NKHE,
     '            NNB,NONY,NPNE,NPNY,nr,NVHE,NVHP(1,1,1,nr),NXI,
     '            ny,ny2,NYNP(1,1,1,1,0,1,nr),NYNY,CYNY,RATIO,
     '            FIX(1,1,nx_opt),*9999)
C MPN 3Jul2003 OLD
C                CALL GETEQVNONY(IDO,INP,NBH,NBJ,NENP(1,0,nr),NKB,NKHE,
C     '            NNB,NONY,NPNE,NPNY,nr,NVHE,NVHP(1,1,1,nr),NXI,
C     '            ny,ny2,NYNP,NYNY,CYNY,RATIO,
C     '            FIX(1,1,nx_opt),*9999)
                IF(ny2.NE.ny) THEN
C                 There is an equivalent mesh degree of freedom for which an
C                 no as already been assigned.
                  DONE=.TRUE.
                  IF(ny2.EQ.0) THEN !dof not used
C                   Set to zero.
                    FIX(ny,1,nx_opt)=.TRUE.
                    nk=NPNY(1,ny,0)
                    nv=NPNY(2,ny,0)
                    nj=NPNY(3,ny,0)
                    np=NPNY(4,ny,0)
                    XP(nk,nv,nj,np)=0.0d0
                  ELSE IF(FIX(ny2,1,nx_opt)) THEN
C                   Set dof to same value as equivalent dof.
                    FIX(ny,1,nx_opt)=.TRUE.
                    nk=NPNY(1,ny,0)
                    nv=NPNY(2,ny,0)
                    nv2=NPNY(2,ny2,0)
                    nj=NPNY(3,ny,0)
                    np=NPNY(4,ny,0)
                    XP(nk,nv,nj,np)=RATIO*XP(nk,nv2,nj,np)
                  ELSE
                    NOYT=NONY(0,ny2,2)
                    NONY(0,ny,2)=NOYT !CYNO(0) is already 0
                    DO noy=1,NOYT
                      no2=NONY(noy,ny2,2)
                      NONY(noy,ny,2)=no2
                      COY=RATIO*CONY(noy,ny2,2)
                      CONY(noy,ny,2)=COY
                      IF(no2.LT.NOPM) THEN
                        nyo=NYNO(0,no2,2)+1
                        NYNO(0,no2,2)=nyo
                        IF(nyo.LE.NYOM) THEN
                          NYNO(nyo,no2,2)=ny
                          CYNO(nyo,no2,2)=COY
                        ELSE
                          INCR_NYOM=.TRUE.
                        ENDIF
                      ENDIF
                    ENDDO ! noy
                  ENDIF !ny=0/FIX
                ENDIF !ny2.NE.ny
              ENDIF !coincident versions

C GM 19/3/96 Check for derivative
              IF(.NOT.DONE.AND.KTYP1B.EQ.2) THEN !derivative as angle
                nk=NPNY(1,ny,0)
                IF((nk.EQ.2.OR.nk.EQ.3)) THEN
                  DONE=.TRUE.
                  nv=NPNY(2,ny,0)
                  nj=NPNY(3,ny,0)
                  np=NPNY(4,ny,0)
                  FOUND=.FALSE.
                  !loop over from nj to find number of free components
                  DO nj2=nj+1,NJT
                    nh=NH_LOC(nj2,nx_opt)
                    ny2=NYNP(nk,nv,nh,np,2,1,nr)
                    IF(.NOT.FIX(ny2,1,nx_opt)) THEN !another component
                      FOUND=.TRUE.
                    ENDIF
                  ENDDO
                  IF(FOUND) THEN !add to optimisation
                    no=no+1
                    IF(no.LE.NOPM) THEN
                      PAOPTY(no)=1 !Parameter is a geometric dof
                    ELSE
                      INCR_NOPM=.TRUE.
                    ENDIF
                    !tell all ny they use this no
                    DO nj2=1,NJT
                      nh=NH_LOC(nj2,nx_opt)
                      ny2=NYNP(nk,nv,nh,np,2,1,nr)
                      IF(.NOT.FIX(ny2,1,nx_opt)) THEN !another component
                        NONY(0,ny2,2)=NONY(0,ny2,2)+1
                        IF(NONY(0,ny2,2).LE.NOYM) THEN
                          NONY(NONY(0,ny2,2),ny2,2)=no
                          CONY(NONY(0,ny2,2),ny2,2)=0.0d0
                        ELSE
                          INCR_NOYM=.TRUE.
                        ENDIF
                        !tell this no it uses all ny
                        IF(no.LE.NOOPM) THEN
                          NYNO(0,no,2)=NYNO(0,no,2)+1
                          IF(NYNO(0,no,2).LE.NYOM) THEN
                            NYNO(NYNO(0,no,2),no,2)=ny2
                            CYNO(0,no,2)=0.0d0
                            CYNO(NYNO(0,no,2),no,2)=0.0d0
                          ELSE
                            INCR_NYOM=.TRUE.
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDIF
                ENDIF !nk
              ENDIF !derivative as angle

              IF(.NOT.DONE) THEN
                no=no+1
                IF(no.LE.NOPM) THEN
                  PAOPTY(no)=1 !Parameter is a geometric dof
                ELSE
                  INCR_NOPM=.TRUE.
                ENDIF
!                PAOPTY(no)=1 !Parameter is a geometric dof
                NONY(0,ny,2)=1
                NONY(1,ny,2)=no
                CONY(1,ny,2)=1.0d0
                IF(no.LE.NOOPM) THEN
                  NYNO(0,no,2)=1
                  NYNO(1,no,2)=ny
                  CYNO(0,no,2)=0.0d0
                  CYNO(1,no,2)=1.0d0
                ENDIF
              ENDIF !.NOT.DONE
            ENDIF !.NOT.FIX
          ENDDO !no_nynr

          CALL ASSERT(.NOT.INCR_NOPM,'>>Increase NOPM',ERROR,*9999)
          CALL ASSERT(.NOT.INCR_NOYM,'>>Increase NOYM',ERROR,*9999)
          CALL ASSERT(.NOT.INCR_NYOM,'>>Increase NYOM',ERROR,*9999)
          NOT(2,1,nr,nx_opt)=no
          IF(G1SCALING) THEN
            DO nl=1,NLT
              IF(JTYP2B.EQ.1.AND.NPL(4,0,nl).LT.0) THEN
C               Line is mapped to another line so don't include it
                NONL(nl)=0
              ELSE
C               See if the nodes at the end of the line are in the fit
                NODEIN1=.FALSE.
                DERIVIN1=.FALSE.
                NODEIN2=.FALSE.
                DERIVIN2=.FALSE.
                np=NPL(2,1,nl)
                nk=NPL(4,1,nl)
C               Check if the node is in the optimisation problem
                DO nj=1,NJT
                  DO nv=1,NVHP(nj,np,1,nr)
                    ny1=NYNP(1,nv,nj,np,2,1,nr)
                    IF(.NOT.FIX(ny1,1,nx_opt)) THEN
                      NODEIN1=.TRUE.
                    ENDIF
                    IF(nk.NE.0) THEN !Line has a derivative term
                      ny2=NYNP(nk,nv,nj,np,2,1,nr)
                      IF(.NOT.FIX(ny2,1,nx_opt)) THEN
                        DERIVIN1=.TRUE.
                      ENDIF
                    ENDIF
                  ENDDO !nv
                ENDDO !nj
                np=NPL(3,1,nl)
                nk=NPL(5,1,nl)
C               Check if the node is in the optimisation problem
                DO nj=1,NJT
                  DO nv=1,NVHP(nj,np,1,nr)
                    ny1=NYNP(1,nv,nj,np,2,1,nr)
                    IF(.NOT.FIX(ny1,1,nx_opt)) THEN
                      NODEIN2=.TRUE.
                    ENDIF
                    IF(nk.NE.0) THEN !Line has a derivative term
                      ny2=NYNP(nk,nv,nj,np,2,1,nr)
                      IF(.NOT.FIX(ny2,1,nx_opt)) THEN
                        DERIVIN2=.TRUE.
                      ENDIF
                    ENDIF
                  ENDDO !nv
                ENDDO !nj
                IF(NODEIN1.OR.NODEIN2.OR.DERIVIN1.OR.DERIVIN2) THEN !Line nl is to be optimised
                  no=no+1
                  IF(no.LE.NOPM) THEN
                    PAOPTY(no)=2 !Parameter is a line
                  ENDIF
                  NYNO(0,no,2)=0 !no is not coupled to a ny but a nl
                  CYNO(0,no,2)=0.0d0
                  IF(no.LE.NOPM) THEN
                    NONL(nl)=no
                    NLNO(no)=nl
                  ENDIF
                ELSE
                  NONL(nl)=0
                ENDIF
              ENDIF
            ENDDO !End of nl loop
          ENDIF
          NTOPTI=no
          CALL ASSERT(NTOPTI.LE.NOPM,'>>Increase NOPM',ERROR,*9999)

C*** Calculate the constraints
C*** SEN 20/02/03 Commented out since the MinLSSQP optimiser doesn't have
C         nonlinear constraints at the moment. When they are added to the
C         code then we can reable the following.

          nocont=0
C*** SEN 20/02/03 Start of commenting out.
C         IF(KTYP1B.EQ.1) THEN !Derivative constraints
C           DO nonode=1,NPNODE(0,nr)
C             np=NPNODE(nonode,nr)
CC            Find number of versions and derivatives at the node
C             nj=NJ_LOC(NJL_GEOM,1,nr)
C             NKTOT=NKH(nj,np,1)
C             NVTOT=NVHP(nj,np,1,nr)
C             DO njj2=2,NJ_LOC(NJL_GEOM,0,nr)
C               nj=NJ_LOC(NJL_GEOM,njj2,nr)
C               IF(NKH(nj,np,1).NE.NKTOT.AND.NKTOT.NE.0) THEN
C                 WRITE(OP_STRING,
C    '              '('' Differing numbers of derivatives at node '''
C    '              //',I3)') np
C                 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                 NKTOT=0
C               ENDIF
C               IF(NVHP(nj,np,1,nr).NE.NVTOT.AND.NVTOT.NE.0) THEN
C                 WRITE(OP_STRING,
C    '              '('' Differing numbers of versions at node '''
C    '              //',I3)') np
C                 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                 NVTOT=0
C               ENDIF
C             ENDDO
C             DO nk=2,NKTOT
C               IF(nk.EQ.2.OR.nk.EQ.3.OR.nk.EQ.5) THEN !first deriv
C                 INOPTI=.FALSE.
C                 DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
C                   nj=NJ_LOC(NJL_GEOM,njj2,nr)
C                   DO nv=1,NVHP(nj,np,1,nr)
C                     ny=NYNP(nk,nv,nj,np,0,1,nr)
C                     IF(ny.GT.0) THEN
C                       no=NONY(1,ny,2)
C                       IF(no.GT.0) THEN
CC                        Only need one constraint if mesh dof map to
CC                        one solution dof.
C                         IF(NYNO(1,no,2).EQ.ny) INOPTI=.TRUE.
C                       ENDIF !no.GT.0
C                     ENDIF !ny.GT.0
C                   ENDDO !nv
C                 ENDDO !nj
C                 IF(INOPTI) THEN
C                   nocont=nocont+1
C                   IF(nocont.LE.NCOM) THEN
C                     PMIN(nocont+NTOPTI)=1.0d0
C                     PMAX(nocont+NTOPTI)=1.0d0
C                   ENDIF
C                 ENDIF !INOPTI
C               ENDIF !nk
C             ENDDO !nk
C           ENDDO !nonode (np)
C         ENDIF !Derivative constraints
C         IF(G1SCALING) THEN
C           DO nl=1,NLT
C             IF(NONL(nl).GT.0) THEN !The line is an optimising variable
C               nocont=nocont+1
C               IF(nocont.LE.NCOM) THEN
C                 PMIN(nocont+NTOPTI)=0.0d0
C                 PMAX(nocont+NTOPTI)=0.0d0
C               ENDIF
C             ENDIF
C           ENDDO !nl
C         ENDIF !g1scaling
C*** SEN 20/02/03 End of commenting out.
          NTCNTR=nocont
          CALL ASSERT(NTCNTR.LE.NCOM,'>> Increase NCOM',ERROR,*9999)

C*** Calculate the residuals

          ndr=0
          DO nd=1,NDT
            IF(LD(nd).NE.0) THEN
              ndr=ndr+1
              LDR(ndr)=nd
            ENDIF
          ENDDO
          LDR(0)=ndr
          IF(KTYP29.EQ.2.OR.KTYP29B.EQ.3) THEN
            nores=ndr !residual for each projection
          ELSE IF(KTYP29B.EQ.2) THEN
            nores=1 !one residual
          ELSE !KTYP29B.EQ.1
            nores=NJ_LOC(NJL_GEOM,0,nr)*ndr !resid. for each component
          ENDIF
          IF(KTYP12.EQ.1) THEN !Sobolev smoothing
            nores=nores+1
          ENDIF
          NT_RES=nores
          CALL ASSERT(NT_RES.LE.NREM,'>> Increase NREM',ERROR,*9999)




        ELSE IF(KTYP27.EQ.12) THEN !Activation times

C*** Calculate the optimisation variables
          CALL ASSERT(NYOM.GE.1,'>>Increase NYOM',ERROR,*9999)
          CALL ASSERT(NRCM.GE.2,'>>Increase NRCM',ERROR,*9999)
          INCR_NOPM=.FALSE.
          no=0
          DO no_nynr=1,NYNR(0,0,1,nr,nx_sol)
            ny=NYNR(no_nynr,0,1,nr,nx_sol) !Heart nys
            CONY(0,ny,2)=0.0d0
            IF(.NOT.FIX(ny,1,nx_sol)) THEN !mesh dof is in the fit
              no=no+1
              IF(no.LE.NOPM) THEN
                PAOPTY(no)=1 !Parameter is a geometric dof
              ELSE
                INCR_NOPM=.TRUE.
              ENDIF
              NONY(0,ny,2)=1
              NONY(1,ny,2)=no
              CONY(1,ny,2)=1.0d0
              IF(no.LE.NOOPM) THEN
                NYNO(0,no,2)=1
                NYNO(1,no,2)=ny
                CYNO(0,no,2)=0.0d0
                CYNO(1,no,2)=1.0d0
              ENDIF
            ENDIF !fix
          ENDDO !no_nynr

          NTOPTI=no

          IF(ACTN_OPTI_WAVE_PARAM) THEN
            IF(DABS(ACTN_MIN(3)).GT.ZERO_TOL.OR.
     '        DABS(ACTN_MAX(3)).GT.ZERO_TOL) THEN !Optimise Tmem jump
              NTOPTI=NTOPTI+1
              PAOPTY(NTOPTI)=2
            ENDIF
            IF(DABS(ACTN_MIN(4)).GT.ZERO_TOL.OR.
     '        DABS(ACTN_MAX(4)).GT.ZERO_TOL) THEN !Optimise wave width
              NTOPTI=NTOPTI+1
              PAOPTY(NTOPTI)=3
            ENDIF
          ENDIF

          NOT(1,1,nr,nx_opt)=NTOPTI
          NOT(1,1,nr,nx_sol)=no !is this correct ??
          NOT(2,1,nr,nx_opt)=NTOPTI
          NOT(2,1,nr,nx_sol)=no !is this correct ??

          CALL ASSERT(.NOT.INCR_NOPM,'>>Increase NOPM',ERROR,*9999)
          CALL ASSERT(no.LE.NOOPM,'>>Increase NOOPM',ERROR,*9999)


C*** Calculate the constraints
          NTCNTR=0  !No nonlinear constraints (but each activation
                    !time is bounded


C LKC 11-JUL-2002 The number of residuals is also == to ISIZE_PHI(1)

C*** Calculate the residuals
          IF(ACTN_IREGULARISE.LE.1) THEN
            NT_RES=NPLIST4(0)
          ELSE IF(ACTN_IREGULARISE.EQ.2) THEN
            NT_RES=NPLIST4(0)+1
C          ELSE IF(ACTN_IREGULARISE.EQ.3) THEN
C            NT_RES=NPLIST4(0)+1
          ELSE
            ERROR='>>Update code: Unknown residuals'
            GOTO 9999
          ENDIF
          !Number of outer surface nodes which have recordings.
          CALL ASSERT(NT_RES.LE.NREM,'>>Increase NREM',ERROR,*9999)
          CALL ASSERT(NT_RES.LE.NY_TRANSFER_M,
     '      '>>Increase NY_TRANSFER_M',ERROR,*9999)

        ELSE IF(KTYP27.EQ.13) THEN ! Dipole optimisation

C!!! this is specific for a single dipole
          NTOPTI=6  !variables
          NTCNTR=0  !No nonlinear constraints
          nss=1
C LKC 11-JUL-2002 Generalise for optimising against potential and
C     magnetic fields specified by KTYP27B
C          NT_RES=ISIZE_MFI(1,nss) !residuals


C Note: nr=1 (default) is currently also the electrode nr
          IF(KTYP27B.EQ.1) THEN !magnetic electrodes
            NT_RES=ISIZE_MFI(1,nss)
          ELSEIF(KTYP27B.EQ.2) THEN !potential electrodes
            NT_RES=ISIZE_PHI(1)
          ELSEIF(KTYP27B.EQ.3) THEN !both mag and potl.
            NT_RES=ISIZE_MFI(1,nss)+ISIZE_PHI(1)
          ELSE
            ERROR='Unknown KTYP27B - Code needs updating'
            GOTO 9999
          ENDIF

          CALL ASSERT(NT_RES.GE.1,'>> No residuals in problem',
     '      ERROR,*9999)

          CALL ASSERT(NTOPTI.LE.NOPM,'>>Increase NOPM',ERROR,*9999)
          CALL ASSERT(NT_RES.LE.NREM,'>>Increase NREM',ERROR,*9999)


        ENDIF !ktyp27

      ENDIF

      CALL EXITS('GLOBALO')
      RETURN
 9999 CALL ERRORS('GLOBALO',ERROR)
      CALL EXITS('GLOBALO')
      RETURN 1
      END
CC AJPe


