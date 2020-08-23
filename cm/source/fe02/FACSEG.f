      SUBROUTINE FACSEG(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,NKJE,
     '  NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,NVJF,SE,SF,XP,
     '  ERROR,*)

C#### Subroutine: FACSEG
C###  Description:
C###    FACSEG defines the various parameters associated with global
C###    face segments nf=1,NFT and element sides nf=1,NFE(nb).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),
     '  NFFACE(0:NF_R_M,NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NKEF(0:4,16,6,NBFM),NLF(4,NFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM)
      REAL*8 SE(NSM,NBFM,NEM),SF(NSM,NBFM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IDO1,IDO2,IDO3,II(5),inp1,inp2,J,M,MB,mb2,mf,MM,
     '  MP,MP1(4),MP2(4),N,n1,n2,naf,nb,nbe,ne,nee,nef,nf,ni,ni1,ni2,
     '  ni3,nicollapse,NITB,nj,njj1,njj2,nk,NKJF(NKM,NNM,NJM),NKTOT,nl,
     '  NM(3),NM2,NM3,nn,nne,NNTOT,noelem,np,nr,NUMCOLLAPSED,NUMNODES,
     '  NUMNODES2,NUMNODES3,SIMP3D_FACEB
      CHARACTER ERROR_STRING*64
      LOGICAL ALL,COLLAPSED,COLLAPSED1,COLLAPSED2,FACENODE,FOUND

      DATA II/1,2,3,1,2/

      CALL ENTERS('FACSEG',*9999)

      DO nb=1,NBFT

C*** Calculate NNF

C
C     OR 11-09-2006
C
C     Initializing variable NNF(0,6,NBFM) for auxiliary basis functions
C
        IF (NBC(nb).EQ.0) THEN  ! auxiliary basis fubction
          DO nf=1,6             ! Fine as long as the 2nd component of NNF
                                ! is fixed to a constant value
            NNF(0,nf,nb) = 0
          ENDDO                 ! nf
! RGB 11/12/97 adding 3D simplex elements

        ELSEIF(NIT(nb).EQ.3.AND.IBT(1,1,nb).EQ.3) THEN !3D and simplex
          NFE(nb) = 4
          DO nf = 1,NFE(nb)
            NNF(0,nf,nb) = 3
          ENDDO !nf
          NNF(1,1,nb)=1 ! No Xi coordinates for unstructured meshes
          NNF(2,1,nb)=1
          NNF(3,1,nb)=4
          NNF(4,1,nb)=3
          NNF(1,2,nb)=1 ! No Xi coordinates for unstructured meshes
          NNF(2,2,nb)=2
          NNF(3,2,nb)=1
          NNF(4,2,nb)=4
          NNF(1,3,nb)=1 ! No Xi coordinates for unstructured meshes
          NNF(2,3,nb)=3
          NNF(3,3,nb)=2
          NNF(4,3,nb)=1
          NNF(1,4,nb)=1 ! No Xi coordinates for unstructured meshes
          NNF(2,4,nb)=4
          NNF(3,4,nb)=3
          NNF(4,4,nb)=2
          DO nf=1,NNT(nb)
            DO nn=1,NNF(0,1,nb)
              NKEF(0,nn,nf,nb)=1
              DO nk=1,NKT(1,nb)
                NKEF(nk,nn,nf,nb)=nk
              ENDDO
            ENDDO
          ENDDO
        ELSEIF(NIT(nb).EQ.2.AND.IBT(1,1,nb).EQ.3) THEN !2D simplex
          SIMP3D_FACEB=nb ! Simplex face basis for 3D
          NFE(nb)=1
          DO nf=1,NFE(nb)
            NNF(0,nf,nb)=3
          ENDDO !nf
          NNF(1,1,nb)=1       ! No Xi coordinates for unstructured meshes
          NNF(2,1,nb)=1
          NNF(3,1,nb)=2
          NNF(4,1,nb)=3
          DO nf=1,NNT(nb)
            DO nn=1,NNF(0,1,nb)
              NKEF(0,nn,nf,nb)=1
              DO nk=1,NKT(1,nb)
                NKEF(nk,nn,nf,nb)=nk
              ENDDO
            ENDDO
          ENDDO
        ELSE
     '  IF(NNT(nb).GT.0.AND.NKT(0,nb).GT.0.AND.IBT(1,NIT(nb),nb).NE.9)
     '    THEN
C         if there are nodes and derivatives in the basis calc. NNF
          NITB=NIT(nb)
          IF(NITB.EQ.1) THEN
            NFE(nb)=0
          ELSE IF(NITB.EQ.2) THEN
            NNF(0,1,nb)=NNT(nb)
            NNF(1,1,nb)=3
            DO nn=1,NNT(nb)
              NNF(1+nn,1,nb)=nn
              NKEF(0,nn,1,nb)=NKT(nn,nb)
              DO nk=1,NKT(nn,nb)
                NKEF(nk,nn,1,nb)=nk
              ENDDO !nk
            ENDDO !nn
            NFE(nb)=1
          ELSE IF(NITB.EQ.3) THEN

C CPB 27/7/95 Generalising code for sector basis functions.

C           Find the maximum extents of this basis
            DO ni=1,3
              NM(ni)=0
              DO nn=1,NNT(nb)
                IF(INP(nn,ni,nb).GT.NM(ni)) NM(ni)=INP(nn,ni,nb)
              ENDDO !nn
            ENDDO !ni
            nf=1
C           Loop over the xi directions (ni1) normal to the face
            DO ni1=1,NITB
C             Find the xi directions that lie in the face
              ni2=II(ni1+1)
              ni3=II(ni1+2)
              COLLAPSED=IBT(1,ni1,nb).EQ.5.OR.IBT(1,ni1,nb).EQ.6
              DO n1=1,NM(ni1),NM(ni1)-1
                NNTOT=0
                DO nn=1,NNT(nb)
                  FACENODE=.FALSE.
                  IF(COLLAPSED) THEN
                    IF(IBT(1,ni1,nb).EQ.5) THEN
                      IF(INP(nn,ni1,nb).EQ.n1.OR.
     '                  INP(nn,IBT(3,ni1,nb),nb).EQ.1)
     '                  FACENODE=.TRUE.
                    ELSE
                      IF(INP(nn,ni1,nb).EQ.n1.OR.
     '                  INP(nn,IBT(3,ni1,nb),nb).EQ.
     '                  NM(IBT(3,ni1,nb)))
     '                  FACENODE=.TRUE.
                    ENDIF
                  ELSE
                    IF(INP(nn,ni1,nb).EQ.n1) FACENODE=.TRUE.
                  ENDIF
                  IF(FACENODE) THEN
C                   Found a face node
                    NNTOT=NNTOT+1
                    NNF(1+NNTOT,nf,nb)=nn
                  ENDIF
                ENDDO !nn
C               Check that the face is not just a line in the collapsed
C               case
                NUMCOLLAPSED=0
                nicollapse=1
                IF(IBT(1,ni2,nb).EQ.5.OR.IBT(1,ni2,nb).EQ.6) THEN
                  NUMCOLLAPSED=1
                  nicollapse=ni2
                  IF(IBT(2,ni2,nb).EQ.4) THEN
                    NUMNODES2=2
                  ELSE
                    NUMNODES2=IBT(2,ni2,nb)+1
                  ENDIF
                ELSE
                  IF(IBT(1,ni2,nb).EQ.2) THEN
                    NUMNODES2=2
                  ELSE
                    NUMNODES2=IBT(2,ni2,nb)+1
                  ENDIF
                ENDIF
                IF(IBT(1,ni3,nb).EQ.5.OR.IBT(1,ni3,nb).EQ.6) THEN
                  NUMCOLLAPSED=NUMCOLLAPSED+1
                  nicollapse=ni3
                  IF(IBT(2,ni3,nb).EQ.4) THEN
                    NUMNODES3=2
                  ELSE
                    NUMNODES3=IBT(2,ni3,nb)+1
                  ENDIF
                ELSE
                  IF(IBT(1,ni3,nb).EQ.2) THEN
                    NUMNODES3=2
                  ELSE
                    NUMNODES3=IBT(2,ni3,nb)+1
                  ENDIF
                ENDIF
                IF(NUMCOLLAPSED.NE.0.AND.IBT(3,nicollapse,nb).NE.ni1)
     '            THEN
                  IF(nicollapse.EQ.ni2) THEN
                    NUMNODES=(NUMNODES2-1)*NUMNODES3+1
                  ELSE
                    NUMNODES=(NUMNODES3-1)*NUMNODES2+1
                  ENDIF
                ELSE
                  NUMNODES=NUMNODES2*NUMNODES3
                ENDIF
                IF(NNTOT.EQ.NUMNODES) THEN
                  NNF(1,nf,nb)=ni1
C                 Order the local face nodes correctly
                  DO nn=1,NNTOT
                    NM2=16
                    NM3=16
                    M=nn
                    DO N=nn,NNTOT
                      IF(INP(NNF(1+N,nf,nb),ni3,nb).LT.NM3) THEN
                        NM2=INP(NNF(1+N,nf,nb),ni2,nb)
                        NM3=INP(NNF(1+N,nf,nb),ni3,nb)
                        M=N
                      ELSE IF(INP(NNF(1+N,nf,nb),ni3,nb).EQ.NM3) THEN
                        IF(INP(NNF(1+N,nf,nb),ni2,nb).LT.NM2) THEN
                          NM2=INP(NNF(1+N,nf,nb),ni2,nb)
                          M=N
                        ENDIF
                      ENDIF
                    ENDDO
C                   exchange the nn and m components
                    n2=NNF(1+nn,nf,nb)
                    NNF(1+nn,nf,nb)=NNF(1+M,nf,nb)
                    NNF(1+M,nf,nb)=n2
                  ENDDO !nn
                  NNF(0,nf,nb)=NNTOT
C                 Determine derivatives in the face (if any)
                  DO n=1,NNF(0,nf,nb)
                    nn=NNF(1+n,nf,nb)
                    NKTOT=1
                    NKEF(1,n,nf,nb)=1
C cpb 16/6/98 Fixing NKEF for tricubic Hermite elements
C                   First find ni2 derivatives
                    DO nk=1,NKT(nn,nb)
                      IDO1=IDO(nk,nn,ni1,nb)
                      IDO2=IDO(nk,nn,ni2,nb)
                      IDO3=IDO(nk,nn,ni3,nb)
                      IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.1) THEN
                        NKTOT=NKTOT+1
                        NKEF(NKTOT,n,nf,nb)=nk
                      ENDIF
                    ENDDO !nk
C                   Next find ni3 derivatives
                    DO nk=1,NKT(nn,nb)
                      IDO1=IDO(nk,nn,ni1,nb)
                      IDO2=IDO(nk,nn,ni2,nb)
                      IDO3=IDO(nk,nn,ni3,nb)
                      IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.2) THEN
                        NKTOT=NKTOT+1
                        NKEF(NKTOT,n,nf,nb)=nk
                      ENDIF
                    ENDDO !nk
C                   Finally find ni2 and ni3 derivatives
                    DO nk=1,NKT(nn,nb)
                      IDO1=IDO(nk,nn,ni1,nb)
                      IDO2=IDO(nk,nn,ni2,nb)
                      IDO3=IDO(nk,nn,ni3,nb)
                      IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.2) THEN
                        NKTOT=NKTOT+1
                        NKEF(NKTOT,n,nf,nb)=nk
                      ENDIF
                    ENDDO !nk
                    NKEF(0,n,nf,nb)=NKTOT
                  ENDDO !n
                  nf=nf+1
                ENDIF
              ENDDO !n1
            ENDDO !ni1
            NFE(nb)=nf-1
          ENDIF
        ENDIF
      ENDDO !nb

C*** Calculate NFF and NPF

      NFT=0
      DO nr=1,NRT
        NFFACE(0,nr)=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NJ_LOC(NJL_GEOM,0,nr).GE.2.AND.NIT(NBJ(1,ne)).EQ.2) THEN
C           The element is a 2D element in 2/3D space so the face is
C           just the element itself.
            NFT=NFT+1
            NFFACE(0,nr)=NFFACE(0,nr)+1
            IF(NFT.LE.NFM.AND.NFFACE(0,nr).LE.NF_R_M) THEN
              NFFACE(NFFACE(0,nr),nr)=NFT
              nb=NBJ(1,ne)
              DO ni=1,2
                NPF(2*(ni-1)+1,NFT)=ni
                IF(IBT(1,ni,nb).EQ.1) THEN
                  NPF(2*(ni-1)+2,NFT)=IBT(2,ni,nb)
                ELSE IF(IBT(1,ni,nb).EQ.2) THEN
                  IF(IBT(2,ni,nb).EQ.1) THEN
                    NPF(2*(ni-1)+2,NFT)=4
                  ELSE
                    NPF(2*(ni-1)+2,NFT)=4+IBT(2,ni,nb)
                  ENDIF
                ELSE IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
                  NPF(2*(ni-1)+2,NFT)=IBT(2,ni,nb)
                ENDIF
              ENDDO !ni
              NFF(1,ne)=NFT
              NPF(5,NFT)=1
              NPF(6,NFT)=ne
              NPF(7,NFT)=0
              NPF(8,NFT)=1
              NPF(9,NFT)=0
              DO njj1=1,3 !geometry/fibres/field
                DO njj2=1,NJ_LOC(njj1,0,nr)
                  nj=NJ_LOC(njj1,njj2,nr)
                  NBJF(nj,NFT)=NBJ(nj,ne)
                ENDDO !njj2
              ENDDO !njj1
            ENDIF !NFT <= NFM
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' Number of faces NFT='',I6)') NFT
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3.AND.NIT(NBJ(1,ne)).EQ.3)
     '        THEN
C           3D element in 3D space.
            MB=NBJ(1,ne) !is basis type # of parent elem for nj=1
            DO mf=1,NFE(MB)
C             3D simplex elements - RGB
              IF(IBT(1,1,MB).EQ.3) THEN !Simplex
                nf=0
                FOUND=.FALSE.
                DO WHILE(.NOT.FOUND.AND.nf.LT.MIN(NFT,NFM))
                  nf=nf+1
                  IF(NPF(5,nf).EQ.1) THEN ! Check for 1 element only
                    ALL=.TRUE.
                    nee=NPF(6,nf)
                    nef=NPF(8,nf)
                    nbe=NBJ(1,nee)
                    nn=0
                    DO WHILE(nn.LT.NNF(0,nef,nbe).AND.ALL)
                      nn=nn+1
                      nne=NNF(1+nn,nef,nbe)
                      np=NPNE(nne,nbe,nee)
                      MM=NNF(1+nn,mf,MB)
                      MP=NPNE(MM,MB,ne)
                      IF(mp.NE.np) ALL=.FALSE.
                    ENDDO !nn
                    IF(ALL) THEN
                      NFF(mf,ne)=nf
                      NPF(5,nf)=2
                      NPF(7,nf)=ne
                      NPF(9,nf)=mf
                      IF(DOP) THEN
                        WRITE(OP_STRING,'('' ne='',I5,'' mf='',I2,'
     '                    //''' nft='',I5,'' NPF:'',9I6)')
     '                    ne,mf,nf,(NPF(j,nf),j=1,9)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                      FOUND=.TRUE.
                      IF(NRE(nee).NE.nr) THEN
                        NFFACE(0,nr)=NFFACE(0,nr)+1
                        IF(NFFACE(0,nr).LE.NF_R_M) THEN
                          NFFACE(NFFACE(0,nr),nr)=nf
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !nf
                IF(.NOT.FOUND) THEN
                  NFT=NFT+1
                  NFFACE(0,nr)=NFFACE(0,nr)+1
                  IF(NFT.LE.NFM.AND.NFFACE(0,nr).LE.NF_R_M) THEN
                    NFFACE(NFFACE(0,nr),nr)=NFT
                    NFF(mf,ne)=NFT
                    NPF(1,NFT)=1 !No global Xi coordinates for
                    NPF(3,NFT)=2 !unstructured meshes
                    !NPF(2,NFT)=1
                    !NPF(4,NFT)=1
                    IF(IBT(1,1,MB).EQ.1) THEN
                      NPF(2,NFT)=IBT(2,1,MB)
                    ELSE IF(IBT(1,1,MB).EQ.2) THEN
                      IF(IBT(2,1,MB).EQ.1) THEN
                        NPF(2,NFT)=4
                      ELSE
                        NPF(2,NFT)=4+IBT(2,1,MB)
                      ENDIF
                    ELSE IF(IBT(1,1,MB).EQ.5.OR.IBT(1,1,MB).EQ.6)
     '                  THEN
                      NPF(2,NFT)=IBT(2,1,MB)
                    ENDIF
                    IF(IBT(1,1,MB).EQ.1) THEN
                      NPF(4,NFT)=IBT(2,1,MB)
                    ELSE IF(IBT(1,1,MB).EQ.2) THEN
                      IF(IBT(2,1,MB).EQ.1) THEN
                        NPF(4,NFT)=4
                      ELSE
                        NPF(4,NFT)=4+IBT(2,1,MB)
                      ENDIF
                    ELSE IF(IBT(1,1,MB).EQ.5.OR.IBT(1,1,MB).EQ.6)
     '                  THEN
                      NPF(4,NFT)=IBT(2,1,MB)
                    ENDIF
                    NPF(5,NFT)=1
                    NPF(6,NFT)=ne
                    NPF(7,NFT)=0
                    NPF(8,NFT)=mf
                    NPF(9,NFT)=0
                    DO njj1=1,3 !geometry/fibres/field
                      DO njj2=1,NJ_LOC(njj1,0,nr)
                        nj=NJ_LOC(njj1,njj2,nr)
                        NBJF(nj,NFT)=SIMP3D_FACEB
                      ENDDO !njj2
                    ENDDO !njj1
                  ENDIF
                ENDIF
              ELSE !Not 3D simplex
                ni1=NNF(1,mf,MB)
                ni2=II(ni1+1)
                ni3=II(ni1+2)
C             Check to see that the face actually has is a surface by
C             checking for repeated nodes.
C             If the coordinate system is not rect. Cart. then some
C             funny transformations may be done on theta, so repeated
C             nodes do not imply constant theta.  Only assume collapsed
C             if the repeated node is at the origin.
                DO inp1=1,4
                  MP1(inp1)=0
                  MP2(inp1)=0
                ENDDO
                COLLAPSED1=.TRUE.
                COLLAPSED2=.TRUE.
                nn=0
                DO WHILE((COLLAPSED1.OR.COLLAPSED2).
     '            AND.nn.LT.NNF(0,mf,MB))
                  nn=nn+1
                  MM=NNF(1+nn,mf,MB)
                  MP=NPNE(MM,MB,ne)
                  IF(COLLAPSED1) THEN
                    inp1=INP(MM,ni2,MB) !node index in 1st xi dirn.
                    IF(MP1(inp1).EQ.0) THEN !first node found
                      MP1(inp1)=MP
                    ELSE IF(MP1(inp1).NE.MP) THEN !more than 1 node
                      COLLAPSED1=.FALSE.
                    ELSE IF(ITYP10(nr).GE.2) THEN !not rect. Cart.
C                     Only collapsed if r=0 (cylindrical or spherical),
C                     or lambda=0 or mu=0 (prolate or oblate).
                      IF(XP(1,1,NJ_LOC(NJL_GEOM,1,nr),MP).NE.0.0d0.AND
     '                  .(ITYP10(nr).LE.3.OR.
     '                  XP(1,1,NJ_LOC(NJL_GEOM,2,nr),MP).NE.0.0d0)) THEN
                        COLLAPSED1=.FALSE.
                      ENDIF
                    ENDIF
                  ENDIF
                  IF(COLLAPSED2) THEN
                    inp2=INP(MM,ni3,MB) !node index in 2nd xi dirn.
                    IF(MP2(inp2).EQ.0) THEN
                      MP2(inp2)=MP
                    ELSE IF(MP2(inp2).NE.MP) THEN
                      COLLAPSED2=.FALSE.
                    ELSE IF(ITYP10(nr).GE.2) THEN !not rect. Cart.
C                     Only collapsed if r=0 (cylindrical or spherical),
C                     or lambda=0 or mu=0 (prolate or oblate).
                      IF(XP(1,1,NJ_LOC(NJL_GEOM,1,nr),MP).NE.0.0d0.AND
     '                  .(ITYP10(nr).LE.3.OR.
     '                  XP(1,1,NJ_LOC(NJL_GEOM,2,nr),MP).NE.0.0d0)) THEN
                        COLLAPSED2=.FALSE.
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !nn
                IF(COLLAPSED1.OR.COLLAPSED2) THEN
                  NFF(mf,ne)=0 !face has no surface area
                ELSE
C               Try to find a global face matching current element face
                  nf=0
                  FOUND=.FALSE.
                  DO WHILE(.NOT.FOUND.AND.nf.LT.MIN(NFT,NFM))
                    nf=nf+1
                    IF(NPF(5,nf).EQ.1) THEN !faces adjoining only 1 elem
                      nb=NBJF(1,nf) !is basis type # of face for nj=1
                      IF(NNT(nb).EQ.NNF(0,mf,MB)) THEN
                        ALL=.TRUE.
                        nee=NPF(6,nf)
                        nef=NPF(8,nf)
                        nbe=NBJ(1,nee)
C                   Old 7Nov97
C                    CALL CALC_FACE_INFORMATION_IND(NBJ(1,nee),
C     '                NBJF(1,nf),nef,NKE(1,1,1,nee),NKEF,NKF,NNF,
C     '                NPNE(1,1,nee),NPNF,nr,NVJE(1,1,1,nee),NVJF,
C     '                SE(1,1,nee),SF,ERROR,*9999)
                        nn=0
                        DO WHILE(nn.LT.NNT(nb).AND.ALL)
                          nn=nn+1
C                    np=NPNF(nn,nb)
                          nne=NNF(1+nn,nef,nbe)
                          np=NPNE(nne,nbe,nee)
                          MM=NNF(1+nn,mf,MB)
                          MP=NPNE(MM,MB,ne)
                          IF(mp.NE.np) ALL=.FALSE.
                        ENDDO
                        IF(ALL) THEN
                          NFF(mf,ne)=nf
                          NPF(5,nf)=2
                          NPF(7,nf)=ne
                          NPF(9,nf)=mf
                          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                            WRITE(OP_STRING,'('' ne='',I5,'' mf='',I2,'
     '                        //''' nft='',I5,'' NPF:'',9I6)')
     '                        ne,mf,nf,(NPF(j,nf),j=1,9)
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                          ENDIF
                          FOUND=.TRUE.
                          IF(NRE(nee).NE.nr) THEN
                            NFFACE(0,nr)=NFFACE(0,nr)+1
                            IF(NFFACE(0,nr).LE.NF_R_M) THEN
                              NFFACE(NFFACE(0,nr),nr)=nf
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO !nf
                  IF(.NOT.FOUND) THEN
                    NFT=NFT+1
                    NFFACE(0,nr)=NFFACE(0,nr)+1
                    IF(NFT.LE.NFM.AND.NFFACE(0,nr).LE.NF_R_M) THEN
                      NFFACE(NFFACE(0,nr),nr)=NFT
                      NFF(mf,ne)=NFT
C                   Set face interpolation types
                      NPF(1,NFT)=ni2
                      IF(IBT(1,ni2,MB).EQ.1) THEN
                        NPF(2,NFT)=IBT(2,ni2,MB)
                      ELSE IF(IBT(1,ni2,MB).EQ.2) THEN
                        IF(IBT(2,ni2,MB).EQ.1) THEN
                          NPF(2,NFT)=4
                        ELSE
                          NPF(2,NFT)=4+IBT(2,ni2,MB)
                        ENDIF
                      ELSE IF(IBT(1,ni2,MB).EQ.5.OR.IBT(1,ni2,MB).EQ.6)
     '                    THEN
                        NPF(2,NFT)=IBT(2,ni2,MB)
                      ENDIF
                      NPF(3,NFT)=ni3
                      IF(IBT(1,ni3,MB).EQ.1) THEN
                        NPF(4,NFT)=IBT(2,ni3,MB)
                      ELSE IF(IBT(1,ni3,MB).EQ.2) THEN
                        IF(IBT(2,ni3,MB).EQ.1) THEN
                          NPF(4,NFT)=4
                        ELSE
                          NPF(4,NFT)=4+IBT(2,ni3,MB)
                        ENDIF
                      ELSE IF(IBT(1,ni3,MB).EQ.5.OR.IBT(1,ni3,MB).EQ.6)
     '                    THEN
                        NPF(4,NFT)=IBT(2,ni3,MB)
                      ENDIF
                      DO njj1=1,3
                        DO njj2=1,NJ_LOC(njj1,0,nr)
                          nj=NJ_LOC(njj1,njj2,nr)
                          mb2=NBJ(nj,ne) !is basis # of parent elem for nj
                          IF(mb2.NE.0) THEN
                            CALL FIND_FACE_BASIS(IBT,mb2,nb,mf,NNF,
     '                        ERROR,*9999)
                            IF(nb.EQ.0) THEN
                              WRITE(ERROR,'('' >>No face basis found'
     '                          //' for face '',I5,'' of basis '',I5)'
     '                          ) NFT,mb2
                              GOTO 9999
                            ENDIF
C                         Have found basis nb that matches the face basis
                            NBJF(nj,NFT)=nb
                            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                          call mp_setlock()
                              WRITE(OP_STRING,
     '                          '('' Found face basis nb='','//
     '                          'I2,'' NNT(nb)='',I2,
     '                          '//''' NKT(nn=0..,nb)='',10I2)') nb,
     '                          NNT(nb),(NKT(nn,nb),nn=0,NNT(nb))
                              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                          call mp_unsetlock()
                            ENDIF
                          ELSE !mb2=0 so no basis defined for this nj
                            NBJF(nj,NFT)=0
                          ENDIF
                        ENDDO !njj2
                      ENDDO !njj1
                      NPF(5,NFT)=1
                      NPF(6,NFT)=ne
                      NPF(7,NFT)=0
                      NPF(8,NFT)=mf
                      NPF(9,NFT)=0
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                        WRITE(OP_STRING,'('' ne='',I5,'' mf='',I2,'
     '                    //''' nft='',I5,'' NPF:'',9I6)') ne,mf,NFT,
     '                    (NPF(j,NFT),j=1,9)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                      ENDIF
                    ENDIF !NFT <= NFM
                  ENDIF !not found
                ENDIF !face collapsed
              ENDIF !3D simplex
            ENDDO !mf
          ENDIF ! 2d/3d elements
        ENDDO !noelem
        
        WRITE(ERROR_STRING,'(''>> Increase NF_R_M, try '',I6)') 
     '    NFFACE(0,nr)
        CALL ASSERT(NFFACE(0,nr).LE.NF_R_M,ERROR_STRING,ERROR,*9999)

      ENDDO !nr


      CALL ASSERT(NFT.LE.NFM,'>>Increase NFM',ERROR,*9999)

C     Set up NLF(naf,nf) = line #s of local arcs nnf of face nf

      DO nf=1,NFT
        ne=NPF(6,nf)
        nef=NPF(8,nf)
        nr=NRE(ne)
        CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf),nef,
     '    NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,nr,
     '    NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
        nb=NBJF(1,nf) !basis fn for first geometric variable
        DO naf=1,NLE(nb)
          nl=0
          FOUND=.FALSE.
          DO WHILE(.NOT.FOUND.AND.nl.LT.NLT)
            nl=nl+1
            FOUND=.TRUE.
            IF(NIT(nb).EQ.2.AND.IBT(1,1,nb).EQ.3) THEN ! 2D Simplex
              IF(NPNF(NNL(1,naf,nb),nb).NE.NPL(2,1,nl).AND.
     '           NPNF(NNL(2,naf,nb),nb).NE.NPL(2,1,nl))
     '           FOUND=.FALSE.
              IF(NPNF(NNL(1,naf,nb),nb).NE.NPL(3,1,nl).AND.
     '           NPNF(NNL(2,naf,nb),nb).NE.NPL(3,1,nl))
     '           FOUND=.FALSE.
            ELSE
              DO n1=1,NNL(0,naf,nb)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
               WRITE(OP_STRING,'(/'' nf='',I5,'' NAF='',I2,'' nl='',I5,'
     '              //''' N1='',I1,'' NNL='',I2,'' NPL(N1+1,nj=1,nl)='',
     '              I2,'' nb='',I2)') nf,naf,nl,n1,NNL(n1,naf,nb),
     '              NPL(n1+1,1,nl),nb
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
                ENDIF
                IF(NNL(n1,naf,nb).NE.0) THEN
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                    WRITE(OP_STRING,'('' NPNF='',I5)')
     '                NPNF(NNL(n1,naf,nb),nb)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                  ENDIF
                  IF(NPNF(NNL(n1,naf,nb),nb).NE.NPL(n1+1,1,nl))
     '              FOUND=.FALSE.
                ENDIF
              ENDDO !n1
            ENDIF
          ENDDO !nl
          IF(FOUND) THEN
            NLF(naf,nf)=nl
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(/'' For face nf='',I5,'
     '          //''' local arc NAF='',I2,'
     '          //''' is global arc NLF(naf,nf)='',I5)') nf,naf,nl
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ELSE
            NLF(naf,nf)=0
          ENDIF
        ENDDO !naf
      ENDDO !nf

      CALL EXITS('FACSEG')
      RETURN
 9999 CALL ERRORS('FACSEG',ERROR)
      CALL EXITS('FACSEG')
      RETURN 1
      END


