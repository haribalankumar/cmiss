      SUBROUTINE IPMESH8_DYNAM(AMAP,AREANUM,BMAX,ELEMLIST,GEOMFRAC,
     '  IBT,IDO,INP,MAX_BRANCH,MAX_RAN_PTS,NBJ,NEELEM,NKJE,NKJ,
     '  NPF,NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,
     '  POINTS,POINTCOORD,PTAR,PWEIGHT,RANDOM_COORD,SE,START,
     '  STARTNODE,TCX,TMAX,WEIGHTS,XA,XE,XJPOWER,XP,ERROR,*)

C#### Subroutine: IPMESH8_DYNAM
C###  Description:
C###    defines a tree using the Monte-Carlo method. It is possible
C###    to grow a stand alone tree or grow a tree over a host FE
C###    mesh.
C***  Martin Buist December 1996

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'

!  Parameter list
      INTEGER BMAX,IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),MAX_BRANCH,MAX_RAN_PTS,NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),TMAX
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!  Local variables
      INTEGER AMAP(MAX_BRANCH),AREANUM(MAX_RAN_PTS),ESTART,
     '  ELEMLIST(MAX_BRANCH,4),GENERATIONS,ICHAR,INFO,LINETYPE,nb,nb1,
     '  nb2,nba,nbg,ne,ne1,nee,nj,nk,nlms,nn,nnp,NOQUES,nonode,noelem,
     '  NOELEMS(0:6),np,npp,np1,NPSTART,NPFIN,NPFINISH,nr2,
     '  nrc,NRANDOM,ns,NTREES,ntt,NUMAREAS,NUMELEMI,NUMPTS,
     '  PTAR(MAX_BRANCH),STARTNODE(TMAX,2),TCX(0:NP_R_M,2,0:BMAX),
     '  TOTALAREAS,NJT_TEMP
      REAL*8 BLENGTH,BSIZE(3),COEF(2),COFM(3),
     '  GEOMFRAC(MAX_BRANCH,2,3),HOLD,NORML(4),NUMELEM,
     '  POINT1(3),POINT2(3),POINT3(3),POINT11(2),POINT22(2),
     '  POINTCOORD(MAX_BRANCH,3),PWEIGHT(MAX_RAN_PTS),PXI,
     '  RANDOM_COORD(MAX_RAN_PTS,NJT),S1,S2,SBSIZE,START(TMAX,3),
     '  WEIGHTS(MAX_BRANCH),X(3),XI(3),XJPOWER(MAX_BRANCH,3),Z(3)
      LOGICAL BRANCH,ELEM,ELEMSMOOTH,FILEIP,GEOM,SAMETYPE,
     '  SAMENUMXI,SURF_TREE,POINTS(MAX_BRANCH)

      CALL ENTERS('IPMESH8_DYNAM',*9999)

      CALL ASSERT((NJT.EQ.2).OR.(NJT.EQ.3),
     '  '>>ERROR can only grow in 2 or 3 dimensions',ERROR,*9999)

C Initialisation
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      nb=1
      BLENGTH=0.65d0
      SBSIZE=0.5d0
      DO nj=1,3
        BSIZE(nj)=0.0d0
      ENDDO !nj

      CALL TREE_START(GEOM,IBT,IDO,INP,NBJ,NEELEM,
     '  NKJE,NPF,NP_INTERFACE,NPNE,nr2,
     '  NTREES,NVJE,SE,START,STARTNODE,SURF_TREE,
     '  TMAX,XA,XE,XP,ERROR,*9999)

      IF(SURF_TREE) THEN
        NJT_TEMP=NJT-1
      ELSE
        NJT_TEMP=NJT
      ENDIF

      DO ntt=1,NTREES
        CALL TREE_GEOM(ELEM,ELEMLIST,ELEMSMOOTH,ESTART,
     '    GENERATIONS,GEOMFRAC,MAX_BRANCH,NBJ,NOELEMS,
     '    NXI,ntt,POINTS,WEIGHTS,XJPOWER,ERROR,*9999)

C Enter the basis function type for the trees
        IDEFLT(1)=1
        FORMAT='($,'' Enter the tree basis function type [1]:'',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=nb
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) nb=IDATA(1)

        CALL TREE_POINTS(ELEM,ELEMLIST,GEOMFRAC,MAX_BRANCH,MAX_RAN_PTS,
     '    NOELEMS,NRANDOM,ntt,POINTS,PWEIGHT,RANDOM_COORD,XJPOWER,
     '    WEIGHTS,ERROR,*9999)

C Enter the fraction of the distance to the centre of mass the tree
C will grow to at each generation
        RDEFLT(1)=0.65d0
        FORMAT='($,'' Enter the branch length ratio [0.65]: '',F8.6)'
        IF(IOTYPE.EQ.3) RDATA(1)=BLENGTH
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) BLENGTH=RDATA(1)

C Enter the smallest branch which can be made in real space
        RDEFLT(1)=0.5d0
        FORMAT='($,'' Enter smallest allowable branch [0.5]: '',F8.6)'
        IF(IOTYPE.EQ.3) RDATA(1)=SBSIZE
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SBSIZE=RDATA(1)

C Set up the initial states of the branch arrays
        DO nrc=1,NRANDOM
          AREANUM(nrc)=1  !all points in area 1
        ENDDO !nrc

        DO nba=1,MAX_BRANCH
          AMAP(nba)=1
          PTAR(nba)=1  !all branches in area 1
          DO nj=1,NJT_TEMP
            POINTCOORD(nba,nj)=0.0d0
          ENDDO !nj
        ENDDO !nba

        DO nba=0,NP_R_M
          DO nbg=0,BMAX
            TCX(nba,1,nbg)=0
            TCX(nba,2,nbg)=0
          ENDDO !nbg
        ENDDO !nba

        NUMPTS=1
        NUMAREAS=1
        TOTALAREAS=1

C Set the start of the current tree
        IF(ELEM) THEN
          IF(GEOM) THEN
C            DO nj=1,NJT_TEMP
C              POINTCOORD(1,nj)=START(ntt,nj)+
C     '          ELEMLIST(STARTNODE(ntt,2),nj+1)
C            ENDDO !nj
            DO ne=1,NOELEMS(0)
              IF(ELEMLIST(ne,1).EQ.STARTNODE(ntt,2)) THEN
                DO nj=1,NJT_TEMP
                  POINTCOORD(1,nj)=START(ntt,nj)+
     '              ELEMLIST(ne,nj+1)
                ENDDO !nj
              ENDIF
            ENDDO !ne
          ELSE
            DO ne=1,NOELEMS(0)
              IF(ELEMLIST(ne,1).EQ.ESTART) THEN
                DO nj=1,NJT_TEMP
                  POINTCOORD(1,nj)=START(ntt,nj)+
     '              ELEMLIST(ne,nj+1)
                ENDDO !nj
              ENDIF
            ENDDO !ne
          ENDIF
        ELSE
          DO nj=1,NJT_TEMP
            POINTCOORD(1,nj)=START(ntt,nj)
          ENDDO !nj
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'('' START (TRANS) '',F8.4,F8.4,F8.4)')
     '      POINTCOORD(1,1),POINTCOORD(1,2),POINTCOORD(1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

C Loop over the number of generations required
        DO nbg=1,GENERATIONS
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' '
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Generation number '',I4)') nbg
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*) ' '
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          DO nba=1,NUMAREAS
            BRANCH=.TRUE.
            DO nj=1,NJT_TEMP
              COFM(nj)=0.0d0
            ENDDO !nj

            HOLD=0.0d0
            DO nrc=1,NRANDOM
              IF(AREANUM(nrc).EQ.nba) THEN
                DO nj=1,NJT_TEMP
                  COFM(nj)=COFM(nj)+(RANDOM_COORD(nrc,nj)*
     '              PWEIGHT(nrc))
                ENDDO !nj
                HOLD=HOLD+PWEIGHT(nrc)
              ENDIF
            ENDDO !nrc

            IF(HOLD.GT.ZERO_TOL) THEN
              DO nj=1,NJT_TEMP
                COFM(nj)=COFM(nj)/HOLD
              ENDDO !nj
            ELSE
              WRITE(OP_STRING,'(''>>WARNING: no point point found. '')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO nj=1,NJT_TEMP
                COFM(nj)=0.0d0
              ENDDO !nj
              BRANCH=.FALSE.
            ENDIF

            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' Center of mass '',F12.4,F12.4,F12.4)')
     '          COFM(1),COFM(2),COFM(3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(NJT_TEMP.EQ.3) THEN
              DO nj=1,NJT_TEMP
                POINT1(nj)=COFM(nj)
                POINT2(nj)=POINTCOORD(PTAR(nba),nj)
              ENDDO !nj
              IF(nbg.EQ.1) THEN !branch to split along shortest axis
                DO nj=1,NJT_TEMP
                  POINT3(nj)=COFM(nj)
                ENDDO !nj
                IF((NOELEMS(1).LE.NOELEMS(2)).AND.
     '            (NOELEMS(1).LE.NOELEMS(3))) THEN
                  POINT3(1)=POINT3(1)+1.0d0
                ELSE IF((NOELEMS(2).LE.NOELEMS(3)).AND.
     '            (NOELEMS(2).LE.NOELEMS(1))) THEN
                  POINT3(2)=POINT3(2)+1.0d0
                ELSE
                  POINT3(3)=POINT3(3)+1.0d0
                ENDIF
              ELSE
                DO nj=1,NJT_TEMP
                  POINT3(nj)=POINTCOORD(TCX(PTAR(nba),1,0),nj)
                ENDDO !nj
              ENDIF
              CALL PLANE_FROM_3_PTS(NORML,2,POINT1,POINT2,POINT3,
     '          ERROR,*9999)
            ELSE
              DO nj=1,NJT_TEMP
                POINT11(nj)=COFM(nj)
                POINT22(nj)=POINTCOORD(PTAR(nba),nj)
              ENDDO !nj
              CALL LINE_FROM_2_PTS(COEF,POINT11,POINT22,LINETYPE,
     '          ERROR,*9999)
            ENDIF

            TOTALAREAS=TOTALAREAS+1
            NUMPTS=NUMPTS+1

C Create the new branch
            DO nj=1,NJT_TEMP
              POINTCOORD(NUMPTS,nj)=(POINTCOORD(PTAR(nba),nj)*
     '          (1.0d0-BLENGTH))+(COFM(nj)*BLENGTH)
              BSIZE(nj)=POINTCOORD(NUMPTS,nj)-POINTCOORD(PTAR(nba),nj)
            ENDDO !nj

            BSIZE(1)=DSQRT(BSIZE(1)**2+BSIZE(2)**2+BSIZE(3)**2)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Branch size '',F12.4)') BSIZE(1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
C Check to see the branch length is above the allowable tolerance
            IF(BSIZE(1).LT.SBSIZE) THEN
              BRANCH=.FALSE.
              NUMPTS=NUMPTS-1
              TOTALAREAS=TOTALAREAS-1
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Branch discarded'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF

C Update connectivity array, area number arrays
            IF(BRANCH) THEN
              IF(NJT_TEMP.EQ.3) THEN
                DO nrc=1,NRANDOM
                  IF(AREANUM(nrc).EQ.nba) THEN
                    IF((RANDOM_COORD(nrc,1)*NORML(1))
     '                +(RANDOM_COORD(nrc,2)*NORML(2))
     '                +(RANDOM_COORD(nrc,3)*NORML(3))
     '                +NORML(4).GE.0.0d0) THEN
                      AREANUM(nrc)=TOTALAREAS
                    ELSE
                      AREANUM(nrc)=nba
                    ENDIF
                  ENDIF
                ENDDO !nrc
              ELSE
                DO nrc=1,NRANDOM
                  IF(AREANUM(nrc).EQ.nba) THEN
                    IF(LINETYPE.EQ.1) THEN
                      IF(RANDOM_COORD(nrc,2).GE.
     '                  ((RANDOM_COORD(nrc,1)*COEF(1))+COEF(2))) THEN
                        AREANUM(nrc)=TOTALAREAS
                      ELSE
                        AREANUM(nrc)=nba
                      ENDIF
                    ELSE
                      IF(RANDOM_COORD(nrc,1).GE.COEF(1)) THEN
                        AREANUM(nrc)=TOTALAREAS
                      ELSE
                        AREANUM(nrc)=nba
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !nrc
              ENDIF

              AMAP(nba)=TOTALAREAS
              AMAP(TOTALAREAS)=nba

              IF(ELEMSMOOTH) THEN
C Make all the elements around the same length
                CALL ASSERT(SBSIZE.GT.ZERO_TOL,
     '            '>>ERROR Can not use smoothing with zero branch',
     '            ERROR,*9999)
                NPFIN=NUMPTS
                NUMELEM=DNINT(BSIZE(1)/SBSIZE)
                NUMELEMI=INT(NUMELEM)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Split into # elements '',I4)')
     '              NUMELEMI
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                DO nj=1,NJT_TEMP
                  POINTCOORD(NUMPTS+NUMELEMI-1,nj)=
     '              POINTCOORD(NUMPTS,nj)
                ENDDO !nj
                DO nlms=1,(NUMELEMI-1)
                  S1=DBLE(NUMELEMI-nlms)
                  S2=DBLE(nlms)
                  DO nj=1,NJT_TEMP
                    POINTCOORD(NPFIN,nj)=((S1/NUMELEM)*
     '                POINTCOORD(PTAR(nba),nj))
     '                +((S2/NUMELEM)*POINTCOORD(NUMPTS+NUMELEMI-1,nj))
                  ENDDO !nj
                  NPFIN=NPFIN+1
                ENDDO !nlms
                TCX(NUMPTS,1,0)=PTAR(nba)
                TCX(PTAR(nba),2,0)=TCX(PTAR(nba),2,0)+1
                TCX(PTAR(nba),2,TCX(PTAR(nba),2,0))=NUMPTS
                IF(NUMELEM.GT.1.5d0) THEN
                  DO nlms=NUMPTS+1,NPFIN
                    TCX(nlms,1,0)=TCX(nlms-1,1,0)
                    TCX(nlms-1,2,0)=TCX(nlms-1,2,0)+1
                    TCX(nlms-1,2,TCX(nlms-1,2,0))=nlms
                  ENDDO !nlms
                ENDIF
                NUMPTS=NPFIN
                PTAR(nba)=NUMPTS
                PTAR(TOTALAREAS)=NUMPTS
              ELSE
C Make one element for each branch
                TCX(NUMPTS,1,0)=PTAR(nba)
                TCX(PTAR(nba),2,0)=TCX(PTAR(nba),2,0)+1
                TCX(PTAR(nba),2,TCX(PTAR(nba),2,0))=NUMPTS
                PTAR(nba)=NUMPTS
                PTAR(TOTALAREAS)=NUMPTS
              ENDIF !smoothing
            ENDIF !branch

          ENDDO !nba
          NUMAREAS=TOTALAREAS
        ENDDO !nbg

        IF(ELEM) THEN
          IF(DOP) THEN
            DO np=1,NUMPTS
              WRITE(OP_STRING,'('' NP '',I4,F12.6,F12.6,F12.6)')
     '          np,POINTCOORD(np,1),POINTCOORD(np,2),POINTCOORD(np,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !np
          ENDIF
          DO np=1,NUMPTS
            IF(SURF_TREE) POINTCOORD(np,3)=0.0d0
            ne=ELEMLIST(1,1)
            nb1=NBJ(1,ne)
            DO ne1=1,NOELEMS(0)
              DO nj=1,NJT
                IF(POINTCOORD(np,nj).GT.1.0d0) THEN
                  ne=NXI(nj,1,ne)
                  POINTCOORD(np,nj)=POINTCOORD(np,nj)-1.0d0
                ELSE IF(POINTCOORD(np,nj).LT.0.0d0) THEN
                  ne=NXI(-nj,1,ne)
                  POINTCOORD(np,nj)=POINTCOORD(np,nj)+1.0d0
                ENDIF
              ENDDO !nj
            ENDDO !ne1
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '        nr2,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,
     '        ERROR,*9999)
            DO nj=1,NJT
              XI(nj)=POINTCOORD(np,nj)
            ENDDO !nj
            DO nj=1,NJT
              nb2=NBJ(nj,ne)
              X(nj)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),INP(1,1,nb1),
     '          nb2,1,XI,XE(1,nj))
            ENDDO !nj
            CALL XZ(ITYP10(nr2),X,Z)
            DO nj=1,NJT
              POINTCOORD(np,nj)=Z(nj)
            ENDDO !nj
          ENDDO !np
          IF(DOP) THEN
            DO np=1,NUMPTS
              WRITE(OP_STRING,'('' NP '',I4,F12.6,F12.6,F12.6)')
     '          np,POINTCOORD(np,1),POINTCOORD(np,2),POINTCOORD(np,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !np
          ENDIF
        ENDIF

C Write the tree into nodes and elements in CMISS format
        NPSTART=NPT(0)+1
        NPFINISH=NPT(0)+NUMPTS

        IF(GEOM) THEN
          DO np=NPSTART,NPFINISH-1
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            DO nj=1,NJT
              XP(1,1,nj,np)=POINTCOORD((np-NPT(0))+1,nj)
            ENDDO !nj
          ENDDO !np
        ELSE
C Nodal coordinates
          DO np=NPSTART,NPFINISH
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            DO nj=1,NJT
              XP(1,1,nj,np)=POINTCOORD((np-NPT(0)),nj)
            ENDDO !nj
          ENDDO !np
        ENDIF

        IF(GEOM) THEN
          ne1=1
        ELSE
          ne1=0
        ENDIF

C Element numbers
        ne=NET(0)
        DO np=1,NUMPTS
          IF(GEOM.AND.(np.EQ.1)) THEN
            IF(TCX(np,2,0).GT.0) THEN
              DO nba=1,TCX(np,2,0)  !loop over number of children
                ne=ne+1
                NPNE(1,nb,ne)=STARTNODE(ntt,1)
                NPNE(2,nb,ne)=TCX(np,2,nba)+NPT(0)-ne1
              ENDDO !nba
            ENDIF
          ELSE
            IF(TCX(np,2,0).GT.0) THEN
              DO nba=1,TCX(np,2,0)  !loop over number of children
                ne=ne+1
                NPNE(1,nb,ne)=np+NPT(0)-ne1
                NPNE(2,nb,ne)=TCX(np,2,nba)+NPT(0)-ne1
              ENDDO !nba
            ENDIF
          ENDIF
        ENDDO !np

        IF(GEOM) THEN
          NPFINISH=NPFINISH-1
        ENDIF

C Sizing arrays
        ne1=NEELEM(0,nr)+1
        np1=NPNODE(0,nr)+1
        NEELEM(0,nr)=ne-NET(0)+NEELEM(0,nr)
        nee=NET(0)
        NET(nr)=ne                  !highest element# in region nr
        NET(0)=ne                   !highest element# in all regions
                                    !number elements in region nr
        NEELEM(0,0)=ne              !number elements in all regions
        npp=NPT(0)
        NPT(nr)=NPFINISH            !highest node# in region nr
        NPT(0)=NPFINISH             !highest node# in all regions
        NPNODE(0,nr)=NPNODE(0,nr)+NPFINISH-NPSTART+1
                                  !number nodes in region nr
        NPNODE(0,0)=NPFINISH        !number nodes in all regions

C Other necessary arrays
        DO noelem=ne1,NEELEM(0,nr)
          ne=noelem+nee-ne1+1
          NEELEM(noelem,nr)=NE
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            NBJ(nj,ne)=nb
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                NKJE(nk,nn,nj,ne)=nk
              ENDDO !nk
            ENDDO !nn
          ENDDO
          NRE(ne)=nr
        ENDDO

        DO nonode=np1,NPNODE(0,nr)
          np=nonode+npp-np1+1
          NPNODE(nonode,nr)=np
          DO nj=1,NJT
            NKJ(nj,np)=NKT(0,nb) !one version per nj per node
            NVJP(nj,np)=1
          ENDDO !nj
        ENDDO !nonode (np)

        DO nb1=1,NBFT
          DO noelem=ne1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ns=1,NST(nb1)+NAT(nb1)
              SE(ns,nb1,ne)=1.0D0
            ENDDO !ns
            DO nn=1,NNT(nb1)
C KAT 23Feb01: now handled by NKJE above
C              DO nk=1,NKT(nn,nb1)
C                NKE(nk,nn,nb1,ne)=NK
C              ENDDO !nk
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                NVJE(nn,nb1,nj,ne)=1 !version one of nn,nj in elem ne
              ENDDO !nj
            ENDDO !nn
C Update alternate bases
            SAMETYPE=.FALSE.
            SAMENUMXI=.FALSE.
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              IF(NBC(nb1).EQ.NBC(NBJ(nj,ne)).OR.NBC(nb1).EQ.7)
     '          SAMETYPE=.TRUE. !Same basis type or extended basis
              IF(NIT(nb1).EQ.NIT(NBJ(nj,ne))) SAMENUMXI=.TRUE.
            ENDDO
            IF(NNT(nb1).GT.0.AND.SAMETYPE.AND.SAMENUMXI.AND.
     '        nb1.NE.nb) THEN
              DO nnp=1,8
                NPNE(nnp,nb1,ne)=NPNE(nnp,nb,ne)
              ENDDO
            ENDIF
          ENDDO !noelem (ne)
        ENDDO !nb1
      ENDDO !ntt

      CALL EXITS('IPMESH8_DYNAM')
      RETURN
 9999 CALL ERRORS('IPMESH8_DYNAM',ERROR)
      CALL EXITS('IPMESH8_DYNAM')
      RETURN 1
      END


