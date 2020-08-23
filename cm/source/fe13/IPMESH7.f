      SUBROUTINE IPMESH7(NBJ,NEELEM,NENP,NKJE,NKJ,
     '  NPNE,NPNODE,NPQ,NQLIST,nr,NRE,NVJE,NVJP,NXQ,
     '  SE,XP,XQ,ERROR,*)

C#### Subroutine: IPMESH7
C###  Description:
C###    IPMESH7 defines mesh parameters for grid point mesh.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh01.cmn'
      INCLUDE 'mxch.inc'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),NPQ(NQM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NQLIST(0:NQM),
     '  nr,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,
     '  nb,nb1,nc,ne,nj,nk,nn,nnp,ns,noelem,
     '  nonode,NOQUES,np,nq,nq1,nq2,nq3,nq4,nqN,nqS,nq_BRANCH,
     '  n1list
      LOGICAL BRANCH,FILEIP,INLIST,NEIGHBOUR_IN_LIST,ON_SECOND_BRANCH,
     '  SAMENUMXI,SAMETYPE

      CALL ENTERS('IPMESH7',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
! Initialising BRANCH MLB 4/4/97
      BRANCH=.FALSE.

      FORMAT='($,'' Enter basis function# for mesh [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=MESH1_nb
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) MESH1_nb=IDATA(1)
      nb=MESH1_nb

! Store grid# corresponding to node#
      noelem=0
      nonode=0
      ne = NEELEM(0,0) !initial element#
      np = NPNODE(0,0) !initial node#
      nq = 1           !initial grid#
      nqS= nq          !stores grid# at start of row
      NQLIST(0)=0      !initialise #nq's in NQLIST
      ON_SECOND_BRANCH=.FALSE.

! Main grid point loop
      DO WHILE(nqS.GT.0)
        nonode=nonode+1 !increment node counter
        np=np+1         !increment node#
        CALL ASSERT(np.LE.NPM,'>>Too many nodes. Increase NPM',
     '    ERROR,*9999)
C GMH 8/1/97 Update cmgui link
        CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
        NPNODE(nonode,nr)=np
        NPQ(nq)=np           !store node# at current grid#
        NQLIST(0)=NQLIST(0)+1
        NQLIST(NQLIST(0))=nq !store nq in current list
        nq4=nq
        nq3=NXQ(-1,1,nq4,1) !is neighbour in -ve Xi1
        nq2=NXQ(-2,1,nq4,1) !is neighbour in -ve Xi2
        nq1=NXQ(-1,1,nq2,1) !is neighbour in -ve Xi1 & -ve Xi2

        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' NQLIST(0)='',I6)') NQLIST(0)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' np='',I6,11X,'' nq1..='',4I6,'
     '      //''' *** nq'')') np,nq1,nq2,nq3,nq4
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !dop

! Define new element if all nodes are defined
        IF(nq1.GT.0.AND.nq2.GT.0.AND.nq3.GT.0.AND.nq4.GT.0) THEN
          IF(NPQ(nq1).GT.0.AND.NPQ(nq2).GT.0.AND
     '      .NPQ(nq3).GT.0.AND.NPQ(nq4).GT.0) THEN
            noelem=noelem+1 !increment element counter
            ne=ne+1         !increment element#
            CALL ASSERT(ne.LE.NEM,'>>Too many elements. Increase NEM',
     '        ERROR,*9999)
            NEELEM(noelem,nr)=ne
            NPNE(1,nb,ne)=NPQ(nq1)  !1st node of element ne
            NPNE(2,nb,ne)=NPQ(nq2)  !2nd node of element ne
            NPNE(3,nb,ne)=NPQ(nq3)  !3rd node of element ne
            NPNE(4,nb,ne)=NPQ(nq4)  !4th node of element ne
            IF(DOP) THEN
              WRITE(OP_STRING,'('' ne='',I6,'' NPNE(1..4,nb,ne)='',4I6,'
     '          //''' *** new element'' )') ne,(NPNE(nn,nb,ne),nn=1,4)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF !dop
          ENDIF !npq
        ENDIF !nq1..4

! Define coords of new node
        DO nj=1,NJT
          XP(1,1,nj,np)=XQ(nj,nq)
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,'('' XP(1,1,nj,'',I5,''): '',3E12.3)')
     '      np,(XP(1,1,nj,np),nj=1,NJT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !dop

! Find neighbouring grid point
        nqN=NXQ(1,1,nq,1) !is neighbour in +ve Xi1

! Check for wrapped line (i.e. nqN in list of current grid pts)
        NEIGHBOUR_IN_LIST=INLIST(nqN,NQLIST(1),NQLIST(0),n1list)

! If wrapped, define new element (if all nodes OK) at end of row
        IF(NEIGHBOUR_IN_LIST.AND..NOT.ON_SECOND_BRANCH) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,'('' neighbour is in list'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !dop
          nq4=nqN
          nq3=NXQ(-1,1,nq4,1) !is neighbour in -ve Xi1
          nq2=NXQ(-2,1,nq4,1) !is neighbour in -ve Xi2
          nq1=NXQ(-2,1,nq3,1) !is neighbour in -ve Xi1 & -ve Xi2
          IF(nq1.GT.0.AND.nq2.GT.0.AND.nq3.GT.0.AND.nq4.GT.0) THEN
            IF(NPQ(nq1).GT.0.AND.NPQ(nq2).GT.0.AND
     '        .NPQ(nq3).GT.0.AND.NPQ(nq4).GT.0) THEN
              noelem=noelem+1 !increment element counter
              ne=ne+1         !increment element#
              NEELEM(noelem,nr)=ne
              NPNE(1,nb,ne)=NPQ(nq1)  !1st node of element ne
              NPNE(2,nb,ne)=NPQ(nq2)  !2nd node of element ne
              NPNE(3,nb,ne)=NPQ(nq3)  !3rd node of element ne
              NPNE(4,nb,ne)=NPQ(nq4)  !4th node of element ne
              IF(DOP) THEN
                WRITE(OP_STRING,'( '' np='',I6,10X,'' nq1..4='',4I6,'
     '            //''' *** nq'')') np,nq1,nq2,nq3,nq4
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' ne='',I6,'' NPNE(1..4,nb,ne)='','
     '            //'4I6,'' *** new element at row end'' )')
     '            ne,(NPNE(nn,nb,ne),nn=1,4)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF !dop
            ENDIF !npq
          ENDIF !nq1..4,npq
        ENDIF !in_list

! Move to nqN if nq is not a bdry pt and mesh is not wrapped
        IF(nqN.GT.0.AND..NOT.NEIGHBOUR_IN_LIST) THEN !move to nqN
          IF(NXQ(1,0,nq,1).GT.1) THEN
            BRANCH=.TRUE. !flag branch in elements
            nq_BRANCH=nq  !store grid pt at branch
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Branch at nq='',I6)') nq
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF !dop
          ENDIF !nxq
          nq=nqN

! Else move to another branch or start of next row
        ELSE
          IF(BRANCH) THEN
            ON_SECOND_BRANCH=.TRUE.
            nq=NXQ(1,2,nq_BRANCH,1) !is 2nd neighbour in +ve Xi1
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Go to 2nd branch at nq='',I6)') nq
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF !dop
          ELSE
            ON_SECOND_BRANCH=.FALSE.
            nq=NXQ(2,1,nqS,1) !is grid# above nqS
            nqS=nq            !reset start nq
            NQLIST(0)=0       !reinitialise #nq's in NQLIST
          ENDIF !branch
          BRANCH=.FALSE.      !reset branch
        ENDIF !nqN>0

      ENDDO !while nqs>0

! Fill out various element and node arrays
      NEELEM(0,nr)=noelem !#elements in region nr
      NPNODE(0,nr)=nonode !#nodes in region nr

      NET(nr)=ne !highest element# in region nr
      NPT(nr)=np !highest node# in region nr
      NET(0) =ne !highest element# in all regions
      NPT(0) =np !highest node# in all regions
      NEELEM(0,0)=NEELEM(0,0)+NEELEM(0,nr) !#elements in all regions
      NPNODE(0,0)=NPNODE(0,0)+NPNODE(0,nr) !#nodes in all r

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          NBJ(nj,ne)=MESH1_nb
          DO nn=1,NNT(nb)
            DO nk=1,NKT(nn,nb)
              NKJE(nk,nn,nj,ne)=nk
            ENDDO !nk
          ENDDO !nn
        ENDDO
        NRE(ne)=nr
      ENDDO !noelem

      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        DO nj=1,NJT
          NKJ(nj,np)=NKT(0,MESH1_nb)
          DO nc=1,NCM
            NVJP(nj,np)=1 !one version per nj per node
          ENDDO !nc
        ENDDO !nj
      ENDDO !nonode (np)

      DO nb1=1,NBFT
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO ns=1,NST(nb1)+NAT(nb1)
            SE(ns,nb1,ne)=1.0D0
          ENDDO !ns
          DO nn=1,NNT(nb1)
C KAT 23Feb01: now handled by NKJE above
C            DO nk=1,NKT(nn,nb1)
C              NKE(nk,nn,nb1,ne)=NK
C            ENDDO !nk
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              NVJE(nn,nb1,nj,ne)=1 !version one of nn,nj in elem ne
            ENDDO !nj
          ENDDO !nn

! Update alternate bases
          SAMETYPE=.FALSE.
          SAMENUMXI=.FALSE.
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            IF(NBC(nb1).EQ.NBC(NBJ(nj,ne)).OR.NBC(nb1).EQ.7)
     '        SAMETYPE=.TRUE. !Same basis type or extended basis
            IF(NIT(nb1).EQ.NIT(NBJ(nj,ne))) SAMENUMXI=.TRUE.
          ENDDO
          IF(NNT(nb1).GT.0.AND.SAMETYPE.AND.SAMENUMXI.AND.
     '      nb1.NE.nb) THEN
            DO nnp=1,8
              NPNE(nnp,nb1,ne)=NPNE(nnp,nb,ne)
            ENDDO
          ENDIF
        ENDDO !noelem (ne)
      ENDDO !nb1

      CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

      CALL EXITS('IPMESH7')
      RETURN
 9999 CALL ERRORS('IPMESH7',ERROR)
      CALL EXITS('IPMESH7')
      RETURN 1
      END


