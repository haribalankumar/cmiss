      SUBROUTINE CHMESH(IBT,IDO,INP,NBJF,NBJ,NEELEM,NEL,NELIST,NENP,NFF,
     '  NFFACE,NGAP,NKJE,NKEF,NKJ,NLF,NLL,NLLINE,NNF,NNL,
     '  NP_INTERFACE,NPF,NPL,NPLIST,NPNE,NPNF,NPNODE,NRE,NRLIST,NUNK,
     '  NVJE,NVJF,NVJL,NVJP,DF,DL,PAOPTI,PG,RG,SE,SF,WG,XA,XE,XG,XIG,XP,
     '  STRING,ERROR,*)

C#### Subroutine: CHMESH
C###  Description:
C###    CHMESH allows the user to change aspects of the mesh.
C###    Changing the surface shape of a mesh was developed as part of a
C###    Year 4 project in 1995 as an attempt to customise the torso mesh
C###    to individual people.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'chmesh0.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),
     '  NFFACE(0:NF_R_M,NRM),NGAP(NIM,NBM),NP_INTERFACE(0:NPM,0:3),
     '  NKEF(0:4,16,6,NBFM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),NRE(NEM),NUNK(NKM,NJM,NPM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),
     '  NVJP(NJM,NPM)
      REAL*8 DF(NFM),DL(3,NLM),PAOPTI(*),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INCREASE_NVM,INPI,INPT(3),N3CO,n,nb,nb_new,
     '  nb_old,ne,ni,NITB,nk,NKJCHEK,NKJNEED,NKSTART,nj,njj,
     '  NJL_TYPE,nn,NNTB,no_ne,no_np,no_nrlist,np,NPTNEW,NPTOLD,nr,
     '  ns,NTRL,nu,nv
      REAL*8 DEPTH,LLOOSE_TOL,DIFF,REAL_TEMP(1),TOTAL,WIDTH,X(11),XI(3),
     '  Z
      LOGICAL ALL_REGIONS,BASIS,CBBREV,CHANGE,EXISTS,NEWVERSION,USED,
     '  UVERSIONS
!     Functions
      INTEGER IFROMC
      REAL*8 PXI

      CALL ENTERS('CHMESH',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM change mesh
C###  Parameter:      <measurement Z>
C###  Parameter:      <change/nochange>[change]
C###  Parameter:      <region (#s/all) [1]>
C###  Description:
C###    To change the shape of the mesh to customise a generic mesh for
C###    specific cases. measurement gives a data set for a specific case.

        OP_STRING(1)=STRING(1:IEND)//' from optimiser'
        OP_STRING(2)=BLANK(1:IEND)//'<measurement Z>'
        OP_STRING(3)=BLANK(1:15)//'<change/nochange>[change]'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change mesh from optimiser
C###  Parameter:      <region (#s/all) [1]>
C###  Description:
C###    Uses an optimiser to fit a generic mesh to a specific set of
C###    data points.

        OP_STRING(1)=STRING(1:IEND)//' from optimiser'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change mesh basis BASIS#
C###  Parameter:      <(geometry/fibre/field) #>[field 1]
C###    Specify the variable.
C###  Parameter:      <element (#s/all)[all]>
C###    Specifies the elements to be changed.
C###  Parameter:      <region (#s/all) [1]>
C###    Specify the region of the mesh to be updated.
C###  Parameter:      <nodes_from NODE#[1]>
C###    Specifies the node number assigned to the next node created.  By
C###    default the node number will be the next number after the
C###    largest in the file.
C###  Description:
C###    Changes the basis used in the mesh for a specified variable,
C###    creating necessary new nodes.

        OP_STRING(1)=STRING(1:IEND)//' basis #'
        OP_STRING(2)=BLANK(1:15)//'<(geometry/fibre/field) #>[field 1]'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change mesh unique_versions
C###  Parameter:      <(geometry/fibre/field) #>[field 1]
C###    Specify the variable.
C###  Parameter:      <region (#s/all) [1]>
C###    Specify the region of the mesh to be updated.
C###  Description:
C###    Changes the mesh for a specified variable so that there is a
C###    unique version of each node for each element.  This enables the
C###    C1 continuity restriction for cubic Hermite elements to be
C###    relaxed.

        OP_STRING(1)=STRING(1:IEND)//' unique_versions'
        OP_STRING(2)=BLANK(1:15)//'<(geometry/fibre/field) #>[field 1]'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHMESH',ERROR,*9999)
      ELSE
        CALL ASSERT(NET(0).GT.0,'>>Read in a mesh first',ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        UVERSIONS=.FALSE.
        BASIS=.FALSE.
        IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
          CHANGE=.TRUE.
          DO n=1,NPC
            POLY_COEFFS(n)=PAOPTI(n)
          ENDDO
          DO n=1,NCC
            COS_COEFFS(n)=PAOPTI(NPC+n)
          ENDDO
        ELSE IF(CBBREV(CO,'UNIQUE_VERSIONS',2,noco+1,NTCO,N3CO)) THEN
          UVERSIONS=.TRUE.
        ELSE IF(CBBREV(CO,'BASIS',2,noco+1,NTCO,N3CO)) THEN
          BASIS=.TRUE.
          nb_new=IFROMC(CO(N3CO+1))
          CALL ASSERT(nb_new.GT.0.AND.nb_new.LE.NBT,'Basis not defined',
     '      ERROR,*9999)
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
          CALL ASSERT(NELIST(0).GT.0,
     '      '>> No Elements to refine',ERROR,*9999)
          IF(CBBREV(CO,'NODES_FROM',2,noco+1,NTCO,N3CO)) THEN
            NPTNEW=IFROMC(CO(N3CO+1))-1
          ELSE
            NPTNEW=NPT(0)
          ENDIF
        ELSE
          IF(CBBREV(CO,'NOCHANGE',3,noco+1,NTCO,N3CO)) THEN
            CHANGE=.FALSE.
          ELSE
            CHANGE=.TRUE.
          ENDIF
        ENDIF
        IF(BASIS.OR.UVERSIONS) THEN
          IF(CBBREV(CO,'GEOMETRY',2,noco+1,NTCO,N3CO)) THEN
            CALL ASSERT(CALL_ELEM,'Geometry mesh not defined',
     '        ERROR,*9999)
            NJL_TYPE=NJL_GEOM
            njj=IFROMC(CO(N3CO+1))
          ELSE IF(CBBREV(CO,'FIBRE',2,noco+1,NTCO,N3CO)) THEN
            CALL ASSERT(CALL_ELFB,'Fibre mesh not defined',ERROR,*9999)
            NJL_TYPE=NJL_FIBR
            njj=IFROMC(CO(N3CO+1))
          ELSE IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN
            CALL ASSERT(CALL_ELFD,'Field mesh not defined',ERROR,*9999)
            NJL_TYPE=NJL_FIEL
            njj=IFROMC(CO(N3CO+1))
          ELSE
            CALL ASSERT(CALL_ELFD,'Field mesh not defined',ERROR,*9999)
            NJL_TYPE=NJL_FIEL
            njj=1
          ENDIF
        ENDIF

        IF(BASIS) THEN
          LLOOSE_TOL=DSQRT(LOOSE_TOL)
          NITB=NIT(nb_new)
          NNTB=NNT(nb_new)
C         Find total number of nodes in each direction
          DO ni=1,NITB
            IF(IBT(1,ni,nb_new).EQ.1) THEN !Lagrange
              INPT(ni)=IBT(2,ni,nb_new)
            ELSE IF(IBT(1,ni,nb_new).EQ.2) THEN !Hermite
              INPT(ni)=1
            ENDIF
          ENDDO !ni
          CALL ASSERT(NBI(nb_new).LT.5.OR.NBI(nb_new).GT.7,
     '      'Use unit scale factors',ERROR,*9999)

          DO no_ne=1,NELIST(0)   !Loop over elements
            ne=NELIST(no_ne)
            nr=NRE(ne) !is region# for current element
            nj=NJ_LOC(NJL_TYPE,njj,nr)
            nb_old=NBJ(nj,ne)
            CALL ASSERT(NIT(nb_old).EQ.NITB,
     '        'Bases have differing numbers of xi-directions',
     '        ERROR,*9999)

            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '        nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)

C KAT 24Jul99: Simple algorithm that looks for each node it needs and
C              creates new ones if necessary.
            DO nn=1,NNTB
C             Find xi-coordinates of the node
              DO ni=1,NITB
                INPI=INP(nn,ni,nb_new)-1
                XI(ni)=DBLE(INPI)/INPT(ni)
              ENDDO !ni
C             Look for a node at this point
              CALL REFINE_FINDNODE(IBT,IDO,INP,nb_old,NBJ,ne,NKJ,np,
     '          NPLIST,NPNODE,NPTNEW,NPTOLD,nr,NUNK,NVJP,XE,XI,XP,
     '          EXISTS,ERROR,*9999)

C KAT 8Oct99: Initiallizing in REFINE_FINDNODE instead
C              IF(.NOT.EXISTS) THEN
CC               Initialize for all nj
C                DO njj1=1,3
C                  DO njj2=1,NJ_LOC(njj1,0,nr)
C                    nj1=NJ_LOC(njj1,njj2,nr)
C                    IF(nj1.LE.NJT) THEN
C                      NKJ(nj1,np)=1
C                      NUNK(1,nj1,np)=1
CC                     NVJP and XP(1,1,nj1,np) already set in REFINE_FINDNODE.
C                    ELSE
C                      NKJ(nj1,np)=0
C                      NVJP(nj1,np)=0
C                    ENDIF
C                  ENDDO !njj2
C                ENDDO !njj1
C              ENDIF !not EXISTS

C             Number of derivatives needed at the node
              NKJNEED=NKT(nn,nb_new)
C             Check how many derivatives already at the node are used
              nk=NKJ(nj,np)
              USED=.FALSE. !for now
              DO WHILE(nk.GT.0.AND..NOT.USED)
                IF(NUNK(nk,nj,np).NE.0) THEN
                  USED=.TRUE.
                ELSE
                  nk=nk-1
                ENDIF
              ENDDO
C             Number of derivatives that can be compared with existing node
              NKJCHEK=MIN(nk,NKJNEED)
C             One more than number of derivs already at node
              NKSTART=nk+1
              IF(NKJNEED.GE.NKSTART) THEN
                CALL ASSERT(NKJNEED.LE.NKM,'Increase NKM',ERROR,*9999)
                NKJ(nj,np)=NKJNEED
              ENDIF !NKJNEED.GE.NKSTART
              DO nk=1,NKJNEED
                nu=IDO(nk,nn,0,nb_new)
                NUNK(nk,nj,np)=nu
                X(nu)=PXI(IBT(1,1,nb_old),IDO(1,1,0,nb_old),
     '            INP(1,1,nb_old),nb_old,nu,XI,XE(1,nj))
C               If an angle is 2*pi set it to zero
                IF((nj.EQ.2.AND.ITYP10(nr).GE.2).OR.
     '            (nj.EQ.3.AND.ITYP10(nr).GE.3)) THEN
                  IF(DABS(X(1)-2.0d0*PI).LT.LOOSE_TOL) X(1)=0.0d0
                ENDIF
              ENDDO !nk
C             Determine the version to use by checking if a version with the
C             appropriate information already exists.
              NEWVERSION=.TRUE. !assume a new version unless we find a match
              nv=0
              DO WHILE(nv.LT.NVJP(nj,np).AND.NEWVERSION)
                nv=nv+1
                NEWVERSION=.FALSE. !assume node matches until it doesn't
                nk=0
                DO WHILE(nk.LT.NKJCHEK.AND..NOT.NEWVERSION)
                  nk=nk+1
                  nu=NUNK(nk,nj,np)
                  DIFF=DABS(X(nu)-XP(nk,nv,nj,np))
                  IF(DIFF.GT.LLOOSE_TOL*(DABS(X(nu))+LLOOSE_TOL))
     '              NEWVERSION=.TRUE. !doesn't match
                ENDDO !nk
              ENDDO !nv
              IF(NEWVERSION) THEN
                nv=nv+1
                NKSTART=1
                CALL ASSERT(nv.LE.NVM,' >>ERROR: Increase NVM',ERROR,
     '            *9999)
                NVJP(nj,np)=nv
              ENDIF !NEWVERSION
C             Include new derivatives
              DO nk=NKSTART,NKJNEED
                nu=NUNK(nk,nj,np)
                XP(nk,nv,nj,np)=X(nu)
              ENDDO !nk

              NPNE(nn,nb_new,ne)=np
              NVJE(nn,nb_new,nj,ne)=nv
              DO nk=1,NKJNEED
                NKJE(nk,nn,nj,ne)=nk
              ENDDO

            ENDDO !nn

            DO ns=1,NST(nb_new)
              SE(ns,nb_new,ne)=0d0
            ENDDO !se

            NBJ(nj,ne)=nb_new

          ENDDO !elements

! Update interface info
          CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)

          DO nr=1,NRT
            CALL ASSERT(NPNODE(0,nr).LE.NP_R_M,'>>Increase NP_R_M',
     '        ERROR,*9999)
          ENDDO !nr
          CALL ASSERT(NPNODE(0,0).LE.NPM,'>>Increase NPM',
     '      ERROR,*9999)

c cpb 20/7/95 NENP needs to be calculated before the lines are
c calculated in LINSEG
          CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

          IF(NBI(nb_new).EQ.2.OR.NBI(nb_new).EQ.3) THEN
! scale factors based on specified element or global derivs
C ***           Transfer SE to DL
            CALL LINSEG(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '        NLLINE,NNL,NPL,NPNE,NPNODE,NVJE,NVJL,XP,ERROR,*9999)
            CALL SEDL(IBT,IDO,nb_new,NEELEM,NLL,NNL,NPL,DL,SE,
     '        ERROR,*9999)
          ENDIF

          CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '      NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '      DL,SE,XP,ERROR,*9999)
          CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '      NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,
     '      NRE,NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,
     '      XE,XG,XP,ERROR,*9999)

          CALL_EQUA=.FALSE.
          CALL_FIT=.FALSE.
          CALL_OPTI=.FALSE.
          CALL_INIT=.FALSE. !Initial conditions need to be set at new
          CALL_SOLV=.FALSE.  !nodes

        ELSE IF(UVERSIONS) THEN !unique versions
C KAT 18Jul99:  At present this is a simple algorithm that copies
C               version 1 to all other versions.
C         Set all nodes to zero versions
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            nj=NJ_LOC(NJL_TYPE,njj,nr)
            DO no_np=1,NPNODE(0,nr)
              np=NPNODE(no_np,nr)
              NVJP(nj,np)=0
            ENDDO !np
          ENDDO !nr
          INCREASE_NVM=0
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            nj=NJ_LOC(NJL_TYPE,njj,nr)
            DO no_ne=1,NEELEM(0,nr)
              ne=NEELEM(no_ne,nr)
              nb=NBJ(nj,ne)
              DO nn=1,NNT(nb)
                np=NPNE(nn,nb,ne)
                nv=NVJP(nj,np)+1
                NVJP(nj,np)=nv
                IF(nv.LE.NVM) THEN
                  NVJE(nn,nb,nj,ne)=nv
                  IF(nv.GT.1) THEN
                    DO nk=1,NKT(0,nb)
                      XP(nk,nv,nj,np)=XP(nk,1,nj,np)
                    ENDDO !nk
                  ENDIF !nv>1
                ELSE IF(nv.GT.INCREASE_NVM) THEN
                  INCREASE_NVM=nv
                ENDIF !nv<=NVM
              ENDDO !nn
            ENDDO !ne
          ENDDO !nr

          IF(INCREASE_NVM.NE.0) THEN
            IEND=0
            CALL APPENDC(IEND,'Increase NVM to ',ERROR)
            CALL APPENDI(IEND,INCREASE_NVM,ERROR)
            GOTO 9999
          ENDIF !INCREASE_NVM

C!!!      need to update NVJL

        ELSE

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            IF(CHANGE) THEN !Assume mesh has been changed
              IF(NJT.EQ.2) THEN
                CALL CHMESH_CALC_ARCLENGTHS(IBT,IDO,INP,NBJ,NEELEM,NGAP,
     '            NPNE,nr,DL,SE,TOTAL,WG,XIG,XP,0.0d0,ERROR,*9999)
              ENDIF
              CALL CHMESH_ALTER(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,
     '          NP_INTERFACE,NPL,NPNE,NPNODE,nr,NVJL,DL,SE,TOTAL,XP,
     '          ERROR,*9999)
C   Mesh altered.  Update scale factors.
              CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '          NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '          DL,SE,XP,ERROR,*9999)
            ENDIF

C Calculate a circumferential measurement (if requested)
            IF(CBBREV(CO,'MEASUREMENT',2,noco+1,NTCO,N3CO)) THEN
C 2-MAY-1998 LKC need to pass in an array
C            CALL PARSRL(CO(N3CO+1),1,NTRL,Z,ERROR,*9999)
              CALL PARSRL(CO(N3CO+1),1,NTRL,REAL_TEMP,ERROR,*9999)
              Z=REAL_TEMP(1)
              CALL CHMESH_CALC_ARCLENGTHS(IBT,IDO,INP,NBJ,NEELEM,NGAP,
     '          NPNE,nr,DL,SE,TOTAL,WG,XIG,XP,Z,ERROR,*9999)
C 14-OCT-2002 LKC need to supply a format for this output due
C    to diffent defaults for different OS's
C              WRITE(OP_STRING,*)' MEASUREMENT =',TOTAL
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

              WRITE(OP_STRING,'('' MEASUREMENT = '',E12.5)') TOTAL
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

            ENDIF
            IF(CBBREV(CO,'TWODMEAS',2,noco+1,NTCO,N3CO)) THEN
C 2-MAY-1998 LKC need to pass in an array
C            CALL PARSRL(CO(N3CO+1),1,NTRL,Z,ERROR,*9999)
              CALL PARSRL(CO(N3CO+1),1,NTRL,REAL_TEMP,ERROR,*9999)
              Z=REAL_TEMP(1)
              CALL MESHXY(IBT,IDO,INP,NBJ,nr,NEELEM,NPNE,DEPTH,SE,
     '          XP,WIDTH,Z,ERROR,*9999)
              WRITE(OP_STRING,*)'x = ',WIDTH,' y= ',DEPTH
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO

        ENDIF !UVERSIONS
      ENDIF


      CALL EXITS('CHMESH')
      RETURN
 9999 CALL ERRORS('CHMESH',ERROR)
      CALL EXITS('CHMESH')
      RETURN 1
      END


