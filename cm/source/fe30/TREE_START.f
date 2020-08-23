      SUBROUTINE TREE_START(GEOM,IBT,IDO,INP,NBJ,NEELEM,
     '  NKJE,NPF,NP_INTERFACE,NPNE,nr2,NTREES,NVJE,SE,
     '  START,STARTNODE,SURF_TREE,TMAX,XA,XE,XP,ERROR,*)

C#### Subroutine: TREE_START
C###  Description:
C###    TREE_START defines the starting positions for the trees
C###    grown by IPMESH8.
C###  Created by Martin Buist December 1996

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'

!     Parameter list
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  nr2,nr3,ntt,NTREES,
     '  NVJE(NNM,NBFM,NJM,NEM),TMAX
      REAL*8 SE(NSM,NBFM,NEM),START(TMAX,3),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL SURF_TREE
!     Local variables
      INTEGER IBEG,ICHAR,IEND,INFO,LD,ne,NITB,ni,
     '  nj,noelem,NOQUES,STARTNODE(TMAX,2)
      REAL*8 XI(3),XP1(3),Z(3)
      CHARACTER CHAR1*10
      LOGICAL FILEIP,FOUND,GEOM

      CALL ENTERS('TREE_START',*9999)

C Initialisation
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      ne=0
      nr2=1
      NTREES=1
      DO ntt=1,8
        STARTNODE(ntt,1)=0
        STARTNODE(ntt,2)=0
        DO nj=1,NJT
          START(ntt,nj)=0.0d0
        ENDDO !nj
      ENDDO !ntt
      DO nj=1,NJT
        XP1(nj)=0.0d0
        XI(nj)=0.0d0
      ENDDO !nj

      SURF_TREE=.TRUE.
      IF(IOTYPE.EQ.3) THEN
        IF(SURF_TREE) THEN
          ADATA(1)='Y'
        ELSE
          ADATA(1)='N'
        ENDIF
      ENDIF
      FORMAT='($,'' Grow on a 2d surface in 3d space [Y]? '',A)'
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        IF(ADATA(1).EQ.'Y') THEN
          SURF_TREE=.TRUE.
        ELSE
          SURF_TREE=.FALSE.
        ENDIF
      ENDIF


      IDEFLT(1)=1
      FORMAT='($,'' How many trees are there to grow [1] '',I2)'
      IF(IOTYPE.EQ.3) IDATA(1)=NTREES
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NTREES=IDATA(1)
      CALL ASSERT(NTREES.LE.8,'>>ERROR maximum number of trees is 8',
     '  ERROR,*9999)
      CALL ASSERT(NTREES.GE.1,'>>ERROR minimum number of trees is 1',
     '  ERROR,*9999)

      IDEFLT(1)=1
      FORMAT='($,'' Enter the region to grow trees within [1]: '',I2)'
      IF(IOTYPE.EQ.3) IDATA(1)=nr2
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NRM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) nr2=IDATA(1)

      DO ntt=1,NTREES
        WRITE(CHAR1,'(I1)') ntt
        CALL STRING_TRIM(CHAR1,IBEG,IEND)
        GEOM=.TRUE.
        IF(IOTYPE.EQ.3) THEN
          IF(GEOM) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        FORMAT='($,'' Do you want to grow from an existing '
     '    // 'element for tree '//CHAR1(IBEG:IEND)//' [Y]? '',A)'
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            GEOM=.TRUE.
          ELSE
            GEOM=.FALSE.
          ENDIF
        ENDIF

        IF(GEOM) THEN
          IDEFLT(1)=1
          FORMAT='($,'' Enter node number to grow from [1] '',I4)'
          IF(IOTYPE.EQ.3) IDATA(1)=STARTNODE(ntt,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) STARTNODE(ntt,1)=IDATA(1)
          CALL ASSERT(NP_INTERFACE(STARTNODE(ntt,1),0).GT.0,
     '      '>>ERROR no region found for the node',ERROR,*9999)
          nr3=NP_INTERFACE(STARTNODE(ntt,1),1)
          DO nj=1,NJT
            XP1(nj)=XP(1,1,nj,STARTNODE(ntt,1))
          ENDDO !nj
          IF(ITYP10(nr3).NE.1) THEN !need rc for dexi_point
            CALL COORD(ITYP10(nr3),1,XP1,Z,ERROR,*9999)
            DO nj=1,NJT
              XP1(nj)=Z(nj)
            ENDDO
          ENDIF

          LD=0
          FOUND=.FALSE.
          noelem=1
          NITB=NIT(NBJ(1,NEELEM(1,nr2)))
          DO WHILE((.NOT.FOUND).AND.(noelem.LE.NEELEM(0,nr2)))
            ne=NEELEM(noelem,nr2)
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '        nr2,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,
     '        ERROR,*9999)
            DO ni=1,NITB !initialising XI
              XI(ni)=0.5d0
            ENDDO
            CALL DEXI_POINT(IBT,IDO,INP,LD,NBJ,ne,NITB,nr2,0.d0,XE,XI,
     '        XI,XP1,.FALSE.,ERROR,*9999)
            IF(LD.NE.0) THEN
              FOUND=.TRUE.
            ELSE
              noelem=noelem+1
            ENDIF
          ENDDO

          CALL ASSERT(LD.NE.0,'>>ERROR element not found for node',
     '      ERROR,*9999)

C          CALL XCOORD(IBT,IDO,INP,NBJ,ne,NEELEM,NKE,NPF,
C     '      NPNE,NPNODE,NVJE,SE,LOOSE_TOL,XA,XE,XI,XP,XP1,ERROR,
C     '      *9999)
C          CALL ASSERT(ne.GT.0,
C     '      '>>ERROR start must lie within an element',ERROR,*9999)

          STARTNODE(ntt,2)=LD
          DO nj=1,NJT
            START(ntt,nj)=XI(nj)
          ENDDO !nj

        ELSE  ! not geom

          DO nj=1,NJT
            RDEFLT(1)=0.0d0
            WRITE(CHAR1,'(I1)') nj
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            FORMAT='($,'' Enter starting coordinate '
     '        //CHAR1(IBEG:IEND)// ' [0]:'',F8.6)'
            IF(IOTYPE.EQ.3) RDATA(1)=START(ntt,nj)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) START(ntt,nj)=RDATA(1)
          ENDDO !nj

        ENDIF

        IF(DOP) THEN
          WRITE(OP_STRING,'('' For tree number '',I4)') ntt
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Starting node number is '',I4)')
     '      STARTNODE(ntt,1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Starting element number is '',I4)')
     '      STARTNODE(ntt,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Xi coords '',F6.4,F6.4,F6.4)')
     '      START(ntt,1),START(ntt,2),START(ntt,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO !ntt

      CALL EXITS('TREE_START')
      RETURN
 9999 CALL ERRORS('TREE_START',ERROR)
      CALL EXITS('TREE_START')
      RETURN 1
      END



