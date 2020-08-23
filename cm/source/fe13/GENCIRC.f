      SUBROUTINE GENCIRC(MAX_STRAHLER,nb,NBJ,NEELEM,NPNE,NPNODE,
     '   nr,NRE,NXI,ORDER,ANG,CE,nhx,XP,ERROR,*)


C#### Subroutine: GENCIRC
C###  Description:
C###    GENCIRC generates pulmonary venous and arterial trees.
C###    Created by KSB and SNH May 2001

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'mesh00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter List
      INTEGER MAX_STRAHLER,nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nhx,nr,NRE(NEM),nr0,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),ORDER
      REAL*8 ANG,angle,CE(NMM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER factor(2),ne,ne0,nj,no_branch,noelem,nonode,np,np0,ord
      REAL*8 DIAMETER(11),LENGTH(11)

      CALL ENTERS('GENCIRC',*9999)
      IF (CIRC_MODEL(nhx).EQ.1) THEN !Strahler ordering
        DO nr0=1,MAX_STRAHLER,1
          IF (nhx.EQ.1) THEN     !Arterial tree
            DIAMETER(nr0)=HUANG_ARTERY_DIAM(nr0)
            LENGTH(nr0)=HUANG_ARTERY_LENGTH(nr0)
          ELSEIF (nhx.EQ.2) THEN  !Venous tree
            DIAMETER(nr0)=HUANG_VEIN_DIAM(nr0)
            LENGTH(nr0)=HUANG_VEIN_LENGTH(nr0)
            N_STRAHLER_ORD(nr0)=N_STRAHLER_ORDVEN(nr0)
          ENDIF !nhx.EQ.1
        ENDDO !nr0
        np=NPT(0)+1                      !Initial node
        ne = NET(0)+1
        DO nj=1,4
          XP(1,1,nj,np)=0.d0      !Starting position.
          XP(2,1,nj,np)=0.d0
        ENDDO !nj
        noelem=1
        nonode=1
        NPNE(1,nb,ne)=np            !First element is between nodes
        NPNE(2,nb,ne)=np+1            !1 and 2.
        NPNODE(1,nr)=np
        NPNODE(2,nr)=np+1
        NEELEM(noelem,nr)=ne
        NBJ(1,ne)=nb
        NBJ(2,ne)=nb
        NBJ(3,ne)=nb
        ord=MAX_STRAHLER       !Start off at the largest order.
        CE(6,NE)=DBLE(ord)
        CE(9,ne)=LENGTH(11)   !Initial branch length
        CE(10,ne)=CE(9,ne)     !Current branch length
        CE(11,ne)=CE(10,ne)    !Branch length at previous time-step
        CE(12,ne)=(PI*DIAMETER(11)**2.d0)/4.0d0 !Current branch cross-sectional area
        CE(14,ne)=CE(12,ne)*CE(10,ne) !Current branch volume
        CE(13,ne)=CE(14,ne)    !Initial branch volume
        NXI(-1,0,ne)=0         !No parent branch until coupled to heart
        np=np+1
        DO nj=1,4              !Derivatives not defined
          XP(1,1,nj,(np-1))=0.d0
          XP(1,1,nj,np)=0.d0
        ENDDO !nj
        XP(1,1,3,np)=XP(1,1,3,(np-1))+CE(10,ne)    !Move along z axis.
        N_STRAHLER_ORD(ord)=N_STRAHLER_ORD(ord)-1
        ne0=NEELEM(1,nr)-1  !Initialising parent element
        factor(1)=1          !Used to branch left (N=1)
        factor(2)=-1         !or right (N=2)
        DO WHILE(ord.GE.ORDER) !While there are still elements to add.
          ne0=ne0+1 !Parent element
          np0=NPNE(2,nb,ne0) !End node of parent.
          DO no_branch=1,2   !Set up 2 daughters, find details
            noelem=noelem+1
            nonode=nonode+1
            NPNODE(nonode,nr)=np
            ne=ne+1          !New element
            np=np+1          !New node
            NEELEM(noelem,nr)=ne
            NPNE(1,nb,ne)=np0
            NPNE(2,nb,ne)=np
            NBJ(1,ne)=nb
            NBJ(2,ne)=nb
            NBJ(3,ne)=nb
            NXI(-1,0,ne)=1   !One parent
            NXI(-1,1,ne)=ne0
            NXI(1,0,ne0)=no_branch   !Parent has two daughters
            NXI(1,no_branch,ne0)=ne
            DO nj=1,4  !Derivatives not defined
              XP(1,1,nj,np)=0.d0
            ENDDO !nj
            IF ((no_branch.EQ.2)
     '        .AND.(CE(6,ne-1).EQ.CE(6,ne0)))THEN
              ord=ord-1  !If the other daughter = the parent
            ELSEIF (N_STRAHLER_ORD(ord).EQ.0) THEN
              ord=ord-1
            ELSEIF(ord.ne.MAX_STRAHLER) THEN
              IF((NO_BRANCH.EQ.1)  !If there are still higher orders to be used
     '          .and.(N_STRAHLER_ORD(ord+1).NE.0)) ord=ord+1
            ENDIF
            CE(6,ne)=DBLE(ord) !Strahler Order
            N_STRAHLER_ORD(ord)=N_STRAHLER_ORD(ord)-1
            CE(9,ne)=LENGTH(ord) !Branch lengths
            CE(10,ne)=CE(9,ne)
            CE(11,ne)=CE(10,ne) !XP(*,*,4,*) stores angle w.r.t. z axis.
            angle=ANG*(PI/180)  !converting degrees to radians
            XP(1,1,4,np)=XP(1,1,4,np0)+factor(no_branch)*angle
            XP(1,1,3,np)=XP(1,1,3,np0)+CE(10,ne)*COS(XP(1,1,4,np))
            XP(1,1,2,np)=XP(1,1,2,np0)+CE(10,ne)*SIN(XP(1,1,4,np))
            CE(12,ne)=(PI*DIAMETER(ord)**2.d0)/4.0d0 !x-sec
            CE(14,ne)=CE(12,ne)*CE(10,ne) !Branch volumes
            CE(13,ne)=CE(14,ne)
            NRE(ne)=nr
          ENDDO !NO_BRANCH
        ENDDO !ord.ne.0
        nonode=nonode+1
        NPNODE(nonode,nr)=np
        NET(nr)=ne      !highest element # in region nr
        NET(0)=ne       !highest element # in any region
        NEELEM(0,nr)=noelem !total # of elements in reion nr
        NEELEM(0,0)=ne  !total # of elements in all regions
        NPNODE(0,nr)=nonode !total # of nodes in region nr
        NPNODE(0,0)=np  !total # of nodes in all regions
        NPT(nr)=np      !highest node in region nr
        NPT(0)=np       !highest node in all regions
      ENDIF ! CIRC_MODEL(nhx).EQ.1
C OPCIRC has now been removed May 2003. May also remove GENCIRC
c do not currently use - will maybe use as a test mesh for flow solution
C      CALL OPCIRC(MAX_STRAHLER,NEELEM,nr,N_STRAHLER_ORD,NXI,ord,CE,
C     '  DIAMETER,ERROR,*9999) !Output stats on mesh
      
      CALL EXITS('GENCIRC')
      RETURN
 9999 CALL ERRORS('GENCIRC',ERROR)
      CALL EXITS('GENCIRC')
      RETURN 1
      END


