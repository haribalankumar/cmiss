      SUBROUTINE VORONOI(IBT,nb_voronoi,N_BDRY,N_IBDRY,N_INTNL,NBJ,
     '  NEELEM,NENFVC,NENP,NFVC,NKJ,NKJE,NMAX_GENERATION,NODENVC,
     '  NODENVCB,NP_INTERFACE,NPLIST,NPNE,NPNODE,nr,NRE,nr_host,
     '  nr_target,NVCNODE,NVJE,NVJP,NXI,SE,VC,VC_INIT,XNFV,XP,ZA,
     '  CONVERT_MESH,ERROR,*)

C#### Subroutine: VORONOI
C###  Description:
C###    Voronoi routine.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),nb_voronoi,N_BDRY,N_IBDRY,N_INTNL,
     &  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NENFVC(0:NFVCM,NFVM),
     &  NENP(NPM,0:NEPM,0:NRM),NFVC(2,0:NFVCM,NVCM),NKJ(NJM,NPM),
     &  NKJE(NKM,NNM,NJM,NEM),NMAX_GENERATION,NODENVC(NVCM),
     &  NODENVCB(NVCBM),NP_INTERFACE(0:NPM,0:3),NPLIST(0:NPM),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),nr_host,
     &  nr_target,NVCNODE(2,NP_R_M),NVJE(NNM,NBFM,NJM,NEM),
     &  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),VC(0:NVCM),VC_INIT(2,NVCM),
     '  XNFV(-(NJM+1):NJM,NFVM),XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM)
      LOGICAL CONVERT_MESH
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER bnvc,n,nonode,novc,np,np_to_find,nvc

      CALL ENTERS('VORONOI',*9999)

C***  Voronoi cells are created from Delaunay triangulation
C     ..Initialisation
      DO nonode=1,NPNODE(0,nr)
        NVCNODE(TYPE,nonode)=0
      ENDDO !nonode
      VC(0)=0.d0
      DO novc=1,NVCM
        VC(novc)=0.d0
      ENDDO !novc

      DO n=1,N_BDRY
        np_to_find=NPLIST(n)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          IF(np.EQ.np_to_find) THEN !found it
            CALL ASSERT(NVCNODE(TYPE,nonode).EQ.0,
     '        '>>Node already specified',ERROR,*9999)
            NVCNODE(TYPE,nonode)=BOUNDARY
            GOTO 30
          ENDIF
        ENDDO !nonode
 30     CONTINUE
      ENDDO !n

      DO n=1,N_IBDRY
        np_to_find=NPLIST(N_BDRY+n)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          IF(np.EQ.np_to_find) THEN !found it
            CALL ASSERT(NVCNODE(TYPE,nonode).EQ.0,
     '        '>>Node already specified',ERROR,*9999)
            NVCNODE(TYPE,nonode)=INTBOUN
            GOTO 60
          ENDIF
        ENDDO !nonode
 60     CONTINUE
      ENDDO !n

      DO n=1,N_INTNL
        np_to_find=NPLIST(N_BDRY+N_IBDRY+n)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          IF(np.EQ.np_to_find) THEN !found it
            CALL ASSERT(NVCNODE(TYPE,nonode).EQ.0,
     '        '>>Node already specified',ERROR,*9999)
            NVCNODE(TYPE,nonode)=INTERNAL
            GOTO 90
          ENDIF
        ENDDO !nonode
 90     CONTINUE
      ENDDO !n
      NVCT=N_IBDRY+N_INTNL

      DO nonode=1,NPNODE(0,nr)
        CALL ASSERT(NVCNODE(TYPE,nonode).NE.0,
     '    '>> all nodes need to be specified a type',ERROR,*9999)
      ENDDO

      nvc=1
      bnvc=1
      DO nonode=1,NPNODE(0,nr)
        IF(NVCNODE(TYPE,nonode).EQ.INTBOUN.OR.
     '    NVCNODE(TYPE,nonode).EQ.INTERNAL) THEN
          NVCNODE(MAP,nonode)=nvc
          CALL ASSERT(nvc.LE.NVCM,'>> Increase NVCM',ERROR,*9999)
          NODENVC(nvc)=nonode
          nvc=nvc+1
        ELSEIF(NVCNODE(TYPE,nonode).EQ.BOUNDARY) THEN
          NVCNODE(MAP,nonode)=bnvc
          NODENVCB(bnvc)=nonode
          bnvc=bnvc+1
          CALL ASSERT(bnvc.LE.NVCBM,'>> Increase NVCBM',ERROR,*9999)
        ELSE
          ERROR='Incorrect type found in Voronoi node'
          GOTO 9999
        ENDIF
      ENDDO

C     ..Calculate the Voronoi mesh arrays
      CALL CALC_VORO(NBJ,NENFVC,NENP,NFVC,NODENVC,NPLIST,NPNE,NPNODE,
     '  nr,NXI,VC,VC_INIT,XNFV,XP,ZA,ERROR,*9999)

C SMAR009 19/01/99 removed NVCNODE,
      IF(CONVERT_MESH)THEN
C       IF(DEL_RET.EQ.2)THEN
C... check basis funciton is 2-D, linear (uses 4 nodes)
        CALL ASSERT(NIT(nb_voronoi).EQ.2,'>>basis function must be 2-D',
     '    ERROR,*9999)
        CALL ASSERT(IBT(1,1,nb_voronoi).EQ.1,
     '    '>>basis function must be linear',ERROR,*9999)
C       Code under development for generating alveolar mesh.
        CALL CALC_VORO_ALVEOLI(nb_voronoi,NBJ,NEELEM,NENFVC,NENP,NFVC,
     '    NKJ,NKJE,NMAX_GENERATION,NODENVC,NP_INTERFACE,NPNE,NPNODE,nr,
     '    NRE,nr_host,nr_target,NVCNODE,NVJE,NVJP,SE,XP,ZA,ERROR,*9999)
      ENDIF

C set up NFVCL
c      DO nvc=1,NVCT
c        nonode=NODENVC(nvc) !node corresponding to cell
c        np=NPNODE(nonode,nr)
c        nfvl=0
c        nfv=0
c        DO nonode2=1,NFVC(1,0,nvc) !for each connected node
c ! find common elements and set up mapping for faces
c          nonode3=NFVC(1,nonode2,nvc) !local connected node #
c          np2=NPNODE(nonode3,nr)
c          NELST(0)=0 !find list of common elements (->NELST)
c          DO noelem=1,NENP(np,0,2) !for all ne connected to np
c            ne=NENP(np,noelem,2) !global element #
c            DO noelem2=1,NENP(np2,0,2) !for all ne connected to np2
c              ne2=NENP(np2,noelem2,2) !global element #
c              IF(ne2.EQ.ne)THEN !is common to both nodes
c                NELST(0)=NELST(0)+1
c                NELST(NELST(0))=ne
c              ENDIF !ne2.eq.ne
c            ENDDO !noelem2
c          ENDDO !noelem
c          nfv=nfv+1 !increment the number of faces
c          NFVCL(0,nfv)=NELST(0) !# of Voronoi vertices in face nfv

c        ENDDO !nonode2
c      ENDDO !nonode

      CALL EXITS('VORONOI')
      RETURN
 9999 CALL ERRORS('VORONOI',ERROR)
      CALL EXITS('VORONOI')
      RETURN 1
      END


