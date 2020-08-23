      SUBROUTINE CALC_VORO(NBJ,NENFVC,NENP,NFVC,NODENVC,NPLIST,
     '  NPNE,NPNODE,nr,NXI,VC,VC_INIT,XNFV,XP,ZA,ERROR,*)

C#### Subroutine: CALC_VORO
C###  Description:
C###    CALC_VORO calculates the voronoi mesh

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'mxch.inc'
c      INCLUDE 'cmiss$reference:voro00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NENFVC(0:NFVCM,NFVM),NENP(NPM,0:NEPM,0:NRM),
     '  NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 VC(0:NVCM),VC_INIT(2,NVCM),XNFV(-(NJM+1):NJM,NFVM),
     '  XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nvc,nonode,np,noelem,ne,nb,nn,ADJNODE,cnonode,nfvl,nfv,
     '  cnp
      LOGICAL FOUND,ADD1,ADD2

      CALL ENTERS('CALC_VORO',*9999)
C     CALL ASSERT(CALL_RECONNECT,'>>Need to reconnect Delaunay'//
C     '  ' mesh before calculating Voronoi diagram',
C     '  ERROR,*9999)
C     ..Get the nonode numbers for each np
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        NPLIST(np)=nonode
      ENDDO
      
C     ..Initialisation
      nfv=1
      
C     ..Internal nodes only
      DO nvc=1,NVCT
        
C       Initialisation
        NFVC(1,0,nvc)=0
        NFVC(2,0,nvc)=0
        VC(nvc)=0.d0
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)

C       ..Loop over the elements surrounding this node to get the
C         nodes adjoining node np
        DO noelem=1,NENP(np,0,nr)
          ne=NENP(np,noelem,nr)

C         ..Loop over the element's local nodes to get the surrounding
C           nodes
          nb=NBJ(1,ne)
          DO nn=1,NNT(nb)
            ADJNODE=NPNE(nn,nb,ne)

            IF(ADJNODE.NE.np) THEN

C             ..Find the corresponding nonode of ADJNODE
              cnonode=NPLIST(ADJNODE)

C             ..Make sure that we haven't done this node already
              ADD1=.TRUE.
              FOUND=.FALSE.
              nfvl=1
              DO WHILE(.NOT.FOUND.AND.nfvl.LE.NFVC(1,0,nvc))
                cnp=NPNODE(NFVC(1,nfvl,nvc),nr)
                IF(ADJNODE.EQ.cnp) THEN
                  ADD1=.FALSE.
                  FOUND=.TRUE.
                ELSE
                  nfvl=nfvl+1
                ENDIF
              ENDDO
            ELSE
              ADD1=.FALSE.
            ENDIF

            IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2.AND.ADD1) THEN

C               ..Calculate face information & volume contribution
                CALL CALC_VORO_FACE2D(ADJNODE,NBJ,NENP,np,NPNE,nr,
     '          VC(nvc),XNFV(-(NJM+1),nfv),XP,ZA,ADD2,ERROR,*9999)

            ELSEIF(NJ_LOC(NJL_GEOM,0,nr).EQ.3.AND.ADD1) THEN

C               ..Calculate face information & volume contribution
              CALL CALC_VORO_FACE3D(ADJNODE,NBJ,NENFVC,NENP,nfv,np,NPNE,
     '          nr,NXI,VC(nvc),XNFV(-(NJM+1),nfv),XP,ZA,ADD2,ERROR,
     '          *9999)

            ENDIF
            
            IF(ADD1.AND.ADD2) THEN
              NFVC(1,0,nvc)=NFVC(1,0,nvc)+1
              NFVC(2,0,nvc)=NFVC(2,0,nvc)+1
              CALL ASSERT(NFVC(1,0,nvc).LE.NFVCM,
     '          '>>Increase NFVCM',ERROR,*9999)
              NFVC(1,NFVC(1,0,nvc),nvc)=cnonode
              NFVC(2,NFVC(2,0,nvc),nvc)=nfv
              nfv=nfv+1

C             ..Increase NFVCM 'cos NFVM not in ippara file, its
C               calculated from NVCM*(NFVCM/2)
              CALL ASSERT(nfv.LE.NFVM,'>>Increase NFVCM',ERROR,*9999)
            ENDIF
          ENDDO !nn
        ENDDO !noelem
        VC(nvc)=DABS(VC(nvc))
      ENDDO !nvc
C     ..Calculate total number of faces
      NFVT=nfv-1

C     ..Calculate the total volume of the mesh
      VC(0)=0.d0
      DO nvc=1,NVCT
        VC(0)=VC(0)+VC(nvc)
        VC_INIT(1,nvc)=VC(nvc)
        VC_INIT(2,nvc)=VC(nvc)
      ENDDO

C     ..Ensure that reconnect is called before this routine is
C      CALL_RECONNECT=.FALSE.

      CALL EXITS('CALC_VORO')
      RETURN
 9999 CALL ERRORS('CALC_VORO',ERROR)
      CALL EXITS('CALC_VORO')
      RETURN 1
      END



