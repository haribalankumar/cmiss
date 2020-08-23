      SUBROUTINE MVINTNOD(NFVC,NODENVC,NPNODE,nr,NVCNODE,XNFV,XP,
     '  ERROR,*)

C#### Subroutine: MVINTNOD
C###  Description:
C###    It performs one sweep through the nodes, and is called from the
C###    outside numerous times to smooth the mesh.
C###    The internal nodes of a Voronoi mesh are moved
C###    according to a face weighted smoothing process. Care needs to
C###    be taken here that there is enough a small enough time step
C###    between moving the boundary & internal boundary nodes at each
C###    iteration, otherwise volume distribution may lag boundary motion


      IMPLICIT NONE
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVCNODE(2,NP_R_M)
      REAL*8 XNFV(-(NJM+1):NJM,NFVM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,nvc,nfvl,np,nj,cnp,nfv,nonode,cnonode,njj
      REAL*8 OLDPOS(3),POSADD(3),CELLAREA

      CALL ENTERS('MVINTNOD',*9999)

      DO I=1,SMOOTHIT
        DO nvc=1,NVCT
          nonode=NODENVC(nvc)

C         ..Boundary and internal boundary nodes have been moved,
C         we only wish to move the internal nodes.
          IF(NVCNODE(TYPE,nonode).EQ.BOUNDARY) THEN
            np=NPNODE(nonode,nr)

C           ..Get the original position, and zero the face
C           weighted  addition to position
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              OLDPOS(njj)=XP(1,1,nj,np)
              POSADD(njj)=0.d0
            ENDDO

            CELLAREA=0.d0

C           ..Compute the face weighted addition to position, and
C           the total surface area of the cell
            DO nfvl=1,NFVC(1,0,nvc)
              cnonode=NFVC(1,nfvl,nvc)
              cnp=NPNODE(cnonode,nr)
              nfv=NFVC(2,nfvl,nvc)
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj,nr)
                POSADD(njj)=POSADD(njj)+XP(1,1,nj,cnp)*XNFV(FAREA,
     '            nfv)
              ENDDO
              CELLAREA=CELLAREA+XNFV(FAREA,nfv)
            ENDDO

C           ..Store new node positions, and calculate the velocity
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)

C             ..new position
              XP(1,1,nj,np)=POSADD(njj)/CELLAREA

C             ..node velocity
              XP(1,2,nj,np)=XP(1,2,nj,np)+(XP(1,1,nj,np)
     '          -OLDPOS(njj))/DT
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      CALL EXITS('MVINTNOD')
      RETURN
 9999 CALL ERRORS('MVINTNOD',ERROR)
      CALL EXITS('MVINTNOD')
      RETURN 1
      END


