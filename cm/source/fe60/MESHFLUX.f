      SUBROUTINE MESHFLUX(NFVC,NODENVC,NPNODE,nr,XNFV,XP,ZNFVMSH,
     '  ERROR,*)

C#### Subroutine: MESHFLUX
C###  Description:
C###    <HTML>
C###    The velocities of the nodes are used to estimate the scalar
C###    quantity meshflux, defined as the rate of volume swept by the
C###    face. Thus it will be the normal component of the face velocity,
C###    which can be obtained by a simple average, multiplied by the
C###    face area
C###    <PRE>
C###
C###    Meshflux(i,j)=0.5d0*Area(i,j)(NodeVel_(i) + NodeVel_(j)).N_(i,j)
C###
C###    .OR. centroid based mesh flux (More Accurate)
C###
C###    Meshflux(i,j) = Area(i,j)/d(i,j)*
C###                    [(X_(j)- c_(i,j)).NodeVel_(j) +
C###                     (c_(i,j) - X_(i)).NodeVel_(i)]
C###
C###    Where Area    = Area of Voronoi face
C###          NodeVel = Velocity of the node
C###          N       = unit outward normal from node i
C###          X       = node position
C###          c       = face centroid
C###          d       = internodal distance
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'voro00.inc'
      INCLUDE 'voro01.inc'
!     Parameter List
      INTEGER NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),
     '  NPNODE(0:NP_R_M,0:NRM),nr
      REAL*8 XNFV(-(NJM+1):NJM,NFVM),XP(NKM,NVM,NJM,NPM),ZNFVMSH(NFVM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER np,nfvl,cnp,nfv,nj,njj,nonode,cnonode,nvc
      REAL*8 VEL(3),DDOT,XJ_CIJ(3),CIJ_XI(3),TERM

      CALL ENTERS('MESHFLUX',*9999)

      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)

C       ..Loop over adjoining nodes
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)
          cnp=NPNODE(cnonode,nr)
          nfv=NFVC(2,nfvl,nvc)

          IF(CENTMESH) THEN

C           ..First make the node position minus the centroid part
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              XJ_CIJ(njj)=XP(1,1,nj,cnp)-XNFV(CENTINDX(njj),nfv)
              CIJ_XI(njj)=XNFV(CENTINDX(njj),nfv)-XP(1,1,nj,np)
            ENDDO

C           ..Now get the term in square brackets
            TERM=0.d0
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              TERM=TERM+XJ_CIJ(njj)*XP(1,2,nj,cnp)+
     '          CIJ_XI(njj)*XP(1,2,nj,np)
            ENDDO

C           ..Multiply by area and inv dist to get the centroid based
C           mesh flux
            ZNFVMSH(nfv)=TERM*XNFV(FAREA,nfv)*XNFV(IDIST,nfv)

          ELSE

C           ..Add the two adjoining velocities together to get a total
C           ..Note velocities stored in version 2 part of XP
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              VEL(njj)=0.5d0*(XP(1,2,nj,np)+XP(1,2,nj,cnp))
            ENDDO

C           ..Face velocity = average of 2 adjoining nodes
C           ..Mesh flux     = rate of volume swept by the face
C           = Speed of face x Area of face x direction
C           that face is moving in

            ZNFVMSH(nfv)=XNFV(FAREA,nfv)*
     '        DDOT(NJ_LOC(NJL_GEOM,0,nr),VEL,1,XNFV(1,nfv),1)
          ENDIF
        ENDDO
      ENDDO

      CALL EXITS('MESHFLUX')
      RETURN
 9999 CALL ERRORS('MESHFLUX',ERROR)
      CALL EXITS('MESHFLUX')
      RETURN 1
      END


