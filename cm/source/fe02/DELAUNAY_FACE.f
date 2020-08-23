      LOGICAL FUNCTION DELAUNAY_FACE(NBJ,FACE,ne1,ne2,NPNE,nr,NVJE,
     '  XP,ZA)

C#### Subroutine: DELAUNAY_FACE
C###  Description:
C###    Checks to see whether the face of the tetrahedron or triangle
C###    is Delaunay by checking to see the node in ne1 corresponding
C###    (opposite) to FACE is within the circumsphere of ne2,
C###    then the edge is NOT Delaunay.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      REAL*8 SMALL
      PARAMETER(SMALL=1.d-12)
!     Parameter List
      INTEGER NBJ(NJM,NEM),FACE,ne1,ne2,NPNE(NNM,NBFM,NEM),nr,
     '  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM)
!     Local Variables
      INTEGER ne1_node,ne1_versn,nb1,njj,nj
      REAL*8 RSQ,DIST

      nb1=NBJ(1,ne1)
      ne1_node=NPNE(FACE,nb1,ne1)
      ne1_versn=NVJE(FACE,nb1,1,ne1)
      RSQ=0.d0
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        nj=NJ_LOC(NJL_GEOM,njj,nr)
        DIST=ZA(1,njj,1,ne2)-XP(1,ne1_versn,nj,ne1_node)
        RSQ=RSQ+DIST**2
      ENDDO
      IF((ZA(1,NJ_LOC(NJL_GEOM,0,nr)+1,1,ne2)-SMALL).GT.RSQ) THEN
        DELAUNAY_FACE=.FALSE.
      ELSE
        DELAUNAY_FACE=.TRUE.
      ENDIF

      RETURN
      END

