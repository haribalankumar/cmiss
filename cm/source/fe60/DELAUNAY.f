      SUBROUTINE DELAUNAY(nb,NBJ,NEELEM,NENP,NKJE,NPLIST,NPNE,
     '  NPNODE,nr,NRE,NVJE,N_BDRY,N_IBDRY,N_INTNL,SE,
     '  XP,INTERNAL,ERROR,*)

C#### Subroutine: DELAUNAY
C###  Description:
C###    Generates a 3-dimensional Delaunay mesh from
C###    boundary information.
C###
C MHT 6-12-01 Chris Were's GENMESH program, converted from C to
C   Fortran by John Bodley.  Will remain commented out until converted
C     to proper CMISS format.
      
C#### Variable: NODE(LDNODE,3)
C###  Type: REAL*8
C###  Set_up: MAKE_TETRAHEDRAL
C###  Description:
C###    NODE Contains the coordinates of the tetrahedral vertices

C#### Variable: CIRC(LDCIRC,3)
C###  Type: REAL*8
C###  Set_up:
C###  Description:
C###    CIRC Contains the n-dimensional circumcentre of the tetrahedral

C#### Variable: TETRA(0:LDTETRA,4)
C###  Type: INTEGER
C###  Set_up:
C###  Description:
C###    TETRA Contains tetrahedral structural information
C###
C###         TETRA(0,HEAD) Delaunay mesh head tetrahedral
C###         TETRA(0,TAIL) Delaunay mesh tail tetrahedral
C###         TETRA(0,NEXT) Pointer to next tetrahedral
C###
C###         HEAD : 1 - Pointer to tetrahedral head
C###         TAIL : 2 - Pointer to tetrahedral tail
C###         PREV : 3 - Relative pointer to previous tetrahedral
C###         NEXT : 4 - Relative pointer to next tetrahedral

C#### Variable: ELEM(LDELEM,4)
C###  Type: INTEGER
C###  Set_up: MAKE_TETRAHEDRAL
C###  Description:
C###    ELEM Contains the global node number for each tetrahedral vertex

C#### Variable: FACE(LDFACE,NFACES)
C###  Type: INTEGER
C###  Set_up:
C###  Description:
C###    FACE Contains a pointer to the neighbouring tetrahedral on a
C###    given face

C#### Variable: FOUND(0:N_GM)
C###  Type: INTEGER
C###  Set_up: FOUND_LIST
C###  Description:
C###   FOUND Contains a list of pointers to the tetrahedral

C**** Variable: INTERNAL
C***  Type: LOGICAL
C***  Set_up:
C***  Description:
C***    INTERNAL TRUE if additional internal nodes are allowed to be
C***    added 

C#### Variable: INODE(N_GM)
C###  Type: INTEGER
C###  Set_up:
C###  Description:
C###   INODE Contains boundary node indices

C#### Variable: ITYPE(N_GM)
C###  Type: INTEGER
C###  Set_up:
C###  Description:
C###    ITYPE Contains boundary node type

C**** Variable: BC(LDBC,3)
C***  Type: REAL*8
C***  Set_up:
C***  Description:
C***    BC Contains boundary information -> Looks as though this is not used

      IMPLICIT NONE
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'genmesh.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),N_BDRY,N_IBDRY,N_INTNL
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      LOGICAL INTERNAL
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER BDRY(NBOUNDM),ELEM(LDELEM,4),FACE(LDFACE,NFACES),
     '  FOUND(0:N_GM),I,INODE(N_GM),IREGION(0:N_GM),ITYPE(N_GM),J,
     '  TETRA(0:LDTETRA,4)
      REAL*8 BC(LDBC,3),CIRC(LDCIRC,3),GAPSQ,NODE(LDNODE,3),
     '  RSQ(N_GM)
      LOGICAL LGAPSQ,PLANE

      CALL ENTERS('DELAUNAY',*9999)

      IF(NJT.EQ.2)THEN
        NPTS=3
      ELSE IF(NJT.EQ.3)THEN
        NPTS=4
      ENDIF
      GAPSQ=1.d0 !should be controlled from command line
      DO i=0,LDTETRA
        DO j=1,4
          TETRA(i,j)=0
        ENDDO !j
      ENDDO !i
 ! Calculate the 'sets' of B-IB nodes.  i.e. corners will have >1 B node
      CALL MAKE_SETS_DELAUNAY(BDRY,NPLIST,N_BDRY,N_IBDRY,XP,ERROR,*9999)
      ! Read in the initial geometry
      CALL READ_GEOMETRY(BDRY,INODE,IREGION,ITYPE,N_BDRY,N_IBDRY,
     '  N_INTNL,NPLIST,BC,NODE,XP,LGAPSQ,ERROR,*9999)
      ! Make the initial tetrahedra
      CALL MAKE_TETRAHEDRAL(ELEM,FACE,IREGION,TETRA,CIRC,NODE,RSQ,
     '  PLANE,ERROR,*9999)
      DO J=1,3 ! Regions - excluding the initial tetrahedral
        DO I=1,IREGION(0)
          IF(IREGION(I).EQ.J) THEN
            CALL ASSEMBLE_VORO(ELEM,FACE,FOUND,I,TETRA,CIRC,NODE,RSQ,
     '        ERROR,*9999)
            CALL ADD_NODE(ELEM,FACE,FOUND,I,TETRA,CIRC,NODE,RSQ,
     '        PLANE,ERROR,*9999)
            IF(PLANE) THEN
              CALL FIX_TETRAHEDRAL(ELEM,FACE,TETRA,CIRC,NODE,RSQ,
     '          PLANE,ERROR,*9999)
              CALL CHECK_INTEGRITY(ELEM,FACE,TETRA,CIRC,NODE,RSQ,
     '          ERROR,*9999)
            ENDIF !PLANE
          ENDIF !IREGION
        ENDDO !I
      ENDDO !J
      ! Peel temporary tetrahedral
      CALL PEEL_AWAY(ELEM,FACE,IREGION,TETRA,ERROR,*9999)
      IF(INTERNAL.AND.LGAPSQ) THEN
 ! Calculate internal points
        CALL INTERNAL_NODES(ELEM,FACE,FOUND,IREGION,TETRA,CIRC, GAPSQ,
     &    NODE,RSQ,ERROR,*9999)
      ENDIF !INTERNAL
      ! Write out the Delaunay mesh
      CALL DELAUNAY_NE(nb,NBJ,NEELEM,NENP,NKJE,NPNE,NPNODE,nr,NRE,NVJE,
     '  ELEM,TETRA,SE,ERROR,*9999)

      CALL EXITS('DELAUNAY')
      RETURN
 9999 CALL ERRORS('DELAUNAY',ERROR)
      CALL EXITS('DELAUNAY')
      RETURN 1
      END


