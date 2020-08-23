      SUBROUTINE FACE_INT_PREP(NEA,NEAXT,NF_EA,NHE,NHF,NIEF,
     '  NPF,nr,nx,INTEGRATE,ERROR,*)

C#### Subroutine: FACE_INT_PREP
C###  Description:
C###    Finds element connectivity and xi-direction information for a
C###    face and determines whether an integral needs to be evaluated.
C#### Variable: NIEF(0:ni_f)
C###  Type: INTEGER
C###  Set_up: FACE_INT_PREP
C###  Description:
C###    NIEF(ni_f) is the element xi direction corresponding to face
C###    xi-direction ni.
C###    NIEF(0) is the out-of-face element xi-direction.
C#### Variable: NEAXT
C###  Type: INTEGER
C###  Set_up: FACE_INT_PREP
C###  Description:
C###    NEAXT is the number of elements adjacent to the face. (<=2)
C#### Variable: NEA(neax)
C###  Type: INTEGER
C###  Set_up: FACE_INT_PREP
C###  Description:
C###    NEA(neax) are the adjacent elements to the face.
C#### Variable: NF_EA(neax)
C###  Type: INTEGER
C###  Set_up: FACE_INT_PREP
C###  Description:
C###    NF_EA(neax) are the local face numbers in the adjacent elements.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NEA(2),NEAXT,NF_EA(2),NHE(NEM),
     '  NHF,NIEF(0:2),NPF(9),nr,nx
      LOGICAL INTEGRATE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NJF

      CALL ENTERS('FACE_INT_PREP',*9999)

      NJF=NJ_LOC(NJL_GEOM,0,nr)
      CALL ASSERT(NJF.EQ.3,'>>Only 3D face integrals are implemented',
     '  ERROR,*9999)
      IF(NJF.EQ.1) THEN
        NIEF(0)=1
C        ne=NEP
      ELSE IF(NJF.EQ.2) THEN
C        NIEF(1)=NPL(1)
C        NIEF(0)=3-NIEF(1)
C        NEA(1)=NEL(1)
      ELSE !IF(NJF.EQ.3) THEN
        NIEF(1)=NPF(1) !first in-face xi direction
        NIEF(2)=NPF(3) !second in-face xi direction
        NIEF(0)=6-NIEF(1)-NIEF(2) !out-of-face xi direction
        NEAXT=NPF(5) !number of adjacent elements
        NEA(1)=NPF(6) !adjacent element numbers
        NEA(2)=NPF(7)
        NF_EA(1)=NPF(8) !local face number in element ne1
        NF_EA(2)=NPF(9) !local face number in element ne2
      ENDIF !NJF

      NHF=NHE(NEA(1)) !num. dependent variables
      IF(ITYP15(nr,nx).EQ.3) THEN !Derivative Discontinuity Term
        INTEGRATE=.TRUE.
      ELSE !ITYP15(nr,nx)=1or2 !Flux Term
C       If this is an interior face flux term then fluxes on either side
C       of the face will cancel unless derivatives are discontinuous.
        INTEGRATE=NEAXT.NE.2 !for now
C KAT 25May99: Weights are now zero where derivatives are continuous.
C        IF(.NOT.INTEGRATE) THEN !2 adjacent elements
C          DO nhx=1,NHF
C            nh=NH_LOC(nhx,nx)
C            nb_f=NBHF(nh)
C            DO neax=1,2 !loop over adjacent elements
C              nb_e(neax)=NBH(nh,1,NEA(neax))
C            ENDDO !neax
C            DO nn_f=1,NNT(nb_f)
C              DO neax=1,2 !loop over adjacent elements
C                nn_e=NNF(1+nn_f,NF_EA(neax),nb_e(neax))
C                nv(neax)=NVHE(nn_e,nb_e(neax),nh,NEA(neax))
C              ENDDO !neax
C              IF(nv(1).NE.nv(2)) INTEGRATE=.TRUE.
C            ENDDO !nn_f
C          ENDDO !nhx
C        ENDIF !NEAXT=2
      ENDIF

      CALL EXITS('FACE_INT_PREP')
      RETURN
 9999 CALL ERRORS('FACE_INT_PREP',ERROR)
      CALL EXITS('FACE_INT_PREP')
      RETURN 1
      END


