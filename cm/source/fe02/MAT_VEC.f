      SUBROUTINE MAT_VEC(NITB,nr,A_VECTOR,B_VECTOR,C_VECTOR,
     '  DXRCXI,XG,ERROR,*)

C#### Subroutine: MAT_VEC
C###  Description:
C###    MAT_VEC calculates direction cosines of undeformed
C###    material vectors.
C###    This routine assumes XG contains microstructural orientations,
C###    and DXRCDXI contains derivatives of rc coords wrt Xi.

C     MPN rewritten 25-Apr-96

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),
     '  DXRCXI(3,3),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj
      REAL*8 FIBRE_ORIENT(3,3)

      CALL ENTERS('MAT_VEC',*9999)

C     Compute components of undeformed orthonormal fibre reference
C     vectors at Gauss pt NG wrt rc coord system
      CALL FIBRE_REF_VECS(NITB,nr,FIBRE_ORIENT(1,1),
     '  FIBRE_ORIENT(1,2),FIBRE_ORIENT(1,3),DXRCXI,ERROR,*9999)

C     Rotate fibre reference vectors into true material fibre vectors
      CALL MAT_VEC_ROTATE(NITB,nr,FIBRE_ORIENT,XG,ERROR,*9999)
      DO nj=1,3
        A_VECTOR(nj)=FIBRE_ORIENT(nj,1)
        B_VECTOR(nj)=FIBRE_ORIENT(nj,2)
        C_VECTOR(nj)=FIBRE_ORIENT(nj,3)
      ENDDO

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Undeformed material vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    a_vector    b_vector    c_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nj=1,3
          WRITE(OP_STRING,'(1X,3D12.3)')
     '      A_VECTOR(nj),B_VECTOR(nj),C_VECTOR(nj)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nj
CC$      call mp_unsetlock()
      ENDIF

C 21/2/97 LC archived section :  old MPN 25-Apr-96

CC**** Delta is angle between h (wall normal) and v (normal to wall in
CC**** yz-plane).
CC**** Phi is angle between v and r.
CC**** Note: If Beta>0 A_VECTOR is calculated with Beta terms included.


      CALL EXITS('MAT_VEC')
      RETURN
 9999 CALL ERRORS('MAT_VEC',ERROR)
      CALL EXITS('MAT_VEC')
      RETURN 1
      END


