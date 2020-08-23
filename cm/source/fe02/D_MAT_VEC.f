      SUBROUTINE D_MAT_VEC(NITB,nr,DXRCXI,D2XRCXI,
     '  D_FIBRE_ORIENT,FIBRE_ORIENT,XG,ERROR,*)

C#### Subroutine: D_MAT_VEC
C###  Description:
C###    D_MAT_VEC calculates direction cosines of undeformed
C###    material vectors and their derivatives wrt local element xi
C###    coordinates.
C###    This routine assumes XG contains microstructural orientations,
C###    and DXRCXI and D2XRCXI contains derivatives of rc coords wrt xi.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 DXRCXI(3,3),D2XRCXI(3,3),
     '  D_FIBRE_ORIENT(3,3,3),FIBRE_ORIENT(3,3),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni,njj
      REAL*8 D_F_VECTOR(3,3),D_G_VECTOR(3,3),D_H_VECTOR(3,3)

      CALL ENTERS('D_MAT_VEC',*9999)

C     Compute components of undeformed orthonormal fibre reference
C     vectors at Gauss pt NG wrt rc coord system
      CALL D_FIBRE_REF_VECS(NITB,nr,DXRCXI,D2XRCXI,
     '  D_F_VECTOR,D_G_VECTOR,D_H_VECTOR,
     '  FIBRE_ORIENT(1,1),FIBRE_ORIENT(1,2),FIBRE_ORIENT(1,3),
     '  ERROR,*9999)
      DO ni=1,NITB
        DO njj=1,3
          D_FIBRE_ORIENT(njj,1,ni)=D_F_VECTOR(njj,ni)
          D_FIBRE_ORIENT(njj,2,ni)=D_G_VECTOR(njj,ni)
          D_FIBRE_ORIENT(njj,3,ni)=D_H_VECTOR(njj,ni)
        ENDDO
      ENDDO

C     Rotate fibre reference vectors into true material fibre vectors
      CALL D_MAT_VEC_ROTATE(NITB,nr,D_FIBRE_ORIENT,
     '  FIBRE_ORIENT,XG,ERROR,*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Undeformed material vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    a_vector    b_vector    c_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO njj=1,3
          WRITE(OP_STRING,'(1X,3D12.3)')
     '      FIBRE_ORIENT(njj,1),FIBRE_ORIENT(njj,2),FIBRE_ORIENT(njj,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nj
CC$      call mp_unsetlock()
      ENDIF

C 21/2/97 LC archived section :  old MPN 25-Apr-96

CC**** Delta is angle between h (wall normal) and v (normal to wall in
CC**** yz-plane).
CC**** Phi is angle between v and r.
CC**** Note: If Beta>0 A_VECTOR is calculated with Beta terms included.


      CALL EXITS('D_MAT_VEC')
      RETURN
 9999 CALL ERRORS('D_MAT_VEC',ERROR)
      CALL EXITS('D_MAT_VEC')
      RETURN 1
      END


