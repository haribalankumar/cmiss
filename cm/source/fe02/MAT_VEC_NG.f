      SUBROUTINE MAT_VEC_NG(NITB,nr,A_VECTOR,B_VECTOR,C_VECTOR,XG,
     '  ERROR,*)

C#### Subroutine: MAT_VEC_NG
C###  Description:
C###    MAT_VEC_NG calculates direction cosines of undeformed
C###    material vectors at Gauss point ng.
C###    This routine assumes XG contains Gauss pt coordinates,
C###    microstructural orientations, and derivatives wrt Xi.

C     MPN rewritten 25-Apr-96

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 DXRCXI(3,3)

      CALL ENTERS('MAT_VEC_NG',*9999)

C     Compute derivatives of rc coords wrt Xi coords
      CALL DXRCDXI(NITB,nr,DXRCXI,XG,ERROR,*9999)

C     Calculate direction cosines
      CALL MAT_VEC(NITB,nr,A_VECTOR,B_VECTOR,C_VECTOR,DXRCXI,XG,
     '  ERROR,*9999)

      CALL EXITS('MAT_VEC_NG')
      RETURN
 9999 CALL ERRORS('MAT_VEC_NG',ERROR)
      CALL EXITS('MAT_VEC_NG')
      RETURN 1
      END


