      SUBROUTINE D_FIBRE_REF_VECS(NITB,nr,DXRCXI,D2XRCXI,
     '  D_F_VECTOR,D_G_VECTOR,D_H_VECTOR,F_VECTOR,G_VECTOR,H_VECTOR,
     '  ERROR,*)

C#### Subroutine: D_FIBRE_REF_VECS
C###  Description:
C###    D_FIBRE_REF_VECS calculates direction cosines of undeformed
C###    fibre reference vectors and their derivatives wrt xi coords.
C###    F_VECTOR coincides with the Xi1 base vector.
C###    G_VECTOR lies in the Xi1-Xi2 plane and is perpendicular to Xi1.
C###    H_VECTOR is perpendicular to the Xi1-Xi2 plane.
C###    This routines assumes that DXRCXI and D2XRCXI contain the first
C###    and second derivatives of rectangular Cartesian coordinates
C###    wrt xi coords.

C     MPN 25-Apr-96: written for fibre/imbric/sheet angle reference

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 DXRCXI(3,3),D2XRCXI(3,3,3),
     '  D_F_VECTOR(3,3),D_G_VECTOR(3,3),D_H_VECTOR(3,3),
     '  F_VECTOR(3),G_VECTOR(3),H_VECTOR(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni,njj

      CALL ENTERS('D_FIBRE_REF_VECS',*9999)

C     Initialise all vectors
      DO njj=1,3
        F_VECTOR(njj)=0.0d0
        G_VECTOR(njj)=0.0d0
        H_VECTOR(njj)=0.0d0
        DO ni=1,NITB
          D_F_VECTOR(njj,ni)=0.0d0
          D_G_VECTOR(njj,ni)=0.0d0
          D_H_VECTOR(njj,ni)=0.0d0
        ENDDO !ni
      ENDDO !nj1

C     F_VECTOR is the normalised undeformed Xi1 base vector
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        F_VECTOR(njj)=DXRCXI(njj,1)
        DO ni=1,NITB
          D_F_VECTOR(njj,ni)=D2XRCXI(njj,ni,1)
        ENDDO !ni
      ENDDO !njj
      CALL D_NORMALISE(NITB,D_F_VECTOR,F_VECTOR,ERROR,*9999)

      IF(NITB.GE.2) THEN !2D or 3D
C       H_VECTOR is the undeformed Xi1-Xi2 plane normal
C CS 4/9/97 for 2D case ensure that H_VECTOR is the same for all elems.
        IF(NJT.NE.2) THEN
          CALL D_CROSS(NITB,DXRCXI(1,1),DXRCXI(1,2),H_VECTOR,
     '      D2XRCXI(1,1,1),D2XRCXI(1,1,2),D_H_VECTOR)
          CALL D_NORMALISE(NITB,D_H_VECTOR,H_VECTOR,ERROR,*9999)
        ELSE
          H_VECTOR(1)=0.0d0
          H_VECTOR(2)=0.0d0
          H_VECTOR(3)=1.0d0
        ENDIF
C       G_VECTOR lies in the undeformed Xi1-Xi2 plane and is
C       normal to F_VECTOR
        CALL D_CROSS(NITB,H_VECTOR,F_VECTOR,G_VECTOR,
     '    D_H_VECTOR,D_F_VECTOR,D_G_VECTOR)
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Undeformed fibre reference vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    f_vector    g_vector    h_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO njj=1,3
          WRITE(OP_STRING,'(1X,3D12.3)')
     '      F_VECTOR(njj),G_VECTOR(njj),H_VECTOR(njj)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !ni
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('D_FIBRE_REF_VECS')
      RETURN
 9999 CALL ERRORS('D_FIBRE_REF_VECS',ERROR)
      CALL EXITS('D_FIBRE_REF_VECS')
      RETURN 1
      END


