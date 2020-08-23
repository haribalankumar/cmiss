      SUBROUTINE FIBRE_REF_VECS(NITB,nr,F_VECTOR,G_VECTOR,H_VECTOR,
     '  dXRC_dXI,ERROR,*)

C#### Subroutine: FIBRE_REF_VECS
C###  Description:
C###    FIBRE_REF_VECS calculates direction cosines of undeformed
C###    fibre reference vectors at Gauss point ng.
C###    F_VECTOR coincides with the Xi1 base vector.
C###    G_VECTOR lies in the Xi1-Xi2 plane and is perpendicular to Xi1.
C###    H_VECTOR is perpendicular to the Xi1-Xi2 plane.
C###    This routines assumes that dXRC_dXI contains the derivatives of
C###    rectangular Cartesian coordinates w.r.t. to Xi coordinates.

C     MPN 25-Apr-96: written for fibre/imbric/sheet angle reference

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 F_VECTOR(3),G_VECTOR(3),H_VECTOR(3),dXRC_dXI(3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj
      LOGICAL ZERO_F_VECTOR

      CALL ENTERS('FIBRE_REF_VECS',*9999)

C     Initialise all vectors
      DO nj=1,3
        F_VECTOR(nj)=0.0d0
        G_VECTOR(nj)=0.0d0
        H_VECTOR(nj)=0.0d0
      ENDDO !nj1

C     F_VECTOR is the normalised undeformed Xi1 base vector
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        F_VECTOR(nj)=dXRC_dXI(nj,1)
      ENDDO !nj
      CALL NORMALISE(3,F_VECTOR,ERROR,*9999)

      IF(NITB.GE.2) THEN !2D or 3D
C news MPN 29June2000: check for zero length F_VECTOR
C                      and correct if necesary
        ZERO_F_VECTOR=.TRUE.
        DO nj=1,3
          IF(DABS(F_VECTOR(nj)).GT.ZERO_TOL) ZERO_F_VECTOR=.FALSE.
        ENDDO !nj
        IF(ZERO_F_VECTOR) THEN
C         ...so set G_VECTOR to be the normalised undef Xi2 base vector
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            G_VECTOR(nj)=dXRC_dXI(nj,2)
          ENDDO !nj1
          CALL NORMALISE(3,G_VECTOR,ERROR,*9999)
C         ...then F_VECTOR is the undeformed Xi2-Xi3 plane normal
          CALL CROSS(dXRC_dXI(1,2),dXRC_dXI(1,3),F_VECTOR)
          CALL NORMALISE(3,F_VECTOR,ERROR,*9999)
C         ...and H_VECTOR lies in the undeformed Xi2-Xi3 plane and is
C         normal to G_VECTOR
          CALL CROSS(F_VECTOR,G_VECTOR,H_VECTOR)
        ELSE !F_VECTOR is not all zero
C         H_VECTOR is the undeformed Xi1-Xi2 plane normal
C CS 4/9/97 for 2D case ensure that H_VECTOR is the same for all elems.
          IF(NJT.NE.2) THEN
            CALL CROSS(dXRC_dXI(1,1),dXRC_dXI(1,2),H_VECTOR)
            CALL NORMALISE(3,H_VECTOR,ERROR,*9999)
          ELSE
            H_VECTOR(1)=0.0d0
            H_VECTOR(2)=0.0d0
            H_VECTOR(3)=1.0d0
          ENDIF
C         G_VECTOR lies in the undeformed Xi1-Xi2 plane and is
C         normal to F_VECTOR
          CALL CROSS(H_VECTOR,F_VECTOR,G_VECTOR)
        ENDIF !ZERO_F_VECTOR
C newe MPN 29June2000
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Undeformed fibre reference vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    f_vector    g_vector    h_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nj=1,3
          WRITE(OP_STRING,'(1X,3D12.3)')
     '      F_VECTOR(nj),G_VECTOR(nj),H_VECTOR(nj)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !ni
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('FIBRE_REF_VECS')
      RETURN
 9999 CALL ERRORS('FIBRE_REF_VECS',ERROR)
      CALL EXITS('FIBRE_REF_VECS')
      RETURN 1
      END


