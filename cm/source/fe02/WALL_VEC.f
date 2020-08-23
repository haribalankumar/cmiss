      SUBROUTINE WALL_VEC(NITB,nr,W1_VECTOR,W2_VECTOR,W3_VECTOR,
     '  DXRCXI,XG,ERROR,*)

C#### Subroutine: WALL_VEC
C###  Description:
C###    WALL_VEC calculates direction cosines of undeformed
C###    wall vectors at Gauss point ng.
C###    W1_VECTOR coincides with the circumferential dirn base vector.
C###    W2_VECTOR lies in the Xi1-Xi2 plane, perpendicular to W1_VECTOR.
C###    W3_VECTOR is perpendicular to the Xi1-Xi2 plane.
C###    This routine assumes that XG contains the Gauss pt
C###    position/derivs wrt Xi and DXRCDXI contains derivatives
C###    of rc coords wrt Xi.

C     MPN 16-Jun-97: written for LIST STRAIN/STRESS output
C     THIS ROUTINE HAS NOT BEEN FULLY DEBUGGED.
C     CS 16-Feb-01 debugged some more, now works for rc coords

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 DXRCXI(3,3),W1_VECTOR(3),W2_VECTOR(3),W3_VECTOR(3),
     '  XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj1,njj_theta !SMAR009 22/12/98 ,nj2
      REAL*8 DZX,X(3)       !dXrc_dXref,

      CALL ENTERS('WALL_VEC',*9999)

C     Initialise all vectors and
C     put undef coords into X for DZX function call below
      DO nj1=1,3
        W1_VECTOR(nj1)=0.0d0
        W2_VECTOR(nj1)=0.0d0
        W3_VECTOR(nj1)=0.0d0
        IF(nj1.LE.NJ_LOC(NJL_GEOM,0,nr)) THEN
          X(nj1)=XG(nj1,1)
        ELSE
          X(nj1)=0.0d0
        ENDIF
      ENDDO !nj1

C     W1_VECTOR is the normalised undeformed circumferential base vector
C     CS 16/2/2001 fixing rc case
      IF(ITYP10(nr).EQ.1) THEN !rc
        DO nj1=1,NJ_LOC(NJL_GEOM,0,nr)
          W1_VECTOR(nj1)=DXRCXI(nj1,1)
        ENDDO !nj1
C!!!    MPN March2010, following comments by me stuff up, does not work
C!!!    for rc with long-axis = X        
C!!!    MPN Feb2004: long axis of LV is assumed to be the z-axis and 
C!!!    circumferential dirn is in X-Y plane (thus set Z component to zero)
C!!!    This could eventually be generalised for a specified long axis
C        W1_VECTOR(3)=0.0d0
      ELSE IF((ITYP10(nr).EQ.2).OR.(ITYP10(nr).EQ.3)) THEN !cylindrical/spherical polar
        njj_theta=2
      ELSE IF(ITYP10(nr).GE.4) THEN !prolate/oblate sph.
        njj_theta=3
      ENDIF
      IF(ITYP10(nr).NE.1) THEN
        DO nj1=1,NJ_LOC(NJL_GEOM,0,nr)
          W1_VECTOR(nj1)=DZX(ITYP10(nr),nj1,njj_theta,X)
        ENDDO !nj1
      ENDIF
      CALL NORMALISE(3,W1_VECTOR,ERROR,*9999)

      IF(NITB.GE.2) THEN !2D or 3D
C       W3_VECTOR is the undeformed Xi1-Xi2 plane normal
        CALL CROSS(DXRCXI(1,1),DXRCXI(1,2),W3_VECTOR)
        CALL NORMALISE(3,W3_VECTOR,ERROR,*9999)
C       W2_VECTOR lies in the undeformed Xi1-Xi2 plane and is
C       normal to W1_VECTOR
        CALL CROSS(W3_VECTOR,W1_VECTOR,W2_VECTOR)
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Undeformed wall vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''   w1_vector   w2_vector   w3_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nj1=1,3
          WRITE(OP_STRING,'(1X,3D12.3)')
     '      W1_VECTOR(nj1),W2_VECTOR(nj1),W3_VECTOR(nj1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nj1
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('WALL_VEC')
      RETURN
 9999 CALL ERRORS('WALL_VEC',ERROR)
      CALL EXITS('WALL_VEC')
      RETURN 1
      END


