      SUBROUTINE MAT_VEC_ROTATE(NITB,nr,FIBRE_ORIENT,XG,ERROR,*)

C#### Subroutine: MAT_VEC_ROTATE
C###  Description:
C###    MAT_VEC_ROTATE rotates the fibre reference vector orientations
C###    into the true material fibre directions using the fibre, sheet
C###    and sheet-normal angles stored in XG.

C     MPN 25-Apr-96: written to reflect Laminar Structure
C                    of the Heart paper

C#### Comment: FIBRE IMBRICATION and SHEET ANGLES
C###  Description:
C###    Eqtn 5 of "Laminar Structure of the Heart" by
C###    Le Grice, Hunter and Smaill defines how the fibrous-sheet
C###    material vectors can be expressed in terms of fibre reference
C###    coordinates. The base vectors of this coord system,
C###    F_VECTOR, G_VECTOR and H_VECTOR, are orthonormal and are
C###    computed from the XI base vectors (which are not orthonormal
C###    in general). Components of the fibrous-sheet material vectors
C###    are computed with respect to global cartesian coordinates and
C###    stored in A_VECTOR, B_VECTOR and C_VECTOR.
C###    Three successive coordinate rotations (the order IS important)
C###    are required to transform components . Firstly, the
C###    fibrous-sheet material vectors are rotated by GAMA
C###    about the fibre axis until the new sheet
C###    axis lies in the Xi1-Xi2 plane. Secondly, the resulting vectors
C###    are rotated by BETA about the new sheet axis until the new
C###    fibre axis also lies in the Xi1-Xi2 plane. Lastly, the
C###    resulting vectors are rotated by ALFA about the new sheet-normal
C###    axis (which now coincides with the Xi1-Xi2 plane normal) until
C###    the new fibre axis coincides with the Xi1 base vector.

C**** ALFA is fibre angle.
C**** BETA is imbrication angle.
C**** GAMA is sheet angle.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 XG(NJM,NUM),FIBRE_ORIENT(3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj
      REAL*8 ALFA,BETA,GAMA

      CALL ENTERS('MAT_VEC_ROTATE',*9999)

C     To compute undeformed material vectors from the fibre reference
C     apply the inverse of the transformation described above.
      IF(NJ_LOC(NJL_FIBR,0,nr).GE.1) THEN !fibre angles defined
        nj=NJ_LOC(NJL_FIBR,1,nr)
        ALFA=XG(nj,1) !fibre angle subtends the Xi1 base vector
C                     !and fibre axis in the Xi1-Xi2 plane.
C                     !ie) angle that rotates the Xi1 base vector
C                     !(about the sheet-normal axis) into the
C                     !direction of the fibre axis given that
C                     !the sheet axis lies in the Xi1-Xi2 plane
C                     !and the sheet-normal coincides with the
C                     !Xi1-Xi2 plane normal.
C       Rotate fibre ref vectors about the sheet-normal axis by ALFA
        CALL ROT_COORDSYS(3,ALFA,FIBRE_ORIENT,ERROR,*9999)
      ELSE
        IF(DOP.AND.NITB.GE.2) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' >>>WARNING: Assuming fibre '
     '      //'angle is zero'')')
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ENDIF

      IF(NJ_LOC(NJL_FIBR,0,nr).GE.2) THEN !imbrication angles defined
        nj=NJ_LOC(NJL_FIBR,2,nr)
        BETA=XG(nj,1) !imbrication angle subtends the fibre axis
C                     !and the Xi1-Xi2 plane (about the sheet axis)
C                     !given that the sheet axis lies in
C                     !the Xi1-Xi2 plane
C       Rotate resulting vectors about the new sheet axis by BETA
        CALL ROT_COORDSYS(2,BETA,FIBRE_ORIENT,ERROR,*9999)
      ELSE
        IF(DOP.AND.NITB.EQ.3) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' >>>WARNING: Assuming imbrication '
     '      //'angle is zero'')')
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ENDIF

      IF(NJ_LOC(NJL_FIBR,0,nr).GE.3) THEN !sheet angles defined
        nj=NJ_LOC(NJL_FIBR,3,nr)
C    CS 29/1/2001 Removed -PI/2 rotation
C        GAMA=XG(nj,1) !sheet angle-PI/2 subtends sheet axis and the
        GAMA=XG(nj,1) !sheet angle subtends sheet axis and the
C                     !Xi1-Xi2 plane about the true fibre axis
C    CS 29/1/2001 Removed -PI/2 rotation
C!!! TEMPORARY: MPN 5May96
C!!!      NOTE: this is temporary until we redefine the way the sheet
C!!!      angle is defined. Better to rotate by GAMA above and not
C!!!      GAMA-PI/2. This means the sheet angles for the full
C!!!      heart mesh will need transforming.
C!!!      See also below.

C         Rotate resulting vectors about the fibre axis by GAMA-PI/2
C          CALL ROT_COORDSYS(1,GAMA-PI/2.0d0,FIBRE_ORIENT,ERROR,*9999)

C!!! see TEMPORARY comment above
C!!! Eventually replace above call with this routine call
C    CS 30/1/2001 done
CC         Rotate resulting vectors about the fibre axis by GAMA
          CALL ROT_COORDSYS(1,GAMA,FIBRE_ORIENT,ERROR,*9999)

      ELSE
        IF(DOP.AND.NITB.EQ.3) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' >>>WARNING: Assuming sheet '
     '      //'angle is zero'')')
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
C!!! see TEMPORARY comment above
C!!! Eventually remove this routine call
C        IF(NITB.EQ.3) THEN !only rotate for 3D elements
C          CALL ROT_COORDSYS(1,-PI/2.0d0,FIBRE_ORIENT,ERROR,*9999)
C        ENDIF
      ENDIF

      CALL EXITS('MAT_VEC_ROTATE')
      RETURN
 9999 CALL ERRORS('MAT_VEC_ROTATE',ERROR)
      CALL EXITS('MAT_VEC_ROTATE')
      RETURN 1
      END


