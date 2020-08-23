      SUBROUTINE D_MAT_VEC_ROTATE(NITB,nr,D_FIBRE_ORIENT,FIBRE_ORIENT,
     '  XG,ERROR,*)

C#### Subroutine: D_MAT_VEC_ROTATE
C###  Description:
C###    MAT_VEC_ROTATE rotates the fibre reference vector orientations
C###    and their derivatives wrt xi into the true material fibre
C###    directions and their derivatives using the fibre, sheet and
C###    sheet-normal angles stored in XG.
C###  See-Also: FIBRE IMBRICATION and SHEET ANGLES

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 XG(NJM,NUM),D_FIBRE_ORIENT(3,3,3),FIBRE_ORIENT(3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni,nj,NU1(3)
      REAL*8 ALFA,BETA,GAMA,D_ALFA(3),D_BETA(3),D_GAMA(3)

      DATA NU1/2,4,7/

      CALL ENTERS('D_MAT_VEC_ROTATE',*9999)

C***  See MAT_VEC_ROTATE for an explanation of rotation angles

      IF(NJ_LOC(NJL_FIBR,0,nr).GE.1) THEN !fibre angles defined
        nj=NJ_LOC(NJL_FIBR,1,nr)
        ALFA=XG(nj,1)
        DO ni=1,NITB
          D_ALFA(ni)=XG(nj,NU1(ni))
        ENDDO
C       Rotate fibre ref vectors about the sheet-normal axis by ALFA
        CALL D_ROT_COORDSYS(3,NITB,ALFA,FIBRE_ORIENT,
     '    D_ALFA,D_FIBRE_ORIENT,ERROR,*9999)
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
        BETA=XG(nj,1)
        DO ni=1,NITB
          D_BETA(ni)=XG(nj,NU1(ni))
        ENDDO
C       Rotate resulting vectors about the new sheet axis by BETA
        CALL D_ROT_COORDSYS(2,NITB,BETA,FIBRE_ORIENT,
     '    D_BETA,D_FIBRE_ORIENT,ERROR,*9999)
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
        GAMA=XG(nj,1)
        DO ni=1,NITB
          D_GAMA(ni)=XG(nj,NU1(ni))
        ENDDO
C    CS 29/1/2001 Removing -PI/2 rotation as described below
C!!! TEMPORARY: MPN 5May96
C!!!      NOTE: this is temporary until we redefine the way the sheet
C!!!      angle is defined. Better to rotate by GAMA above and not
C!!!      GAMA-PI/2. This means the sheet angles for the full
C!!!      heart mesh will need transforming.
C!!!      See also below.

C       Rotate resulting vectors about the fibre axis by GAMA-PI/2
C        CALL D_ROT_COORDSYS(1,NITB,GAMA-PI/2.0d0,FIBRE_ORIENT,
C     '    D_GAMA,D_FIBRE_ORIENT,ERROR,*9999)
C       Rotate resulting vectors about the fibre axis by GAMA
        CALL D_ROT_COORDSYS(1,NITB,GAMA,FIBRE_ORIENT,
     '    D_GAMA,D_FIBRE_ORIENT,ERROR,*9999)
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
C    removed see comment above
C!!! see TEMPORARY comment above
C!!! Eventually remove this routine call
C        IF(NITB.EQ.3) THEN !only rotate for 3D elements
C          DO ni=1,NITB
C            D_GAMA(ni)=0.0d0
C          ENDDO
C          CALL D_ROT_COORDSYS(1,NITB,-PI/2.0d0,FIBRE_ORIENT,
C     '      D_GAMA,D_FIBRE_ORIENT,ERROR,*9999)
C        ENDIF
      ENDIF

      CALL EXITS('D_MAT_VEC_ROTATE')
      RETURN
 9999 CALL ERRORS('D_MAT_VEC_ROTATE',ERROR)
      CALL EXITS('D_MAT_VEC_ROTATE')
      RETURN 1
      END


