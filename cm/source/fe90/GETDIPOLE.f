      SUBROUTINE GETDIPOLE(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,n,
     '  nr,CENTRE,DIPOLE_CEN,DIPOLE_DIR,DIRECTION,TIME,ERROR,*)

C#### Subroutine: GETDIPOLE
C###  Description:
C###    GETDIPOLE returns the centre and direction of dipole n on
C###    region nr at a particular time.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM),n,nr
      REAL*8 CENTRE(*),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM),DIRECTION(*),TIME
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nj,nt
      REAL*8 TIMEPOS

      CALL ENTERS('GETDIPOLE',*9999)

C LKC 2-FEB-1999
      CALL ASSERT(USE_DIPOLE.GE.1,'Set USE_DIPOLE to 1',
     '  ERROR,*9999)

      IF(DIPOLE_CEN_NTIME(n,nr).GT.0) THEN
        CALL ASSERT(TIME.LE.DIPOLE_CEN(4,DIPOLE_CEN_NTIME(n,nr),n,nr)
     '    .AND.TIME.GE.DIPOLE_CEN(4,0,n,nr),
     '    '>>Time is outside dipole path',ERROR,*9999)
        nt=0
        DO WHILE(nt.LT.DIPOLE_CEN_NTIME(n,nr).AND.TIME.GT.
     '    DIPOLE_CEN(4,nt+1,n,nr))
          nt=nt+1
        ENDDO
        CALL ASSERT(nt.LT.DIPOLE_CEN_NTIME(n,nr),
     '    '>>Could not find dipole',ERROR,*9999)
        IF(DABS(DIPOLE_CEN(4,nt+1,n,nr)-DIPOLE_CEN(4,nt,n,nr)).GT.
     '    ZERO_TOL) THEN
          TIMEPOS=(TIME-DIPOLE_CEN(4,nt,n,nr))/
     '      (DIPOLE_CEN(4,nt+1,n,nr)-DIPOLE_CEN(4,nt,n,nr))
        ELSE
          TIMEPOS=TIME
        ENDIF
        DO nj=1,NJT
          CENTRE(nj)=(1.0d0-TIMEPOS)*DIPOLE_CEN(nj,nt,n,nr)+
     '      TIMEPOS*DIPOLE_CEN(nj,nt+1,n,nr)
        ENDDO !nj
      ELSE
        DO nj=1,NJT
          CENTRE(nj)=DIPOLE_CEN(nj,0,n,nr)
        ENDDO !nj
      ENDIF
      IF(DIPOLE_DIR_NTIME(n,nr).GT.0) THEN
        CALL ASSERT(TIME.LE.DIPOLE_DIR(4,DIPOLE_DIR_NTIME(n,nr),n,nr)
     '    .AND.TIME.GE.DIPOLE_DIR(4,0,n,nr),
     '    '>>Time is outside dipole path',ERROR,*9999)
        nt=0
        DO WHILE(nt.LT.DIPOLE_DIR_NTIME(n,nr).AND.TIME.GT.
     '    DIPOLE_DIR(4,nt+1,n,nr))
          nt=nt+1
        ENDDO
        CALL ASSERT(nt.LT.DIPOLE_DIR_NTIME(n,nr),
     '    '>>Could not find dipole',ERROR,*9999)
        IF(DABS(DIPOLE_DIR(4,nt+1,n,nr)-DIPOLE_DIR(4,nt,n,nr)).GT.
     '    ZERO_TOL) THEN
          TIMEPOS=(TIME-DIPOLE_DIR(4,nt,n,nr))/
     '      (DIPOLE_DIR(4,nt+1,n,nr)-DIPOLE_DIR(4,nt,n,nr))
        ELSE
          TIMEPOS=TIME
        ENDIF
        DO nj=1,NJT
          DIRECTION(nj)=(1.0d0-TIMEPOS)*DIPOLE_DIR(nj,nt,n,nr)+
     '      TIMEPOS*DIPOLE_DIR(nj,nt+1,n,nr)
        ENDDO !nj
      ELSE
        DO nj=1,NJT
          DIRECTION(nj)=DIPOLE_DIR(nj,0,n,nr)
        ENDDO !nj
      ENDIF

      CALL EXITS('GETDIPOLE')
      RETURN
9999  CALL ERRORS('GETDIPOLE',ERROR)
      CALL EXITS('GETDIPOLE')
      RETURN 1
      END


