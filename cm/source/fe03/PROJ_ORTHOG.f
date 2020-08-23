      SUBROUTINE PROJ_ORTHOG(IBT,IDO,INP,NBJ,DISTANCE,
     '  XE,XI,ZD,FOUND,ERROR,*)

C#### Subroutine: PROJ_ORTHOG
C###  Description:
C###    <HTML>
C###    Finds the closest point on a given element to a given
C###    world coordinate (ZD).  Starts from the given XI position.
C###    If FOUND is false then the projection from the element to the
C###    coordinate must be orthogonal.  If a suitable projection is
C###    found, FOUND is set to true. The square of the scalar distance
C###    between the coordinate and the point is returned in DISTANCE.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM)
      REAL*8 DISTANCE,XE(NSM,NJM),XI(NIM),ZD(NJM)
      LOGICAL FOUND
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IT,ITMAX,NITB
      REAL*8 USER_TOL
      LOGICAL FOUND_E

      CALL ENTERS('PROJ_ORTHOG',*9999)
      !get deformed/undeformed into XE
      ITMAX=10 ! max # iterations to use
      NITB=NIT(NBJ(1)) ! find the dimension of the element
C     Loop over dimension of elements
      IF(NITB.EQ.1) THEN ! 1D elements
        IF(ITYP10(1).EQ.1) THEN !rect cartesian
          CALL CLOS11(IBT,IDO,INP,IT,ITMAX,NBJ,
     '      DISTANCE,XE,XI(1),ZD,FOUND,ERROR,*9999)
        ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar
        ELSE IF(ITYP10(1).EQ.3) THEN !spherical polar
        ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
        ELSE IF(ITYP10(1).EQ.5) THEN !oblate spheroidal
        ENDIF
        IF(.NOT.FOUND.AND.XI(1).GE.0.0D0.AND.XI(1).LE.1.0D0) THEN
          FOUND=.TRUE.
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' Xi= '',F7.3,'' SQ='',E12.3)') XI(1),DISTANCE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDIF
        IF(IT.GE.ITMAX) THEN
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' WARNING: Convergence not reached in CLOS1X'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDIF ! it.ge.itmax

      ELSE IF(NITB.EQ.2) THEN !2D elements
        IF(ITYP10(1).EQ.1) THEN !rectanglar cartesian
          CALL CLOS21(IBT,IDO,INP,IT,ITMAX,NBJ,
     '      DISTANCE,XE,XI,ZD,FOUND,ERROR,*9999)
        ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar
        ELSE IF(ITYP10(1).EQ.3) THEN !spherical polar
        ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
          CALL CLOS24(IBT,IDO,INP,IT,ITMAX,NBJ,
     '      DISTANCE,XE,XI,ZD,ERROR,*9999)
        ELSE IF(ITYP10(1).EQ.5) THEN !oblate spheroidal
        ENDIF
        IF(.NOT.FOUND.AND.XI(1).GE.0.0D0.AND.XI(1).LE.1.0D0.AND
     '    .XI(2).GE.0.0D0.AND.XI(2).LE.1.0D0) THEN
          FOUND=.TRUE.
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' Xi: '',2F7.3,'' SQ='',E12.3)') XI(1),XI(2),DISTANCE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDIF
        IF(IT.GE.ITMAX) THEN
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' WARNING: Convergence not reached in CLOS2X'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDIF ! it.ge.itmax

      ELSE IF(NITB.EQ.3) THEN !3D elements

C 24/2/97 LC archived section :  GMH 29/10/95 Not being done correctly

        IF(ITYP10(1).EQ.1) THEN !rectanglar cartesian
          FOUND_E=.FALSE.
          USER_TOL=LOOSE_TOL
          CALL CLOS31(IBT,IDO,INP,IT,ITMAX,NBJ,
     '      DISTANCE,USER_TOL,XE,XI,ZD,FOUND_E,FOUND,ERROR,*9999)
C new MPN 7Nov97: the following should be uncommented when
C                  CLOS3D is debugged.
CC new MPN 5Nov97: new CLOS3D routine for general curvilinear coords
C        ELSE
C          CALL CLOS3D(ITYP10(nr),IBT,IDO,INP,IT,ITMAX,NBJ,nr,
C     '      DISTANCE,XE,XI,ZD,ERROR,*9999)
C old MPN 7Nov97: the following should be commented when
C                 CLOS3D is debugged.
        ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar
        ELSE IF(ITYP10(1).EQ.3) THEN !spherical polar
        ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
        ELSE IF(ITYP10(1).EQ.5) THEN !oblate spheroidal
C end old
        ENDIF
        IF(.NOT.FOUND.AND.XI(1).GE.0.0D0.AND.XI(1).LE.1.0D0.AND
     '    .XI(2).GE.0.0D0.AND.XI(2).LE.1.0D0.AND
     '    .XI(3).GE.0.0D0.AND.XI(3).LE.1.0D0) THEN
          FOUND=.TRUE.
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' Xi: '',3F7.3,'' SQ='',E12.3)') XI(1),XI(2),
     '        XI(3),DISTANCE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDIF
        IF(IT.GE.ITMAX) THEN
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' WARNING: Convergence not reached in CLOS3X'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDIF ! it.ge.itmax
      ENDIF !nitb=1,2,3

      CALL EXITS('PROJ_ORTHOG')
      RETURN
 9999 CALL ERRORS('PROJ_ORTHOG',ERROR)
      CALL EXITS('PROJ_ORTHOG')
      RETURN 1
      END


