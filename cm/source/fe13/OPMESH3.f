      SUBROUTINE OPMESH3(ERROR,*)

C#### Subroutine: OPMESH3
C###  Description:
C###    OPMESH3 outputs circular or eccentric spheres mesh.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mesh03.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nocircle,nosphere,nr
      LOGICAL BEMREGION,FEMREGION

      CALL ENTERS('OPMESH3',*9999)

      WRITE(OP_STRING,'(/'' Number of regions = '',I2)') MESH3_NR
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      BEMREGION=.FALSE.
      FEMREGION=.FALSE.
      DO nr=1,MESH3_NR
        IF(MESH3_TYPE(nr).EQ.1) THEN
          BEMREGION=.TRUE.
          IF(nr.EQ.1) THEN
            IF(MESH3_BEMTYPE.EQ.1) THEN
              IF(NJT.EQ.2) THEN
                WRITE(OP_STRING,'('' Region '',I2,'' is a '
     '            //'circular BE region'')') nr
              ELSE IF(NJT.EQ.3) THEN
                WRITE(OP_STRING,'('' Region '',I2,'' is a '
     '            //'spherical BE region'')') nr
              ENDIF
            ELSE IF(MESH3_BEMTYPE.EQ.2) THEN
              WRITE(OP_STRING,'('' Region '',I2,'' is an '
     '          //'annular BE region'')') nr
            ENDIF
          ELSE
            WRITE(OP_STRING,'('' Region '',I2,'' is a BE region'')')
     '        nr
          ENDIF
        ELSE IF(MESH3_TYPE(nr).EQ.2) THEN
          FEMREGION=.TRUE.
          WRITE(OP_STRING,'('' Region '',I2,'' is a FE region'')') nr
        ENDIF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO !nr

      IF(NJT.EQ.2)THEN
        WRITE(OP_STRING,'('' Number of circles = '',I3)') NSPHERES
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(BEMREGION) THEN
          WRITE(OP_STRING,'( '' BE basis function number '','
     '      //'''is '',I2)') MESH3_NB(1,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(FEMREGION) THEN
          WRITE(OP_STRING,'( '' FE basis function number '','
     '      //'''is '',I2)') MESH3_NB(1,2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        DO nocircle=1,NSPHERES
          WRITE(OP_STRING,'(/'' Radius of circle '',I3,'' = '',D12.4)')
     '      nocircle,MESH3_RAD(nocircle)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Number of elements around '','
     '      //'''circle '',I3,'' = '',I3)')
     '      nocircle,MESH3_S(nocircle,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Conductivity inside circle '','
     '      //'I3,'' = '',D12.4)') nocircle,SIGMA(nocircle)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ELSEIF (NJT.EQ.3)THEN
        WRITE(OP_STRING,'('' Number of spheres = '',I3)') NSPHERES
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(BEMREGION) THEN
          WRITE(OP_STRING,'( '' Main (central) BE basis function '
     '      //'number is '',I2)') MESH3_NB(1,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Bottom (apex node 1) BE sector basis '
     '      //'function number is '',I2)') MESH3_NB(2,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Top (apex node 3) BE sector basis '
     '      //'function number is '',I2)') MESH3_NB(3,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(FEMREGION) THEN
          WRITE(OP_STRING,'( '' Main (central) BFE basis function '
     '      //'number is '',I2)') MESH3_NB(1,2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Bottom (apex node 1) FE sector basis '
     '      //'function number is '',I2)') MESH3_NB(2,2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Top (apex node 3) FE sector basis '
     '      //'function number is '',I2)') MESH3_NB(3,2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        DO nosphere=1,NSPHERES
          WRITE(OP_STRING,'(/'' Radius of sphere '',I2,'' = '',D12.4)')
     '      nosphere,MESH3_RAD(nosphere)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Number of elements around '','
     '      //'''sphere '',I2,'' (in the theta direction) = '',I3)')
     '      nosphere,MESH3_S(nosphere,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Number of elements in the azimuthal '','
     '      //'''direction of sphere '',I2,'' = '',I3)')
     '      nosphere,MESH3_S(nosphere,2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Conductivity inside sphere '','
     '      //'I2,'' = '',D12.4)') nosphere,SIGMA(nosphere)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF !njt loop

      CALL EXITS('OPMESH3')
      RETURN
 9999 CALL ERRORS('OPMESH3',ERROR)
      CALL EXITS('OPMESH3')
      RETURN 1
      END


