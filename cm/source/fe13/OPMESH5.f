      SUBROUTINE OPMESH5(ERROR,*)

C#### Subroutine: OPMESH5
C###  Description:
C###    OPMESH5 outputs open cylindrical mesh.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mesh05.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nocylinder

      CALL ENTERS('OPMESH5',*9999)
      IF(NJT.EQ.2)THEN
        ERROR='>> NO 2d cylindrical mesh set up'
      ELSE IF(NJT.EQ.3)THEN
        WRITE(OP_STRING,'(/'' Number of cylinders =  '',I2)')NCYLINDERS
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' Basis function# of bicubic hermite '','
     '    //'''is '',I2)') MESH5_NB
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Height of the cylinders = '',D12.4)')
     '    HEIGHT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nocylinder=1,NCYLINDERS
          WRITE(OP_STRING,'(/'' Radius of cylinder '',I2,'' = '''
     '      //',D12.4)') nocylinder,MESH5_RAD(nocylinder)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Number of elements around '','
     '      //'''cylinder '',I2,'' = '',I3)')
     '      nocylinder,MESH5_S(nocylinder,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Number of elements up '','
     '      //'''cylinder '',I2,'' = '',I3)')
     '      nocylinder,MESH5_S(nocylinder,2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF !njt loop

      CALL EXITS('OPMESH5')
      RETURN
 9999 CALL ERRORS('OPMESH5',ERROR)
      CALL EXITS('OPMESH5')
      RETURN 1
      END


