      SUBROUTINE OPCUST(ERROR,*)
C SMAR009 18/01/99 removed nr,
C#### Subroutine: OPCUST
C###  Description:
C###    OPCUST outputs customisation parameters for customising a mesh.

      IMPLICIT NONE
!     SMAR009 18/12/98 INCLUDE 'cmiss$reference:b00.cmn'
!     SMAR009 18/12/98 INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'chmesh0.cmn'
!     SMAR009 18/12/98 INCLUDE 'cmiss$reference:geom00.cmn'
!     SMAR009 18/12/98 INCLUDE 'cmiss$reference:inout00.cmn'
!     Parameter List
C SMAR009 18/01/99 removed INTEGER nr
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n

      CALL ENTERS('OPCUST',*9999)
      IF(CUSTOMISATION_TYPE.EQ.5) THEN
        WRITE(OP_STRING,'(/'' Polynomial and cosine functions used'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The number of polynomial coefficients '
     '  //'used ='',I2)')NPC
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The coefficient values are '',20D12.4)')
     '  (POLY_COEFFS(n),n=1,NPC)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The number of cosine coefficients '
     '  //'used ='',I2)')NCC
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The coefficient values are '',20D12.4)')
     '  (COS_COEFFS(n),n=1,NCC)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSEIF(CUSTOMISATION_TYPE.EQ.6) THEN
        WRITE(OP_STRING,'(/'' 2-d measurments used'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The number of 2-d measurements '
     '  //'used ='',I2)')NTDM
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO n=1,NTDM
          WRITE(OP_STRING,'('' Measurment number '',I3)')n
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' The measurement in the x direction '
     '    //'is '',D12.4)')WDMEASURE(1,n)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' The measurement in the y direction '
     '    //'is '',D12.4)')WDMEASURE(2,n)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' The height of the measurements '
     '    //'is '',D12.4,''/'')')WDMEASURE(3,n)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !n
        IF(SCLTYPE.EQ.1) THEN
          WRITE(OP_STRING,'('' Simple length scaling used'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Total length ='',I2)')TORSO_LENGTHS(0,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSEIF(SCLTYPE.EQ.2) THEN
          WRITE(OP_STRING,'('' Variable length scaling used'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Number of length measurements '
     '    //'used = '',I3)')NTL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO n=1,NTL
            WRITE(OP_STRING,'('' Measurment number '',I3)')n
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Height on generic and corresponding''
     '        //'' actual height is '',2(D12.4,2X) )')
     '        TORSO_LENGTHS(1,n),TORSO_LENGTHS(2,n)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !n
        ELSE
          ERROR='>> No such length scaling'
          GOTO 9999
        ENDIF
      ELSE
        ERROR='>>No such customisation type'
        GOTO 9999
      ENDIF !customisation type

      CALL EXITS('OPCUST')
      RETURN
 9999 CALL ERRORS('OPCUST',ERROR)
      CALL EXITS('OPCUST')
      RETURN 1
      END


