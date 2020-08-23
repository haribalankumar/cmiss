      SUBROUTINE OPOBJE(ERROR,*)

C#### Subroutine: OPOBJE
C###  Description:
C###    OPOBJE outputs objects.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'obje00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,ND1,ND2,noobje,nopart,NOSEGM,NTPART

      CALL ENTERS('OPOBJE',*9999)

      DO noobje=1,NTOBJE
        NOSEGM=NSOBJE(1,noobje) !is segment number of complete object
        NTPART=NSOBJE(2,noobje) !is number of parts to object
C        CALL STRING_TRIM(OBJECT_NAME,IBEG,IEND)
        CALL STRING_TRIM(OBJECT_NAME(noobje),IBEG,IEND)
        WRITE(OP_STRING,'(/'' Object '',I2,'' ='',A,'' (segment '','
     '    //'''number '',I3,'//''') has '',I2,'' parts '')')
     '    noobje,OBJECT_NAME(noobje)(IBEG:IEND),NOSEGM,NTPART
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nopart=1,NTPART
          ND1=NSOBJE(2+nopart,noobje) !is 1st nd number of object part
          ND2=NSOBJE(3+nopart,noobje) !is 2nd nd number of object part
          WRITE(OP_STRING,'('' Part number '',I2,'' starts at nd='','
     '      //'I3,'//''' and ends at nd='',I3)') nopart,ND1,ND2
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO

      CALL EXITS('OPOBJE')
      RETURN
 9999 CALL ERRORS('OPOBJE',ERROR)
      CALL EXITS('OPOBJE')
      RETURN 1
      END


