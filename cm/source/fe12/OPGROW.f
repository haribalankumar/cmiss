      SUBROUTINE OPGROW(ERROR,*)

C#### Subroutine: OPGROW
C###  Description:
C###    OPGROW outputs growth law parameters.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'grow00.cmn'
      INCLUDE 'ktyp60.cmn'
!     Parameter List
      CHARACTER ERROR*(*)

      CALL ENTERS('OPGROW',*9999)
      IF(KTYP60.EQ.1) THEN
        WRITE(OP_STRING,'('' Growth law:'','
     '    //''' Density = k* strain energy density'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Growth proportionality constant is '',D12.4)') GROW2
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ELSE IF(KTYP60.EQ.2) THEN
        WRITE(OP_STRING,'('' Growth law:'','
     '    //''' Density = constant + k* strain energy density'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Growth density constant         is '',D12.4)') GROW1
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Growth proportionality constant is '',D12.4)') GROW2
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ELSE IF(KTYP60.EQ.3) THEN !Carter growth law

      ELSE IF(KTYP60.EQ.4) THEN !Huiskes growth law

      ELSE IF(KTYP60.EQ.5) THEN !Tree growth law #1
        WRITE(OP_STRING,'('' Growth law:'','
     '    //''' Length change = constant - k* oxygen gradient'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Growth length change constant   is '',D12.4)') GROW1
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Growth proportionality constant is '',D12.4)') GROW2
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('OPGROW')
      RETURN
 9999 CALL ERRORS('OPGROW',ERROR)
      CALL EXITS('OPGROW')
      RETURN 1
      END


