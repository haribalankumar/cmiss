C**** MPA4000 dummy routines

      SUBROUTINE MBCLS(A)
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      DIMENSION A(*)
      CHARACTER ERROR*10
      WRITE(OP_STRING(1),*) ' >>Link with MAP 4000 library'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
9999  RETURN
      END

      SUBROUTINE MBOPN(A)
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      DIMENSION A(*)
      CHARACTER ERROR*10
      WRITE(OP_STRING(1),*) ' >>Link with MAP 4000 library'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
9999  RETURN
      END
