C**** CMISS Module FD16: Dummy Coronary angiography routines

      SUBROUTINE ANGIO(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module FE16: need sub ANGIO')
      RETURN 1
      END

      SUBROUTINE LIRECO(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module FE16: need sub LIRECO')
      RETURN 1
      END
