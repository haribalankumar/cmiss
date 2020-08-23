C**** CMISS Module FZLAPACK.F: Dummy LAPACK routines

      SUBROUTINE DP0TRS(A,*)
      DIMENSION A(*)
      WRITE(*,*) '>>Link with Lapack library'
      RETURN 1
      END

      SUBROUTINE DTRS(A,*)
      DIMENSION A(*)
      WRITE(*,*) '>>Link with BLAS library'
      RETURN 1
      END

