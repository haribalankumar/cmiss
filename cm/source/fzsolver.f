C**** CMISS Module fzsolver.f: Dummy routines for the solver library interface

C       SUBROUTINE ALLOC_SOLVER()
C       CALL FLAG_ERROR(0,'Link with Solver library')
C       RETURN 1
C       END

      SUBROUTINE CONVERT_SPARSE_5TO1()
      CALL FLAG_ERROR(0,
     &  'Link with Solver library: need sub CONVERT_SPARSE_5TO1')
      RETURN 1
      END

      SUBROUTINE FREE_SOLVER()
      RETURN
      END

      SUBROUTINE IPSOLU_SOLVER()
      CALL FLAG_ERROR(0,'Link with Solver library')
      CALL ERRORIN('IPSOLU_SOLVER')
      RETURN 1
      END

      SUBROUTINE OPSOLU_SOLVER()
      CALL FLAG_ERROR(0,'Link with Solver library')
      CALL ERRORIN('OPSOLU_SOLVER')
      RETURN 1
      END

      SUBROUTINE SOLVE_SYSTEM()
      CALL FLAG_ERROR(0,'Link with Solver library')
      CALL ERRORIN('SOLVE_SYSTEM')
      RETURN 1
      END
