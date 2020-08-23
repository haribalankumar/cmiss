C**** CMISS Module FZAMG.F: Dummy AMG routines

      SUBROUTINE AMG_SET_PARMS(*)
      CALL FLAG_ERROR(0,'cm is not built with amg1r6')
      CALL ERRORIN('AMG_SET_PARMS')
      RETURN 1
      END

      SUBROUTINE AMG_COPY_MATRIX(*)
      CALL FLAG_ERROR(0,'cm is not built with amg1r6 library')
      CALL ERRORIN('AMG_COPY_MATRIX')
      RETURN 1
      END

      SUBROUTINE AMG_WRAPPER(*)
      CALL FLAG_ERROR(0,'cm is not built with amg1r6')
      CALL ERRORIN('AMG_WRAPPER')
      RETURN 1
      END
