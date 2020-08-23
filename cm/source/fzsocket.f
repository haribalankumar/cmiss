C**** CMISS Module FZSOCKET.F: Dummy SOCKET routines

      INTEGER FUNCTION FSKC2F(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with Socket library')
      FSKC2F=0
      RETURN
      END

      INTEGER FUNCTION FSKCLOSE(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with Socket library')
      FSKCLOSE=0
      RETURN
      END

      INTEGER FUNCTION FSKCONNECT(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with Socket library')
      FSKCONNECT=0
      RETURN
      END

      INTEGER FUNCTION FSKF2C(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with Socket library')
      FSKF2C=0
      RETURN
      END

      INTEGER FUNCTION FSKLEN(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with Socket library')
      FSKLEN=0
      RETURN
      END

      INTEGER FUNCTION FSKREAD(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with Socket library')
      FSKREAD=0
      RETURN
      END

      INTEGER FUNCTION FSKSELECT(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with Socket library')
      FSKSELECT=0
      RETURN
      END

      INTEGER FUNCTION FSKWRITE(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with Socket library')
      FSKWRITE=0
      RETURN
      END
