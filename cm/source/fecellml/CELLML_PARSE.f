      SUBROUTINE CELLML_PARSE(variant,ERROR_CODE)

C#### Subroutine: CELLML_PARSE
C###  Description:
C###  <p>CELLML_PARSE takes a variant number and parses the URI given by
C###  that variant's entry in the CELLML_URIS array. ERROR_CODE is
C###  returned with a non-zero value if an error occurs while parsing
C###  the file. CELLML_INITIALISE must be called before a call to this
C###  routine.</p><p>Error codes: <ol>
C###  <li>Zero length URI specified for given variant;</li>
C###  </ol></p>

      IMPLICIT NONE

      INCLUDE 'cellml.cmn'

      !Parameter list
      INTEGER variant,ERROR_CODE
      !Local variables
      INTEGER IBEG,IEND,LENGTH

      !Initialise error code
      ERROR_CODE = 0

      CALL STRING_TRIM(CELLML_URIS(variant),IBEG,IEND)
      CALL CellMLProcessorParse(CELLML_URIS(variant)(IBEG:IEND),
     '  CELLML_MODELS(variant),ERROR_CODE)

      RETURN
      END


