      SUBROUTINE FIELDML_WRITE_FIELD_REFS_CELLML(IFILE,NUM_FIELDS,
     '  FIELD_NAMES)

C#### Subroutine: FIELDML_WRITE_FIELD_REFS_CELLML
C###  Description:
C###  For each of the FIELD_NAMES writes out the corresponding FieldML
C###  field_ref element for a single component real field.

      IMPLICIT NONE
      
      INCLUDE 'cell02.cmn'

      !Parameter list
      INTEGER IFILE,NUM_FIELDS
      CHARACTER FIELD_NAMES(*)*(CELL_NAME_LENGTH)
      !Local variables
      INTEGER field,IBEG,IEND

      DO field=1,NUM_FIELDS
        CALL STRING_TRIM(FIELD_NAMES(field),IBEG,IEND)
        WRITE(IFILE,'(''    <field_ref ref="'
     '    //FIELD_NAMES(field)(IBEG:IEND)//'">'')')
        WRITE(IFILE,'(''      <component_ref ref="value">'')')
        WRITE(IFILE,'(''        <label name="value"/>'')')
        WRITE(IFILE,'(''      </component_ref>'')')
        WRITE(IFILE,'(''    </field_ref>'')')
      ENDDO !field=1,NUM_FIELDS
      
      RETURN
      END
      

