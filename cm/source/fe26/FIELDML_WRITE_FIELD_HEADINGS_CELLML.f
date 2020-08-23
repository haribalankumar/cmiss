      SUBROUTINE FIELDML_WR_FLD_HEADINGS_CELLML(IFILE,NUM_FIELDS,
     '  FIELD_NAMES)

C#### Subroutine: FIELDML_WR_FLD_HEADINGS_CELLML
C###  Description:
C###  For each of the FIELD_NAMES writes out the corresponding FieldML
C###  field element for a single component real field.

      IMPLICIT NONE
      
      INCLUDE 'cell02.cmn'

      !Parameter list
      INTEGER IFILE,NUM_FIELDS
      CHARACTER FIELD_NAMES(*)*(CELL_NAME_LENGTH)
      !Local variables
      INTEGER field,IBEG,IEND

      DO field=1,NUM_FIELDS
        WRITE(IFILE,'(''  <field '')')
        CALL STRING_TRIM(FIELD_NAMES(field),IBEG,IEND)
        WRITE(IFILE,'(''    name="'//FIELD_NAMES(field)(IBEG:IEND)//
     '    '"'')')
        WRITE(IFILE,'(''    value_type="real"'')')
        WRITE(IFILE,'(''    coordinate_system='//
     '    '"rectangular cartesian">'')')
        WRITE(IFILE,'(''    <component name="value"/>'')')
        WRITE(IFILE,'(''  </field>'')')
      ENDDO !field=1,NUM_FIELDS
      
      RETURN
      END
      

