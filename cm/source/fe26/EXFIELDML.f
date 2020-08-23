      SUBROUTINE EXFIELDML(LIST,NQLIST,NXLIST,XQ,YQS,
     '  CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,CELL_ICQS_NAMES,
     '  CELL_RCQS_VALUE,CELL_RCQS_SPATIAL,CELL_RCQS_NAMES,
     '  CELL_YQS_VALUE,CELL_YQS_SPATIAL,CELL_YQS_NAMES,ICQS_SPATIAL,
     '  IICQS_SPATIAL,IRCQS_SPATIAL,RCQS_SPATIAL,STRING,ERROR,*)

C#### Subroutine: EXFIELDML
C###  Description:
C###  EXFIELDML exports to FieldML. Initially just a test for exporting
C###  cellular grid fields as nodes.

C *** Created by DPN 07-03-2003 (based on EXPOIN)
C     Just trying out exporting cellular fields as node fields in
C     FieldML - mainly because the current methods of exporting these
C     would need updating to cope with CellML, so instead I'll begin
C     experimenting with FieldML.

      IMPLICIT none
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'grou00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER LIST(0:NLISTM),NQLIST(0:NQM),NXLIST(0:NXM),
     '  CELL_ICQS_SPATIAL(NQIM,NQVM),
     '  CELL_ICQS_VALUE(NQIM,NQVM),CELL_RCQS_SPATIAL(NQRM,NQVM),
     '  CELL_YQS_SPATIAL(NIQSM,NQVM),ICQS_SPATIAL(NQISVM,NQM),
     '  IICQS_SPATIAL(0:NQISVM,NQVM),IRCQS_SPATIAL(0:NQRSVM,NQVM)
      REAL*8 XQ(NJM,NQM),YQS(NIQSM,NQM),
     '  CELL_RCQS_VALUE(NQRM,NQVM),CELL_YQS_VALUE(NIQSM,NQVM),
     '  RCQS_SPATIAL(NQRSVM,NQM)
      CHARACTER ERROR*(*),
     '  CELL_ICQS_NAMES(NQIM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH),
     '  STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IBEG3,IBEG4,IEND,IEND1,IEND2,IEND3,IEND4,
     '  IFROMC,ILISTMBR,index,LOCAL_ICQS(NQIM),N3CO,nc,niqV,nj,NOLIST,
     &  nq,nr,NUM_COMPS,NUMFIELDS,NUM_FIELD_LIST,nx,nxc,OFFSET,
     &  VAL_INDEX,nqsv,j,k,variant
      REAL*8 FIELD_LIST(99),LOCAL_RCQS(NQRM)
      CHARACTER CHAR*30,CHAR1*5,CHAR2*5,CHAR3*5,COMPONENT_NAME*80,
     &  FIELD_NAME*80,GROUP_NAME*80,FILE*100,LABEL*30,OUTPUT*11,
     &  FIELDML_TYPE*50
      LOGICAL CBBREV

      CALL ENTERS('EXFIELDML',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM export fieldml<;FILENAME[default]> cellml
C###  Description:
C###    Export grid point cellular fields as FieldML node fields.
C###  Parameter:      <number LIST[all]>
C###    Specify the grid points to export.
C###  Parameter:      <as NAME[fields]>
C###    The regionML group name.
C###  Parameter:      <region #[1]>
C###    Specify the region number.
C###  Parameter:      <offset OFFSET[10000]>
C###    Add OFFSET to node numbers (names).
C###  Parameter:      <class>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']> cellml'
        OP_STRING(2)=BLANK(1:15)//'<number LIST[all]>'
        OP_STRING(3)=BLANK(1:15)//'<as NAME[fields]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<offset OFFSET[10000]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXFIELDML',ERROR,*9999)
      ELSE

        nc=1 !Temporary MPN 12-Nov-94

        IF(CBBREV(CO,'CELLML',2,noco+1,NTCO,N3CO)) THEN
          FIELDML_TYPE='CELLML'
        ELSE
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

C news MPN 9-11-95: OFFSET option to add const to all points
        IF(CBBREV(CO,'OFFSET',2,noco+1,NTCO,N3CO)) THEN
          OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          OFFSET=10000
        ENDIF

        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        IF(CBBREV(CO,'NUMBER' ,2,noco+1,NTCO,N3CO).OR
     '    .CBBREV(CO,'ELEMENT',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NQM,NQLIST(0),NQLIST(1),ERROR,*9999)
        ELSE
          IF(FIELDML_TYPE(1:6).EQ.'CELLML') THEN
            NQLIST(0)=NQT
            DO nq=1,NQT
              NQLIST(nq)=nq
            ENDDO
          ENDIF
        ENDIF
        CALL ASSERT(NQLIST(0).LE.NQM,'>> Increase NQM',ERROR,
     '    *9999)

        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          GROUP_NAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          GROUP_NAME='fields'
        ENDIF

        
        CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(NQLIST(0).GT.0) THEN
          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.fml','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)

          !write the XML header
          WRITE(IFILE,'(''<?xml version="1.0" '//
     '      'encoding="iso-8859-1"?>'')')
          WRITE(IFILE,'(''<!--Generated by CMISS(cm)-->'')')
          !The regionML wrapper to get a group name
          WRITE(IFILE,'(''<regionml>'')')
          CALL STRING_TRIM(GROUP_NAME,IBEG,IEND)
          WRITE(IFILE,'(''<group name="'//GROUP_NAME(IBEG:IEND)//'">''
     &      )')
          !open the fieldml element
          WRITE(IFILE,'(''<fieldml'')')
          WRITE(IFILE,
     '      '(''  xmlns="http://www.physiome.org.nz/fieldml/0.1#"'')')
          WRITE(IFILE,
     '      '(''  xmlns:fieldml='//
     '      '"http://www.physiome.org.nz/fieldml/0.1#">'')')
          
          IF(FIELDML_TYPE(1:6).EQ.'CELLML') THEN
            !Define the coordinate field
            WRITE(IFILE,'(''  <field '')')
            WRITE(IFILE,'(''    name="coordinates"'')')
            WRITE(IFILE,'(''    value_type="real"'')')
            WRITE(IFILE,'(''    coordinate_system='//
     '        '"rectangular cartesian">'')')
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              IF(nj.EQ.1) THEN
                WRITE(IFILE,'(''    <component name="x"/>'')')
              ELSEIF(nj.EQ.2) THEN
                WRITE(IFILE,'(''    <component name="y"/>'')')
              ELSEIF(nj.EQ.3) THEN
                WRITE(IFILE,'(''    <component name="z"/>'')')
              ENDIF
            ENDDO
            WRITE(IFILE,'(''  </field>'')')
            !Define the cell type (variant) field
            WRITE(IFILE,'(''  <field '')')
            WRITE(IFILE,'(''    name="cell_type"'')')
            WRITE(IFILE,'(''    value_type="integer"'')')
            WRITE(IFILE,'(''    coordinate_system='//
     '        '"rectangular cartesian">'')')
            WRITE(IFILE,'(''    <component name="type"/>'')')
            WRITE(IFILE,'(''  </field>'')')
C           Define the cellular fields for each variant - for now, we
C           don't care if some of the fields in different variants have
C           the same name as they will simply appear as the same field
C           in CMGUI...which is generally what we will want. Might want
C           to make the fields unique across variants one day.
            DO variant=1,CELL_NUM_VARIANTS
              NUM_FIELD_LIST=CELL_NUM_STATE(variant)
     '          +CELL_NUM_DERIVED(variant)
              CALL FIELDML_WR_FLD_HEADINGS_CELLML(IFILE,
     '          NUM_FIELD_LIST,CELL_YQS_NAMES(1,variant))
              NUM_FIELD_LIST=CELL_NUM_PARAMETERS(variant)
     '          +CELL_NUM_PROTOCOL(variant)
              CALL FIELDML_WR_FLD_HEADINGS_CELLML(IFILE,
     '          NUM_FIELD_LIST,CELL_RCQS_NAMES(1,variant))
            ENDDO !variant=1,CELL_NUM_VARIANTS
            !Define the coordinate field labels template
            WRITE(IFILE,'(''  <labels_template '//
     '        'name="coordinates_template">'')')
            WRITE(IFILE,'(''    <field_ref ref="coordinates">'')')
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              IF(nj.EQ.1) THEN
                WRITE(IFILE,'(''      <component_ref ref="x">'')')
                WRITE(IFILE,'(''        <label name="value"/>'')')
                WRITE(IFILE,'(''      </component_ref>'')')
              ELSEIF(nj.EQ.2) THEN
                WRITE(IFILE,'(''      <component_ref ref="y">'')')
                WRITE(IFILE,'(''        <label name="value"/>'')')
                WRITE(IFILE,'(''      </component_ref>'')')
              ELSEIF(nj.EQ.3) THEN
                WRITE(IFILE,'(''      <component_ref ref="z">'')')
                WRITE(IFILE,'(''        <label name="value"/>'')')
                WRITE(IFILE,'(''      </component_ref>'')')
              ENDIF
            ENDDO
            WRITE(IFILE,'(''    </field_ref>'')')
            WRITE(IFILE,'(''  </labels_template>'')')
            !Define the cell_type field labels template
            WRITE(IFILE,'(''  <labels_template '//
     '        'name="cell_type_template">'')')
            WRITE(IFILE,'(''    <field_ref ref="cell_type">'')')
            WRITE(IFILE,'(''      <component_ref ref="type">'')')
            WRITE(IFILE,'(''        <label name="value"/>'')')
            WRITE(IFILE,'(''      </component_ref>'')')
            WRITE(IFILE,'(''    </field_ref>'')')
            WRITE(IFILE,'(''  </labels_template>'')')
            !Define the cellular field template for each variant
            DO variant=1,CELL_NUM_VARIANTS
              WRITE(CHAR,'(I12)') variant
              CALL STRING_TRIM(CHAR,IBEG,IEND)
              WRITE(FIELD_NAME,'(''cellml_template_'//
     '          CHAR(IBEG:IEND)//''')')
              CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
              WRITE(IFILE,'(''  <labels_template '//
     '          'name="'//FIELD_NAME(IBEG:IEND)//'">'')')
              NUM_FIELD_LIST=CELL_NUM_STATE(variant)
     '          +CELL_NUM_DERIVED(variant)
              CALL FIELDML_WRITE_FIELD_REFS_CELLML(IFILE,
     '          NUM_FIELD_LIST,CELL_YQS_NAMES(1,variant))
              NUM_FIELD_LIST=CELL_NUM_PARAMETERS(variant)
     '          +CELL_NUM_PROTOCOL(variant)
              CALL FIELDML_WRITE_FIELD_REFS_CELLML(IFILE,
     '          NUM_FIELD_LIST,CELL_RCQS_NAMES(1,variant))
              WRITE(IFILE,'(''  </labels_template>'')')
            ENDDO !variant=1,CELL_NUM_VARIANTS
            !loop through the grid points writing out the nodes
            DO NOLIST=1,NQLIST(0)
              nq=NQLIST(NOLIST)
              WRITE(CHAR,'(I12)') nq+OFFSET
              CALL STRING_TRIM(CHAR,IBEG,IEND)
              WRITE(IFILE,'(''  <node name="'//CHAR(IBEG:IEND)//'">'')')
              !Define the coordinate field values
              WRITE(IFILE,'(''    <assign_labels template_name='//
     '          '"coordinates_template">'')')
              NUM_FIELD_LIST=0
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                NUM_FIELD_LIST=NUM_FIELD_LIST+1
                FIELD_LIST(NUM_FIELD_LIST)=XQ(nj,nq)
              ENDDO
              WRITE(IFILE,'(6X,3E13.5)') (FIELD_LIST(nj),nj=1,
     '          NUM_FIELD_LIST)
              WRITE(IFILE,'(''    </assign_labels>'')')
              !Define the cell_type field value
              WRITE(IFILE,'(''    <assign_labels template_name='//
     '          '"cell_type_template">'')')
              WRITE(IFILE,'(6X,I12)') ICQS_SPATIAL(1,nq)
              WRITE(IFILE,'(''    </assign_labels>'')')
              !Define the cellular field values
              variant = ICQS_SPATIAL(1,nq)
              WRITE(CHAR,'(I12)') variant
              CALL STRING_TRIM(CHAR,IBEG,IEND)
              WRITE(FIELD_NAME,'(''cellml_template_'//
     '          CHAR(IBEG:IEND)//''')')
              CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
              WRITE(IFILE,'(''    <assign_labels template_name="'
     '          //FIELD_NAME(IBEG:IEND)//'">'')')
              NUM_FIELD_LIST=CELL_NUM_STATE(variant)
     '          +CELL_NUM_DERIVED(variant)
              WRITE(IFILE,'(6X,6E13.5)') (YQS(j,nq),j=1,NUM_FIELD_LIST)
              CALL CELL_ASSIGN_SPATIAL(CELL_ICQS_VALUE,LOCAL_ICQS,
     '          ICQS_SPATIAL,IICQS_SPATIAL,IRCQS_SPATIAL,nq,
     '          CELL_RCQS_VALUE,LOCAL_RCQS,RCQS_SPATIAL,ERROR,*9999)
              NUM_FIELD_LIST=CELL_NUM_PARAMETERS(variant)
     '          +CELL_NUM_PROTOCOL(variant)
              WRITE(IFILE,'(6X,6E13.5)') (LOCAL_RCQS(j),j=1,
     '          NUM_FIELD_LIST)
              WRITE(IFILE,'(''    </assign_labels>'')')
              !and close the node element
              WRITE(IFILE,'(''  </node>'')')
            ENDDO
          ENDIF !IF(FIELDML_TYPE(1:6).EQ.'CELLML')

          !close the fieldml element
          WRITE(IFILE,'(''</fieldml>'')')
          !close the group element
          WRITE(IFILE,'(''</group>'')')
          !close the regionml element
          WRITE(IFILE,'(''</regionml>'')')
          
          CALL CLOSEF(IFILE,ERROR,*9999)
        ENDIF !IF(NQLIST(0).GT.0)

      ENDIF

      CALL EXITS('EXFIELDML')
      RETURN
 9999 CALL ERRORS('EXFIELDML',ERROR)
      CALL EXITS('EXFIELDML')
      RETURN 1
      END


