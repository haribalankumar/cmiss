      SUBROUTINE EXPORT_FIELD_HEADING_CELL(IFILE,nb,ne,
     '  NFIELDT,NQS,NQXI,
     '  CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,CELL_ICQS_NAMES,
     '  CELL_RCQS_VALUE,CELL_RCQS_SPATIAL,CELL_RCQS_NAMES,
     '  CELL_YQS_VALUE,CELL_YQS_SPATIAL,CELL_YQS_NAMES,
     '  ICQS_SPATIAL,IICQS_SPATIAL,IRCQS_SPATIAL,RCQS_SPATIAL,
     '  ERROR,*)

C#### Subroutine: EXPORT_FIELD_HEADING_CELL
C###  Description:
C###    Writes out the field headings when exporting cellular
C###    material parameters.
C *** Created - DPN 28 September 1999

C *** DPN 28 September 1999 - to start with, only worry about
C ***   spatially varying REAL parameters.

C *** DPN 16 November 1999 - now worry about exporting all parameters
C ***   and variables

      IMPLICIT NONE
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

      !Parameters
      INTEGER IFILE,nb,ne,NFIELDT,NQS(NEQM),NQXI(0:NIM,NQSCM),
     '  CELL_ICQS_SPATIAL(NQIM,NQVM),CELL_ICQS_VALUE(NQIM,NQVM),
     '  CELL_RCQS_SPATIAL(NQRM,NQVM),CELL_YQS_SPATIAL(NIQSM,NQVM),
     '  ICQS_SPATIAL(NQISVM,NQM),
     '  IICQS_SPATIAL(0:NQISVM,NQVM),
     '  IRCQS_SPATIAL(0:NQRSVM,NQVM)
      REAL*8 CELL_RCQS_VALUE(NQRM,NQVM),CELL_YQS_VALUE(NIQSM,NQVM),
     '  RCQS_SPATIAL(NQRSVM,NQM)
      CHARACTER ERROR*(*),
     '  CELL_ICQS_NAMES(NQIM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH)
      !Local variables
      INTEGER IBEG,IEND,IBEG1,IEND1,IBEG2,IEND2,IBEG3,IEND3,field,
     '  nqv,nqsc,k
      CHARACTER CHAR*30,BASES*54,XI*54
      LOGICAL CELL_FIELD

      CALL ENTERS('EXPORT_FIELD_HEADING_CELL',*9999)

C *** Set the number of fields
C       (model ID + cell type + # variables + # integer parameters
C       less the cell variant + # real parameters)
      NFIELDT=1+1+CELL_NUM_STATE(1)+(NQIT-1)+NQRT

C KAT 23/3/00: moving to EXELEM as this is not part of the field heading
CC *** NJH_FIELD_BASE_LIST(0) is the number of fields
C      NJH_FIELD_BASE_LIST(0)=NFIELDT

C *** Write out the number of fields
      WRITE(IFILE,'( '' #Fields='',I12)') NFIELDT

C *** Set the bases string for use with the fields and the xi string
      nqsc=NQS(ne)
      IF(nb.NE.0) THEN
        IF(NIT(nb).EQ.1) THEN
          BASES='l.Lagrange'
          WRITE(XI,'(3X,''#xi1='',I5)') NQXI(1,nqsc)-1
        ELSE IF(NIT(nb).EQ.2) THEN
          BASES='l.Lagrange*l.Lagrange'
          WRITE(XI,'(3X,''#xi1='',I5,'
     '      //''', #xi2='',I5)') NQXI(1,nqsc)-1,
     '      NQXI(2,nqsc)-1
        ELSE IF(NIT(nb).EQ.3) THEN
          BASES='l.Lagrange*l.Lagrange*l.Lagrange'
          WRITE(XI,'(3X,''#xi1='',I5,'
     '      //''', #xi2='',I5,'', #xi3='',I5)')
     '      NQXI(1,nqsc)-1,NQXI(2,nqsc)-1,
     '      NQXI(3,nqsc)-1
        ENDIF
      ELSE
        ERROR='>>Invalid basis number'
        GOTO 9999
      ENDIF
      CALL STRING_TRIM(BASES,IBEG2,IEND2)
      CALL STRING_TRIM(XI,IBEG3,IEND3)

C *** Loop through the fields
      DO field=1,NFIELDT

C KAT 23/3/00: moving to EXELEM as this is not part of the field heading
C        !used to loop through the variables to be exported
C        NJH_FIELD_BASE_LIST(field)=field

        IF(field.EQ.1) THEN
C ***     Write out the header for the model ID
          WRITE(IFILE,'(1X,I1,'') model_id, field, constant, '
     '      //'string, #Components=1'')') 1
          WRITE(IFILE,'(1X,''  ID.'')')
        ELSEIF(field.EQ.2) THEN
C ***     Write out the header for the cell type field
          WRITE(IFILE,'(1X,I1,'') cell_type, field, integer, '
     '      //'#Components=1'')') 2
          WRITE(IFILE,'(3X,''variant.  '',A,'', no modify, grid '
     '      //'based.'')') BASES(IBEG2:IEND2)
          WRITE(IFILE,'(3X,A)') XI(IBEG3:IEND3)
        ELSEIF(field.LT.NQIT+2) THEN
C ***     Write out the headers for the integer parameter fields
          nqv=field-1 !skip variant number
          CALL STRING_TRIM(CELL_ICQS_NAMES(nqv,1),IBEG,IEND)
          WRITE(CHAR,'(I12)') field
          CALL STRING_TRIM(CHAR,IBEG1,IEND1)
          CELL_FIELD=.FALSE.
          DO k=1,IICQS_SPATIAL(0,1)
            IF(IICQS_SPATIAL(k,1).EQ.nqv) THEN
              CELL_FIELD=.TRUE.
            ENDIF
          ENDDO !k
          IF(CELL_FIELD) THEN
C ***       normal element based field
            WRITE(IFILE,'(1X,A,'') '',A,'', field, integer, '
     '        //'#Components=1'')') CHAR(IBEG1:IEND1),
     '        CELL_ICQS_NAMES(nqv,1)(IBEG:IEND)
            WRITE(IFILE,'(3X,''value.  '',A,'', no modify, grid '
     '        //'based.'')') BASES(IBEG2:IEND2)
            WRITE(IFILE,'(3X,A)') XI(IBEG3:IEND3)
          ELSE
C ***       indexed element field
            WRITE(IFILE,'(1X,A,'') '',A,'', field, indexed, '
     '        //'Index_field=cell_type, #Values='',I12,'', '
     '        //'integer, #Components=1'')')
     '        CHAR(IBEG1:IEND1),
     '        CELL_ICQS_NAMES(nqv,1)(IBEG:IEND),
     '        CELL_NUM_VARIANTS
            WRITE(IFILE,'(1X,''  value.'')')
          ENDIF
        ELSEIF(field.LT.NQRT+NQIT+2) THEN
C ***     Write out the headers for the real parameter fields
          nqv=field-1-NQIT
          CALL STRING_TRIM(CELL_RCQS_NAMES(nqv,1),IBEG,IEND)
          WRITE(CHAR,'(I12)') field
          CALL STRING_TRIM(CHAR,IBEG1,IEND1)
          CELL_FIELD=.FALSE.
          DO k=1,IRCQS_SPATIAL(0,1)
            IF(IRCQS_SPATIAL(k,1).EQ.nqv) THEN
              CELL_FIELD=.TRUE.
            ENDIF
          ENDDO !k
          IF(CELL_FIELD) THEN
C ***       normal element based field
            WRITE(IFILE,'(1X,A,'') '',A,'', field, real, '
     '        //'#Components=1'')') CHAR(IBEG1:IEND1),
     '        CELL_RCQS_NAMES(nqv,1)(IBEG:IEND)
            WRITE(IFILE,'(3X,''value.  '',A,'', no modify, grid '
     '        //'based.'')') BASES(IBEG2:IEND2)
            WRITE(IFILE,'(3X,A)') XI(IBEG3:IEND3)
          ELSE
C ***       indexed element field
            WRITE(IFILE,'(1X,A,'') '',A,'', field, indexed, '
     '        //'Index_field=cell_type, #Values='',I12,'', '
     '        //'real, #Components=1'')')
     '        CHAR(IBEG1:IEND1),
     '        CELL_RCQS_NAMES(nqv,1)(IBEG:IEND),
     '        CELL_NUM_VARIANTS
            WRITE(IFILE,'(1X,''  value.'')')
          ENDIF
        ELSE
C ***     Write out the headers for the state variables - always normal
C ***     element based fields
          nqv=field-1-NQIT-NQRT
          CALL STRING_TRIM(CELL_YQS_NAMES(CELL_STATE_OFFSET(1)-1+nqv,1),
     '      IBEG,IEND)
          WRITE(CHAR,'(I12)') field
          CALL STRING_TRIM(CHAR,IBEG1,IEND1)
          WRITE(IFILE,'(1X,A,'') '',A,'', field, real, '
     '      //'#Components=1'')') CHAR(IBEG1:IEND1),
     '      CELL_YQS_NAMES(CELL_STATE_OFFSET(1)-1+nqv,1)(IBEG:IEND)
          WRITE(IFILE,'(3X,''value.  '',A,'', no modify, grid '
     '      //'based.'')') BASES(IBEG2:IEND2)
          WRITE(IFILE,'(3X,A)') XI(IBEG3:IEND3)
        ENDIF

      ENDDO !nqv

C *** Now write out the values for the indexed fields
      WRITE(IFILE,'(1X,''Values:'')')
C *** The model ID
      WRITE(IFILE,'(1X,''MODEL_ID'')')
C *** The integer parameters
      DO nqv=2,NQIT !skip variant number
        CELL_FIELD=.FALSE.
        DO k=1,IICQS_SPATIAL(0,1)
          IF(IICQS_SPATIAL(k,1).EQ.nqv) THEN
            CELL_FIELD=.TRUE.
          ENDIF
        ENDDO !k
        IF(CELL_FIELD) THEN
          !normal node field
          ! ** do nothing **
        ELSE
          !indexed field - write out the variant values
          WRITE(IFILE,'(I12)')
     '      (CELL_ICQS_VALUE(nqv,k),k=1,CELL_NUM_VARIANTS)
        ENDIF
      ENDDO !nqv
C *** The real parameters
      DO nqv=1,NQRT
        CELL_FIELD=.FALSE.
        DO k=1,IRCQS_SPATIAL(0,1)
          IF(IRCQS_SPATIAL(k,1).EQ.nqv) THEN
            CELL_FIELD=.TRUE.
          ENDIF
        ENDDO !k
        IF(CELL_FIELD) THEN
          !normal node field
          ! ** do nothing **
        ELSE
          !indexed field
          WRITE(IFILE,'(9E13.5)')
     '      (CELL_RCQS_VALUE(nqv,k),k=1,CELL_NUM_VARIANTS)
        ENDIF
      ENDDO !nqv

      CALL EXITS('EXPORT_FIELD_HEADING_CELL')
      RETURN

 9999 CALL ERRORS('EXPORT_FIELD_HEADING_CELL',ERROR)
      CALL EXITS('EXPORT_FIELD_HEADING_CELL')
      RETURN 1

      END


