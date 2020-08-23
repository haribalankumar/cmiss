      SUBROUTINE EXPORT_FIELD_HEADING_PROP(FIELD_BASE_TYPE,IFILE,nb,
     '  NFIELDT,nhx,NM_LIST,no_nhx,FIELD_NAME,NAME,ERROR,*)

C#### Subroutine: EXPORT_FIELD_HEADING_PROP
C###  Description:
C###    EXPORT_FIELD_HEADING_PROP writes the field heading when
C###    exporting element based properties from CE (EXELEM).

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER FIELD_BASE_TYPE,IFILE,nb,NFIELDT,nhx,
     '  NM_LIST(0:30),no_nhx
      LOGICAL NAME
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG1,IBEG2,IBEG3,IEND1,IEND2,IEND3,IFIELD,NCOMPT,
     '  NCOORD,NFIELD
      CHARACTER BASES*54,BASE_TYPE*50,COMPONENT_NAME*80,
     '  COORDINATE_SYSTEM(6)*21,FIELD_NAME*80,FIELD_TYPE(3)*10,
     '  PROPERTY_11_03(21)*13
      DATA COORDINATE_SYSTEM/
     '  'rectangular cartesian','cylindrical polar    ',
     '  'spherical polar      ','prolate spheroidal   ',
     '  'oblate spheroidal    ','fibre                '/
      DATA FIELD_TYPE/
     '  'coordinate','anatomical','field     '/
      DATA PROPERTY_11_03/
     '  'length1      ','length2      ','length3      ',
     '  'xsec_area    ','hyd_diameter ','bvolume1     ',
     '  'bvolume2     ','vol_ratio    ','RBC_flow     ',
     '  'app_visc     ','hematocrit   ','major_axis   ',
     '  'minor_axis   ','Rfactor      ','resistance   ',
     '  'hematcrit0   ','flow         ','transit_t    ',
     '  'length_a0    ','length_C0    ','probability_Q'/

      CALL ENTERS('EXPORT_FIELD_HEADING_PROP',*9999)

      FIELD_BASE_TYPE=1
      NFIELDT=NM_LIST(0)
      IF(.NOT.NAME)THEN
c      IF(ITYP2(nr,nx).EQ.11)THEN !pulmonary transport
        CALL STRING_TRIM(PROPERTY_11_03(nhx),IBEG1,IEND1)
        FIELD_NAME=PROPERTY_11_03(nhx)(IBEG1:IEND1)
        COMPONENT_NAME=PROPERTY_11_03(nhx)(IBEG1:IEND1)
      ELSE
        COMPONENT_NAME=FIELD_NAME
      ENDIF
        NFIELD=no_nhx !for list of fields in .exelem file
        NCOMPT=1 !# of components in field, don't change
        IFIELD=3 !type is field, this mustn't change
        NCOORD=1 !just assuming rectangular cartesian
c      ENDIF

C**   write the number of fields
      CALL STRING_TRIM(FIELD_NAME,IBEG1,IEND1)
      CALL STRING_TRIM(FIELD_TYPE(IFIELD),IBEG2,IEND2)
      CALL STRING_TRIM(COORDINATE_SYSTEM(NCOORD),IBEG3,IEND3)
      WRITE(IFILE,'( I2,'') '',A,'', '',A,'', '',A,'
     '  //''', #Components='',I1)') NFIELD, FIELD_NAME(IBEG1:IEND1),
     '  FIELD_TYPE(IFIELD)(IBEG2:IEND2),
     '  COORDINATE_SYSTEM(NCOORD)(IBEG3:IEND3),NCOMPT
C**   write element field component heading
      CALL STRING_TRIM(COMPONENT_NAME,IBEG1,IEND1)
      BASE_TYPE='grid based'
      IF(nb.NE.0) THEN
        IF(NIT(nb).EQ.1) THEN
          BASES='l.Lagrange'
        ELSE IF(NIT(nb).EQ.2) THEN
          BASES='l.Lagrange*l.Lagrange'
        ELSE IF(NIT(nb).EQ.3) THEN
          BASES='l.Lagrange*l.Lagrange*l.Lagrange'
        ENDIF
      ELSE
        ERROR='>>Invalid basis number'
        GOTO 9999
      ENDIF
      CALL STRING_TRIM(BASES,IBEG2,IEND2)
      CALL STRING_TRIM(BASE_TYPE,IBEG3,IEND3)
      WRITE(IFILE,'(3X,A,''.  '',A,'', no modify, '','
     '  //'A,''.'')') COMPONENT_NAME(IBEG1:IEND1),
     '  BASES(IBEG2:IEND2),BASE_TYPE(IBEG3:IEND3)

      CALL EXITS('EXPORT_FIELD_HEADING_PROP')
      RETURN
 9999 CALL ERRORS('EXPORT_FIELD_HEADING_PROP',ERROR)
      CALL EXITS('EXPORT_FIELD_HEADING_PROP')
      RETURN 1
      END



