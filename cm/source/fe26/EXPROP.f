      SUBROUTINE EXPROP(NBJ,NEELEM,NELIST,NRLIST,NXLIST,CE,STRING,
     '  XAB,ERROR,*)

C#### Subroutine: EXPROP
C###  Description:
C###    EXPROP exports element based properties from the CE array.

C****   ELEM_NAME  is '1..#elements' or group name
C****   ELEM_TOTAL is total #elements in exported list

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NRLIST(0:NRM),NXLIST(0:NXM)
      REAL*8 CE(NMM,NEM,NXM),XAB(NORM,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER FIELD_BASE_TYPE,IBEG,
     '  IBEG1,IBEG2,IBEG3,IBEG4,IEND,IEND1,IEND2,IEND3,IEND4,IFROMC,
     '  N3CO,nb,ne,NFIELDT,nm,NM_LIST(0:30),NO_NELIST,no_nm,
     '  nr,nx,nxc,offset_elem
      CHARACTER CHAR1*5,CHAR2*5,CHAR3*5,CHAR4*20,ELEM_NAME*50,
     '  FIELD_NAME*80,FILE*200
      LOGICAL ALL_REGIONS,AUTONAME,CBBREV,NAME,ELEM_FIELD

      CALL ENTERS('EXPROP',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        WRITE(CHAR2,'(I5)') NEELEM(1,1)
        WRITE(CHAR3,'(I5)') NEELEM(NEELEM(0,1),1)
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
        CHAR4=CHAR2(IBEG2:IEND2)//'..'//CHAR3(IBEG3:IEND3)
        CALL STRING_TRIM(CHAR4,IBEG4,IEND4)

C---------------------------------------------------------------------

C#### Command: FEM export properties<;FILENAME[default]>
C###  Parameter:    <elements (GROUP/#s/all)[all]>
C###    Specify either element group, element numbers or all elements
C###    to be exported.
C###  Parameter:    <region (#s/all)[1]>
C###    Limit to elements of specified regions.
C###  Parameter:      <as NAME[0...0]>
C###    Name the element group with a character name.
C###  Parameter:      <name NAME[fieldname]>
C###    Override the default field name.
C###  Parameter:      <offset_elem OFFSET[0]>
C####   Add OFFSET to element numbers.
C###  Parameter:      <offset_node OFFSET[offset_elem]>
C###    Add OFFSET to node numbers
C###  Parameter:      <index #s[1]>
C###    For `field' specify the index in the YP, YQ or YQS dependent
C###    variable array for the desired aspect of the solution.

        OP_STRING(1)=STRING(1:IEND)
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<element (GROUP/#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<as NAME['//CHAR4(IBEG4:IEND4)//']>'
        OP_STRING(5)=BLANK(1:15)//'<offset_elem OFFSET[0]>'
        OP_STRING(6)=BLANK(1:15)//'<offset_node OFFSET[offset_elem]>'
        OP_STRING(7)=BLANK(1:15)//'<index #s[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXPROP',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
        CALL ASSERT(NELIST(0).GT.0,'>> No Elements to export ',
     '    ERROR,*9999)
C MHT 8-Sep-2004 adding classes
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'INDEX',3,noco+1,NTCO,N3CO)) THEN
          CDATA(1)='ELEMENTS' !not really, but don't use list anyway
          CALL PARSILG(NM_LIST,NMM,CDATA(1),CO(N3CO+1),ERROR,*9999)
        ELSE
          NM_LIST(0)=1
          NM_LIST(1)=1
        ENDIF
        IF(CBBREV(CO,'AUTONAME',4,noco+1,NTCO,N3CO)) THEN
          AUTONAME=.TRUE.
        ELSE
          AUTONAME=.FALSE.
        ENDIF
c        nx=1
        IF(CBBREV(CO,'OFFSET_ELEM',8,noco+1,NTCO,N3CO)) THEN
          offset_elem=IFROMC(CO(N3CO+1))
        ELSE
          offset_elem=0
        ENDIF
        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          ELEM_NAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          IF(AUTONAME) THEN
            WRITE(CHAR1,'(I5)') nr
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            ELEM_NAME='region_'//CHAR1(IBEG1:IEND1)
          ELSE !autoname
            WRITE(CHAR1,'(I5)') NELIST(1)+offset_elem
            WRITE(CHAR2,'(I5)') NELIST(NELIST(0))+offset_elem
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
            ELEM_NAME=CHAR1(IBEG1:IEND1)//'..'//CHAR2(IBEG2:IEND2)
          ENDIF !autoname
        ENDIF
        CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(CBBREV(CO,'NAME',1,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          FIELD_NAME=CO(N3CO+1)(IBEG1:IEND1)
          NAME=.TRUE.
        ELSE
          NAME=.FALSE.
        ENDIF
        IF(CBBREV(CO,'FIELD',4,noco+1,NTCO,N3CO)) THEN
          ELEM_FIELD=.TRUE.
        ELSE
          ELEM_FIELD=.FALSE.
        ENDIF
        IF(NELIST(0).GT.0) THEN
          NFIELDT=NM_LIST(0)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.exelem','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          CALL STRING_TRIM(ELEM_NAME,IBEG,IEND)
          WRITE(IFILE,'( '' Group name: '',A)') ELEM_NAME(IBEG:IEND)
          WRITE(IFILE,'( '' Shape.  Dimension='',I1)') 
     &      NIT(NBJ(1,NEELEM(1,nr)))
          WRITE(IFILE,'( '' #Scale factor sets=0'')')
          WRITE(IFILE,'( '' #Nodes=0'')')
          WRITE(IFILE,'( '' #Fields='',I1 )') NFIELDT
          DO no_nm=1,NM_LIST(0)
            nm=NM_LIST(no_nm)
            ne=NEELEM(1,nr) !first element in region
            nb=NBJ(1,ne)
            CALL EXPORT_FIELD_HEADING_PROP(FIELD_BASE_TYPE,IFILE,nb,
     '        NFIELDT,nm,NM_LIST,no_nm,FIELD_NAME,NAME,
     '        ERROR,*9999)
            IF(NIT(nb).EQ.1) THEN
              WRITE(IFILE,'(3X,''#xi1=1'')')
            ELSE IF(NIT(nb).EQ.2) THEN
              WRITE(IFILE,'(3X,''#xi1=1  #xi2=1'')')
            ELSE IF(NIT(nb).EQ.3) THEN
              WRITE(IFILE,'(3X,''#xi1=1  #xi2=1  #xi3=1'')')
            ENDIF
          ENDDO !no_nm(nm)
          DO no_nelist=1,NELIST(0)
            ne=NELIST(no_nelist)
            WRITE(IFILE,'(1X,''Element: '',I12,'' 0 0'' )')
     '        ne+offset_elem
            WRITE(IFILE,'(3X,''Values:'')')
            DO no_nm=1,NM_LIST(0)
              nm=NM_LIST(no_nm)
C              IF(ITYP3(nr,nx).EQ.3.AND.nm.EQ.17) !increases flow value
C     '          CE(nm,ne,nx)=CE(nm,ne,nx)*1000.d0 !for viewing
              IF(.NOT.ELEM_FIELD)THEN
                IF(NIT(NBJ(1,ne)).EQ.1) THEN !1d
                  WRITE(IFILE,'(3X,2(E18.8))') CE(nm,ne,nx),CE(nm,ne,nx)
                ELSE IF(NIT(NBJ(1,ne)).EQ.2) THEN !2d
                  WRITE(IFILE,'(3X,4(E18.8))') CE(nm,ne,nx),
     '              CE(nm,ne,nx),CE(nm,ne,nx),CE(nm,ne,nx)
                ELSE !3d
                  WRITE(IFILE,'(3X,4(E18.8))') CE(nm,ne,nx),
     '              CE(nm,ne,nx),CE(nm,ne,nx),CE(nm,ne,nx)
                  WRITE(IFILE,'(3X,4(E18.8))') CE(nm,ne,nx),
     '              CE(nm,ne,nx),CE(nm,ne,nx),CE(nm,ne,nx)
                ENDIF
              ELSE IF(ELEM_FIELD)THEN  !Write XAB elem field
                IF(NIT(NBJ(1,ne)).EQ.1) THEN !1d
                  WRITE(IFILE,'(3X,2(E18.8))') XAB(nm,ne),XAB(nm,ne)
                ELSE IF(NIT(NBJ(1,ne)).EQ.2) THEN !2d
                  WRITE(IFILE,'(3X,4(E18.8))') XAB(nm,ne),
     '              XAB(nm,ne),XAB(nm,ne),XAB(nm,ne)
                ELSE !3d
                  WRITE(IFILE,'(3X,4(E18.8))') XAB(nm,ne),
     '              XAB(nm,ne),XAB(nm,ne),XAB(nm,ne)
                  WRITE(IFILE,'(3X,4(E18.8))') XAB(nm,ne),
     '              XAB(nm,ne),XAB(nm,ne),XAB(nm,ne)
                ENDIF
              ENDIF
            ENDDO !no_nm(nm)
          ENDDO !no_nelist(ne)
        ENDIF
        CALL CLOSEF(IFILE,ERROR,*9999)
      ENDIF

      CALL EXITS('EXPROP')
      RETURN
 9999 CALL ERRORS('EXPROP',ERROR)
      CALL EXITS('EXPROP')
      RETURN 1
      END


