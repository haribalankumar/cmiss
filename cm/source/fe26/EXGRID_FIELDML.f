      SUBROUTINE EXGRID_FIELDML(NLATNE,NQNLAT,NQS,NQSCNB,NQXI,
     &  NRLIST,offset_elem,offset_node,XQ,YQ,
     &  ACTIVATION_TIME,ELASTIC_TUBE,PATH,POTENTIAL,PRESSURE,RADIUS,
     &  VELOCITY,ERROR,*)

C#### Subroutine: EXGRID_FIELDML
C###  Description:
C###    EXGRID_FIELDML creates an fieldml file containing lattice-based
C###    grid nodes and connectivity. Elements are exported as (2/4/8)
C###    noded elements with linear basis functions.
C***  Created by Greg Sands, Sept 2003

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER NLATNE(NEQM+1),NQNLAT(NEQM*NQEM),NQS(NEQM),
     &  NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),
     &  offset_elem,offset_node
      REAL*8 XQ(NJM,NQM),YQ(NYQM,NIQM,NAM,NXM)
      LOGICAL ACTIVATION_TIME,ELASTIC_TUBE,PATH,POTENTIAL,PRESSURE,
     &  RADIUS,VELOCITY
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER ELEMENT,IBEG,IEND,nb,ne,
     &  NITB,nj,NJQ,nlat,no_nrlist,nq,nqq,nr,nx,SCHEME

      CHARACTER BasisName*8,CHAR1*10,CoordName*14,ElementName*10,
     &  ElNodeName*15,ElNodeValName*14,InterpName*23,NodeName*8

      CALL ENTERS('EXGRID_FIELDML',*9999)

      DO no_nrlist=1,NRLIST(0)
        nr=NRLIST(no_nrlist)
        SCHEME=NQS(1)
        nb=NQSCNB(SCHEME)
        NITB=NQXI(0,SCHEME)
        NJQ=NJ_LOC(NJL_GEOM,0,nr)

        WRITE(BasisName,'(''"Basis'',I1,''"'')') nr
        WRITE(CoordName,'(''"Coordinates'',I1,''"'')') nr
        WRITE(ElementName,'(''"Element'',I1,''"'')') nr
        WRITE(ElNodeName,'(''"ElementNodes'',I1,''"'')') nr
        WRITE(ElNodeValName,'(''"NodalValues'',I1,''"'')') nr
        WRITE(InterpName,'(''"ElementInterpolation'',I1,''"'')') nr
        WRITE(NodeName,'(''"Nodes'',I1,''"'')') nr

C***    Write out field definition
        WRITE(IFILE,'(/,2X,''<field name='',A,'' value_type="real"'
     &    //' coordinate_system="rectangular cartesian">'')')
     &    CoordName
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          IF(nj.EQ.1) WRITE(IFILE,'(4X,''<component name="x" />'')')
          IF(nj.EQ.2) WRITE(IFILE,'(4X,''<component name="y" />'')')
          IF(nj.EQ.3) WRITE(IFILE,'(4X,''<component name="z" />'')')
        ENDDO
        WRITE(IFILE,'(2X,''</field>'')')

C***    Write out node template
        WRITE(IFILE,'(/,2X,''<labels_template name='',A,''>'')')
     &    NodeName
        WRITE(IFILE,'(4X,''<field_ref ref='',A,''>'')')
     &    CoordName
        DO nj=1,NJQ
          IF(nj.EQ.1) WRITE(IFILE,'(6X,''<component_ref ref="x">'')')
          IF(nj.EQ.2) WRITE(IFILE,'(6X,''<component_ref ref="y">'')')
          IF(nj.EQ.3) WRITE(IFILE,'(6X,''<component_ref ref="z">'')')
          WRITE(IFILE,'(8X,''<label name="value" />'')')
          WRITE(IFILE,'(6X,''</component_ref>'')')
        ENDDO
        WRITE(IFILE,'(4X,''</field_ref>'')')
        WRITE(IFILE,'(2X,''</labels_template>'')')

C***    Write out basis mapping
        WRITE(IFILE,'(/,2X,''<mapping name='',A,'
     &    //''' basis="l.Lagrange*l.Lagrange*l.Lagrange">'')')
     &    BasisName
        WRITE(IFILE,'(4X,''<coefficients>'')')
        WRITE(IFILE,'(6X,''<element_lookup>'')')
        WRITE(IFILE,'(8X,''<field_lookup>'')')
        WRITE(IFILE,'(10X,''<component_lookup>'')')
        WRITE(IFILE,'(12X,''<label_lookup indices='',A,'' />'')')
     &    ElNodeValName
        WRITE(IFILE,'(10X,''</component_lookup>'')')
        WRITE(IFILE,'(8X,''</field_lookup>'')')
        WRITE(IFILE,'(6X,''</element_lookup>'')')
        WRITE(IFILE,'(4X,''</coefficients>'')')
        WRITE(IFILE,'(2X,''</mapping>'')')

C***    Write out element template
        WRITE(IFILE,'(/,2X,''<labels_template name='',A,''>'')')
     &    ElementName
        DO nj=1,NJQ
          WRITE(IFILE,'(4X,''<node_lookup>'')')
          WRITE(IFILE,'(6X,''<node_index>'')')
          WRITE(IFILE,'(8X,''<element_lookup>'')')
          WRITE(IFILE,'(10X,''<label_lookup indices='',A,''>'')')
     &      ElNodeName
          WRITE(IFILE,'(12X,''<label_lookup indices="'',I1,''" />'')')nj
          WRITE(IFILE,'(10X,''</label_lookup>'')')
          WRITE(IFILE,'(8X,''</element_lookup>'')')
          WRITE(IFILE,'(6X,''</node_index>'')')
          WRITE(IFILE,'(6X,''<field_lookup>'')')
          WRITE(IFILE,'(8X,''<component_lookup>'')')
          WRITE(IFILE,'(10X,''<label_lookup indices="value" />'')')
          WRITE(IFILE,'(8X,''</component_lookup>'')')
          WRITE(IFILE,'(6X,''</field_lookup>'')')
          WRITE(IFILE,'(4X,''</node_lookup>'')')
        ENDDO
        WRITE(IFILE,'(2X,''</labels_template>'')')

C***    Write out element interpolation
        WRITE(IFILE,'(/,2X,''<element_interpolation name='',A,''>'')')
     &    InterpName
        WRITE(IFILE,'(4X,''<field_ref ref='',A,''>'')')
     &    CoordName
        DO nj=1,NJQ
          IF(nj.EQ.1) WRITE(IFILE,'(6X,''<component_ref ref="x">'')')
          IF(nj.EQ.2) WRITE(IFILE,'(6X,''<component_ref ref="y">'')')
          IF(nj.EQ.3) WRITE(IFILE,'(6X,''<component_ref ref="z">'')')
          WRITE(IFILE,'(8X,''<mapping_ref ref='',A,'' />'')')
     &      BasisName
          WRITE(IFILE,'(8X,''<label name='',A,''>'')')
     &      ElNodeValName
          WRITE(IFILE,'(10X,''<labels_template_ref ref='',A,'' />'')')
     &      ElementName
          WRITE(IFILE,'(8X,''</label>'')')
          WRITE(IFILE,'(6X,''</component_ref>'')')
        ENDDO
        WRITE(IFILE,'(4X,''</field_ref>'')')
        WRITE(IFILE,'(2X,''</element_interpolation>'')')

C***    Write out nodes
        DO nq=NQR(1,nr),NQR(2,nr)
          WRITE(CHAR1,'(I8)') nq+offset_node
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
          WRITE(IFILE,'(/,2X,''<node name="'',A,''">'')')
     &      CHAR1(IBEG:IEND)
          WRITE(IFILE,'(4X,''<assign_labels template_name='',A,''>'')')
     &      NodeName
          WRITE(IFILE,'(6X,5F24.12)') (XQ(nj,nq),nj=1,NJQ)
          WRITE(IFILE,'(4X,''</assign_labels>'')')
          WRITE(IFILE,'(2X,''</node>'')')
        ENDDO

C***    Write out elements
        DO ne=1,NEQM
          WRITE(CHAR1,'(I8)') ne+offset_elem
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
          WRITE(IFILE,'(/,2X,''<element name="'',A,''" shape='
     &      //'"line*line*line">'')') CHAR1(IBEG:IEND)
          WRITE(IFILE,'(4X,''<element_interpolation_ref ref='',A,'
     &      //''' />'')')
     &      InterpName
          WRITE(IFILE,'(4X,''<label name='',A,''>'')')
     &      ElNodeName
          WRITE(IFILE,'(6X,8I10)')
     &       (NQNLAT(nlat)+offset_node,nlat=NLATNE(ne),NLATNE(ne+1)-1)
          WRITE(IFILE,'(4X,''</label>'')')
          WRITE(IFILE,'(2X,''</element>'')')
        ENDDO

      ENDDO

      CALL EXITS('EXGRID_FIELDML')
      RETURN
 9999 CALL ERRORS('EXGRID_FIELDML',ERROR)
      CALL EXITS('EXGRID_FIELDML')
      RETURN 1
      END
