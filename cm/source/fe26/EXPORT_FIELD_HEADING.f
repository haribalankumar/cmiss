      SUBROUTINE EXPORT_FIELD_HEADING(ELEMENT_DIMENSION,
     '  FIELD_BASE_TYPE,IBT,IFILE,IVALUE,iy,nb,NBH,
     '  nc,NDERIV,ne,NFIELDT,nhx,NHE,NHP,nj,
     '  NJLIST,np,nr,NVERSIONS,NW,nx,SPECIAL_BASIS_FLAG,
     '  AUTONAME,DATAFILE,DEFFIBS,FIELD_EX_TYPE,INPUT_FIELD_NAME,
     '  SET_FIELD_NAME,GRID_NUMBERS,ERROR,*)

C#### Subroutine: EXPORT_FIELD_HEADING
C###  Description:
C###    EXPORT_FIELD_HEADING writes the field heading when exporting
C###    nodes (EXNODE) or elements (EXELEM).
C     SPECIAL_BASIS_FLAG is zero except for the Hermite simplicies or
C     sector elements, where
C       1 is Hermite simplex with apex at node 1
C       2 is Hermite simplex with apex at node 3
C       3 is sector element
C       4 other bases with cross derivs set to zero

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
      INTEGER ELEMENT_DIMENSION,FIELD_BASE_TYPE,IBT(3,NIM,NBFM),
     '  IFILE,IVALUE,iy,
     '  nb,NBH(NHM,NCM,NEM),
     '  nc,NDERIV,ne,NFIELDT,
     '  nhx,NHE(NEM),NHP(NPM),nj,
     '  np,nr,
     '  NVERSIONS,
     '  NW(NEM,3),nx,SPECIAL_BASIS_FLAG
      CHARACTER INPUT_FIELD_NAME*50,FIELD_EX_TYPE*8,ERROR*(*)
      LOGICAL AUTONAME,DATAFILE,DEFFIBS,GRID_NUMBERS,
     '  NJLIST(3),SET_FIELD_NAME
!     Local Variables
      INTEGER CLEN,
     '  IB2,IBEG1,IBEG2,IBEG3,IBEG5,IE2,IEND1,IEND2,IEND3,
     '  IEND5,IFIELD,INTSTR(1024),ISOCKET,
     '  nb_aux,NCOMP,NCOMPT,NCOORD,
     '  NFIELD,nh_aux,nhx2,nk
      CHARACTER BASES*54,BASE_TYPE*50,CHAR2*2,COMPONENT_NAME*80,
     '  COORDINATE_COMPONENTS(3,6)*17,COORDINATE_SYSTEM(6)*21,
     '  DERIVATIVES(7)*15,DERIVATIVE_STRING*100,
     '  FIELD_NAME*80,FIELD_TYPE(3)*10,ISIZE*1,
     '  TEMP_FIELD_NAME*80
      LOGICAL FIBRES_DEFINED,
     '  FIELDS_DEFINED,WRITE_HEADING
      DATA COORDINATE_SYSTEM/
     '  'rectangular cartesian',
     '  'cylindrical polar    ',
     '  'spherical polar      ',
     '  'prolate spheroidal   ',
     '  'oblate spheroidal    ',
     '  'fibre                '/
      DATA COORDINATE_COMPONENTS/
     '  'x                ','y                ','z                ',
     '  'r                ','theta            ','z                ',
     '  'r                ','theta            ','phi              ',
     '  'lambda           ','mu               ','theta            ',
     '  'lambda           ','mu               ','theta            ',
     '  'fibre angle      ','imbrication angle','sheet angle      '/
      DATA FIELD_TYPE/
     '  'coordinate',
     '  'anatomical',
     '  'field     '/
      DATA DERIVATIVES/
     '  'd/ds1,','d/ds2,','d2/ds1ds2,','d/ds3,','d2/ds1ds3,',
     '  'd2/ds2ds3,','d3/ds1ds2ds3,'/

      CALL ENTERS('EXPORT_FIELD_HEADING',*9999)

C LKC 15-NOV-97 Changing the way the total number of field is
C calculated as it is wrong for nodes which are located at region
C interfaces and also have different field components in
C different regions
C
C      FIBRES_DEFINED=NJ_LOC(NJL_FIBR,0,nr).GT.0
C      FIELDS_DEFINED=NJ_LOC(NJL_FIEL,0,nr).GT.0
C LKC 20-DEC-97 new way of determining if fibres/fields have been
C     defined
      FIBRES_DEFINED=NJLIST(2)
      FIELDS_DEFINED=NJLIST(3)

      FIELD_BASE_TYPE=1
      IF(nj.GT.0) THEN
C**     non-solution field (geometric, anatomical, fitted, ...)
        NFIELDT=1
        IF(FIBRES_DEFINED) NFIELDT=NFIELDT+1
        IF(FIELDS_DEFINED) NFIELDT=NFIELDT+1
C DB 30-JUL-99 Exports grid point number with geometry
        IF((NE.GT.0).AND.(NQT.GT.0).AND.GRID_NUMBERS) NFIELDT=NFIELDT+1

        IF(NJL_GEOM.EQ.NJ_TYPE(nj,1)) THEN
          NFIELD=1
          NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
          FIELD_NAME='coordinates'
          IFIELD=1
          NCOORD=ITYP10(nr)
          NCOMP=NJ_TYPE(nj,2)
          COMPONENT_NAME=COORDINATE_COMPONENTS(NCOMP,NCOORD)
        ELSE IF(NJL_FIBR.EQ.NJ_TYPE(nj,1)) THEN
          NFIELD=2
C         Avoid writing out indices if not required
          IF(NFIELD.GT.NFIELDT) SPECIAL_BASIS_FLAG=-1
          NCOMPT=NJ_LOC(NJL_FIBR,0,nr)
          FIELD_NAME='fibres'
          IFIELD=2
          NCOORD=6
          NCOMP=NJ_TYPE(nj,2)
          COMPONENT_NAME=COORDINATE_COMPONENTS(NCOMP,NCOORD)
        ELSE IF(NJL_FIEL.EQ.NJ_TYPE(nj,1)) THEN
          IF(FIBRES_DEFINED) THEN
            NFIELD=3
          ELSE
            NFIELD=2
          ENDIF
C         Avoid writing out indices if not required
          IF(NFIELD.GT.NFIELDT) SPECIAL_BASIS_FLAG=-1
          NCOMPT=NJ_LOC(NJL_FIEL,0,nr)
          IFIELD=3
          NCOMP=NJ_TYPE(nj,2)

          IF ((ITYP3(nr,nx).EQ.1).AND.(ITYP2(nr,nx).EQ.5)
     '      .AND.(ITYP4(nr,nx).EQ.3)) THEN
C flow in elastic tube
            FIELD_NAME='radius'
            NCOORD=ITYP10(nr)
            WRITE(COMPONENT_NAME,'(I2)') NCOMP

C!!! CS 8/4/2003 time integration & cellular based modelling with
C!!! fields defined does not necessarily mean the fields are radii!
C          ELSEIF ((ITYP5(nr,nx).EQ.2).AND.(ITYP2(nr,nx).EQ.9).AND.
C     '      CALL_FIEL.AND.CALL_ELFB) THEN
CC purkinje fibre tree with radius defined
C            FIELD_NAME='radius'
C            NCOORD=ITYP10(nr)
C            COMPONENT_NAME='r'

          ELSE
            FIELD_NAME='general'
            NCOORD=1
            WRITE(COMPONENT_NAME,'(I2)') NCOMP
          ENDIF
        ELSE
          ERROR='>>Invalid NJ number'
          GOTO 9999
        ENDIF

      ELSE !solution field

        IF(NX_TYPE(nx).EQ.NX_SOLVE) THEN
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
          IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '       ITYP4(nr,nx).EQ.7) THEN !Collocation
            FIELD_BASE_TYPE=2 !Grid Based
C KAT 23/3/00: moving to EXELEM as this is not part of the field heading
C            NJH_FIELD_BASE_LIST(0)=NJH_FIELD_BASE_LIST(0)+1
C            NJH_FIELD_BASE_LIST(NJH_FIELD_BASE_LIST(0))=nhx
          ENDIF
          IF(ITYP5(nr,nx).EQ.1) THEN !Static analysis
            IF(ITYP2(nr,nx).EQ.1) THEN !Linear elasticity
C GMH 12/12/95 Do everything the same as finite elasticity
C              Export as a coordinate field of 3 components
C              rather than 3 standard fields
              IF(np.GT.0) THEN
                NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                NFIELDT=NHP(np)-NCOMPT+1
              ELSE
                NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                NFIELDT=NHE(ne)-NCOMPT+1
              ENDIF
              IF(nhx.LE.NCOMPT) THEN
                NFIELD=1
                NCOMP=nhx
                FIELD_NAME='deformed'
                IFIELD=1
                NCOORD=ITYP10(nr)
                COMPONENT_NAME=COORDINATE_COMPONENTS(NCOMP,NCOORD)
              ELSE
                NFIELD=nhx-NCOMPT+1
                NCOMPT=1
                NCOMP=1
                WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME=FIELD_NAME
              ENDIF
            ELSE IF(ITYP2(nr,nx).EQ.2.OR.
     '        ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
C             Finite elasticity or fluid mechanics/const vol constraint
              IF(np.GT.0) THEN
                NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                NFIELDT=NHP(np)-NCOMPT+1
              ELSE
                NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                NFIELDT=NHE(ne)-NCOMPT+1
C!!! DB.  Temporary. Don't export auxiliary only basis
                DO nhx2=NJ_LOC(NJL_GEOM,0,nr)+1,NHE(ne)
                  nh_aux=NH_LOC(nhx2,nx)
                  nb_aux=NBH(nh_aux,1,ne)
                  IF(NBC(nb_aux).EQ.0) NFIELDT=NFIELDT-1
                ENDDO
              ENDIF
              IF(FIBRES_DEFINED) NFIELDT=NFIELDT+1
              IF(DEFFIBS) THEN !deformed fibre angle fields
                NFIELD=NFIELDT
                NCOMPT=NJ_LOC(NJL_FIBR,0,nr)
                FIELD_NAME='deformed_fibres'
                IFIELD=2
                NCOORD=6
                NCOMP=NJ_TYPE(nhx,2)
                COMPONENT_NAME=COORDINATE_COMPONENTS(NCOMP,NCOORD)
              ELSE
                IF(nhx.LE.NCOMPT) THEN !geom cmpts
                  NFIELD=1
                  NCOMP=nhx
                  FIELD_NAME='deformed'
                  IFIELD=1
                  NCOORD=ITYP10(nr)
                  COMPONENT_NAME=COORDINATE_COMPONENTS(NCOMP,NCOORD)
                ELSE !hydrostatic pressure field
                  NFIELD=nhx-NCOMPT+1
                  IF(np.EQ.0) THEN
                    DO nhx2=NJ_LOC(NJL_GEOM,0,nr)+1,nhx
                      nh_aux=NH_LOC(nhx2,nx)
                      nb_aux=NBH(nh_aux,1,ne)
                      IF(NBC(nb_aux).EQ.0) NFIELD=NFIELD-1
                    ENDDO !nhx2
                  ENDIF
                  NCOMPT=1
                  NCOMP=1
                  WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                  IFIELD=3
                  NCOORD=1
                  COMPONENT_NAME=FIELD_NAME
                ENDIF
              ENDIF
            ELSE IF((ITYP2(nr,nx).EQ.3).OR. !Laplace equation
     '          (ITYP2(nr,nx).EQ.4).OR. !Helmholtz equation
     '          (ITYP2(nr,nx).EQ.5).OR. !Poisson equation
     '          (ITYP2(nr,nx).EQ.6).OR. !Linear 2nd order elliptic
     '          (ITYP2(nr,nx).EQ.7).OR. !Biharmonic equation
     '          (ITYP2(nr,nx).EQ.9)) THEN !Oxygen transport
              IF(np.GT.0) THEN
                NFIELDT=NHP(np)
              ELSE
                NFIELDT=NHE(ne)
              ENDIF
              NFIELD=nhx
              IF(nhx.EQ.1) THEN
                NCOMPT=1
                NCOMP=1
                IF(nc.EQ.1) THEN
                  FIELD_NAME='potential'
                ELSE
                  FIELD_NAME='flux'
                ENDIF
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME='potential'
              ELSE
                NCOMPT=1
                NCOMP=1
                WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME=FIELD_NAME
              ENDIF
            ELSE IF(ITYP2(nr,nx).EQ.8.AND. !Fluid mechanics
     '        ITYP3(nr,nx).NE.3) THEN !not const vol constraint (above)
C???????????to be done
              IF(np.GT.0) THEN
                NFIELDT=NHP(np)
              ELSE
                NFIELDT=NHE(ne)
              ENDIF
              NFIELD=nhx
              NCOMPT=1
              NCOMP=1
              WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
              IFIELD=3
              NCOORD=1
              COMPONENT_NAME=FIELD_NAME
            ENDIF
          ELSE IF(ITYP5(nr,nx).EQ.2) THEN !Time integration
            IF(ITYP2(nr,nx).EQ.1) THEN !Linear elasticity
              IF((NW(ne,1).EQ.1).OR. !Truss
     '          (NW(ne,1).EQ.2).OR. !Batten
     '          (NW(ne,1).EQ.11).OR. !Plane stress
     '          (NW(ne,1).EQ.12)) THEN !Plane strain
                IF(np.GT.0) THEN
                  NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                  NFIELDT=NHP(np)-NCOMPT+1
                ELSE
                  NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                  NFIELDT=NHE(ne)-NCOMPT+1
                ENDIF
                IF(nhx.LE.NCOMPT) THEN
                  NFIELD=1
                  NCOMP=nhx
                  FIELD_NAME='displacement'
                  IFIELD=3
                  NCOORD=ITYP10(nr)
                  COMPONENT_NAME=COORDINATE_COMPONENTS(NCOMP,NCOORD)
                ELSE
                  NFIELD=nhx-NCOMPT+1
                  NCOMPT=1
                  NCOMP=1
                  WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                  IFIELD=3
                  NCOORD=1
                  COMPONENT_NAME=FIELD_NAME
                ENDIF
              ELSE IF((NW(ne,1).EQ.3).OR. !Beam (Kirchhoff)
     '            (NW(ne,1).EQ.6).OR. !Plate (Kirchhoff)
     '            (NW(ne,1).EQ.7)) THEN !Shell
                IF(np.GT.0) THEN
                  NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                  NFIELDT=NHP(np)-NCOMPT+1
                ELSE
                  NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                  NFIELDT=NHE(ne)-NCOMPT+1
                ENDIF
                IF(nhx.LE.NCOMPT) THEN
                  NFIELD=1
                  NCOMP=nhx
                  FIELD_NAME='displacement'
                  IFIELD=3
                  NCOORD=ITYP10(nr)
                  COMPONENT_NAME=COORDINATE_COMPONENTS(NCOMP,NCOORD)
                ELSE
C???????????????to be done
                  NFIELD=nhx-NCOMPT+1
                  NCOMPT=1
                  NCOMP=1
                  WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                  IFIELD=3
                  NCOORD=1
                  COMPONENT_NAME=FIELD_NAME
                ENDIF
              ELSE
C?????????????to be done
                IF(np.GT.0) THEN
                  NFIELDT=NHP(np)
                ELSE
                  NFIELDT=NHE(ne)
                ENDIF
                NFIELD=nhx
                NCOMPT=1
                NCOMP=1
                WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME=FIELD_NAME
              ENDIF
            ELSE IF(ITYP2(nr,nx).EQ.2) THEN !Finite elasticity
              IF(np.GT.0) THEN
                NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                NFIELDT=NHP(np)-NCOMPT+1
              ELSE
                NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                NFIELDT=NHE(ne)-NCOMPT+1
              ENDIF
              IF(nhx.LE.NCOMPT) THEN
                NFIELD=1
                NCOMP=nhx
                FIELD_NAME='deformed'
                IFIELD=1
                NCOORD=ITYP10(nr)
                COMPONENT_NAME=COORDINATE_COMPONENTS(NCOMP,NCOORD)
              ELSE
                NFIELD=nhx-NCOMPT+1
                NCOMPT=1
                NCOMP=1
                WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME=FIELD_NAME
              ENDIF
            ELSE IF(ITYP2(nr,nx).EQ.3) THEN !Advection-diffusion
              IF(np.GT.0) THEN
                NFIELDT=NHP(np)
              ELSE
                NFIELDT=NHE(ne)
              ENDIF
              NFIELD=nhx
              IF(nhx.EQ.1) THEN
                NCOMPT=1
                NCOMP=1
                FIELD_NAME='concentration'
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME='concentration'
              ELSE
                NCOMPT=1
                NCOMP=1
                WRITE (FIELD_NAME,'(''nh_'',I1)') NHX
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME=FIELD_NAME
              ENDIF
            ELSE IF(ITYP2(nr,nx).EQ.4) THEN !Wave equation
              IF(np.GT.0) THEN
                NFIELDT=NHP(np)
              ELSE
                NFIELDT=NHE(ne)
              ENDIF
              NFIELD=nhx
              IF(nhx.EQ.1) THEN
                NCOMPT=1
                NCOMP=1
                FIELD_NAME='potential'
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME='potential'
              ELSE
                NCOMPT=1
                NCOMP=1
                WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME=FIELD_NAME
              ENDIF
            ELSE IF(ITYP2(nr,nx).EQ.5) THEN !Navier-Stokes equations
              IF(np.GT.0) THEN
                NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                NFIELDT=NHP(np)-NCOMPT+1
              ELSE
                NCOMPT=NJ_LOC(NJL_GEOM,0,nr)
                NFIELDT=NHE(ne)-NCOMPT+1
              ENDIF
              IF(ITYP3(nr,nx).EQ.2)THEN !lung gas transport
                NFIELD=1
                NFIELDT=1
                NCOMPT=1
                NCOMP=1
                IF(FIELD_EX_TYPE.EQ.'RADIUS')THEN
                  FIELD_NAME='radius'
                ELSE
                  FIELD_NAME='concentration'
                ENDIF
                COMPONENT_NAME=FIELD_NAME
                IFIELD=3
                NCOORD=1
              ELSE
                IF(nhx.LE.NCOMPT) THEN
                  NFIELD=1
                  NCOMP=nhx
                  FIELD_NAME='velocity'
                  IFIELD=3
                  NCOORD=ITYP10(nr)
                  COMPONENT_NAME=COORDINATE_COMPONENTS(NCOMP,NCOORD)
                ELSE
                  NFIELD=nhx-NCOMPT+1
                  NCOMPT=1
                  NCOMP=1
                  WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                  IFIELD=3
                  NCOORD=1
                  COMPONENT_NAME=FIELD_NAME
                ENDIF
              ENDIF
            ELSE IF(ITYP2(nr,nx).EQ.6) THEN !Bio-heat transfer
              IF(np.GT.0) THEN
                NFIELDT=NHP(np)
              ELSE
                NFIELDT=NHE(ne)
              ENDIF
              NFIELD=nhx
              IF(nhx.EQ.1) THEN
                NCOMPT=1
                NCOMP=1
                FIELD_NAME='concentration'
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME='concentration'
              ELSE
                NCOMPT=1
                NCOMP=1
                WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME=FIELD_NAME
              ENDIF
            ELSE IF(ITYP2(nr,nx).EQ.7) THEN !*Maxwell equations
C???????????to be done
              IF(np.GT.0) THEN
                NFIELDT=NHP(np)
              ELSE
                NFIELDT=NHE(ne)
              ENDIF
              NFIELD=nhx
              NCOMPT=1
              NCOMP=1
              WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
              IFIELD=3
              NCOORD=1
              COMPONENT_NAME=FIELD_NAME
            ELSE IF(ITYP2(nr,nx).EQ.8) THEN !Huygens activation
              ERROR='>>Not implemented'
              GOTO 9999
            ELSE IF(ITYP2(nr,nx).EQ.9) THEN !cellular based modelling
              IF(ITYP19(nr,nx).EQ.1) THEN !electrical
                IF(ITYP3(nr,nx).EQ.1) THEN !Cubic - norecovery
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='potential'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='potential'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ELSE IF(ITYP3(nr,nx).EQ.2) THEN !FitzHugh-Nagumo equations
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='activation'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='activation'
                  ELSE IF(nhx.EQ.2) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='recovery'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='recovery'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ELSE IF(ITYP3(nr,nx).EQ.3) THEN !van Capelle-Durrer
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='potential'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='potential'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''ionic concentration'',I1)')
     '                nhx-1
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ELSE IF(ITYP3(nr,nx).EQ.4) THEN !Beeler-Reuter
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='potential'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='potential'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''ionic concentration'',I1)')
     '                nhx-1
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ELSE IF(ITYP3(nr,nx).EQ.5) THEN !Jafri-Rice-Winslow
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='potential'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='potential'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''ionic concentration'',I1)')
     '                nhx-1
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ELSE IF(ITYP3(nr,nx).EQ.6) THEN !Luo-Rudy
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='potential'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='potential'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''ionic concentration'',I1)')
     '                nhx-1
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ELSE IF(ITYP3(nr,nx).EQ.7) THEN !DFN
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='potential'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='potential'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''ionic concentration'',I1)')
     '                nhx-1
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ELSE IF(ITYP3(nr,nx).EQ.8) THEN !Noble98
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='potential'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='potential'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''ionic concentration'',I1)')
     '                nhx-1
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ELSE IF(ITYP3(nr,nx).EQ.9) THEN !Hodgkin-Huxley
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='potential'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='potential'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''ionic concentration'',I1)')
     '                nhx-1
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ELSE IF(ITYP3(nr,nx).EQ.10) THEN !User defined
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='potential'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='potential'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''ionic concentration'',I1)')
     '                nhx-1
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ELSE
                  ERROR='>>Ionic current model not implemented'
                  GOTO 9999
                ENDIF
                IF(NFIELDT.EQ.0) NFIELDT=1 !NHE/NHP not set in bidomain
              ELSEIF(ITYP19(nr,nx).EQ.7) THEN !coupled
                IF(ITYP3(nr,nx).LE.2) THEN !N98/LR and HMT
                  IF(np.GT.0) THEN
                    NFIELDT=NHP(np)
                  ELSE
                    NFIELDT=NHE(ne)
                  ENDIF
                  NFIELD=nhx
                  IF(nhx.EQ.1) THEN
                    NCOMPT=1
                    NCOMP=1
                    FIELD_NAME='potential'
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME='potential'
                  ELSE
                    NCOMPT=1
                    NCOMP=1
                    WRITE (FIELD_NAME,'(''ionic concentration'',I1)')
     '                nhx-1
                    IFIELD=3
                    NCOORD=1
                    COMPONENT_NAME=FIELD_NAME
                  ENDIF
                ENDIF
              ENDIF !ITYP19(nr,nx)
            ELSE IF(ITYP2(nr,nx).EQ.10) THEN !Oxygen transport
              IF(np.GT.0) THEN
                NFIELDT=NHP(np)
              ELSE
                NFIELDT=NHE(ne)
              ENDIF
              NFIELD=nhx
              NCOMPT=1
              NCOMP=1
              WRITE (FIELD_NAME,'(''concentration'',I1)') nhx
              IFIELD=3
              NCOORD=1
              COMPONENT_NAME=FIELD_NAME
            ELSE IF(ITYP2(nr,nx).EQ.11) THEN !Pulmonary transport
              IF(np.GT.0) THEN
                NFIELDT=NHP(np) !number of fields
              ELSE
                NFIELDT=NHE(ne) !number of fields
              ENDIF
              NFIELD=nhx
              NCOMPT=1
              NCOMP=1
              IF(ITYP3(nr,nx).EQ.1)THEN !gas mixing
                FIELD_NAME='concentration'
                COMPONENT_NAME='concentration'
              ELSE IF(ITYP3(nr,nx).EQ.2)THEN !water+heat
                IF(nhx.EQ.1)THEN
                  FIELD_NAME='temperature'
                  COMPONENT_NAME='temperature'
                ELSE
                  FIELD_NAME='humidity'
                  COMPONENT_NAME='humidity'
                ENDIF
              ELSE IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.4.
     '          OR.ITYP3(nr,nx).EQ.6) THEN !capillary blood flow, P-F-R, art-cap-ven blood flow
                IF(nhx.EQ.1) THEN
                  FIELD_NAME='pressure'
                  COMPONENT_NAME='pressure'
C!!!  Need to add auxillary field name 'flow'
                ENDIF
              ELSE IF(ITYP3(nr,nx).EQ.5)THEN !MCT
                FIELD_NAME='mucus_area'
                COMPONENT_NAME='mucus_area'
              ENDIF !ITYP3
              IFIELD=3
              NCOORD=1
            ENDIF !ITYP2
          ELSE IF(ITYP5(nr,nx).EQ.3) THEN !Modal analysis
C?????????to be done
            IF(np.GT.0) THEN
              NFIELDT=NHP(np)
            ELSE
              NFIELDT=NHE(ne)
            ENDIF
            NFIELD=nhx
            NCOMPT=1
            NCOMP=1
            WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
            IFIELD=3
            NCOORD=1
            COMPONENT_NAME=FIELD_NAME
c cpb 1/5/95 Replacing Fourier analysis with Quasi-static analysis
c        ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Fourier analysis
cC?????????to be done
c          IF(np.GT.0) THEN
c            NFIELDT=NHP(np)
c          ELSE
c            NFIELDT=NHE(ne)
c          ENDIF
c          NFIELD=nhx
c          NCOMPT=1
c          NCOMP=1
c          WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
c          IFIELD=3
c          NCOORD=1
c          COMPONENT_NAME=FIELD_NAME

          ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Quasi-static analysis
            IF(ITYP2(nr,nx).EQ.3) THEN !Laplace's equation
              IF(np.GT.0) THEN
                NFIELDT=NHP(np)
              ELSE
                NFIELDT=NHE(ne)
              ENDIF
              NFIELD=nhx
              IF(nhx.EQ.1) THEN
                NCOMPT=1
                NCOMP=1
                FIELD_NAME='potential'
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME='potential'
              ELSE
                NCOMPT=1
                NCOMP=1
                WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
                IFIELD=3
                NCOORD=1
                COMPONENT_NAME=FIELD_NAME
              ENDIF
            ENDIF
          ELSE IF(ITYP5(nr,nx).EQ.5) THEN !Wavefront analysis
            NFIELDT=1
            NFIELD=1
            IF(np.GT.0) THEN
              NCOMPT=NHP(np)
            ELSE
              NCOMPT=NHE(ne)
            ENDIF
            NCOMP=nhx
            WRITE (FIELD_NAME,'(''solution'')')
            IFIELD=3
            NCOORD=1
            WRITE (COMPONENT_NAME,'(''nh_'',I1)') nhx
          ELSE IF(ITYP5(nr,nx).EQ.6) THEN !Buckling analysis
C?????????to be done
            IF(np.GT.0) THEN
              NFIELDT=NHP(np)
            ELSE
              NFIELDT=NHE(ne)
            ENDIF
            NFIELD=nhx
            NCOMPT=1
            NCOMP=1
            WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
            IFIELD=3
            NCOORD=1
            COMPONENT_NAME=FIELD_NAME
          ENDIF
        ELSE IF(NX_TYPE(nx).EQ.NX_FIT) THEN
          IF(KTYP8.EQ.4) THEN
            IF(np.GT.0) THEN
              NFIELDT=NHP(np)
            ELSE
              NFIELDT=NHE(ne)
            ENDIF
            NFIELD=nhx
            IF(nhx.EQ.1) THEN
              NCOMPT=1
              NCOMP=1
              FIELD_NAME='potential'
              IFIELD=3
              NCOORD=1
              COMPONENT_NAME='potential'
            ELSE
              NCOMPT=1
              NCOMP=1
              WRITE (FIELD_NAME,'(''nh_'',I1)') nhx
              IFIELD=3
              NCOORD=1
              COMPONENT_NAME=FIELD_NAME
            ENDIF
          ELSE IF(KTYP8.EQ.3) THEN
            IF(np.GT.0) THEN
              NFIELDT=NHP(np)
            ELSE
              NFIELDT=NHE(ne)
            ENDIF
            NFIELD=nhx
            NCOMPT=1
            NCOMP=1
            FIELD_NAME='field1'
            IFIELD=3
            NCOORD=1
            COMPONENT_NAME='field1'
          ELSE
            ERROR='>>This fit type is not implemented'
            GOTO 9999
          ENDIF
        ELSE
          ERROR='>>This nx type is not implemented'
          GOTO 9999
        ENDIF
      ENDIF
      IF((NJ.EQ.1).OR.(nhx.EQ.1).AND..NOT.DEFFIBS) THEN
C**     write the number of fields
        IF(DATAFILE) THEN
          WRITE(IFILE,'( '' #Fields='',I1 )') NFIELDT
        ELSE
          IF(FSKWRITE(NFIELDT,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
        ENDIF
      ENDIF
C!!! DB.  Temporary.  Stops auxilliary basis export
      IF(NFIELD.LE.NFIELDT) THEN
        IF(nb.EQ.0) THEN
          WRITE_HEADING=.TRUE.
        ELSE
          IF(NBC(nb).NE.0) THEN
            WRITE_HEADING=.TRUE.
          ELSE
            WRITE_HEADING=.FALSE.
          ENDIF
        ENDIF
      ELSE
        WRITE_HEADING=.FALSE.
      ENDIF
C GMH 9/2/97 This is a patch to make this (old) code compatible with
C            the cmgui link.  This code is only enabled by use of the
C             <autoname> flag.  The field names should be consistent,
C            with the exception of analytic and flux flags.  These
C            are treated as special cases here.
      IF(AUTONAME) THEN
C       Check for a scalar field, and name the component 'value'
        IF((NCOMPT.EQ.1).AND.(NCOORD.EQ.1)) THEN !1 compo, cartesian
          COMPONENT_NAME='value'
        ENDIF !scalar
        IF(iy.EQ.7) THEN !analytic
          TEMP_FIELD_NAME=FIELD_NAME
          CALL STRING_TRIM(TEMP_FIELD_NAME,IBEG1,IEND1)
          FIELD_NAME='analytic_'//TEMP_FIELD_NAME(IBEG1:IEND1)
        ENDIF
        CALL STRING_TRIM(FIELD_NAME,IBEG1,IEND1)
        WRITE (FIELD_NAME(IEND1+1:),
     '    '(''('',I1,''_'',I1,''_'',I1,''_'',I1,'')'')')
     '    nx,nr,nc,iy
      ENDIF
      IF(WRITE_HEADING) THEN
        IF(NCOMP.EQ.1) THEN
C**       write field heading
!NEWS AJP 31/5 Adding user-specified field names
          IF(SET_FIELD_NAME) THEN
C AJS MAR-2007 Initialise FIELD_NAME as blank (prevent bad output)
            FIELD_NAME(1:80)=' '
            FIELD_NAME(1:50)=INPUT_FIELD_NAME(1:50)
          ENDIF
!newe
          CALL STRING_TRIM(FIELD_NAME,IBEG1,IEND1)
          IF(DATAFILE) THEN
            WRITE(CHAR2,'(I2)') NCOMPT
            CALL STRING_TRIM(CHAR2,IB2,IE2)
            CALL STRING_TRIM(FIELD_TYPE(IFIELD),IBEG2,IEND2)
            CALL STRING_TRIM(COORDINATE_SYSTEM(NCOORD),IBEG3,IEND3)
            IF((NCOORD.EQ.4).OR.(NCOORD.EQ.5)) THEN
              WRITE(IFILE,'( I2,'') '',A,'', '',A,'', '',A,'
     '          //''', focus='',E12.4,'', #Components='',A)') NFIELD,
     '          FIELD_NAME(IBEG1:IEND1),FIELD_TYPE(IFIELD)(IBEG2:IEND2),
     '          COORDINATE_SYSTEM(NCOORD)(IBEG3:IEND3),FOCUS,
     '          CHAR2(IB2:IE2)
            ELSE
              WRITE(IFILE,'( I2,'') '',A,'', '',A,'', '',A,'
     '          //''', #Components='',A)') NFIELD,
     '          FIELD_NAME(IBEG1:IEND1),FIELD_TYPE(IFIELD)(IBEG2:IEND2),
     '          COORDINATE_SYSTEM(NCOORD)(IBEG3:IEND3),CHAR2(IB2:IE2)
            ENDIF
          ELSE
            CLEN=FSKLEN(FIELD_NAME(IBEG1:IEND1))
            CALL FSKF2C(FIELD_NAME(IBEG1:IEND1),CLEN,INTSTR)
            IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
            IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1) GOTO 9999
            IF(FSKWRITE(IFIELD,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
            IF(FSKWRITE(NCOORD,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
C?????????DB.  Sort out focus
            IF((NCOORD.EQ.4).OR.(NCOORD.EQ.5)) THEN
              IF(FSKWRITE(1.D0,SK_DOUBLE_FLOAT,1,CONNID2).EQ.-1)
     '          GOTO 9999
            ENDIF
            IF(FSKWRITE(NCOMPT,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
          ENDIF
        ENDIF
        IF(NP.GT.0) THEN
C**       write node field component heading
          CALL STRING_TRIM(COMPONENT_NAME,IBEG1,IEND1)
          IF(DATAFILE) THEN

C new CS 3/9/98 Added derivative string output
            IF(NDERIV.GT.0) THEN

              DERIVATIVE_STRING='('
              CALL STRING_TRIM(DERIVATIVE_STRING,IBEG5,IEND5)

              DO nk=1,NDERIV
                  DERIVATIVE_STRING=DERIVATIVE_STRING(IBEG5:IEND5)
     '              //DERIVATIVES(nk)
                  CALL STRING_TRIM(DERIVATIVE_STRING,IBEG5,IEND5)
              ENDDO !nk

              DERIVATIVE_STRING=DERIVATIVE_STRING(IBEG5:IEND5-1)//')'
              CALL STRING_TRIM(DERIVATIVE_STRING,IBEG5,IEND5)

            ELSE
              DERIVATIVE_STRING=' '
              CALL STRING_TRIM(DERIVATIVE_STRING,IBEG5,IEND5)
            ENDIF

            IF(IVALUE.LT.100) THEN
              ISIZE='2'
            ELSE IF (IVALUE.LT.1000) THEN
              ISIZE='3'
            ELSE
              ISIZE='4'
            ENDIF
            IF(NVERSIONS.GT.1) THEN
              WRITE(IFILE,'(3X,A,''.  Value index='',I'//ISIZE//
     '          ','', #Derivatives='',I2,'' '',A,'', #Versions='',I2)')
     '          COMPONENT_NAME(IBEG1:IEND1),IVALUE,NDERIV,
     '          DERIVATIVE_STRING(IBEG5:IEND5),NVERSIONS
            ELSE
              WRITE(IFILE,'(3X,A,''.  Value index='',I'//ISIZE//
     '          ','', #Derivatives='',I2,'' '',A)')
     '          COMPONENT_NAME(IBEG1:IEND1),IVALUE,NDERIV,
     '          DERIVATIVE_STRING(IBEG5:IEND5)
            ENDIF
          ELSE
            CLEN=FSKLEN(COMPONENT_NAME(IBEG1:IEND1))
            CALL FSKF2C(COMPONENT_NAME(IBEG1:IEND1),CLEN,INTSTR)
            IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
            IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1) GOTO 9999
            IF(FSKWRITE(IVALUE,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
            IF(FSKWRITE(NDERIV,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
            IF(FSKWRITE(NVERSIONS,SK_LONG_INT,1,CONNID2).EQ.-1)
     '        GOTO 9999
          ENDIF
        ELSE
C**       write element field component heading
C**       write coordinate field component
          CALL STRING_TRIM(COMPONENT_NAME,IBEG1,IEND1)
          IF(.NOT.DATAFILE) THEN
            CLEN=FSKLEN(COMPONENT_NAME(IBEG1:IEND1))
            CALL FSKF2C(COMPONENT_NAME(IBEG1:IEND1),CLEN,INTSTR)
            IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
            IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1) GOTO 9999
          ENDIF
          IF(FIELD_BASE_TYPE.EQ.1) THEN !Node Based
            BASE_TYPE='standard node based'
            CALL CALCULATE_BASIS_STRING(IBT,nb,DATAFILE,BASES,
     '        SPECIAL_BASIS_FLAG,ERROR,*9999)
          ELSE IF(FIELD_BASE_TYPE.EQ.2) THEN ! Grid Based
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
          ELSE
            ERROR='>>Invalid field base type'
            GOTO 9999
          ENDIF
          IF(DATAFILE) THEN
            CALL STRING_TRIM(BASES,IBEG2,IEND2)
            CALL STRING_TRIM(BASE_TYPE,IBEG3,IEND3)
          ENDIF
          IF((ELEMENT_DIMENSION.EQ.3).AND.(NCOMP.EQ.2).AND.
     '      ((NCOORD.EQ.2).OR.(NCOORD.EQ.3))) THEN
            IF(DATAFILE) THEN
              WRITE(IFILE,'(3X,A,''.  '',A,'', increasing in xi1, '','
     '          //'A,''.'')') COMPONENT_NAME(IBEG1:IEND1),
     '          BASES(IBEG2:IEND2),BASE_TYPE(IBEG3:IEND3)
            ELSE
              ISOCKET=1
              IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '          GOTO 9999
            ENDIF
C?????????DB.  phi for spherical ?
          ELSE
            IF(((ELEMENT_DIMENSION.EQ.3).OR.
     '        (ELEMENT_DIMENSION.EQ.1)).AND.
     '        (NCOMP.EQ.3).AND.
     '        ((NCOORD.EQ.4).OR.(NCOORD.EQ.5))) THEN
              IF(DATAFILE) THEN
                WRITE(IFILE,'(3X,A,''.  '',A,'', decreasing in xi1, '','
     '            //'A,''.'')') COMPONENT_NAME(IBEG1:IEND1),
     '            BASES(IBEG2:IEND2),BASE_TYPE(IBEG3:IEND3)
              ELSE
                ISOCKET=-1
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
              ENDIF
            ELSE
              IF(DATAFILE) THEN
                WRITE(IFILE,'(3X,A,''.  '',A,'', no modify, '','
     '            //'A,''.'')') COMPONENT_NAME(IBEG1:IEND1),
     '            BASES(IBEG2:IEND2),BASE_TYPE(IBEG3:IEND3)
              ELSE
                ISOCKET=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('EXPORT_FIELD_HEADING')
      RETURN
 9999 CALL ERRORS('EXPORT_FIELD_HEADING',ERROR)
      CALL EXITS('EXPORT_FIELD_HEADING')
      RETURN 1
      END


