      SUBROUTINE EXGAUS(NBJ,NEELEM,NRLIST,NXLIST,CGE,
     '  XIG,YG,STRING,ERROR,*)

C#### Subroutine: EXGAUS
C###  Description:
C###    EXGAUS exports Gauss points and constitutive parameters

C***  Created by CS 1 March 2000. Adapted from EXDATA.

      IMPLICIT none
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NRLIST(0:NRM),
     &  NXLIST(0:NXM)
      REAL*8 CGE(NMM,NGM,NEM,NXM),XIG(NIM,NGM,NBM),
     &  YG(NIYGM,NGM,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,IFROMC,il,
     &  N3CO,nb,ne,ng,ni,niyg,nr,node,noelem,nx,OFFSET,
     &  offset_elem
      CHARACTER FILE*100,POINT_NAME*50
      LOGICAL ALL_REGIONS,CBBREV,PARAMS,EXPORT_YG

      CALL ENTERS('EXGAUS',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM export gauss<;FILENAME[default]>
C###  Parameter:      <offset OFFSET[0]>
C###    Add OFFSET to Gauss point numbers
C###  Parameter:      <as NAME[gauss_points]>
C###    Add OFFSET to element numbers
C###  Parameter:      <as NAME[gauss_points]>
C###    Name of the Gauss point set.
C###  Parameter:      <parameters>
C###    Outputs Guass point constitutive parameters.
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###    Limit to Gauss points of specified regions.
C###  Description:
C###    Exports Gauss points and constitutive parameters.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']>'
        OP_STRING(2)=BLANK(1:15)//'<parameters>'
        OP_STRING(3)=BLANK(1:15)//'<offset OFFSET[0]>'
        OP_STRING(4)=BLANK(1:15)//'<offset_elem OFFSET[0]>'
        OP_STRING(5)=BLANK(1:15)//'<as NAME[gauss_points]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM export gauss<;FILENAME[default]> yg
C###  Parameter:      <offset OFFSET[0]>
C###    Add OFFSET to Gauss point numbers
C###  Parameter:      <offset_elem OFFSET[0]>
C###    Add OFFSET to element numbers
C###  Parameter:      <as NAME[gauss_points]>
C###    Name of the Gauss point set.
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###    Limit to Gauss points of specified regions.
C###  Description:
C###    Exports Gauss points and constitutive parameters.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     &    //']> yg'
        OP_STRING(2)=BLANK(1:15)//'<offset OFFSET[0]>'
        OP_STRING(3)=BLANK(1:15)//'<offset_elem OFFSET[0]>'
        OP_STRING(4)=BLANK(1:15)//'<as NAME[gauss_points]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXPOIN',ERROR,*9999)
      ELSE

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nx=NXLIST(1)

        IF(CBBREV(CO,'OFFSET_ELEM',8,noco+1,NTCO,N3CO)) THEN
          offset_elem=IFROMC(CO(N3CO+1))
        ELSE
          offset_elem=0
        ENDIF

        IF(CBBREV(CO,'OFFSET',2,noco+1,NTCO,N3CO)) THEN
          OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          OFFSET=0
        ENDIF

        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          POINT_NAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          POINT_NAME='gauss_points'
        ENDIF

        IF(CBBREV(CO,'PARAMETERS',1,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          PARAMS=.TRUE.
        ELSE
          PARAMS=.FALSE.
        ENDIF

        IF(CBBREV(CO,'YG',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          EXPORT_YG=.TRUE.
        ELSE
          EXPORT_YG=.FALSE.
        ENDIF

        CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        CALL STRING_TRIM(FILE,IBEG,IEND)
        CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.exdata','NEW',
     '    'SEQUEN','FORMATTED',132,ERROR,*9999)

C       Write the group name
        CALL STRING_TRIM(POINT_NAME,IBEG,IEND)
        WRITE(IFILE,'( '' Group name: '',A)')
     '    POINT_NAME(IBEG:IEND)
C       Write the number of fields
        IF(.NOT.EXPORT_YG) THEN
          IF(.NOT.PARAMS) THEN
            WRITE(IFILE,'(1X,''#Fields=1'')')
          ELSE
            WRITE(IFILE,'(1X,''#Fields='',I3)') (ILT(1,nr,nx)+1)
          ENDIF
        ELSE
          WRITE(IFILE,'(1X,''#Fields= 7'')') ! 6 stress/strain components
        ENDIF
C       Write the field heading
        WRITE(IFILE,'(1X,''1) element_xi, field, '
     '    //'element_xi, #Components=1'')')
        WRITE(IFILE,'(1X,''  1.  Value index= 1, '
     '    //'#Derivatives=0'')')
        IF(PARAMS) THEN
          DO il=1,ILT(1,nr,nx)!loop over all constitutive parameters
            WRITE(IFILE,'(1X,I3,'') '',I3,'', coordinate, '
     '        //'rectangular cartesian, #Components=1'')') (il+1),il
            WRITE(IFILE,'(1X,''  1.  Value index= 1, '
     '        //'#Derivatives=0'')')
          ENDDO ! il
        ELSE IF (EXPORT_YG) THEN
          DO niyg=1,6 !loop over 6 stress/strain components
            WRITE(IFILE,'(I2,'') yg'',I1,'', field, '
     '        //'real, #Components=1'')')
     '        (niyg+1),niyg
            WRITE(IFILE,'(1X,''  1.  Value index= 1, '
     '        //'#Derivatives=0'')')
          ENDDO ! il

        ENDIF
C       Loop over elements and Gauss points writing each
C       Gauss point as an embedded node
        node=0+OFFSET
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          DO ng=1,NGT(nb)
            node=node+1
            WRITE(IFILE,'(1X,''Node: '',I9)') node
            WRITE(IFILE,'(1X,''E '',I9,I2,3E13.5)')
     '        ne+offset_elem,NIT(nb),(XIG(ni,ng,nb),ni=1,NIT(nb))
            IF(PARAMS) THEN
              DO il=1,ILT(1,nr,nx) !loop over constitutive parameters
                WRITE(IFILE,'(1X,E13.5)') CGE(il,ng,ne,nx)
              ENDDO ! il
            ELSE IF (EXPORT_YG) THEN
              DO niyg=1,6 !loop over 6 stress/strain components
                WRITE(IFILE,'(1X,E13.5)') YG(niyg,ng,ne)
              ENDDO ! niyg
            ENDIF
          ENDDO !ng
        ENDDO !noelem (ne)
        CALL CLOSEF(IFILE,ERROR,*9999)
      ENDIF

      CALL EXITS('EXGAUS')
      RETURN
 9999 CALL ERRORS('EXGAUS',ERROR)
      CALL EXITS('EXGAUS')
      RETURN 1
      END


C*** 8-JUN-1999 LKC EXGEOM -- all map3d code -- Archived
C*** 8-JUN-1999 LKC EXGEOM_NE -- all map3d code -- Archived


