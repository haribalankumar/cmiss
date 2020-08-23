      SUBROUTINE UPFG(UPDATE,CP,NBH,NEELEM,NKH,NKJ,NPLIST,NPNODE,NRLIST,
     '  NVHP,NVJP,NYNP,OPERATION,PART2,XAB,XP,YG,YP,ZD,STRING,ERROR,*)

C#### Subroutine: UPFG
C###  Description:
C###    This subroutine is used by upgeom, upfiel, upsolu
C###    and upmate to substitute, add or multiply the value and
C###    derivatives between: XP,YP,CP; XE,YE,CE; XQ,YQ,CQ; XG,YG,CG.
C###    This subroutine has been developed to provide a consistent
C###    framework for the common task of operating on data between major
C###    arrays in cmiss. Programmers in the past have written code to
C###    operate on data between arrays which is specific to their demands,
C###    e.g. move the solution into a field. This results in a lot of
C###    very specific code that is difficult to extend. With a little
C###    planning addition, subtraction, multiplication, and division
C###    could've been included. Also specific code lead to inconsistant
C###    command syntax to do similar task, e.g. fem upd field solution
C###    and fem upd solution from field. This framework will provide a
C###    standard syntax.
C###    IMPLEMENTED ([OP] = substitute/add/subtract/multiply/divide):
C###    XP/YP/CP = XP/YP/CP [OP] XP/YP/CP/constant
C###    YG = YG [OP] YG/constant
C###    TO BE IMPLEMENTED:
C###    YG = XG/YG/CG [OP] XG/CG/constant
C###    XE/YE/CE = XE/YE/CE [OP] XE/YE/CE/constant
C###    XQ/YQ/CQ = XQ/YQ/CQ [OP] XQ/YQ/CQ/constant
C###


      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
!     Parameter List
      INTEGER INDICES(10,2),PART2,NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),NPLIST(0:NPM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJP(NJM,NPM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8  CP(NMM,NPM,NXM),XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),
     &  YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),ZD(NJM,NDM)
      CHARACTER OPERATION*16,STRING*(MXCH),UPDATE*16,ERROR*(*)
!      LOGICAL
!     Local Variables
      INTEGER i,IBEG,IEND,IFROMC,ILT_TOT2,N3CO,nc,NELISTL(0:NEM),
     '  nhc1,nhc2,nk,nmc1,nmc2,nonode,nonr,NPLIST_LOCAL(0:1000),
     '  nr,NREMAIN,NSE,NTIMES,NUMVALUES,nx1,nx2,nxc1,nxc2
      REAL*8 CONST,RFROMC,SCALE
      LOGICAL ALL_REGIONS,CBBREV
      CHARACTER ATPOINTS*16,FROM*16,TO*16

      CALL ENTERS('UPFG',*9999)

      IF(UPDATE(1:4).EQ.'HELP') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C---------------------------------------------------------------------

C#### Command: FEM update geometry
C###  Parameter:      ->    (x/y/z)
C###    Specify the material parameter number to update.
C###  Parameter:      ->    <region #[1]>
C###    Specify the region number to update.
C###  Parameter:      <nodes/*elements/*gauss_pts/*grid_pts[nodes]>
C###    Specify at which points to update
C###  Parameter:      <all,value,dxi1,dxi2,dxi1dxi2,dxi3,
C###                     dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>
C###    Specify which value or derivative to use
C###  Parameter:      substitute/add/substract/multiply/divide
C###    Specify how to update
C###  Parameter:      ->    <constant VALUE#>
C###    Specify to update from a constant and the value of the constant
C###  Parameter:      ->    <geometry (x/y/z)>
C###    Specify to update from a geometric field and the dimension
C###  Parameter:      ->    <field FIELD#>
C###    Specify to update from a field and the field number
C###  Parameter:      ->    <solution>
C###    Specify to update from a solution field
C###  Parameter:      ->        <niy #[1]>
C###    Specify the niy index in the solution field
C###  Parameter:      ->        <depvar #[1]>
C###    Specify the dependent variable number to update from
C###  Parameter:      ->        <class #[1]>
C###    Specify the class number to update from
C###  Parameter:      ->    <*material parameter #>
C###    Specify to update from a material field and the parameter number
C###  Parameter:      ->        <depvar #[1]>
C###    Specify the dependent variable number
C###  Parameter:      ->        <class #[1]>
C###    Specify the class number
C###  Parameter:      <all,value,dxi1,dxi2,dxi1dxi2,dxi3,
C###                     dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>
C###    Specify which value or derivative to use
C###  Parameter:      <scale_factor VALUE[1.0]>
C###    Specify what to scale the FROM array before updating the TO array
C###  Description: Updates the specified geometric field with
C###    a constant, a geometric field, a field, a solution or a material
C###    field. These updating array can substitute, add, subtract, multiply
C###    or divide the field being updated. The array that is
C###    updating the geometric field can be scaled by a scale factor
C###    before the updating operation.
C###    The update can occur on the node, element, gauss or grid points
C###    parameter array.
C###    Use 'fem list material' the list the material parameters
C###    to operate on.
C###    NOTE: A number of these features have not been implemented. These
C###    features will be implements as they are needed.The '*' indicate
C###    the features which have not been implemented.
C###

C---------------------------------------------------------------------

C#### Command: FEM update fibre
C###  Parameter:      ->    FIBRE#
C###    Specify the dependent variable number to update.
C###  Parameter:      ->    <region #[1]>
C###    Specify the class number to update.
C###  Parameter:      <nodes/*elements/*gauss_pts/*grid_pts[nodes]>
C###    Specify at which points to update
C###  Parameter:      <all,value,dxi1,dxi2,dxi1dxi2,dxi3,
C###                     dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>
C###    Specify which value or derivative to use
C###  Parameter:      substitute/add/substract/multiply/divide
C###    Specify how to update
C###  Parameter:      ->    <constant VALUE#>
C###    Specify to update from a constant and the value of the constant
C###  Parameter:      ->    <geometry (x/y/z)>
C###    Specify to update from a geometric field and the dimension
C###  Parameter:      ->    <fibre FIBRE#>
C###    Specify to update from a fibre field
C###  Parameter:      ->    <field FIELD#>
C###    Specify to update from a field and the field number
C###  Parameter:      ->    <solution>
C###    Specify to update from a solution field
C###  Parameter:      ->        <niy #[1]>
C###    Specify the niy index in the solution field
C###  Parameter:      ->        <depvar #[1]>
C###    Specify the dependent variable number to update from
C###  Parameter:      ->        <class #[1]>
C###    Specify the class number to update from
C###  Parameter:      ->    <*material parameter #>
C###    Specify to update from a material field and the parameter number
C###  Parameter:      ->        <depvar #[1]>
C###    Specify the dependent variable number
C###  Parameter:      ->        <class #[1]>
C###    Specify the class number
C###  Parameter:      <all,value,dxi1,dxi2,dxi1dxi2,dxi3,
C###                     dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>
C###    Specify which value or derivative to use
C###  Parameter:      <scale_factor VALUE[1.0]>
C###    Specify what to scale the FROM array before updating the TO array
C###  Description: Updates the specified fibre field with
C###    a constant, a geometric field, a fibre, a field, 
C###    a solution or a material
C###    field. These array can substitute, add, subtract, multiply or divide
C###    the field that is being updated. The array that is
C###    updating the field can be scaled by a scale factor before the
C###    updating operation.
C###    The update can occur at the node, element, gauss points
C###    or grid points.
C###    Use 'fem list material' the list the material parameters
C###    to operate on.
C###    NOTE: A number of these features have not been implemented. These
C###    features will be implements as they are needed.The '*' indicate
C###    the features which have not been implemented.
C###

C---------------------------------------------------------------------

C#### Command: FEM update field
C###  Parameter:      ->    FIELD#
C###    Specify the dependent variable number to update.
C###  Parameter:      ->    <region #[1]>
C###    Specify the class number to update.
C###  Parameter:      <nodes/*elements/*gauss_pts/*grid_pts[nodes]>
C###    Specify at which points to update
C###  Parameter:      <all,value,dxi1,dxi2,dxi1dxi2,dxi3,
C###                     dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>
C###    Specify which value or derivative to use
C###  Parameter:      substitute/add/substract/multiply/divide
C###    Specify how to update
C###  Parameter:      ->    <constant VALUE#>
C###    Specify to update from a constant and the value of the constant
C###  Parameter:      ->    <geometry (x/y/z)>
C###    Specify to update from a geometric field and the dimension
C###  Parameter:      ->    <field FIELD#>
C###    Specify to update from a field and the field number
C###  Parameter:      ->    <solution>
C###    Specify to update from a solution field
C###  Parameter:      ->        <niy #[1]>
C###    Specify the niy index in the solution field
C###  Parameter:      ->        <depvar #[1]>
C###    Specify the dependent variable number to update from
C###  Parameter:      ->        <class #[1]>
C###    Specify the class number to update from
C###  Parameter:      ->    <*material parameter #>
C###    Specify to update from a material field and the parameter number
C###  Parameter:      ->        <depvar #[1]>
C###    Specify the dependent variable number
C###  Parameter:      ->        <class #[1]>
C###    Specify the class number
C###  Parameter:      <all,value,dxi1,dxi2,dxi1dxi2,dxi3,
C###                     dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>
C###    Specify which value or derivative to use
C###  Parameter:      <scale_factor VALUE[1.0]>
C###    Specify what to scale the FROM array before updating the TO array
C###  Description: Updates the specified field with
C###    a constant, a geometric field, a field, a solution or a material
C###    field. These array can substitute, add, subtract, multiply or divide
C###    the field that is being updated. The array that is
C###    updating the field can be scaled by a scale factor before the
C###    updating operation.
C###    The update can occur at the node, element, gauss points
C###    or grid points.
C###    Use 'fem list material' the list the material parameters
C###    to operate on.
C###    NOTE: A number of these features have not been implemented. These
C###    features will be implements as they are needed.The '*' indicate
C###    the features which have not been implemented.
C###

C---------------------------------------------------------------------

C#### Command: FEM update solution
C###  Parameter:      ->    <niy #[1]>
C###    Specify the niy index to update.
C###  Parameter:      ->    <depvar #[1]>
C###    Specify the dependent variable number to update.
C###  Parameter:      ->    <region #[1]>
C###    Specify the region number to update.
C###  Parameter:      ->    <class #[1]>
C###    Specify the class number to update.
C###  Parameter:      <nodes/*elements/*gauss_pts/*grid_pts[nodes]>
C###    Specify at which points to update
C###  Parameter:      <all,value,dxi1,dxi2,dxi1dxi2,dxi3,
C###                     dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>
C###    Specify which value or derivative to use
C###  Parameter:      substitute/add/substract/multiply/divide
C###    Specify how to update
C###  Parameter:      ->    <constant VALUE#>
C###    Specify to update from a constant and the value of the constant
C###  Parameter:      ->    <geometry (x/y/z)>
C###    Specify to update from a geometric field and the dimension
C###  Parameter:      ->    <fibre FIBRE#>
C###    Specify to update from a fibre field
C###  Parameter:      ->    <field FIELD#>
C###    Specify to update from a field and the field number
C###  Parameter:      ->    <solution>
C###    Specify to update from a solution field
C###  Parameter:      ->        <niy #[1]>
C###    Specify the niy index in the solution field
C###  Parameter:      ->        <depvar #[1]>
C###    Specify the dependent variable number to update from
C###  Parameter:      ->        <class #[1]>
C###    Specify the class number to update from
C###  Parameter:      ->    <*material parameter #>
C###    Specify to update from a material field and the parameter number
C###  Parameter:      ->        <depvar #[1]>
C###    Specify the dependent variable number
C###  Parameter:      ->        <class #[1]>
C###    Specify the class number
C###  Parameter:      <all,value,dxi1,dxi2,dxi1dxi2,dxi3,
C###                     dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>
C###    Specify which value or derivative to use
C###  Parameter:      <scale_factor VALUE[1.0]>
C###    Specify what to scale the FROM array before updating the TO array
C###  Description: Updates the specified solution field with
C###    a constant, a geometric field, a field, a solution or a material
C###    field.  These updating array can substitute, add, subtract, multiply
C###    or divide the field being updated. The array that is
C###    updating the solution field can be scaled by a scale factor
C###    before the updating operation.
C###    The update can occur on the node, element, gauss or grid points
C###    parameter array.
C###    Use 'fem list material' the list the material parameters
C###    to operate on.
C###    NOTE: A number of these features have not been implemented. These
C###    features will be implements as they are needed.The '*' indicate
C###    the features which have not been implemented.
C###

C---------------------------------------------------------------------

C#### Command: FEM update material
C###  Parameter:      ->    parameter #
C###    Specify the material parameter number to update.
C###  Parameter:      ->    <depvar #[1]>
C###    Specify the dependent variable number to update.
C###  Parameter:      ->    <region #[1]>
C###    Specify the region number to update.
C###  Parameter:      ->    <class #[1]>
C###    Specify the class number to update.
C###  Parameter:      <nodes/*elements/*gauss_pts/*grid_pts[nodes]>
C###    Specify at which points to update
C###  Parameter:      <all,value,dxi1,dxi2,dxi1dxi2,dxi3,
C###                     dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>
C###    Specify which value or derivative to use
C###  Parameter:      substitute/add/substract/multiply/divide
C###    Specify how to update
C###  Parameter:      ->    <constant VALUE#>
C###    Specify to update from a constant and the value of the constant
C###  Parameter:      ->    <geometry (x/y/z)>
C###    Specify to update from a geometric field and the dimension
C###  Parameter:      ->    <fibre FIBRE#>
C###    Specify to update from a fibre field
C###  Parameter:      ->    <field FIELD#>
C###    Specify to update from a field and the field number
C###  Parameter:      ->    <solution>
C###    Specify to update from a solution field
C###  Parameter:      ->        <niy #[1]>
C###    Specify the niy index in the solution field
C###  Parameter:      ->        <depvar #[1]>
C###    Specify the dependent variable number to update from
C###  Parameter:      ->        <class #[1]>
C###    Specify the class number to update from
C###  Parameter:      ->    <*material parameter #>
C###    Specify to update from a material field and the parameter number
C###  Parameter:      ->        <depvar #[1]>
C###    Specify the dependent variable number
C###  Parameter:      ->        <class #[1]>
C###    Specify the class number
C###  Parameter:      <all,value,dxi1,dxi2,dxi1dxi2,dxi3,
C###                     dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>
C###    Specify which value or derivative to use
C###  Parameter:      <scale_factor VALUE[1.0]>
C###    Specify what to scale the FROM array before updating the TO array
C###  Description: Updates the specified material field with
C###    a constant, a geometric field, a field, a solution or a material
C###    field.  These updating array can substitute, add, subtract, multiply
C###    or divide the field being updated. The array that is
C###    updating the material field can be scaled by a scale factor
C###    before the updating operation.
C###    The update can occur on the node, element, gauss or grid points
C###    parameter array.
C###    Use 'fem list material' the list the material parameters
C###    to operate on.
C###    NOTE: A number of these features have not been implemented. These
C###    features will be implements as they are needed.The '*' indicate
C###    the features which have not been implemented.
C###


        IF(UPDATE(5:12).EQ.'GEOMETRY') THEN
          OP_STRING(1)=STRING(1:IEND)//'  (x/y/z)'
          OP_STRING(2)=BLANK(1:IEND)//'    <region #[1]>'
          NSE=2
        ELSEIF(UPDATE(5:9).EQ.'FIBRE') THEN
          OP_STRING(1)=STRING(1:IEND)//' FIBRE#'
          OP_STRING(2)=BLANK(1:IEND)//'    <region #[1]>'
          NSE=2
        ELSEIF(UPDATE(5:9).EQ.'FIELD') THEN
          OP_STRING(1)=STRING(1:IEND)//' FIELD#'
          OP_STRING(2)=BLANK(1:IEND)//'    <region #[1]>'
          NSE=2
        ELSEIF(UPDATE(5:12).EQ.'SOLUTION') THEN
          OP_STRING(1)=STRING(1:IEND)//' '
          OP_STRING(2)=BLANK(1:IEND)//'    <niy #[1]>'
          OP_STRING(3)=BLANK(1:IEND)//'    <depvar #[1]>'
          OP_STRING(4)=BLANK(1:IEND)//'    <disp/force [disp]>'
          OP_STRING(5)=BLANK(1:IEND)//'    <region #[1]>'
          OP_STRING(6)=BLANK(1:IEND)//'    <class #[1]>'
          NSE=5
        ELSEIF(UPDATE(5:12).EQ.'MATERIAL') THEN
          OP_STRING(1)=STRING(1:IEND)//' parameter #'
          OP_STRING(2)=BLANK(1:IEND)//'    <depvar #[1]>'
          OP_STRING(3)=BLANK(1:IEND)//'    <region #[1]>'
          OP_STRING(4)=BLANK(1:IEND)//'    <class #[1]>'
          OP_STRING(5)=BLANK(1:IEND)//'    <for solve/optimisation/'
     '      //'fit>'
          NSE=5        
        ENDIF

C        OP_STRING(NSE+1)=BLANK(1:IEND)//'<at (nodes/*elements/'
C     '    //'~gauss_pts/*grid_pts)[nodes]>'
C DMAL 01 JULY 2003 This ^ would cause a failed compilation under irix
C hence I added a concantination between the /* to fix it.
C I don't understand why, it didn't give a problem on other platforms.
        OP_STRING(NSE+1)=BLANK(1:IEND)//'<nodes/*elements/'
     &    //'~gauss_pts/'//'*grid_pts[nodes]>'
        OP_STRING(NSE+2)=BLANK(1:IEND)//'<all,value,dxi1,dxi2,'
     &    //'dxi1dxi2,dxi3,dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>'
        OP_STRING(NSE+3)=BLANK(1:IEND)//'substitute/add/substract/'
     &    //'multiply/divide'
        OP_STRING(NSE+4)=BLANK(1:IEND)//'    <constant VALUE#>'
        OP_STRING(NSE+5)=BLANK(1:IEND)//'    <geometry (x/y/z)>'
        OP_STRING(NSE+6)=BLANK(1:IEND)//'    <*fibre FIBRE#>'
        OP_STRING(NSE+7)=BLANK(1:IEND)//'    <field FIELD#>'
        OP_STRING(NSE+8)=BLANK(1:IEND)//'    <solution>'
        OP_STRING(NSE+9)=BLANK(1:IEND)//'         <niy #[1]>'
        OP_STRING(NSE+10)=BLANK(1:IEND)//'         <depvar #[1]>'
        OP_STRING(NSE+11)=BLANK(1:IEND)//'         <disp/force [disp]>'
        OP_STRING(NSE+12)=BLANK(1:IEND)//'         <class #[1]>'
        OP_STRING(NSE+13)=BLANK(1:IEND)//'    <material parameter #>'
        OP_STRING(NSE+14)=BLANK(1:IEND)//'         <depvar #[1]>'
        OP_STRING(NSE+15)=BLANK(1:IEND)//'         <class #[1]>'
        OP_STRING(NSE+16)=BLANK(1:IEND)//'         <for solve/'
     &    //'optimisation/fit>'
        OP_STRING(NSE+17)=BLANK(1:IEND)//'<all,value,dxi1,dxi2,'
     &    //'dxi1dxi2,dxi3,dxi1dxi3,dxi2dxi3,dxi1dxi2dxi3[all]>'
        OP_STRING(NSE+18)=BLANK(1:IEND)//'<scale_factor VALUE#[1.0]>'
        OP_STRING(NSE+19)=BLANK(1:IEND)//'~ - partially implemented'
        OP_STRING(NSE+20)=BLANK(1:IEND)//'* - not implemented'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE
C       INDICES stores nm,niy,nh,nj,nx,nr in that order
C       INDICES(i,1)=TO array indices
C       INDICES(i,2)=FROM array indices
        DO i=1,6
          INDICES(i,1)=0
          INDICES(i,2)=0
        ENDDO

        nc=1 ! set as default for updating geometry, changes to 2 for force
        
       ! Determine which points to update at
        IF(CBBREV(CO,'NODES',2,noco+1,PART2,N3CO))THEN
          ATPOINTS='NODES'
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,
     &      *9999)
C MHT   09Apr03 for updating fields at node groups.
          CALL PARSE_NODES(NPNODE,NPLIST,N3CO-1,NRLIST,PART2-2,CO,
     '      ERROR,*9999)

        ELSEIF(CBBREV(CO,'ELEMENTS',2,noco+1,PART2,N3CO))THEN
          ATPOINTS='ELEMENTS'
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,
     '      *9999)
          CALL PARSE_ELEMENTS(NEELEM,NELISTL,n3co,NRLIST,NTCO,CO,ERROR,
     '      *9999)
        ELSEIF(CBBREV(CO,'GAUSS_PTS',2,noco+1,PART2,N3CO))THEN
          ATPOINTS='GAUSS_PTS'
        ELSEIF(CBBREV(CO,'GRID_PTS',2,noco+1,PART2,N3CO))THEN
          ATPOINTS='GRID_PTS'
          ERROR='>> ''AT GRID_PTS'' not implemented'
          GOTO 9999
        ELSEIF(CBBREV(CO,'DATA',2,noco+1,PART2,N3CO))THEN
          ATPOINTS='DATA_PTS'
        ELSE
          ATPOINTS='NODES'
          NPLIST(0)=0
          DO nonr=1,NRLIST(0)
            nr=NRLIST(nonr)
            DO nonode=1,NPNODE(0,nr)
              NPLIST(nonode)=NPNODE(nonode,nr)
            ENDDO !nonode
            NPLIST(0)=NPLIST(0)+NPNODE(0,nr)
          ENDDO !nonr
        ENDIF

        IF(CBBREV(CO,'ALL',3,noco+1,PART2,N3CO)) THEN
          nk=0
        ELSEIF(CBBREV(CO,'VALUE',3,noco+1,PART2,N3CO)) THEN
          nk=1
        ELSEIF(CBBREV(CO,'DXI1',4,noco+1,PART2,N3CO)) THEN
          nk=2
        ELSEIF(CBBREV(CO,'DXI2',4,noco+1,PART2,N3CO)) THEN
          nk=3
        ELSEIF(CBBREV(CO,'DXI1DXI2',8,noco+1,PART2,N3CO)) THEN
          nk=4
        ELSEIF(CBBREV(CO,'DXI3',4,noco+1,PART2,N3CO)) THEN
          nk=5
        ELSEIF(CBBREV(CO,'DXI1DXI3',8,noco+1,PART2,N3CO)) THEN
          nk=6
        ELSEIF(CBBREV(CO,'DXI2DXI3',8,noco+1,PART2,N3CO)) THEN
          nk=7
        ELSEIF(CBBREV(CO,'DXI1DXI2DXI3',12,noco+1,PART2,N3CO)) THEN
          nk=8
        ELSE
          nk=0
        ENDIF
        INDICES(7,1)=nk

        IF(CBBREV(CO,'ALL',3,PART2,NTCO,N3CO)) THEN
          nk=0
        ELSEIF(CBBREV(CO,'VALUE',3,PART2,NTCO,N3CO)) THEN
          nk=1
        ELSEIF(CBBREV(CO,'DXI1',4,PART2,NTCO,N3CO)) THEN
          nk=2
        ELSEIF(CBBREV(CO,'DXI2',4,PART2,NTCO,N3CO)) THEN
          nk=3
        ELSEIF(CBBREV(CO,'DXI1DXI2',8,PART2,NTCO,N3CO)) THEN
          nk=4
        ELSEIF(CBBREV(CO,'DXI3',4,PART2,NTCO,N3CO)) THEN
          nk=5
        ELSEIF(CBBREV(CO,'DXI1DXI3',8,PART2,NTCO,N3CO)) THEN
          nk=6
        ELSEIF(CBBREV(CO,'DXI2DXI3',8,PART2,NTCO,N3CO)) THEN
          nk=7
        ELSEIF(CBBREV(CO,'DXI1DXI2DXI3',12,PART2,NTCO,N3CO)) THEN
          nk=8
        ELSE
          nk=0
        ENDIF
        INDICES(7,2)=nk
        
        IF((INDICES(7,1).EQ.0).AND.(INDICES(7,2).GT.0))THEN
          ERROR='>>ERROR: No TO value/derivative index defined'
          GOTO 9999
        ENDIF
        IF((INDICES(7,2).EQ.0).AND.(INDICES(7,1).GT.0))THEN
          ERROR='>>ERROR: No FROM value/derivative index defined'
          GOTO 9999
        ENDIF
        
        ! Get TO array indices
        IF(CBBREV(CO,'REGION',2,noco+1,PART2,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF
        INDICES(6,1)=nr

        IF(DOP)THEN
          WRITE(OP_STRING,'(/'' Updating : '',A16)') UPDATE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      
        IF(ATPOINTS(1:5).EQ.'NODES')THEN
          IF(UPDATE(1:8).EQ.'GEOMETRY')THEN
            TO='XP'
            IF(CBBREV(CO,'X',1,noco+1,noco+2,N3CO)) THEN
              INDICES(4,1)=NJ_LOC(NJL_GEOM,1,nr)
            ELSEIF(CBBREV(CO,'Y',1,noco+1,noco+2,N3CO)) THEN
              INDICES(4,1)=NJ_LOC(NJL_GEOM,2,nr)
            ELSEIF(CBBREV(CO,'Z',1,noco+1,noco+2,N3CO)) THEN
              INDICES(4,1)=NJ_LOC(NJL_GEOM,3,nr)
            ELSE
              ERROR='>> Unidentified (x/y/z) field in command'
              GOTO 9999
            ENDIF
          ELSEIF(UPDATE(1:5).EQ.'FIBRE')THEN
            TO='XP'
            INDICES(4,1)=NJ_LOC(NJL_FIBR,IFROMC(CO(noco+1)),nr)
          ELSEIF(UPDATE(1:5).EQ.'FIELD')THEN
            TO='XP'
            INDICES(4,1)=NJ_LOC(NJL_FIEL,IFROMC(CO(noco+1)),nr)
          ELSEIF(UPDATE(1:8).EQ.'SOLUTION')THEN
            TO='YP'
            IF(CBBREV(CO,'NIY',2,noco+1,PART2,N3CO)) THEN
              INDICES(2,1)=IFROMC(CO(N3CO+1))
            ELSE
              INDICES(2,1)=1
            ENDIF
            IF(CBBREV(CO,'DEPVAR',2,noco+1,PART2,N3CO)) THEN
              nhc1=IFROMC(CO(N3CO+1))
            ELSE
              nhc1=1
            ENDIF
            IF(CBBREV(CO,'DISPLACEMENT',4,noco+1,PART2,N3CO)) THEN
              nc=1 !update the geometry solution
            ELSEIF(CBBREV(CO,'FORCE',4,noco+1,PART2,N3CO)) THEN
              nc=2 !update the force solution
            ENDIF
            IF(CBBREV(CO,'CLASS',2,noco+1,PART2,N3CO)) THEN
              nxc1=IFROMC(CO(N3CO+1))
            ELSE
              nxc1=1
            ENDIF
            CALL NX_LOC(NX_INQUIRE,nxc1,nx1,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx1.NE.0,'Invalid class',ERROR,*9999)
            INDICES(5,1)=nx1
            INDICES(3,1)=NH_LOC(nhc1,nx1)
          ELSEIF(UPDATE(1:8).EQ.'MATERIAL')THEN
            TO='CP'
            IF(CBBREV(CO,'CLASS',2,noco+1,PART2,N3CO)) THEN
              nxc1=IFROMC(CO(N3CO+1))
            ELSE
              nxc1=1
            ENDIF
            IF(CBBREV(CO,'FOR',3,noco+1,PART2,N3CO)) THEN
              IF(CBBREV(CO,'SOLVE',3,noco+1,PART2,N3CO)) THEN
                CALL NX_LOC(NX_INQUIRE,nxc1,nx1,NX_SOLVE,ERROR,*9999)
              ELSEIF(CBBREV(CO,'FIT',3,noco+1,PART2,N3CO)) THEN
                CALL NX_LOC(NX_INQUIRE,nxc1,nx1,NX_FIT,ERROR,*9999)
              ELSEIF(CBBREV(CO,'OPTIMISATION',3,noco+1,PART2,N3CO)) THEN
                CALL NX_LOC(NX_INQUIRE,nxc1,nx1,NX_OPTI,ERROR,*9999)
              ENDIF
            ELSE
              CALL NX_LOC(NX_INQUIRE,nxc1,nx1,NX_SOLVE,ERROR,*9999)
            ENDIF
            CALL ASSERT(nx1.NE.0,'Invalid class',ERROR,*9999)
            INDICES(5,1)=nx1
            IF(CBBREV(CO,'DEPVAR',2,noco+1,PART2,N3CO)) THEN
              nhc1=IFROMC(CO(N3CO+1))
            ELSE
              nhc1=1
            ENDIF
            INDICES(3,1)=NH_LOC(nhc1,nx1)
            ILT_TOT2=ILT(1,nr,nx1)
            IF((ITYP5(nr,nx1).EQ.2).AND.(ITYP2(nr,nx1).EQ.3))
     '        ILT_TOT2=ILT(1,nr,nx1)/KTYP3A(nx1)
            IF(CBBREV(CO,'PARAMETER',2,noco+1,PART2,N3CO)) THEN
              nmc1=IFROMC(CO(N3CO+1))
            ELSE
              nmc1=1
            ENDIF
            INDICES(1,1)=(nhc1-1)*ILT_TOT2+nmc1
          ENDIF

          ! Get FROM array
          IF(CBBREV(CO,'CONSTANT',2,PART2,NTCO,N3CO))THEN
            FROM='CONSTANT'
            CONST=RFROMC(CO(N3CO+1))
            DO i=1,6
              INDICES(i,2)=INDICES(i,1)
            ENDDO
          ELSEIF(CBBREV(CO,'GEOMETRY',2,PART2,NTCO,N3CO))THEN
            FROM='XP'
            IF(CBBREV(CO,'X',1,PART2+1,PART2+2,N3CO)) THEN
              INDICES(4,2)=NJ_LOC(NJL_GEOM,1,nr)
            ELSEIF(CBBREV(CO,'Y',1,PART2+1,PART2+2,N3CO)) THEN
              INDICES(4,2)=NJ_LOC(NJL_GEOM,2,nr)
            ELSEIF(CBBREV(CO,'Z',1,PART2+1,PART2+2,N3CO)) THEN
              INDICES(4,2)=NJ_LOC(NJL_GEOM,3,nr)
            ELSE
              ERROR='>> Unidentified (x/y/z) field in command'
              GOTO 9999
            ENDIF
          ELSEIF(CBBREV(CO,'FIBRE',2,PART2,NTCO,N3CO))THEN
            FROM='XP'
            INDICES(4,2)=NJ_LOC(NJL_FIBR,IFROMC(CO(N3CO+1)),nr)
          ELSEIF(CBBREV(CO,'FIELD',2,PART2,NTCO,N3CO))THEN
            FROM='XP'
            INDICES(4,2)=NJ_LOC(NJL_FIEL,IFROMC(CO(N3CO+1)),nr)
          ELSEIF(CBBREV(CO,'SOLUTION',2,PART2,NTCO,N3CO))THEN
            FROM='YP'
            IF(CBBREV(CO,'NIY',2,PART2,NTCO,N3CO)) THEN
              INDICES(2,2)=IFROMC(CO(N3CO+1))
            ELSE
              INDICES(2,2)=1
            ENDIF
            IF(CBBREV(CO,'DEPVAR',2,PART2,NTCO,N3CO)) THEN
              nhc2=IFROMC(CO(N3CO+1))
            ELSE
              nhc2=1
            ENDIF
            ! GR The displacement option (nc=1) is already the default so
            ! there is no need to test for it.
            IF(CBBREV(CO,'FORCE',4,PART2,NTCO,N3CO)) THEN
              nc=2 !update the force solution
            !ELSE IF(CBBREV(CO,'DISPLACEMENT',4,PART2,NTCO,N3CO)) THEN
            !  nc=1 !update the geometry solution
            ENDIF
            IF(CBBREV(CO,'CLASS',2,PART2,NTCO,N3CO)) THEN
              nxc2=IFROMC(CO(N3CO+1))
            ELSE
              nxc2=1
            ENDIF
            CALL NX_LOC(NX_INQUIRE,nxc2,nx2,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx2.NE.0,'Invalid class',ERROR,*9999)
            INDICES(5,2)=nx2
            INDICES(3,2)=NH_LOC(nhc2,nx2)
          ELSEIF(CBBREV(CO,'MATERIAL',2,PART2,NTCO,N3CO))THEN
            FROM='CP'
c DMAL 01 JULY 2003 nmc2 used before being set.
c            INDICES(1,2)=(nh2-1)*ILT_TOT2+nmc2
            IF(CBBREV(CO,'CLASS',2,PART2,NTCO,N3CO)) THEN
              nxc2=IFROMC(CO(N3CO+1))
            ELSE
              nxc2=1
            ENDIF
            IF(CBBREV(CO,'FOR',3,PART2,NTCO,N3CO)) THEN
              IF(CBBREV(CO,'SOLVE',3,PART2,NTCO,N3CO)) THEN
                CALL NX_LOC(NX_INQUIRE,nxc2,nx2,NX_SOLVE,ERROR,*9999)
              ELSEIF(CBBREV(CO,'FIT',3,PART2,NTCO,N3CO)) THEN
                CALL NX_LOC(NX_INQUIRE,nxc2,nx2,NX_FIT,ERROR,*9999)
              ELSEIF(CBBREV(CO,'OPTIMISATION',3,PART2,NTCO,N3CO)) THEN
                CALL NX_LOC(NX_INQUIRE,nxc2,nx2,NX_OPTI,ERROR,*9999)
              ENDIF
            ELSE
              CALL NX_LOC(NX_INQUIRE,nxc2,nx2,NX_SOLVE,ERROR,*9999)
            ENDIF
            CALL ASSERT(nx2.NE.0,'Invalid class',ERROR,*9999)
            INDICES(5,2)=nx2
            IF(CBBREV(CO,'DEPVAR',2,PART2,NTCO,N3CO)) THEN
              nhc2=IFROMC(CO(N3CO+1))
            ELSE
              nhc2=1
            ENDIF
            INDICES(3,2)=NH_LOC(nhc2,nx2)
            ILT_TOT2=ILT(1,nr,nx2)
            IF((ITYP5(nr,nx2).EQ.2).AND.(ITYP2(nr,nx2).EQ.3))
     '        ILT_TOT2=ILT(1,nr,nx2)/KTYP3A(nx2)
            IF(CBBREV(CO,'PARAMETER',2,PART2,NTCO,N3CO)) THEN
              nmc2=IFROMC(CO(N3CO+1))
            ELSE
              nmc2=1
            ENDIF
            INDICES(1,2)=(nhc2-1)*ILT_TOT2+nmc2
          ELSE
            ERROR='>> Not implemented'
            GOTO 9999
          ENDIF
          
          NUMVALUES=NPLIST(0)
          
        ELSE IF(ATPOINTS(1:8).EQ.'DATA_PTS')THEN
          nk=0
          INDICES(7,1)=nk
          INDICES(7,2)=nk
          ! Get TO array
          IF(UPDATE(1:5).EQ.'FIELD')THEN
            TO='ZD'
c            INDICES(4,1)=NJ_LOC(NJL_FIEL,IFROMC(CO(noco+1)),nr)
            INDICES(4,1)=IFROMC(CO(noco+1))
          ENDIF
          ! Get FROM array
          IF(CBBREV(CO,'CONSTANT',2,PART2,NTCO,N3CO))THEN
            FROM='CONSTANT'
            CONST=RFROMC(CO(N3CO+1))
            DO i=1,6
              INDICES(i,2)=INDICES(i,1)
            ENDDO
          ELSEIF(CBBREV(CO,'FIELD',2,PART2,NTCO,N3CO))THEN
            FROM='ZD'
            INDICES(4,2)=NJ_LOC(NJL_FIEL,IFROMC(CO(N3CO+1)),nr)
          ENDIF

          NUMVALUES=NDT

        ELSEIF(ATPOINTS(1:9).EQ.'GAUSS_PTS')THEN
          ! Get TO array indices
          IF(UPDATE(1:8).EQ.'SOLUTION')THEN
            TO='YG'
            IF(CBBREV(CO,'NIY',2,noco+1,PART2,N3CO)) THEN
              INDICES(2,1)=IFROMC(CO(N3CO+1))
            ELSE
              INDICES(2,1)=1
            ENDIF
            IF(CBBREV(CO,'DEPVAR',2,noco+1,PART2,N3CO)) THEN
              nhc1=IFROMC(CO(N3CO+1))
            ELSE
              nhc1=1
            ENDIF
            IF(CBBREV(CO,'CLASS',2,noco+1,PART2,N3CO)) THEN
              nxc1=IFROMC(CO(N3CO+1))
            ELSE
              nxc1=1
            ENDIF
            CALL NX_LOC(NX_INQUIRE,nxc1,nx1,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx1.NE.0,'Invalid class',ERROR,*9999)
            INDICES(5,1)=nx1
            INDICES(3,1)=NH_LOC(nhc1,nx1)
          ELSE
            ERROR='>> Not implemented'
            GOTO 9999
          ENDIF

          ! Get FROM array indices
          IF(CBBREV(CO,'CONSTANT',2,PART2,NTCO,N3CO))THEN
            FROM='CONSTANT'
            CONST=RFROMC(CO(N3CO+1))
            DO i=1,6
              INDICES(i,2)=INDICES(i,1)
            ENDDO
          ELSEIF(CBBREV(CO,'SOLUTION',2,PART2,NTCO,N3CO))THEN
            FROM='YG'
            IF(CBBREV(CO,'NIY',2,PART2,NTCO,N3CO)) THEN
              INDICES(2,2)=IFROMC(CO(N3CO+1))
            ELSE
              INDICES(2,2)=1
            ENDIF
            IF(CBBREV(CO,'DEPVAR',2,PART2,NTCO,N3CO)) THEN
              nhc2=IFROMC(CO(N3CO+1))
            ELSE
              nhc2=1
            ENDIF
            IF(CBBREV(CO,'CLASS',2,PART2,NTCO,N3CO)) THEN
              nxc2=IFROMC(CO(N3CO+1))
            ELSE
              nxc2=1
            ENDIF
            CALL NX_LOC(NX_INQUIRE,nxc2,nx2,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx2.NE.0,'Invalid class',ERROR,*9999)
            INDICES(5,2)=nx2
            INDICES(3,2)=NH_LOC(nhc2,nx2)
          ELSE
            ERROR='>> Not implemented'
            GOTO 9999
          ENDIF
          NUMVALUES=NEM*NGM
          
        ELSE IF(ATPOINTS(1:8).EQ.'ELEMENTS')THEN
C Get TO array
          IF(UPDATE(1:5).EQ.'FIELD')THEN
            TO='XAB'
            INDICES(4,1)=NEJ_LOC(IFROMC(CO(noco+1)),nr)
          ENDIF
C Get FROM array
          IF(CBBREV(CO,'CONSTANT',2,PART2,NTCO,N3CO))THEN
            FROM='CONSTANT'
            CONST=RFROMC(CO(N3CO+1))
            DO i=1,6
              INDICES(i,2)=INDICES(i,1)
            ENDDO
          ELSEIF(CBBREV(CO,'FIELD',2,PART2,NTCO,N3CO))THEN
            FROM='XAB'
            INDICES(4,2)=NEJ_LOC(IFROMC(CO(N3CO+1)),nr)
          ELSE
            ERROR='>> Not implemented'
            GOTO 9999
          ENDIF
          
          NUMVALUES=NELISTL(0)
          
        ENDIF
        INDICES(6,2)=INDICES(6,1)

        ! Get scale factor value
        IF(CBBREV(CO,'SCALE_FACTOR',2,PART2,NTCO,N3CO))THEN
          SCALE=RFROMC(CO(N3CO+1))
        ELSE
          SCALE=1.0d0
        ENDIF

        IF(DOP)THEN
          WRITE(OP_STRING,'(/'' TO Array   : '',A16)') TO
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' FROM Array : '',A16)') FROM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' INDICES   [TO   FROM]:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' nm =   '',2I6)') INDICES(1,1),
     '      INDICES(1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' niy=   '',2I6)') INDICES(2,1),
     '      INDICES(2,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' nh =   '',2I6)') INDICES(3,1),
     '      INDICES(3,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' nj =   '',2I6)') INDICES(4,1),
     '      INDICES(4,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' nx =   '',2I6)') INDICES(5,1),
     '      INDICES(5,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' nr =   '',2I6)') INDICES(6,1),
     '      INDICES(6,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' nk =   '',2I6)') INDICES(7,1),
     '      INDICES(7,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(ATPOINTS.eq.'NODES')THEN !ARC ONLY IF ITS A NODE GROUP 10-11-11
C MHT 15-07-05 breaking up into smaller groups of nodes.
C              to avoid memory allocation problems
        NTIMES=0
        NREMAIN=NPLIST(0)
        DO WHILE(NREMAIN.GT.0)
          NPLIST_LOCAL(0)=MIN(NREMAIN,1000)
          DO i=1,NPLIST_LOCAL(0)
            NPLIST_LOCAL(i)=NPLIST(i+NTIMES*1000)
          ENDDO !i
          NREMAIN=NREMAIN-NPLIST_LOCAL(0)
          NUMVALUES=NPLIST_LOCAL(0)
          CALL UPFG_OPERATE(ATPOINTS,CONST,CP,FROM,INDICES,NBH,nc,
     &      NEELEM,NELISTL,NKH,NKJ,NPLIST_LOCAL,NUMVALUES,NVHP,NVJP,
     &      NYNP,OPERATION,SCALE,TO,XAB,XP,YG,YP,ZD,ERROR,*9999)
          NTIMES=NTIMES+1
        ENDDO !WHILE
        ELSE
          CALL UPFG_OPERATE(ATPOINTS,CONST,CP,FROM,INDICES,NBH,nc,
     &      NEELEM,NELISTL,NKH,NKJ,NPLIST_LOCAL,NUMVALUES,NVHP,NVJP,
     &      NYNP,OPERATION,SCALE,TO,XAB,XP,YG,YP,ZD,ERROR,*9999)
        ENDIF      
      ENDIF

      CALL EXITS('UPFG')
      RETURN
 9999 CALL ERRORS('UPFG',ERROR)
      CALL EXITS('UPFG')
      RETURN 1
      END


