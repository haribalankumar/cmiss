      SUBROUTINE UPSOLU(CP,IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,
     '  NEP,NHE,NHP,NHQ,NKHE,NKH,NKJ,NONY,NPF,NPLIST,NPLIST3,NPLIST4,
     '  NPNE,NPNODE,NPNY,NQNY,NRLIST,NRLIST2,NVHE,NVHP,NVJE,NVJP,NW,
     '  NXLIST,NYNE,NYNO,NYNP,NYNR,NYQNR,CURVCORRECT,PHI,PHI_H,SE,XAB,
     '  XIP,XP,YG,YP,YQ,YQS,ZA,ZD,ZE,ZG,ZP,STRING,ERROR,*)
C SMAR009 18/01/99 removed NP_INTERFACE,

C#### Subroutine: UPSOLU
C###  Description:
C###    <html><pre>
C###    EITHER:
C###    Updates YP array from GEOMETRY or FIELD variables.
C###    Substitutes/adds/subtracts geometry or field vars to current
C###    solution. This is useful for FE40 problems going from
C###    displacements to final geometry.
C###    OR:
C###    Transfers fitted deformed fibre angles into the YP arrray.
C###    OR:
C###    If the COUPLED option is chosen the solution stored in YP(ny,1)
C###    is copied from the ny's of the SOURCE_REGION into all ny1s
C###    that are coupled to ny during the global coupled solution.
C###    Requires coupled solution to be set up in define solve.
C###    OR:
C###    If the CAVITY_REFERENCE option is chosen YP(ny,1) is copied into
C###    YP(ny,10) which is used as the reference state for FE50 Cavity
C###    problems.
C###    OR:
C###    If the CONVERGED_REFERENCE option is chosen YP(ny,1) is copied
C###    into YP(ny,11) which is used as the reference state for FE50
C###    active tension problems.
C###    OR:
C###    If the HISTORY option is chosen either PHI or PHI_H is copied
C###    into YP(ny,1,nx) and saved to a history file.
C###    </pre></html>

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NEP(NPM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJ(NJM,NPM),NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NPF(9,NFM),NPLIST(0:NPM),NPLIST3(0:NPM),NPLIST4(0:NPM),
     '  NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NQNY(2,NYQM,0:NRCM,NXM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
C SMAR009 18/01/99 removed NP_INTERFACE(0:NPM,0:3),
      REAL*8 CP(NMM,NPM,NXM),CURVCORRECT(2,2,NNM,NEM),
     '  PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  SE(NSM,NBFM,NEM),XAB(NORM,NEM),XIP(NIM,NPM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     '  YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),ZE(NSM,NHM),
     &  ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER i,IBEG,IBEG1,IEND,IEND1,IFROMC,IY,irow,j,N3CO,na,nb,nc,ne,
     '  nh,nh1,nh2,nh_num,nh_update,NHLIST(0:9),nhx,niq,
     '  niy_f,niy_t,nj,njj,NIQLIST(0:1),NIQSLIST(0:1),NIYLIST(0:16),
     '  NJJLIST(0:9),njj_num,nk,nn,no,noelem,nolist,nonode,nonr,
     '  no_nrlist,no_nynr,notime,np,npp,np2,np_coronary,NP_average,nq,
     '  nr,nr1,nr_coronary,nr_host,nr_source,NUMTIMEDATA,nv,
     '  nx,nxc,nx_f,nx_t,ny,ny1,ny2,nyo,PART2
      REAL*8 DXIX(3,3),SUM,TIME,XI(3),YPMAX(16),YPMIN(16),Z(3),ZRC(3),
     '  ZRC_AVGE(3)
      LOGICAL ABBREV,ALL_REGIONS,AVGENODE,CBBREV,DEFFIBS,ENDFILE,
     '  EXCLUDE,NBMATCH,NKMATCH,NVMATCH,OPPHI,YPDATA,YQDATA,YQSDATA
      CHARACTER FILEFORMAT*6,OUTFILE*100,OPERATION*16,OPERATION2*10,
     '  TYPE*20,UPDATE*16

      CALL ENTERS('UPSOLU',*9999)

      CALL STRING_TRIM(FILE00,IBEG1,IEND1)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

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
C###  Parameter:      <at (nodes/*elements/*gauss_pts/*grid_pts[nodes]>
C###    Specify at which points to update
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

          UPDATE='HELPSOLUTION'
          CALL UPFG(UPDATE,%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     '      %VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     '      %VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     &      STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update solution geometry
C###  Parameter:      <(substitute/add/subtract)[substitute]>
C###   Specify the method of modifing the solution
C###  Parameter:      <YP_index IY#[1]>
C###   Specify the index of the YP array to update
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <using (fit/solve)[solve]>
C###    Specify the method of update
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    This command will modify solution values by adding or
C###    subtracting various quantities. The option add geometry
C###    is useful for fe40 problems to go from displacements to
C###    final geometry.

        OP_STRING(1)=STRING(1:IEND)//' geometry'
        OP_STRING(2)=BLANK(1:15)//
     '    '<substitute/add/subtract[substitute]>'
        OP_STRING(3)=BLANK(1:15)//'<YP_index IY#[1]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update solution field #s
C###  Parameter:      <(substitute/add/subtract)[substitute]>
C###   Specify the method of modifing the solution
C###  Parameter:      <YP_index IY#[1]>
C###   Specify the index of the YP array to update
C###  Parameter:      <nh nh#s>
C###   Specify the dependent variable numbers to update. If not
C###   specified the fields specified are mapped to the ddependent
C###   variables starting from the first dependent variable.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <using (fit/solve)[solve]>
C###    Specify the method of update
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    This command will modify solution values by adding or
C###    subtracting various quantities. The option add geometry
C###    is useful for fe40 problems to go from displacements to
C###    final geometry.

        OP_STRING(1)=STRING(1:IEND)//' field #s'
        OP_STRING(2)=BLANK(1:15)//
     '    '<substitute/add/subtract[substitute]>'
        OP_STRING(3)=BLANK(1:15)//'<YP_index IY#[1]>'
        OP_STRING(4)=BLANK(1:15)//'<nh nh#s>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update solution deformed_fibres
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Transfers fitted deformed fibre angles into the YP arrray
C###

        OP_STRING(1)=STRING(1:IEND)//' deformed_fibres'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update solution coupled
C###  Parameter:      source_region #
C###   Specify the source region
C###  Parameter:      <using (fit/solve)[solve]>
C###   Specify wiether a fit or solve is used
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    If the COUPLED option is chosen the solution stored in YP(ny,1)
C###    is copied from the ny's of the SOURCE_REGION into all ny1s
C###    that are coupled to ny during the global coupled solution.
C###    Requires coupled solution to be set up in define solve.

        OP_STRING(1)=STRING(1:IEND)//' coupled'
        OP_STRING(2)=BLANK(1:15)//'source_region #'
        OP_STRING(3)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update solution cavity_reference
C###  Parameter:      <average NODE_NUM#[1]
C###    Specify the calculation of the average NH# from
C###    the given list of deformed nodes.
C###  Parameter:         <in NH#[1]>
C###    Specify the NH# to update into
C###  Parameter:         <node (#s/GROUP/all)[all]>>
C###    Specify the node numbers
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <using (fit/solve)[solve]>
C###   Specify wiether fit or solve is used
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    If the CAVITY_REFERENCE option is chosen YP(ny,1) is copied into
C###    YP(ny,10) which is used as the reference state for FE50 Cavity
C###    problems.
C###    The 'average NODE_NUM' option is used to calculate the average
C###    NH# rectangular cartesian coordinate from the given list of
C###    deformed nodes, and then change the deformed coordinates for
C###    node NODE_NUM so that it is positioned at that averaged NH#
C###    coordinate. This is  for use in the computation of deformed
C###    cavity element volumes. The node list and NODE_NUM must
C###    belong to the specified region.

        OP_STRING(1)=STRING(1:IEND)//' cavity_reference'
        OP_STRING(2)=BLANK(1:15)//'<average NODE_NUM#[1]'
        OP_STRING(3)=BLANK(1:15)//'   <in NH#[1]>'
        OP_STRING(4)=BLANK(1:15)//'   <node (#s/GROUP/all)[all]>>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C *** DPN 28 April 2001 - Adding a update command to enable the storing
C ***   of converged solutions for use in active tension simulations

C#### Command: FEM update solution converged_reference
C###  Parameter:      <save/restore [save]>
C###    Specify whether to save YP(ny,1) to YP(ny,11) or to restore
C###    YP(ny,11) to YP(ny,1)
C###  Parameter:      <region (#s/all)[all]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    If the CONVERGED_REFERENCE option is specified, then a reference
C###    copy of YP(ny,1) is either saved to or restored from YP(ny,11).
C###    When restoring the saved copy, ZP and ZA are also restored

        OP_STRING(1)=STRING(1:IEND)//' converged_reference'
        OP_STRING(2)=BLANK(1:15)//'<save/restore [save]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update solution coronary_mesh
C###  Parameter:      <update_region (#s/all)[2]>
C###    Specify the region of the coronary mesh
C###  Parameter:      <host_region (#s/all)[1]>
C###    Specify the host region for the coronary mesh
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' coronary_mesh'
        OP_STRING(2)=BLANK(1:15)//'<update_region (#s/all)[2]>'
        OP_STRING(3)=BLANK(1:15)//'<host_region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C----------------------------------------------------------------------

C#### Command: FEM update solution dependent
C###  Parameter:      <from_class #[1]>
C###    Specify the class to update from
C###  Parameter:      <to_class #[2]>
C###    Specify the class to update
C###  Parameter:      <from_index #[1>
C###    Specify the index to update from
C###  Parameter:      <to_index #[1]>
C###    Specify the index to update to
C###  Parameter:      <region #[1]>
C###    Specify the region numbers to update.
C###  Description:
C###    Update nodal solutions from 1 class/region to
C###    another class/region, index location to another
C###    index location in YP array

        OP_STRING(1)=STRING(1:IEND)//' dependent'
        OP_STRING(2)=BLANK(1:15)//'<from_class #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<to_class #[2]>'
        OP_STRING(4)=BLANK(1:15)//'<from_index #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<to_index #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C----------------------------------------------------------------------

C#### Command: FEM update solution grid
C###  Parameter:      <from_class #[1]>
C###    Specify the class to update from
C###  Parameter:      <to_class #[2]>
C###    Specify the class to update
C###  Parameter:      <region #[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <level #[1]>
C###    Specify the grid level to update
C###  Parameter:      <niq_index #[1]>
C###    Specify the niq index to update.
C###  Description:
C###    Update grid point solutions from one class/region
C###    to another class.

        OP_STRING(1)=STRING(1:IEND)//' grid'
        OP_STRING(2)=BLANK(1:15)//'<from_class #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<to_class #[2]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<level #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<niq_index #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C----------------------------------------------------------------------

C#### Command: FEM update solution history
C###  Parameter:      <outfile FILENAME[$current]>
C###    Specify the history filename (binhis or iphist) to output the
C###    solution matrix to.
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as binary or ascii.
C###  Parameter:      <from (PHI/PHI_H)[PHI]>
C###    Specify whether the solution matrix is PHI or PHI_H.
C###  Description:
C###    Creates a history file containing either the PHI(nytr,nts) or
C###    PHI_H(nytr2,nts) matrix. The history file contains only the
C###    outer surface region for PHI or the first surface region for
C###    PHI_H defined by the single layer transfer matrix.

        OP_STRING(1)=STRING(1:IEND)//' history'
        OP_STRING(2)=BLANK(1:15)//'<outfile FILENAME['
     '    //FILE00(IBEG1:IEND1)//'_new]>'
        OP_STRING(3)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(4)=BLANK(1:15)//'<(from PHI/PHI_H)[PHI]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C----------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPSOLU',ERROR,*9999)
      ELSE
        TYPE=' '
        IF(CBBREV(CO,'SUBSTITUTE',4,noco+1,NTCO,N3CO)) THEN
          OPERATION='SUBSTITUTE'
          PART2=N3CO
          IF(.NOT.(CBBREV(CO,'GEOMETRY',2,noco+1,PART2,N3CO).OR.
     '       CBBREV(CO,'FIELD',2,noco+1,PART2,N3CO))) TYPE='OPERATION'
        ELSEIF(CBBREV(CO,'ADD',2,noco+1,NTCO,N3CO)) THEN
          OPERATION='ADD'
          PART2=N3CO
          IF(.NOT.(CBBREV(CO,'GEOMETRY',2,noco+1,PART2,N3CO).OR.
     '       CBBREV(CO,'FIELD',2,noco+1,PART2,N3CO))) TYPE='OPERATION'
        ELSEIF(CBBREV(CO,'SUBTRACT',4,noco+1,NTCO,N3CO)) THEN
          OPERATION='SUBTRACT'
          PART2=N3CO
          IF(.NOT.(CBBREV(CO,'GEOMETRY',2,noco+1,PART2,N3CO).OR.
     '       CBBREV(CO,'FIELD',2,noco+1,PART2,N3CO))) TYPE='OPERATION'
        ELSEIF(CBBREV(CO,'MULTIPLY',2,noco+1,NTCO,N3CO)) THEN
          OPERATION='MULTIPLY'
          PART2=N3CO
          IF(.NOT.(CBBREV(CO,'GEOMETRY',2,noco+1,PART2,N3CO).OR.
     '       CBBREV(CO,'FIELD',2,noco+1,PART2,N3CO))) TYPE='OPERATION'
        ELSEIF(CBBREV(CO,'DIVIDE',2,noco+1,NTCO,N3CO)) THEN
          OPERATION='DIVIDE'
          PART2=N3CO
          IF(.NOT.(CBBREV(CO,'GEOMETRY',2,noco+1,PART2,N3CO).OR.
     '       CBBREV(CO,'FIELD',2,noco+1,PART2,N3CO))) TYPE='OPERATION'
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

!Initialise DXIX
        DO i=1,3
          DO j=1,3
            DXIX(i,j)=0.0d0
          ENDDO
        ENDDO

        IF(TYPE(1:9).EQ.'OPERATION')THEN
          ! DO NOTHING : All updating done in subroutine UPFG
        ELSEIF((CBBREV(CO,'DEPENDENT',3,noco+1,NTCO,N3CO)).OR.
     '    (CBBREV(CO,'GRID',3,noco+1,NTCO,N3CO))) THEN
          IF(CBBREV(CO,'FROM_CLASS',6,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_f,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_f.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)

          IF(CBBREV(CO,'TO_CLASS',4,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=2
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_t,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_t.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)


C LKC 26-JAN-2000 Generalise for niy index locations
          IF(CBBREV(CO,'FROM_INDEX',6,noco+1,NTCO,N3CO)) THEN
            niy_f=IFROMC(CO(N3CO+1))
          ELSE
            niy_f=1
          ENDIF

          IF(CBBREV(CO,'TO_INDEX',6,noco+1,NTCO,N3CO)) THEN
            niy_t=IFROMC(CO(N3CO+1))
          ELSE
            niy_t=1
          ENDIF

        ELSE
          IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
              CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '          ERROR,*9999)
            ELSE
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this solve class',ERROR,*9999)
            ENDIF
          ELSE
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '        ERROR,*9999)
          ENDIF
        ENDIF

        DEFFIBS=.FALSE.
        IF(TYPE(1:9).EQ.'OPERATION') THEN
          ! DO NOTHING: Upadating done in subroutine UPFG
        ELSEIF(CBBREV(CO,'GEOMETRY',1,noco+1,NTCO,n3co))THEN
          TYPE='GEOMETRY'
C         Assign NHLIST(0)=0 since we are not specifying any nh in this command
          NHLIST(0)=0
C         Put njj #'s for geometry into NJJLIST
          NJJLIST(0)=NJT
          DO njj=1,NJJLIST(0)
            NJJLIST(njj)=njj
          ENDDO !njj
          IF(CBBREV(CO,'YP_INDEX',1,noco+1,NTCO,n3co)) THEN
            IY=IFROMC(CO(n3co+1))
            CALL ASSERT(IY.GT.0,'>>IY must be > 0',ERROR,*9999)
          ELSE
            IY=1
          ENDIF
        ELSEIF(CBBREV(CO,'FIELD',1,noco+1,NTCO,n3co).AND.
     '    TYPE(1:9).NE.'OPERATION') THEN
          TYPE='FIELD'
C         Put list of field #'s into NJJLIST
          CALL PARSIL(CO(n3co+1),9,NJJLIST(0),NJJLIST(1),ERROR,*9999)
          IF(NJJLIST(0).EQ.0) THEN
            CO(noco+1)='?'
            GO TO 1
          ENDIF
          DO njj=1,NJJLIST(0)
            CALL ASSERT(NJJLIST(njj).LE.NJ_LOC(NJL_FIEL,0,0),
     '        '>>Too few field variables defined',ERROR,*9999)
          ENDDO !njj
          IF(CBBREV(CO,'YP_INDEX',1,noco+1,NTCO,n3co)) THEN
            IY=IFROMC(CO(n3co+1))
            CALL ASSERT(IY.GT.0,'>>IY must be > 0',ERROR,*9999)
          ELSE
            IY=1
          ENDIF
          IF(CBBREV(CO,'NH',1,noco+1,NTCO,n3co)) THEN
            CALL PARSIL(CO(n3co+1),9,NHLIST(0),NHLIST(1),ERROR,*9999)
            IF(NHLIST(0).EQ.0) THEN
              CO(noco+1)='?'
              GO TO 1
            ENDIF
            CALL ASSERT(NJJLIST(0).EQ.NHLIST(0),'>>Number of fields must
     '        equal the number of dependent variables',ERROR,*9999)
          ELSE
            NHLIST(0)=0
          ENDIF
        ELSE IF(CBBREV(CO,'DEFORMED_FIBRES',1,noco+1,NTCO,n3co)) THEN
          DEFFIBS=.TRUE.
C         Put njj #'s for fibres into NJJLIST
          CALL ASSERT(NJ_LOC(NJL_FIBR,0,0).GT.0,'>>Fibres not defined',
     '      ERROR,*9999)
          NJJLIST(0)=NJ_LOC(NJL_FIBR,0,0)
          DO njj=1,NJJLIST(0)
            NJJLIST(njj)=njj
          ENDDO !njj
          IY=5 !temporary location for deformed fibre angles
        ELSE IF(CBBREV(CO,'CAVITY_REFERENCE',2,noco+1,NTCO,n3co)) THEN
          TYPE='CAVITY_REFERENCE'
          IY=1
          IF(CBBREV(CO,'AVERAGE',1,noco+1,NTCO,N3CO)) THEN
            NP_average=IFROMC(CO(N3CO+1))
            CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
            CALL ASSERT(NPLIST(0).GT.0,'>>No nodes in list!',
     '        ERROR,*9999)
            AVGENODE=.TRUE.
            IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
              nhx=IFROMC(CO(N3CO+1))
              nh_update=NH_LOC(nhx,nx)
              CALL ASSERT(nh_update.GT.0,'>>nh # not defined',
     '          ERROR,*9999)
            ELSE
              nh_update=1
            ENDIF
          ELSE
            AVGENODE=.FALSE.
          ENDIF
C *** DPN 28 April 2001 - Adding command to save/restore a copy of
C ***   YP(ny,1)
        ELSE IF(CBBREV(CO,'CONVERGED_REFERENCE',2,noco+1,NTCO,n3co))
     '      THEN
          TYPE='CONVERGED_REFERENCE'
          CALL ASSERT(NIYM.GE.11,'>>NIYM must be at least 11',ERROR,
     '      *9999)
          IF(CBBREV(CO,'SAVE',4,noco+1,NTCO,N3CO)) THEN
            OPERATION2='SAVE'
          ELSE IF(CBBREV(CO,'RESTORE',4,noco+1,NTCO,N3CO)) THEN
            OPERATION2='RESTORE'
          ELSE
            OPERATION2='SAVE'
          ENDIF
        ELSE IF(CBBREV(CO,'COUPLED',2,noco+1,NTCO,n3co)) THEN
          CALL ASSERT(IS_COUPLED(nx),
     '      '>>Define solve for coupled problem',ERROR,*9999)
          TYPE='COUPLED'
          IF(CBBREV(CO,'SOURCE_REGION',1,noco+1,NTCO,n3co)) THEN
            nr_source=IFROMC(CO(n3co+1))
            CALL ASSERT(nr_source.GT.0,'>>Source region cannot be '
     '        //'global coupled region (nr=0)',ERROR,*9999)
          ELSE
            CO(noco+1)='?'
            GO TO 1
          ENDIF
        ELSE IF(CBBREV(CO,'CORONARY_MESH',2,noco+1,NTCO,n3co)) THEN
          TYPE='CORONARY_MESH'
          IF(CBBREV(CO,'UPDATE_REGION',2,noco+1,NTCO,n3co)) THEN
            nr_coronary=IFROMC(CO(n3co+1))
          ELSE
            nr_coronary=2
          ENDIF
          IF(CBBREV(CO,'HOST_REGION',2,noco+1,NTCO,n3co)) THEN
            nr_host=IFROMC(CO(n3co+1))
            CALL ASSERT(nr_coronary.GT.0,'>>coronary region must be '
     '        //'greater than zero',ERROR,*9999)
          ELSE
            nr_host=1
          ENDIF
          CALL ASSERT(nr_coronary.EQ.nr_host,'>>coronary region'
     '       //'must be different from host region',ERROR,*9999)

        ELSE IF(CBBREV(CO,'DEPENDENT',3,noco+1,NTCO,N3CO)) THEN
          TYPE='DEPENDENT'
          nr=NRLIST(1)
        ELSE IF(CBBREV(CO,'GRID',3,noco+1,NTCO,N3CO).AND.
     '    TYPE(1:9).NE.'OPERATION') THEN
          TYPE='GRID'
          nr=NRLIST(1)
        ELSE IF(CBBREV(CO,'HISTORY',2,noco+1,NTCO,N3CO)) THEN
          TYPE='HISTORY'
          IF(CBBREV(CO,'OUTFILE',2,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            OUTFILE=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            OUTFILE=FILE00(IBEG1:IEND1)//'_new'
          ENDIF
          CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)

          IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF

          IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'PHI_H',4)) THEN
              CALL ASSERT(EVALUATE_INVERSE,'>>Evaluate PHI_H first',
     '          ERROR,*9999)
              OPPHI=.FALSE.
            ELSE
              CALL ASSERT(EVALUATE_PHI,'>>Evaluate PHI first',
     '          ERROR,*9999)
              OPPHI=.TRUE.
            ENDIF
          ELSE
            CALL ASSERT(EVALUATE_PHI,'>>Evaluate PHI first',ERROR,*9999)
            OPPHI=.TRUE.
          ENDIF
        ELSEIF(TYPE(1:9).EQ.'OPERATION')THEN
        ELSE
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        IF((.NOT.DEFFIBS).AND.(TYPE(1:9).NE.'OPERATION')) THEN
          IF(CBBREV(CO,'SUBSTITUTE',2,noco+1,NTCO,N3CO)) THEN
            OPERATION='SUBSTITUTE'
          ELSE IF(CBBREV(CO,'ADD',2,noco+1,NTCO,N3CO)) THEN
            OPERATION='ADD'
          ELSE IF(CBBREV(CO,'SUBTRACT',2,noco+1,NTCO,N3CO)) THEN
            OPERATION='SUBTRACT'
          ELSE
            OPERATION='SUBSTITUTE'
          ENDIF
        ELSEIF(DEFFIBS.AND.(TYPE(1:9).NE.'OPERATION')) THEN
          OPERATION='SUBSTITUTE'
        ENDIF

        IF(TYPE(1:9).EQ.'OPERATION')THEN

          UPDATE='SOLUTION'
          CALL UPFG(UPDATE,CP,NBH,NEELEM,NKH,NKJ,NPLIST,NPNODE,
     '      NRLIST,NVHP,NVJP,NYNP,OPERATION,PART2,XAB,XP,YG,YP,ZD,
     '      STRING,ERROR,*9999)

        ELSEIF(TYPE(1:7).EQ.'COUPLED') THEN
C         Update solution in specified region(s) from ny's that
C         are coupled to the source region
          DO no_nynr=1,NYNR(0,0,1,nr_source,nx)
C           loop over global ny variables of source region
            ny=NYNR(no_nynr,0,1,nr_source,nx) !global ny var number
            IF(NONY(0,ny,2,0,nx).GT.0) THEN !ny is solved for globally
              CALL ASSERT(NONY(0,ny,2,0,nx).EQ.1,'>> Cant handle '
     '          //'couplings with 1 ny to many no.s',ERROR,*9999)
              no=NONY(1,ny,2,0,nx) !solution var no for global region
              DO nyo=1,NYNO(0,no,2,0,nx)
C               loop over ny variables attached to global no variable
                ny1=NYNO(nyo,no,2,0,nx)
                nr1=NPNY(6,ny1,0,nx)
C new MPN 4Jul2003: don't want to copy YP value if ny in same region as ny1
                IF(ny1.NE.ny.AND.nr1.NE.nr_source) THEN
C old                IF(ny1.NE.ny) THEN
                  IF(NPNY(0,ny1,0,nx).EQ.1) THEN
                    np=NPNY(4,ny1,0,nx)
C GMH 8/1/97 Update cmgui link
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                  ENDIF
                  YP(ny1,1,nx)=YP(ny,1,nx)
                ENDIF
              ENDDO !nyo (ny1)
            ENDIF
          ENDDO !no_nynr (ny)

        ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN

C$OMP PARALLEL DO
C$OMP&PRIVATE(EXCLUDE,nc,nh1,nh2,nk,np,np2,npp,nv,ny1,ny2)
C$OMP&SHARED(niy_f,niy_t,NKH,NPNODE,nr,nx_f,nx_t,NYNP,YP)
          DO npp=1,NPNODE(0,nr)
            np=NPNODE(npp,nr)

            EXCLUDE=.FALSE.
            DO np2=1,CPLST(0,1)
              IF(np.EQ.CPLST(np2,1)) EXCLUDE=.TRUE.
            ENDDO

            !find ny's for each nr/nx and copy across
            nv=1
            nh1=NH_LOC(1,nx_f)
            nh2=NH_LOC(1,nx_t)
            DO nc=1,NCT(nr,nx_f)
              IF((nc.GT.1).AND.EXCLUDE) THEN
                !don't update - explicit no flux bc has been applied
              ELSE
                DO nk=1,NKH(nh1,np,nc,nr)
                  ny1=NYNP(nk,nv,nh1,np,0,nc,nr)
                  ny2=NYNP(nk,nv,nh2,np,0,nc,nr)
C LKC 26-JAN-2000 generalise for niy_f and niy_t
C                  YP(ny2,1,nx_t)=YP(ny1,1,nx_f)
                  YP(ny2,niy_t,nx_t)=YP(ny1,niy_f,nx_f)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
C$OMP END PARALLEL DO

        ELSE IF(TYPE(1:4).EQ.'GRID') THEN

          IF(CBBREV(CO,'LEVEL',3,noco+1,NTCO,N3CO)) THEN
            na=IFROMC(CO(N3CO+1))
          ELSE
            na=1
          ENDIF
          IF(CBBREV(CO,'NIQ_INDEX',3,noco+1,NTCO,N3CO)) THEN
            niq=IFROMC(CO(N3CO+1))
          ELSE
            niq=1
          ENDIF

          DO nq=NQR(1,nr),NQR(2,nr)
            YQ(nq,niq,na,nx_t)=YQ(nq,niq,na,nx_f)
          ENDDO !nq

        ELSE IF(TYPE(1:8).EQ.'HISTORY') THEN

          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)

          na=1
          NIYLIST(0)=1
          NIYLIST(1)=1
          NIQLIST(0)=0
          NIQSLIST(0)=0
          YPDATA=.TRUE.
          YQDATA=.FALSE.
          YQSDATA=.FALSE.

          IF(OPPHI) THEN
            ! ouput PHI array to history file
C***        PHI is defined as the second or outer region of the single layer
C***        (epicardial to torso) transfer matrix
            nr=TRSF_NR_OUTER
            NRLIST(0)=1
            NRLIST(1)=nr
            DO nonr=0,NRLIST(0)
              NRLIST2(nonr)=NRLIST(nonr)
            ENDDO

C***        Initialise YP(ny,1,nx)
            DO nc=1,NCT(nr,nx)
              DO no_nynr=1,NYNR(0,0,nc,nr,nx)
                ny=NYNR(no_nynr,0,nc,nr,nx)
                YP(ny,1,nx)=0.0d0
              ENDDO
            ENDDO

C***        Open up a new history file for writing
            CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '        nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '        YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',FILEFORMAT,OUTFILE,
     '        'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
            DO notime=1,NTST
              TIME=DBLE(notime)
              irow=0
              DO nolist=1,NPLIST4(0)
                np=NPLIST4(nolist)
                DO nhx=1,NHP(np,nr,nx)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,1,nr) !nc=1
                    CALL ASSERT(MAX(NKH(nh,np,1,nr)-KTYP93(1,nr),1)
     '                .EQ.1,'>>Only valid for nk = 1',ERROR,*9999)
                    irow=irow+1
                    ny=NYNP(1,nv,nh,np,0,1,nr) !nk=nc=1
                    YP(ny,1,nx)=PHI(irow,notime)
                  ENDDO !nv
                ENDDO !nhx
              ENDDO !nolist
              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '          NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',
     '          FILEFORMAT,OUTFILE,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,
     '          YQDATA,YQSDATA,ERROR,*9999)
            ENDDO !notime
          ELSE
            ! ouput PHI_H array to history file
C***        PHI_H is defined as the first region of the single layer
C***        (epicardial to torso) transfer matrix
            nr=TRSF_NR_FIRST
            NRLIST(0)=1
            NRLIST(1)=nr
            DO nonr=0,NRLIST(0)
              NRLIST2(nonr)=NRLIST(nonr)
            ENDDO

C***        Initialise YP(ny,1,nx)
            DO nc=1,NCT(nr,nx)
              DO no_nynr=1,NYNR(0,0,nc,nr,nx)
                ny=NYNR(no_nynr,0,nc,nr,nx)
                YP(ny,1,nx)=0.0d0
              ENDDO
            ENDDO

C***        Open up a new history file for writing
            CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '        nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '        YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',FILEFORMAT,OUTFILE,
     '        'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
            DO notime=1,NTST
              TIME=DBLE(notime)
              irow=0
              DO nolist=1,NPLIST3(0)
                np=NPLIST3(nolist)
                DO nhx=1,NHP(np,nr,nx)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,1,nr) !nc=1
                    CALL ASSERT(MAX(NKH(nh,np,1,nr)-KTYP93(1,nr),1)
     '                .EQ.1,'>>Only valid for nk = 1',ERROR,*9999)
                    irow=irow+1
                    ny=NYNP(1,nv,nh,np,0,1,nr) !nk=nc=1
                    YP(ny,1,nx)=PHI_H(irow,notime)
                  ENDDO !nv
                ENDDO !nhx
              ENDDO !nolist
              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '          NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',
     '          FILEFORMAT,OUTFILE,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,
     '          YQDATA,YQSDATA,ERROR,*9999)
            ENDDO !notime
          ENDIF

C***      Close the history file
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '      NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '      OUTFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

C *** DPN 28 April 2001 - Command to save/restore YP(ny,1)

        ELSE IF(TYPE(1:19).EQ.'CONVERGED_REFERENCE') THEN

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO no_nynr=1,NYNR(0,0,1,nr,nx)
C             loop over global ny variables of source region
              ny=NYNR(no_nynr,0,1,nr,nx) !global ny var number
              IF(OPERATION2(1:4).EQ.'SAVE') THEN
                YP(ny,11,nx) = YP(ny,1,nx)
              ELSE IF(OPERATION2(1:7).EQ.'RESTORE') THEN
                YP(ny,1,nx) = YP(ny,11,nx)
              ENDIF
            ENDDO !no_nynr (ny)
          ENDDO !no_nrlist=1,NRLIST(0)
          ! Update ZP and ZA
          CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr,nx),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,
     '      ERROR,*9999)

        ELSE !not COUPLED
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            IF(DEFFIBS) THEN
C             Initialise ZP,ZA
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nhx=1,NHP(np,nr,nx)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,1,nr) !nc=1
                    DO nk=1,NKH(nh,np,1,nr) !nc=1
                      ZP(nk,nv,nh,np,1)=0.0d0 !nc=1
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nh
              ENDDO !nonode (np)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                DO nhx=1,NH_LOC(0,nx)
                  nh=NH_LOC(nhx,nx)
                  DO na=1,NAT(NBH(nh,1,ne)) !nc=1
                    ZA(na,nh,1,ne)=0.0d0 !nc=1
                  ENDDO !na
                ENDDO !nh
              ENDDO !noelem (ne)
            ELSE
              CALL YPZP(IY,NBH,NEELEM,NHE,NHP(1,nr,nx),NKH(1,1,1,nr),
     '          NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,
     '          ERROR,*9999)
            ENDIF
            IF(TYPE(1:16).EQ.'CAVITY_REFERENCE') THEN
              IF(AVGENODE) THEN
C new MPN 17Apr97: average RC coords instead of curv. coords.
C               Calc RC coords of node to be changed
                DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                  nh=NH_LOC(nhx,nx)
                  Z(nhx)=ZP(1,1,nh,NP_average,1)
                ENDDO !nhx (nh)
                CALL XZ(ITYP11(nr),Z,ZRC_AVGE)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' NP_average='',I5,'' Z: '','
     '              //'3D12.4,'' ZRC_AVGE: '',3D12.4)') NP_average,
     '              (Z(nhx),nhx=1,NJ_LOC(NJL_GEOM,0,nr)),
     '              (ZRC_AVGE(nhx),nhx=1,NJ_LOC(NJL_GEOM,0,nr))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF !DOP
C               Calc the average nh_update RC coord for nodes in NPLIST
                SUM=0.0d0
                DO nonode=1,NPLIST(0)
                  np=NPLIST(nonode)
                  DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                    nh=NH_LOC(nhx,nx)
                    Z(nhx)=ZP(1,1,nh,np,1)
                  ENDDO !nhx (nh)
C                 Convert curv. coords to rc coords
                  CALL XZ(ITYP11(nr),Z,ZRC)
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,'('' np='',I5,'' ZRC: '',3D12.4)')
     '                np,(ZRC(nhx),nhx=1,NJ_LOC(NJL_GEOM,0,nr))
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF !DOP
                  SUM=SUM+ZRC(nh_update)
                ENDDO !nonode (np)
                ZRC_AVGE(nh_update)=SUM/DBLE(NPLIST(0))
C               Convert averaged RC coords to curv. coords
                CALL ZX(ITYP11(nr),ZRC_AVGE,Z)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' Avge  curv coords Z: '',3D12.4)')
     '              (Z(nhx),nhx=1,NJ_LOC(NJL_GEOM,0,nr))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF !DOP
C               Put averaged curv. coords into node NP_average
C               (unless there is more that one version)
                DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                  nh=NH_LOC(nhx,nx)
                  IF(NVHP(nh,NP_average,1,nr).EQ.1)
     '              ZP(1,1,nh,NP_average,1)=Z(nhx)
                ENDDO !nhx
C old MPN 17Apr97
CC               Calc the average nh_update value for nodes in NPLIST
C                SUM=0.0d0
C                DO nonode=1,NPLIST(0)
C                  np=NPLIST(nonode)
C                  SUM=SUM+ZP(1,1,nh_update,np,1)
C                ENDDO !np
C                ZP(1,1,nh_update,NP_average,1)=SUM/DBLE(NPLIST(0))
C end old
C               replace averaged ZP into YP(ny,1)
                CALL ZPYP(1,NBH,NEELEM,NHE,NHP(1,nr,nx),NKH(1,1,1,nr),
     '            NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,
     '            ERROR,*9999)
              ENDIF !AVGENODE
C             Copy YP(ny,1) into YP(ny,10) for FE50 Cavity problems
              CALL ZPYP(10,NBH,NEELEM,NHE,NHP(1,nr,nx),NKH(1,1,1,nr),
     '          NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,
     '          ERROR,*9999)
            ELSE IF(TYPE(1:13).EQ.'CORONARY_MESH') THEN
              CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr_host,nx),
     '          NKH(1,1,1,nr_host),
     '          NPNODE,nr_host,NVHP(1,1,1,nr_host),nx,NYNE,
     '          NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
              DO nonode=1,NPNODE(0,nr_coronary)
                np_coronary=NPNODE(nonode,nr_coronary)
                ne=NEP(np_coronary)
                CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),NKHE(1,1,1,ne),
     '            NPF(1,1),NPNE(1,1,ne),nr_host,NVHE(1,1,1,ne),
     '            NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '            ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
                DO nj=1,NJT
                  xi(nj)=XIP(nj,nonode)
                ENDDO
                CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '            NHE(ne,nx),nr_host,nx,DXIX,ZE,ZG,XI,ERROR,*9999)
              ENDDO
              DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                nh=NH_LOC(nhx,nx)
                ZP(1,1,nh,np_coronary,1)=ZG(1,nhx)
              ENDDO
              CALL ZPYP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr_coronary,nx),
     '          NKH(1,1,1,nr_coronary),
     '          NPNODE,nr_coronary,NVHP(1,1,1,nr_coronary),
     '          nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
            ELSE !other UPSOLU TYPEs
              nc=1 !only ever doing displacement (yet...)
C GMH 13/2/97 avoid unused error (should propagate through following)
              nc=nc
              NBMATCH=.TRUE.
              NKMATCH=.TRUE.
              NVMATCH=.TRUE.
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
C news MPN 11-Jul-96: extended to handle list of njj's for deffibs
                DO njj_num=1,NJJLIST(0)
                  njj=NJJLIST(njj_num)
                  IF(NHLIST(0).GT.0) THEN
                    nh_num=NHLIST(njj_num)
                    nh=NH_LOC(nh_num,nx)
                  ELSE IF(NHLIST(0).EQ.0) THEN
                    nh=NH_LOC(njj_num,nx)
                  ENDIF
                  IF(DEFFIBS) THEN
                    nj=NJ_LOC(NJL_FIEL,njj,nr)
                    DO nv=1,NVJP(nj,np)
                      DO nk=1,NKJ(nj,np)
                        ZP(nk,nv,nh,np,1)=XP(nk,nv,nj,np) !nc=1
                      ENDDO !nk
                    ENDDO !nv
                  ELSE
                    IF(TYPE(1:8).EQ.'GEOMETRY') THEN
                      nj=NJ_LOC(NJL_GEOM,njj,nr)
                    ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
                      nj=NJ_LOC(NJL_FIEL,njj,nr)
                    ENDIF
                    IF(NVJP(nj,np).NE.NVHP(nh,np,1,nr)) NVMATCH=.FALSE.
                    IF(NKJ(nj,np).NE.NKH(nh,np,1,nr)) NKMATCH=.FALSE.
                    DO nv=1,NVHP(nh,np,1,nr) !nc=1
                      DO nk=1,NKH(nh,np,1,nr) !nc=1
                        IF(OPERATION(1:10).EQ.'SUBSTITUTE') THEN
                          ZP(nk,nv,nh,np,1)=XP(nk,nv,nj,np)
                        ELSE IF(OPERATION(1:3).EQ.'ADD') THEN
                          ZP(nk,nv,nh,np,1)=
     '                      ZP(nk,nv,nh,np,1)+XP(nk,nv,nj,np)
                        ELSE IF(OPERATION(1:8).EQ.'SUBTRACT') THEN
                          ZP(nk,nv,nh,np,1)=
     '                      ZP(nk,nv,nh,np,1)-XP(nk,nv,nj,np)
                        ENDIF
                      ENDDO !nk
                    ENDDO !nv
                  ENDIF
                ENDDO !njj_num
C old
C                DO nhx=1,NHP(np,nr,nx)
C                  nh=NH_LOC(nhx,nx)
C                  IF(TYPE(1:8).EQ.'GEOMETRY') THEN
C                    nj=NJ_LOC(NJL_GEOM,nhx,nr)
C                  ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
C                    nj=NJ_LOC(NJL_FIEL,nhx,nr)
C                  ENDIF
C                  DO nv=1,NVHP(nh,np,1,nr) !nc=1
C                    DO nk=1,NKH(nh,np,1,nr) !nc=1
C                      IF(OPERATION(1:10).EQ.'SUBSTITUTE') THEN
C                        ZP(nk,nv,nh,np,1)=XP(nk,nv,nj,np)
C                      ELSE IF(OPERATION(1:3).EQ.'ADD') THEN
C                        ZP(nk,nv,nh,np,1)=ZP(nk,nv,nh,np,1)+
C     '                    XP(nk,nv,nj,np)
C                      ELSE IF(OPERATION(1:8).EQ.'SUBTRACT') THEN
C                        ZP(nk,nv,nh,np,1)=ZP(nk,nv,nh,np,1)-
C     '                    XP(nk,nv,nj,np)
C                      ENDIF
C                    ENDDO !nk
C                  ENDDO !nv
C                ENDDO !nh
              ENDDO !nonode (np)
              CALL ZPYP(IY,NBH,NEELEM,NHE,NHP(1,nr,nx),NKH(1,1,1,nr),
     '          NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '          ERROR,*9999)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                nb=NBH(nh,nc,ne)
                IF(NBJ(nj,ne).NE.nb) NBMATCH=.FALSE.
                DO nn=1,NNT(nb)
                  IF(NVJE(nn,nb,nj,ne).NE.NVHE(nn,nb,nh,ne))
     '              NVMATCH=.FALSE.
                ENDDO !nn
              ENDDO !noelem
              IF(.NOT.NBMATCH) THEN
                WRITE(OP_STRING,
     '            '('' >>WARNING: Bases not consistent'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ENDIF
              IF(.NOT.NVMATCH) THEN
                WRITE(OP_STRING,
     '            '('' >>WARNING: Versions not consistent'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ENDIF
              IF(.NOT.NKMATCH) THEN
                WRITE(OP_STRING,
     '            '('' >>WARNING: Derivatives not consistent'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF !TYPE.EQ.'CAVITY_REFERENCE'
          ENDDO !nr
        ENDIF
      ENDIF

      CALL EXITS('UPSOLU')
      RETURN
 9999 CALL ERRORS('UPSOLU',ERROR)
      CALL EXITS('UPSOLU')
      RETURN 1
      END


