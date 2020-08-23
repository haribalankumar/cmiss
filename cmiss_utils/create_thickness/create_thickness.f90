! -*-f90-*-
PROGRAM Create_thickness

  ! Created: Geoff Carden 4/99
  ! Modified: GPC 1/02
  ! Compile: f90 -n32 -c create_thickness.o -g create_thickness.f90
  ! Link: f90 -n32 -o create_thickness create_thickness.o constants.o

  ! this program reads in existing node and element files for a
  ! membrane of a bone and assigns a thickness to the bone. The
  ! program asks for the name or the files to be altered and adds the
  ! extra nodes and changes the element file to suit the thickness
  ! measurement.

  ! The program ultimately writes out the new ipnode and ipelem files
  ! so the 3_D basis functions need to be st up in CMISS prior to
  ! reading in the new files. Linear and cubic hermite basis functions
  ! (without versions) can be handled

  ! The program also generates .ipslic and .ipthic files which
  ! specify the nodes in each node slice and the thickness of the
  ! bone at each node.

  USE constants
  USE fileio
  USE thickness
  USE slicefile
  USE ipelemfile
  USE ipnodefile
  USE createnodes
  IMPLICIT NONE

  ! Node data
  INTEGER(I4)::NumNode
  INTEGER(I4)::NumCoord=3
  INTEGER(I4),ALLOCATABLE::NumDeriv(:)
  INTEGER(I4)::MaxNumDeriv
  INTEGER(I4)::MaxNumVersions
  INTEGER(I4),ALLOCATABLE::NodeVersion(:,:)
  INTEGER(I4),ALLOCATABLE::NodeIndex(:)
  INTEGER(I4),ALLOCATABLE::NodeSliceCount(:)
  REAL(dp),ALLOCATABLE::NodeValue(:,:,:,:)
  ! NodeValue(:,:,:,:)
  !   The first  dimension is the number of nodes
  !   The second dimension is the number of coordinates
  !   The third  dimension is the number of versions
  !   The fourth dimension is the position and the value of the derivatives 

  ! New data
  INTEGER(I4)::NumNode_new
  INTEGER(I4)::NumNodeVers
  INTEGER(I4)::MaxNewNode
  INTEGER(I4)::MaxNumNodeVers
  INTEGER(I4),ALLOCATABLE::NodeVersion_new(:,:)
  INTEGER(I4),ALLOCATABLE::NodeIndex_new(:)
  REAL(dp),ALLOCATABLE::NodeValue_new(:,:,:,:)


  ! Element data
  INTEGER(I4)::NumElem                            ! How many of the buggers we have
  INTEGER(I4),ALLOCATABLE::ElemIndex(:)           ! Index # of each element
  INTEGER(I4),ALLOCATABLE::ElemNode(:,:)          ! List of nodes in each element
  INTEGER(I4),ALLOCATABLE::ElemNode_new(:,:)      ! New values created
  INTEGER(I4),ALLOCATABLE::ElemVersion(:,:,:)     ! Versions of element nodes
  INTEGER(I4),ALLOCATABLE::ElemVersion_new(:,:,:) ! New versions numbers


  ! Slice data
  INTEGER(I4)::NumSlice                           !
  INTEGER(I4)::MaxNodeInSlice                     ! Maximum of nodes in any slice
  INTEGER(I4),ALLOCATABLE::SliceIndex(:)          ! Index # of slice
  INTEGER(I4),ALLOCATABLE::SliceNodes(:,:)        ! List of nodes in slice
  INTEGER(I4),ALLOCATABLE::NumNodeInSlice(:)      ! Number in each slice


  ! Redistribution data
  INTEGER(I4)::NumRedist                          ! Number of redistribution splines
  CHARACTER(LEN=FILE_LEN)::Redist(MaxNumRedist)   ! Character strings of redist data
  LOGICAL::do_redist                              ! Flag to do redist


  ! Projection data
  REAL(dp),ALLOCATABLE::ThickArray(:)
  ! Stores the thickness at each of the original nodes
  INTEGER(I4)::Node_Offset
  ! Set offset for the new node numbers. The default is 0
  INTEGER(I4),ALLOCATABLE::NodeProjectIndex(:,:)


  ! Filenames
  CHARACTER(LEN=FILE_LEN)::DefaultName
  CHARACTER(LEN=FILE_LEN)::SliceFileName
  CHARACTER(LEN=FILE_LEN)::ThickFileName
  CHARACTER(LEN=FILE_LEN)::NodeFileName
  CHARACTER(LEN=FILE_LEN)::ElemFileName
  CHARACTER(LEN=FILE_LEN_3D)::OutputFileName
  CHARACTER(LEN=FULL_FILE_LEN)::FullSliceFileName
  CHARACTER(LEN=FULL_FILE_LEN)::FullThickFileName
  CHARACTER(LEN=FULL_FILE_LEN)::FullNodeFileName
  CHARACTER(LEN=FULL_FILE_LEN)::FullElemFileName
  CHARACTER(LEN=FULL_FILE_LEN_3D)::FullNodeFileName3D
  CHARACTER(LEN=FULL_FILE_LEN_3D)::FullElemFileName3D


  ! IO unit numbers
  INTEGER(I4),PARAMETER::UnitSF=20   ! Slice File
  INTEGER(I4),PARAMETER::UnitIP=22   ! IPNode file
  INTEGER(I4),PARAMETER::UnitPM=40   !
  INTEGER(I4),PARAMETER::UnitTF=50   ! Thickness file
  INTEGER(I4),PARAMETER::UnitElem=30 ! IPElem file


  ! Local variables
  REAL(dp)::Bone_Thickness
  INTEGER(I4)::i,j,k,value_read,int_read,def_int,ierr,len,itmp
  CHARACTER(LEN=3)::filetype,slice_read,thick_read,output_read

  REAL(DP)::fudge_alpha
  REAL(DP)::ratio_value
  INTEGER(I4):: new_basis
  INTEGER(I4):: node_offset
  LOGICAL:: prompt_slice
  LOGICAL:: prompt_thick
  LOGICAL:: read_thick
  LOGICAL:: ratio_thick
  LOGICAL:: one_node
  LOGICAL:: new_versions
  LOGICAL:: verbose

  ! External functions
  INTEGER, EXTERNAL::iargc


  ! Set defaults
  DefaultName='filename'
  SliceFileName=''
  ThickFileName=''
  NodeFileName=''
  ElemFileName=''
  OutputFileName=''

  Redist=''
  NumRedist=0
  do_redist=.false.

  verbose=.false.
  one_node=.true.
  new_versions=.true.
  ratio_thick=.true.
  ratio_value=1.0
  prompt_slice=.false.
  prompt_thick=.false.
  read_thick=.false.
  node_offset=0
  new_basis=0
  fudge_alpha=1.0


  ! Parse command line
  IF (iargc().gt.0) THEN
    i=1
    DO WHILE (i.le.iargc())
      call getarg(i,buff)

      ! Filenames
      IF (buff(1:5).eq.'-file') THEN
        i=i+1
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,DefaultName)
      ELSE IF (buff(1:6).eq.'-slice') THEN
        i=i+1
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,SliceFileName)
      ELSE IF (buff(1:6).eq.'-thick') THEN
        i=i+1
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,ThickFileName)
      ELSE IF (buff(1:5).eq.'-node') THEN
        i=i+1
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,NodeFileName)
      ELSE IF (buff(1:5).eq.'-elem') THEN
        i=i+1
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,ElemFileName)
      ELSE IF (buff(1:7).eq.'-output') THEN
        i=i+1
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,OutputFileName)

      ! Operations
      ELSE IF (buff.eq.'-prompt_slice') THEN
        prompt_slice=.true.
      ELSE IF (buff.eq.'-read_slice' .or. buff.eq.'-prompt_slice...not') THEN
        prompt_slice=.false.

      ELSE IF (buff.eq.'-prompt_thick') THEN
        prompt_thick=.true.
      ELSE IF (buff.eq.'-read_thick' .or. buff.eq.'-prompt_thick...not') THEN
        prompt_thick=.false.

      ELSE IF (buff.eq.'-no_new_versions') THEN
        new_versions=.false.
      ELSE IF (buff.eq.'-new_versions') THEN
        new_versions=.true.

      ELSE IF (buff.eq.'-one_node...not') THEN
        one_node=.false.
      ELSE IF (buff.eq.'-one_node') THEN
        one_node=.true.

      ELSE IF (buff.eq.'-ratio_thick...not') THEN
        ratio_thick=.false.
      ELSE IF (buff.eq.'-ratio_thick') THEN
        ratio_thick=.true.

      ELSE IF (buff.eq.'-add_basis') THEN
        i=i+1
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,buff)
        READ(buff,*,iostat=ierr) new_basis
        IF (ierr.ne.0) CALL Usage()

      ELSE IF (buff.eq.'-offset') THEN
        i=i+1
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,buff)
        READ(buff,*,iostat=ierr) node_offset
        IF (ierr.ne.0) CALL Usage()

      ELSE IF (buff.eq.'-fudge') THEN
        i=i+1
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,buff)
        READ(buff,*,iostat=ierr) fudge_alpha
        IF (ierr.ne.0) CALL Usage()

      ELSE IF (buff.eq.'-ratio') THEN
        i=i+1
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,buff)
        READ(buff,*,iostat=ierr) ratio_value
        IF (ierr.ne.0) CALL Usage()
        prompt_thick=.false.
        read_thick=.false.

      ELSE IF (buff.eq.'-redist') THEN
        i=i+1
        NumRedist=NumRedist+1
        IF (NumRedist.gt.MaxNumRedist) THEN
          WRITE(stderr,*) 'create_thickness: NumRedist > MaxNumRedist (',&
            & NumRedist,' > ',MaxNumRedist,')'
          WRITE(stderr,*) '  Could be a good idea to increase&
            & MaxNumRedist in constants.f90'
          STOP
        ENDIF
        IF (i.gt.iargc()) CALL Usage()
        call getarg(i,Redist(NumRedist))
        do_redist=.true.

      ! Other args
      ELSE IF (buff.eq.'-verb') THEN
        verbose=.true.
      ELSE IF (buff.eq.'-help') THEN
        CALL Usage()
        STOP
      ELSE
        CALL Usage()
        STOP
      ENDIF

      i=i+1
    ENDDO

    ! Set filenames
    IF (SliceFileName.eq.'')  SliceFileName=DefaultName
    IF (ThickFileName.eq.'')  ThickFileName=DefaultName
    IF (NodeFileName.eq.'')   NodeFileName=DefaultName
    IF (ElemFileName.eq.'')   ElemFileName=DefaultName
    IF (OutputFileName.eq.'') OutputFileName=TRIM(DefaultName)//'3D'
    CALL CleanFileName( FullSliceFileName,SliceFileName,'.ipslic' )
    CALL CleanFileName( FullThickFileName,ThickFileName,'.ipthic' )
    CALL CleanFileName( FullNodeFileName,NodeFileName,'.ipnode' )
    CALL CleanFileName( FullElemFileName,ElemFileName,'.ipelem' )
    CALL CleanFileName( FullNodeFileName3D,OutputFileName,'.ipnode' )
    CALL CleanFileName( FullElemFileName3D,OutputFileName,'.ipelem' )


  ! Get the info from the user the old painful way
  ELSE
    WRITE(stdout,'(A)',advance='no') ' > Read existing slice file <y,n> [n]? '
    READ(stdin,'(A)') buff
    buff=ADJUSTL(buff)
    IF (buff.eq.'y' .or. buff.eq.'Y') prompt_slice=.false.

    ! Get the default file name/slice file name
    WRITE(stdout,'(A)',advance='no') ' Slice file name: '
    READ(stdin,'(A)') buff
    buff=ADJUSTL(buff)
    IF (buff.ne.' ') THEN
      SliceFileName=buff
    ELSE
      SliceFileName=DefaultName
    ENDIF
    CALL CleanFileName( FullSliceFileName,SliceFileName,'.ipslic' )
    DefaultName=SliceFileName

    ! Get the ipnode file name
    WRITE(stdout,'(3A)',advance='no') ' Input node file name [',&
      & TRIM(DefaultName),']: '
    READ(stdin,'(A)') buff
    buff=ADJUSTL(buff)
    IF (buff.ne.' ') THEN
      NodeFileName=buff
    ELSE
      NodeFileName=DefaultName
    ENDIF
    CALL CleanFileName( FullNodeFileName,NodeFileName,'.ipnode' )
    DefaultName=NodeFileName

    ! Get the ipelem file name
    WRITE(stdout,'(3A)',advance='no') ' Input element file name [',&
      & TRIM(DefaultName),']: '
    READ(stdin,'(A)') buff
    buff=ADJUSTL(buff)
    IF (buff.ne.' ') THEN
      ElemFileName=buff
    ELSE
      ElemFileName=DefaultName
    ENDIF
    CALL CleanFileName( FullElemFileName,ElemFileName,'.ipelem' )
    DefaultName=ElemFileName

    ! Get the output file names
    OutputFileName=TRIM(DefaultName)//'3D'
    WRITE(stdout,'(4A)',advance='no') ' Enter the name for the 3D',&
      & ' Node and Element Output files [',TRIM(OutputFileName),']: '
    READ(stdin,'(A)') buff
    buff=ADJUSTL(buff)
    IF (buff.ne.' ') THEN
      OutputFileName=buff
    ELSE
      OutputFileName=TRIM(DefaultName)//'3D'
    ENDIF
    CALL CleanFileName( FullNodeFileName3D,OutputFileName,'.ipnode' )
    CALL CleanFileName( FullElemFileName3D,OutputFileName,'.ipelem' )

    ! Get node offset
    WRITE(stdout,'(2A)',advance='no') ' Set the node offset for the ',&
      & 'new nodes [0]: '
    READ(stdin,'(A)') buff
    IF (buff.ne.'') THEN
      READ(buff,*,iostat=ierr) itmp
      IF (ierr.eq.0) node_offset=itmp
    ENDIF

    ! Get thickness settings
    one_node=.true.
    WRITE(stdout,'(A)',advance='no') ' > Only add one node per slice <y,n> [y]? '
    READ(stdin,'(A)') buff
    buff=ADJUSTL(buff)
    IF (buff.eq.'n' .or. buff.eq.'N') one_node=.false.

    ! Thickness data is only needed if we are processing slice by slice
    IF (.not.one_node) THEN
      prompt_thick=.false.
      WRITE(stdout,'(A)',advance='no') ' > Read existing thickness file <y,n> [n]? '
      READ(stdin,'(A)') buff
      buff=ADJUSTL(buff)
      IF (buff.ne.'y' .and. buff.ne.'Y') prompt_thick=.true.

      WRITE(stdout,'(3A)',advance='no') ' Enter the thickness file name [',&
        & TRIM(DefaultName),']: '
      READ(stdin,'(A)') buff
      buff=ADJUSTL(buff)
      IF (buff.ne.' ') THEN
        ThickFileName=buff
      ELSE
        ThickFileName=DefaultName
      ENDIF
      CALL CleanFileName( FullThickFileName,ThickFileName,'.ipthic' )
      DefaultName=ThickFileName
    ENDIF
  ENDIF


  !------------------------------------------------
  ! Run the program
  !------------------------------------------------


  !------------------------------------------------
  ! Read in the slice data
  !------------------------------------------------
  IF (prompt_slice) THEN
    CALL PromptSliceDim( NumSlice,MaxNodeInSlice,verbose )
  ELSE
    WRITE(stdout,'(3A)') 'Reading in Slice File ''',&
      & TRIM(FullSliceFileName),''''
    CALL ReadSliceDim( NumSlice,MaxNodeInSlice,UnitSF, &
      & FullSliceFileName,verbose )
    WRITE(stdout,'(2(A,I4))') ' Slice dimensions: ',&
      & NumSlice,' x ',MaxNodeInSlice
  ENDIF

  ALLOCATE(SliceIndex(NumSlice),                &
    &      NumNodeInSlice(NumSlice),            &
    &      SliceNodes(NumSlice,MaxNodeInSlice), &
    &      stat=ierr)
  IF (ierr.ne.0) THEN
    WRITE(stderr,*) 'create_thickness: error allocating memory for&
      & slice data'
    STOP
  ENDIF

  IF (prompt_slice) THEN
    CALL PromptSliceFile( NumSlice,MaxNodeInSlice,SliceIndex,&
      & NumNodeInSlice,SliceNodes,verbose )

    CALL WriteSliceFile( NumSlice,SliceIndex,NumNodeInSlice,&
      & SliceNodes,UnitSF,FullSliceFileName,verbose )
  ELSE
    CALL ReadSliceFile( NumSlice,MaxNodeInSlice,SliceIndex,&
      & NumNodeInSlice,SliceNodes,UnitSF,FullSliceFileName,verbose )
  ENDIF


  !------------------------------------------------
  ! Read in the ipnode data.
  !------------------------------------------------
  WRITE(stdout,'(3A)') 'Reading in Node File ''',TRIM(FullNodeFileName),''''
  CALL ReadIPNodeDim( NumNode,NumCoord,MaxNumDeriv,MaxNumVersions, &
    & UnitIP,FullNodeFileName,verbose )
  MaxNumVersions=MAX(3*MaxNumVersions,32)

  ALLOCATE(NumDeriv(NumCoord), &
         & NodeIndex(NumNode), &
         & NodeSliceCount(NumNode), &
         & NodeVersion(NumNode,NumCoord), &
         & NodeValue(NumNode,NumCoord,MaxNumVersions,MaxNumDeriv+1), &
         & stat=ierr)
  IF (ierr.ne.0) THEN
    WRITE(stderr,*) 'create_thickness: Error allocating memory for&
      & reading in node data'
    STOP
  ENDIF

  CALL ReadIPNodeFile( NumNode,NumCoord,MaxNumDeriv,MaxNumVersions, &
    & NumDeriv,NodeIndex,NodeVersion,NodeValue, &
    & UnitIP,FullNodeFileName,verbose )
  WRITE(stdout,'(A,I4)') ' Number of nodes: ',NumNode

  ! Count the number of slices that each node appears in.
  CALL CountMentionsOfNodesInSlices( NumNode ,NumSlice, &
    & NodeIndex, NodeSliceCount, NumNodeInSlice, SliceNodes, verbose  )


  ! -----------------------------------------------
  ! Create the thickness data (from file or stdin)
  ! -----------------------------------------------
  ALLOCATE(ThickArray(NumNode))
  IF (.not.one_node) THEN
    ! Read thickness values from prompting
    IF (prompt_thick) THEN
      OPEN(unit=unittf,file=FullThickFilename,status='replace')    
      CALL PromptThickness(ThickArray,NumNode,NodeIndex,&
        & NumNodeInSlice,SliceNodes,NumSlice,ratio_thick,unittf)
      CLOSE(unittf)
    ! Read values from a thickness file
    ELSE IF (read_thick) THEN
      OPEN(unit=unittf,file=FullThickFilename,status='old',iostat&
        & =ierr)
      IF(ierr.NE.0) THEN
        WRITE(stderr,'(3a)') 'Cannot open "',TRIM(FullThickFilename)&
          & ,'"'
        STOP
      ENDIF
      CALL ReadThickness(ThickArray,NumNode,NodeIndex,ratio_thick&
        & ,unittf)
      WRITE(stdout,'(">> Thickness File read")')
      WRITE(stdout,'()')
      CLOSE(unittf)
    ! Use values from command line
    ELSE
      ThickArray(:)=ratio_value
    ENDIF
  ENDIF


  !------------------------------------------------
  ! Read in the element data.
  !------------------------------------------------
  WRITE(stdout,'(3A)') 'Reading in Element File ''',&
    & TRIM(FullElemFileName),''''
  CALL ReadElemDim( NumElem,UnitElem,FullElemFileName,verbose )
  ALLOCATE(ElemIndex(NumElem),                                  &
    &      ElemNode(NumElem,NumNodePerElem_2D),                 &
    &      ElemVersion(NumElem,NumNodePerElem_2D,4),            &
    &      stat=ierr)
  IF (ierr.ne.0) THEN
    WRITE(stderr,*) 'create_thickness: Error allocating element memory'
    STOP
  ENDIF
  CALL ReadElemFile( NumElem,NumNode,NumCoord,NumNodePerElem_2D,&
    & ElemIndex,ElemNode,ElemVersion,NodeIndex,NodeVersion, &
    & UnitElem,FullElemFileName,Verbose )
  WRITE(stdout,'(A,I4)') ' Number of elements: ',NumElem


  !------------------------------------------------------------------
  ! Allocate space for new nodes and elements
  !------------------------------------------------------------------
  MaxNewNode=3*NumNode      ! These values are arbitary. If we run
  MaxNumNodeVers=32*NumNode ! into trouble, increase the values

  ALLOCATE(NodeValue_New(MaxNewNode,NumCoord,MaxNumVersions,MaxNumDeriv+1),&
    &      NodeIndex_New(MaxNewNode),                    &
    &      NodeVersion_New(MaxNewNode,NumCoord),         &
    &      NodeProjectIndex(MaxNumNodeVers,4),           &
    &      ElemNode_new(NumElem,NumNodePerElem_3D),      &
    &      ElemVersion_new(NumElem,NumNodePerElem_3D,NumCoord+1), &
    &      stat=ierr)
  IF (ierr.ne.0) THEN
    WRITE(stderr,*) 'create_thickness: Error allocating memory for&
      & generation of new nodes'
    STOP
  ENDIF

  !------------------------------------------------------------------
  ! Creat the new nodes
  !------------------------------------------------------------------
  WRITE(stdout,'(A)') 'Creating New Nodes'
  CALL NewNodes( NumNode,NumCoord,MaxNumVersions,MaxNodeInSlice,  &
    & NumDeriv,Node_Offset,NumNode_New,MaxNewNode,                &
    & NodeValue,NodeIndex,NodeVersion,                            &
    & NodeValue_New,NodeIndex_New,NodeVersion_New,                &
    & NumSlice,NumNodeInSlice,SliceNodes, ThickArray,             &
    & NodeProjectIndex, NumNodeVers, MaxNumNodeVers,              &
    & Fudge_alpha,Ratio_thick,One_node,new_versions,verbose )
  WRITE(stdout,'(A,I4)') ' Number of nodes: ',NumNode_New

  !------------------------------------------------------------------
  ! Creat the new elements
  !------------------------------------------------------------------
  WRITE(stdout,'(A)') 'Creating New Elements'
  CALL NewElems( NumElem,NumNode_new,NumCoord,NumNodeVers, &
    & NumNodePerElem_2D, NumNodePerElem_3D,                &
    & ElemNode,ElemVersion,ElemNode_new,ElemVersion_new,   &
    & NumSlice,NodeSliceCount,NumNodeInSlice,SliceNodes,   &
    & NodeProjectIndex,NodeIndex_New,NodeVersion_New,      &
    & verbose )
  WRITE(stdout,'(A,I4)') ' Number of elements: ',NumElem

  !------------------------------------------------------------------
  ! Redistribute points
  !------------------------------------------------------------------
  IF (do_redist) THEN
    CALL RedistPoints( NumRedist,Redist, &
    & NumNode_new,NumCoord,NodeIndex_new,NodeVersion_new,NodeValue_new )
  ENDIF

  !------------------------------------------------------------------
  ! Write out the new 3D node file
  !------------------------------------------------------------------
  WRITE(stdout,'(3A)') 'Writing out new Node File ''',&
    & TRIM(FullNodeFileName3D),''''
  CALL WriteIPNodeFile(NumNode_New,NumCoord,NumDeriv,&
    & NodeIndex_New,NodeVersion_New,NodeValue_New,UnitIP,&
    & FullNodeFileName3D,OutputFileName,verbose)

  !------------------------------------------------------------------
  ! Write out the new 3D elem file
  !------------------------------------------------------------------
  WRITE(stdout,'(3A)') 'Writing out new Element File ''',&
    & TRIM(FullElemFileName3D),''''
  CALL WriteElemFile( NumElem,NumCoord,NumNodePerElem_3D, &
    & ElemIndex,ElemNode_new,ElemVersion_new,new_basis,   &
    & UnitElem,FullElemFileName3D,OutputFileName,verbose )


  !------------------------------------------------------------------
  ! Clean up -- no real need to do this
  !------------------------------------------------------------------
  DEALLOCATE(ElemIndex, &
     & ElemNode,        &
     & ElemNode_new,    &
     & NodeIndex,       &
     & NodeSliceCount,  &
     & NodeValue,       &
     & NodeVersion,     &
     & NodeIndex_new,   &
     & NodeValue_new,   &
     & NodeVersion_new, &
     & ElemVersion,     &
     & ElemVersion_new, &
     & ThickArray,      &
     & stat=ierr)

  STOP
END PROGRAM Create_thickness
