! -*-f90-*-
MODULE slicefile

	USE constants
	USE fileio
	IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------

  SUBROUTINE PromptSliceDim( NumSlice,MaxNodeInSlice,verbose )

    ! Prompt the user for the dimensions of the slice information

    IMPLICIT NONE
    ! Parameters
    INTEGER(I4)::NumSlice
    INTEGER(I4)::MaxNodeInSlice
    LOGICAL::verbose
    ! Local variables
    INTEGER(I4)::ierr


    WRITE(stdout,'(A)',advance='no') ' Number of Slices: '
    READ(stdin,*,iostat=ierr) NumSlice
    IF (ierr.ne.0) THEN
      WRITE(stderr,*) 'PromptSliceDim: error reading dimensions'
      STOP
    ENDIF

    MaxNodeInSlice=MaximumNodesInSlice

    RETURN
  END SUBROUTINE PromptSliceDim

  !----------------------------------------------------------------

  SUBROUTINE PromptSliceFile( NumSlice,MaxNodeInSlice,SliceIndex,&
    & NumNodeInSlice,SliceNodes,verbose )

    ! Prompts the nodal slice file

    IMPLICIT NONE
    ! Subroutine args
    INTEGER(I4)::NumSlice
    INTEGER(I4)::MaxNodeInSlice
    INTEGER(I4)::SliceIndex(:)
    INTEGER(I4)::NumNodeInSlice(:)
    INTEGER(I4)::SliceNodes(:,:)
    LOGICAL::verbose
    ! Local variables
    INTEGER(I4)::i,j,ierr,nlow,ntoken,int_read,def_num,counter
    INTEGER(I4)::NumNodeInSliceFile


    SliceIndex(:)=0
    NumNodeInSlice(:)=0
    SliceNodes(:,:)=0
    counter=1
    def_num=1

    NumNodeInSliceFile=0
    DO i=1,NumSlice
      WRITE(stdout,*)

      WRITE(stdout,'(A,I5,A)',advance='no') ' Slice number [',counter,']: '
      READ(stdin,*) int_read
      SliceIndex(i)=counter
      IF (int_read.ne.0) SliceIndex(i)=int_read
      counter=SliceIndex(i)+1

      WRITE(stdout,'(A,I5,A)',advance='no') ' Number of nodes for the&
        & slice [',def_num,']: '
      READ(stdin,*) int_read
      NumNodeInSlice(i)=def_num
      IF (int_read.ne.0) THEN
        NumNodeInSlice(i)=int_read
        def_num=int_read
      ENDIF

      IF (NumNodeInSlice(i).gt.MaxNodeInSlice) THEN
        WRITE(stderr,*) 'PromptSliceFile: too many nodes in slice;',&
          & NumNodeInSlice(i),'>',MaxNodeInSlice
        STOP
      ENDIF

      ! Read node numbers for slice
      nlow=0
      WRITE(stdout,'(A)',advance='no') ' Slice Nodes: '
      DO WHILE (nlow.lt.NumNodeInSlice(i))
        READ(stdin,'(A)') buff
        ntoken=StrTokCnt(buff)
        READ(buff,*,iostat=ierr) (SliceNodes(i,j),j=nlow+1,nlow+ntoken)
        IF (ierr.ne.0) THEN
          WRITE(stderr,*) 'PromptSliceFile: error reading line'
          STOP
        ENDIF
        nlow=nlow+ntoken
      ENDDO

      NumNodeInSliceFile=NumNodeInSliceFile + NumNodeInSlice(i)
    END DO
    WRITE(stdout,'(A,I6,A)') 'PromptSliceFile: Read  (',&
      & NumNodeInSliceFile,' nodes)'
    WRITE(stdout,'()') 

    RETURN
  END SUBROUTINE PromptSliceFile

  !-----------------------------------------------------------------

	SUBROUTINE ReadSliceDim( NumSlice,MaxNodeInSlice,Unit,FileName,verbose )

    ! Read the slice file to get the number of slices and the maximum
    ! number of nodes per slice

    IMPLICIT NONE
    ! Parameters
    INTEGER(I4)::NumSlice
    INTEGER(I4)::MaxNodeInSlice
    INTEGER(I4)::Unit
    CHARACTER(LEN=*)::FileName
    LOGICAL::verbose
    ! Local variables
    INTEGER(I4)::i,k,ierr,nlow,ntoken
    INTEGER(I4)::NumNodesInSlice


    OPEN(unit=Unit,file=FileName,status='old',iostat=ierr)
	  IF (ierr.ne.0) THEN
      WRITE(stderr,'(3A)') 'ReadSliceDim: Can''t open file "',TRIM(FileName),'"'
      STOP
	  ENDIF
    CALL GetLine( Unit,buff,verbose,ierr )
    READ(buff(scan(buff,':')+1:),*) NumSlice

    MaxNodeInSlice=0
    DO i=1,NumSlice
      ! Slice Number
      CALL GetLine( Unit,buff,verbose,ierr )

      ! Number of nodes in slice
      CALL GetLine( Unit,buff,verbose,ierr )
      READ(buff(scan(buff,':')+1:),*) NumNodesInSlice
      MaxNodeInSlice=MAX(MaxNodeInSlice,NumNodesInSlice)
      IF (verbose) WRITE(stdout,*) NumNodesInSlice,MaxNodeInSlice

      ! Nodes
      nlow=0
      DO WHILE (nlow.lt.NumNodesInSlice)
        CALL GetLine( Unit,buff,verbose,ierr )
        buff=buff(scan(buff,':')+1:)
        ntoken=StrTokCnt(buff)
        nlow=nlow+ntoken
      ENDDO

    END DO

    CLOSE(unit=Unit)
    RETURN
  END SUBROUTINE ReadSliceDim

  !-----------------------------------------------------------------

	SUBROUTINE ReadSliceFile( NumSlice,MaxNodeInSlice,SliceIndex,&
    & NumNodeInSlice,SliceNodes,Unit,FileName,verbose )

    ! Read a node slice file. ReadSliceDim() should previously have
    ! been called to get the dimensions of the arrays.

    IMPLICIT NONE
    ! Parameters
    INTEGER(I4)::NumSlice
    INTEGER(I4)::MaxNodeInSlice
    INTEGER(I4)::SliceIndex(:)
    INTEGER(I4)::NumNodeInSlice(:)
    INTEGER(I4)::SliceNodes(:,:)
    INTEGER(I4)::Unit
    CHARACTER(LEN=*)::FileName
    LOGICAL::verbose
    ! Local variables
    INTEGER(I4)::i,j,k,num,ierr,nlow,ntoken
    INTEGER(I4)::NumNodeInSliceFile


    OPEN(unit=Unit,file=FileName,status='old',iostat=ierr)
	  IF (ierr.ne.0) THEN
      WRITE(stderr,*) 'ReadSliceFile: Can''t open file "',TRIM(FileName),'"'
      STOP
	  ENDIF
    CALL GetLine( Unit,buff,verbose,ierr )
    READ(buff(scan(buff,':')+1:),*) num
    IF (num.ne.NumSlice) THEN
      WRITE(stderr,*) 'ReadSliceFile: NumSlice dimenions missmatch:',&
        & num,NumSlice
      STOP
    ENDIF


    SliceIndex(:)=0
    NumNodeInSlice(:)=0
    SliceNodes(:,:)=0
    NumNodeInSliceFile=0

    DO i=1,NumSlice
      CALL GetLine( Unit,buff,verbose,ierr )
      READ(buff(scan(buff,':')+1:),*) SliceIndex(i) 
      IF (verbose) WRITE(stdout,*) SliceIndex(i)

      CALL GetLine( Unit,buff,verbose,ierr )
      READ(buff(scan(buff,':')+1:),*) NumNodeInSlice(i)
      IF (verbose) WRITE(stdout,*) NumNodeInSlice(i)

      IF (NumNodeInSlice(i).gt.MaxNodeInSlice) THEN
        WRITE(stderr,*) 'ReadSliceFile: too many nodes in slice;',&
          & NumNodeInSlice(i),'>',MaxNodeInSlice
        STOP
      ENDIF

      nlow=0
      DO WHILE (nlow.lt.NumNodeInSlice(i))
        CALL GetLine( Unit,buff,verbose,ierr )
        buff=buff(scan(buff,':')+1:)
        ntoken=StrTokCnt(buff)
        READ(buff,*,iostat=ierr) (SliceNodes(i,j),j=nlow+1,nlow+ntoken)
        IF (ierr.ne.0) THEN
          WRITE(stderr,*) 'ReadSliceFile: error reading slice file'
          STOP
        ENDIF
        nlow=nlow+ntoken
      ENDDO
      IF (verbose) WRITE(stdout,*) (SliceNodes(i,j),j=1,&
        & NumNodeInSlice(i))

      NumNodeInSliceFile=NumNodeInSliceFile + NumNodeInSlice(i)
    END DO

    IF (verbose) THEN
      WRITE(stdout,'(A,I6,A)') 'ReadSliceFile: Slice file read  (',&
        & NumNodeInSliceFile,' nodes)'
      WRITE(stdout,*)
    ENDIF

    CLOSE(unit=Unit)
    RETURN
  END SUBROUTINE ReadSliceFile

  !-----------------------------------------------------------------

	SUBROUTINE WriteSliceFile( NumSlice,SliceIndex,NumNodeInSlice,&
    & SliceNodes,Unit,FileName,verbose )

    ! Write a node slice file.

    IMPLICIT NONE
    ! Parameters
    INTEGER(I4)::NumSlice
    INTEGER(I4)::SliceIndex(:)
    INTEGER(I4)::NumNodeInSlice(:)
    INTEGER(I4)::SliceNodes(:,:)
    LOGICAL::old_style_slice_files
    INTEGER(I4)::Unit
    CHARACTER(LEN=*)::FileName
    LOGICAL::verbose
    ! Local variables
    INTEGER(I4)::i,j,ierr
    INTEGER(I4)::NumNodeInSliceFile


    OPEN(unit=Unit,file=FileName,status='unknown',iostat=ierr)
	  IF (ierr.ne.0) THEN
      WRITE(stderr,*) 'WriteSliceFile: Can''t open file "',TRIM(FileName),'"'
      STOP
	  ENDIF

    WRITE(Unit,'(A,I5)') ' Number of Slices: ',NumSlice
    WRITE(Unit,*)

    NumNodeInSliceFile=0
    DO i=1,NumSlice
      WRITE(Unit,'(2(A,I5))') ' Slice number [',i,']: ',SliceIndex(i)
      WRITE(Unit,'(2(A,I5))') ' Number of nodes for the slice [',&
        & NumNodeInSlice(i),']: ',NumNodeInSlice(i)
      WRITE(Unit,'(A,200(1X,I6))') ' Slice Nodes: ',&
        & (SliceNodes(i,j),j=1,NumNodeInSlice(i))
      WRITE(Unit,*)

      NumNodeInSliceFile=NumNodeInSliceFile+NumNodeInSlice(i)
    ENDDO

    IF (verbose) THEN
      WRITE(stdout,'(A,I6,A)') 'WriteSliceFile: Slice file wrote  (',&
        & NumNodeInSliceFile,' nodes)'
      WRITE(stdout,*)
    ENDIF

    CLOSE(unit=Unit)
    RETURN
  END SUBROUTINE WriteSliceFile

  !-----------------------------------------------------------------

END MODULE slicefile
