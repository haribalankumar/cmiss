! -*-f90-*-
MODULE ipelemfile

	USE constants
  USE fileio
	IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------

  SUBROUTINE ReadElemDim( NumElem,Unit,FileName,verbose )

    IMPLICIT NONE
    ! Routine arguments
    INTEGER(I4)::NumElem
    ! File Arguments
    INTEGER(I4)::Unit
    CHARACTER(LEN=*)::FileName
    LOGICAL::verbose
    ! Local Variables
    INTEGER::ierr

    ! Read and Write the new element file 
    OPEN(unit=Unit,file=FileName,status='old',iostat=ierr)
    IF (ierr.NE.0) THEN
      WRITE(stderr,'(3a)') 'ReadElemDim: Can''t open file "',&
        & TRIM(FileName),'"'
      STOP
    ENDIF

    ! First read number of elements (on the third non-blank line)
    CALL GetLine( unit,buff,verbose,ierr )
    CALL GetLine( unit,buff,verbose,ierr )
    CALL GetLine( unit,buff,verbose,ierr )
    READ(buff(scan(buff,':')+1:),*) NumElem
    IF (verbose) WRITE(stdout,*) NumElem

    CLOSE(unit=Unit)
    RETURN
  END SUBROUTINE ReadElemDim

  !----------------------------------------------------------------

  SUBROUTINE ReadElemFile( NumElem,NumNode,NumCoord,NumNodePerElem,&
    & ElemIndex,ElemNode,ElemVersion,NodeIndex,NodeVersion, &
    & Unit,FileName,Verbose )

    IMPLICIT NONE
    ! Data args
    INTEGER(I4)::NumElem
    INTEGER(I4)::NumNode
    INTEGER(I4)::NumCoord
    INTEGER(I4)::NumNodePerElem
    INTEGER(I4)::ElemIndex(:)
    INTEGER(I4)::ElemNode(:,:)
    INTEGER(I4)::ElemVersion(:,:,:)
    INTEGER(I4)::NodeIndex(:)
    INTEGER(I4)::NodeVersion(:,:)
    ! File Arguments
    INTEGER(I4)::Unit
    CHARACTER(LEN=*)::FileName
    LOGICAL::verbose
    ! Local Variables
    INTEGER::i,j,k,l,ierr


    ! Initialise values to zero
    ElemNode=0
    ElemVersion(:,:,1:NumCoord)=1
    ElemVersion(:,:,NumCoord+1)=0
    k=1

    ! Read and Write the new element file 
    OPEN(unit=Unit,file=FileName,status='old',iostat=ierr)
    IF (ierr.NE.0) THEN
      WRITE(stderr,'(3a)') 'ReadElemFile: Can''t open file "',&
        & TRIM(FileName),'"'
      STOP
    ENDIF

    ! First read number of elements (on the third non-blank line)
    CALL GetLine( unit,buff,verbose,ierr )
    CALL GetLine( unit,buff,verbose,ierr )
    CALL GetLine( unit,buff,verbose,ierr )
    READ(buff(scan(buff,':')+1:),*) NumElem
    IF (verbose) WRITE(stdout,*) NumElem

    ! Read in the data for each element
    DO i=1,NumElem
      CALL GetLine( unit,buff,verbose,ierr )
      READ(buff(scan(buff,':')+1:),*) ElemIndex(i)
      IF (verbose) WRITE(stdout,*) ElemIndex(i)

      DO j=1,4 ! These are 4 lines not needed
        CALL GetLine( unit,buff,verbose,ierr )
      ENDDO

      CALL GetLine( unit,buff,verbose,ierr )
      READ(buff(scan(buff,':')+1:),*) (ElemNode(i,j),&
        & j=1,NumNodePerElem_2D)
      IF (verbose) WRITE(stdout,*) (ElemNode(i,j),&
        & j=1,NumNodePerElem_2D)

      ! Compare the elemnt with the existing node information to see
      ! if the nodes in the element have versions
      DO j=1,NumNodePerElem_2D

        ! Find the nodal_index
        DO l=1,NumNode
          IF (NodeIndex(l).eq.ElemNode(i,j)) GOTO 100
        ENDDO
        WRITE(stderr,*) 'ReadElemFile: Node ',ElemNode(i,j),&
          & 'not found'
        STOP
100     CONTINUE

        ! Read the values in
        DO k=1,NumCoord
          IF (NodeVersion(l,k).gt.1) THEN
            ElemVersion(i,j,NumCoord+1)=1
            CALL GetLine( unit,buff,verbose,ierr )
            READ(buff(scan(buff,':')+1:),*,iostat=ierr) ElemVersion(i,j,k)
            IF (ierr.ne.0) THEN
              WRITE(stderr,'(2A,/,3A)') 'ReadElemFile: error reading ',&
                & 'element versions from line','''',TRIM(buff),''''
              STOP
            ENDIF
            IF (verbose) WRITE(stdout,*) ElemVersion(i,j,k)
          ELSE
            ElemVersion(i,j,NumCoord+1)=0
            ElemVersion(i,j,k)=1
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    CLOSE(unit=Unit)
    RETURN
  END SUBROUTINE ReadElemFile

  !----------------------------------------------------------------

  SUBROUTINE WriteElemFile(NumElem,NumCoord,NumNodePerElem,&
    & ElemIndex,ElemNodes,ElemVersion,New_Basis,Unit,FileName,Heading,Out)

    IMPLICIT NONE
    ! Data values
    INTEGER(I4)::NumElem
    INTEGER(I4)::NumCoord
    INTEGER(I4)::NumNodePerElem
    INTEGER(I4)::ElemIndex(:)
    INTEGER(I4)::ElemNodes(:,:)
    INTEGER(I4)::ElemVersion(:,:,:)
    INTEGER(I4)::New_Basis
    ! File Arguments
    INTEGER(I4)::Unit
    CHARACTER(LEN=*)::FileName
    CHARACTER(LEN=*)::Heading
    LOGICAL::Out
    ! Local Variables
    INTEGER(I4)::i,j,k,ierr
    INTEGER(I4)::occurance(NumNodePerElem)

    ! Write the new element file 
    OPEN(unit=Unit,file=FileName,status='unknown',iostat=ierr)
    IF (ierr.NE.0) THEN
      WRITE(stderr,'(3a)') 'WriteElemFile: Can''t open file "',&
        & TRIM(FileName),'"'
      STOP
    ENDIF

    ! Write the header information
    WRITE(Unit,'(A)') ' CMISS Version 1.21 ipelem File Version 2'
    WRITE(Unit,'(2A)') ' Heading: ',TRIM(ADJUSTL(Heading))
    WRITE(Unit,*)

    WRITE(Unit,'(A,I6)') ' The number of elements is [1]: ',NumElem
    WRITE(Unit,*)

    ! Write out the elements
    DO i=1,NumElem

      ! Element info
      WRITE(Unit,'(2(A,I6))') ' Element number [',i,']: ',ElemIndex(i)
      WRITE(Unit,'(2(A,I1))') ' The number of geometric Xj-coordinates&
        & is [',NumCoord,']: ',NumCoord
      DO j=1,NumCoord
        WRITE(Unit,'(2(A,I2))') ' The basis function type for&
          & geometric variable',j,' is [1]: ',1
      ENDDO

      ! Node numbers for element
      WRITE(Unit,'(A,I2,A,I5,10(X,I5)))') ' Enter the',NumNodePerElem,' global&
        & numbers for basis 1: ',(ElemNodes(i,j),j=1,NumNodePerElem)
      ! Version information for nodes; first count the number
      ! of times any node occurs in the element.
      occurance(:)=1
      DO j=2,NumNodePerElem
        DO k=j-1,1,-1
          IF (ElemNodes(i,k).eq.ElemNodes(i,j)) THEN
            occurance(j)=occurance(k)+1
            GOTO 100
          ENDIF
        ENDDO
100     CONTINUE
      ENDDO

      ! Write out the version numbers. Do these numbers have to be
      ! grouped Do these values have to be grouped so that all
      ! versions of one node are in the one clump?
      DO j=1,NumNodePerElem
        IF (ElemVersion(i,j,NumCoord+1).eq.1) THEN
          DO k=1,NumCoord
            WRITE(Unit,'(A,I3,A,I5,A,I1,A,I2)' ) ' The version number&
              & for occurrence',occurance(j),' of node ',ElemNodes(i&
              & ,j),', njj=',k,' is [ 1]: ',ElemVersion(i,j,k)
          ENDDO
        ENDIF
      ENDDO

      ! Write out addition basis ordering
      IF(New_Basis.gt.0) THEN
        WRITE(Unit,'(2(A,I2),A,I5,10(X,I5)))') ' Enter the',NumNodePerElem,' numbers&
          & for basis',New_Basis,' [prev]: ',(ElemNodes(i,j),j=1,NumNodePerElem)
      ENDIF


      WRITE(Unit,*)
    ENDDO

    RETURN
  END SUBROUTINE WriteElemFile

  !----------------------------------------------------------------
  
END MODULE ipelemfile

