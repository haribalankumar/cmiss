! -*-f90-*-
MODULE ipnodefile

	USE constants
  USE fileio
	IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------

  SUBROUTINE ReadIPNodeDim( NumNode,NumCoord,MaxNumDeriv, &
    & MaxNumVersions,Unit,FileName,Out )

    ! Reads dimensions of the node data contained in an ipnode file

    IMPLICIT NONE
    ! Data arguments
    INTEGER(I4)::NumNode
    INTEGER(I4)::NumCoord
    INTEGER(I4)::MaxNumDeriv
    INTEGER(I4)::MaxNumVersions
    ! File arguments
    INTEGER(I4)::Unit
    CHARACTER(LEN=*)::FileName
    LOGICAL::Out
    ! Internal variables
    INTEGER(I4)::i,j,k,v,ierr,idum,node_version
    INTEGER(I4), ALLOCATABLE::num_deriv(:)
    LOGICAL, ALLOCATABLE::versions(:)


    MaxNumDeriv=0
    MaxNumVersions=0

    OPEN(unit=Unit,file=FileName,status='old',iostat=ierr)
    IF (ierr.ne.0) THEN
      WRITE(stderr,*) 'ReadIPNodeDim: Cannot find ipnode file  "'&
        & ,TRIM(FileName),'"'
      STOP
    ENDIF

    CALL GetLine( Unit,buff,out,ierr )
    CALL GetLine( Unit,buff,out,ierr )

    CALL GetLine( Unit,buff,out,ierr )
    READ(buff(scan(buff,':')+1:),*) NumNode

    CALL GetLine( Unit,buff,out,ierr )
    READ(buff(scan(buff,':')+1:),*) NumCoord

    ! Alloc and initialise
    ALLOCATE(versions(NumCoord), &
           & num_deriv(NumCoord),&
           & stat=ierr)
    IF (ierr.ne.0) THEN
      WRITE(stderr,*) 'ReadIPNodeDim: alloc failure'
      STOP
    ENDIF
    num_deriv=0
    versions=.FALSE.

    ! Read in version flags
    DO i=1,NumCoord
      CALL GetLine( Unit,buff,out,ierr )
      buff=ADJUSTL(buff(scan(buff,'?')+1:))

      IF (TRIM(buff).eq.'y' .or. TRIM(buff).eq.'Y') THEN
        versions(i)=.TRUE.
      ENDIF
    ENDDO

    ! Read in derivative info
    DO i=1,NumCoord
      CALL GetLine( Unit,buff,out,ierr )
      READ(buff(scan(buff,':')+1:),*) num_deriv(i)
    ENDDO
    MaxNumDeriv=MAXVAL(num_deriv)

    ! Read the rest of the file
    DO i=1,NumNode

      CALL GetLine( unit,buff,out,ierr )

      ! Loop over the coordinates
      DO j=1,NumCoord

        ! If we have versions
        IF (versions(j)) THEN
          CALL GetLine( unit,buff,out,ierr )
          READ(buff(scan(buff,':')+1:),*) node_version
        ELSE
          node_version=1
        ENDIF
        MaxNumVersions=MAX(MaxNumVersions,node_version)

        DO v=1,node_version
          IF (node_version.gt.1) THEN
            CALL GetLine( unit,buff,out,ierr )
          ENDIF

          CALL GetLine( unit,buff,out,ierr )

          IF (num_deriv(j).eq.0) THEN
          ELSE IF (num_deriv(j).eq.1) THEN
            CALL GetLine( unit,buff,out,ierr )

          ELSE IF (num_deriv(j).eq.3) THEN
            DO k=2,4
              CALL GetLine( unit,buff,out,ierr )
            ENDDO

          ELSE
            WRITE(stderr,*) 'ReadIPNodeDim: ',num_deriv(j),&
              & ' is a strange number of derivatives'
            STOP
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    CLOSE(unit=Unit)

    DEALLOCATE(versions,num_deriv,stat=ierr)
    IF (ierr.ne.0) THEN
      WRITE(stderr,*) 'ReadIPNodeDim: dealloc failure'
      STOP
    ENDIF

    RETURN
  END SUBROUTINE ReadIPNodeDim

  !----------------------------------------------------------------

  SUBROUTINE ReadIPNodeFile( NumNode,NumCoord,MaxNumDeriv, &
    & MaxNumVersions,NumDeriv,NodeIndex,NodeVersion,NodeValue, &
    & Unit,FileName,Out )

    ! Reads nodal data from an ipnode file. Note that ReadIPNodeDim
    ! should have been called first to find the dimensions of the data.

    IMPLICIT NONE
    ! Data arguments
    INTEGER(I4)::NumNode
    INTEGER(I4)::NumCoord
    INTEGER(I4)::MaxNumDeriv
    INTEGER(I4)::MaxNumVersions
    INTEGER(I4)::NumDeriv(:)
    INTEGER(I4)::NodeIndex(:)
    INTEGER(I4)::NodeVersion(:,:)
    REAL(dp)::NodeValue(:,:,:,:)
    ! File arguments
    INTEGER(I4)::Unit
    CHARACTER(LEN=*)::FileName
    LOGICAL::Out
    ! Internal variables
    INTEGER(I4)::i,j,k,v,ierr,idum
    LOGICAL, ALLOCATABLE::versions(:)


    OPEN(unit=Unit,file=FileName,status='old',iostat=ierr)
    IF (ierr.ne.0) THEN
      WRITE(stderr,*) 'ReadIPNodeFile: Cannot find ipnode file  "'&
        & ,TRIM(FileName),'"'
      STOP
    ENDIF

    CALL GetLine( Unit,buff,out,ierr )
    CALL GetLine( Unit,buff,out,ierr )

    CALL GetLine( Unit,buff,out,ierr )
    READ(buff(scan(buff,':')+1:),*) NumNode

    CALL GetLine( Unit,buff,out,ierr )
    READ(buff(scan(buff,':')+1:),*) NumCoord

    ! Alloc and initialise
    ALLOCATE(versions(NumCoord),&
           & stat=ierr)
    IF (ierr.ne.0) THEN
      WRITE(stderr,*) 'ReadIPNodeFile: alloc failure'
      STOP
    ENDIF
    versions=.FALSE.
    NumDeriv=0

    ! Read in version flags
    DO i=1,NumCoord
      CALL GetLine( Unit,buff,out,ierr )
      buff=ADJUSTL(buff(scan(buff,'?')+1:))

      IF (TRIM(buff).eq.'y' .or. TRIM(buff).eq.'Y') THEN
        Versions(i)=.TRUE.
      ENDIF
    ENDDO

    ! Read in derivative info
    DO i=1,NumCoord
      CALL GetLine( Unit,buff,out,ierr )
      READ(buff(scan(buff,':')+1:),*) NumDeriv(i)
    ENDDO
    IF (MAXVAL(NumDeriv).gt.MaxNumDeriv) THEN
      WRITE(stderr,*) 'ReadIPNodeFile: Max num derivative too large;',&
        & MAXVAL(NumDeriv),'>',MaxNumDeriv
      STOP
    ENDIF

    ! Read the rest of the file
    DO i=1,NumNode

      CALL GetLine( unit,buff,out,ierr )
      READ(buff(scan(buff,':')+1:),*) NodeIndex(i)

      ! Loop over the coordinates
      DO j=1,NumCoord

        ! If we have versions
        IF (Versions(j)) THEN
          CALL GetLine( unit,buff,out,ierr )
          READ(buff(scan(buff,':')+1:),*) NodeVersion(i,j)
        ELSE
          NodeVersion(i,j)=1
        ENDIF
        IF (NodeVersion(i,j).gt.MaxNumVersions) THEN
          WRITE(stderr,*) 'ReadIPNodeFile: exceeded maximum'// &
            & ' number of versions. ',NodeVersion(i,j),'>',&
            & MaxNumVersions
          STOP
        ENDIF

        DO v=1,NodeVersion(i,j)
          IF (NodeVersion(i,j).gt.1) THEN
            CALL GetLine( unit,buff,out,ierr )
          ENDIF

          CALL GetLine( unit,buff,out,ierr )
          READ(buff(scan(buff,':')+1:),*) NodeValue(i,j,v,1)
          IF (out) WRITE(stdout,*) NodeValue(i,j,v,1)

          IF (NumDeriv(j).eq.0) THEN
          ELSE IF (NumDeriv(j).eq.1) THEN
            CALL GetLine( unit,buff,out,ierr )
            READ(buff(scan(buff,':')+1:),*) NodeValue(i,j,v,2)
            IF (out) WRITE(stdout,*) NodeValue(i,j,v,2)

          ELSE IF (NumDeriv(j).eq.3) THEN
            DO k=2,4
              CALL GetLine( unit,buff,out,ierr )
              READ(buff(scan(buff,':')+1:),*) NodeValue(i,j,v,k)
              IF (out) WRITE(stdout,*) NodeValue(i,j,v,k)
            ENDDO

          ELSE
            WRITE(stderr,*) 'ReadIPNodeFile: ',NumDeriv(j),&
              & ' is a strange number of derivatives'
            STOP
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    CLOSE(unit=Unit)

    DEALLOCATE(versions,stat=ierr)
    IF (ierr.ne.0) THEN
      WRITE(stderr,*) 'ReadIPNodeFile: dealloc failure'
      STOP
    ENDIF

    RETURN
  END SUBROUTINE ReadIPNodeFile

  !----------------------------------------------------------------

  SUBROUTINE WriteIPNodeFile(NumNode,NumCoord,NumDeriv,&
    & NodeIndex,NodeVersion,NodeValue,Unit,FileName,Heading,Out)

    ! Writes nodal values to an ipnode file.

    IMPLICIT NONE
    ! Data values
    INTEGER(I4)::NumNode
    INTEGER(I4)::NumCoord
    INTEGER(I4)::NumDeriv(:)
    INTEGER(I4)::NodeIndex(:)
    INTEGER(I4)::NodeVersion(:,:)
    REAL(dp)::NodeValue(:,:,:,:)
    ! File arguments
    INTEGER(I4)::Unit
    CHARACTER(LEN=*)::FileName
    CHARACTER(LEN=*)::Heading
    LOGICAL::Out

    ! Internal Parameters
    LOGICAL,PARAMETER::debugnode=.false.
    ! Internal variables
    INTEGER(I4)::i,j,k,v,ierr,idum


    OPEN(unit=Unit,file=FileName,status='unknown',iostat=ierr)
    IF (ierr.ne.0) THEN
      WRITE(stderr,*) 'WriteIPNodeFile: Can''t open ipnode file  "'&
        & ,TRIM(FileName),'"'
      STOP
    ENDIF

    WRITE(Unit,'(A)')  ' CMISS Version 1.21 ipnode File Version 2'
    WRITE(Unit,'(2A)') ' Heading: ',TRIM(ADJUSTL(Heading))
    WRITE(Unit,*)

    WRITE(Unit,'(2(A,I5))')' The number of nodes is [',NumNode,']: ',NumNode
    WRITE(Unit,'(2(A,I1))')' Number of coordinates [',NumCoord,']: ',NumCoord
    DO i=1,NumCoord
      WRITE(Unit,'(A,I1,A))')' Do you want prompting for different&
        & versions of nj=',i,' [N]? Y'
    ENDDO
    DO i=1,NumCoord
      WRITE(Unit,'(A,I1,A,I1)')' The number of derivatives for&
        & coordinate ',i,' is [0]: ',NumDeriv(i)
    ENDDO
    WRITE(Unit,*)

    DO i=1,NumNode
      WRITE(Unit,'(2(A,I5))')' Node number [',NodeIndex(i),']: '&
        & ,NodeIndex(i)

      DO j=1,NumCoord
        WRITE(Unit,'(A,I1,A,I2)')' The number of versions for nj=',j&
          & ,' is [1]: ',NodeVersion(i,j)

        DO k=1,NodeVersion(i,j)
          IF (NodeVersion(i,j).gt.1) THEN
            WRITE(Unit,'(A,I2,A)')' For version number',k,':'
          ENDIF

          WRITE(Unit,'(A,I1,A,E12.5,A,E25.17)')' The Xj(',j,') coordinate&
            & is [',0.0D0,']: ',NodeValue(i,j,k,1)
          IF (NumDeriv(j).eq.3) THEN
            WRITE(Unit,'(A,E12.5,A,E25.17))')' The derivative wrt&
              & direction 1 is [',0.0D0,']: '&
              & ,NodeValue(i,j,k,2)
            WRITE(Unit,'(A,E12.5,A,E25.17))')' The derivative wrt&
              & direction 2 is [',0.0D0,']: '&
              & ,NodeValue(i,j,k,3)
            WRITE(Unit,'(A,E12.5,A,E25.17))')' The derivative wrt&
              & directions 1 & 2 is [',0.0D0,']: '&
              & ,NodeValue(i,j,k,4)
          ENDIF
        ENDDO
      ENDDO

      WRITE(Unit,*)
    ENDDO
    CLOSE(unit=Unit)

    RETURN
  END SUBROUTINE WriteIPNodeFile

  !----------------------------------------------------------------

END MODULE ipnodefile
