! -*-f90-*-
MODULE thickness

	USE constants
	USE fileio
	IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------

  SUBROUTINE PromptThickness(ThickArray,NumNodes,Nodal_index,&
    &NumNodes_Slice,Nodes_Slice,NumSlices,ratio_log,Unit)

    ! This subroutine asks the user whether constant thickness is
    ! wanted for the 3D mesh. If so a thickness file is written by
    ! the routine WriteThickConstant. Otherwise the thickness for
    ! every node is prompted with some shortcuts to speed up the
    ! process.

    IMPLICIT NONE

    INTEGER(i4)::NumNodes,NumSlices,Unit
    REAL(dp)::ThickArray(:),SliceThick(NumSlices)
    INTEGER(i4)::Nodal_index(:)
    INTEGER(i4)::NumNodes_Slice(:)
    INTEGER(i4)::Nodes_Slice(:,:)
    LOGICAL::ratio_log
    ! local variables
    INTEGER(i4)::Enter_opt,Enter_opt2,limit
    INTEGER(i4)::coun,i,j,k,number,def_int,counter,ierr
    REAL(dp)::def,Real_read,Bone_thickness
    LOGICAL::error


    !WRITE(unit,'(" Thickness file",/)')
    error=.true.
    DO WHILE ( error)
      def_int=1
      WRITE(stdout,'()')
      WRITE(stdout,'(" Enter the option for specifying the Bone thickness &
        &[1]")')
      WRITE(stdout,'(" Constant Spatially  <1>")')
      WRITE(stdout,'(" By Slice            <2>")')
      WRITE(stdout,'(" By Node             <3>")')
      READ(stdin,'(A)') Buff
      READ(Buff,*,iostat=ierr) Enter_opt
      IF (Buff.eq.'' .or. ierr.ne.0) Enter_opt=def_int
      IF (Enter_opt.lt.1 .or. Enter_opt.gt.3) Enter_opt=1

      WRITE(stdout,'()')
      WRITE(stdout,'(" How is thickness calculated? [1]")')
      WRITE(stdout,'(" As ratio of distance to centre      <1>")')
      WRITE(stdout,'(" As projection by physical quantity  <2>")')
      READ(stdin,'(A)') Buff
      READ(Buff,'(I32)',iostat=ierr) Enter_opt2
      IF (Buff.eq.'' .or. ierr.ne.0) Enter_opt2=def_int
      IF (Enter_opt2.ne.1 .and. Enter_opt2.ne.2) Enter_opt2=1

      ratio_log = Enter_opt2.eq.1

      error=.False.
      def=1.0
      ! Create a header for the thickness file
      WRITE(unit,'(" Thickness file")')
      IF (ratio_log) THEN
        WRITE(unit,'(" Thickness based on projection by specified & 
          &ratio of distance to centre")')
        WRITE(unit,'()')
      ELSE
        WRITE(unit,'(" Thickness based on projection by specified &
          &quantity")')
        WRITE(unit,'()')
      END IF

      IF (Enter_opt.eq.1) THEN
        WRITE(stdout,'(" The thickness is [",E12.4,"]: ")',advance='no') Def
        READ(stdin,'(E10.4)') real_read
        Bone_thickness=Def
        IF (ABS(real_read).gt.Zero_tolerance) Bone_thickness=real_read
        ! Default cond.
        DO i=1,NumNodes
          ThickArray(i)=Bone_thickness
          WRITE(unit,'(" Thickness at node ",I5,": ",E14.5)')&
            & Nodal_index(i),ThickArray(i)
        END DO
        error=.false.
      ELSE IF (Enter_opt.eq.2) THEN
        k=1
        DO i=1,NumSlices
          WRITE(stdout,'(" Thickness at the nodes in slice",I5," is [",&
            &E9.4,"] : ")', advance='no') i,Def
          READ(stdin,'(E12.4)') real_read
          SliceThick(i)=Def
          IF (ABS(real_read).gt.Zero_tolerance) THEN
            SliceThick(i)=real_read 
            Def=real_read ! Updates the default value
          END IF
          ! cond.
          DO j=1,NumNodes_Slice(i)
            number=Nodes_Slice(i,j) ! This is the node from the slice
            coun=0
            DO WHILE (nodal_index(k) .ne. number)
              coun=coun+1
              k=k+1
              IF (k.eq.NumNodes+1) k=1
              IF (coun.eq.NumNodes+3) THEN
                WRITE(stdout,'(" ERROR: Tried to read in a node number &
                  &that does not exist")')
                STOP
              ENDIF
            END DO
            ! Loop Finds where the particular number is stored in the
            ! node array. May not be consecutively numbered whereas
            ! they are consecutively stored here. Note: if the arrays
            ! are consecutively numbered this routine works very fast.

            ! Store the node index into slice_index so that it can be
            ! easily found later without going though all the nodes
            ! again.

            ! Note: if more then one slice contains the same node and
            ! the thickness is different the thickness will be set to
            ! the value in the later slice. It is not possible to
            ! have different versions for the same node in thickness
            ! at this stage.
            thickArray(k)=SliceThick(i)
          ENDDO
        ENDDO
        DO i=1,NumNodes
          WRITE(unit,'(" Thickness at node ",i5,": ",E14.5)')& 
            & Nodal_Index(i),ThickArray(i)
        END DO
        error=.false.

      ELSE IF (Enter_opt.eq.3) THEN
        counter=0
        WRITE(stdout,'(" To accept the default for the current node: <&
          &Enter>")') 
        WRITE(stdout,'(" To accept the current default for the remaining& 
          & nodes: <-1.> at any thickness prompt.")')
        WRITE(stdout,'(" To accept the current default for the next 10 &
          &nodes: <-10.> at any thickness prompt.")')
        main:DO WHILE (i.le.Numnodes)
          WRITE(stdout,'(" Thickness at node ",i5," is [",E12.4,&
            &"] : ")',advance='no') Nodal_index(i),def
          counter=counter+1
          IF (counter.eq.10) THEN
            counter=0
            WRITE(stdout,*) ! reinitialises the record
          ENDIF
          READ(stdin,'(E12.5)') real_read
          ThickArray(i)=def
          IF (NINT(real_read) .eq. -1) THEN
            DO j=i,NumNodes
              ThickArray(j)=Def
              WRITE(unit,'(" Thickness at node ",I5,": ",E14.5)')&
                & Nodal_index(j),ThickArray(j)
            ENDDO
            EXIT main
          ELSE IF (NINT(real_read) .eq. -10) THEN
            IF (i+9 .le. Numnodes) THEN
              limit=i+9
            ELSE
              limit=Numnodes
            ENDIF
            DO j=i,limit
              ThickArray(j)=Def
              WRITE(unit,'(" Thickness at node ",I5,": ",E14.5)')&
                & Nodal_index(j),ThickArray(j)
              i=i+1
            ENDDO
          ELSE IF (real_read.gt.Zero_tolerance) THEN
            ThickArray(i)=real_read 
            Def=real_read
          ENDIF
          IF (NINT(real_read).ne.-1 .and. NINT(real_read).ne.-10) THEN
            WRITE(unit,'(" Thickness at node ",I5,": ",E14.5)')&
              & Nodal_index(i),ThickArray(i)
            i=i+1
          ENDIF
        ENDDO main
        error=.false.
      ELSE
        WRITE(stdout,'(" Incorrect option entered ")') 
        error=.true.
      END IF
    ENDDO
    WRITE(stdout,'(">> New thickness file written")')

    RETURN
  END SUBROUTINE PromptThickness

  !----------------------------------------------------------------

  SUBROUTINE ReadThickness(ThickArray,NumNodes,Nodal_index,ratio_log,unit)

    ! This subroutine reads the thickness from the .ipthic file

    IMPLICIT NONE

    INTEGER(i4)::NumNodes,Unit
    REAL(dp)::Thickarray(:)
    INTEGER(i4)::Nodal_index(:)
    LOGICAL::ratio_log
    !local variables
    INTEGER(i4)::i,j,k,nodeNumRead
    CHARACTER(LEN=5)::ratio_query

    ratio_log=.false.

    READ(unit,*)
    READ(unit,'(44X,A5)',advance='no') ratio_query
    READ(unit,*)
    IF (ratio_query.eq.'ratio') ratio_log=.true.
    READ(unit,*)
    DO i=1,NumNodes
      READ(unit,'(19X,I5)',advance='no') NodeNumRead
      IF (NodeNumRead.eq.Nodal_index(i)) THEN
        READ(unit,'(2X,E14.5)') ThickArray(i)
      ELSE
        WRITE(stdout,'(" ERROR: File out of order or incomplete")')
        WRITE(stdout,'(" At node",I5," of the input file")') NodeNumRead
        STOP
      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE ReadThickness

  !-----------------------------------------------------------------

END MODULE thickness

