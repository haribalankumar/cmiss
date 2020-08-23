! -*-f90-*-
MODULE createnodes

  USE constants
  USE spline
  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------

  SUBROUTINE NewElems( NumElem, NumNode, NumCoord, NumNodeVers,  &
    & NumNodePerElem_2D, NumNodePerElem_3D,                      &
    & ElemNode_old,ElemVersion_old,ElemNode_new,ElemVersion_new, &
    & NumSlice,NodeSliceCount,NumNodeInSlice,SliceNodes,         &
    & NodeProjectIndex,NodeIndex,NodeVersion,verbose  )

    IMPLICIT none
    ! Arguments
    INTEGER(I4)::NumElem
    INTEGER(I4)::NumNode
    INTEGER(I4)::NumCoord
    INTEGER(I4)::NumSlice
    INTEGER(I4)::NumNodeVers
    INTEGER(I4)::NumNodePerElem_2D
    INTEGER(I4)::NumNodePerElem_3D
    ! Data arrays
    INTEGER(I4)::Elemnode_old(:,:)
    INTEGER(I4)::ElemVersion_old(:,:,:)
    INTEGER(I4)::Elemnode_new(:,:)
    INTEGER(I4)::ElemVersion_new(:,:,:)
    ! Node information
    INTEGER(I4)::NodeProjectIndex(:,:)
    INTEGER(I4)::NodeIndex(:)
    INTEGER(I4)::NodeVersion(:,:)
    ! Slice Info
    INTEGER(I4)::NodeSliceCount(:)
    INTEGER(I4)::NumNodeInSlice(:)
    INTEGER(I4)::SliceNodes(:,:)
    ! Other
    LOGICAL::Verbose
    ! Local variables
    INTEGER(I4)::i,j,k,l,num_vers,tmp_node
    INTEGER(I4)::old_node,old_vers
    INTEGER(I4)::new_node,new_vers


    ! Copy across old values
    DO i=1,NumElem
      DO j=1,NumNodePerElem_2D
        ! Copy old nodes
        Elemnode_new(i,j)=Elemnode_old(i,j)

        ! Copy old versions
        DO k=1,NumCoord+1
          ElemVersion_new(i,j,k)=ElemVersion_old(i,j,k)
        ENDDO
      ENDDO
    ENDDO

    ! Set the new values
    DO i=1,NumElem
      ! Set new elements
      DO j=1,NumNodePerElem_2D
        ! Find the node/version
        old_node=Elemnode_old(i,j)
        old_vers=ElemVersion_old(i,j,1)

        ! Cope with cases where more than one slice per node
        CALL FindNewProjectionPoint( tmp_node,old_node,NumNode,NumSlice,    &
    &     NodeIndex,NodeSliceCount,NumNodeInSlice,SliceNodes,Elemnode_old(i,:), &
    &     verbose  )
        IF (tmp_node.ne.old_node) THEN
          old_node=tmp_node
          old_vers=1
        ENDIF

        DO k=1,NumNodeVers
          IF (NodeProjectIndex(k,1).eq.old_node .and. &
            & NodeProjectIndex(k,2).eq.old_vers) GOTO 100
        ENDDO
        WRITE(stderr,*) 'NewElems: Can''t find node/version ',old_node,old_vers
        STOP
100     CONTINUE
        new_node=NodeProjectIndex(k,3)
        new_vers=NodeProjectIndex(k,4)

        DO k=1,NumNode
          IF (NodeIndex(k).eq.new_node) GOTO 200
        ENDDO
        WRITE(stderr,*) 'NewElems: Can''t find node ',new_node
        STOP
200     CONTINUE
        num_vers=NodeVersion(k,1)

        ! Set the new nodes
        Elemnode_new(i,j+NumNodePerElem_2D)=new_node
        ! Set the new versions
        DO l=1,NumCoord
          ElemVersion_new(i,j+NumNodePerElem_2D,l)=new_vers
        ENDDO

        IF (num_vers.gt.1) THEN
          ElemVersion_new(i,j+NumNodePerElem_2D,NumCoord+1)=1
        ELSE
          ElemVersion_new(i,j+NumNodePerElem_2D,NumCoord+1)=0
        ENDIF

      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE NewElems

  !----------------------------------------------------------------

  SUBROUTINE NewNodes( NumNode,NumCoord,MaxVersion,MaxSlice, &
    & NumDeriv,node_offset,NumNode_new,MaxNewNode,           &
    & NodeValue_old,NodeIndex_old,NodeVersion_old,           &
    & NodeValue_new,NodeIndex_new,NodeVersion_new,           &
    & NumSlice,NumNodeInSlice,SliceNodes,node_thick,         &
    & NodeProjectIndex,NumNodeVers,MaxNumNodeVers,           &
    & fudge_alpha,ratio_thick,one_node,new_versions,verbose )

    IMPLICIT NONE
    INTEGER(I4)::NumNode,NumCoord,NumDeriv(3)
    INTEGER(I4)::MaxVersion,MaxSlice,MaxNewNode
    ! Old data
    INTEGER(I4)::NodeIndex_old(:)
    INTEGER(I4)::NodeVersion_old(:,:)
    REAL(dp)::NodeValue_old(:,:,:,:)
    ! New data
    INTEGER(I4)::NodeIndex_new(:)
    INTEGER(I4)::NodeVersion_new(:,:)
    REAL(dp)::NodeValue_new(:,:,:,:)
    INTEGER(I4)::NodeProjectIndex(:,:)
    INTEGER(I4)::NumNodeVers,MaxNumNodeVers
    ! Slice info
    INTEGER(I4)::NumSlice
    INTEGER(I4)::NumNodeInSlice(:)
    INTEGER(I4)::SliceNodes(:,:)
    REAL(dp)::node_thick(:)
    ! Other crap
    INTEGER(I4)::node_offset
    INTEGER(I4)::NumNode_new
    REAL(DP)::fudge_alpha
    LOGICAL::ratio_thick
    LOGICAL::one_node
    LOGICAL::new_versions
    LOGICAL::verbose
    ! Internal variables
    INTEGER(I4)::i,j,k,l,m,n,o,index(MaxSlice)
    INTEGER(I4)::maxnodenum,offset,num_ver_tmp
    INTEGER(I4)::nminus,nplus,npoint,num_point
    REAL(dp)::centre(NumCoord),bone_centre(NumCoord),ALPHA
    LOGICAL::pre,post


    IF (one_node) THEN
      WRITE(stderr,*) 'Generating One node per slice'
    ELSE
      WRITE(stderr,*) 'Generating N nodes per slice'
    ENDIF

    maxnodenum=maxval(NodeIndex_old(1:NumNode))
    offset=maxnodenum+node_offset

    ! Set to zero
    NumNodeVers=0
    NodeProjectIndex=0

    NodeIndex_new=0
    NodeVersion_new=0
    NodeValue_new=0

    ! Copy across old node values
    DO i=1,NumNode
      IF (verbose) WRITE(stderr,*) 'Have node ',NodeIndex_old(i),' in set,'
      NodeIndex_new(i)=NodeIndex_old(i)
      DO j=1,NumCoord
        NodeVersion_new(i,j)=NodeVersion_old(i,j)
        DO k=1,NodeVersion_new(i,j)
          NodeValue_new(i,j,k,:)=NodeValue_old(i,j,k,:)
        ENDDO
      ENDDO
    ENDDO
    IF (verbose) WRITE(stderr,*) 'Have ',NumNode,' nodes in set'

    n=NumNode
    o=0

    !---------------------------------------------
    ! Calculate the centre of the bone
    !---------------------------------------------
    num_point=0
    bone_centre=0.0
    DO i=1,NumSlice
      DO j=1,NumNodeInSlice(i)
        DO k=1,NumCoord
          bone_centre(k)=bone_centre(k)+NodeValue_old(m,k,1,1)  
          num_point = num_point+1
        ENDDO
      ENDDO
    ENDDO
    DO k=1,NumCoord
      bone_centre(k)=bone_centre(k)/float(num_point)
    ENDDO

    !---------------------------------------------
    ! For each slice, creat the new internal nodes
    !---------------------------------------------
    DO i=1,NumSlice
      pre  = .false.
      post = .false.

      ! Calculate the point which we will project to.

      ! If we have more than 2 nodes in a slice, then use the
      ! geometric centre of the slice
      IF (NumNodeInSlice(i).gt.2) THEN
        centre=0.0
        DO j=1,NumNodeInSlice(i)
          ! Find the node
          DO m = 1,NumNode
            IF (SliceNodes(i,j).eq.NodeIndex_old(m)) GOTO 100
          ENDDO
          WRITE(stderr,*) 'NewNodes: Can''t find node ',SliceNodes(i,j)
          STOP

100       CONTINUE
          index(j)=m

          DO k=1,NumCoord
            centre(k)=centre(k)+NodeValue_old(m,k,1,1)
          ENDDO
        ENDDO

        DO k=1,NumCoord
          centre(k) = centre(k)/float(NumNodeInSlice(i))
        ENDDO

      ! If there are only one or two points in slice we project to
      ! the centre of a neighbouring slice.
      ELSE

        ! Find neighbouring slice we can project to: should project to
        ! closest slice rather than next in list.
        DO j=1,NumSlice-1
          nplus  = i+j
          nminus = i-j
          IF (nplus.le.NumSlice) THEN
            IF (NumNodeInSlice(nplus).gt.2) THEN
              npoint=nplus
              post=.true.
              GOTO 200
            ENDIF
          ENDIF
          IF (nminus.ge.1) THEN
            IF (NumNodeInSlice(nminus).gt.2) THEN
              npoint=nminus
              pre=.true.
              GOTO 200
            ENDIF
          ENDIF
        ENDDO
        WRITE(stderr,*) 'NewNodes: Can''t find slice with more than two nodes'
        STOP

200     CONTINUE
        ! Calculate centre
        centre=0.0
        DO j=1,NumNodeInSlice(npoint)
          ! Find the node
          DO m = 1,NumNode
            IF (SliceNodes(npoint,j).eq.NodeIndex_old(m)) GOTO 300
          ENDDO
          WRITE(stderr,*) 'NewNodes: Can''t find node ',SliceNodes(npoint,j)
          STOP

300       CONTINUE
          DO k=1,NumCoord
            centre(k)=centre(k)+NodeValue_old(m,k,1,1)
          ENDDO
        ENDDO
        DO k=1,NumCoord
          centre(k) = centre(k)/float(NumNodeInSlice(npoint))
        ENDDO

        ! Load index
        DO j=1,NumNodeInSlice(i)
          ! Find the node
          DO m = 1,NumNode
            IF (SliceNodes(i,j).eq.NodeIndex_old(m)) GOTO 400
          ENDDO
          WRITE(stderr,*) 'NewNodes: Can''t find node ',SliceNodes(i,j)
          STOP
400       CONTINUE
          index(j)=m
        ENDDO

      ENDIF

      ! Scale the centre position -- value is now supplied by
      ! a command line arg.
      ! This is not recomended -- instead use the redistribution
      ! using a spline.
      ! ALPHA=0.95
      ALPHA=fudge_alpha
      DO k=1,NumCoord
        centre(k)=(1.0d0-ALPHA)*bone_centre(k) + (ALPHA)*centre(k)
      ENDDO

      !------------------------------
      ! One central point
      !------------------------------
      IF (one_node) THEN

        IF (.not.pre) THEN
          n=n+1
          o=o+1
        ENDIF
        IF (verbose) THEN
          IF (post) THEN
            WRITE(stderr,*) 'Adding node ',o+offset,' to set (will be appended to)'
          ELSE IF (PRE) THEN
            WRITE(stderr,*) 'Appending to node ',o+offset
          ELSE
            WRITE(stderr,*) 'Adding node ',o+offset,' to set'
          ENDIF
        ENDIF

        IF (n.gt.MaxNewNode) THEN
          WRITE(stderr,*) 'Exceeded MaxNewNode ; ',n,MaxNewNode
          STOP
        ENDIF

        ! Create the new point
        NodeIndex_new(n)=o+offset

        ! Add versions to the new points, if flagged to do so
        IF (new_versions) THEN
          DO m=1,NumNodeInSlice(i)

            DO j=1,NumCoord
              num_ver_tmp=NodeVersion_new(n,j)+NodeVersion_old(index(m),j)
              IF (num_ver_tmp.gt.MaxVersion) THEN
                WRITE(stderr,*) 'Exceeded MaxVersion ; ',n,num_ver_tmp,MaxVersion
                STOP
              ENDIF

              DO k=1,NodeVersion_old(index(m),j)
                l=NodeVersion_new(n,j)+k
                NodeValue_new(n,j,l,1)=centre(j)
                !NodeValue_new(n,j,l,2:)=NodeValue_old(index(m),j,k,2:)
                NodeValue_new(n,j,l,2:)=0.0

                IF (j.eq.1) THEN
                  ! Set the correspondence of the index of the surface
                  ! and projected points.
                  NumNodeVers=NumNodeVers+1
                  IF (NumNodeVers.gt.MaxNumNodeVers) THEN
                    WRITE(stderr,*) 'NewNodes: Not enough space in NodeProjectIndex ;',NumNodeVers
                    STOP
                  ENDIF
                  NodeProjectIndex(NumNodeVers,1)=NodeIndex_old(index(m))
                  NodeProjectIndex(NumNodeVers,2)=k
                  NodeProjectIndex(NumNodeVers,3)=NodeIndex_new(n)
                  NodeProjectIndex(NumNodeVers,4)=NodeVersion_new(n,j)+k
                  !print '(10(1x,i5))', NumNodeVers,NodeIndex_old(index(j))&
                  !  & ,NodeIndex_new(n),k,NodeVersion_new(n,j)+k
                ENDIF
              ENDDO
              NodeVersion_new(n,j)=num_ver_tmp
            ENDDO
          ENDDO

        ! No versions to be added
        ELSE
          DO j=1,NumCoord
            NodeVersion_new(n,j)=1
            NodeValue_new(n,j,1,1)=centre(j)
            NodeValue_new(n,j,l,2:)=0.0
          ENDDO

          ! Set the correspondence of the index of the surface
          ! and projected points.
          DO m=1,NumNodeInSlice(i)
            DO k=1,NodeVersion_old(index(m),1)
              NumNodeVers=NumNodeVers+1
              NodeProjectIndex(NumNodeVers,1)=NodeIndex_old(index(m))
              NodeProjectIndex(NumNodeVers,2)=k
              NodeProjectIndex(NumNodeVers,3)=NodeIndex_new(n)
              NodeProjectIndex(NumNodeVers,4)=1
              !print '(10(1x,i5))', NumNodeVers,NodeIndex_old(index(j))&
              !  & ,NodeIndex_new(n),k,NodeVersion_new(n,j)+k
            ENDDO
          ENDDO
        ENDIF

        IF (post.and.(i.ne.NumSlice)) THEN
          n=n-1
          o=o-1
        ENDIF

      !------------------------------
      ! Project towards the centre
      !------------------------------
      ELSE
        ! Project the points toward the centre
        DO j=1,NumNodeInSlice(i)
          n=n+1
          o=o+1
          IF (verbose) WRITE(stderr,*) 'Adding node ',o,' to set'
          IF (n.gt.MaxNewNode) THEN
            WRITE(stderr,*) 'NewNodes: Exceeded MaxNewNode ; ',n,MaxNewNode
            STOP
          ENDIF

          IF (ratio_thick) THEN
            ALPHA=MIN(1.0d0,node_thick(i))
          ELSE
            ALPHA=0.0d0
            DO k=1,NumCoord
              ALPHA=ALPHA+(centre(k)-NodeValue_old(index(j),k,1,1))**2
            ENDDO
            IF (SQRT(ALPHA).ne.0.0d0) THEN
              ALPHA=MIN(1.0d0,node_thick(i)/SQRT(ALPHA))
            ENDIF
          ENDIF

          ! Create the new point
          NodeIndex_new(n)=o+offset
          NodeVersion_new(n,:)=NodeVersion_old(index(j),:)
          NodeValue_new(n,:,:,:)=0
          NodeValue_new(n,:,:,2:)=NodeValue_old(index(j),:,:,2:)

          DO k=1,NumCoord
            DO l=1,NodeVersion_new(n,k)
              NodeValue_new(n,k,l,1)&
                & = (1.0d0-ALPHA)*NodeValue_old(index(j),k,l,1) &
                & + ALPHA*centre(k)

              IF (k.eq.1) THEN
                ! Set the correspondence of the index of the surface
                ! and projected points.
                NumNodeVers=NumNodeVers+1
                IF (NumNodeVers.gt.MaxNumNodeVers) THEN
                  WRITE(stderr,*) 'NewNodes: Not enough space in NodeProjectIndex ;',NumNodeVers
                  STOP
                ENDIF
                NodeProjectIndex(NumNodeVers,1)=NodeIndex_old(index(j))
                NodeProjectIndex(NumNodeVers,2)=l
                NodeProjectIndex(NumNodeVers,3)=NodeIndex_new(n)
                NodeProjectIndex(NumNodeVers,4)=l
                !print '(10(1x,i5))', NumNodeVers,NodeIndex_old(index(j))&
                !  & ,NodeIndex_new(n),l,l
              ENDIF
            ENDDO
          ENDDO

        ENDDO


      ENDIF
    ENDDO

    NumNode_new=n

    !DO i = 1,NumNode_new
    !  print *, i,NodeIndex_new(i),NodeVersion_new(i,:)      
    !ENDDO

    RETURN
  END SUBROUTINE NewNodes

  !----------------------------------------------------------------

  SUBROUTINE RedistPoints( NumRedist,Redist, &
    & NumNode,NumCoord,NodeIndex,NodeVersion,NodeValue )

    ! Redistribute points along a specified spline

    IMPLICIT NONE
    ! Arguments
    INTEGER(I4)::NumRedist
    CHARACTER(LEN=*)::Redist(:)
    INTEGER(I4)::NumNode
    INTEGER(I4)::NumCoord
    INTEGER(I4)::NodeIndex(:)
    INTEGER(I4)::NodeVersion(:,:)
    REAL(DP)::NodeValue(:,:,:,:)
    ! Local veriables
    INTEGER(I4)::i,j,k,l,idx,npoints,nspline,ierr,jerr
    INTEGER(I4),ALLOCATABLE::node_points(:),node_spline(:)
    REAL(DP),ALLOCATABLE::x_points(:,:),x_spline(:,:)
    REAL(DP)::ALPHA

    DO i=1,NumRedist
      ! Count the number of variables
      k=scan(Redist(i),':')
      npoints=1
      DO j=1,k-1
        IF (Redist(i)(j:j).eq.',') npoints=npoints+1
      ENDDO
      nspline=1
      DO j=k+1,LEN(Redist(i))
        IF (Redist(i)(j:j).eq.',') nspline=nspline+1
      ENDDO
      ALLOCATE(node_points(npoints), &
        &      node_spline(nspline), &
        &      x_points(npoints,NumCoord), &
        &      x_spline(nspline,NumCoord), &
        &      stat=ierr)
      IF (ierr.ne.0) THEN
        WRITE(stderr,*) 'RedistPoints: alloc failure'
        STOP
      ENDIF
      ! Read them in
      READ(Redist(i)(:k-1),*,iostat=ierr) (node_points(j),j=1,npoints)
      READ(Redist(i)(k+1:),*,iostat=jerr) (node_spline(j),j=1,nspline)
      IF (ierr.ne.0 .or. jerr.ne.0) THEN
        WRITE(stderr,*) 'RedistPoints: error reading points from ''',&
          & TRIM(Redist(i)),''''
        STOP
      ENDIF
      WRITE(*,*) 'Redistributing: ',i,' of ',NumRedist
      WRITE(*,*) '   Points to move:         ',npoints,':',node_points
      WRITE(*,*) '   Points defining spline: ',nspline,':',node_spline
      WRITE(*,*)

      ! Find the node coordinates
      DO j=1,nspline
        DO k=1,NumNode
          IF (NodeIndex(k).eq.node_spline(j)) GOTO 100
        ENDDO
        WRITE(stderr,*) 'RedistPoints[1]: unable to find node ',node_spline(j)
        STOP
100     CONTINUE
        idx=k
        DO k=1,NumCoord
          x_spline(j,k)=NodeValue(idx,k,1,1)
        ENDDO

        IF(node_points(1).ne.node_spline(1) .or. &
          &node_points(npoints).ne.node_spline(nspline)) THEN
          WRITE(stderr,*) 'Error!'
          WRITE(stderr,*) 'RedistPoints: the strange spline'//&
            & ' definition syntax requires that the start and end'
          WRITE(stderr,*) '   points of the spline be the first'//&
            & ' and last points of the list of points to move.'
          WRITE(stderr,*) '   Sorry about that.'
          WRITE(stderr,*)
          WRITE(stderr,*) '   Try something like'
          WRITE(stderr,*)
          WRITE(stderr,'(a)',advance='no') '      -redist '
          IF(node_points(1).ne.node_spline(1)) THEN
            WRITE(stderr,'(i3)',advance='no') node_spline(1)
          ENDIF
          DO k=1,npoints
            WRITE(stderr,'(i3)',advance='no') node_points(k)
          ENDDO
          IF(node_points(npoints).ne.node_spline(nspline)) THEN
            WRITE(stderr,'(i3)',advance='no') node_spline(nspline)
          ENDIF
          WRITE(stderr,'(a)',advance='no') ':'
          DO k=1,nspline
            WRITE(stderr,'(i3)',advance='no') node_spline(k)
          ENDDO
          STOP
        ENDIF
        !WRITE(*,*) 'coords: ',j,x_spline(j,:)
      ENDDO

      ! Redistribute them
      ! Needs to be modified for more than 1 point in a spline
      CALL InterpLinear( x_spline,x_points,nspline,npoints,NumCoord )
      !CALL InterpCubic( x_spline,x_points,nspline,npoints,NumCoord )

      ! Put back in the NodeValue array
      DO j=1,npoints
        DO k=1,NumNode
          IF (NodeIndex(k).eq.node_points(j)) GOTO 200
        ENDDO
        WRITE(stderr,*) 'RedistPoints[2]: unable to find node ',node_points(j)
        STOP
200     CONTINUE
        idx=k
        DO k=1,NumCoord
          DO l=1,NodeVersion(idx,k)
            NodeValue(idx,k,l,1)=x_points(j,k)
          ENDDO
        ENDDO
        !WRITE(*,*) 'points: ',j,x_points(j,:)
      ENDDO

      ! Clean up
      DEALLOCATE(node_points, &
        &        node_spline, &
        &        x_points,    &
        &        x_spline,    &
        &        stat=ierr)
    ENDDO

    RETURN
  END SUBROUTINE RedistPoints

  !----------------------------------------------------------------

  SUBROUTINE CountMentionsOfNodesInSlices( NumNode,NumSlice, &
    & NodeIndex,NodeSliceCount,NumNodeInSlice,SliceNodes,verbose  )

    ! Count the number of slices that each node appears in

    IMPLICIT none
    ! Arguments
    INTEGER(I4)::NumNode
    INTEGER(I4)::NumSlice
    ! Node Data
    INTEGER(I4)::NodeIndex(:)
    INTEGER(I4)::NodeSliceCount(:)
    ! Slice data
    INTEGER(I4)::NumNodeInSlice(:)
    INTEGER(I4)::SliceNodes(:,:)
    ! Other
    LOGICAL::Verbose
    ! Local variables
    INTEGER(I4)::i,j,k,node

    NodeSliceCount(:)=0

    DO i=1,NumSlice
      DO j=1,NumNodeInSlice(i)
        node=SliceNodes(i,j)

        DO k=1,NumNode
          IF(NodeIndex(k).eq.node) THEN
            NodeSliceCount(k)=NodeSliceCount(k)+1
            GOTO 100
          ENDIF
        ENDDO
        WRITE(stderr,*) 'CountMentionsOfNodesInSlices: Can''t find node',node
        STOP
100     CONTINUE
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE CountMentionsOfNodesInSlices

  !----------------------------------------------------------------

  SUBROUTINE FindNewProjectionPoint( NewNode,Node,NumNode,NumSlice, &
    & NodeIndex,NodeSliceCount,NumNodeInSlice,SliceNodes,Elemnode,  &
    & verbose  )

    ! Check to see if a node is in more than one slice. If
    ! it is then find the correct projection point.

    IMPLICIT none
    ! Arguments
    INTEGER(I4)::Node
    INTEGER(I4)::NumNode
    INTEGER(I4)::NumSlice
    INTEGER(I4)::NewNode
    ! Node Data
    INTEGER(I4)::NodeIndex(:)
    INTEGER(I4)::NodeSliceCount(:)
    ! Slice data
    INTEGER(I4)::NumNodeInSlice(:)
    INTEGER(I4)::SliceNodes(:,:)
    INTEGER(I4)::Elemnode(:)
    ! Other
    LOGICAL::Verbose
    ! Local variables
    INTEGER(I4)::i,j,k,l


    NewNode=Node

    ! Find the node, check if it is mentioned once only
    DO i=1,NumNode
      IF (NodeIndex(i).eq.Node) THEN
        GOTO 100
      ENDIF
    ENDDO
    WRITE(stderr,*) 'FindNewProjectionPoint: Can''t find node ',Node
    STOP
100 CONTINUE

    IF (NodeSliceCount(i).le.1) THEN
      RETURN
    ENDIF

    ! Otherwise, find slice that the node is in
    DO i=1,NumSlice
      DO j=1,NumNodeInSlice(i)
        IF (SliceNodes(i,j).eq.Node) THEN
          ! We are in this slice
          DO k=1,NumNodeInSlice(i)
            IF (k.ne.j) THEN
              ! Check over the list of nodes in the element
              DO l=1,NumNodePerElem_2D
                IF (Elemnode(l).eq.SliceNodes(i,k)) THEN
                  NewNode=Elemnode(l)
                  IF(verbose) THEN
                    WRITE(stderr,*) 'FindNewProjectionPoint: node ',node,' -> ',NewNode
                  ENDIF
                  RETURN
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    WRITE(stderr,*) 'FindNewProjectionPoint: Can''t find node in any slice',Node

    RETURN
  END SUBROUTINE FindNewProjectionPoint

  !----------------------------------------------------------------

END MODULE createnodes
